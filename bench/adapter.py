"""ctypes adapter around the Fortran GJK-EPA / QuickHull C ABI.

Conventions (must match src/collision/gjkepa.f90):
    Minkowski support:  S(d) = S_A(d) - S_B(-d)
    EPA normal n is the outward normal of the closest face of A-B
    Witness:  A = support(A, n),  B = support(B, -n)
    To separate: translate B by sign(n·(A-B)) * n * depth
    depth >= 0; colli_type 0=separated, 1=penetrating, 2=face contact
    TOL_FF_ is a support-feature band for classification, not GJK epsilon
"""

from __future__ import annotations

import os
import threading
import time
from ctypes import (
    CDLL,
    POINTER,
    byref,
    c_double,
    c_int,
)
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

GJK_TRACE_MAX = 64
EPA_TRACE_MAX = 128
QUICKHULL_MAX_FACES = 4096

GJK_REASON = {
    1: "continue",
    2: "support_plane_miss",
    3: "duplicate_support",
    4: "origin_in_tetra",
    5: "origin_close",
    6: "max_iter",
}
EPA_REASON = {
    1: "expand",
    2: "converged",
    3: "no_expand",
    4: "max_iter",
    5: "fail",
}

_INFO = {
    0: "ok",
    1: "bad_input",
    2: "gjk_iter_limit",
    3: "epa_fail",
}

_HULL_LOCK = threading.Lock()
_LIB: CDLL | None = None
_LIB_PATH: Path | None = None


def _as_c_points(verts: np.ndarray) -> tuple[np.ndarray, int]:
    arr = np.ascontiguousarray(np.asarray(verts, dtype=np.float64))
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"expected (n,3) vertices, got {arr.shape}")
    return arr, int(arr.shape[0])


def _bind_signatures(lib: CDLL) -> None:
    lib.gjkepa_query_c.argtypes = [
        POINTER(c_double),
        c_int,
        POINTER(c_double),
        c_int,
        c_int,
        c_double,
        c_int,
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
    ]
    lib.quickhull_c.argtypes = [
        POINTER(c_double),
        c_int,
        c_int,
        POINTER(c_double),
        POINTER(c_int),
        POINTER(c_int),
    ]


def _span_degenerate(verts: np.ndarray) -> bool:
    if verts.shape[0] < 4:
        return True
    centered = verts - verts.mean(axis=0)
    try:
        singular = np.linalg.svd(centered, compute_uv=False)
    except np.linalg.LinAlgError:
        return True
    return float(singular[-1]) < 1e-8 * max(float(singular[0]), 1e-12)


def load_library(path: Path) -> CDLL:
    global _LIB, _LIB_PATH
    path = path.resolve()
    if _LIB is not None and _LIB_PATH == path:
        return _LIB
    lib_dir = str(path.parent)
    fortran_hint = os.environ.get("PATH", "")
    os.environ["PATH"] = lib_dir + os.pathsep + fortran_hint
    if os.name == "nt":
        os.add_dll_directory(lib_dir)
        # gfortran runtime often sits next to the compiler
        for extra in fortran_hint.split(os.pathsep):
            if extra and Path(extra).is_dir():
                try:
                    os.add_dll_directory(extra)
                except OSError:
                    pass
    lib = CDLL(str(path))
    lib.gjkepa_abi_version.restype = c_int
    lib.gjkepa_abi_version.argtypes = []
    ver = int(lib.gjkepa_abi_version())
    if ver != 1:
        raise RuntimeError(f"unexpected gjkepa C ABI {ver}")
    lib.gjkepa_query_c.restype = None
    lib.quickhull_c.restype = None
    _bind_signatures(lib)
    _LIB = lib
    _LIB_PATH = path
    return lib


@dataclass
class GjkStep:
    direction: np.ndarray
    support: np.ndarray
    closest: np.ndarray
    nsimp: int
    reason: str


@dataclass
class EpaStep:
    normal: np.ndarray
    dist: float
    support: np.ndarray
    nfaces: int
    nverts: int
    reason: str


@dataclass
class QueryResult:
    collision: bool
    colli_type: int
    info: int
    info_name: str
    gjk_iters: int
    epa_iters: int
    epa_info: int
    nearest_a: np.ndarray
    nearest_b: np.ndarray
    normal: np.ndarray
    point: np.ndarray
    depth: float
    algo_s: float
    transfer_s: float
    gjk_trace: list[GjkStep] = field(default_factory=list)
    epa_trace: list[EpaStep] = field(default_factory=list)
    status: str = "separated"
    warning: str = ""
    envelope_skip: bool = False

    def resolve_b_translation(self) -> np.ndarray:
        """World translation of B that this adapter treats as depenetration.

        EPA n is outward on the closest A-B face. Witnesses are
        A=support(A,n), B=support(B,-n). The sign is taken from
        n · (A-B) so we do not guess using body centroids.
        Moving B by +sign(n, A-B) * n * depth walks the pair toward contact.
        """
        if not self.collision or self.depth <= 0.0:
            return np.zeros(3)
        n = np.asarray(self.normal, dtype=np.float64)
        ln = float(np.linalg.norm(n))
        if ln < 1e-15:
            return np.zeros(3)
        n = n / ln
        sep = np.asarray(self.nearest_a, dtype=np.float64) - np.asarray(self.nearest_b, dtype=np.float64)
        sign = 1.0 if float(np.dot(n, sep)) >= 0.0 else -1.0
        return n * float(self.depth) * sign


def classify_status(result: QueryResult, scale: float, contact_rel: float = 1e-4) -> QueryResult:
    contact_eps = max(1e-9, contact_rel * max(scale, 1e-12))
    if result.info == 1:
        result.status = "failed"
        result.warning = "invalid input (need at least one (n,3) vertex array)"
        return result
    if result.info == 3:
        result.status = "failed"
        result.warning = "EPA failed to seed or expand; GJK hit is not treated as a clean miss"
        return result
    if result.info == 2:
        result.warning = "GJK reached the iteration cap; result may be incomplete"
    if result.gjk_iters == 0 and not result.collision and result.info == 0:
        result.envelope_skip = True
        result.status = "separated"
        return result
    if not result.collision:
        result.status = "separated"
        return result
    if result.colli_type == 2 or result.depth <= contact_eps:
        result.status = "contact"
        return result
    result.status = "penetrating"
    return result


def query(
    lib: CDLL,
    verts_a: np.ndarray,
    verts_b: np.ndarray,
    *,
    version: int = 2,
    tol_ff: float = 1.0,
    want_trace: bool = False,
    scale: float | None = None,
) -> QueryResult:
    """Call the Fortran algorithm. Vertices must already be in world coordinates."""
    a, n1 = _as_c_points(verts_a)
    b, n2 = _as_c_points(verts_b)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        res = QueryResult(
            collision=False,
            colli_type=0,
            info=1,
            info_name="bad_input",
            gjk_iters=0,
            epa_iters=0,
            epa_info=0,
            nearest_a=np.zeros(3),
            nearest_b=np.zeros(3),
            normal=np.zeros(3),
            point=np.zeros(3),
            depth=0.0,
            algo_s=0.0,
            transfer_s=0.0,
        )
        return classify_status(res, 1.0)

    collision = c_int()
    colli_type = c_int()
    info = c_int()
    gjk_iters = c_int()
    epa_iters = c_int()
    epa_info = c_int()
    nearest_a = (c_double * 3)()
    nearest_b = (c_double * 3)()
    normal = (c_double * 3)()
    point = (c_double * 3)()
    depth = c_double()
    elapsed = c_double()
    n_gjk = c_int()
    n_epa = c_int()
    gjk_dir = (c_double * (3 * GJK_TRACE_MAX))()
    gjk_support = (c_double * (3 * GJK_TRACE_MAX))()
    gjk_closest = (c_double * (3 * GJK_TRACE_MAX))()
    gjk_nsimp = (c_int * GJK_TRACE_MAX)()
    gjk_reason = (c_int * GJK_TRACE_MAX)()
    epa_nml = (c_double * (3 * EPA_TRACE_MAX))()
    epa_dist = (c_double * EPA_TRACE_MAX)()
    epa_support = (c_double * (3 * EPA_TRACE_MAX))()
    epa_nfaces = (c_int * EPA_TRACE_MAX)()
    epa_nverts = (c_int * EPA_TRACE_MAX)()
    epa_reason = (c_int * EPA_TRACE_MAX)()

    t0 = time.perf_counter()
    lib.gjkepa_query_c(
        a.ctypes.data_as(POINTER(c_double)),
        c_int(n1),
        b.ctypes.data_as(POINTER(c_double)),
        c_int(n2),
        c_int(version),
        c_double(tol_ff),
        c_int(1 if want_trace else 0),
        byref(collision),
        byref(colli_type),
        byref(info),
        byref(gjk_iters),
        byref(epa_iters),
        byref(epa_info),
        nearest_a,
        nearest_b,
        normal,
        point,
        byref(depth),
        byref(elapsed),
        gjk_dir,
        gjk_support,
        gjk_closest,
        gjk_nsimp,
        gjk_reason,
        byref(n_gjk),
        epa_nml,
        epa_dist,
        epa_support,
        epa_nfaces,
        epa_nverts,
        epa_reason,
        byref(n_epa),
    )
    transfer_s = time.perf_counter() - t0

    gtrace: list[GjkStep] = []
    etrace: list[EpaStep] = []
    if want_trace:
        for i in range(int(n_gjk.value)):
            gtrace.append(
                GjkStep(
                    direction=np.array([gjk_dir[3 * i], gjk_dir[3 * i + 1], gjk_dir[3 * i + 2]]),
                    support=np.array([gjk_support[3 * i], gjk_support[3 * i + 1], gjk_support[3 * i + 2]]),
                    closest=np.array([gjk_closest[3 * i], gjk_closest[3 * i + 1], gjk_closest[3 * i + 2]]),
                    nsimp=int(gjk_nsimp[i]),
                    reason=GJK_REASON.get(int(gjk_reason[i]), str(int(gjk_reason[i]))),
                )
            )
        for i in range(int(n_epa.value)):
            etrace.append(
                EpaStep(
                    normal=np.array([epa_nml[3 * i], epa_nml[3 * i + 1], epa_nml[3 * i + 2]]),
                    dist=float(epa_dist[i]),
                    support=np.array([epa_support[3 * i], epa_support[3 * i + 1], epa_support[3 * i + 2]]),
                    nfaces=int(epa_nfaces[i]),
                    nverts=int(epa_nverts[i]),
                    reason=EPA_REASON.get(int(epa_reason[i]), str(int(epa_reason[i]))),
                )
            )

    sc = float(scale) if scale is not None else float(
        max(np.ptp(a, axis=0).max(initial=0.0), np.ptp(b, axis=0).max(initial=0.0), 1e-12)
    )
    res = QueryResult(
        collision=bool(collision.value),
        colli_type=int(colli_type.value),
        info=int(info.value),
        info_name=_INFO.get(int(info.value), str(int(info.value))),
        gjk_iters=int(gjk_iters.value),
        epa_iters=int(epa_iters.value),
        epa_info=int(epa_info.value),
        nearest_a=np.array([nearest_a[0], nearest_a[1], nearest_a[2]], dtype=np.float64),
        nearest_b=np.array([nearest_b[0], nearest_b[1], nearest_b[2]], dtype=np.float64),
        normal=np.array([normal[0], normal[1], normal[2]], dtype=np.float64),
        point=np.array([point[0], point[1], point[2]], dtype=np.float64),
        depth=float(depth.value),
        algo_s=float(elapsed.value),
        transfer_s=transfer_s,
        gjk_trace=gtrace,
        epa_trace=etrace,
    )
    res = classify_status(res, sc)
    if _span_degenerate(a) or _span_degenerate(b):
        res.status = "failed"
        extra = "input is not a full-dimensional convex body"
        res.warning = f"{res.warning}; {extra}" if res.warning else extra
    return res


def convex_hull_mesh(lib: CDLL, points: np.ndarray) -> tuple[np.ndarray, int]:
    """Return (nfaces, 3, 3) triangles from the Fortran QuickHull, plus info."""
    pts, n = _as_c_points(points)
    faces = np.zeros((3, 3, QUICKHULL_MAX_FACES), dtype=np.float64, order="F")
    nfaces = c_int()
    info = c_int()
    with _HULL_LOCK:
        lib.quickhull_c(
            pts.ctypes.data_as(POINTER(c_double)),
            c_int(n),
            c_int(QUICKHULL_MAX_FACES),
            faces.ctypes.data_as(POINTER(c_double)),
            byref(nfaces),
            byref(info),
        )
    nf = int(nfaces.value)
    mesh = np.transpose(faces[:, :, :nf], (2, 1, 0)).copy() if nf > 0 else np.zeros((0, 3, 3))
    return mesh, int(info.value)


def unique_vertices(mesh: np.ndarray, atol: float = 1e-10) -> np.ndarray:
    if mesh.size == 0:
        return np.zeros((0, 3))
    pts = mesh.reshape(-1, 3)
    out: list[np.ndarray] = []
    for p in pts:
        if not any(np.linalg.norm(p - q) <= atol for q in out):
            out.append(p)
    return np.vstack(out) if out else np.zeros((0, 3))
