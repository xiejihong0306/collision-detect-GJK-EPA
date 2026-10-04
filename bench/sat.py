"""Independent SAT reference for convex polyhedra given triangle meshes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SatResult:
    intersect: bool
    uncertain: bool
    reason: str
    mtv: np.ndarray | None
    mtv_depth: float | None
    scale: float
    tol: float


def _unit(v: np.ndarray, eps: float) -> np.ndarray | None:
    n = float(np.linalg.norm(v))
    if n < eps:
        return None
    return v / n


def _edges_from_faces(faces: np.ndarray) -> np.ndarray:
    edges = []
    for tri in faces:
        for i in range(3):
            e = tri[(i + 1) % 3] - tri[i]
            edges.append(e)
    return np.asarray(edges, dtype=np.float64)


def _face_normals(faces: np.ndarray, eps: float) -> list[np.ndarray]:
    out: list[np.ndarray] = []
    for tri in faces:
        n = _unit(np.cross(tri[1] - tri[0], tri[2] - tri[0]), eps)
        if n is not None:
            out.append(n)
    return out


def _project(verts: np.ndarray, axis: np.ndarray) -> tuple[float, float]:
    d = verts @ axis
    return float(d.min()), float(d.max())


def _axis_seen(axis: np.ndarray, used: list[np.ndarray], eps: float) -> bool:
    for u in used:
        if abs(float(np.dot(u, axis))) > 1.0 - 1e-8:
            return True
    return False


def sat_check(
    verts_a: np.ndarray,
    faces_a: np.ndarray,
    verts_b: np.ndarray,
    faces_b: np.ndarray,
    *,
    rel_tol: float = 1e-6,
) -> SatResult:
    """Boolean SAT plus a translation MTV that treats containment correctly.

    Per axis the exit distance is min(Bmax-Amin, Amax-Bmin), not the raw
    interval-intersection length. Containment therefore uses the thinner
    remaining wall, not the contained body's projected width.
    """
    va = np.asarray(verts_a, dtype=np.float64)
    vb = np.asarray(verts_b, dtype=np.float64)
    fa = np.asarray(faces_a, dtype=np.float64)
    fb = np.asarray(faces_b, dtype=np.float64)
    if va.size == 0 or vb.size == 0 or fa.size == 0 or fb.size == 0:
        return SatResult(False, True, "empty_geometry", None, None, 0.0, 0.0)

    scale = float(max(np.linalg.norm(np.ptp(va, axis=0)), np.linalg.norm(np.ptp(vb, axis=0)), 1e-12))
    tol = max(1e-12, rel_tol * scale)
    eps_axis = max(1e-12, 1e-8 * scale)

    axes: list[np.ndarray] = []
    for n in _face_normals(fa, eps_axis) + _face_normals(fb, eps_axis):
        if not _axis_seen(n, axes, eps_axis):
            axes.append(n)

    ea = _edges_from_faces(fa)
    eb = _edges_from_faces(fb)
    skipped_cross = 0
    for u in ea:
        un = _unit(u, eps_axis)
        if un is None:
            skipped_cross += 1
            continue
        for v in eb:
            vn = _unit(v, eps_axis)
            if vn is None:
                skipped_cross += 1
                continue
            c = _unit(np.cross(un, vn), eps_axis)
            if c is None:
                skipped_cross += 1
                continue
            if not _axis_seen(c, axes, eps_axis):
                axes.append(c)

    best_depth = np.inf
    best_axis = None
    for axis in axes:
        amin, amax = _project(va, axis)
        bmin, bmax = _project(vb, axis)
        if amax < bmin - tol or bmax < amin - tol:
            return SatResult(False, False, "separated", None, None, scale, tol)
        # One-sided exits; valid for overlap and for A-contains-B / B-contains-A.
        exit_pos = amax - bmin
        exit_neg = bmax - amin
        if exit_pos <= exit_neg:
            depth = exit_pos
            direction = axis
        else:
            depth = exit_neg
            direction = -axis
        if depth < best_depth:
            best_depth = depth
            best_axis = direction

    uncertain = skipped_cross > 0 and not np.isfinite(best_depth)
    if best_axis is None:
        return SatResult(True, True, "no_usable_axis", None, None, scale, tol)
    return SatResult(True, uncertain, "overlap", best_axis, float(best_depth), scale, tol)
