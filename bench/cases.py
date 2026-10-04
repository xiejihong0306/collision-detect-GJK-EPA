"""Named convex-pair scenes for the bench and regression suite."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .adapter import convex_hull_mesh, unique_vertices
from .geom import Pose, box_faces, box_vertices, compose_rotation, skinny_vertices, tetra_faces, tetra_vertices


@dataclass
class Body:
    name: str
    local_verts: np.ndarray
    local_faces: np.ndarray
    pose: Pose = field(default_factory=Pose)

    def world_verts(self) -> np.ndarray:
        return self.pose.apply(self.local_verts)

    def world_faces(self) -> np.ndarray:
        return np.array([self.pose.apply(tri) for tri in self.local_faces], dtype=np.float64)


@dataclass
class Case:
    name: str
    title: str
    body_a: Body
    body_b: Body
    expect: str
    expect_depth: float | None = None
    expect_normal_axis: int | None = None
    feature_tol: float = 1.0
    notes: str = ""


def _box(name: str, size=(1.0, 1.0, 1.0), t=(0.0, 0.0, 0.0), R=None) -> Body:
    pose = Pose(np.asarray(t, dtype=np.float64), np.eye(3) if R is None else np.asarray(R, dtype=np.float64))
    return Body(name, box_vertices(size), box_faces(size), pose)


def _tet(name: str, t=(0.0, 0.0, 0.0), R=None) -> Body:
    pose = Pose(np.asarray(t, dtype=np.float64), np.eye(3) if R is None else np.asarray(R, dtype=np.float64))
    return Body(name, tetra_vertices(), tetra_faces(), pose)


def _rotated(R_axis, deg: float) -> np.ndarray:
    return compose_rotation(np.eye(3), np.asarray(R_axis, dtype=np.float64), np.deg2rad(deg))


def catalog() -> dict[str, Case]:
    cases: list[Case] = [
        Case("separated", "Clearly separated cubes", _box("A"), _box("B", t=(2.0, 0.0, 0.0)), "separated"),
        Case("approach", "Cubes 0.05 apart", _box("A"), _box("B", t=(1.05, 0.0, 0.0)), "separated"),
        Case("face_touch", "Face contact", _box("A"), _box("B", t=(1.0, 0.0, 0.0)), "contact", 0.0, 0),
        Case("shallow", "Shallow overlap 0.05", _box("A"), _box("B", t=(0.95, 0.0, 0.0)), "penetrating", 0.05, 0),
        Case("overlap", "Overlap 0.5", _box("A"), _box("B", t=(0.5, 0.0, 0.0)), "penetrating", 0.5, 0),
        Case("deep", "Deep overlap 0.9", _box("A"), _box("B", t=(0.1, 0.0, 0.0)), "penetrating", 0.9, 0),
        Case(
            "contained",
            "Small cube inside large cube",
            _box("A", (2.0, 2.0, 2.0)),
            _box("B", (0.4, 0.4, 0.4), t=(0.8, 0.8, 0.8)),
            "penetrating",
        ),
        Case("coincident", "Identical coincident cubes", _box("A"), _box("B"), "penetrating", 1.0),
        Case(
            "vertex_face",
            "Tetra vertex into cube top",
            _box("A"),
            _tet("B", t=(0.3, 0.3, 0.7)),
            "penetrating",
        ),
        Case(
            "edge_edge",
            "Crossing-edge tetrahedra",
            Body(
                "A",
                (va := np.array([[0.0, -0.2, 0.0], [0.0, 0.2, 0.0], [-0.4, 0.0, -0.4], [-0.4, 0.0, 0.4]])),
                tetra_faces(va),
            ),
            Body(
                "B",
                (vb := np.array([[-0.2, 0.0, 0.05], [0.2, 0.0, 0.05], [0.0, -0.4, 0.45], [0.0, 0.4, 0.45]])),
                tetra_faces(vb),
            ),
            "penetrating",
        ),
        Case("face_face", "Matching cube faces", _box("A"), _box("B", t=(1.0, 0.2, 0.2)), "contact", 0.0, 0),
        Case(
            "parallel_faces",
            "Parallel-face boxes",
            _box("A", (1.2, 0.8, 0.6)),
            _box("B", (1.0, 0.7, 0.5), t=(1.2, 0.05, 0.05)),
            "contact",
            0.0,
            0,
        ),
        Case(
            "near_coplanar",
            "Nearly coplanar overlap",
            _box("A"),
            _box("B", t=(0.999, 0.0, 0.0)),
            "penetrating",
            0.001,
            0,
        ),
        Case("box", "Long box vs cube", _box("A"), _box("B", (2.0, 0.4, 0.6), t=(0.7, 0.2, 0.1)), "penetrating"),
        Case("tetra", "Two tetrahedra", _tet("A"), _tet("B", t=(0.25, 0.1, 0.1)), "penetrating"),
        Case(
            "skinny",
            "Thin rod through cube",
            _box("A"),
            Body("B", skinny_vertices(), box_faces((0.08, 0.08, 2.4)), Pose(np.array([0.46, 0.46, -0.7]))),
            "penetrating",
        ),
        Case(
            "rotated",
            "Rotated overlapping cubes",
            _box("A"),
            _box("B", t=(0.4, 0.1, 0.0), R=_rotated((0.0, 0.0, 1.0), 18.0)),
            "penetrating",
        ),
        Case("tiny_gap", "Gap 1e-5", _box("A"), _box("B", t=(1.0 + 1e-5, 0.0, 0.0)), "separated"),
        Case("tiny_pen", "Penetration 1e-5", _box("A"), _box("B", t=(1.0 - 1e-5, 0.0, 0.0)), "penetrating", 1e-5, 0),
        Case(
            "scale_mismatch",
            "Unit vs 20x cube",
            _box("A"),
            _box("B", (20.0, 20.0, 20.0), t=(0.5, -9.5, -9.5)),
            "penetrating",
        ),
        Case(
            "large_offset",
            "Small feature far from origin",
            _box("A", t=(1.0e4, 0.0, 0.0)),
            _box("B", (0.2, 0.2, 0.2), t=(1.0e4 + 0.85, 0.4, 0.4)),
            "penetrating",
        ),
        Case(
            "duplicates",
            "Cube with repeated vertices",
            Body("A", np.vstack([box_vertices(), box_vertices()[:3]]), box_faces(), Pose()),
            _box("B", t=(0.5, 0.0, 0.0)),
            "penetrating",
            0.5,
            0,
        ),
        Case(
            "degenerate",
            "Colinear 'tetra' (must fail or refuse)",
            Body(
                "A",
                np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                np.zeros((0, 3, 3)),
            ),
            _box("B"),
            "failed",
            notes="Input is not a volume. SAT is uncertain; algorithm must not silently look separated.",
        ),
    ]
    return {c.name: c for c in cases}


def default_case() -> Case:
    return catalog()["overlap"]


def random_convex(lib, rng: np.random.Generator, n: int = 12, scale: float = 1.0) -> Body:
    raw = rng.normal(size=(n, 3))
    raw = raw / np.maximum(np.linalg.norm(raw, axis=1, keepdims=True), 1e-12)
    raw *= scale * (0.4 + 0.6 * rng.random((n, 1)))
    mesh, info = convex_hull_mesh(lib, raw)
    if info != 0 or mesh.shape[0] < 4:
        # Fallback: a tetra that is always valid.
        return _tet("rand")
    verts = unique_vertices(mesh)
    return Body("rand", verts, mesh, Pose())


def random_pose(rng: np.random.Generator, translate_scale: float = 1.5) -> Pose:
    axis = rng.normal(size=3)
    ang = float(rng.uniform(0.0, 2.0 * np.pi))
    t = rng.uniform(-translate_scale, translate_scale, size=3)
    return Pose(t, compose_rotation(np.eye(3), axis, ang))
