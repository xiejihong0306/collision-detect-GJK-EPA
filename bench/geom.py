"""Rigid poses and primitive convex meshes. No collision math."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def rotation_matrix(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    k = np.asarray(axis, dtype=np.float64)
    n = float(np.linalg.norm(k))
    if n < 1e-15:
        return np.eye(3)
    k = k / n
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    kx, ky, kz = k
    K = np.array([[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]])
    return np.eye(3) * c + (1.0 - c) * np.outer(k, k) + s * K


def compose_rotation(R: np.ndarray, axis: np.ndarray, angle_rad: float) -> np.ndarray:
    return rotation_matrix(axis, angle_rad) @ np.asarray(R, dtype=np.float64)


@dataclass
class Pose:
    translation: np.ndarray = field(default_factory=lambda: np.zeros(3))
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))

    def copy(self) -> Pose:
        return Pose(self.translation.copy(), self.rotation.copy())

    def apply(self, local: np.ndarray) -> np.ndarray:
        pts = np.asarray(local, dtype=np.float64)
        return pts @ self.rotation.T + self.translation

    def as_dict(self) -> dict:
        return {
            "translation": self.translation.tolist(),
            "rotation": self.rotation.tolist(),
        }

    @staticmethod
    def from_dict(data: dict) -> Pose:
        return Pose(
            translation=np.asarray(data["translation"], dtype=np.float64),
            rotation=np.asarray(data["rotation"], dtype=np.float64),
        )


def box_vertices(size: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> np.ndarray:
    sx, sy, sz = size
    pts = []
    for z in (0.0, sz):
        for y in (0.0, sy):
            for x in (0.0, sx):
                pts.append([x, y, z])
    return np.asarray(pts, dtype=np.float64)


def tetra_vertices() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )


def skinny_vertices() -> np.ndarray:
    return box_vertices((0.08, 0.08, 2.4))


def characteristic_scale(verts_a: np.ndarray, verts_b: np.ndarray) -> float:
    def ext(v: np.ndarray) -> float:
        if v.size == 0:
            return 0.0
        return float(np.linalg.norm(np.ptp(v, axis=0)))

    return max(ext(verts_a), ext(verts_b), 1e-12)


def box_faces(size: tuple[float, float, float] = (1.0, 1.0, 1.0)) -> np.ndarray:
    """12 triangles for an axis-aligned box with min corner at origin."""
    v = box_vertices(size)
    # vertex index: ix + 2*iy + 4*iz
    tris = [
        (0, 2, 3),
        (0, 3, 1),
        (4, 5, 7),
        (4, 7, 6),
        (0, 1, 5),
        (0, 5, 4),
        (2, 6, 7),
        (2, 7, 3),
        (0, 4, 6),
        (0, 6, 2),
        (1, 3, 7),
        (1, 7, 5),
    ]
    return np.array([[v[i], v[j], v[k]] for i, j, k in tris], dtype=np.float64)


def tetra_faces(verts: np.ndarray | None = None) -> np.ndarray:
    v = tetra_vertices() if verts is None else np.asarray(verts, dtype=np.float64)
    tris = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    return np.array([[v[i], v[j], v[k]] for i, j, k in tris], dtype=np.float64)
