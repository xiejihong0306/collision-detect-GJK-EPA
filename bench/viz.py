"""PyVista world view + optional Minkowski pane. Discrete queries only."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from ctypes import CDLL

from .adapter import QueryResult
from .cases import Case, catalog
from .geom import characteristic_scale, compose_rotation
from .snapshot import save_snapshot
from .verify import run_query, run_sat

STATUS_COLOR = {
    "separated": "#2ca02c",
    "contact": "#ffbf00",
    "penetrating": "#d62728",
    "failed": "#7b2cbf",
}


def _faces_to_poly(faces: np.ndarray):
    import pyvista as pv

    if faces.size == 0:
        return pv.PolyData(np.zeros((0, 3)))
    pts = np.asarray(faces, dtype=np.float64).reshape(-1, 3)
    n = faces.shape[0]
    cells = np.hstack([np.full((n, 1), 3), np.arange(3 * n).reshape(n, 3)]).astype(np.int64)
    return pv.PolyData(pts, cells)


def _hud(case: Case, q: QueryResult, sat, mode: str, constraint: str | None, paused: bool, traces: bool, draw_s: float) -> str:
    n = q.normal
    ntxt = "unavailable" if q.status == "failed" and q.info == 1 else f"({n[0]:+.4f},{n[1]:+.4f},{n[2]:+.4f})"
    sat_txt = "unavailable (no faces)" if sat is None else (
        f"{'overlap' if sat.intersect else 'separate'}" + (" uncertain" if sat.uncertain else "")
    )
    gjk = "envelope_skip" if q.envelope_skip else (f"iters={q.gjk_iters}" if q.gjk_iters else "unavailable")
    epa = f"iters={q.epa_iters} info={q.epa_info}" if q.collision else "not run"
    if q.info == 3:
        epa = f"FAIL info={q.epa_info}"
    t = case.body_b.pose.translation
    pause = "PAUSED" if paused else "LIVE"
    axis = constraint or "free"
    warn = q.warning or ""
    return (
        f"{case.name}  {pause}  mode={mode}  axis={axis}  traces={'on' if traces else 'off'}\n"
        f"status={q.status.upper()}  type={q.colli_type}  info={q.info_name}\n"
        f"GJK {gjk}   EPA {epa}\n"
        f"depth={q.depth:.6g}  normal={ntxt}\n"
        f"A={np.array2string(q.nearest_a, precision=4)}  B={np.array2string(q.nearest_b, precision=4)}\n"
        f"point={np.array2string(q.point, precision=4)}\n"
        f"SAT {sat_txt}   scale={characteristic_scale(case.body_a.world_verts(), case.body_b.world_verts()):.4g}  tol_ff={case.feature_tol}\n"
        f"time algo={q.algo_s*1e3:.3f}ms  call={q.transfer_s*1e3:.3f}ms  draw={draw_s*1e3:.3f}ms\n"
        f"B pose t=({t[0]:+.3f},{t[1]:+.3f},{t[2]:+.3f})\n"
        f"{warn}\n"
        f"C camera  T move  R rotate  X/Y/Z constrain  [/] case  0 reset\n"
        f"wheel: T/C nudge B, R spin B   P pause  N step  S save  D traces"
    )


class Workbench:
    def __init__(self, lib: CDLL, case: Case, *, root: Path, lib_path: Path, artifacts: Path):
        self.lib = lib
        self.root = root
        self.lib_path = lib_path
        self.artifacts = artifacts
        self.cases = catalog()
        self.names = list(self.cases)
        self.case = case
        self.initial_b = case.body_b.pose.copy()
        self.mode = "translate"
        self.constraint: str | None = None
        self.paused = False
        self.traces = False
        self.q: QueryResult | None = None
        self.sat = None
        self.scene_ver = 0
        self.draw_s = 0.0
        self._drag = False
        self._grab = np.zeros(3)
        self._plane_p = np.zeros(3)
        self._plane_n = np.array([0.0, 0.0, 1.0])
        self._mink_cam_ready = False
        self.plotter = None

    def _refresh_query(self) -> None:
        self.scene_ver += 1
        ver = self.scene_ver
        self.q = run_query(self.lib, self.case, want_trace=self.traces)
        if ver != self.scene_ver:
            return
        if self.case.body_a.local_faces.size and self.case.body_b.local_faces.size:
            self.sat = run_sat(self.case)
        else:
            self.sat = None

    def _minkowski_poly(self):
        import pyvista as pv

        a = self.case.body_a.world_verts()
        b = self.case.body_b.world_verts()
        pts = (a[:, None, :] - b[None, :, :]).reshape(-1, 3)
        if pts.shape[0] > 400:
            idx = np.linspace(0, pts.shape[0] - 1, 400).astype(int)
            pts = pts[idx]
        cloud = pv.PolyData(pts)
        try:
            return cloud.delaunay_3d().extract_surface()
        except Exception:
            return cloud

    def _sync_meshes(self) -> None:
        import pyvista as pv

        t0 = time.perf_counter()
        fa = self.case.body_a.world_faces()
        fb = self.case.body_b.world_faces()
        color = STATUS_COLOR.get(self.q.status if self.q else "separated", "#888888")
        self.plotter.subplot(0, 0)
        if fa.size:
            self.plotter.add_mesh(
                _faces_to_poly(fa),
                name="body_a",
                color="#4c78a8",
                opacity=0.35,
                show_edges=True,
                edge_color="black",
                reset_camera=False,
            )
        if fb.size:
            self.plotter.add_mesh(
                _faces_to_poly(fb),
                name="body_b",
                color=color,
                opacity=0.45,
                show_edges=True,
                edge_color="black",
                reset_camera=False,
            )
        if self.q and self.q.collision:
            n = self.q.normal
            ln = float(np.linalg.norm(n))
            if ln > 1e-12:
                arrow = pv.Arrow(
                    start=self.q.point,
                    direction=n / ln,
                    scale=max(0.25, self.q.depth if self.q.depth > 1e-6 else 0.25),
                )
                self.plotter.add_mesh(arrow, name="normal", color="magenta", reset_camera=False)
            line = pv.Line(self.q.nearest_a, self.q.nearest_b)
            self.plotter.add_mesh(line, name="witness", color="white", line_width=3, reset_camera=False)
            preview = self.case.body_b.world_faces() + self.q.resolve_b_translation()
            if preview.size:
                self.plotter.add_mesh(
                    _faces_to_poly(preview),
                    name="preview",
                    color="white",
                    opacity=0.12,
                    show_edges=True,
                    reset_camera=False,
                )
        else:
            for name in ("normal", "witness", "preview"):
                try:
                    self.plotter.remove_actor(name, render=False)
                except Exception:
                    pass

        self.plotter.subplot(0, 1)
        mink = self._minkowski_poly()
        if mink.n_points > 0:
            try:
                self.plotter.remove_actor("minkowski_hint", render=False)
            except Exception:
                pass
            self.plotter.add_mesh(mink, name="minkowski", color="#f58518", opacity=0.25, show_edges=True, reset_camera=False)
            if not getattr(self, "_mink_cam_ready", False):
                self.plotter.reset_camera()
                self._mink_cam_ready = True
        else:
            try:
                self.plotter.remove_actor("minkowski", render=False)
            except Exception:
                pass
            self.plotter.add_text("Minkowski A-B (unavailable)", name="minkowski_hint", font_size=10)
        if self.traces and self.q and self.q.gjk_trace:
            last = self.q.gjk_trace[-1]
            self.plotter.add_mesh(pv.Sphere(radius=0.03, center=(0, 0, 0)), name="origin", color="red", reset_camera=False)
            self.plotter.add_mesh(
                pv.Line(np.zeros(3), last.support),
                name="gjk_support",
                color="black",
                reset_camera=False,
            )
        self.plotter.subplot(0, 0)
        self.plotter.add_text(
            _hud(self.case, self.q, self.sat, self.mode, self.constraint, self.paused, self.traces, self.draw_s),
            name="hud",
            font_size=9,
            color=STATUS_COLOR.get(self.q.status if self.q else "separated", "white"),
        )
        self.draw_s = time.perf_counter() - t0

    def _set_case(self, name: str) -> None:
        self.case = catalog()[name]
        self.initial_b = self.case.body_b.pose.copy()
        self._mink_cam_ready = False
        self._refresh_query()
        self._sync_meshes()
        self.plotter.subplot(0, 0)
        self.plotter.reset_camera()

    def _cycle_case(self, step: int) -> None:
        i = self.names.index(self.case.name) if self.case.name in self.names else 0
        self._set_case(self.names[(i + step) % len(self.names)])

    def _ray(self, x: int, y: int):
        import vtk

        ren = self.plotter.renderers[0]
        picker = vtk.vtkWorldPointPicker()
        picker.Pick(x, y, 0, ren)
        p0 = np.array(picker.GetPickPosition(), dtype=np.float64)
        picker.Pick(x, y, 1, ren)
        p1 = np.array(picker.GetPickPosition(), dtype=np.float64)
        d = p1 - p0
        n = float(np.linalg.norm(d))
        if n < 1e-12:
            return p0, np.array(ren.GetActiveCamera().GetDirectionOfProjection())
        return p0, d / n

    def _intersect_plane(self, origin, direction):
        den = float(np.dot(direction, self._plane_n))
        if abs(den) < 1e-12:
            return None
        t = float(np.dot(self._plane_p - origin, self._plane_n) / den)
        return origin + t * direction

    def _on_press(self, obj, _evt) -> None:
        if self.mode == "camera":
            return
        x, y = obj.GetEventPosition()
        origin, direction = self._ray(x, y)
        cam = self.plotter.renderers[0].GetActiveCamera()
        self._plane_n = np.array(cam.GetViewPlaneNormal(), dtype=np.float64)
        self._plane_p = self.case.body_b.pose.translation.copy()
        hit = self._intersect_plane(origin, direction)
        if hit is None:
            return
        self._drag = True
        self._grab = hit - self.case.body_b.pose.translation

    def _on_move(self, obj, _evt) -> None:
        if not self._drag:
            return
        x, y = obj.GetEventPosition()
        origin, direction = self._ray(x, y)
        hit = self._intersect_plane(origin, direction)
        if hit is None:
            return
        if self.mode == "translate":
            delta = hit - self._grab - self.case.body_b.pose.translation
            if self.constraint == "x":
                delta = np.array([delta[0], 0.0, 0.0])
            elif self.constraint == "y":
                delta = np.array([0.0, delta[1], 0.0])
            elif self.constraint == "z":
                delta = np.array([0.0, 0.0, delta[2]])
            self.case.body_b.pose.translation = self.case.body_b.pose.translation + delta
        elif self.mode == "rotate":
            axis = {"x": np.array([1.0, 0.0, 0.0]), "y": np.array([0.0, 1.0, 0.0]), "z": np.array([0.0, 0.0, 1.0])}.get(
                self.constraint, self._plane_n
            )
            self.case.body_b.pose.rotation = compose_rotation(self.case.body_b.pose.rotation, axis, 0.02 * (hit[0] - self._grab[0] - self.case.body_b.pose.translation[0]))
        if not self.paused:
            self._refresh_query()
        self._sync_meshes()
        self.plotter.render()

    def _on_release(self, _obj, _evt) -> None:
        self._drag = False

    def _wheel_axis(self) -> np.ndarray:
        cam = self.plotter.renderers[0].GetActiveCamera()
        named = {
            "x": np.array([1.0, 0.0, 0.0]),
            "y": np.array([0.0, 1.0, 0.0]),
            "z": np.array([0.0, 0.0, 1.0]),
        }
        if self.constraint in named:
            return named[self.constraint]
        if self.mode == "rotate":
            n = np.array(cam.GetViewPlaneNormal(), dtype=np.float64)
        else:
            n = -np.array(cam.GetDirectionOfProjection(), dtype=np.float64)
        ln = float(np.linalg.norm(n))
        return n / ln if ln > 1e-12 else np.array([0.0, 0.0, 1.0])

    def _nudge_pose(self, sign: int) -> None:
        axis = self._wheel_axis()
        scale = characteristic_scale(self.case.body_a.world_verts(), self.case.body_b.world_verts())
        if self.mode == "rotate":
            self.case.body_b.pose.rotation = compose_rotation(
                self.case.body_b.pose.rotation, axis, sign * np.deg2rad(8.0)
            )
        else:
            self.case.body_b.pose.translation = self.case.body_b.pose.translation + axis * sign * 0.04 * max(scale, 0.25)
        if not self.paused:
            self._refresh_query()
        self._sync_meshes()
        self.plotter.render()

    def _on_wheel_forward(self, _obj, _evt) -> None:
        self._nudge_pose(+1)

    def _on_wheel_backward(self, _obj, _evt) -> None:
        self._nudge_pose(-1)

    def run(self) -> None:
        import pyvista as pv

        self._refresh_query()
        self.plotter = pv.Plotter(shape=(1, 2), title="GJK-EPA verification bench")
        self.plotter.subplot(0, 0)
        self.plotter.add_axes()
        self.plotter.show_bounds(grid="back", location="outer", color="gray")
        self.plotter.add_text("world", font_size=10, position="upper_left")
        self.plotter.subplot(0, 1)
        self.plotter.add_axes()
        self.plotter.add_text("Minkowski A-B (debug hull / traces)", font_size=10)
        self.plotter.subplot(0, 0)
        self._sync_meshes()
        self.plotter.reset_camera()

        iren = self.plotter.iren.interactor
        style = iren.GetInteractorStyle()
        style.OnMouseWheelForward = lambda *a: None
        style.OnMouseWheelBackward = lambda *a: None
        iren.AddObserver("LeftButtonPressEvent", self._on_press, 1.0)
        iren.AddObserver("MouseMoveEvent", self._on_move, 1.0)
        iren.AddObserver("LeftButtonReleaseEvent", self._on_release, 1.0)
        iren.AddObserver("MouseWheelForwardEvent", self._on_wheel_forward, 1.0)
        iren.AddObserver("MouseWheelBackwardEvent", self._on_wheel_backward, 1.0)

        def bind(key, fn):
            self.plotter.add_key_event(key, fn)

        bind("c", lambda: setattr(self, "mode", "camera") or self._sync_meshes() or self.plotter.render())
        bind("t", lambda: setattr(self, "mode", "translate") or self._sync_meshes() or self.plotter.render())
        bind("r", lambda: setattr(self, "mode", "rotate") or self._sync_meshes() or self.plotter.render())
        bind("x", lambda: setattr(self, "constraint", None if self.constraint == "x" else "x") or self._sync_meshes() or self.plotter.render())
        bind("y", lambda: setattr(self, "constraint", None if self.constraint == "y" else "y") or self._sync_meshes() or self.plotter.render())
        bind("z", lambda: setattr(self, "constraint", None if self.constraint == "z" else "z") or self._sync_meshes() or self.plotter.render())
        bind("bracketright", lambda: self._cycle_case(1))
        bind("bracketleft", lambda: self._cycle_case(-1))
        bind("0", lambda: (setattr(self.case.body_b, "pose", self.initial_b.copy()), self._refresh_query(), self._sync_meshes(), self.plotter.render()))
        bind("p", lambda: setattr(self, "paused", not self.paused) or self._sync_meshes() or self.plotter.render())
        bind("n", lambda: (self._refresh_query(), self._sync_meshes(), self.plotter.render()))
        bind("d", lambda: (setattr(self, "traces", not self.traces), self._refresh_query(), self._sync_meshes(), self.plotter.render()))

        def _save():
            path = self.artifacts / f"scene_{self.case.name}.json"
            save_snapshot(
                path,
                root=self.root,
                case=self.case,
                query=self.q,
                sat=self.sat,
                seed=None,
                lib_path=self.lib_path,
                feature_tol=self.case.feature_tol,
            )
            print(f"saved {path}")

        bind("s", _save)
        print("Discrete collision queries only — fast dragging can skip intermediate contacts.")
        self.plotter.show()
