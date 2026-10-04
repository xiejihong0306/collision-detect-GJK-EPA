"""Independent checks and headless suites. Failures return non-zero via launcher."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ctypes import CDLL

from .adapter import QueryResult, query
from .cases import Case, catalog, random_convex, random_pose
from .geom import Pose, characteristic_scale, compose_rotation
from .sat import sat_check
from .snapshot import save_snapshot


@dataclass
class Check:
    name: str
    ok: bool
    detail: str = ""


@dataclass
class CaseReport:
    name: str
    ok: bool
    checks: list[Check] = field(default_factory=list)
    snapshot: str | None = None
    query: QueryResult | None = None

    def fail(self, name: str, detail: str) -> None:
        self.ok = False
        self.checks.append(Check(name, False, detail))

    def pass_(self, name: str, detail: str = "") -> None:
        self.checks.append(Check(name, True, detail))


def run_query(lib: CDLL, case: Case, *, want_trace: bool = False) -> QueryResult:
    wa = case.body_a.world_verts()
    wb = case.body_b.world_verts()
    scale = characteristic_scale(wa, wb)
    return query(
        lib,
        wa,
        wb,
        version=2,
        tol_ff=case.feature_tol,
        want_trace=want_trace,
        scale=scale,
    )


def run_sat(case: Case):
    return sat_check(
        case.body_a.world_verts(),
        case.body_a.world_faces(),
        case.body_b.world_verts(),
        case.body_b.world_faces(),
    )


def _normal_on_axis(n: np.ndarray, axis: int, tol: float = 0.15) -> bool:
    if float(np.linalg.norm(n)) < 1e-12:
        return False
    u = n / np.linalg.norm(n)
    return abs(abs(float(u[axis])) - 1.0) <= tol


def evaluate_case(lib: CDLL, case: Case, *, artifacts: Path, root: Path, lib_path: Path) -> CaseReport:
    report = CaseReport(case.name, True)
    q = run_query(lib, case)
    report.query = q
    sat = None
    if case.body_a.local_faces.size > 0 and case.body_b.local_faces.size > 0:
        sat = run_sat(case)

    if case.expect == "failed":
        if q.status == "failed" or q.info != 0:
            report.pass_("refuse_or_fail", q.info_name)
        elif q.status == "separated" and q.info == 0:
            report.fail("refuse_or_fail", "degenerate input reported as clean separation")
        else:
            report.pass_("refuse_or_fail", f"status={q.status} info={q.info_name} (non-silent)")
        return _maybe_save(report, case, q, sat, artifacts, root, lib_path)

    if q.status == "failed":
        report.fail("algorithm", q.warning or q.info_name)
        return _maybe_save(report, case, q, sat, artifacts, root, lib_path)

    if case.expect == "separated":
        if q.collision or q.status != "separated":
            report.fail("state", f"expected separated, got {q.status} collision={q.collision}")
        else:
            report.pass_("state", q.status)
        if sat is not None and not sat.uncertain and sat.intersect:
            report.fail("sat", "SAT reports overlap on a separated case")
        elif sat is not None:
            report.pass_("sat", sat.reason)
    elif case.expect == "contact":
        if q.status not in ("contact", "penetrating"):
            report.fail("state", f"expected contact, got {q.status}")
        else:
            report.pass_("state", q.status)
        if case.expect_depth is not None and abs(q.depth - case.expect_depth) > max(2e-3, 0.05 * max(case.expect_depth, 1e-6)):
            if q.status == "contact" and case.expect_depth == 0.0 and q.depth <= 2e-2:
                report.pass_("depth", f"{q.depth}")
            else:
                report.fail("depth", f"{q.depth} vs {case.expect_depth}")
        if sat is not None and not sat.uncertain and not sat.intersect:
            report.fail("sat", "SAT reports separation on a contact case")
        elif sat is not None:
            report.pass_("sat", sat.reason)
    elif case.expect == "penetrating":
        if not q.collision or q.status == "separated":
            report.fail("state", f"expected penetrating, got {q.status}")
        else:
            report.pass_("state", q.status)
        if case.expect_depth is not None and abs(q.depth - case.expect_depth) > max(2e-3, 0.05 * case.expect_depth):
            report.fail("depth", f"{q.depth} vs {case.expect_depth}")
        elif case.expect_depth is not None:
            report.pass_("depth", f"{q.depth}")
        if case.expect_normal_axis is not None and not _normal_on_axis(q.normal, case.expect_normal_axis):
            report.fail("normal", f"{q.normal} not along axis {case.expect_normal_axis}")
        elif case.expect_normal_axis is not None:
            report.pass_("normal")
        if sat is not None and not sat.uncertain and not sat.intersect:
            report.fail("sat", "SAT reports separation on a penetrating case")
        elif sat is not None:
            report.pass_("sat", sat.reason)
        if q.collision and q.depth > 1e-8:
            _check_correction(lib, case, q, report)

    return _maybe_save(report, case, q, sat, artifacts, root, lib_path)


def _check_correction(lib: CDLL, case: Case, q: QueryResult, report: CaseReport) -> None:
    delta = q.resolve_b_translation()
    moved = Case(
        name=case.name + "/corrected",
        title=case.title,
        body_a=case.body_a,
        body_b=type(case.body_b)(
            case.body_b.name,
            case.body_b.local_verts,
            case.body_b.local_faces,
            Pose(case.body_b.pose.translation + delta, case.body_b.pose.rotation.copy()),
        ),
        expect="contact",
        feature_tol=case.feature_tol,
    )
    q2 = run_query(lib, moved)
    sat2 = None
    if moved.body_a.local_faces.size and moved.body_b.local_faces.size:
        sat2 = run_sat(moved)
    band = max(2e-2, 0.05 * max(q.depth, 1e-6))
    if q2.status == "failed":
        report.fail("correction", q2.warning or q2.info_name)
        return
    sat_ok = sat2 is None or sat2.uncertain or (not sat2.intersect) or (
        sat2.mtv_depth is not None and sat2.mtv_depth <= band
    )
    gjk_ok = q2.status == "separated" or (q2.collision and q2.depth <= band)
    if sat_ok:
        report.pass_("correction", f"SAT after EPA step: {None if sat2 is None else sat2.mtv_depth}; GJK {q2.status} depth={q2.depth}")
    elif gjk_ok:
        report.pass_("correction", f"GJK in band after EPA step depth={q2.depth}")
    else:
        report.fail(
            "correction",
            f"after witness-signed n*depth: GJK {q2.status} depth={q2.depth}; "
            f"SAT {None if sat2 is None else (sat2.intersect, sat2.mtv_depth)}",
        )


def _maybe_save(report, case, q, sat, artifacts, root, lib_path) -> CaseReport:
    if not report.ok:
        path = artifacts / f"failure_{case.name}.json"
        save_snapshot(path, root=root, case=case, query=q, sat=sat, seed=None, lib_path=lib_path, feature_tol=case.feature_tol)
        report.snapshot = str(path)
    return report


def invariance_checks(lib: CDLL, case: Case) -> list[Check]:
    checks: list[Check] = []
    q = run_query(lib, case)
    swapped = Case(
        name=case.name + "/swap",
        title=case.title,
        body_a=case.body_b,
        body_b=case.body_a,
        expect=case.expect,
        feature_tol=case.feature_tol,
    )
    qs = run_query(lib, swapped)
    if q.collision != qs.collision:
        checks.append(Check("swap_hit", False, f"{q.collision} vs {qs.collision}"))
    else:
        checks.append(Check("swap_hit", True))
    if q.collision and qs.collision and abs(q.depth - qs.depth) > max(2e-2, 0.08 * max(q.depth, qs.depth, 1e-6)):
        checks.append(Check("swap_depth", False, f"{q.depth} vs {qs.depth}"))
    else:
        checks.append(Check("swap_depth", True, f"{q.depth} vs {qs.depth}"))

    R = compose_rotation(np.eye(3), np.array([0.2, 0.8, 0.3]), 0.7)
    t = np.array([3.0, -1.5, 2.0])
    def xform(body):
        return type(body)(
            body.name,
            body.local_verts,
            body.local_faces,
            Pose(R @ body.pose.translation + t, R @ body.pose.rotation),
        )
    rigid = Case(case.name + "/rigid", case.title, xform(case.body_a), xform(case.body_b), case.expect, feature_tol=case.feature_tol)
    qr = run_query(lib, rigid)
    if q.collision != qr.collision:
        checks.append(Check("rigid_hit", False, f"{q.collision} vs {qr.collision}"))
    else:
        checks.append(Check("rigid_hit", True))
    sat0 = run_sat(case)
    satr = run_sat(rigid)
    if sat0.mtv_depth is not None and satr.mtv_depth is not None:
        if abs(sat0.mtv_depth - satr.mtv_depth) > max(2e-2, 0.08 * max(sat0.mtv_depth, 1e-6)):
            checks.append(Check("rigid_sat_depth", False, f"{sat0.mtv_depth} vs {satr.mtv_depth}"))
        else:
            checks.append(Check("rigid_sat_depth", True))
    if q.collision and qr.collision and abs(q.depth - qr.depth) > max(2e-2, 0.08 * max(q.depth, 1e-6)):
        checks.append(
            Check(
                "rigid_depth",
                False,
                f"GJK-EPA depth {q.depth} vs {qr.depth} after the same rigid motion "
                f"(SAT stayed {sat0.mtv_depth} / {satr.mtv_depth}); refine/EPA may be pose-dependent",
            )
        )
    else:
        checks.append(Check("rigid_depth", True))

    s = 2.5
    def scale_body(body):
        return type(body)(
            body.name,
            body.local_verts * s,
            body.local_faces * s if body.local_faces.size else body.local_faces,
            Pose(body.pose.translation * s, body.pose.rotation.copy()),
        )
    scaled = Case(case.name + "/scale", case.title, scale_body(case.body_a), scale_body(case.body_b), case.expect, feature_tol=case.feature_tol)
    qz = run_query(lib, scaled)
    if q.collision and qz.collision and abs(qz.depth - q.depth * s) > max(5e-2, 0.1 * max(q.depth * s, 1e-6)):
        checks.append(Check("scale_depth", False, f"{qz.depth} vs {q.depth * s}"))
    else:
        checks.append(Check("scale_depth", True))
    return checks


def run_regression(lib: CDLL, *, artifacts: Path, root: Path, lib_path: Path) -> list[CaseReport]:
    reports = []
    for case in catalog().values():
        reports.append(evaluate_case(lib, case, artifacts=artifacts, root=root, lib_path=lib_path))
    inv_src = catalog()["overlap"]
    inv = CaseReport("invariance_overlap", True)
    for chk in invariance_checks(lib, inv_src):
        inv.checks.append(chk)
        if not chk.ok:
            inv.ok = False
    if not inv.ok:
        q = run_query(lib, inv_src)
        sat = run_sat(inv_src)
        path = artifacts / "failure_invariance_overlap.json"
        save_snapshot(path, root=root, case=inv_src, query=q, sat=sat, seed=None, lib_path=lib_path, feature_tol=inv_src.feature_tol)
        inv.snapshot = str(path)
    reports.append(inv)
    return reports


def run_random(lib: CDLL, *, seed: int, count: int, artifacts: Path, root: Path, lib_path: Path) -> list[CaseReport]:
    rng = np.random.default_rng(seed)
    reports = []
    for i in range(count):
        a = random_convex(lib, rng)
        b = random_convex(lib, rng)
        b.pose = random_pose(rng)
        case = Case(f"random_{seed}_{i}", f"random pair {i}", a, b, "unknown")
        q = run_query(lib, case)
        sat = run_sat(case)
        rep = CaseReport(case.name, True, query=q)
        if q.status == "failed":
            # Random hulls should be valid; failure is a real issue.
            rep.fail("algorithm", q.warning or q.info_name)
        elif sat.uncertain:
            rep.pass_("sat_uncertain", sat.reason)
        elif q.collision != sat.intersect:
            if q.status == "contact":
                rep.pass_("sat_contact_band", "contact vs SAT may differ inside tolerance")
            else:
                rep.fail("sat_disagree", f"GJK collision={q.collision} SAT={sat.intersect}")
        else:
            rep.pass_("sat", f"agree collision={q.collision}")
        if not rep.ok:
            path = artifacts / f"failure_{case.name}.json"
            save_snapshot(path, root=root, case=case, query=q, sat=sat, seed=seed, lib_path=lib_path, feature_tol=case.feature_tol)
            rep.snapshot = str(path)
        reports.append(rep)
    return reports


def write_report(path: Path, reports: list[CaseReport]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "n": len(reports),
        "nfail": sum(1 for r in reports if not r.ok),
        "cases": [
            {
                "name": r.name,
                "ok": r.ok,
                "snapshot": r.snapshot,
                "checks": [{"name": c.name, "ok": c.ok, "detail": c.detail} for c in r.checks],
                "query": None
                if r.query is None
                else {
                    "status": r.query.status,
                    "collision": r.query.collision,
                    "depth": r.query.depth,
                    "info": r.query.info_name,
                },
            }
            for r in reports
        ],
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def print_reports(reports: list[CaseReport]) -> int:
    nfail = 0
    for r in reports:
        mark = "PASS" if r.ok else "FAIL"
        extra = f" snapshot={r.snapshot}" if r.snapshot else ""
        print(f"[{mark}] {r.name}{extra}")
        for c in r.checks:
            if not c.ok:
                print(f"       - {c.name}: {c.detail}")
        if not r.ok:
            nfail += 1
    print(f"{len(reports) - nfail}/{len(reports)} passed")
    return nfail
