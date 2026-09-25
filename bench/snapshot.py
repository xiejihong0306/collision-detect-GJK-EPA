"""Save and load failure / scene snapshots with actual geometry."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import __version__
from .adapter import QueryResult
from .cases import Body, Case
from .geom import Pose
from .sat import SatResult


def _git_rev(root: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
        )
        return out.stdout.strip() if out.returncode == 0 else "unavailable"
    except OSError:
        return "unavailable"


def _arr(x) -> list:
    return np.asarray(x, dtype=np.float64).tolist()


def dump_query(q: QueryResult) -> dict:
    return {
        "collision": q.collision,
        "colli_type": q.colli_type,
        "info": q.info,
        "info_name": q.info_name,
        "status": q.status,
        "warning": q.warning,
        "gjk_iters": q.gjk_iters,
        "epa_iters": q.epa_iters,
        "epa_info": q.epa_info,
        "nearest_a": _arr(q.nearest_a),
        "nearest_b": _arr(q.nearest_b),
        "normal": _arr(q.normal),
        "point": _arr(q.point),
        "depth": q.depth,
        "algo_s": q.algo_s,
        "transfer_s": q.transfer_s,
        "envelope_skip": q.envelope_skip,
        "gjk_trace": [
            {
                "direction": _arr(s.direction),
                "support": _arr(s.support),
                "closest": _arr(s.closest),
                "nsimp": s.nsimp,
                "reason": s.reason,
            }
            for s in q.gjk_trace
        ],
        "epa_trace": [
            {
                "normal": _arr(s.normal),
                "dist": s.dist,
                "support": _arr(s.support),
                "nfaces": s.nfaces,
                "nverts": s.nverts,
                "reason": s.reason,
            }
            for s in q.epa_trace
        ],
    }


def dump_sat(s: SatResult | None) -> dict | None:
    if s is None:
        return None
    return {
        "intersect": s.intersect,
        "uncertain": s.uncertain,
        "reason": s.reason,
        "mtv": None if s.mtv is None else _arr(s.mtv),
        "mtv_depth": s.mtv_depth,
        "scale": s.scale,
        "tol": s.tol,
    }


def save_snapshot(
    path: Path,
    *,
    root: Path,
    case: Case,
    query: QueryResult,
    sat: SatResult | None,
    seed: int | None,
    lib_path: Path,
    feature_tol: float,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format": "gjkepa-bench-snapshot-v1",
        "saved_at": datetime.now(timezone.utc).isoformat(),
        "bench_version": __version__,
        "git": _git_rev(root),
        "library": str(lib_path),
        "seed": seed,
        "feature_tol": feature_tol,
        "case": {
            "name": case.name,
            "title": case.title,
            "expect": case.expect,
            "notes": case.notes,
            "a": {
                "name": case.body_a.name,
                "local_verts": _arr(case.body_a.local_verts),
                "local_faces": _arr(case.body_a.local_faces),
                "pose": case.body_a.pose.as_dict(),
            },
            "b": {
                "name": case.body_b.name,
                "local_verts": _arr(case.body_b.local_verts),
                "local_faces": _arr(case.body_b.local_faces),
                "pose": case.body_b.pose.as_dict(),
            },
        },
        "algorithm": dump_query(query),
        "sat": dump_sat(sat),
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def load_snapshot(path: Path) -> Case:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("format") != "gjkepa-bench-snapshot-v1":
        raise ValueError(f"unsupported snapshot format: {data.get('format')}")
    c = data["case"]

    def body(key: str) -> Body:
        d = c[key]
        return Body(
            name=d["name"],
            local_verts=np.asarray(d["local_verts"], dtype=np.float64),
            local_faces=np.asarray(d["local_faces"], dtype=np.float64),
            pose=Pose.from_dict(d["pose"]),
        )

    return Case(
        name=c.get("name", path.stem),
        title=c.get("title", path.stem),
        body_a=body("a"),
        body_b=body("b"),
        expect=c.get("expect", "unknown"),
        feature_tol=float(data.get("feature_tol", 1.0)),
        notes=c.get("notes", ""),
    )
