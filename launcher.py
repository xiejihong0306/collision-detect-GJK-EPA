#!/usr/bin/env python3
"""Sole user entry for the GJK-EPA verification bench."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Bootstrap uses only the standard library so it can create .venv first.
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from bench.bootstrap import prepare  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="launcher.py",
        description="3D interactive / headless verification bench for this repository's GJK-EPA.",
    )
    p.add_argument("--case", help="named case (see --list-cases)")
    p.add_argument("--list-cases", action="store_true", help="print case names and exit")
    p.add_argument("--headless", action="store_true", help="no render window; run a test suite")
    p.add_argument("--suite", choices=("regression", "random", "all"), default="regression")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--count", type=int, default=100)
    p.add_argument("--replay", type=Path, help="snapshot JSON written by this bench")
    p.add_argument("--traces", action="store_true", help="request GJK/EPA iteration traces")
    p.add_argument("--artifacts", type=Path, default=ROOT / "artifacts")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    lib_path = prepare(ROOT, headless=args.headless or args.list_cases)
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true" if args.headless else "false")

    from bench.adapter import load_library
    from bench.cases import catalog, default_case
    from bench.snapshot import load_snapshot
    from bench.verify import print_reports, run_query, run_random, run_regression, write_report

    lib = load_library(lib_path)
    artifacts = Path(args.artifacts)
    artifacts.mkdir(parents=True, exist_ok=True)

    if args.list_cases:
        for name, case in catalog().items():
            print(f"{name:16}  {case.title}  expect={case.expect}")
        return 0

    if args.replay is not None:
        case = load_snapshot(Path(args.replay))
        q = run_query(lib, case, want_trace=args.traces)
        print(f"replay {args.replay}")
        print(f"  status={q.status} collision={q.collision} type={q.colli_type} info={q.info_name}")
        print(f"  depth={q.depth} normal={q.normal} iters GJK={q.gjk_iters} EPA={q.epa_iters}")
        if q.warning:
            print(f"  warning={q.warning}")
        if args.headless:
            return 0 if q.status != "failed" or case.expect == "failed" else 1
        from bench.viz import Workbench

        Workbench(lib, case, root=ROOT, lib_path=lib_path, artifacts=artifacts).run()
        return 0

    if args.headless:
        reports = []
        if args.suite in ("regression", "all"):
            reports.extend(run_regression(lib, artifacts=artifacts, root=ROOT, lib_path=lib_path))
        if args.suite in ("random", "all"):
            reports.extend(
                run_random(
                    lib,
                    seed=args.seed,
                    count=args.count,
                    artifacts=artifacts,
                    root=ROOT,
                    lib_path=lib_path,
                )
            )
        write_report(artifacts / "last_report.json", reports)
        nfail = print_reports(reports)
        print(f"report: {artifacts / 'last_report.json'}")
        return 1 if nfail else 0

    cases = catalog()
    if args.case:
        if args.case not in cases:
            print(f"unknown case {args.case!r}. Use --list-cases.", file=sys.stderr)
            return 2
        case = cases[args.case]
    else:
        case = default_case()
    from bench.viz import Workbench

    Workbench(lib, case, root=ROOT, lib_path=lib_path, artifacts=artifacts).run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
