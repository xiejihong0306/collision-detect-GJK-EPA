# GJK-EPA verification bench

Interactive and headless checks for **this repository's** Fortran GJK-EPA.
Results come from `GJKEPA_query` through a C ABI. SAT is an independent
reference, not a stand-in for the algorithm.

This is a discrete-query bench. Fast dragging can skip intermediate contacts;
there is no continuous collision detection.

## Start

From the repository root, with no pre-build step:

```text
python launcher.py
python launcher.py --case face_touch
python launcher.py --headless --suite regression
python launcher.py --headless --suite random --seed 42 --count 1000
python launcher.py --replay artifacts/failure_xxx.json
```

`launcher.py` is the only user entry. It uses `.venv`, installs
`bench/requirements.txt` there if needed, configures CMake, and incrementally
builds the `gjkepa_c` shared library. A failed build aborts and does not load
an older binary.

## Keys (interactive)

| Key | Action |
| --- | --- |
| C | camera (left-drag rotates the view) |
| T | translate B |
| R | rotate B |
| X Y Z | constrain translation / rotation axis (toggle) |
| [ ] | previous / next case |
| 0 | reset B pose |
| P | pause automatic queries |
| N | single-step query |
| D | GJK/EPA traces + Minkowski pane content |
| S | save snapshot under `artifacts/` |
| 滚轮 | 改 B 的位姿（不缩放窗口） |

The mouse wheel never zooms the camera. In translate or camera mode it
nudges B along the view axis (or X/Y/Z if constrained). In rotate mode it
spins B about that axis. Middle-mouse pans, right-mouse still zooms.
Left-drag moves B in translate/rotate modes on the current view plane.

Colors: green separated, yellow contact / tolerance band, red penetrating,
purple failed or undecidable. Failures are never painted as a clean miss.

## Adapter conventions

```
S(d) = S_A(d) - S_B(-d)
```

EPA `n` is the outward normal of the closest face of `A-B`. Witnesses are
`support(A, n)` and `support(B, -n)`. Depenetration moves **B** by
`sign(n · (A-B)) * n * depth` (witness sign, not centroid sign).

`TOL_FF` is the original feature-classification band (default `1.0` to match
existing Fortran tests), not the GJK/EPA geometric epsilon.

## Assumptions

- Vertices are `(n,3)` in the same length units as the original library.
- Random clouds are passed through this repo's QuickHull before use.
- QuickHull still uses `SAVE` internals; hull calls are serialized.
- Iteration traces are optional diagnostics and default off.
- SAT MTV uses one-sided exit distances so containment is not confused with
  projected overlap length.
- Symmetric or coincident poses may have non-unique normals; tests allow that.

## Known algorithm issues (not hidden)

`invariance_overlap` currently **fails** and is supposed to: after the same
rigid motion of both cubes, SAT depth stays 0.5 while GJK–EPA depth jumps
from 0.5 to about 0.945. Snapshot: `artifacts/failure_invariance_overlap.json`.
Cause: when the GJK simplex is not a tetrahedron, depth comes from
`refine_penetration` sampling world-axis / diagonal directions, which is not
rotation-invariant. Do not “fix” this by loosening the test.

## Layout

```
launcher.py
bench/bootstrap.py     venv + CMake
bench/adapter.py       ctypes + conventions
bench/geom.py          poses and primitives
bench/sat.py           independent SAT
bench/cases.py         named scenes
bench/verify.py        suites
bench/snapshot.py      save / replay
bench/viz.py           PyVista
src/collision/gjkepa_capi.f90
```
