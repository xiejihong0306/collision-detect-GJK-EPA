This directory keeps the 2023 hand-written GJK-EPA module (`GCLIB_GJKEPA.f90`).
It is not compiled. The active implementation lives in `src/collision/gjkepa.f90`
and uses horizon EPA instead of rebuilding the polytope with QuickHull each step.
