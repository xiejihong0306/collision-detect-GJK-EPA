# GJK-EPA

Fortran library for 3D convex-body collision detection
([Gilbert–Johnson–Keerthi](https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm)
+ Expanding Polytope) and 3D convex hull meshing (QuickHull).

三维凸体碰撞检测（GJK + EPA）与凸包三角剖分（QuickHull）的 Fortran 实现。

## Algorithms

### GJK — does the Minkowski difference contain the origin?

```
S(d) = S_A(d) - S_B(-d)

iterate a simplex toward 0 in A - B
0 in simplex  <=>  A intersects B
```

### EPA — how deep, and in which direction?

Start from the GJK tetrahedron. Repeatedly expand the face closest to the
origin along its outward normal `n` until the support mapping adds no depth:

```
w = S(n)
if  n · w - d  <  eps
    penetration depth = d
    contact normal    = n
else
    delete faces visible from w
    stitch new triangles along the horizon
```

QuickHull is a first-class preprocessing tool: raw points → convex triangle
mesh → feed the vertices to `GJKEPA`.

## Build

Requires CMake ≥ 3.20 and a Fortran 2008 compiler (`gfortran` or Intel `ifx`/`ifort`).
OpenMP is used if the compiler provides it.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The example program prints a cube-vs-rotated-cube query:

```bash
./build/bin/example_two_polyhedra
```

On Windows the binaries are under `build\bin\`.

## API

```fortran
use GCLIB_GJKEPA, only: GJKEPA
use GCLIB_QuickHull, only: QuickHull

call GJKEPA(version, tol_ff, p1, p2, collision, colli_type, &
            nearest, normal, point, depth)
! p1, p2     : (n,3) convex-hull vertices
! colli_type : 0 separated / 1 penetrating / 2 face contact
! version    : 1 coarse contact point, 2 feature classification

call QuickHull(points, hull_mesh, info)
! hull_mesh  : (nfaces,3,3) triangle mesh
! info       : 0 success, nonzero failure (never pauses / stops)
```

## Layout

```
src/geom/geom3d.f90         shared 3D primitives
src/collision/gjkepa.f90    GJK + horizon EPA
src/hull/quickhull.f90      QuickHull (from the former quick-hull-mesh repo)
src/legacy/                 original 2023 GJK-EPA source, not compiled
examples/two_polyhedra.f90
tests/
```

## License

MIT. Copyright (c) 2023-2026 XIE Jihong / Jairox.
