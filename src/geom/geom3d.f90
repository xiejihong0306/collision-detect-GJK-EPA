! geom3d — shared 3D primitives for GJK / EPA / QuickHull
!
!   a x b = | i     j     k    |
!           | a1    a2    a3   |
!           | b1    b2    b3   |
!
!   signed plane distance:  d = n · (x - p0),   n = (p1-p0) x (p2-p1) / |...|
!
module geom3d
    implicit none
    private

    integer, parameter, public :: rk = selected_real_kind(15, 307)
    integer, parameter, public :: ik = selected_int_kind(9)

    real(rk), parameter, public :: GEOM_EPS = 1.0e-12_rk
    real(rk), parameter, public :: GEOM_EPS8 = 1.0e-8_rk

    public :: cross3, unit3, plane_normal, signed_dist_plane
    public :: support_vertex, minkowski_support, collect_supports
    public :: origin_in_tetra, point_in_triangle
    public :: unique_mesh_vertices
    public :: foot_point_line, feet_segment_segment
    public :: polygon_centroid, sort_coplanar_ccw
    public :: barycentric_triangle, same_dir
    public :: rotate_axis_angle

contains

    pure function cross3(a, b) result(c)
        real(rk), intent(in) :: a(3), b(3)
        real(rk) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), &
             a(3)*b(1) - a(1)*b(3), &
             a(1)*b(2) - a(2)*b(1)]
    end function cross3

    pure function unit3(v) result(u)
        real(rk), intent(in) :: v(3)
        real(rk) :: u(3), n
        n = norm2(v)
        if (n < GEOM_EPS) then
            u = 0.0_rk
        else
            u = v / n
        end if
    end function unit3

    pure function same_dir(a, b) result(yes)
        real(rk), intent(in) :: a(3), b(3)
        logical :: yes
        yes = dot_product(a, b) > 0.0_rk
    end function same_dir

    ! tri(3,3): three vertices, each row xyz
    pure function plane_normal(tri) result(n)
        real(rk), intent(in) :: tri(3, 3)
        real(rk) :: n(3)
        n = unit3(cross3(tri(2, :) - tri(1, :), tri(3, :) - tri(2, :)))
    end function plane_normal

    pure function signed_dist_plane(p, tri) result(d)
        real(rk), intent(in) :: p(3), tri(3, 3)
        real(rk) :: d, n(3)
        n = plane_normal(tri)
        d = dot_product(p - tri(1, :), n)
    end function signed_dist_plane

    function support_vertex(pts, dir) result(s)
        real(rk), intent(in) :: pts(:, :), dir(3)
        real(rk) :: s(3)
        integer :: i, imax
        real(rk) :: best, d
        imax = 1
        best = -huge(1.0_rk)
        do i = 1, size(pts, 1)
            d = dot_product(dir, pts(i, :))
            if (d > best) then
                best = d
                imax = i
            end if
        end do
        s = pts(imax, :)
    end function support_vertex

    ! Minkowski difference support:
    !   S_{A-B}(d) = S_A(d) - S_B(-d)
    function minkowski_support(p1, p2, dir) result(s)
        real(rk), intent(in) :: p1(:, :), p2(:, :), dir(3)
        real(rk) :: s(3)
        s = support_vertex(p1, dir) - support_vertex(p2, -dir)
    end function minkowski_support

    subroutine collect_supports(pts, dir, tol, out)
        real(rk), intent(in) :: pts(:, :), dir(3), tol
        real(rk), allocatable, intent(out) :: out(:, :)
        integer :: i, n, c
        real(rk) :: best, d
        real(rk), allocatable :: tmp(:, :)
        n = size(pts, 1)
        best = -huge(1.0_rk)
        do i = 1, n
            d = dot_product(dir, pts(i, :))
            if (d > best) best = d
        end do
        allocate(tmp(n, 3))
        c = 0
        do i = 1, n
            d = dot_product(dir, pts(i, :))
            if (d > best - tol) then
                c = c + 1
                tmp(c, :) = pts(i, :)
            end if
        end do
        allocate(out(c, 3))
        if (c > 0) out = tmp(1:c, :)
        deallocate(tmp)
    end subroutine collect_supports

    ! Origin in tetrahedron ABCD via consistent signed volumes:
    !   V(A,B,C,D) = (B-A) · ((C-A) x (D-A))
    ! Origin inside iff the four tets (O,BCD), (A,O,CD), (A,B,O,D), (A,B,C,O)
    ! have the same sign as V(A,B,C,D) (including faces: |V| < eps).
    pure function origin_in_tetra(tet) result(inside)
        real(rk), intent(in) :: tet(4, 3)
        logical :: inside
        real(rk) :: a(3), b(3), c(3), d(3), o(3)
        real(rk) :: v0, v1, v2, v3, v4
        a = tet(1, :); b = tet(2, :); c = tet(3, :); d = tet(4, :)
        o = 0.0_rk
        v0 = scalar_triple(b - a, c - a, d - a)
        if (abs(v0) < GEOM_EPS) then
            inside = .false.
            return
        end if
        v1 = scalar_triple(b - o, c - o, d - o)
        v2 = scalar_triple(o - a, c - a, d - a)
        v3 = scalar_triple(b - a, o - a, d - a)
        v4 = scalar_triple(b - a, c - a, o - a)
        inside = (v1*v0 >= -GEOM_EPS8) .and. (v2*v0 >= -GEOM_EPS8) .and. &
                 (v3*v0 >= -GEOM_EPS8) .and. (v4*v0 >= -GEOM_EPS8)
    end function origin_in_tetra

    pure function scalar_triple(u, v, w) result(s)
        real(rk), intent(in) :: u(3), v(3), w(3)
        real(rk) :: s
        s = dot_product(u, cross3(v, w))
    end function scalar_triple

    ! Point-in-triangle in 3D, including boundary. Projects to the most
    ! stable 2D plane (largest |n_i|).
    pure function point_in_triangle(tri, p) result(inside)
        real(rk), intent(in) :: tri(3, 3), p(3)
        logical :: inside
        real(rk) :: n(3), u(2), v(2), w(2), uu, uv, vv, wu, wv, den, s, t
        integer :: k
        n = plane_normal(tri)
        k = maxloc(abs(n), dim=1)
        u = drop_axis(tri(2, :) - tri(1, :), k)
        v = drop_axis(tri(3, :) - tri(1, :), k)
        w = drop_axis(p - tri(1, :), k)
        uu = dot_product(u, u)
        uv = dot_product(u, v)
        vv = dot_product(v, v)
        wu = dot_product(w, u)
        wv = dot_product(w, v)
        den = uu*vv - uv*uv
        if (abs(den) < GEOM_EPS) then
            inside = .false.
            return
        end if
        s = (vv*wu - uv*wv) / den
        t = (uu*wv - uv*wu) / den
        inside = (s >= -GEOM_EPS8) .and. (t >= -GEOM_EPS8) .and. (s + t <= 1.0_rk + GEOM_EPS8)
    end function point_in_triangle

    pure function drop_axis(v, k) result(q)
        real(rk), intent(in) :: v(3)
        integer, intent(in) :: k
        real(rk) :: q(2)
        select case (k)
        case (1)
            q = v(2:3)
        case (2)
            q = [v(1), v(3)]
        case default
            q = v(1:2)
        end select
    end function drop_axis

    ! Unique vertices of a triangular mesh mesh(nfaces, 3, 3)
    subroutine unique_mesh_vertices(mesh, verts, info)
        real(rk), intent(in) :: mesh(:, :, :)
        real(rk), allocatable, intent(out) :: verts(:, :)
        integer, intent(out) :: info
        integer :: i, j, k, nf, nv, cap
        real(rk), allocatable :: buf(:, :)
        real(rk) :: p(3)
        logical :: found
        info = 0
        if (size(mesh, 2) /= 3 .or. size(mesh, 3) /= 3) then
            info = 1
            return
        end if
        nf = size(mesh, 1)
        if (nf < 1) then
            info = 1
            return
        end if
        cap = max(8, nf)
        allocate(buf(cap, 3))
        nv = 0
        do i = 1, nf
            do j = 1, 3
                p = mesh(i, j, :)
                found = .false.
                do k = 1, nv
                    if (norm2(p - buf(k, :)) < GEOM_EPS) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) then
                    if (nv == cap) then
                        call grow_buf(buf, cap)
                    end if
                    nv = nv + 1
                    buf(nv, :) = p
                end if
            end do
        end do
        allocate(verts(nv, 3))
        verts = buf(1:nv, :)
        deallocate(buf)
    end subroutine unique_mesh_vertices

    subroutine grow_buf(buf, cap)
        real(rk), allocatable, intent(inout) :: buf(:, :)
        integer, intent(inout) :: cap
        real(rk), allocatable :: tmp(:, :)
        integer :: n
        n = size(buf, 1)
        allocate(tmp(n, 3))
        tmp = buf
        deallocate(buf)
        cap = cap * 2
        allocate(buf(cap, 3))
        buf(1:n, :) = tmp
        deallocate(tmp)
    end subroutine grow_buf

    pure function foot_point_line(p, line) result(f)
        real(rk), intent(in) :: p(3), line(2, 3)
        real(rk) :: f(3), d(3), t, dn
        d = line(2, :) - line(1, :)
        dn = dot_product(d, d)
        if (dn < GEOM_EPS) then
            f = line(1, :)
            return
        end if
        t = dot_product(p - line(1, :), d) / dn
        f = line(1, :) + t * d
    end function foot_point_line

    ! Closest points on two (infinite) lines. Parallel: mid of line1 + foot.
    function feet_segment_segment(l1, l2) result(feet)
        real(rk), intent(in) :: l1(2, 3), l2(2, 3)
        real(rk) :: feet(2, 3)
        real(rk) :: p1(3), q1(3), p2(3), q2(3)
        real(rk) :: d1(3), d2(3), r(3)
        real(rk) :: a, b, c, e, f, den, s, t
        p1 = l1(1, :); q1 = l1(2, :)
        p2 = l2(1, :); q2 = l2(2, :)
        d1 = q1 - p1
        d2 = q2 - p2
        r = p1 - p2
        a = dot_product(d1, d1)
        b = dot_product(d1, d2)
        c = dot_product(d1, r)
        e = dot_product(d2, d2)
        f = dot_product(d2, r)
        den = a*e - b*b
        if (abs(den) < GEOM_EPS) then
            feet(1, :) = 0.5_rk * (p1 + q1)
            feet(2, :) = foot_point_line(feet(1, :), l2)
        else
            s = (b*f - c*e) / den
            t = (a*f - b*c) / den
            feet(1, :) = p1 + s * d1
            feet(2, :) = p2 + t * d2
        end if
    end function feet_segment_segment

    pure function polygon_centroid(pts) result(c)
        real(rk), intent(in) :: pts(:, :)
        real(rk) :: c(3)
        integer :: i, n
        n = size(pts, 1)
        c = 0.0_rk
        if (n < 1) return
        do i = 1, 3
            c(i) = sum(pts(:, i)) / real(n, rk)
        end do
    end function polygon_centroid

    function sort_coplanar_ccw(pts) result(ordered)
        real(rk), intent(in) :: pts(:, :)
        real(rk) :: ordered(size(pts, 1), 3)
        real(rk) :: c(3), n(3), v1(3), v2(3), ang, best
        integer :: i, j, npts, idx
        logical :: used(size(pts, 1))
        npts = size(pts, 1)
        ordered = 0.0_rk
        if (npts < 1) return
        if (npts == 1) then
            ordered(1, :) = pts(1, :)
            return
        end if
        c = polygon_centroid(pts)
        if (npts >= 3) then
            n = unit3(cross3(pts(2, :) - pts(1, :), pts(3, :) - pts(1, :)))
        else
            n = [0.0_rk, 0.0_rk, 1.0_rk]
        end if
        if (norm2(n) < GEOM_EPS) n = [0.0_rk, 0.0_rk, 1.0_rk]
        used = .false.
        ordered(1, :) = pts(1, :)
        used(1) = .true.
        do i = 2, npts
            best = huge(1.0_rk)
            idx = -1
            do j = 1, npts
                if (used(j)) cycle
                v1 = pts(j, :) - c
                v2 = ordered(i - 1, :) - c
                ang = atan2(dot_product(n, cross3(v2, v1)), dot_product(v1, v2))
                ang = modulo(ang + 2.0_rk*acos(-1.0_rk), 2.0_rk*acos(-1.0_rk))
                if (ang < best) then
                    best = ang
                    idx = j
                end if
            end do
            if (idx < 0) exit
            ordered(i, :) = pts(idx, :)
            used(idx) = .true.
        end do
    end function sort_coplanar_ccw

    ! Barycentric coords of p on triangle ABC.  p = u A + v B + w C, u+v+w=1
    pure subroutine barycentric_triangle(tri, p, uvw)
        real(rk), intent(in) :: tri(3, 3), p(3)
        real(rk), intent(out) :: uvw(3)
        real(rk) :: v0(3), v1(3), v2(3), d00, d01, d11, d20, d21, den
        v0 = tri(2, :) - tri(1, :)
        v1 = tri(3, :) - tri(1, :)
        v2 = p - tri(1, :)
        d00 = dot_product(v0, v0)
        d01 = dot_product(v0, v1)
        d11 = dot_product(v1, v1)
        d20 = dot_product(v2, v0)
        d21 = dot_product(v2, v1)
        den = d00*d11 - d01*d01
        if (abs(den) < GEOM_EPS) then
            uvw = [1.0_rk, 0.0_rk, 0.0_rk]
            return
        end if
        uvw(2) = (d11*d20 - d01*d21) / den
        uvw(3) = (d00*d21 - d01*d20) / den
        uvw(1) = 1.0_rk - uvw(2) - uvw(3)
    end subroutine barycentric_triangle

    ! Rodrigues: R = I cosθ + (k k^T)(1-cosθ) + K sinθ
    pure function rotate_axis_angle(p, axis, angle_deg) result(q)
        real(rk), intent(in) :: p(3), axis(3), angle_deg
        real(rk) :: q(3), k(3), ang, c, s
        k = unit3(axis)
        if (norm2(k) < GEOM_EPS) then
            q = p
            return
        end if
        ang = angle_deg * acos(-1.0_rk) / 180.0_rk
        c = cos(ang)
        s = sin(ang)
        q = p*c + cross3(k, p)*s + k*dot_product(k, p)*(1.0_rk - c)
    end function rotate_axis_angle

end module geom3d
