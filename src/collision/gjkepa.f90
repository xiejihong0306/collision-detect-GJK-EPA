! GJK-EPA collision query
!
! GJK (Gilbert-Johnson-Keerthi):
!   S(d) = S_A(d) - S_B(-d)
!   iterate a simplex toward the origin in A-B
!   origin in simplex  <=>  A and B intersect
!
! EPA (Expanding Polytope):
!   start from the GJK tetrahedron
!   repeatedly expand the face closest to the origin
!   along its outward normal until S(n) adds no new depth:
!       n · w - d  <  eps     =>     depth = d,  normal = n
!
module GCLIB_GJKEPA
    use geom3d
    implicit none
    private

    integer, parameter :: GJK_MAX_ITER = 64
    integer, parameter :: EPA_MAX_ITER = 128
    integer, parameter :: EPA_MAX_VERTS = 256
    integer, parameter :: EPA_MAX_FACES = 512
    integer, parameter :: EPA_MAX_EDGES = 1024
    integer, parameter, public :: GJK_TRACE_MAX = 64
    integer, parameter, public :: EPA_TRACE_MAX = 128

    public :: GJKEPA, GJKEPA_query

    ! Diagnostic-only records. Filling them must not change collision results.
    type, public :: gjk_step_t
        real(rk) :: dir(3) = 0.0_rk
        real(rk) :: support(3) = 0.0_rk
        real(rk) :: closest(3) = 0.0_rk
        integer :: nsimp = 0
        integer :: reason = 0
    end type gjk_step_t

    type, public :: epa_step_t
        real(rk) :: normal(3) = 0.0_rk
        real(rk) :: dist = 0.0_rk
        real(rk) :: support(3) = 0.0_rk
        integer :: nfaces = 0
        integer :: nverts = 0
        integer :: reason = 0
    end type epa_step_t

    type :: epa_face
        integer :: v(3) = 0
        real(rk) :: n(3) = 0.0_rk
        real(rk) :: dist = 0.0_rk
        logical :: alive = .false.
        logical :: visible = .false.
    end type epa_face

contains

    subroutine GJKEPA(version_, TOL_FF_, &
                      p1_, p2_, collision_, colliType_, &
                      nearest_points_, collision_normal_, collision_point_, penetration_depth_)
        implicit none
        integer, intent(in) :: version_
        real(rk), intent(in) :: TOL_FF_
        real(rk), intent(in) :: p1_(:, :), p2_(:, :)
        logical, intent(out) :: collision_
        integer, intent(out) :: colliType_
        real(rk), intent(out) :: nearest_points_(2, 3)
        real(rk), intent(out) :: collision_normal_(3)
        real(rk), intent(out) :: collision_point_(3)
        real(rk), intent(out) :: penetration_depth_
        integer :: info, gjk_iters, epa_iters, epa_info

        call GJKEPA_query(version_, TOL_FF_, p1_, p2_, collision_, colliType_, &
                          nearest_points_, collision_normal_, collision_point_, &
                          penetration_depth_, info, gjk_iters, epa_iters, epa_info)
    end subroutine GJKEPA

    ! Same collision results as GJKEPA, plus status / iteration counts.
    ! info: 0 ok, 1 bad input, 2 GJK hit iter limit, 3 EPA failed to seed/expand
    subroutine GJKEPA_query(version_, TOL_FF_, p1_, p2_, collision_, colliType_, &
                            nearest_points_, collision_normal_, collision_point_, &
                            penetration_depth_, info_, gjk_iters_, epa_iters_, epa_info_, &
                            gjk_trace, n_gjk_trace, epa_trace, n_epa_trace)
        implicit none
        integer, intent(in) :: version_
        real(rk), intent(in) :: TOL_FF_
        real(rk), intent(in) :: p1_(:, :), p2_(:, :)
        logical, intent(out) :: collision_
        integer, intent(out) :: colliType_
        real(rk), intent(out) :: nearest_points_(2, 3)
        real(rk), intent(out) :: collision_normal_(3)
        real(rk), intent(out) :: collision_point_(3)
        real(rk), intent(out) :: penetration_depth_
        integer, intent(out) :: info_, gjk_iters_, epa_iters_, epa_info_
        type(gjk_step_t), optional, intent(inout) :: gjk_trace(:)
        integer, optional, intent(out) :: n_gjk_trace
        type(epa_step_t), optional, intent(inout) :: epa_trace(:)
        integer, optional, intent(out) :: n_epa_trace

        real(rk) :: simplex(4, 3), dir(3), epa_n(3), epa_d
        integer :: nsimp, info
        logical :: hit, gjk_capped

        collision_ = .false.
        colliType_ = 0
        nearest_points_ = 0.0_rk
        collision_normal_ = 0.0_rk
        collision_point_ = 0.0_rk
        penetration_depth_ = 0.0_rk
        info_ = 0
        gjk_iters_ = 0
        epa_iters_ = 0
        epa_info_ = 0
        if (present(n_gjk_trace)) n_gjk_trace = 0
        if (present(n_epa_trace)) n_epa_trace = 0

        if (size(p1_, 2) /= 3 .or. size(p2_, 2) /= 3) then
            info_ = 1
            return
        end if
        if (size(p1_, 1) < 1 .or. size(p2_, 1) < 1) then
            info_ = 1
            return
        end if

        if (.not. sphere_envelope_hit(p1_, p2_)) return

        call run_gjk(p1_, p2_, hit, simplex, nsimp, dir, gjk_iters_, gjk_capped, &
                     gjk_trace, n_gjk_trace)
        if (gjk_capped) info_ = 2
        if (.not. hit) then
            if (simplex_touches_origin(simplex, nsimp, collision_normal_, penetration_depth_)) then
                collision_ = .true.
                nearest_points_ = witness_points(p1_, p2_, collision_normal_)
                collision_point_ = collision_point_from_supports(version_, p1_, p2_, collision_normal_)
                colliType_ = classify_collision(p1_, p2_, collision_normal_, TOL_FF_)
            else
                nearest_points_(1, :) = support_vertex(p1_, dir)
                nearest_points_(2, :) = support_vertex(p2_, -dir)
            end if
            return
        end if

        collision_ = .true.
        penetration_depth_ = huge(1.0_rk)
        call refine_penetration(p1_, p2_, collision_normal_, penetration_depth_)
        if (nsimp >= 4) then
            call run_epa(p1_, p2_, simplex, nsimp, info, epa_n, epa_d, epa_iters_, &
                         epa_trace, n_epa_trace)
            epa_info_ = info
            if (info == 0 .and. epa_d > 1.0e-9_rk .and. epa_d < penetration_depth_) then
                penetration_depth_ = epa_d
                collision_normal_ = epa_n
            else if (info /= 0) then
                info_ = 3
            end if
        end if
        if (penetration_depth_ > 1.0e8_rk) penetration_depth_ = 0.0_rk
        if (norm2(collision_normal_) < GEOM_EPS) collision_normal_ = [0.0_rk, 0.0_rk, 1.0_rk]

        nearest_points_ = witness_points(p1_, p2_, collision_normal_)
        collision_point_ = collision_point_from_supports(version_, p1_, p2_, collision_normal_)
        colliType_ = classify_collision(p1_, p2_, collision_normal_, TOL_FF_)
    end subroutine GJKEPA_query

    logical function sphere_envelope_hit(p1, p2) result(hit)
        real(rk), intent(in) :: p1(:, :), p2(:, :)
        real(rk) :: c1(3), c2(3), r1, r2
        integer :: i
        real(rk), parameter :: pad = 1.0_rk
        c1 = polygon_centroid(p1)
        c2 = polygon_centroid(p2)
        r1 = 0.0_rk
        r2 = 0.0_rk
        do i = 1, size(p1, 1)
            r1 = max(r1, norm2(p1(i, :) - c1))
        end do
        do i = 1, size(p2, 1)
            r2 = max(r2, norm2(p2(i, :) - c2))
        end do
        hit = norm2(c1 - c2) <= r1 + r2 + pad
    end function sphere_envelope_hit

    !----------------------------------------------------------------
    ! GJK
    !----------------------------------------------------------------
    subroutine run_gjk(p1, p2, hit, simplex, nsimp, dir, niter, capped, trace, ntrace)
        real(rk), intent(in) :: p1(:, :), p2(:, :)
        logical, intent(out) :: hit
        real(rk), intent(out) :: simplex(4, 3)
        integer, intent(out) :: nsimp
        real(rk), intent(out) :: dir(3)
        integer, intent(out) :: niter
        logical, intent(out) :: capped
        type(gjk_step_t), optional, intent(inout) :: trace(:)
        integer, optional, intent(out) :: ntrace
        real(rk) :: w(3), c1(3), c2(3), v(3), vn2
        integer :: iter, reason

        hit = .false.
        capped = .false.
        niter = 0
        nsimp = 0
        simplex = 0.0_rk
        if (present(ntrace)) ntrace = 0
        c1 = polygon_centroid(p1)
        c2 = polygon_centroid(p2)
        dir = c1 - c2
        if (norm2(dir) < GEOM_EPS) dir = [1.0_rk, 0.0_rk, 0.0_rk]

        do iter = 1, GJK_MAX_ITER
            niter = iter
            w = minkowski_support(p1, p2, dir)
            reason = 1
            ! If the support plane does not pass the origin, bodies are separate.
            if (dot_product(w, dir) < -1.0e-9_rk) then
                hit = .false.
                reason = 2
                call record_gjk_step(trace, ntrace, dir, w, w, nsimp, reason)
                return
            end if
            if (nsimp >= 1) then
                if (already_vertex(simplex, nsimp, w)) then
                    call closest_point_on_simplex(simplex, nsimp, v)
                    hit = (norm2(v) <= 1.0e-6_rk)
                    reason = 3
                    call record_gjk_step(trace, ntrace, dir, w, v, nsimp, reason)
                    return
                end if
            end if
            nsimp = min(nsimp + 1, 4)
            simplex(nsimp, :) = w
            if (nsimp == 4 .and. origin_in_tetra(simplex)) then
                hit = .true.
                reason = 4
                call record_gjk_step(trace, ntrace, dir, w, [0.0_rk, 0.0_rk, 0.0_rk], nsimp, reason)
                return
            end if
            call reduce_simplex_toward_origin(simplex, nsimp, v)
            vn2 = dot_product(v, v)
            if (vn2 <= 1.0e-12_rk) then
                hit = .true.
                reason = 5
                call record_gjk_step(trace, ntrace, dir, w, v, nsimp, reason)
                return
            end if
            call record_gjk_step(trace, ntrace, dir, w, v, nsimp, reason)
            dir = -v
        end do
        call closest_point_on_simplex(simplex, nsimp, v)
        hit = (norm2(v) <= 1.0e-6_rk)
        capped = .true.
        call record_gjk_step(trace, ntrace, dir, w, v, nsimp, 6)
    end subroutine run_gjk

    subroutine record_gjk_step(trace, ntrace, dir, support, closest, nsimp, reason)
        type(gjk_step_t), optional, intent(inout) :: trace(:)
        integer, optional, intent(inout) :: ntrace
        real(rk), intent(in) :: dir(3), support(3), closest(3)
        integer, intent(in) :: nsimp, reason
        integer :: k
        if (.not. present(trace) .or. .not. present(ntrace)) return
        if (ntrace >= size(trace)) return
        ntrace = ntrace + 1
        k = ntrace
        trace(k)%dir = dir
        trace(k)%support = support
        trace(k)%closest = closest
        trace(k)%nsimp = nsimp
        trace(k)%reason = reason
    end subroutine record_gjk_step

    subroutine reduce_simplex_toward_origin(s, n, closest)
        real(rk), intent(inout) :: s(4, 3)
        integer, intent(inout) :: n
        real(rk), intent(out) :: closest(3)
        real(rk) :: p(3), best_p(3), best_d, d
        real(rk) :: cand(4, 3)
        integer :: best_n, i, j, k, m

        if (n <= 3) then
            call closest_point_on_simplex(s, n, closest)
            return
        end if

        if (origin_in_tetra(s)) then
            closest = 0.0_rk
            return
        end if

        best_d = huge(1.0_rk)
        best_n = 3
        best_p = s(1, :)
        ! drop one vertex at a time; keep the triangle closest to the origin
        do i = 1, 4
            m = 0
            do j = 1, 4
                if (j == i) cycle
                m = m + 1
                cand(m, :) = s(j, :)
            end do
            call closest_point_on_simplex(cand, 3, p)
            d = norm2(p)
            if (d < best_d) then
                best_d = d
                best_p = p
                best_n = 3
                do k = 1, 3
                    s(k, :) = cand(k, :)
                end do
            end if
        end do
        n = best_n
        closest = best_p
    end subroutine reduce_simplex_toward_origin

    subroutine refine_penetration(p1, p2, nml, depth)
        real(rk), intent(in) :: p1(:, :), p2(:, :)
        real(rk), intent(inout) :: nml(3), depth
        real(rk) :: dirs(3, 12), w(3), d, u(3), best, best_n(3)
        integer :: i
        dirs(:, 1) = [1.0_rk, 0.0_rk, 0.0_rk]
        dirs(:, 2) = [-1.0_rk, 0.0_rk, 0.0_rk]
        dirs(:, 3) = [0.0_rk, 1.0_rk, 0.0_rk]
        dirs(:, 4) = [0.0_rk, -1.0_rk, 0.0_rk]
        dirs(:, 5) = [0.0_rk, 0.0_rk, 1.0_rk]
        dirs(:, 6) = [0.0_rk, 0.0_rk, -1.0_rk]
        dirs(:, 7) = unit3([1.0_rk, 1.0_rk, 1.0_rk])
        dirs(:, 8) = unit3([1.0_rk, 1.0_rk, -1.0_rk])
        dirs(:, 9) = unit3([1.0_rk, -1.0_rk, 1.0_rk])
        dirs(:, 10) = unit3([-1.0_rk, 1.0_rk, 1.0_rk])
        dirs(:, 11) = unit3([1.0_rk, -1.0_rk, -1.0_rk])
        dirs(:, 12) = unit3([-1.0_rk, -1.0_rk, 1.0_rk])
        best = huge(1.0_rk)
        best_n = nml
        do i = 1, 12
            u = unit3(dirs(:, i))
            if (norm2(u) < GEOM_EPS) cycle
            w = minkowski_support(p1, p2, u)
            d = dot_product(u, w)
            if (d >= -1.0e-12_rk .and. d < best) then
                best = max(0.0_rk, d)
                best_n = u
            end if
        end do
        if (best < huge(1.0_rk)/2.0_rk) then
            depth = best
            nml = best_n
        end if
    end subroutine refine_penetration

    logical function simplex_touches_origin(s, n, nml, depth) result(yes)
        real(rk), intent(in) :: s(4, 3)
        integer, intent(in) :: n
        real(rk), intent(out) :: nml(3), depth
        real(rk) :: p(3), d
        nml = 0.0_rk
        depth = 0.0_rk
        yes = .false.
        if (n < 1) return
        call closest_point_on_simplex(s, n, p)
        d = norm2(p)
        if (d <= 1.0e-6_rk) then
            yes = .true.
            depth = 0.0_rk
            nml = unit3(p)
            if (norm2(nml) < GEOM_EPS .and. n >= 3) nml = plane_normal(s(1:3, :))
            if (norm2(nml) < GEOM_EPS) nml = [1.0_rk, 0.0_rk, 0.0_rk]
        end if
    end function simplex_touches_origin

    subroutine closest_point_on_simplex(s, n, p)
        real(rk), intent(in) :: s(4, 3)
        integer, intent(in) :: n
        real(rk), intent(out) :: p(3)
        real(rk) :: a(3), b(3), ab(3), t, uvw(3), tri(3, 3)
        select case (n)
        case (1)
            p = s(1, :)
        case (2)
            a = s(1, :); b = s(2, :)
            ab = b - a
            t = dot_product(-a, ab) / max(dot_product(ab, ab), GEOM_EPS)
            t = min(1.0_rk, max(0.0_rk, t))
            p = a + t * ab
        case default
            tri = s(1:3, :)
            call barycentric_triangle(tri, [0.0_rk, 0.0_rk, 0.0_rk], uvw)
            uvw = max(uvw, 0.0_rk)
            if (sum(uvw) < GEOM_EPS) then
                p = s(1, :)
            else
                uvw = uvw / sum(uvw)
                p = uvw(1)*s(1, :) + uvw(2)*s(2, :) + uvw(3)*s(3, :)
            end if
        end select
    end subroutine closest_point_on_simplex

    !----------------------------------------------------------------
    ! EPA — horizon expansion
    !----------------------------------------------------------------
    subroutine run_epa(p1, p2, simplex, nsimp, info, normal, depth, niter, trace, ntrace)
        real(rk), intent(in) :: p1(:, :), p2(:, :)
        real(rk), intent(in) :: simplex(4, 3)
        integer, intent(in) :: nsimp
        integer, intent(out) :: info
        real(rk), intent(out) :: normal(3), depth
        integer, intent(out) :: niter
        type(epa_step_t), optional, intent(inout) :: trace(:)
        integer, optional, intent(out) :: ntrace

        real(rk) :: verts(EPA_MAX_VERTS, 3)
        type(epa_face) :: faces(EPA_MAX_FACES)
        integer :: nv, nf, iter, imin
        real(rk) :: w(3)
        logical :: expanded

        info = 0
        niter = 0
        normal = 0.0_rk
        depth = 0.0_rk
        if (present(ntrace)) ntrace = 0

        call seed_polytope(p1, p2, simplex, nsimp, verts, nv, faces, nf, info)
        if (info /= 0 .or. nf < 4) then
            if (nsimp >= 3) then
                normal = unit3(plane_normal(simplex(1:3, :)))
                depth = abs(dot_product(normal, simplex(1, :)))
                info = 0
                call record_epa_step(trace, ntrace, normal, depth, simplex(1, :), nf, nv, 2)
            else
                info = 1
                call record_epa_step(trace, ntrace, normal, depth, [0.0_rk, 0.0_rk, 0.0_rk], nf, nv, 5)
            end if
            return
        end if

        do iter = 1, EPA_MAX_ITER
            niter = iter
            imin = closest_alive_face(faces, nf)
            if (imin < 1) then
                info = 2
                call record_epa_step(trace, ntrace, normal, depth, [0.0_rk, 0.0_rk, 0.0_rk], nf, nv, 5)
                return
            end if
            w = minkowski_support(p1, p2, faces(imin)%n)
            if (dot_product(w, faces(imin)%n) - faces(imin)%dist < 1.0e-6_rk) then
                normal = faces(imin)%n
                depth = max(0.0_rk, faces(imin)%dist)
                call record_epa_step(trace, ntrace, normal, depth, w, nf, nv, 2)
                return
            end if
            call expand_polytope(w, verts, nv, faces, nf, expanded)
            if (.not. expanded) then
                normal = faces(imin)%n
                depth = max(0.0_rk, faces(imin)%dist)
                call record_epa_step(trace, ntrace, normal, depth, w, nf, nv, 3)
                return
            end if
            call record_epa_step(trace, ntrace, faces(imin)%n, faces(imin)%dist, w, nf, nv, 1)
        end do

        imin = closest_alive_face(faces, nf)
        if (imin > 0) then
            normal = faces(imin)%n
            depth = max(0.0_rk, faces(imin)%dist)
            info = 0
            call record_epa_step(trace, ntrace, normal, depth, [0.0_rk, 0.0_rk, 0.0_rk], nf, nv, 4)
        else
            info = 3
            call record_epa_step(trace, ntrace, normal, depth, [0.0_rk, 0.0_rk, 0.0_rk], nf, nv, 5)
        end if
    end subroutine run_epa

    subroutine record_epa_step(trace, ntrace, nml, dist, support, nfaces, nverts, reason)
        type(epa_step_t), optional, intent(inout) :: trace(:)
        integer, optional, intent(inout) :: ntrace
        real(rk), intent(in) :: nml(3), dist, support(3)
        integer, intent(in) :: nfaces, nverts, reason
        integer :: k
        if (.not. present(trace) .or. .not. present(ntrace)) return
        if (ntrace >= size(trace)) return
        ntrace = ntrace + 1
        k = ntrace
        trace(k)%normal = nml
        trace(k)%dist = dist
        trace(k)%support = support
        trace(k)%nfaces = nfaces
        trace(k)%nverts = nverts
        trace(k)%reason = reason
    end subroutine record_epa_step

    subroutine seed_polytope(p1, p2, simplex, nsimp, verts, nv, faces, nf, info)
        real(rk), intent(in) :: p1(:, :), p2(:, :), simplex(4, 3)
        integer, intent(in) :: nsimp
        real(rk), intent(out) :: verts(EPA_MAX_VERTS, 3)
        type(epa_face), intent(out) :: faces(EPA_MAX_FACES)
        integer, intent(out) :: nv, nf, info
        real(rk) :: tet(4, 3), a(3), b(3), c(3), n(3)
        integer :: k

        info = 0
        verts = 0.0_rk
        nv = 0
        nf = 0
        faces = epa_face()

        tet = simplex
        if (nsimp < 4) then
            if (nsimp < 3) then
                info = 1
                return
            end if
            a = simplex(1, :); b = simplex(2, :); c = simplex(3, :)
            n = unit3(cross3(b - a, c - a))
            if (norm2(n) < GEOM_EPS) then
                info = 1
                return
            end if
            tet(1, :) = a; tet(2, :) = b; tet(3, :) = c
            tet(4, :) = minkowski_support(p1, p2, n)
            if (.not. origin_in_tetra(tet)) then
                tet(4, :) = minkowski_support(p1, p2, -n)
            end if
            if (abs(signed_dist_plane(tet(4, :), tet(1:3, :))) < GEOM_EPS8) then
                info = 1
                return
            end if
        end if

        if (.not. origin_in_tetra(tet)) then
            n = unit3(cross3(tet(2, :) - tet(1, :), tet(3, :) - tet(1, :)))
            tet(4, :) = minkowski_support(p1, p2, -n)
        end if

        do k = 1, 4
            nv = nv + 1
            verts(nv, :) = tet(k, :)
        end do
        ! Faces: ABC, ACD, ADB, BCD  (winding later fixed by origin)
        call add_face(faces, nf, verts, 1, 2, 3)
        call add_face(faces, nf, verts, 1, 3, 4)
        call add_face(faces, nf, verts, 1, 4, 2)
        call add_face(faces, nf, verts, 2, 4, 3)
        if (nf < 4) info = 1
    end subroutine seed_polytope

    subroutine add_face(faces, nf, verts, i1, i2, i3)
        type(epa_face), intent(inout) :: faces(:)
        integer, intent(inout) :: nf
        real(rk), intent(in) :: verts(:, :)
        integer, intent(in) :: i1, i2, i3
        real(rk) :: e1(3), e2(3), n(3), d
        integer :: ia, ib, ic
        if (nf >= EPA_MAX_FACES) return
        ia = i1; ib = i2; ic = i3
        e1 = verts(ib, :) - verts(ia, :)
        e2 = verts(ic, :) - verts(ia, :)
        n = unit3(cross3(e1, e2))
        if (norm2(n) < GEOM_EPS) return
        d = dot_product(n, verts(ia, :))
        if (d < 0.0_rk) then
            n = -n
            d = -d
            ib = i3
            ic = i2
        end if
        nf = nf + 1
        faces(nf)%v = [ia, ib, ic]
        faces(nf)%n = n
        faces(nf)%dist = d
        faces(nf)%alive = .true.
        faces(nf)%visible = .false.
    end subroutine add_face

    integer function closest_alive_face(faces, nf) result(imin)
        type(epa_face), intent(in) :: faces(:)
        integer, intent(in) :: nf
        integer :: i
        real(rk) :: best
        imin = -1
        best = huge(1.0_rk)
        do i = 1, nf
            if (faces(i)%alive .and. faces(i)%dist < best) then
                best = faces(i)%dist
                imin = i
            end if
        end do
    end function closest_alive_face

    logical function already_vertex(verts, nv, w) result(yes)
        real(rk), intent(in) :: verts(:, :), w(3)
        integer, intent(in) :: nv
        yes = vertex_index(verts, nv, w) > 0
    end function already_vertex

    integer function vertex_index(verts, nv, w) result(idx)
        real(rk), intent(in) :: verts(:, :), w(3)
        integer, intent(in) :: nv
        integer :: i
        idx = -1
        do i = 1, nv
            if (norm2(verts(i, :) - w) < GEOM_EPS8) then
                idx = i
                return
            end if
        end do
    end function vertex_index

    subroutine expand_polytope(w, verts, nv, faces, nf, expanded)
        real(rk), intent(in) :: w(3)
        real(rk), intent(inout) :: verts(:, :)
        integer, intent(inout) :: nv
        type(epa_face), intent(inout) :: faces(:)
        integer, intent(inout) :: nf
        logical, intent(out) :: expanded
        integer :: i, e1, e2, nvis, nhor
        integer :: horizon(EPA_MAX_EDGES, 2)
        integer :: iw

        expanded = .false.
        nvis = 0
        do i = 1, nf
            faces(i)%visible = .false.
            if (.not. faces(i)%alive) cycle
            if (dot_product(faces(i)%n, w) > faces(i)%dist + GEOM_EPS8) then
                faces(i)%visible = .true.
                nvis = nvis + 1
            end if
        end do
        if (nvis == 0) return

        call collect_horizon(faces, nf, horizon, nhor)
        if (nhor < 3) return

        iw = vertex_index(verts, nv, w)
        if (iw < 1) then
            if (nv >= EPA_MAX_VERTS) return
            nv = nv + 1
            iw = nv
            verts(iw, :) = w
        end if

        do i = 1, nf
            if (faces(i)%visible) faces(i)%alive = .false.
        end do
        do i = 1, nhor
            e1 = horizon(i, 1)
            e2 = horizon(i, 2)
            call add_face(faces, nf, verts, iw, e1, e2)
        end do
        expanded = .true.
    end subroutine expand_polytope

    subroutine collect_horizon(faces, nf, horizon, nhor)
        type(epa_face), intent(in) :: faces(:)
        integer, intent(in) :: nf
        integer, intent(out) :: horizon(:, :)
        integer, intent(out) :: nhor
        integer :: i, e, a, b
        integer :: ev(3)
        nhor = 0
        do i = 1, nf
            if (.not. faces(i)%alive .or. .not. faces(i)%visible) cycle
            ev = faces(i)%v
            do e = 1, 3
                a = ev(e)
                b = ev(modulo(e, 3) + 1)
                if (.not. edge_shared_with_visible(faces, nf, i, a, b)) then
                    if (nhor >= EPA_MAX_EDGES) return
                    nhor = nhor + 1
                    horizon(nhor, 1) = a
                    horizon(nhor, 2) = b
                end if
            end do
        end do
    end subroutine collect_horizon

    logical function edge_shared_with_visible(faces, nf, self, a, b) result(yes)
        type(epa_face), intent(in) :: faces(:)
        integer, intent(in) :: nf, self, a, b
        integer :: j
        yes = .false.
        do j = 1, nf
            if (j == self) cycle
            if (.not. faces(j)%alive .or. .not. faces(j)%visible) cycle
            if (face_has_edge(faces(j)%v, a, b)) then
                yes = .true.
                return
            end if
        end do
    end function edge_shared_with_visible

    pure logical function face_has_edge(v, a, b) result(yes)
        integer, intent(in) :: v(3), a, b
        yes = (count_eq(v, a) == 1) .and. (count_eq(v, b) == 1)
    end function face_has_edge

    pure integer function count_eq(v, x) result(c)
        integer, intent(in) :: v(3), x
        integer :: i
        c = 0
        do i = 1, 3
            if (v(i) == x) c = c + 1
        end do
    end function count_eq

    !----------------------------------------------------------------
    ! Witness / contact
    !----------------------------------------------------------------
    function witness_points(p1, p2, nml) result(pts)
        real(rk), intent(in) :: p1(:, :), p2(:, :), nml(3)
        real(rk) :: pts(2, 3)
        pts(1, :) = support_vertex(p1, nml)
        pts(2, :) = support_vertex(p2, -nml)
    end function witness_points

    integer function classify_collision(p1, p2, nml, tol) result(kind)
        real(rk), intent(in) :: p1(:, :), p2(:, :), nml(3), tol
        real(rk), allocatable :: s1(:, :), s2(:, :)
        kind = 1
        call collect_supports(p1, nml, max(tol, 1.0e-6_rk), s1)
        call collect_supports(p2, -nml, max(tol, 1.0e-6_rk), s2)
        if (size(s1, 1) >= 3 .and. size(s2, 1) >= 3) kind = 2
        if (allocated(s1)) deallocate(s1)
        if (allocated(s2)) deallocate(s2)
    end function classify_collision

    function collision_point_from_supports(version, p1, p2, nml) result(pt)
        integer, intent(in) :: version
        real(rk), intent(in) :: p1(:, :), p2(:, :), nml(3)
        real(rk) :: pt(3)
        real(rk), allocatable :: s1(:, :), s2(:, :)
        integer :: n1, n2
        real(rk) :: feet(2, 3), line1(2, 3), line2(2, 3)
        real(rk), parameter :: atol = 1.0e-1_rk

        pt = 0.0_rk
        call collect_supports(p1, nml, atol, s1)
        call collect_supports(p2, -nml, atol, s2)
        n1 = size(s1, 1)
        n2 = size(s2, 1)

        if (version == 1) then
            if (n1 == 1 .and. n2 == 1) then
                pt = 0.5_rk * (s1(1, :) + s2(1, :))
            else if (n1 == 1) then
                pt = s1(1, :)
            else if (n2 == 1) then
                pt = s2(1, :)
            else
                pt = polygon_centroid(s1)
            end if
        else
            ! version 2 (default): feature classification
            if (n1 == 1 .and. n2 == 1) then
                pt = 0.5_rk * (s1(1, :) + s2(1, :))
            else if (n1 == 1) then
                pt = s1(1, :)
            else if (n2 == 1) then
                pt = s2(1, :)
            else if (n1 == 2 .and. n2 == 2) then
                line1 = s1(1:2, :)
                line2 = s2(1:2, :)
                feet = feet_segment_segment(line1, line2)
                pt = 0.5_rk * (feet(1, :) + feet(2, :))
            else if (n1 >= 3 .and. n2 == 2) then
                pt = 0.5_rk * (s2(1, :) + s2(2, :))
            else if (n2 >= 3 .and. n1 == 2) then
                pt = 0.5_rk * (s1(1, :) + s1(2, :))
            else
                pt = polygon_centroid(s1)
            end if
        end if
        if (allocated(s1)) deallocate(s1)
        if (allocated(s2)) deallocate(s2)
    end function collision_point_from_supports

end module GCLIB_GJKEPA
