! C ABI for the verification bench. Calls existing GJKEPA_query / QuickHull.
! Does not reimplement collision geometry.
module gjkepa_capi
    use, intrinsic :: iso_c_binding
    use GCLIB_GJKEPA, only: GJKEPA_query, gjk_step_t, epa_step_t, GJK_TRACE_MAX, EPA_TRACE_MAX
    use GCLIB_QuickHull, only: QuickHull
    implicit none
contains

    subroutine gjkepa_query_c(p1, n1, p2, n2, version, tol_ff, want_trace, &
                              collision, colli_type, info, gjk_iters, epa_iters, epa_info, &
                              nearest_a, nearest_b, normal, point, depth, elapsed_s, &
                              gjk_dir, gjk_support, gjk_closest, gjk_nsimp, gjk_reason, n_gjk, &
                              epa_nml, epa_dist, epa_support, epa_nfaces, epa_nverts, epa_reason, n_epa) &
            bind(c, name="gjkepa_query_c")
        real(c_double), intent(in) :: p1(3, n1), p2(3, n2)
        integer(c_int), intent(in), value :: n1, n2, version, want_trace
        real(c_double), intent(in), value :: tol_ff
        integer(c_int), intent(out) :: collision, colli_type, info, gjk_iters, epa_iters, epa_info
        real(c_double), intent(out) :: nearest_a(3), nearest_b(3), normal(3), point(3)
        real(c_double), intent(out) :: depth, elapsed_s
        real(c_double), intent(out) :: gjk_dir(3, GJK_TRACE_MAX)
        real(c_double), intent(out) :: gjk_support(3, GJK_TRACE_MAX)
        real(c_double), intent(out) :: gjk_closest(3, GJK_TRACE_MAX)
        integer(c_int), intent(out) :: gjk_nsimp(GJK_TRACE_MAX), gjk_reason(GJK_TRACE_MAX), n_gjk
        real(c_double), intent(out) :: epa_nml(3, EPA_TRACE_MAX)
        real(c_double), intent(out) :: epa_dist(EPA_TRACE_MAX)
        real(c_double), intent(out) :: epa_support(3, EPA_TRACE_MAX)
        integer(c_int), intent(out) :: epa_nfaces(EPA_TRACE_MAX), epa_nverts(EPA_TRACE_MAX)
        integer(c_int), intent(out) :: epa_reason(EPA_TRACE_MAX), n_epa

        real(8), allocatable :: a(:, :), b(:, :)
        real(8) :: nearest(2, 3), nml(3), pt(3), dep
        logical :: hit
        integer :: kind, inf, igjk, iepa, ieinfo, i, ng, ne
        integer(8) :: t0, t1, rate
        type(gjk_step_t) :: gtr(GJK_TRACE_MAX)
        type(epa_step_t) :: etr(EPA_TRACE_MAX)

        collision = 0
        colli_type = 0
        info = 1
        gjk_iters = 0
        epa_iters = 0
        epa_info = 0
        nearest_a = 0.0_c_double
        nearest_b = 0.0_c_double
        normal = 0.0_c_double
        point = 0.0_c_double
        depth = 0.0_c_double
        elapsed_s = 0.0_c_double
        n_gjk = 0
        n_epa = 0
        gjk_dir = 0.0_c_double
        gjk_support = 0.0_c_double
        gjk_closest = 0.0_c_double
        gjk_nsimp = 0
        gjk_reason = 0
        epa_nml = 0.0_c_double
        epa_dist = 0.0_c_double
        epa_support = 0.0_c_double
        epa_nfaces = 0
        epa_nverts = 0
        epa_reason = 0

        if (n1 < 1 .or. n2 < 1) return

        allocate(a(n1, 3), b(n2, 3))
        a = transpose(p1)
        b = transpose(p2)

        call system_clock(t0, rate)
        if (want_trace /= 0) then
            call GJKEPA_query(int(version), real(tol_ff, 8), a, b, hit, kind, &
                              nearest, nml, pt, dep, inf, igjk, iepa, ieinfo, &
                              gtr, ng, etr, ne)
        else
            ng = 0
            ne = 0
            call GJKEPA_query(int(version), real(tol_ff, 8), a, b, hit, kind, &
                              nearest, nml, pt, dep, inf, igjk, iepa, ieinfo)
        end if
        call system_clock(t1, rate)

        if (hit) collision = 1
        colli_type = kind
        info = inf
        gjk_iters = igjk
        epa_iters = iepa
        epa_info = ieinfo
        nearest_a = nearest(1, :)
        nearest_b = nearest(2, :)
        normal = nml
        point = pt
        depth = dep
        if (rate > 0) elapsed_s = real(t1 - t0, c_double) / real(rate, c_double)

        if (want_trace /= 0) then
            n_gjk = min(ng, GJK_TRACE_MAX)
            do i = 1, n_gjk
                gjk_dir(:, i) = gtr(i)%dir
                gjk_support(:, i) = gtr(i)%support
                gjk_closest(:, i) = gtr(i)%closest
                gjk_nsimp(i) = gtr(i)%nsimp
                gjk_reason(i) = gtr(i)%reason
            end do
            n_epa = min(ne, EPA_TRACE_MAX)
            do i = 1, n_epa
                epa_nml(:, i) = etr(i)%normal
                epa_dist(i) = etr(i)%dist
                epa_support(:, i) = etr(i)%support
                epa_nfaces(i) = etr(i)%nfaces
                epa_nverts(i) = etr(i)%nverts
                epa_reason(i) = etr(i)%reason
            end do
        end if
        deallocate(a, b)
    end subroutine gjkepa_query_c

    subroutine quickhull_c(pts, n, max_faces, faces, nfaces, info) &
            bind(c, name="quickhull_c")
        real(c_double), intent(in) :: pts(3, n)
        integer(c_int), intent(in), value :: n, max_faces
        real(c_double), intent(out) :: faces(3, 3, max_faces)
        integer(c_int), intent(out) :: nfaces, info
        real(8), allocatable :: p(:, :), mesh(:, :, :)
        integer :: i, nf

        nfaces = 0
        info = 1
        faces = 0.0_c_double
        if (n < 4) return
        allocate(p(n, 3))
        p = transpose(pts)
        call QuickHull(p, mesh, info)
        if (info /= 0 .or. .not. allocated(mesh)) then
            if (info == 0) info = 1
            deallocate(p)
            return
        end if
        nf = min(size(mesh, 1), max_faces)
        nfaces = nf
        do i = 1, nf
            faces(:, :, i) = transpose(mesh(i, :, :))
        end do
        if (size(mesh, 1) > max_faces) info = 2
        deallocate(p, mesh)
    end subroutine quickhull_c

    function gjkepa_abi_version() result(v) bind(c, name="gjkepa_abi_version")
        integer(c_int) :: v
        v = 1
    end function gjkepa_abi_version
end module gjkepa_capi
