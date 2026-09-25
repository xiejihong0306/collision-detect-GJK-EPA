program test_gjkepa
    use GCLIB_GJKEPA, only: GJKEPA
    use geom3d, only: rk, rotate_axis_angle
    implicit none

    integer :: nfail
    nfail = 0

    call case_separated_cubes()
    call case_overlap_cubes()
    call case_deep_overlap()
    call case_face_touch()
    call case_vertex_in_face()
    call case_edge_edge()
    call case_tetrahedra()
    call case_near_miss()
    call case_duplicate_verts()
    call case_rotated_overlap()

    if (nfail == 0) then
        write(*, '(A)') 'test_gjkepa: all passed'
        stop 0
    else
        write(*, '(A,I0,A)') 'test_gjkepa: ', nfail, ' failed'
        stop 1
    end if

contains

    subroutine fail(msg)
        character(*), intent(in) :: msg
        nfail = nfail + 1
        write(*, '(A,A)') 'FAIL: ', msg
    end subroutine fail

    subroutine query(a, b, hit, kind, nml, depth, point)
        real(rk), intent(in) :: a(:, :), b(:, :)
        logical, intent(out) :: hit
        integer, intent(out) :: kind
        real(rk), intent(out) :: nml(3), depth, point(3)
        real(rk) :: nearest(2, 3)
        call GJKEPA(2, 1.0_rk, a, b, hit, kind, nearest, nml, point, depth)
    end subroutine query

    function unit_cube(shift) result(p)
        real(rk), intent(in) :: shift(3)
        real(rk) :: p(8, 3)
        integer :: i, ix, iy, iz
        i = 0
        do iz = 0, 1
            do iy = 0, 1
                do ix = 0, 1
                    i = i + 1
                    p(i, :) = shift + [real(ix, rk), real(iy, rk), real(iz, rk)]
                end do
            end do
        end do
    end function unit_cube

    subroutine case_separated_cubes()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([2.0_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (hit) call fail('separated cubes should not collide')
        if (kind /= 0) call fail('separated cubes type')
    end subroutine case_separated_cubes

    subroutine case_overlap_cubes()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([0.5_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('overlapping cubes should collide')
        if (abs(depth - 0.5_rk) > 2.0e-3_rk) then
            write(*, '(A,F16.8)') '  depth=', depth
            call fail('overlap depth should be ~0.5')
        end if
        if (abs(abs(nml(1)) - 1.0_rk) > 2.0e-2_rk) call fail('overlap normal should be +/-x')
    end subroutine case_overlap_cubes

    subroutine case_deep_overlap()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([0.1_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('deep overlap should collide')
        if (abs(depth - 0.9_rk) > 2.0e-2_rk) then
            write(*, '(A,F16.8)') '  deep depth=', depth
            call fail('deep overlap depth should be ~0.9')
        end if
    end subroutine case_deep_overlap

    subroutine case_face_touch()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([1.0_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('face-touch should be a contact')
        if (depth > 2.0e-2_rk) call fail('face-touch depth should be ~0')
        if (kind /= 2) call fail('face-touch should be type 2')
    end subroutine case_face_touch

    subroutine case_vertex_in_face()
        real(rk) :: cube(8, 3), tet(4, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        cube = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        tet(1, :) = [0.5_rk, 0.5_rk, 0.7_rk]
        tet(2, :) = [1.5_rk, 0.2_rk, 1.4_rk]
        tet(3, :) = [1.5_rk, 0.8_rk, 1.4_rk]
        tet(4, :) = [0.8_rk, 0.5_rk, 1.8_rk]
        call query(cube, tet, hit, kind, nml, depth, point)
        if (.not. hit) call fail('vertex-in-face tetra should collide')
        if (depth <= 0.0_rk) call fail('vertex-in-face should report positive depth')
    end subroutine case_vertex_in_face

    subroutine case_edge_edge()
        real(rk) :: a(4, 3), b(4, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        ! two tetrahedra whose closest features are crossing edges
        a(1, :) = [0.0_rk, -0.2_rk, 0.0_rk]
        a(2, :) = [0.0_rk,  0.2_rk, 0.0_rk]
        a(3, :) = [-0.4_rk, 0.0_rk, -0.4_rk]
        a(4, :) = [-0.4_rk, 0.0_rk,  0.4_rk]
        b(1, :) = [-0.2_rk, 0.0_rk, 0.05_rk]
        b(2, :) = [ 0.2_rk, 0.0_rk, 0.05_rk]
        b(3, :) = [0.0_rk, -0.4_rk, 0.45_rk]
        b(4, :) = [0.0_rk,  0.4_rk, 0.45_rk]
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('crossing-edge tetrahedra should collide')
    end subroutine case_edge_edge

    subroutine case_tetrahedra()
        real(rk) :: a(4, 3), b(4, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a(1, :) = [0.0_rk, 0.0_rk, 0.0_rk]
        a(2, :) = [1.0_rk, 0.0_rk, 0.0_rk]
        a(3, :) = [0.0_rk, 1.0_rk, 0.0_rk]
        a(4, :) = [0.0_rk, 0.0_rk, 1.0_rk]
        b(1, :) = [0.2_rk, 0.2_rk, 0.2_rk]
        b(2, :) = [1.2_rk, 0.2_rk, 0.2_rk]
        b(3, :) = [0.2_rk, 1.2_rk, 0.2_rk]
        b(4, :) = [0.2_rk, 0.2_rk, 1.2_rk]
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('overlapping tetrahedra should collide')
        if (depth <= 0.0_rk) call fail('tetra depth')
    end subroutine case_tetrahedra

    subroutine case_near_miss()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([1.05_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (hit) call fail('1.05-separated cubes should miss')
    end subroutine case_near_miss

    subroutine case_duplicate_verts()
        real(rk) :: a(10, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind
        a(1:8, :) = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        a(9, :) = a(1, :)
        a(10, :) = a(2, :)
        b = unit_cube([0.4_rk, 0.0_rk, 0.0_rk])
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('duplicate-vertex cube should still collide')
    end subroutine case_duplicate_verts

    subroutine case_rotated_overlap()
        real(rk) :: a(8, 3), b(8, 3), nml(3), depth, point(3)
        logical :: hit
        integer :: kind, i
        a = unit_cube([0.0_rk, 0.0_rk, 0.0_rk])
        b = unit_cube([0.35_rk, 0.15_rk, 0.05_rk])
        do i = 1, 8
            b(i, :) = rotate_axis_angle(b(i, :), [0.2_rk, 0.8_rk, 0.3_rk], 27.0_rk)
        end do
        call query(a, b, hit, kind, nml, depth, point)
        if (.not. hit) call fail('rotated overlapping cubes should collide')
        if (norm2(nml) < 0.5_rk) call fail('rotated overlap normal vanished')
        if (depth <= 0.0_rk) call fail('rotated overlap depth')
    end subroutine case_rotated_overlap

end program test_gjkepa
