program two_polyhedra
    use GCLIB_GJKEPA, only: GJKEPA
    use geom3d, only: rk, rotate_axis_angle
    implicit none

    real(rk) :: cube_a(8, 3), cube_b(8, 3)
    logical :: hit
    integer :: kind, i
    real(rk) :: nearest(2, 3), nml(3), point(3), depth

    cube_a = cube_at([0.0_rk, 0.0_rk, 0.0_rk], 1.0_rk)
    cube_b = cube_at([0.4_rk, 0.1_rk, 0.0_rk], 1.0_rk)
    do i = 1, 8
        cube_b(i, :) = rotate_axis_angle(cube_b(i, :), [0.0_rk, 0.0_rk, 1.0_rk], 18.0_rk)
    end do

    call GJKEPA(2, 1.0_rk, cube_a, cube_b, hit, kind, nearest, nml, point, depth)

    write(*, '(A)') 'GJK-EPA example: unit cube vs rotated overlapping cube'
    write(*, '(A,L1)') '  collision        = ', hit
    write(*, '(A,I0)') '  colli_type       = ', kind
    write(*, '(A,3F12.6)') '  normal           = ', nml
    write(*, '(A,F12.6)') '  penetration      = ', depth
    write(*, '(A,3F12.6)') '  contact point    = ', point
    write(*, '(A,3F12.6)') '  nearest on A     = ', nearest(1, :)
    write(*, '(A,3F12.6)') '  nearest on B     = ', nearest(2, :)

contains

    function cube_at(origin, edge) result(p)
        real(rk), intent(in) :: origin(3), edge
        real(rk) :: p(8, 3)
        integer :: i, ix, iy, iz
        i = 0
        do iz = 0, 1
            do iy = 0, 1
                do ix = 0, 1
                    i = i + 1
                    p(i, :) = origin + edge * [real(ix, rk), real(iy, rk), real(iz, rk)]
                end do
            end do
        end do
    end function cube_at

end program two_polyhedra
