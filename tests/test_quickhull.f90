program test_quickhull
    use GCLIB_QuickHull, only: QuickHull
    implicit none

    integer :: nfail
    nfail = 0
    call check_tetra()
    call check_cube()
    call check_duplicates()
    call check_too_few()

    if (nfail == 0) then
        write(*, '(A)') 'test_quickhull: all passed'
        stop 0
    else
        write(*, '(A,I0,A)') 'test_quickhull: ', nfail, ' failed'
        stop 1
    end if

contains

    subroutine fail(msg)
        character(*), intent(in) :: msg
        nfail = nfail + 1
        write(*, '(A,A)') 'FAIL: ', msg
    end subroutine fail

    subroutine check_tetra()
        real(8) :: p(4, 3)
        real(8), allocatable :: mesh(:, :, :)
        integer :: info
        p(1, :) = [0.d0, 0.d0, 0.d0]
        p(2, :) = [1.d0, 0.d0, 0.d0]
        p(3, :) = [0.d0, 1.d0, 0.d0]
        p(4, :) = [0.d0, 0.d0, 1.d0]
        call QuickHull(p, mesh, info)
        if (info /= 0) call fail('tetra info')
        if (.not. allocated(mesh)) then
            call fail('tetra mesh missing')
            return
        end if
        if (size(mesh, 1) /= 4) call fail('tetra should have 4 faces')
        if (size(mesh, 2) /= 3 .or. size(mesh, 3) /= 3) call fail('tetra mesh shape')
    end subroutine check_tetra

    subroutine check_cube()
        real(8) :: p(8, 3)
        real(8), allocatable :: mesh(:, :, :)
        integer :: info, i, ix, iy, iz
        i = 0
        do iz = 0, 1
            do iy = 0, 1
                do ix = 0, 1
                    i = i + 1
                    p(i, :) = [dble(ix), dble(iy), dble(iz)]
                end do
            end do
        end do
        call QuickHull(p, mesh, info)
        if (info /= 0) call fail('cube info')
        if (.not. allocated(mesh)) then
            call fail('cube mesh missing')
            return
        end if
        ! cube convex hull is 12 triangles
        if (size(mesh, 1) /= 12) call fail('cube should have 12 triangles')
    end subroutine check_cube

    subroutine check_duplicates()
        real(8) :: p(6, 3)
        real(8), allocatable :: mesh(:, :, :)
        integer :: info
        p(1, :) = [0.d0, 0.d0, 0.d0]
        p(2, :) = [1.d0, 0.d0, 0.d0]
        p(3, :) = [0.d0, 1.d0, 0.d0]
        p(4, :) = [0.d0, 0.d0, 1.d0]
        p(5, :) = [0.d0, 0.d0, 0.d0]
        p(6, :) = [1.d0, 0.d0, 0.d0]
        call QuickHull(p, mesh, info)
        if (info /= 0) call fail('duplicate-point tetra should still hull')
        if (allocated(mesh)) then
            if (size(mesh, 1) /= 4) call fail('duplicate tetra face count')
        end if
    end subroutine check_duplicates

    subroutine check_too_few()
        real(8) :: p(3, 3)
        real(8), allocatable :: mesh(:, :, :)
        integer :: info
        p(1, :) = [0.d0, 0.d0, 0.d0]
        p(2, :) = [1.d0, 0.d0, 0.d0]
        p(3, :) = [0.d0, 1.d0, 0.d0]
        call QuickHull(p, mesh, info)
        if (info == 0) call fail('three points must fail')
    end subroutine check_too_few

end program test_quickhull
