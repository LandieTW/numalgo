! SLAE - Solution of Linear Algebraic Equations
! Linear equation solution by Gauss-Jordan elimination, equation.
! 
!
!

module linear_solvers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: gauss_jordan

contains
    subroutine gauss_jordan(a, b)

        real(dp), intent(in out) :: a(:, :)
        real(dp), intent(in out) :: b(:, :)

        integer :: n, m
        integer :: i, j, k, l, ll, irow, icol
        integer, allocatable :: indxr(:), indxc(:), ipiv(:)
        real(dp) :: big, dum, pivinv
        
        n = size(a, 1)
        m = size(b, 2)

        if (size(a, 2) /= n) error stop 'gauss_jordan: A must be squared'
        if (size(b, 1) /= n) error stop 'gauss_jordan: B has incompatible lines'

        allocate(ipiv(n), indxr(n), indxc(n))
        ipiv = 0
        indxr = 0
        indxc = 0

        do i = 1, n
            ! Searching the pivot
            big = 0.0_dp
            irow = 0
            icol = 0
            do j = 1, n
                if (ipiv(j) == 1) cycle
                do k = 1, n
                    if (ipiv(k) == 0) then
                        if (abs(a(j,k)) >= big) then
                            big  = abs(a(j,k))
                            irow = j
                            icol = k
                        end if
                    else if (ipiv(k) > 1) then
                        error stop 'gauss_jordan: matriz singular'
                    end if
                end do
            end do
            ipiv(icol) = ipiv(icol) + 1

            ! If necessary, change lines
            if (irow /= icol) then
                call swap_rows(a, irow, icol)
                call swap_rows(b, irow, icol)
            end if

            indxr(i) = irow
            indxc(i) = icol

            if (a(icol,icol) == 0.0_dp) error stop 'gauss_jordan: matriz singular'

            ! Pivot line normalization
            pivinv        = 1.0_dp / a(icol,icol)
            a(icol,icol)  = 1.0_dp
            a(icol,:)     = a(icol,:) * pivinv
            b(icol,:)     = b(icol,:) * pivinv

            ! Reducing other lines
            do ll = 1, n
                if (ll == icol) cycle
                dum        = a(ll,icol)
                a(ll,icol) = 0.0_dp
                a(ll,:)    = a(ll,:) - a(icol,:) * dum
                b(ll,:)    = b(ll,:) - b(icol,:) * dum
            end do
        end do

        ! Interchanging pairs of columns in the inverse
        do l = n, 1, -1
            if (indxr(l) /= indxc(l)) then
                call swap_cols(a, indxr(l), indxc(l))
            end if
        end do

        deallocate(ipiv, indxr, indxc)
    end subroutine gauss_jordan


    pure subroutine swap_rows(mat, i, j)
        real(dp), intent(in out) :: mat(:,:)
        integer,  intent(in)     :: i, j
        real(dp) :: tmp(size(mat,2))
        if (i == j) return
        tmp      = mat(i,:)
        mat(i,:) = mat(j,:)
        mat(j,:) = tmp
    end subroutine swap_rows


    pure subroutine swap_cols(mat, i, j)
        real(dp), intent(in out) :: mat(:,:)
        integer,  intent(in)     :: i, j
        real(dp) :: tmp(size(mat,1))
        if (i == j) return
        tmp      = mat(:,i)
        mat(:,i) = mat(:,j)
        mat(:,j) = tmp
    end subroutine swap_cols


end module linear_solvers

