! SLAE - Solution of Linear Algebraic Equations
! L. U. Crout's Method with partial pivoting
! It initially finds the largest element in each row, and subsequently scales the comparison as if we had 
! initially scaled all the equations to make their maximum coefficient equal to unity
! 

module lu_factorization
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: dp, crout, backsubstitution

    integer,  parameter :: dp = real64
    real(dp), parameter :: min_pivot = 1.0e-20_dp

contains

    ! Given a matrix with physical dimensions, this routine replaces it by the LU decomposition of a rowwise permutation of itself.
    ! This routine is used in combination with backsubstitution to solve linear equations or invert a matrix.
    subroutine crout(a, indx, d, info)
        real(dp), intent(inout) :: a(:, :)
        integer,  intent(out)   :: indx(:)   ! Interchanging rows register
        integer,  intent(out)   :: d         ! Interchanging parity
        integer,  intent(out)   :: info      ! 0 = ok; i > 0 = row i is zero (singular)

        integer  :: i, j, n, imax
        real(dp) :: aamax, dum
        real(dp) :: vv(size(a, 1))
        real(dp) :: tmp_row(size(a, 2))

        n    = size(a, 1)
        d    = 1
        info = 0

        ! Loop over rows to get the implicit scaling information
        do i = 1, n
            aamax = maxval(abs(a(i, :)))
            if (aamax == 0.0_dp) then
                info = i                ! No nonzero largest element
                return
            end if
            vv(i) = 1.0_dp / aamax      ! Save scaling
        end do

        ! Loop over columns
        do j = 1, n

            ! Cases (i < j)
            do i = 1, j - 1
                a(i, j) = a(i, j) - dot_product(a(i, 1:i-1), a(1:i-1, j))
            end do

            ! Cases i >= j: possible pivots
            aamax = 0.0_dp
            imax  = j
            do i = j, n
                a(i, j) = a(i, j) - dot_product(a(i, 1:j-1), a(1:j-1, j))
                dum = vv(i) * abs(a(i, j))
                if (dum >= aamax) then
                    imax  = i
                    aamax = dum
                end if
            end do

            ! Interchanging rows, if necessary
            if (j /= imax) then
                tmp_row  = a(imax, :)
                a(imax, :) = a(j, :)
                a(j, :)    = tmp_row
                d          = -d
                vv(imax)   = vv(j)
            end if
            indx(j) = imax

            ! Nule pivot
            if (a(j, j) == 0.0_dp) a(j, j) = min_pivot

            ! Divide the column below the diagonal by the pivot
            if (j < n) a(j+1:n, j) = a(j+1:n, j) / a(j, j)
        end do
    end subroutine crout



    ! Solves A.x = b using the LU decomposition computed by crout.
    ! It takes into account that b may begin with many zero elements, which makes it efficient for matrix inversion.
    subroutine backsubstitution(a, indx, b)
        real(dp), intent(in)    :: a(:, :)   ! LU matrix from crout
        integer,  intent(in)    :: indx(:)   ! permutation record from crout
        real(dp), intent(inout) :: b(:)

        integer  :: i, ii, ll, n
        real(dp) :: total

        n  = size(b)
        ii = 0

        ! Forward substitution unscrambling the permutation
        do i = 1, n
            ll    = indx(i)
            total = b(ll)
            b(ll) = b(i)
            if (ii /= 0) then
                total = total - dot_product(a(i, ii:i-1), b(ii:i-1))
            else if (total /= 0.0_dp) then
                ii = i
            end if
            b(i) = total
        end do

        ! Back substitution
        do i = n, 1, -1
            b(i) = (b(i) - dot_product(a(i, i+1:n), b(i+1:n))) / a(i, i)
        end do
    end subroutine backsubstitution

end module lu_factorization

