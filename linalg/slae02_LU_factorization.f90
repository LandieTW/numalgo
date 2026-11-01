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
        real(dp), intent(inout) :: A(:, :)      ! Input Matrix A(n, n)
        integer,  intent(out)   :: indx(:)      ! Interchanging rows register (record which row was changed by j row in j step)
        integer,  intent(out)   :: d            ! Interchanging parity or signal
                                                ! Used to calculate the det(A) = d ⋅ ∏U(i,i) - +1 if even number of changes, -1 if odd number
        integer,  intent(out)   :: info         ! 0 = ok; i > 0 = row i is zero (singular)

        ! i, j - row, column counter
        ! n_rows - number of matrix equations
        ! imax - pivot row index (pivot is the max value in j column)
        integer  :: i, j, n_rows, imax
        ! aamax - High value found in Gauss elimination
        ! dum - pivot candidate (temporary value)
        real(dp) :: aamax, dum
        ! vv - scale factors (used to compare normalised pivots)
        real(dp) :: vv(size(A, 1))
        ! tmp_row - temporary row used to change two rows during pivoting
        real(dp) :: tmp_row(size(A, 2))

        n_rows  = size(A, 1)
        d       = 1
        info    = 0

        ! Loop over rows to get the implicit scaling information
        do i = 1, n_rows
            aamax = maxval(abs(A(i, :)))
            if (aamax <= tiny(1.0_dp)) then
                info = i                ! No nonzero largest element
                return
            end if
            vv(i) = 1.0_dp / aamax      ! Save scaling
        end do

        ! Loop over columns
        do j = 1, n_rows

            ! Cases (i < j)
            do i = 1, j - 1
                A(i, j) = A(i, j) - dot_product(A(i, 1:i-1), A(1:i-1, j))
            end do

            ! Cases i >= j: possible pivots
            aamax = 0.0_dp
            imax  = j
            do i = j, n_rows
                A(i, j) = A(i, j) - dot_product(A(i, 1:j-1), A(1:j-1, j))
                dum = vv(i) * abs(A(i, j))
                if (dum >= aamax) then
                    imax  = i
                    aamax = dum
                end if
            end do

            ! Interchanging rows, if necessary
            if (j /= imax) then
                tmp_row  = A(imax, :)
                A(imax, :) = A(j, :)
                A(j, :)    = tmp_row
                d          = -d
                vv(imax)   = vv(j)
            end if
            indx(j) = imax

            ! Nule pivot
            if (abs(A(j, j)) <= min_pivot) A(j, j) = sign(min_pivot, A(j, j))

            ! Divide the column below the diagonal by the pivot
            if (j < n_rows) A(j+1:n_rows, j) = A(j+1:n_rows, j) / A(j, j)
        end do
    end subroutine crout



    ! Solves A.x = b using the LU decomposition computed by crout.
    ! It takes into account that b may begin with many zero elements, which makes it efficient for matrix inversion.
    subroutine backsubstitution(a, indx, b)
        real(dp), intent(in)    :: A(:, :)      ! LU matrix from Crout
        integer,  intent(in)    :: indx(:)      ! permutation record from Crout
        real(dp), intent(inout) :: b(:)         ! independent term vector

        ! p - permutation counter
        ! b_i - match the first index of a non-nule vector b
        ! pr_i - permuted row, used to unshufle the permutation
        ! n_rows - number of matrix equations
        integer  :: p, b_i, pr_i, n_rows
        ! total - temporary accumulator (right-side) - L.y = P.b
        real(dp) :: total

        n_rows  = size(b)
        b_i = 0

        ! Forward substitution unscrambling the permutation
        do p = 1, n_rows
            pr_i    = indx(p)
            total   = b(pr_i)
            b(pr_i) = b(p)
            if (b_i /= 0) then
                total = total - dot_product(A(p, b_i:p-1), b(b_i:p-1))
            else if (abs(total) > tiny(1.0_dp)) then
                b_i = p
            end if
            b(p) = total
        end do

        ! Back substitution
        do p = n_rows, 1, -1
            b(p) = (b(p) - dot_product(A(p, p+1:n_rows), b(p+1:n_rows))) / A(p, p)
        end do
    end subroutine backsubstitution

end module lu_factorization

