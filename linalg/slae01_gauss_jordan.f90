! SLAE - Solution of Linear Algebraic Equations
! Linear equation solution by Gauss-Jordan elimination with full pivoting.
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

        ! A(n, n) - Matrix of coefficients
        ! It's the input, original matrix of the system.
        ! Calculation A.x = B or x = A^-1.B
        real(dp), intent(in out) :: A(:, :)
        ! B(n, m) - Matrix of independent terms
        ! It's the output
        real(dp), intent(in out) :: B(:, :)

        ! n, m - number of lines, columns
        integer :: n, m

        ! i - Main step counter
        ! j, k - Line, column index for pivot searching
        ! l - Final loop counter (interchanges inverse matrix columns)
        ! idx_row_reduced - Index of the lines to be reduced during elimination
        ! pivot_row, pivot_col - Pivot line, column
        integer :: i, j, k, l, idx_row_reduced, pivot_row, pivot_col

        ! save_row, save_col - Saves original line, column that was changed in step i
        ! idx_pivot - Matches if column j is the pivot
        integer, allocatable :: save_row(:), save_col(:), idx_pivot(:)

        ! big - Higher absolute value found on pivot searching.
        !       If big = 0, the matrix is singular
        ! elimination_factor - idx_row_reduced line multiplier relative to pivot line pivot_col.
        !                      Used to subtract elimination_factor x pivot line from line idx_row_reduced
        ! pivinv - (1/pivot) Used to normalize pivot line
        real(dp) :: big, elimination_factor, pivinv
        
        n = size(A, 1)
        m = size(B, 2)

        if (size(A, 2) /= n) error stop 'gauss_jordan: A must be squared'
        if (size(B, 1) /= n) error stop 'gauss_jordan: B has incompatible lines'

        allocate(idx_pivot(n), save_row(n), save_col(n))
        idx_pivot = 0
        save_row = 0
        save_col = 0

        do i = 1, n
            ! Searching the pivot
            big = 0.0_dp
            pivot_row = 0
            pivot_col = 0
            do j = 1, n
                if (idx_pivot(j) == 1) cycle
                do k = 1, n
                    if (idx_pivot(k) == 0) then
                        if (abs(A(j,k)) >= big) then
                            big  = abs(A(j,k))
                            pivot_row = j
                            pivot_col = k
                        end if
                    else if (idx_pivot(k) > 1) then
                        error stop 'gauss_jordan: matriz singular'
                    end if
                end do
            end do
            idx_pivot(pivot_col) = idx_pivot(pivot_col) + 1

            ! If necessary, change lines
            if (pivot_row /= pivot_col) then
                call swap_rows(A, pivot_row, pivot_col)
                call swap_rows(B, pivot_row, pivot_col)
            end if

            save_row(i) = pivot_row
            save_col(i) = pivot_col

            if (A(pivot_col,pivot_col) == 0.0_dp) error stop 'gauss_jordan: matriz singular'

            ! Pivot line normalization
            pivinv        = 1.0_dp / A(pivot_col,pivot_col)
            A(pivot_col,pivot_col)  = 1.0_dp
            A(pivot_col,:)     = A(pivot_col,:) * pivinv
            B(pivot_col,:)     = B(pivot_col,:) * pivinv

            ! Reducing the lines
            do idx_row_reduced = 1, n
                if (idx_row_reduced == pivot_col) cycle
                elimination_factor              = A(idx_row_reduced, pivot_col)
                A(idx_row_reduced, pivot_col)   = 0.0_dp
                A(idx_row_reduced, :)           = A(idx_row_reduced, :) - A(pivot_col, :) * elimination_factor
                B(idx_row_reduced, :)           = B(idx_row_reduced, :) - B(pivot_col, :) * elimination_factor
            end do
        end do

        ! Interchanging pairs of columns in the inverse
        do l = n, 1, -1
            if (save_row(l) /= save_col(l)) then
                call swap_cols(A, save_row(l), save_col(l))
            end if
        end do

        deallocate(idx_pivot, save_row, save_col)
    end subroutine gauss_jordan


    pure subroutine swap_rows(mat, i, j)
        real(dp), intent(in out)    ::  mat(:,:)
        integer,  intent(in)        ::  i, j
        real(dp)                    ::  tmp(size(mat,2))

        if (i == j) return
        tmp      = mat(i,:)
        mat(i,:) = mat(j,:)
        mat(j,:) = tmp

    end subroutine swap_rows


    pure subroutine swap_cols(mat, i, j)
        real(dp), intent(in out)    ::  mat(:,:)
        integer,  intent(in)        ::  i, j
        real(dp)                    ::  tmp(size(mat,1))

        if (i == j) return
        tmp      = mat(:,i)
        mat(:,i) = mat(:,j)
        mat(:,j) = tmp

    end subroutine swap_cols


end module linear_solvers

