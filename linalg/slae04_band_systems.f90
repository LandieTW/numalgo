! SLAE - Solution of Linear Algebraic Equations
! Solves special cases of linear equations systems such as tridiagonal and banded systems,
! tridiagonal systems - matrix has nonzero elements only on the main diagonal plus or minus one column,
! banded systems - matrix has nonzero elements only on the main diagonal plus or minus one diagonal.
!

module special_slae
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    integer,  parameter :: dp = real64
    real(dp), parameter :: TOL = 100.0_dp * epsilon(1.0_dp)     ! Tolerance
    real(dp), parameter :: TINY = 1.0e-20_dp                    ! small value used to substitute pivot in a singular matrix

    public :: tridiagonal
    public :: band_multiply
    public :: band_decomposition
    public :: band_backsubstitution

contains

    ! Solves tridiagonal systems (A.u = r) by Gauss Elimination
    subroutine tridiagonal(a, b, c, r, u)

        ! INPUTS
        ! a - Subdiagonal (elements a(2:n))
        ! b - Main diagonal (elements b(1:n))
        ! c - Superdiagonal (elements c(1:n-1))
        ! r - Independent vector (RHS)
        real(dp), intent(in) :: a(:), b(:), c(:), r(:)
        
        ! OUTPUT
        ! u - Solution vector
        real(dp), intent(out) :: u(:)

        ! LOCALS
        ! j - loop counter
        ! n - number of equations
        ! curr_pivot - Current pivot
        integer :: j, n
        real(dp) :: curr_pivot
        ! Store factors (c(j-1)/curr_pivot) for backsubstitution
        real(dp), allocatable :: gam(:)

        n = size(b)

        if (n < 1) error stop 'tridiagonal: empty system'
        if (abs(b(1)) <= TOL) error stop 'tridiagonal: null pivot'

        allocate(gam(n))

        curr_pivot = b(1)
        u(1) = r(1) / curr_pivot

        ! decomposition and forward substitution
        do j = 2, n
            gam(j) = c(j-1) / curr_pivot
            curr_pivot = b(j) - a(j) * gam(j)
            
            if (abs(curr_pivot) <= TOL) error stop 'tridiagonal: singular matrix'

            u(j) = (r(j) - a(j) * u(j-1)) / curr_pivot
        end do

        ! backsubstitution
        do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
        end do

        deallocate(gam)

        return
    end subroutine tridiagonal


    ! Calculates b = A.x for band matrix A stored in compact format
    subroutine band_multiply(a, x, b, m1, m2, n)

        ! INPUT
        ! a - band matrix stored at a(1:n, 1:m1+m2+1)
        ! x - Input vector
        real(dp), intent(in) :: a(:,:), x(:)
        ! m1, m2 - Number of subdiagonals and superdiagonals
        ! n - Number of equations
        integer, intent(in) :: m1, m2, n

        ! OUTPUT
        ! b - Result vector
        real(dp), intent(out) :: b(:)

        ! LOCALS
        ! i - row counter (equation index)
        ! j - column counter (position within the compact band row)
        ! col_offset - offset between compact column j and global column of A (j + col_offset = global col)
        ! num_rows, num_cols - a matrix dimensions
        ! first_col, last_col - valid compact-column range for row i
        integer :: i, j, col_offset, num_rows, num_cols, first_col, last_col

        num_rows = size(a, 1)
        num_cols = size(a, 2)

        do i = 1, n
            b(i) = 0.0_dp
            ! Global column of A that compact column j represents:  global_col = j + col_offset
            col_offset = i - m1 - 1
            ! Ensure 1 <= j + col_offset <= n
            first_col = max(1, 1 - col_offset)
            last_col  = min(m1 + m2 + 1, n - col_offset)
            do j = first_col, last_col
                b(i) = b(i) + a(i, j) * x(j + col_offset)
            end do
        end do

        return
    end subroutine band_multiply


    ! LU decomposition of a banded matrix A by Crout with partial pivoting.
    subroutine band_decomposition(a, al, indx, d, m1, m2, n)

        ! INPUT / OUTPUT
        ! a - Band matrix
        real(dp), intent(inout) :: a(:, :)
        
        ! OUTPUT
        ! al - L factors (multipliers)
        real(dp), intent(out) :: al(:, :)
        ! indx - Permutation vector (pivots)
        integer, intent(out) :: indx(:)
        ! d - Row changing Signals
        real(dp), intent(out) :: d
        
        ! INPUT
        ! m1, m2 - Band matrix dimensions
        ! n - Number of equations
        integer, intent(in) :: m1, m2, n

        ! LOCALS
        ! i, k - row / elimination-step counters
        ! j - column counter
        ! lower_band_limit - last row index (within current column window) reachable for pivoting/elimination
        ! band_width - total number of diagonals (m1 + m2 + 1)
        ! num_rows, num_cols, num_l_cols - matrix dimensions
        ! pivot_row - row chosen as pivot
        ! pivot_value - candidate/current pivot value
        ! temp - temporary variable (pivot/changing)
        integer :: i, j, k, lower_band_limit, band_width, num_rows, num_cols, num_l_cols, pivot_row
        real(dp) :: pivot_value, temp

        num_rows   = size(a, 1)
        num_cols   = size(a, 2)
        num_l_cols = size(al, 2)
        band_width = m1 + m2 + 1

        ! Moving bands to positioning main diagonal in the first column of the array
        lower_band_limit = m1
        do i = 1, m1
            do j = m1 + 2 - i, band_width
                a(i, j - lower_band_limit) = a(i, j)
            end do
            lower_band_limit = lower_band_limit - 1
            do j = band_width - lower_band_limit, band_width
                a(i, j) = 0.0_dp
            end do
        end do

        ! Main loop: Gauss elimination with partial pivoting
        d = 1.0_dp
        lower_band_limit = m1
        do k = 1, n
            pivot_value = a(k, 1)
            pivot_row   = k
            if (1 < n) lower_band_limit = lower_band_limit + 1

            ! Search the pivot within the reachable band window of the current column
            do j = k + 1, lower_band_limit
                if (abs(a(j, 1)) > abs(pivot_value)) then
                    pivot_value = a(j, 1)
                    pivot_row   = j
                end if
            end do

            indx(k) = pivot_row

            ! singular matrix - Keep going, but set pivot to TINY
            if (abs(pivot_value) <= TINY) then
                a(k, 1)     = TINY
                pivot_value = TINY
            end if

            ! Changing lines
            if (pivot_row /= k) then
                d = -d
                do j = 1, band_width
                    temp             = a(k, j)
                    a(k, j)          = a(pivot_row, j)
                    a(pivot_row, j)  = temp
                end do
            end if

            ! Elimination - the loop must start at k + 1, not at (j + 1) left over from the pivot search
            do i = k + 1, lower_band_limit
                temp        = a(i, 1) / a(k, 1)
                al(k, i - k) = temp
                do j = 2, band_width
                    a(i, j - 1) = a(i, j) - temp * a(k, j)
                end do
                a(i, band_width) = 0.0_dp
            end do
        end do

    end subroutine band_decomposition


    ! Solves system after decomposition, with direct and inverse substitution
    subroutine band_backsubstitution(a, al, indx, b, m1, m2, n)

        ! INPUT
        ! a - Decomposed matrix
        ! al - L factors (multipliers)
        real(dp), intent(in) :: a (:, :), al(:, :)
        ! indx - stores the permutations
        integer, intent(in) :: indx(:)

        ! IN-OUTPUT
        ! b - RHS input, solution output
        real(dp), intent(inout) :: b(:)
        
        ! INPUT
        ! Band parameters (m1, m2 - dimensions; n - number of equations)
        integer, intent(in) :: m1, m2, n

        ! LOCALS
        ! i, k - row / step counters
        ! lower_band_limit - last row index reachable in the current window
        ! band_width - total number of diagonals (m1 + m2 + 1)
        ! num_rows, num_cols, num_l_cols - matrix dimensions
        ! swap_index - row selected for pivoting at step k
        ! temp - temporary variable (pivot/changing)
        integer :: i, k, lower_band_limit, band_width, num_rows, num_cols, num_l_cols, swap_index
        real(dp) :: temp

        num_rows   = size(a, 1)
        num_cols   = size(a, 2)
        num_l_cols = size(al, 2)
        band_width = m1 + m2 + 1

        ! Direct substitution (L . y = P . b)
        lower_band_limit = m1
        do k = 1, n
            swap_index = indx(k)
            if (swap_index /= k) then
                temp        = b(k)
                b(k)        = b(swap_index)
                b(swap_index) = temp
            end if
            if (lower_band_limit < n) lower_band_limit = lower_band_limit + 1
            do i = k + 1, lower_band_limit
                b(i) = b(i) - al(k, i - k) * b(k)
            end do
        end do

        ! Reverse substitution (U . x = y)
        ! For row i, U has at most min(band_width, n - i + 1) entries starting at column i.
        lower_band_limit = 1
        do i = n, 1, -1
            temp = b(i)
            do k = 2, lower_band_limit
                temp = temp - a(i, k) * b(k + i - 1)
            end do
            b(i) = temp / a(i, 1)
            if (lower_band_limit < band_width) lower_band_limit = lower_band_limit + 1
        end do

    end subroutine band_backsubstitution

end module special_slae

