! SLAE - Solution of Linear Algebraic Equations
! Solves a special case of linear equations systems such as tridiagonal systems,
! when matrix has nonzero elements only on the main diagonal plus or minus one column,
! and banded systems, when matrix has nonzero elements only on the main diagonal plus or minus one diagonal.
!

module special_slae
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private

    integer,  parameter :: dp = real64
    real(dp), parameter :: TOL = 100.0_dp * epsilon(1.0_dp)
    real(dp), parameter :: TINY = 1.0e-20_dp

    public :: tridiagonal
    public :: band_multiply
    public :: band_decomposition

contains

    ! Solves for a vector u(1:n) the tridiagonal linear set given by the equation below
    ! [b1 , c1 , 0  , ..., 0     , 0     , 0     ]   [r1    ]   [u1    ]
    ! [a2 , b2 , c2 , ..., 0     , 0     , 0     ]   [r2    ]   [u2    ]
    ! [..., ..., ..., ..., ...   , ...   , ...   ] . [...   ] = [...   ]
    ! [0  , 0  , 0  , ..., a(n-1), b(n-1), c(n-1)]   [r(n-1)]   [u(n-1)]
    ! [0  , 0  , 0  , ..., 0     , a(n)  , b(n)  ]   [r(n)  ]   [u(n)  ]
    subroutine tridiagonal(a, b, c, r, u)
        real(dp), intent(in) :: a(:), b(:), c(:), r(:)
        real(dp), intent(out) :: u(:)

        integer :: j, n
        real(dp) :: bet
        real(dp), allocatable :: gam(:)

        n = size(b)

        if (abs(b(1)) <= TOL) error stop 'tridiagonal: null pivot'

        allocate(gam(n))

        bet = b(1)
        u(1) = r(1) / bet

        ! decomposition and forward substitution
        do j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            
            if (abs(bet) <= TOL) error stop 'tridiagonal: singular matrix'

            u(j) = (r(j) - a(j) * u(j-1)) / bet
        end do

        ! backsubstitution
        do j = n-1 , 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
        end do

        return
    end subroutine tridiagonal

    ! Matrix multiply b = A . x, where A is band diagonal with m1 rows below the diagonal and m2 rows above.
    ! The input vector x and output vector b are stored as x(1: n) and b(1:n), respectively.
    ! The array a(1;n, 1:m1+m2+1) stores A as follows:
    ! The diagonal elements are in a a(1:n, m1+1).
    ! Subdiagonal elements are in a(j:n, 1:m1) (with j > 1 appropriate to teh number of elements on each subdiagonal).
    ! Superdiagonal elements are in a (1:j, m1+2:m1+m2+1) with j < n appropriate to the number of elements on each superdiagonal.
    subroutine band_multiply(a, x, b, m1, m2, n)
        real(dp), intent(in) :: a(:,:), x(:)
        real(dp), intent(out) :: b(:)
        integer, intent(in) :: m1, m2, n

        integer :: i, j, k, mp, np, jmin, jmax

        np = size(a, 1)
        mp = size(a, 2)

        do i = 1, n
            b(i) = 0.0_dp
            k = i - m1 - 1
            jmin = max(1, 1 - k)
            jmax = min(m1 + m2 + 1, n - k)
            do j = jmin, jmax
                b(i) = b(i) + a (i, j) * x(j + k)
            end do
        end do

        return
    end subroutine band_multiply

    ! LU decomposition of a banded matrix A with m1 subdiagonals and m2 superdiagonals,
    ! using Crout's method with partial pivoting.
    subroutine band_decomposition(a, al, indx, d, m1, m2, n)
        real(dp), intent(inout) :: a(:, :)
        real(dp), intent(out) :: al(:, :)
        integer, intent(out) :: indx(:)
        real(dp), intent(out) :: d
        integer, intent(in) :: m1, m2, n


        integer :: i, j, k, l, mm, np, mp, mpl, ipiv
        real(dp) :: dum

        np = size(a, 1)
        mp = size(a, 2)
        mpl = size(al, 2)
        mm = m1 + m2 + 1

        ! Moving bands to positioning main diagonal in the first column of the array
        l = m1
        do i = 1, m1
            do j = m1 + 2 - i, mm
                a(i, j - l) = a(i, j)
            end do
            l = l - 1
            do j = mm - l, mm
                a(i, j) = 0.0_dp
            end do
        end do

        ! Main loop: Gauss elimination with partial pivoting
        d = 1.0_dp
        l = m1
        do k = 1, n
            dum = a(k, 1)
            ipiv = k
            if (1 < n) l = l + 1
            ! Searches the pivot
            do j = k + 1, l
                if (abs(a(j, 1)) > abs(dum)) then
                    dum = a(j, 1)
                    ipiv = j
                end if
            end do
            
            indx(k) = ipiv

            ! singular matrix - Keep going, but set pivot to TINY
            if (abs(dum) <= TINY) a(k, 1) = TINY

            ! Changing lines
            if (ipiv /= k) then
                d = -d
                do j = 1, mm
                    dum = a(k, j)
                    a(k, j) = a(ipiv, j)
                    a(ipiv, j) = dum
                end do
            end if

            ! Elimination
            do i = j + 1, l
                dum = a(i, 1) / a(k, 1)
                al(k, i - k) = dum
                do j = 2, mm
                    a(i, j - 1) = a(i, j) - dum * a(k, j)
                end do
                a(i, mm) = 0.0_dp
            end do
        end do

    end subroutine band_decomposition

    subroutine band_backsubstitution(a, al, indx, b, m1, m2, n)
        real(dp), intent(in) :: a (:, :), al(:, :)
        integer, intent(in) :: indx(:)
        real(dp), intent(inout) :: b(:)
        integer, intent(in) :: m1, m2, n

        integer :: i, k, l, mm, np, mp, mpl
        real(dp) :: dum

        np  = size(a, 1)
        mp  = size(a, 2)
        mpl = size(al, 2)
        mm  = m1 + m2 + 1

        ! Direct subsitution (L . y = P . b)
        l = m1
        do k = 1, n
            i = indx(k)
            if (i /= k) then
                dum   = b(k)
                b(k)  = b(i)
                b(i)  = dum
            end if
            if (l < n) l = l + 1
            do i = k + 1, l
                b(i) = b(i) - al(k, i - k) * b(k)
            end do
        end do

        ! Reverse substitution (U . x = y)
        l = 1
        do i = n, 1, -1
            dum = b(i)
            do k = 2, l
                dum = dum - a(i, k) * b(k + i - 1)
            end do
            b(i) = dum / a(i, 1)
            if (l < mm) l = l + 1
        end do

    end subroutine band_backsubstitution

end module special_slae