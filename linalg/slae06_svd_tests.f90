! SLAE - Solution of Linear Algebraic Equations
! SVD comprises the decomposition of a matrix A(m, n), where m >= n, into U, W and V, such that A = U * W * V^T
!
!
!

program test_svd

    use svd, only: dp, solve_backsubstitution, solve_decomposition

    implicit none

    integer, parameter :: m = 3, n = 3
    integer, parameter :: mp = 3, np = 3
    real(dp) :: A(np, np), W(np), V(np, np), b(mp), x(np)
    real(dp) :: Wmax, Wmin
    integer :: i, j, info

    A = reshape([1.0_dp, 2.0_dp, 3.0_dp, &
                 4.0_dp, 5.0_dp, 6.0_dp, &
                 7.0_dp, 8.0_dp, 10.0_dp ], shape(A))
    b = [1.0_dp, 2.0_dp, 3.0_dp]

    ! A is subscribed with U
    call solve_decomposition(A, m, n, mp, np, W, V, info)

    Wmax = 0.0_dp
    do j = 1, n
        if (W(j) >= Wmax) Wmax = W(j)
    end do

    Wmin = Wmax * 1.0e-6_dp
    do j = 1, n
        if (W(j) <= Wmin) W(j) = 0.0_dp
    end do

    call solve_backsubstitution(A, W, V, m, n, mp, np, b, x)

    print *, "Solution vector x: "
    do i = 1, n
        print *, i, x(i)
    end do

end program test_svd


! Run in the CMD: 
! gfortran -Wall -Wextra -fcheck=all -g slae06_svd.f90 slae06_svd_tests.f90 -o svd_test

