! SLAE - Solution of Linear Algebraic Equations
! Using LU decomposition and backsubstitution routines, it...
!   - finds the inverse of a matrix column by column
!
!

program test_matrix

    use lu_factorization, only: dp, crout, backsubstitution

    implicit none

    integer, parameter :: n = 3
    real(dp) :: a(n, n), a_original(n, n), a_inv(n, n), produto(n, n), det
    integer  :: indx(n), d, info
    integer  :: i, j

    ! Matrix to be inverted
    a = reshape([2.0_dp, 4.0_dp, 1.0_dp, &
                 1.0_dp, 3.0_dp, 2.0_dp, &
                 1.0_dp, 5.0_dp, 4.0_dp], [n, n])
    a_original = a

    ! Set up identity matrix
    a_inv = 0.0_dp
    do i = 1, n
        a_inv(i, i) = 1.0_dp
    end do

    call crout(a, indx, d, info)
    if (info /= 0) error stop 'singular matrix'

    do j = 1, n
        call backsubstitution(a, indx, a_inv(:, j))
    end do

    ! Exact zero test is intentional: 
    ! only a truly null row makes the matrix singular
    print '(a)', '(a) Inverse of A:'
    do i = 1, n
        print '(3f10.4)', a_inv(i, :)
    end do

    ! Check: A . A^-1 should be the identity
    produto = matmul(a_original, a_inv)

    print '(a)', 'A . A^-1 (should be identity):'
    do i = 1, n
        print '(3f10.4)', produto(i, :)
    end do


    det = real(d, dp)
    do j = 1, n
        det = det * a(j, j)
    end do

    print '(a, f10.4)', '(b) Determinant of A = ', det

end program test_matrix


! Run in the CMD: 
! gfortran -Wall -Wextra -fcheck=all -g slae_lu_factorization.f90 slae_lu_factorization_tests.f90 -o slae_lu_factorization_tests

