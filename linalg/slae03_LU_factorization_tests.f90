! SLAE - Solution of Linear Algebraic Equations
! Using LU decomposition and backsubstitution routines, it...
!   - finds the inverse of a matrix column by column
!
!

program test_matrix

    use lu_factorization, only: dp, crout, backsubstitution

    implicit none

    ! n - number of equations
    integer, parameter :: n = 3
    
    ! A - Input Matrix A(n, n)
    ! A_original - A matrix copy
    ! A_inv - A matrix inverted
    ! prod - Calculates the Identity matrix I = A . A^-1
    ! det - A matrix determinant
    real(dp) :: A(n, n), A_original(n, n), A_inv(n, n), prod(n, n), det
    
    ! indx - Interchanging rows register (record which row was changed by j row in j step)
    ! d - Interchanging parity or signal
    !     Used to calculate the det(A) = d ⋅ ∏U(i,i) - +1 if even number of changes, -1 if odd number
    ! info - 0 = ok; i > 0 = row i is zero (singular)
    integer  :: indx(n), d, info
    
    ! i, j - row, column indexes
    integer  :: i, j

    A = reshape([2.0_dp, 4.0_dp, 1.0_dp, &
                 1.0_dp, 3.0_dp, 2.0_dp, &
                 1.0_dp, 5.0_dp, 4.0_dp], [n, n])
    A_original = A

    ! Set up identity matrix
    A_inv = 0.0_dp
    do i = 1, n
        A_inv(i, i) = 1.0_dp
    end do

    call crout(A, indx, d, info)
    if (info /= 0) error stop 'singular matrix'

    do j = 1, n
        call backsubstitution(A, indx, A_inv(:, j))
    end do

    ! Exact zero test is intentional: 
    ! only a truly null row makes the matrix singular
    print '(a)', '(a) Inverse of A:'
    do i = 1, n
        print '(3f10.4)', A_inv(i, :)
    end do

    ! Check: A . A^-1 should be the identity
    prod = matmul(A_original, A_inv)

    print '(a)', 'A . A^-1 (should be identity):'
    do i = 1, n
        print '(3f10.4)', prod(i, :)
    end do


    det = real(d, dp)
    do j = 1, n
        det = det * A(j, j)
    end do

    print '(a, f10.4)', '(b) Determinant of A = ', det

end program test_matrix


! Run in the CMD: 
! gfortran -Wall -Wextra -fcheck=all -g slae02_LU_factorization.f90 slae03_LU_factorization_tests.f90 -o slae03_LU_factorization_tests

