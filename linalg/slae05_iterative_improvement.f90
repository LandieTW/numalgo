! SLAE - Solution of Linear Algebraic Equations
!
!
!
!

module iteractive_methods

    use lu_factorization, only: backsubstitution
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    private
    public :: dp, improve_solution

    integer,  parameter :: dp = real64
    integer, parameter :: MAX_DIM = 500     ! limit system size

contains

    ! Improove iteractively a solution X for A.x = b, considering iteractive iteractive refinement
    subroutine improve_solution(matrix, lu_matrix, n, pivots, rhs, solution)

        integer, intent(in)     :: n
        real(dp), intent(in)    :: matrix(:, :)     ! Original Matrix A(n, n)
        real(dp), intent(in)    :: lu_matrix(:, :)  ! A decomposed in LU
        integer, intent(in)     :: pivots(:)        ! Pivoting indexes Vector
        real(dp), intent(in)    :: rhs(:)           ! Vector b of right side
        real(dp), intent(inout) :: solution(:)      ! aproximated solution

        integer         :: i, j
        real(dp)        :: residual(MAX_DIM)
        real(dp)      :: residual_acc         ! double precision accumulater

        ! Calculates the residual value r = b - A.x (with double precision)
        do i = 1, n
            residual_acc = -dble(rhs(i))
            do j = 1, n
                residual_acc = residual_acc + dble(matrix(i, j)) * dble(solution(j))
            end do
            residual(i) = real(residual_acc)
        end do

        ! Solves LU . δ = r to obtain the error term
        call backsubstitution(lu_matrix, pivots, residual)

        ! Updates solution x ← x - δ
        do i = 1, n
            solution(i) = solution(i) - residual(i)
        end do
    
    end subroutine improve_solution

end module iteractive_methods


! Run in the CMD: 
! gfortran -Wall -Wextra -fcheck=all -g slae02_LU_factorization.f90 slae05_iterative_improvement.f90 -o slae05_iterative_improvement

