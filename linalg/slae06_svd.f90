! SLAE - Solution of Linear Algebraic Equations
! SVD comprises the decomposition of a matrix A(m, n), where m >= n, into U, W and V, such that A = U * W * V^T
!
!
!

module svd

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private
    public :: dp, solve_backsubstitution, pythagorean, solve_decomposition

contains

    ! Solves A.x = B for a vector X, where A is specified by the arrays U, W and V
    subroutine solve_backsubstitution(u, w, v, m, n, mp, np, b, x)

        implicit none
        
        ! U - m x n column-orthogonal matrix (left singular vectors)
        ! W - n singular values (diagonal of the SVD)
        ! V - n x n orthogonal matrix (right singular vectors), with A = U * W * V^T
        real(dp), intent(in) :: u(mp, np), w(np), v(np, np)
        
        ! m, n - Logical dimensions of the matrix A
        ! mp, np - Physical dimensions of the matrix A
        integer, intent(in) :: m, n, mp, np

        ! input right hand side
        real(dp), intent(in) :: b(mp)
        ! output solution vector
        real(dp), intent(out) :: x(np)

        integer :: i, j, jj
        real(dp) :: s, tmp(np)

        ! Calculates U^t * B and divides by Wj
        do j = 1, n
            s = 0.0_dp
            if (w(j) /= 0.0_dp) then
                do i = 1, m
                    s = s + u(i, j) * b(i)
                end do
                s = s / w(j)
            end if
            tmp(j) = s
        end do

        ! Matrix multiply by V to get the solution vector X
        do j = 1, n
            s = 0.0_dp
            do jj = 1, n
                s = s + v(j, jj) * tmp(jj)
            end do
            x(j) = s
        end do
    
    end subroutine solve_backsubstitution

    function pythagorean(a, b) result(pythag)

        real(dp), intent(in) :: a, b
        real(dp) :: absa, absb, pythag

        absa = abs(a)
        absb = abs(b)

        if (absa >= absb) then
            pythag = absa * sqrt(1.0_dp + (absb / absa)**2)
        else
            if (absb == 0.0_dp) then
                pythag = 0.0_dp
            else
                pythag = absb * sqrt(1.0_dp + (absa / absb)**2)
            end if
        end if
    
    end function pythagorean

    ! Calculates the singular value decomposition of a matrix A(m, n), where m >= n, into U, W and V, such that A = U * W * V^T
    subroutine solve_decomposition(a, m, n, mp, np, w, v, info)

        implicit none

        integer, intent(in) :: m, mp, n, np

        integer, intent(out) :: info

        ! Inputed as A, Outputed as U
        ! (A is destroyed in the process)
        real(dp), intent(inout)  :: a(mp, np)
        ! W - singular values
        ! V - orthogonal right hand matrix 
        real(dp), intent(out) :: v(np, np), w(np)

        ! indexes
        integer :: i, its, j, jj, k, l, nm
        
        ! anorm - norm of the matrix A
        ! c, s - cos and sin of the rotation (Given's rotation)
        ! f, g, h - temporary variables (Householder / QR transformation)
        ! x, y, z - temporary variables (QR shift)
        ! tmp_sum - temporary sum for the Householder transformation
        ! scale - scaling factor for the Householder transformation
        ! rv1 - Superdiagonal elements of the bidiagonal matrix
        ! MAX_DIM - maximum dimension of the matrix A
        real(dp) :: anorm, c, s, f, g, h, tmp_sum, scale, x, y, z, rv1(n)

        info = 0
        g = 0.0_dp
        scale = 0.0_dp
        anorm = 0.0_dp

        ! Householder reduction to bidiagonal form
        do i = 1, n
            l = i + 1
            rv1(i) = scale * g
            g = 0.0_dp
            tmp_sum = 0.0_dp
            scale = 0.0_dp
            
            ! Householder at the left side
            if (i <= m) then
                do k = i, m
                    scale = scale + abs(a(k, i))
                end do
                if (scale /= 0.0_dp) then
                    do k = i, m
                        a(k, i) = a(k, i) / scale
                        tmp_sum = tmp_sum + a(k, i)**2
                    end do
                    f = a(i, i)
                    g = -sign(sqrt(tmp_sum), f)
                    h = f * g - tmp_sum
                    a(i, i) = f - g
                    do j = l, n
                        tmp_sum = 0.0_dp
                        do k = i, m
                            tmp_sum = tmp_sum + a(k, i) * a(k, j)
                        end do

                        if (abs(h) > tiny(1.0_dp)) then
                            f = tmp_sum / h
                        else
                            f = 0.0_dp
                        end if

                        do k = i, m
                            a(k, j) = a(k, j) + f * a(k, i)
                        end do
                    end do
                    do k = i, m
                        a(k, i) = a(k, i) * scale
                    end do
                end if
            end if
            w(i) = scale * g
            g = 0.0_dp
            tmp_sum = 0.0_dp
            scale = 0.0_dp

            ! Householder at the right side
            if (i <= m .and. i /= n) then
                do k = l, n
                    scale = scale + abs(a(i, k))
                end do
                if (scale /= 0.0_dp) then
                    do k = l, n
                        a(i, k) = a(i, k) / scale
                        tmp_sum = tmp_sum + a(i, k)**2
                    end do
                    f = a(i, l)
                    g = -sign(sqrt(tmp_sum), f)
                    h = f * g - tmp_sum
                    a(i, l) = f - g
                    do k = l, n
                        rv1(k) = a(i, k) / h
                    end do
                    do j = l, m
                        tmp_sum = 0.0_dp
                        do k = l, n
                            tmp_sum = tmp_sum + a(j, k) * a(i, k)
                        end do
                        do k = l, n
                            a(j, k) = a(j, k) + tmp_sum * rv1(k)
                        end do
                    end do
                    do k = l, n
                        a(i, k) = a(i, k) * scale
                    end do
                end if
            end if
            anorm = max(anorm, abs(w(i)) + abs(rv1(i)))
        end do


        ! V Accumulation (Right-hand transformations)
        g = 0.0_dp
        l = n + 1
        do i = n, 1, -1
            if (i < n) then
                if (g /= 0.0_dp) then
                    do j = l, n
                        v(j, i) = (a(i, j) / a(i, l)) / g
                    end do
                    do j = l, n
                        tmp_sum = 0.0_dp
                        do k = l, n
                            tmp_sum = tmp_sum + a(i, k) * v(k, j)
                        end do
                        do k = l, n
                            v(k, j) = v(k, j) + tmp_sum * v(k, i)
                        end do
                    end do
                end if
                do j = l, n
                    v(i, j) = 0.0_dp
                    v(j, i) = 0.0_dp
                end do
            end if
            v(i, i) = 1.0_dp
            g = rv1(i)
            l = i
        end do

        ! U Accumulation (Left-hand transformations)
        do i = min(m, n), 1, -1
            l = i + 1
            g = w(i)
            do j = l, n
                a(i, j) = 0.0_dp
            end do
            if (g /= 0.0_dp) then
                g = 1.0_dp / g
                do j = l, n
                    tmp_sum = 0.0_dp
                    do k = l, m
                        tmp_sum = tmp_sum + a(k, i) * a(k, j)
                    end do
                    
                    if (abs(a(i, i)) > tiny(1.0_dp)) then
                        f = (tmp_sum / a(i, i)) * g
                    else
                        f = 0.0_dp
                    end if
                    
                    do k = i, m
                        a(k, j) = a(k, j) + f * a(k, i)
                    end do
                end do
                do j = i, m
                    a(j, i) = a(j, i) * g
                end do
            else
                do j = i, m
                    a(j, i) = 0.0_dp
                end do
            end if
            a(i, i) = a(i, i) + 1.0_dp
        end do

        ! Diagonalization of the bidiagonal form (QR iteration with shifts)
        do k = n, 1, -1
            do its = 1, 30

                ! Test for splitting - find where to split the matrix
                do l = k, 1, -1
                    nm = l - 1
                    if ((abs(rv1(l)) + anorm) == anorm) exit
                    if ((abs(w(nm)) + anorm) == anorm) exit
                end do

                ! Given's rotation to eliminate rv1(l)
                c = 0.0_dp
                s = 1.0_dp
                do i = l, k
                    f = s * rv1(i)
                    rv1(i) = c * rv1(i)
                    if ((abs(f) + anorm) == anorm) exit
                    g = w(i)
                    h = pythagorean(f, g)
                    w(i) = h
                    h = 1.0_dp / h
                    c = g * h
                    s = -f * h
                    do j = 1, m
                        y = a(j, nm)
                        z = a(j, i)
                        a(j, nm) = (y * c) + (z * s)
                        a(j, i) = (z * c) - (y * s)
                    end do
                end do

                ! Convergence testing
                z = w(k)
                if (l == k) then
                    if (z < 0.0_dp) then    ! Singular value Force convergence >= 0
                        w(k) = -z
                        do j = 1, n
                            v(j, k) = -v(j, k)
                        end do
                    end if
                    exit
                end if
                if (its == 30) then
                    print *, "No convergence in 30 iterations"
                    info = k
                    return
                end if

                ! Wilkinson's shift and 2nd rotation
                x = w(l)
                nm = k - 1
                y = w(nm)
                g = rv1(nm)
                h = rv1(k)
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0_dp * h * y)
                g = pythagorean(f, 1.0_dp)

                if (abs(x) > tiny(1.0_dp)) then
                    f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x
                else
                    f = 0.0_dp
                end if

                ! Apply the rotation to V and A
                c = 1.0_dp
                s = 1.0_dp
                do j = l, nm
                    i = j + 1
                    g = rv1(i)
                    y = w(i)
                    h = s * g
                    g = c * g
                    z = pythagorean(f, h)
                    rv1(j) = z

                    if (z /= 0.0_dp) then
                        c = f / z
                        s = h / z
                    end if
                    
                        f = (x * c) + (g * s)
                    g = (g * c) - (x * s)
                    h = y * s
                    y = y * c
                    do jj = 1, n
                        x = v(jj, j)
                        z = v(jj, i)
                        v(jj, j) = (x * c) + (z * s)
                        v(jj, i) = (z * c) - (x * s)
                    end do
                    z = pythagorean(f, h)
                    w(j) = z
                    if (z /= 0.0_dp) then
                        z = 1.0_dp / z
                        c = f * z
                        s = h * z
                    end if
                    f = (c * g) + (s * y)
                    x = (c * y) - (s * g)
                    do jj = 1, m
                        y = a(jj, j)
                        z = a(jj, i)
                        a(jj, j) = (y * c) + (z * s)
                        a(jj, i) = (z * c) - (y * s)
                    end do
                end do
                rv1(l) = 0.0_dp
                rv1(k) = f
                w(k) = x
            end do
        end do
    
    end subroutine solve_decomposition

end module svd

