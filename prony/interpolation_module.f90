module cubic_spline_interp
    implicit none
    private
    public :: cspline_interp_complex
!
contains
    ! 主插值子程序
    subroutine cspline_interp_complex(x_old, y_old, x_new, y_new)
        implicit none
        real*8, intent(in) :: x_old(:)        ! 原始横坐标（必须单调递增）
        complex*16, intent(in) :: y_old(:)     ! 原始复数数组
        real*8, intent(in) :: x_new(:)        ! 新横坐标
        complex*16, intent(out) :: y_new(:)    ! 插值结果
        complex*16, allocatable :: y2(:)       ! 二阶导数数组
        integer :: n_old, n_new, i

        n_old = size(x_old)
        n_new = size(x_new)

        ! 检查输入合法性
        if (size(y_old) /= n_old .or. size(y_new) /= n_new) then
            error stop "Error: Input array size mismatch"
        end if

        ! 计算三次样条的二阶导数
        allocate(y2(n_old))
        call calc_spline_coeff(x_old, y_old, y2)

        ! 对每个新点进行插值
        do i = 1, n_new
            y_new(i) = spline_eval(x_old, y_old, y2, x_new(i))
        end do

        deallocate(y2)

    contains
        ! 计算三次样条的二阶导数 y''
        subroutine calc_spline_coeff(x, y, y2)
            implicit none
            real*8, intent(in) :: x(:)
            complex*16, intent(in) :: y(:)
            complex*16, intent(out) :: y2(:)
            complex*16, allocatable :: a(:), b(:), c(:), r(:)
            integer :: n, i

            n = size(x)
            allocate(a(n), b(n), c(n), r(n))

            ! 构建三对角矩阵
            do i = 2, n-1
                a(i) = (x(i) - x(i-1)) / 6.0d0
                b(i) = (x(i+1) - x(i-1)) / 3.0d0
                c(i) = (x(i+1) - x(i)) / 6.0d0
                r(i) = ((y(i+1) - y(i)) / (x(i+1) - x(i))) - &
                       ((y(i) - y(i-1)) / (x(i) - x(i-1)))
            end do

            ! 自然边界条件 (y''=0 在两端)
            a(1) = 0.0d0
            b(1) = 1.0d0
            c(1) = 0.0d0
            r(1) = 0.0d0
                      
            a(n) = 0.0d0
            b(n) = 1.0d0
            c(n) = 0.0d0
            r(n) = 0.0d0

            ! 解三对角方程组（托马斯算法）
            call tridiag_solve(a, b, c, r, y2)

            deallocate(a, b, c, r)
        end subroutine

        ! 托马斯算法解三对角方程组
        subroutine tridiag_solve(a, b, c, r, y)
            implicit none
            complex*16, intent(in) :: a(:), b(:), c(:), r(:)
            complex*16, intent(out) :: y(:)
            complex*16, allocatable :: gam(:)
            complex*16 :: bet
            integer :: n, j

            n = size(a)
            allocate(gam(n))

            bet = b(1)
            y(1) = r(1) / bet

            do j = 2, n
                gam(j) = c(j-1) / bet
                bet = b(j) - a(j) * gam(j)
                y(j) = (r(j) - a(j) * y(j-1)) / bet
            end do

            do j = n-1, 1, -1
                y(j) = y(j) - gam(j+1) * y(j+1)
            end do

            deallocate(gam)
        end subroutine

        ! 计算插值点处的值
        complex*16 function spline_eval(x, y, y2, x_val)
            real*8, intent(in) :: x(:), x_val
            complex*16, intent(in) :: y(:), y2(:)
            integer :: k, khi, klo
            real*8 :: h, a, b

            ! 查找区间
            klo = 1
            khi = size(x)
            do while (khi - klo > 1)
                k = (khi + klo) / 2
                if (x(k) > x_val) then
                    khi = k
                else
                    klo = k
                end if
            end do

            h = x(khi) - x(klo)
            a = (x(khi) - x_val) / h
            b = (x_val - x(klo)) / h

            ! 三次样条公式
            spline_eval = a * y(klo) + b * y(khi) + &
                ((a**3 - a) * y2(klo) + (b**3 - b) *&
                 y2(khi)) * (h**2) / 6.0d0
        end function
    end subroutine cspline_interp_complex
end module cubic_spline_interp
