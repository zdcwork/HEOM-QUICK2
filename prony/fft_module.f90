module fft_module
  contains
    subroutine fftfreq(N, d, freq)
      integer, intent(in) :: N
      real*8, intent(in), optional :: d
      real*8, intent(inout) :: freq(:)
      real*8 :: delta
      integer :: i, cutoff
    
      ! 处理可选参数
      delta = 1.d0
      if (present(d)) delta = d
    
      ! 计算频率分界点
      if (mod(N,2) == 0) then
        cutoff = N/2 - 1
      else
        cutoff = (N-1)/2
      end if
    
      ! 生成频率序列
      do i = 1, N
        if (i-1 <= cutoff) then
          freq(i) = dble(i-1)
        else
          freq(i) = dble(i-1 - N)
        end if
      end do
    
      ! 归一化处理
      freq = freq / (delta*N)
    end subroutine fftfreq
    !
    subroutine fftshift_1d(tt)
      implicit none
      class(*), intent(inout) :: tt(:)
      complex*16, allocatable :: cshifted_tt(:)
      real*8, allocatable :: dshifted_tt(:)
      integer :: n, half
    
      n = size(tt)
      half = n / 2
      !
      select type (tt)
        type is (real*8)
          print *, "Array is of type real*8"
          allocate(dshifted_tt(n))
          dshifted_tt = 0.d0
          ! Swap the first and second halves of the array
          dshifted_tt(1:half) = tt(half+1:n)
          dshifted_tt(half+1:n) = tt(1:half)
          tt = dshifted_tt
        type is (complex*16)
          print *, "Array is of type complex*16"
          allocate(cshifted_tt(n))
          cshifted_tt = 0.d0
          ! Swap the first and second halves of the array
          cshifted_tt(1:half) = tt(half+1:n)
          cshifted_tt(half+1:n) = tt(1:half)
          tt = cshifted_tt
         type is (complex*8)
           print *, "Array is of type complex*8"
           allocate(cshifted_tt(n))
           cshifted_tt = 0.d0
           ! Swap the first and second halves of the array
           cshifted_tt(1:half) = tt(half+1:n)
           cshifted_tt(half+1:n) = tt(1:half)
           tt = cshifted_tt
        class default
          print *, "Array is of an unknown or unsupported type"
          stop
      end select
    end subroutine fftshift_1d
  !
end module fft_module
