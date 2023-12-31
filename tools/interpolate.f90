program interpolate
implicit none

! x1 value1
! x2 value2
! x = ?  given value
!
real*8 :: x1, x2, x, v1, v2, value, dtmp1
!
read(*,*)x1, v1
read(*,*)x2, v2
read(*,*)value
!
!  
! (x - x1) / (value - v1) = (x2 - x1) / (v2 - v1)
!
x = (x2 - x1) / (v2 - v1) * (value - v1) + x1
!
write(*,*)' x = ', x


end program interpolate
