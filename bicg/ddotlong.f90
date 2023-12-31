real*8 function ddotlong(len0, dx, dy)
!
! zx^t * zy (for very long arrays)
!
integer*8, intent(in) :: len0
real*8,    intent(in) :: dx(*), dy(*)
!
integer*8 :: lni
real*8    :: dtmp1
!
dtmp1 = 0.d0
do lni=1,len0
 dtmp1 = dtmp1 + dx(lni) * dy(lni)
end do
ddotlong = dtmp1
end function ddotlong
!
real*8 function cdotlong(len0, cx, cy)
!
integer*8,  intent(in) :: len0
complex*16, intent(in) :: cx(*), cy(*)
!
integer*8 :: lni
real*8    :: dtmp1
!
dtmp1 = 0.d0
do lni=1,len0
  dtmp1 = dtmp1 + dble(cx(lni))*dble(cy(lni)) + dimag(cx(lni))*dimag(cy(lni)) 
end do
cdotlong = dtmp1
end function cdotlong

