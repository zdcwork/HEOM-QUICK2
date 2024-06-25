real*8 function ddotlong(len0, dx, dy)
!
! zx^+ * zy (for very long arrays)
!
integer*8, intent(in) :: len0
real*8, intent(in) :: dx(*), dy(*)
!
integer*8 :: lni
real*8 :: dtmp1
!
dtmp1 = 0.d0
do lni=1,len0
 dtmp1 = dtmp1 + dx(lni) * dy(lni)
end do
ddotlong = dtmp1
end function ddotlong
