complex*16 function zdotclong(len0, zx, zy)
!
! zx^+ * zy (for very long arrays)
!
integer*8, intent(in) :: len0
complex*16, intent(in) :: zx(*), zy(*)
!
integer*8 :: lni
complex*16 ztemp
!
ztemp = dcmplx(0.d0, 0.d0)
do lni=1,len0
 ztemp = ztemp + dconjg(zx(lni)) * zy(lni)
end do
zdotclong = ztemp
end function zdotclong
