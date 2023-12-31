real*8 function dznrm2long(len0, zin)
!
! zin^+ * zin (for very long arrays)
!
integer*8, intent(in) :: len0
complex*16, intent(in) :: zin(*)
!
integer*8 :: lni
real*8 dtemp
!
dtemp = 0.d0
do lni=1,len0
 dtemp = dtemp + dble(dconjg(zin(lni)) * zin(lni))
end do
dznrm2long = dtemp
end function dznrm2long
