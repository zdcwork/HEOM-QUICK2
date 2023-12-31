real*8 function dnormlong(len0, din)
!
! zin^+ * zin (for very long arrays)
!
integer*8, intent(in) :: len0
real*8, intent(in) :: din(*)
!
integer*8 :: lni
real*8 dtemp
!
dtemp = 0.d0
do lni=1,len0
 dtemp = dtemp + din(lni)**2
end do
dnormlong = dtemp
end function dnormlong
