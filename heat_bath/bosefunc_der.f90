subroutine bosefunc_der(npoles, zcoef, zpole, zin, zout)
implicit none
!
! evalute the derivative of
! 1/(1-e^(-z)) = 1/2 + 1/z + sum_j zcoef_j * [1/(z-zpole_j) + 1/(z+zpole_j)]
!
! zin and zout are dimensionless numbers
!
integer,    intent(in)  :: npoles
complex*16, intent(in)  :: zin, zcoef(*), zpole(*)
complex*16, intent(out) :: zout
!
integer :: ni
!
if (cdabs(zin) .le. 1.d-15) then
   write(6,*)
   write(6,*)'bosefunc_der: error! input energy cannot be zero ', zin
   stop
end if
!
zout = -1.d0 / zin**2
!
do ni=1,npoles
   zout = zout - zcoef(ni) * (1.d0 / (zin - zpole(ni))**2 + 1.d0 / (zin + zpole(ni))**2)
end do
!
return
!
end subroutine bosefunc_der
