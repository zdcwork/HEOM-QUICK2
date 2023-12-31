subroutine bosefunc(npoles, zcoef, zpole, zin, zout)
implicit none
!
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
   write(6,*)'bosefunc: error! input energy cannot be zero ', zin
   stop
end if
!
zout = 5.d-1 + 1.d0 / zin
!
do ni=1,npoles
   zout = zout + zcoef(ni) * (1.d0 / (zin - zpole(ni)) + 1.d0 / (zin + zpole(ni)))
end do
!
return
!
end subroutine bosefunc
