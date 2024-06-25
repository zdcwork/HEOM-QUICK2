subroutine calc_spectral_hb(itype, imode, zin, zout)
use matmod
use hbmod
!
! Important:  zin should be in unit of hbar 
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: imode, itype
complex*16, intent(in) :: zin
complex*16, intent(out) :: zout
integer :: ni
!
zout = czero
!
if (itype .eq. 1) then
   do ni=1,ndrude_hb
      zout = zout + zin * dcouple_hb(ni,imode) / dwidth_hb(ni) / (1.d0 + (zin/dwidth_hb(ni))**2)
   end do
else if (itype .eq. 2) then
   do ni=1,ndrude_hb
      zout = zout + zin * dcouple_hb(ni,imode) / dwidth_hb(ni) / (1.d0 + (zin/dwidth_hb(ni))**2)**2
   end do
else if (itype .eq. 3) then
   do ni=1,ndrude_hb
      zout = zout + zin * dcouple_hb(ni,imode) / dwidth_hb(ni) *           &
             ( 1.d0 / (1.d0 + ((zin - dcenter_hb(ni))/dwidth_hb(ni))**2)   &
             + 1.d0 / (1.d0 + ((zin + dcenter_hb(ni))/dwidth_hb(ni))**2) )
   end do
else 
   write(6,*)'calc_spectral_hb: error! unknown itype ', itype
   stop
end if
!
return
end subroutine calc_spectral_hb
