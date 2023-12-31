subroutine refresh_zq_hb
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: imode, icor, idrude
integer :: ni, nj
real*8  :: dtmp1, dtmp2
!
zq_hb(1:nrho,1:nrho,1:nmode_hb) = czero
!
do imode=1,nmode_hb
   do icor=1,ncor_hb
      zq_hb(1:nrho,1:nrho,imode) = zq_hb(1:nrho,1:nrho,imode) + cb_hb(icor,imode) * &
                                   zk_tmp_hb(1:nrho,1:nrho,icor,imode)
   end do
   if (itype_hb .eq. 2) then
      do idrude=1,ndrude_hb
         zq_hb(1:nrho,1:nrho,imode) = zq_hb(1:nrho,1:nrho,imode) + cd_hb(idrude,imode) * &
                                      ztk_tmp_hb(1:nrho,1:nrho,idrude,imode)
      end do
   end if
end do
!
return
end subroutine refresh_zq_hb
