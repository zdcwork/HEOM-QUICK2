subroutine calc_tildeq_hb(ldim,imode,rcgtmp_hb,zout) 
use matmod
use hbmod
!
! Compute the 'ldim'-th element of \tilde{Q} vector
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8,  intent(in)  :: ldim
integer,    intent(in)  :: imode
complex*16, intent(in)  :: rcgtmp_hb(nrho,nrho,*) ! last dimension is nunk*(nbath_hb+1)
complex*16, intent(out) :: zout(nrho,*)
!
integer                 :: ni, ibath, nhalf, inkp, icor
complex*16              :: ctmp1, ctmp2
!
zout(1:nrho,1:nrho) = czero
ibath = (imode - 1) * nkp_hb 
nhalf = nmode_hb * nkp_hb
!
do ni=1,nn_hb
   icor  = ni
   ctmp1 = 2.d0 * eye * dimag(cb_hb(icor,imode))
   inkp  = ibath + ni
   zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + ctmp1 * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
   ctmp1 = 2.d0 * eye * dble(cb_hb(icor,imode)) 
   inkp  = nhalf + ibath + ni
   zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + ctmp1 * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
end do
!
ibath = ibath + nn_hb
if (nk_hb .gt. 0) then  ! only involves itype_hb=3
   icor = npade_hb
   inkp = ibath 
   do ni=1,ndrude_hb*2
      icor = icor + 1
      inkp = inkp + 1
!      cmtmp1(1:nrho,1:nrho) = rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk) + eye * rcgtmp_hb(1:nrho,1:nrho,ldim+(inkp+nhalf)*nunk)
!      cmtmp2(1:nrho,1:nrho) = rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk) - eye * rcgtmp_hb(1:nrho,1:nrho,ldim+(inkp+nhalf)*nunk)
!      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + cb_hb(icor,imode)  * cmtmp1(1:nrho,1:nrho) - &
!                                           dconjg(cb_hb(icor,imode)) * cmtmp2(1:nrho,1:nrho)
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) +              cb_hb(icor,imode)  * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) +        eye * cb_hb(icor,imode)  * rcgtmp_hb(1:nrho,1:nrho,ldim+(inkp+nhalf)*nunk)
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) -       dconjg(cb_hb(icor,imode)) * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + eye * dconjg(cb_hb(icor,imode)) * rcgtmp_hb(1:nrho,1:nrho,ldim+(inkp+nhalf)*nunk)
   end do
end if
!
ibath = ibath + 2 * nk_hb
if (np_hb .gt. 0) then  ! only involves itype_hb=2
   do ni=1,ndrude_hb
      icor = ni
      inkp = ibath + ni
      ctmp1 = 2.d0 * eye * dimag(cd_hb(icor,imode))
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + ctmp1 * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
      ctmp1 = 2.d0 * eye * dble(cd_hb(icor,imode))
      inkp = nhalf + ibath + ni
      zout(1:nrho,1:nrho) = zout(1:nrho,1:nrho) + ctmp1 * rcgtmp_hb(1:nrho,1:nrho,ldim+inkp*nunk)
   end do
end if
!
return
!
end subroutine calc_tildeq_hb
