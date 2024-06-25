subroutine calc_deriv_zk_hb(tt)
!
! purpose: calculate the RHS of EOM for K-terms
!
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,  intent(in) :: tt
integer             :: ni, nj, imode, icor, idrude, iorbs1, ispin
integer*8           :: lni
real*8              :: dtmp0, dtmp1, dtmp2, dtmp3, dtmp4
integer             :: itimes
data itimes /0/
!
if (itimes .eq. 0) then
   write(6,*)'calc_deriv_zk_hb: first entry '
end if
!
if (lsparse) then
   lni = ind_spa(1)
   call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rhotmp_spa(lni), cmtmp1, nrho, nrho, nrho)
   call calchs(cmtmp1)
   call calcdhs(tt, cmtmp1)
else
   call calchs(rhotmp(1,1,1))
   call calcdhs(tt, rhotmp(1,1,1))
end if
cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!
do imode=1,nmode_hb
   cmtmp1(1:nrho,1:nrho) = qa_hb(1:nrho,1:nrho,imode)
!
   do icor=1,ncor_hb
      zk_rhs_hb(1:nrho,1:nrho,icor,imode) = cmtmp1(1:nrho,1:nrho) + cgamma_hb(icor,imode) * zk_tmp_hb(1:nrho,1:nrho,icor,imode)
      call zgemm('n', 'n', nrho, nrho, nrho, -eye, cmtmp3, nrho, zk_tmp_hb(1,1,icor,imode), nrho, cunity, zk_rhs_hb(1,1,icor,imode), nrho)
      call zgemm('n', 'n', nrho, nrho, nrho,  eye, zk_tmp_hb(1,1,icor,imode), nrho, cmtmp3, nrho, cunity, zk_rhs_hb(1,1,icor,imode), nrho)
   end do
!
   if (itype_hb .eq. 2) then
      do idrude=1,ndrude_hb
         icor = idrude + npade_hb
         ztk_rhs_hb(1:nrho,1:nrho,idrude,imode) = zk_tmp_hb(1:nrho,1:nrho,icor,imode) + cgamma_hb(icor,imode) * &
                                                  ztk_tmp_hb(1:nrho,1:nrho,idrude,imode)
         call zgemm('n', 'n', nrho, nrho, nrho, -eye, cmtmp3, nrho, ztk_tmp_hb(1,1,idrude,imode), nrho, &
                    cunity, ztk_rhs_hb(1,1,idrude,imode), nrho)
         call zgemm('n', 'n', nrho, nrho, nrho,  eye, ztk_tmp_hb(1,1,idrude,imode), nrho, cmtmp3, nrho, &
                    cunity, ztk_rhs_hb(1,1,idrude,imode), nrho)
      end do
   end if
end do
!
itimes = itimes + 1
return
end subroutine calc_deriv_zk_hb
