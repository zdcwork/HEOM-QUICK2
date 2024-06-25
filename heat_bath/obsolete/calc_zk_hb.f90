subroutine calc_zk_hb
use matmod
use tmpmatmod
use sparsemod
!
! purpose: compute K and \tilde{K} terms as initial values for 
!          time propagation 
! 
implicit none
include '../include/sizes'
include '../include/common'
!
integer                   :: istat, nwork, info, imode, iorbs, ispin, icor, idrude
integer                   :: ni, nj, ntmp1
integer                   :: nrho2
integer*8                 :: lni
logical                   :: lherm
real*8                    :: dtmp1, dtmp2
real*8,     allocatable   :: rwork(:), dval0(:)
complex*16, allocatable   :: zwork(:,:)
complex*16, allocatable   :: cmat0(:,:), cmat1(:,:)
complex*16, allocatable   :: zmat0(:,:), zmat1(:,:)
!
!write(6,*)
!write(6,*)'calc_zk_hb: entering subroutine'
!call flush(6)
!
if (.not. lhb) then
   write(6,*)
   write(6,*)'calc_zk_hb: error entrance! '
   stop
end if
!
! Q_a = sum_s a_s^+ a_s
!
do imode=1,nmode_hb
   iorbs = imode
   cmtmp1(1:nrho,1:nrho) = czero
   do ispin=1,nspin
      cmtmp2(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,iorbs,ispin) * cunity
      call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, cunity, cmtmp1, nrho)
   end do
   qa_hb(1:nrho,1:nrho,imode) = cmtmp1(1:nrho,1:nrho)
end do
!
allocate(zk_hb(nrho,nrho,ncor_hb,nmode_hb), STAT=istat)
allocate(zk0_hb(nrho,nrho,ncor_hb,nmode_hb), STAT=istat)
allocate(zk_tmp_hb(nrho,nrho,ncor_hb,nmode_hb), STAT=istat)
allocate(zk_rhs_hb(nrho,nrho,ncor_hb,nmode_hb), STAT=istat)
allocate(zq_hb(nrho,nrho,nmode_hb), STAT=istat)
if (itype_hb .eq. 2) then
   allocate(ztk_hb(nrho,nrho,ndrude_hb,nmode_hb), STAT=istat)
   allocate(ztk0_hb(nrho,nrho,ndrude_hb,nmode_hb), STAT=istat)
   allocate(ztk_tmp_hb(nrho,nrho,ndrude_hb,nmode_hb), STAT=istat)
   allocate(ztk_rhs_hb(nrho,nrho,ndrude_hb,nmode_hb), STAT=istat)
end if
zq_hb(1:nrho,1:nrho,1:nmode_hb) = czero
!
! diagonalize system Liouvillian at initial ground state (dhs is not considered)
!
allocate(cmat0(nrho,nrho), cmat1(nrho,nrho), STAT=istat)
if (lsparse) then
   lni = ind_spa(1)
   call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rho_spa(lni), &
                     cmat0, nrho, nrho, nrho)
else 
   cmat0(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
end if
call calchs(cmat0(1,1))
cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
!
nrho2 = nrho**2
nwork = nrho2**2
allocate(zmat0(nrho2,nrho2), zmat1(nrho2,nrho2), dval0(nrho2), STAT=istat)
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), STAT=istat)
!
call zhil2liou('l', 'n', nrho, cmat1, nrho, zmat0, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zmat1, nrho2)
zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
call checkhermicity(zmat0, nrho2, zmat1, nrho2, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)'calc_zk_hb: error! Liouville matrix is not Hermitian '
   stop
end if
!
call zheev('V', 'u', nrho2, zmat0, nrho2, dval0, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'calc_zk_hb: error! diagonalization failed 1 ', info
   stop
end if
!
do imode=1,nmode_hb
   cmtmp1(1:nrho,1:nrho) = qa_hb(1:nrho,1:nrho,imode)
! X^(-1) * Qa
   call zgemv('c', nrho2, nrho2, cunity, zmat0, nrho2, qa_hb(1,1,imode), 1, czero, cmtmp3, 1)
!
   do icor=1,ncor_hb
       ntmp1 = 0
       do nj=1,nrho
          do ni=1,nrho
             ntmp1 = ntmp1 + 1
             cmtmpa(ni,nj) = cmtmp3(ni,nj) * (-1.d0) / ( cgamma_hb(icor,imode) - eye * dval0(ntmp1) )
          end do
       end do
       call zgemv('n', nrho2, nrho2, cunity, zmat0, nrho2, cmtmpa, 1, czero, cmtmp5, 1)
       zk_hb(1:nrho,1:nrho,icor,imode) = cmtmp5(1:nrho,1:nrho)
       zq_hb(1:nrho,1:nrho,imode) = zq_hb(1:nrho,1:nrho,imode) + cb_hb(icor,imode) * cmtmp5(1:nrho,1:nrho)
   end do
!
   if (itype_hb .eq. 2) then
      do idrude=1,ndrude_hb
         icor = idrude + npade_hb
         ntmp1 = 0
         do nj=1,nrho
            do ni=1,nrho
               ntmp1 = ntmp1 + 1
               cmtmpa(ni,nj) = cmtmp3(ni,nj) / ( cgamma_hb(icor,imode) - eye * dval0(ntmp1) )**2
            end do
         end do
         call zgemv('n', nrho2, nrho2, cunity, zmat0, nrho2, cmtmpa, 1, czero, cmtmp5, 1)
         ztk_hb(1:nrho,1:nrho,idrude,imode) = cmtmp5(1:nrho,1:nrho)
!
         zq_hb(1:nrho,1:nrho,imode) = zq_hb(1:nrho,1:nrho,imode) + cd_hb(idrude,imode) * cmtmp5(1:nrho,1:nrho)
      end do
   end if
end do
!
deallocate(cmat0, cmat1, STAT=istat)
deallocate(zmat0, zmat1, dval0, zwork, rwork, STAT=istat)
!
!write(6,*)
!write(6,*)'calc_zk_hb: leaving subroutine'
!call flush(6)
!
return
end subroutine calc_zk_hb
