subroutine calc_qst_hb
use matmod
use hbmod  
use tmpmatmod
use sparsemod
!
! purpose: compute \tilde{Q}_{a} at stationary state
! 
implicit none
include '../include/sizes'
include '../include/common'
!
integer                   :: istat, nwork, info, imode, iorbs, ispin, icor, idrude
integer                   :: ni, nj, ntmp1
integer                   :: isgn
integer                   :: nrho2
integer*8                 :: lni
logical                   :: lherm
real*8                    :: dtmp1
integer,    allocatable   :: ipiv(:)
real*8,     allocatable   :: rwork(:)
complex*16, allocatable   :: zwork(:,:), zval0(:), tmpval(:)
complex*16, allocatable   :: cmat0(:,:), cmat1(:,:)
complex*16, allocatable   :: zmat0(:,:), zmat1(:,:), zmat2(:,:), zmat3(:,:), zmat4(:,:)
complex*16, allocatable   :: zmat9(:,:)
!
!write(6,*)
!write(6,*)'calc_qst_hb: entering subroutine'
!call flush(6)
!
if (.not. lhb) then
   write(6,*)
   write(6,*)'calc_qst_hb: error entrance! '
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
! diagonalize system Liouvillian
!
nrho2 = nrho**2
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
if (igroundsteady .ne. 0) then
   call calcdhs(1.d20, cmat0(1,1))
   cmat1(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
end if
!
nwork = nrho2**2
allocate(ipiv(nrho2), STAT=istat)
allocate(zmat0(nrho2,nrho2), zmat1(nrho2,nrho2), zval0(nrho2), STAT=istat)
allocate(zmat2(nrho2,nrho2), zmat3(nrho2,nrho2), zmat4(nrho2,nrho2), STAT=istat)
allocate(zmat9(nrho2,nrho2), STAT=istat)
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
allocate(zxmat_hb(nrho2,nrho2), dleig_hb(nrho2), STAT=istat)
allocate(zdiag_hb(nrho2,nmode_hb), zdiagp_hb(nrho2,nmode_hb), STAT=istat)
allocate(zdiagc_hb(nrho2,nmode_hb), zdiagpc_hb(nrho2,nmode_hb), STAT=istat)
allocate(zcoef_hb(nmode_hb), zcoefp_hb(nmode_hb), STAT=istat)
allocate(zqst_hb(nrho,nrho,nmode_hb), zqpst_hb(nrho,nrho,nmode_hb), STAT=istat)
allocate(zqcst_hb(nrho,nrho,nmode_hb), zqpcst_hb(nrho,nrho,nmode_hb), STAT=istat)
!
zqst_hb = czero
zqpst_hb = czero
zqcst_hb = czero
zqpcst_hb = czero
!
call zhil2liou('l', 'n', nrho, cmat1, nrho, zmat1, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zmat2, nrho2)
zmat0(1:nrho2,1:nrho2) = -eye * (zmat1(1:nrho2,1:nrho2) - zmat2(1:nrho2,1:nrho2))
zmat9(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2)
!
cmtmp3(1:nrho,1:nrho) = czero  ! Left
cmtmp4(1:nrho,1:nrho) = czero  ! Right
do isgn=1,nsgn
   do ispin=1,nspin
      do iorbs=1,norbs
         cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs,ispin), 0.d0)
         do nj=1,nrho
            do ni=1,nrho
               cmtmp2(ni,nj) = dconjg(cmtmp1(nj,ni))
               cmtmp5(ni,nj) = dconjg(zast_hb(nj,ni,iorbs,ispin,isgn))
            end do
         end do
         if (isgn .eq. 1) then 
            cmtmpa(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho)  ! a^\sigma_{\mu s}
            cmtmpb(1:nrho,1:nrho) = cmtmp1(1:nrho,1:nrho)  ! a^{\bar\sigma}_{\mu s}
         else
            cmtmpa(1:nrho,1:nrho) = cmtmp1(1:nrho,1:nrho)  ! a^\sigma_{\mu s}
            cmtmpb(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho)  ! a^{\bar\sigma}_{\mu s}
         end if
         call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmpa, nrho, zast_hb(1,1,iorbs,ispin,isgn), nrho, cunity, cmtmp3, nrho)
         call zgemm('c', 'n', nrho, nrho, nrho, cunity, zast_hb(1,1,iorbs,ispin,isgn), nrho, cmtmpa, nrho, cunity, cmtmp4, nrho)
!
         call zh2lmlr(nrho, zast_hb(1,1,iorbs,ispin,isgn), nrho, cmtmpb, nrho, zmat1, nrho2)         
         zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) + zmat1(1:nrho2,1:nrho2)
         zmat9(1:nrho2,1:nrho2) = zmat9(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
!
         call zh2lmlr(nrho, cmtmpa, nrho, cmtmp5, nrho, zmat1, nrho2)
         zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) + zmat1(1:nrho2,1:nrho2)
         zmat9(1:nrho2,1:nrho2) = zmat9(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
      end do
   end do
end do
call zhil2liou('l', 'n', nrho, cmtmp3, nrho, zmat1, nrho2)
zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
zmat9(1:nrho2,1:nrho2) = zmat9(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
call zhil2liou('r', 'n', nrho, cmtmp4, nrho, zmat1, nrho2)
zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
zmat9(1:nrho2,1:nrho2) = zmat9(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
!

zmat9(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2)


!
! diagonalize system Liouvillian
!
zmat3(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2)
call zgeev('n', 'v', nrho2, zmat0, nrho2, zval0, zmat1, nrho2, zmat2, nrho2, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error! diagonalization failed 1 ', info
   stop
end if
!
zmat1(1:nrho2,1:nrho2) = zmat2(1:nrho2,1:nrho2)
call zgetrf(nrho2, nrho2, zmat1, nrho2, ipiv, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error message for zgetrf 1 ', info
   stop
end if
call zgetri(nrho2, zmat1, nrho2, ipiv, zwork(1,1), nwork, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error message for zgetri 1 ', info
   stop
end if
!
zmat0(1:nrho2,1:nrho2) = zmat3(1:nrho2,1:nrho2)
do nj=1,nrho2
   do ni=1,nrho2
      zmat3(ni,nj) = zmat1(ni,nj) * zval0(ni)
   end do
end do
zmat4(1:nrho2,1:nrho2) = -zmat0(1:nrho2,1:nrho2)
call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zmat2, nrho2, zmat3, nrho2, cunity, zmat4, nrho2)
!
call cmaxmat(nrho2, nrho2, zmat4, nrho2, dtmp1)
write(6,*)
write(6,*)'calc_qst_hb: check-1 |A-X*Ad*X^(-1)|, max-diff ', dtmp1
call flush(6)
!
! zmat0 --> A
! zmat1 --> X^(-1)
! zmat2 --> X
!
nleig_hb = nrho2
!
zdiag_hb(1:nrho2,1:nmode_hb) = czero
zdiagp_hb(1:nrho2,1:nmode_hb) = czero
zdiagc_hb(1:nrho2,1:nmode_hb) = czero
zdiagpc_hb(1:nrho2,1:nmode_hb) = czero
!
do imode=1,nmode_hb
   do icor=1,ncor_hb
      do ni=1,nleig_hb
         tmpval(ni) = -1.d0 / (cgamma_hb(icor,imode) + zval0(ni)) 
         zdiag_hb(ni,imode) = zdiag_hb(ni,imode) + cb_hb(icor,imode) * tmpval(ni)
      end do
      do ni=1,nleig_hb
         tmpval(ni) = -1.d0 / (dconjg(cgamma_hb(icor,imode)) + zval0(ni)) 
         zdiagc_hb(ni,imode) = zdiagc_hb(ni,imode) + dconjg(cb_hb(icor,imode)) * tmpval(ni)
      end do
   end do
   if (itype_hb .eq. 2) then
      do idrude=1,ndrude_hb
         icor = idrude + npade_hb
         do ni=1,nleig_hb
            tmpval(ni) = 1.d0 / (cgamma_hb(icor,imode) + zval0(ni))**2
            zdiagp_hb(ni,imode) = zdiagp_hb(ni,imode) + cd_hb(idrude,imode) * tmpval(ni)
         end do
         do ni=1,nleig_hb
            tmpval(ni) = 1.d0 / (dconjg(cgamma_hb(icor,imode)) + zval0(ni))**2
            zdiagpc_hb(ni,imode) = zdiagpc_hb(ni,imode) + dconjg(cd_hb(idrude,imode)) * tmpval(ni)
         end do
      end do
   end if
end do
!
do imode=1,nmode_hb
   call zgemv('n', nrho2, nrho2, cunity, zmat1, nrho2, qa_hb(1,1,imode), 1, czero, cmtmp1, 1)
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         ntmp1 = ntmp1 + 1
         cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiag_hb(ntmp1,imode)
      end do
   end do
   call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, czero, zqst_hb(1,1,imode), 1)
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         ntmp1 = ntmp1 + 1
         cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagc_hb(ntmp1,imode)
      end do
   end do
   call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, czero, zqcst_hb(1,1,imode), 1)
   if (itype_hb .eq. 2) then
      call zgemv('n', nrho2, nrho2, cunity, zmat1, nrho2, qa_hb(1,1,imode), 1, czero, cmtmp1, 1)
      ntmp1 = 0
      do nj=1,nrho
         do ni=1,nrho
            ntmp1 = ntmp1 + 1
            cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagp_hb(ntmp1,imode)
         end do
      end do
      call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, cunity, zqst_hb(1,1,imode), 1)
      ntmp1 = 0
      do nj=1,nrho
         do ni=1,nrho
            ntmp1 = ntmp1 + 1
            cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagpc_hb(ntmp1,imode)
         end do
      end do
      call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, cunity, zqcst_hb(1,1,imode), 1)
   end if
end do
!


!
! diagonalize system Liouvillian
!
zmat3(1:nrho2,1:nrho2) = zmat9(1:nrho2,1:nrho2)
call zgeev('n', 'v', nrho2, zmat9, nrho2, zval0, zmat1, nrho2, zmat2, nrho2, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error! diagonalization failed 2 ', info
   stop
end if
!
zmat1(1:nrho2,1:nrho2) = zmat2(1:nrho2,1:nrho2)
call zgetrf(nrho2, nrho2, zmat1, nrho2, ipiv, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error message for zgetrf 2 ', info
   stop
end if
call zgetri(nrho2, zmat1, nrho2, ipiv, zwork(1,1), nwork, info)
if (info .ne. 0) then
   write(6,*)'calc_qst_hb: error message for zgetri 2 ', info
   stop
end if
!
zmat9(1:nrho2,1:nrho2) = zmat3(1:nrho2,1:nrho2)
do nj=1,nrho2
   do ni=1,nrho2
      zmat3(ni,nj) = zmat1(ni,nj) * zval0(ni)
   end do
end do
zmat4(1:nrho2,1:nrho2) = -zmat9(1:nrho2,1:nrho2)
call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zmat2, nrho2, zmat3, nrho2, cunity, zmat4, nrho2)
!
call cmaxmat(nrho2, nrho2, zmat4, nrho2, dtmp1)
write(6,*)
write(6,*)'calc_qst_hb: check-2 |A-X*Ad*X^(-1)|, max-diff ', dtmp1
call flush(6)
!
! zmat9 --> A
! zmat1 --> X^(-1)
! zmat2 --> X
!
nleig_hb = nrho2
!
zdiag_hb(1:nrho2,1:nmode_hb) = czero
zdiagp_hb(1:nrho2,1:nmode_hb) = czero
zdiagc_hb(1:nrho2,1:nmode_hb) = czero
zdiagpc_hb(1:nrho2,1:nmode_hb) = czero
!
do imode=1,nmode_hb
   do icor=1,ncor_hb
      do ni=1,nleig_hb
         tmpval(ni) = -1.d0 / (cgamma_hb(icor,imode) + zval0(ni)) 
         zdiag_hb(ni,imode) = zdiag_hb(ni,imode) + cb_hb(icor,imode) * tmpval(ni)
      end do
      do ni=1,nleig_hb
         tmpval(ni) = -1.d0 / (dconjg(cgamma_hb(icor,imode)) + zval0(ni)) 
         zdiagc_hb(ni,imode) = zdiagc_hb(ni,imode) + dconjg(cb_hb(icor,imode)) * tmpval(ni)
      end do
   end do
   if (itype_hb .eq. 2) then
      do idrude=1,ndrude_hb
         icor = idrude + npade_hb
         do ni=1,nleig_hb
            tmpval(ni) = 1.d0 / (cgamma_hb(icor,imode) + zval0(ni))**2
            zdiagp_hb(ni,imode) = zdiagp_hb(ni,imode) + cd_hb(idrude,imode) * tmpval(ni)
         end do
         do ni=1,nleig_hb
            tmpval(ni) = 1.d0 / (dconjg(cgamma_hb(icor,imode)) + zval0(ni))**2
            zdiagpc_hb(ni,imode) = zdiagpc_hb(ni,imode) + dconjg(cd_hb(idrude,imode)) * tmpval(ni)
         end do
      end do
   end if
end do
!
do imode=1,nmode_hb
   call zgemv('n', nrho2, nrho2, cunity, zmat1, nrho2, qa_hb(1,1,imode), 1, czero, cmtmp1, 1)
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         ntmp1 = ntmp1 + 1
         cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiag_hb(ntmp1,imode)
      end do
   end do
   call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, czero, zqpst_hb(1,1,imode), 1)
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         ntmp1 = ntmp1 + 1
         cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagc_hb(ntmp1,imode)
      end do
   end do
   call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, czero, zqpcst_hb(1,1,imode), 1)
   if (itype_hb .eq. 2) then
      call zgemv('n', nrho2, nrho2, cunity, zmat1, nrho2, qa_hb(1,1,imode), 1, czero, cmtmp1, 1)
      ntmp1 = 0
      do nj=1,nrho
         do ni=1,nrho
            ntmp1 = ntmp1 + 1
            cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagp_hb(ntmp1,imode)
         end do
      end do
      call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, cunity, zqpst_hb(1,1,imode), 1)
      ntmp1 = 0
      do nj=1,nrho
         do ni=1,nrho
            ntmp1 = ntmp1 + 1
            cmtmp3(ni,nj) = cmtmp1(ni,nj) * zdiagpc_hb(ntmp1,imode)
         end do
      end do
      call zgemv('n', nrho2, nrho2, cunity, zmat2, nrho2, cmtmp3, 1, cunity, zqpcst_hb(1,1,imode), 1)
   end if
end do

!
deallocate(cmat0, cmat1, STAT=istat)
deallocate(zmat0, zmat1, zval0, zwork, rwork, tmpval, STAT=istat)
deallocate(zmat2, zmat3, zmat4, ipiv, STAT=istat)
deallocate(zmat9, STAT=istat)
!
!write(6,*)
!write(6,*)'calc_qst_hb: leaving subroutine'
!call flush(6)
!
return
end subroutine calc_qst_hb
