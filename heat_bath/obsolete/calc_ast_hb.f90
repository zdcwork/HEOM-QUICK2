subroutine calc_ast_hb
use matmod
use tmpmatmod
use sparsemod
!
! purpose: compute A^\sigma_{\alpha \mu s} at stationary state
! 
implicit none
include '../include/sizes'
include '../include/common'
!
integer                   :: istat, nwork, info, iorbs, iorbs2, ispin, icor, ialf, isgn
integer                   :: ni, nj, nk, nm, ntmp1, n1, n2
integer                   :: nrho2, norbs2
integer*8                 :: lni
logical                   :: lherm
real*8                    :: dtmp1
real*8,     allocatable   :: rwork(:), dval0(:)
complex*16, allocatable   :: zwork(:,:)
complex*16, allocatable   :: cmat0(:,:), cmat1(:,:)
complex*16, allocatable   :: zmat0(:,:), zmat1(:,:), tmpval(:)
!
!write(6,*)
!write(6,*)'calc_ast_hb: entering subroutine'
!call flush(6)
!
if (.not. lhb) then
   write(6,*)
   write(6,*)'calc_ast_hb: error entrance! '
   stop
end if
!
if (.not. allocated(zast_hb)) then
   allocate(zast_hb(nrho,nrho,norbs,nspin,nsgn), STAT=istat)
end if
zast_hb = czero
!
! diagonalize system Liouvillian
!
nrho2  = nrho**2
norbs2 = norbs**2
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
if (igroundsteady .ne. 0) then
   call calcdhs(1.d20, cmat0(1,1))
   cmat1(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
end if
!
nwork = nrho2**2
allocate(zmat0(nrho2,nrho2), zmat1(nrho2,nrho2), dval0(nrho2), STAT=istat)
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
!
call zhil2liou('l', 'n', nrho, cmat1, nrho, zmat0, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zmat1, nrho2)
zmat0(1:nrho2,1:nrho2) = zmat0(1:nrho2,1:nrho2) - zmat1(1:nrho2,1:nrho2)
call checkhermicity(zmat0, nrho2, zmat1, nrho2, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)'calc_ast_hb: error! Liouville matrix is not Hermitian '
   stop
end if
!
call zheev('V', 'u', nrho2, zmat0, nrho2, dval0, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'calc_ast_hb: error! diagonalization failed 1 ', info
   stop
end if
!
do isgn=1,nsgn
   do ispin=1,nspin
      do iorbs2=1,norbs
         if (isgn .eq. 1) then
            do n2=1,nrho
               do n1=1,nrho
                  cmtmp1(n1,n2) = dcmplx(amsori(n2,n1,iorbs2,ispin), 0.d0)
               end do
            end do
         else
            cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs2,ispin), 0.d0)
         end if
         call zgemv('c', nrho2, nrho2, cunity, zmat0, nrho2, cmtmp1, 1, czero, cmtmp2, 1)
!
         do ialf=1,nalf
            if (nvar2 .eq. norbs2) then
               do iorbs=1,norbs
                  ntmp1 = (iorbs2 - 1) * norbs + iorbs
                  tmpval(1:nrho2) = czero
                  do icor=1,ncor
                     tmpval(1:nrho2) = tmpval(1:nrho2) - cb(ntmp1,ispin,icor,ialf,isgn) /       & 
                                      ( cgama(ntmp1,ispin,icor,ialf,isgn) + eye * dipm(isgn) *  &
                                        eleadinfty(ialf,ispin) - eye * dval0(1:nrho2) )
                  end do
                  ntmp1 = 0
                  do n2=1,nrho
                     do n1=1,nrho
                        ntmp1 = ntmp1 + 1
                        cmtmp3(n1,n2) = cmtmp2(n1,n2) * tmpval(ntmp1)
                     end do
                  end do
                  call zgemv('n', nrho2, nrho2, cunity, zmat0, nrho2, cmtmp3, 1, cunity, zast_hb(1,1,iorbs,ispin,isgn), 1)
               end do
            else
               tmpval(1:nrho2) = czero
               do icor=1,ncor
                  tmpval(1:nrho2) = tmpval(1:nrho2) - cb(iorbs2,ispin,icor,ialf,isgn) /       & 
                                   ( cgama(iorbs2,ispin,icor,ialf,isgn) + eye * dipm(isgn) *  &
                                     eleadinfty(ialf,ispin) - eye * dval0(1:nrho2) )
               end do
               ntmp1 = 0
               do n2=1,nrho
                  do n1=1,nrho
                     ntmp1 = ntmp1 + 1
                     cmtmp3(n1,n2) = cmtmp2(n1,n2) * tmpval(ntmp1)
                  end do
               end do
               call zgemv('n', nrho2, nrho2, cunity, zmat0, nrho2, cmtmp3, 1, cunity, zast_hb(1,1,iorbs2,ispin,isgn), 1)
            end if
         end do  
!
      end do
   end do
end do
!
deallocate(cmat0, cmat1, STAT=istat)
deallocate(zmat0, zmat1, dval0, zwork, rwork, tmpval, STAT=istat)
!
!write(6,*)
!write(6,*)'calc_ast_hb: leaving subroutine'
!call flush(6)
!
return
end subroutine calc_ast_hb
