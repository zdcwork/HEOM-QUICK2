subroutine calc_ax_tfqmr_spa(ldim, rcgtmp, rcgrhs)
!
! purpose : calculate the RHS of EOM for rho 
! note    : matrix is stored in COO sparse format
!    
use matmod
use tmpmatmod
use omp_lib
use matmod_omp
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8,  intent(in)    :: ldim
complex*16, intent(in)    :: rcgtmp(*)
complex*16, intent(out)   :: rcgrhs(*)
integer                   :: itier, itype, dim1, ifactfront, ifactrear, istat
integer                   :: index_nk(MAXTIER)
integer*8                 :: lni, lnj, lnk, lnl, lnm, lnp, lnq
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, jdraw
real*8                    :: cpu1, cpu2
real*8                    :: dtmp1, dtmp2
complex*16                :: cgamma_nk
complex*16                :: ctmp1, ctmp2
character*1               :: transa, transb
integer                   :: nnz, nnz2
integer,    allocatable   :: irow_tmp1(:), icol_tmp1(:)
complex*16, allocatable   :: cval_tmp1(:), cval_tmp2(:), cval_tmp3(:), cval_tmp4(:)
!
integer    :: itimes
data itimes /0/
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'calc_ax_tfqmr_spa: error entry, lsparse = ', lsparse
   stop
end if
if (ldim .ne. lunk_spa) then
   write(6,*)
   write(6,*)'calc_ax_tfqmr_spa: error dimension, ldim != lunk_spa ', ldim, lunk_spa
   stop
end if
if (nthreads_omp .gt. 1) then
   call calc_ax_tfqmr_spa_omp(ldim, rcgtmp, rcgrhs)
   return
end if
if (itimes .eq. 0) then
   write(6,*)
   write(6,*)'calc_ax_tfqmr_spa: first entry'
end if
!
dim1   = nrho**2
itimes = itimes + 1
allocate(irow_tmp1(dim1), icol_tmp1(dim1), STAT=istat)
allocate(cval_tmp1(dim1), cval_tmp2(dim1), STAT=istat)
allocate(cval_tmp3(dim1), cval_tmp4(dim1), STAT=istat)
cval_tmp1(1:dim1) = czero
cval_tmp2(1:dim1) = czero
cval_tmp3(1:dim1) = czero
cval_tmp4(1:dim1) = czero
!
lni = ind_spa(1)
call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rcgtmp(lni), cmtmp1, nrho, nrho, nrho)
if (lhf) then 
   call calchs(cmtmp1)
end if
if (igroundsteady .eq. 0) then
   cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d20, cmtmp1)
   cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
else
   write(6,*)
   write(6,*)'calc_ax_tfqmr_spa: error! unknown igroundsteady ', igroundsteady
   stop
end if
!
!lnl = 1
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
rcgrhs(1:ldim) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
!
! -i (H * rho - rho * H) = -i * H * rho + i * rho * H
lni = ind_spa(1)
call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rcgtmp(lni), irow_spa(lni), icol_spa(lni), &
                 nnz_spa(1), cmtmp3, nrho,  czero, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rcgtmp(lni), irow_spa(lni), icol_spa(lni), &
                 nnz_spa(1), cmtmp3, nrho, cunity, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_extract_coo(nnz_spa(1), irow_spa(lni), icol_spa(lni), rcgrhs(lni),           &
                      cmtmp2, nrho, nrho, nrho)
!
cmtmp2(1:nrho,1:nrho) = czero
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4(1:nrho,1:nrho) = czero
    do ialf=1,nalf
      do imats=1,ncor
        if (lscale) then
           dtmp1 = dbsqrt(iorbs1,ispin,imats,ialf,1)
        else
           dtmp1 = 1.d0
        end if
        call look4drawer(1, ialf, iorbs1, ispin, imats, nn)
        lni = nfirst(2) - 1 + nn
        lnj = ind_spa(lni)
        nnz = nnz_spa(lni)
        call zmat_coo2dns(nnz, irow_spa(lnj), icol_spa(lnj), rcgtmp(lnj), cmtmp1, nrho, nrho, nrho)
        cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + dtmp1 * cmtmp1(1:nrho,1:nrho)
      end do
    end do
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      cmtmp1(ni,nj) = cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
   end do
end do
lni = ind_spa(1)
nnz = nnz_spa(1)
call zmat_extract_coo(nnz, irow_spa(lni), icol_spa(lni), cval_tmp1, cmtmp1, nrho, nrho, nrho)
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lni))
!
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

  do lni=nfirst(itier), nlast(itier)
     if (nnz_spa(lni) .eq. 0) cycle 
!
     cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
     lnl = index_coef_ref(lni)
     do ni=1,itier-1 
        index_nk(ni) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        cgamma_nk = cgamma_nk + cgamma(index_nk(ni))
!
        if (igroundsteady .ne. 0) then
           cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni)) *                         &
                       eleadinfty(ilead(index_nk(ni)), mspin(index_nk(ni)))
        end if
     end do
!
! same itier contribution (N)
!
     lnp = ind_spa(lni)
     nnz = nnz_spa(lni)
     call zmat_add_coo('n', nnz, cgamma_nk, rcgtmp(lnp), cunity, rcgrhs(lnp))
     call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rcgtmp(lnp), irow_spa(lnp), icol_spa(lnp),   &
                      nnz, cmtmp3, nrho,  czero, cmtmp2, nrho,                                &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rcgtmp(lnp), irow_spa(lnp), icol_spa(lnp),   &
                      nnz, cmtmp3, nrho, cunity, cmtmp2, nrho,                                &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
     call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2(1:nrho,1:nrho) = czero
    do ni=1,itier-1
       lnj        = indexcoef(lnl)
       itype      = itypecoef(lnl)
       ifactfront = ifactfcoef(lnl)
       ifactrear  = ifactrcoef(lnl)
       lnl        = lnl + 1
!------
!cycle ! -
!------
      nn   = index_nk(ni)
      lnq  = ind_spa(lnj)
      nnz2 = nnz_spa(lnj) 
      transa = 'n'
      if (mpm(nn) .eq. 1) transa = 'c'
      transb = 'n'
      if (itype .ne. 0) transb = 'c'
!      
      if (.not. offcor) then
          ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
          ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
          call zamsmm_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgtmp(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          call zamsmm_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgtmp(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
      else
          do iorbs2=1,norbs
             ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
             ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
             ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
             call zamsmm_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgtmp(lnq), nnz2,     &
                             irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
             call zamsmm_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgtmp(lnq), nnz2,     &
                             irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          end do
      end if
    end do
!    
    call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
    call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2(1:nrho,1:nrho) = czero
      do isgn=1,nsgn
        transa = 'n'
        if (isgn .ne. 1) transa = 'c'
!               	
        do ispin=1,nspin
          do iorbs1=1,norbs
             cmtmp4(1:nrho,1:nrho) = czero
             cmtmp5(1:nrho,1:nrho) = czero
             do ialf=1,nalf
                do imats=1,ncor
                   lnj        = indexcoef(lnl)
                   itype      = itypecoef(lnl)
                   ifactfront = ifactfcoef(lnl)
                   ifactrear  = ifactrcoef(lnl)
                   lnl        = lnl + 1
! -eye * dtmp1 * op(rhotmp) at right, eye * dtmp2 * op(rhotmp) at left                   
                   dtmp1      = dble(ifactfront)
                   dtmp2      = dble(ifactrear )
                   if (lscale) then
                      dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                      dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                   end if     
                   ctmp1 = dcmplx(0.d0, -dtmp1)
                   ctmp2 = dcmplx(0.d0,  dtmp2)
                   lnq   = ind_spa(lnj)                   
                   nnz2  = nnz_spa(lnj)
                   if (nnz2 .eq. 0) cycle
                   transb = 'n'
                   if (itype .ne. 0) transb = 'c'
                   call zmat_accu_coo2dns(transb, ctmp1, nnz2, irow_spa(lnq), icol_spa(lnq), rcgtmp(lnq),  &
                                          nrho, cmtmp4, nrho)
                   call zmat_accu_coo2dns(transb, ctmp2, nnz2, irow_spa(lnq), icol_spa(lnq), rcgtmp(lnq),  &
                                          nrho, cmtmp5, nrho)
                   !call zamsmm_coo('l', transa, transb, iorbs1, ispin, ctmp1, rcgtmp(lnq), nnz2,    &
                   !                irow_spa(lnq), icol_spa(lnq),                                    &
                   !                cunity, cmtmp2, nrho, cmtmp8, nrho)
                   !call zamsmm_coo('r', transa, transb, iorbs1, ispin, ctmp2, rcgtmp(lnq), nnz2,    &   
                   !                irow_spa(lnq), icol_spa(lnq),                                    &
                   !                cunity, cmtmp2, nrho, cmtmp8, nrho)
                end do
             end do
             call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4, nrho, cunity,  &
                              rcgrhs(lnp), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8, nrho)
             call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5, nrho, cunity,  &
                              rcgrhs(lnp), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8, nrho)
          end do
        end do
      end do
!      call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
!      call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))      
    end if
!
  end do
!
  call cpu_time(cpu2)
!  write(6,100)itier, nlast(itier), cpu2 - cpu1
!  call flush(6)
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
end do 
!
nnz = nnz_spa(1)
lni = ind_spa(1)
rcgrhs(1) = czero
do ni=1,nnz
   if (irow_spa(lni-1+ni) .eq. icol_spa(lni-1+ni)) then
      rcgrhs(lni) = rcgrhs(lni) + rcgtmp(lni-1+ni)
   end if
end do
!
300 continue
!
deallocate(irow_tmp1, icol_tmp1, cval_tmp1, cval_tmp2, cval_tmp3, cval_tmp4, STAT=istat)
!
return
end subroutine calc_ax_tfqmr_spa
