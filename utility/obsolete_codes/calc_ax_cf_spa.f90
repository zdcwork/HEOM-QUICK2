subroutine calc_ax_cf_spa(ldim, rcga, rcgb, rcgrhs, iop)
!
! purpose : calculate the RHS of EOM for rho
!    
! input   : rcga, rcgb, iop
! output  : rcgrhs (iop = 0, rcga has H symmetry)
!                  (iop = 1, rcga has A symmetry)
!
!         
use matmod_omp
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8,  intent(in)  :: ldim
integer,    intent(in)  :: iop
complex*16, intent(in)  :: rcga(*), rcgb(*)
complex*16, intent(out) :: rcgrhs(*)
integer                 :: itier, itype, dim1, ifactfront, ifactrear, istat
integer                 :: index_nk(MAXTIER)
integer*8               :: lni, lnj, lnk, lnl, lnm, lnp, lnq
integer                 :: ni, nj, nk, nl, nm, nn, ntmp3
integer                 :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, jdraw
real*8                  :: cpu1, cpu2
real*8                  :: dtmp1, dtmp2, dtmp3, dtmp4
complex*16              :: cgamma_nk
complex*16              :: zfront, zrear
complex*16              :: ctmp1, ctmp2
character*1             :: transa, transb
integer                 :: nnz, nnz2
integer,    allocatable :: irow_tmp1(:), icol_tmp1(:)
complex*16, allocatable :: cval_tmp1(:), cval_tmp2(:), cval_tmp3(:), cval_tmp4(:)
integer                 :: itimes
data itimes /0/
!
if (iop .ne. 0 .and. iop .ne. 1) then
   write(6,*)'calc_ax_cf_spa: error! unknown iop = ', iop
   stop
end if
if (.not. lsparse) then
   write(6,*)
   write(6,*)'calc_ax_cf_spa: error entry, lsparse = ', lsparse
   stop
end if
if (ldim .ne. lunkcf_spa) then
   write(6,*)
   write(6,*)'calc_ax_cf_spa: error dimension ', ldim, lunkcf_spa
   stop
end if
if (nthreads_omp > 1) then
   call calc_ax_cf_spa_omp(ldim, rcga, rcgb, rcgrhs, iop)
   return 
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
cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!
!lnl = 1
!
if (lfreq) then
   if (lplus) then
      dtmp3 = 1.d0     
   else 
      dtmp3 = -1.d0
   end if
   if (iop .eq. 0) then
      dtmp4 = -1.d0
   else 
      dtmp4 = 1.d0
   end if
   dtmp3 = dtmp3 * dtmp4
end if
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
lni = indcf_spa(1)
nnz = nnzcf_spa(1)
call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rcga(lni), irowcf_spa(lni), icolcf_spa(lni), &
                 nnzcf_spa(1), cmtmp3, nrho,  czero, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rcga(lni), irowcf_spa(lni), icolcf_spa(lni), &
                 nnzcf_spa(1), cmtmp3, nrho, cunity, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_extract_coo(nnz, irowcf_spa(lni), icolcf_spa(lni), rcgrhs(lni),                &
                      cmtmp2, nrho, nrho, nrho)
!
if (lfreq) then
   call zmat_add_coo('n', nnz, dcmplx(dtmp3*dfreq, 0.d0), rcgb(lni), cunity, rcgrhs(lni))
end if
!
!goto 333
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
        lnj = indcf_spa(lni)
        nnz = nnzcf_spa(lni)
        call zmat_coo2dns(nnz, irowcf_spa(lnj), icolcf_spa(lnj), rcgb(lnj), cmtmp1, nrho, nrho, nrho)
        cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + dtmp1 * cmtmp1(1:nrho,1:nrho)
      end do
    end do
    if (iop .eq. 0) then
       call zamsmm1('l', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
       call zamsmm1('r', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
    else 
       call zamsmm1('l', 'n', 'n', iorbs1, ispin,  cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
       call zamsmm1('r', 'n', 'n', iorbs1, ispin,  cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
    end if
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      cmtmp1(ni,nj) = cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
   end do
end do
lni = indcf_spa(1)
nnz = nnzcf_spa(1)
call zmat_extract_coo(nnz, irowcf_spa(lni), icolcf_spa(lni), cval_tmp1, cmtmp1, nrho, nrho, nrho)
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lni))
!
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)
!
  do lni=nfirst(itier), nlast(itier)
     if (nnzcf_spa(lni) .eq. 0) cycle 
!
     lnl = index_coef_ref(lni)
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
     cgamma_nk = czero
     do ni=1,itier-1 
        index_nk(ni) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        cgamma_nk    = cgamma_nk + cgamma(index_nk(ni))
        if (igroundsteady .eq. 1) then
           cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni)) * eleadinfty(ilead(index_nk(ni)), mspin(index_nk(ni)))
        end if
     end do
!
! same itier contribution (N)
!
     lnp = indcf_spa(lni)
     nnz = nnzcf_spa(lni)
     call zmat_add_coo('n', nnz, cgamma_nk, rcga(lnp), cunity, rcgrhs(lnp))
     call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rcga(lnp), irowcf_spa(lnp), icolcf_spa(lnp),   &
                      nnz, cmtmp3, nrho,  czero, cmtmp2, nrho,                                  &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rcga(lnp), irowcf_spa(lnp), icolcf_spa(lnp),   &
                      nnz, cmtmp3, nrho, cunity, cmtmp2, nrho,                                  &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_extract_coo(nnz, irowcf_spa(lnp), icolcf_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
     call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))      
     if (lfreq) then
        call zmat_add_coo('n', nnz, dcmplx(dtmp3*dfreq, 0.d0), rcgb(lnp), cunity, rcgrhs(lnp))
     end if
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
       if (iop .eq. 0) then 
          zfront = -eye
          zrear  = eye 
       else 
          zfront = eye 
          zrear  = -eye
       end if
!------
!cycle ! -
!------
       nn   = index_nk(ni)
       lnq  = indcf_spa(lnj)
       nnz2 = nnzcf_spa(lnj) 
       transa = 'n'
       if (mpm(nn) .eq. 1) transa = 'c'
       transb = 'n'
       if (itype .ne. 0) transb = 'c'
!
       if (.not. offcor) then
          ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * zfront 
          ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * zrear
          call zamsmm_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgb(lnq), nnz2,   &
                          irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          call zamsmm_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgb(lnq), nnz2,   &
                          irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
       else
          do iorbs2=1,norbs
             ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
             ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * zfront
             ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * zrear
             call zamsmm_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgb(lnq), nnz2,       &
                             irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
             call zamsmm_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgb(lnq), nnz2,       &
                             irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          end do
       end if
    end do
    call zmat_extract_coo(nnz, irowcf_spa(lnp), icolcf_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
    call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      !cmtmp2(1:nrho,1:nrho) = czero
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
!
if (lfilter .and. itier .eq. ntier-1) then
   ntmp3 = filtertable(lni)
   if (kpmats(imats) .le. nfilter_long) ntmp3 = ntmp3 + 1
   if (ntmp3 .lt. nfilter_count) cycle
end if
!
                 lnj        = indexcoef(lnl)
                 itype      = itypecoef(lnl)
                 ifactfront = ifactfcoef(lnl)
                 ifactrear  = ifactrcoef(lnl)
                 lnl        = lnl + 1
                 if (iop .eq. 0) then 
                    zfront = -eye
                    zrear  = eye 
                 else 
                    zfront = eye 
                    zrear  = -eye
                 end if
!------
!cycle ! +
!------
                 dtmp1 = dble(ifactfront)
                 dtmp2 = dble(ifactrear )
                 if (lscale) then
                    dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if     
                 ctmp1 = dcmplx(0.d0, -dtmp1) * zfront
                 ctmp2 = dcmplx(0.d0,  dtmp2) * zrear
                 lnq   = indcf_spa(lnj)                   
                 nnz2  = nnzcf_spa(lnj)
                 if (nnz2 .eq. 0) cycle
                 transb = 'n'
                 if (itype .ne. 0) transb = 'c'
                 call zmat_accu_coo2dns(transb, ctmp1, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcgb(lnq),  &
                                        nrho, cmtmp4, nrho)
                 call zmat_accu_coo2dns(transb, ctmp2, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcgb(lnq),  &
                                        nrho, cmtmp5, nrho)
                 !call zamsmm_coo('l', transa, transb, iorbs1, ispin, ctmp1, rcgb(lnq), nnz2,    &
                 !                irowcf_spa(lnq), icolcf_spa(lnq),                              &
                 !                cunity, cmtmp2, nrho, cmtmp8, nrho)
                 !call zamsmm_coo('r', transa, transb, iorbs1, ispin, ctmp2, rcgb(lnq), nnz2,    &   
                 !                irowcf_spa(lnq), icolcf_spa(lnq),                              &
                 !                cunity, cmtmp2, nrho, cmtmp8, nrho)
               end do
             end do
             call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4, nrho, cunity,   &
                              rcgrhs(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8, nrho)
             call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5, nrho, cunity,   &
                              rcgrhs(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8, nrho)
          end do
        end do
      end do
      !call zmat_extract_coo(nnz, irowcf_spa(lnp), icolcf_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
      !call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lnp))      
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
606 continue
!
deallocate(irow_tmp1, icol_tmp1, cval_tmp1, cval_tmp2, cval_tmp3, cval_tmp4, STAT=istat)
!
return
end subroutine calc_ax_cf_spa
