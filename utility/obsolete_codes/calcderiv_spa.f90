subroutine calcderiv_spa(tt)
!
! purpose : calculate the RHS of EOM for rho
!    
! input   : rhotmp (not changed)
! output  : rhorhs (refreshed)
!         
! note    : other arrays such as rho and rho0 should remain intact
!
use matmod
use tmpmatmod
use matmod_omp
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,  intent(in) :: tt
!
integer    :: itier, itype, ifactfront, ifactrear, istat
integer    :: index_nk(MAXTIER)
integer*8  :: lni, lnj, lnk, lnl, lnm, lnp, lnq
integer    :: ni, nj, nk, nl, nm, nn, ntmp3
integer    :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, dim1
real*8     :: cpu1, cpu2, dtmp1, dtmp2, dtmp3, dtmp4
real*8     :: eshift(maxalf,maxspin)
complex*16 :: cgamma_nk, ctmp1, ctmp2, ctmp3, ctmp4
integer    :: itimes
data itimes /0/
!
character*1               :: transa, transb
integer                   :: nnz, nnz2
integer,    allocatable   :: irow_tmp1(:), icol_tmp1(:)
complex*16, allocatable   :: cval_tmp1(:), cval_tmp2(:), cval_tmp3(:), cval_tmp4(:)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'calcderiv_spa: error entry, lsparse = ', lsparse
   stop
end if
if (lpara_omp) then
   call calcderiv_spa_omp(tt)
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
lni = ind_spa(1)
call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rhotmp_spa(lni), cmtmp1, nrho, nrho, nrho)
if (lhf) then
   call calchs(cmtmp1)
end if
call calcdhs(tt, cmtmp1)
cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!
do ispin=1,nspin
  call getengyshift(tt, 1, ispin, eshift(1,ispin))
end do
if (nalf .eq. 2) then
  do ispin=1,nspin
    call getengyshift(tt, 2, ispin, eshift(2,ispin))
  end do
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
rhorhs_spa(1:lunk_spa) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
!
lni = ind_spa(1)
call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rhotmp_spa(lni), irow_spa(lni), icol_spa(lni), &
                 nnz_spa(1), cmtmp3, nrho,  czero, cmtmp2, nrho,                           &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rhotmp_spa(lni), irow_spa(lni), icol_spa(lni), &
                 nnz_spa(1), cmtmp3, nrho, cunity, cmtmp2, nrho,                           &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_extract_coo(nnz_spa(1), irow_spa(lni), icol_spa(lni), rhorhs_spa(lni),           &
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
        call zmat_coo2dns(nnz, irow_spa(lnj), icol_spa(lnj), rhotmp_spa(lnj), cmtmp1, nrho, nrho, nrho)
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
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhorhs_spa(lni))
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)
  do lni=nfirst(itier),nlast(itier)
     if (nnz_spa(lni) .eq. 0) cycle 
!
     cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
     lnl = index_coef_ref(lni)
     do ni=1,itier-1 
        index_nk(ni) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        cgamma_nk    = cgamma_nk + cgamma(index_nk(ni)) + eye * dipm(index_nk(ni)) *  &
                       eshift(ilead(index_nk(ni)), mspin(index_nk(ni)))
     end do
!
! same itier contribution (N)
!
     lnp = ind_spa(lni)
     nnz = nnz_spa(lni)
     call zmat_add_coo('n', nnz, cgamma_nk, rhotmp_spa(lnp), cunity, rhorhs_spa(lnp))
     call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rhotmp_spa(lnp), irow_spa(lnp), icol_spa(lnp),   &
                      nnz, cmtmp3, nrho,  czero, cmtmp2, nrho,                                    &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rhotmp_spa(lnp), irow_spa(lnp), icol_spa(lnp),   &
                      nnz, cmtmp3, nrho, cunity, cmtmp2, nrho,                                    &
                      cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
     call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
     call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhorhs_spa(lnp))
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
!
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
!
          call zamsmm_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rhotmp_spa(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          call zamsmm_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rhotmp_spa(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
!
      else
          do iorbs2=1,norbs
             ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
             ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
             ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
             call zamsmm_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rhotmp_spa(lnq), nnz2,     &
                             irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
             call zamsmm_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rhotmp_spa(lnq), nnz2,     &
                             irow_spa(lnq), icol_spa(lnq), cunity, cmtmp2, nrho, cmtmp8, nrho)
          end do
      end if
    end do
    call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
    call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhorhs_spa(lnp))
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
                 lnl   = lnl + 1
!
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
                 call zmat_accu_coo2dns(transb, ctmp1, nnz2, irow_spa(lnq), icol_spa(lnq), rhotmp_spa(lnq),  &
                                        nrho, cmtmp4, nrho)
                 call zmat_accu_coo2dns(transb, ctmp2, nnz2, irow_spa(lnq), icol_spa(lnq), rhotmp_spa(lnq),  &
                                        nrho, cmtmp5, nrho)
                 !call zamsmm_coo('l', transa, transb, iorbs1, ispin, ctmp1, rhotmp_spa(lnq), nnz2,    &
                 !                irow_spa(lnq), icol_spa(lnq),                                        &
                 !                cunity, cmtmp2, nrho, cmtmp8, nrho)
                 !call zamsmm_coo('r', transa, transb, iorbs1, ispin, ctmp2, rhotmp_spa(lnq), nnz2,    &   
                 !                irow_spa(lnq), icol_spa(lnq),                                        &
                 !                cunity, cmtmp2, nrho, cmtmp8, nrho)
               end do
             end do
             call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4, nrho, cunity,    &
                              rhorhs_spa(lnp), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8, nrho)
             call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5, nrho, cunity,    &
                              rhorhs_spa(lnp), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8, nrho)
          end do
        end do
      end do
      !call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1, cmtmp2, nrho, nrho, nrho)
      !call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhorhs_spa(lnp))      
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
end subroutine calcderiv_spa
