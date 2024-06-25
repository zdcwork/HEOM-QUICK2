subroutine calc_ax_cf_spa_omp(ldim, rcga, rcgb, rhsa, rhsb)
!
! purpose : calculate the RHS of EOM for rho
!    
use matmod
use tmpmatmod
use omp_lib
use matmod_omp
use tmpmatmod_omp
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8,  intent(in)    :: ldim
complex*16, intent(in)    :: rcga(*), rcgb(*)     
complex*16, intent(out)   :: rhsa(*), rhsb(*)
logical                   :: lslow1, lslow2, lfast1, lfast2
integer                   :: itier, itype, dim1, ifactfront, ifactrear, ifactf, ifactr
integer                   :: index_nk(MAXTIER, maxprocs), index_nj(MAXTIER, maxprocs)
integer                   :: index_np(MAXTIER, maxprocs)
integer                   :: nvec_tmp1_omp(maxtier, maxprocs), nvec_tmp2_omp(maxtier, maxprocs)
integer                   :: nvec_tmp3_omp(maxtier, maxprocs), nvec_tmp4_omp(maxtier, maxprocs)
integer*8                 :: lni, lnj, lnk, lnl, lnm, lnp, lnq, lnr
integer                   :: ni, nj, nk, nl, nm, nn, ntmp1, ntmp2, ntmp3
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, jdraw
integer                   :: iopr, isopr
integer                   :: imode, iballs
integer                   :: ndegen, nfast
integer                   :: nnz, nnz2
real*8                    :: cpu1, cpu2, cpu3, cpu4
real*8                    :: dtmp1, dtmp2, dtmp3, dtmp4
real*8                    :: degen
real*8                    :: dsmall_spa
complex*16                :: ctmp1, ctmp2, ctmp3, ctmp4
character*1               :: transa, transb, transc
integer                   :: nrho2, mi, mj, istat, mdim, isgnw, nrho2d
real*8                    :: dfacta, dfactb
complex*16                :: cgamma_nk, cgamma_nj
complex*16                :: zfront, zrear
complex*16                :: zfa2b, zfb2a, zra2b, zrb2a   ! a/b stands for rcga/rcgb
integer,    allocatable   :: irow_tmp1(:), icol_tmp1(:)
complex*16, allocatable   :: cval_tmp1(:)
integer,    allocatable   :: irow_tmp1_omp(:,:), icol_tmp1_omp(:,:)
integer,    allocatable   :: irow_tmp2_omp(:,:), icol_tmp2_omp(:,:)
complex*16, allocatable   :: cval_tmp1_omp(:,:), cval_tmp2_omp(:,:)
integer                   :: nnz_heff
integer,    allocatable   :: irow_heff(:), icol_heff(:)
complex*16, allocatable   :: cval_heff(:)
integer                   :: ithread, icore
integer*8                 :: ifrac
external ifrac
!
integer                   :: itimes
data itimes /0/
!
call cpu_time(cpu3)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'calc_ax_cf_spa_omp: error entry ', lsparse
   stop
end if
if (ldim .ne. lunkcf_spa) then
   write(6,*)
   write(6,*)'calc_ax_cf_spa_omp: error dimension, ldim, lunkcf_spa'
   write(6,*)'calc_ax_cf_spa_omp: ', ldim, lunkcf_spa
   stop
end if
!
if (itimes .eq. 0) then
   write(6,*)'calc_ax_cf_spa_omp: first entry'
end if
if (lproj) then
    dsmall_spa = dpico
else
    dsmall_spa = dsmall
end if
!
dim1   = nrho**2
itimes = itimes + 1
allocate(irow_tmp1(dim1), icol_tmp1(dim1), cval_tmp1(dim1), STAT=istat)
allocate(irow_tmp1_omp(dim1,nthreads_omp), icol_tmp1_omp(dim1,nthreads_omp), STAT=istat)
allocate(irow_tmp2_omp(dim1,nthreads_omp), icol_tmp2_omp(dim1,nthreads_omp), STAT=istat)
allocate(cval_tmp1_omp(dim1,nthreads_omp), cval_tmp2_omp(dim1,nthreads_omp), STAT=istat)
cval_tmp1(1:dim1)                    = czero
cval_tmp1_omp(1:dim1,1:nthreads_omp) = czero
cval_tmp2_omp(1:dim1,1:nthreads_omp) = czero
!
cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!
allocate(irow_heff(dim1), icol_heff(dim1), cval_heff(dim1), STAT=istat)
call zmat_dns2coo(nnz_heff, irow_heff, icol_heff, cval_heff, cmtmp3, nrho, nrho, nrho, dsmall_spa)
!
isgnw = 1
if (lfreq) then
   if (lplus) then
      dtmp3 = 1.d0     
      isgnw = 1
   else 
      dtmp3 = -1.d0
      isgnw = 2
   end if
   dfacta = -dtmp3
   dfactb =  dtmp3
end if
!
if (ltrun_der .or. lad) then
   mdim  = MAXTIER 
   nrho2 = nrho * nrho  
   nrho2d = nrho2 * 2
end if
!
zfb2a = -eye
zrb2a =  eye
zfa2b =  eye
zra2b = -eye
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
rhsa(1:ldim) = czero
rhsb(1:ldim) = czero
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
call zmat_extract_coo(nnzcf_spa(1), irowcf_spa(lni), icolcf_spa(lni), rhsa(lni),         &
                      cmtmp2, nrho, nrho, nrho)
call zmat_zcoomm('r', 'n', 'n', nrho, -eye, rcgb(lni), irowcf_spa(lni), icolcf_spa(lni), &
                 nnzcf_spa(1), cmtmp3, nrho,  czero, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_zcoomm('l', 'n', 'n', nrho,  eye, rcgb(lni), irowcf_spa(lni), icolcf_spa(lni), &
                 nnzcf_spa(1), cmtmp3, nrho, cunity, cmtmp2, nrho,                       &
                 cval_tmp1, irow_tmp1, icol_tmp1, cmtmp7, nrho, cmtmp8, nrho)
call zmat_extract_coo(nnzcf_spa(1), irowcf_spa(lni), icolcf_spa(lni), rhsb(lni),         &
                      cmtmp2, nrho, nrho, nrho)
if (lfreq) then
   call zmat_add_coo('n', nnz, dcmplx(dfacta*dfreq, 0.d0), rcgb(lni), cunity, rhsa(lni))
   call zmat_add_coo('n', nnz, dcmplx(dfactb*dfreq, 0.d0), rcga(lni), cunity, rhsb(lni))
end if
!
!goto 333
cmtmp2(1:nrho,1:nrho) = czero
cmtmp1(1:nrho,1:nrho) = czero
!
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4(1:nrho,1:nrho) = czero
    cmtmp5(1:nrho,1:nrho) = czero
    do ialf=1,nalf
       call look4operator(1, iorbs1, ispin, ialf, iopr)
       if (lscreen .and. jomit(iopr) .eq. 1) cycle
!
!$omp parallel default(shared) private(imats,dtmp1,nn,lni,lnj,nnz,ithread)
      ithread = omp_get_thread_num() + 1
      cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
      cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
!$omp do schedule(dynamic) 
      do imats=1,ncor
        if (lscale) then
           dtmp1 = dbsqrt(iorbs1,ispin,imats,ialf,1)
        else
           dtmp1 = 1.d0
        end if
        call look4drawer(1, ialf, iorbs1, ispin, imats, nn)
        lni = nfirst(2) - 1 + indexdraw(nn)
        lnj = indcf_spa(lni)
        nnz = nnzcf_spa(lni)
        call zmat_accu_coo2dns('n', dcmplx(dtmp1,0.d0), nnz, irowcf_spa(lnj), icolcf_spa(lnj), &
                               rcgb(lnj), nrho, cmtmp4_omp(1,1,ithread), nrho)
        call zmat_accu_coo2dns('n', dcmplx(dtmp1,0.d0), nnz, irowcf_spa(lnj), icolcf_spa(lnj), &
                               rcga(lnj), nrho, cmtmp5_omp(1,1,ithread), nrho)
      end do
!$omp end do
!$omp end parallel
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
         cmtmp5(1:nrho,1:nrho) = cmtmp5(1:nrho,1:nrho) + cmtmp5_omp(1:nrho,1:nrho,icore)
      end do
    end do
!    call zamsmm_spa_omp('l', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho, &
!                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
!    call zamsmm_spa_omp('r', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho, &
!                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
!    call zamsmm_spa_omp('l', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho, &
!                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
!    call zamsmm_spa_omp('r', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho, &
!                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
    if (lprams) then
        call zamsmm1_proj('l', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
        call zamsmm1_proj('r', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
        call zamsmm1_proj('l', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
        call zamsmm1_proj('r', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
    else
        call zamsmm1('l', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
        call zamsmm1('r', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
        call zamsmm1('l', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
        call zamsmm1('r', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
    end if
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      cmtmp7(ni,nj) = cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
      cmtmp8(ni,nj) = cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni))
   end do
end do
lni = indcf_spa(1)
nnz = nnzcf_spa(1)
call zmat_extract_coo(nnz, irowcf_spa(lni), icolcf_spa(lni), cval_tmp1, cmtmp7, nrho, nrho, nrho)
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhsa(lni))
call zmat_extract_coo(nnz, irowcf_spa(lni), icolcf_spa(lni), cval_tmp1, cmtmp8, nrho, nrho, nrho)
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rhsb(lni))
!
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier0
  call cpu_time(cpu1)
  iballs = itier - 1

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn) &
!$omp private(ntmp3,isgn,ispin,iorbs1,iorbs2,ialf,imats,dtmp1,dtmp2,i,j,ithread,ifactf,ifactr,imode)  &
!$omp private(lnp,lnq,nnz,nnz2,transa,transb,transc,ctmp1,ctmp2,zfront,zrear,idraw,cgamma_nj,mi,mj)   & 
!$omp private(iopr,isopr,lnm,lnr,ndegen,degen,nfast,ntmp1,ntmp2,dtmp3,dtmp4,ctmp3,ctmp4,lslow1,lslow2,lfast1,lfast2)
!
  ithread = omp_get_thread_num() + 1
!
!$omp do schedule(dynamic) 
  do lni=nfirst(itier), nlast(itier)
     if (nnzcf_spa(lni) .eq. 0) cycle
!
     cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
     lnl = index_coef_ref(lni)
     if (lsimple) lnm = ioprref(lni)
!
     do ni=1,itier-1 
        index_nk(ni,ithread) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        cgamma_nk = cgamma_nk + cgamma(index_nk(ni,ithread))
        if (igroundsteady .eq. 1) then
           cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                             &
                       eleadinfty(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
        end if
     end do
     !
     if (lad .and. itier .eq. ntier) then
         lslow1 = .true.
         do ni=1,itier-1
            if (mslow(index_nk(ni,ithread)) .ne. 1) then
                lslow1 = .false.
                exit
            end if
         end do
         lfast1 = .true.
         do ni=1,itier-1
            if (mslow(index_nk(ni,ithread)) .eq. 1) then
                lfast1 = .false.
                exit
            end if
         end do
     end if
!
! same itier contribution (N)
!
     lnp = indcf_spa(lni)
     nnz = nnzcf_spa(lni)
     call zmat_add_coo('n', nnz, cgamma_nk, rcga(lnp), cunity, rhsa(lnp))
     call zmat_add_coo('n', nnz, cgamma_nk, rcgb(lnp), cunity, rhsb(lnp))
     if (lfreq) then
        call zmat_add_coo('n', nnz, dcmplx(dfacta*dfreq, 0.d0), rcgb(lnp), cunity, rhsa(lnp))
        call zmat_add_coo('n', nnz, dcmplx(dfactb*dfreq, 0.d0), rcga(lnp), cunity, rhsb(lnp))
     end if
     call zmat_zcoomm2('r', 'n', 'n', nrho, -eye, rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),  &
                       cval_heff, nnz_heff, irow_heff, icol_heff,                                    &
                       cunity, rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),                     &
                       cval_tmp1_omp(1,ithread), irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread), &
                       cval_tmp2_omp(1,ithread), irow_tmp2_omp(1,ithread), icol_tmp2_omp(1,ithread), &
                       cmtmp8_omp(1,1,ithread), nrho)
     call zmat_zcoomm2('l', 'n', 'n', nrho,  eye, rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),  &
                       cval_heff, nnz_heff, irow_heff, icol_heff,                                    &
                       cunity, rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),                     &
                       cval_tmp1_omp(1,ithread), irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread), &
                       cval_tmp2_omp(1,ithread), irow_tmp2_omp(1,ithread), icol_tmp2_omp(1,ithread), &
                       cmtmp8_omp(1,1,ithread), nrho)
     call zmat_zcoomm2('r', 'n', 'n', nrho, -eye, rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),  &
                       cval_heff, nnz_heff, irow_heff, icol_heff,                                    &
                       cunity, rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),                     &
                       cval_tmp1_omp(1,ithread), irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread), &
                       cval_tmp2_omp(1,ithread), irow_tmp2_omp(1,ithread), icol_tmp2_omp(1,ithread), &
                       cmtmp8_omp(1,1,ithread), nrho)
     call zmat_zcoomm2('l', 'n', 'n', nrho,  eye, rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),  &
                       cval_heff, nnz_heff, irow_heff, icol_heff,                                    &
                       cunity, rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),                     &
                       cval_tmp1_omp(1,ithread), irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread), &
                       cval_tmp2_omp(1,ithread), irow_tmp2_omp(1,ithread), icol_tmp2_omp(1,ithread), &
                       cmtmp8_omp(1,1,ithread), nrho)
!
! nearest lower tier contribution (N - 1)
!
    cval_tmp1_omp(1:dim1,ithread) = czero
    cval_tmp2_omp(1:dim1,ithread) = czero
    do ni=1,itier-1
       nn         = index_nk(ni,ithread)
       lnj        = indexcoef(lnl)
       itype      = itypecoef(lnl)
       ifactfront = ifactfcoef(lnl)
       ifactrear  = ifactrcoef(lnl)
       lnl        = lnl + 1
!
       if (psdfff .and. mmfff(nn) .gt. 1) then
           ctmp1 = dble(ifactfront) * ctheta(mspin(nn),mcor(nn),ilead(nn),mpm(nn))
           lnp = indcf_spa(lni)
           nnz = nnzcf_spa(lni)
           lnq = indcf_spa(lnj)
           transa = 'n'
           if (itype .ne. 0) then
               write(6,*)
               write(6,*)'calc_ax_cf_spa_omp: unexpected itype ', itype
               stop
           end if
           call zmat_add_coo(transa, nnz, ctmp1, rcga(lnq), cunity, rhsa(lnp))
           call zmat_add_coo(transa, nnz, ctmp1, rcgb(lnq), cunity, rhsb(lnp))
           cycle
       end if
!
       zfront = -eye
       zrear  =  eye 
!------
!cycle ! -
!------
      lnq  = indcf_spa(lnj)
      nnz2 = nnzcf_spa(lnj) 
      transa = 'n'
      if (mpm(nn) .eq. 1) transa = 'c'
      transb = 'n'
      if (itype .ne. 0) transb = 'c'

      if (.not. offcor) then
         ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * zfront
         ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * zrear
         if (lprams) then
             call zamsmm1_coo_proj('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgb(lnq), nnz2,        &
                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                   irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_coo_proj('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgb(lnq), nnz2,        &
                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                   irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
             call zamsmm1_coo_proj('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcga(lnq), nnz2,        &
                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                   irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_coo_proj('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcga(lnq), nnz2,        &
                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                   irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
         else
             call zamsmm1_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgb(lnq), nnz2,        &
                              irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                              irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgb(lnq), nnz2,        &
                              irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                              irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
             call zamsmm1_coo('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcga(lnq), nnz2,        &
                              irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                              irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_coo('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcga(lnq), nnz2,        &
                              irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                              irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
         end if
!
      else
        do iorbs2=1,norbs
           ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
           ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * zfront
           ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * zrear
           if (lprams) then
               call zamsmm1_coo_proj('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgb(lnq), nnz2,           &
                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                     irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
               call zamsmm1_coo_proj('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgb(lnq), nnz2,           &
                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                     irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)         
               call zamsmm1_coo_proj('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcga(lnq), nnz2,           &
                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                     irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
               call zamsmm1_coo_proj('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcga(lnq), nnz2,           &
                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                     irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)  
           else
               call zamsmm1_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgb(lnq), nnz2,           &
                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
               call zamsmm1_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgb(lnq), nnz2,           &
                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                                irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)         
               call zamsmm1_coo('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcga(lnq), nnz2,           &
                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
               call zamsmm1_coo('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcga(lnq), nnz2,           &
                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cval_tmp2_omp(1,ithread), nnz,   &
                                irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)  
           end if
        end do
      end if
    end do
    call zmat_add_coo('n', nnz, cunity, cval_tmp1_omp(1,ithread), cunity, rhsa(lnp))
    call zmat_add_coo('n', nnz, cunity, cval_tmp2_omp(1,ithread), cunity, rhsb(lnp))
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      do isgn=1,nsgn
        transa = 'n'
        if (isgn .ne. 1) transa = 'c'
!
        do ispin=1,nspin
          orbs0: do iorbs1=1,norbs
             if (lsimple .and. .not. lwalf) then
                iopr = ioprindex(lnm)
                lnm = lnm + 1
                if (iopr .lt. 0) cycle orbs0
             end if
             cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp7_omp(1:nrho,1:nrho,ithread) = czero
             alf0: do ialf=1,nalf
               call look4operator(isgn, iorbs1, ispin, ialf, iopr)
               if (lscreen .and. jomit(iopr) .eq. 1) then
                   cycle alf0
               end if
               if (lsimple .and. lwalf) then
                  iopr = ioprindex(lnm)
                  lnm = lnm + 1
                  if (iopr .lt. 0) cycle alf0
               end if
!
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
                 zfront = -eye
                 zrear  =  eye 
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
                                        nrho, cmtmp4_omp(1,1,ithread), nrho)
                 call zmat_accu_coo2dns(transb, ctmp2, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcgb(lnq),  &
                                        nrho, cmtmp5_omp(1,1,ithread), nrho)
                 call zmat_accu_coo2dns(transb, -ctmp1, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcga(lnq),  &
                                        nrho, cmtmp6_omp(1,1,ithread), nrho)
                 call zmat_accu_coo2dns(transb, -ctmp2, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcga(lnq),  &
                                        nrho, cmtmp7_omp(1,1,ithread), nrho)
               end do
             end do alf0
             if (lprams) then
                 call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4_omp(1,1,ithread), nrho, cunity,   &
                                       rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,   &
                                       rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin, cunity, cmtmp6_omp(1,1,ithread), nrho, cunity,   &
                                       rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin, cunity, cmtmp7_omp(1,1,ithread), nrho, cunity,   &
                                       rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             else
                 call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4_omp(1,1,ithread), nrho, cunity,   &
                                  rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,   &
                                  rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp6_omp(1,1,ithread), nrho, cunity,   &
                                  rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                 call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp7_omp(1,1,ithread), nrho, cunity,   &
                                  rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             end if
          end do orbs0
        end do
      end do
    end if
!
! Adiabatic HEOM part 
!
    if (lad) then
        if (itier .eq. ntier) then
            do isgn=1,nsgn
               transa = 'n'
               if (isgn .ne. 1) transa = 'c'
               do ispin=1,nspin
                  do iorbs1=1,norbs
                     call look4sysopr(isgn, iorbs1, ispin, isopr)
                     cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
                     do ialf=1,nalf
                        iopr = (isopr - 1) * nalf + ialf
                        if (lscreen .and. jomit(iopr) .eq. 1) cycle
                        !
                        loop_draw0: do imats=1,ncor
                           idraw = (iopr - 1) * ncor + imats
                           !
                           if (mslow(idraw) .eq. 1) then
                               lslow2 = .true.
                               lfast2 = .false.
                           else
                               lslow2 = .false.
                               lfast2 = .true. 
                           end if
                           if (lad_fast .and. lfast1 .and. lfast2 .or.           &
                               (.not. lad_fast) .and. lslow1 .and. lslow2 ) then
                               if (itier+1 .gt. ntier_ad) cycle loop_draw0
                           end if
                           !
                           lnj        = indexcoef(lnl)
                           itype      = itypecoef(lnl)
                           ifactfront = ifactfcoef(lnl)
                           ifactrear  = ifactrcoef(lnl)
                           lnl        = lnl + 1
                           lnq        = indcf_spa(lnj) 
                           nnz2       = nnzcf_spa(lnj)
                           !
                           if (lnj .gt. nlast(itier)) then  ! all slow modes (or all fast modes)
                               if (nnz2 .eq. 0) cycle
                               dtmp1 = dble(ifactfront)
                               dtmp2 = dble(ifactrear )
                               ctmp1 = dcmplx(dtmp1, 0.d0)
                               ctmp2 = dcmplx(dtmp2, 0.d0)
                               transb = 'n'
                               if (itype .ne. 0) transb = 'c'
                               call zmat_accu_coo2dns_pair(transb, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcgb(lnq),  &
                                                           ctmp1*zfb2a, cmtmp4_omp(1,1,ithread), nrho,                 &
                                                           ctmp2*zrb2a, cmtmp5_omp(1,1,ithread), nrho)
                               call zmat_accu_coo2dns_pair(transb, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcga(lnq),  &
                                                           ctmp1*zfa2b, cmtmp1_omp(1,1,ithread), nrho,                 &
                                                           ctmp2*zra2b, cmtmp2_omp(1,1,ithread), nrho)
                           !
                           else                                ! adiabatic approximation 
                               !
                               if (ifactfront .eq. 0) cycle    ! higher-tier ADO is zero (contains idential drawers)
!                               
                               call add2list(index_nk(1,ithread), mdim, iballs, index_nj(1,ithread), mdim, idraw, ni)
                               mi    = 0
                               ntmp3 = 0
                               do mj=1,itier
                                  if (ntmp3 .lt. ms2f(index_nj(mj,ithread))) then
                                      ntmp3 = ms2f(index_nj(mj,ithread))
                                      mi = mj
                                  end if
                               end do
                               nfast = mi
                               idraw = index_nj(nfast,ithread)
                               ! analyze degenerate fast modes 
                               ndegen = 1
                               nvec_tmp1_omp(1,ithread) = idraw
                               do mi=1,itier
                                  if (mi .eq. nfast) cycle
                                  ntmp3 = index_nj(mi,ithread)
                                  if (dgama_drawer(ntmp3) .ge. dgama_drawer(idraw)*dratio_fast) then
                                      ndegen = ndegen + 1
                                      nvec_tmp1_omp(ndegen,ithread) = ntmp3
                                  end if
                               end do
                               degen = 1.d0 / dble(ndegen)
                               if (idegen_fast .eq. 0 .and. ndegen .gt. 1) cycle
                               !
                               if (idegen_fast .eq. 2 .and. ndegen .gt. 1 .and. ndegen .le. ndegen_fast) then
!debug
!write(6,*)
!write(6,*)'calc_ax_cf_spa_omp: code not ready yet 1 '
!stop
!debug
                                   lnj        = indexcoef(lnl)
                                   itype      = itypecoef(lnl)
                                   ifactf     = ifactfcoef(lnl)  ! ifactf = fact
                                   ifactr     = ifactrcoef(lnl)  ! ifactr = fact
                                   lnl        = lnl + 1
                                   lnq        = indcf_spa(lnj) 
                                   nnz2       = nnzcf_spa(lnj)
                                   if (nnz2 .eq. 0) cycle
!
                                   call findsubvec(index_nj(1,ithread), itier, nvec_tmp1_omp(1,ithread), ndegen,  &
                                                   index_np(1,ithread), ntmp3, nvec_tmp2_omp(1,ithread))
                                   do ni=1,ndegen
                                      ntmp3 = itier - ni + 1  ! length
                                      ntmp1 = (-1)**(ntmp3 - 1)
                                      nvec_tmp3_omp(ni,ithread) = (-1)**(nvec_tmp2_omp(ni,ithread)-1)  ! gl
                                      nvec_tmp4_omp(ni,ithread) = nvec_tmp3_omp(ni,ithread) * ntmp1    ! gr
                                   end do
                                   lnr = ifrac(ndegen)
                                   dtmp3 = 1.d0 / dble(lnr)
                                   !
                                   do lnr=1,ifrac(ndegen)
                                      dtmp4 = dble(npermu_fact(lnr,ndegen))
                                      ntmp2 = ndegen - npermu(1,lnr,ndegen) + 1
                                      ntmp3 = ndegen - npermu(1,  1,ndegen) + 1
                                      !
                                      ! gl*zgefnl*rcgtmp - gr*rcgtmp*zgefnr
                                      idraw = nvec_tmp1_omp(ntmp2,ithread)
                                      ctmp2 = dcmplx(dble(nvec_tmp3_omp(ntmp3,ithread)), 0.d0)  ! gl(ntmp3)
                                      ctmp1 = dcmplx(dble(nvec_tmp4_omp(ntmp3,ithread)), 0.d0)  ! gr(ntmp3)
                                      !
                                      transb = 'n'
                                      if (itype .ne. 0) transb = 'c'
                                      if (mod(ndegen,2) .eq. 0) then
                                          if (lcop_ad) then
                                              call zlmat_zcoomv4(transb, nrho, ctmp2, nnz_lvl(idraw),                                 &
                                                                 irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                                                 rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 czero, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, -ctmp1, nnz_lvr(idraw),                                &
                                                                 irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                                                 rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 cunity, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, ctmp2, nnz_lvl(idraw),                                 &
                                                                 irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                                                 rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 czero, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, -ctmp1, nnz_lvr(idraw),                                &
                                                                 irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                                                 rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 cunity, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                          else              
                                              call zmat_zcoomm4('l', 'n', transb,  ctmp2, cval_zgfl(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                         &
                                                                rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                                czero, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,             &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('r', 'n', transb, -ctmp1, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                         &
                                                                rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                                cunity, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,            &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('l', 'n', transb,  ctmp2, cval_zgfl(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                         &
                                                                rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                                czero, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,             &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('r', 'n', transb, -ctmp1, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                         &
                                                                rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                                cunity, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,            &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                          end if
                                      else
                                          if (lcop_ad) then
                                              call zlmat_zcoomv4(transb, nrho, ctmp2*zfa2b, nnz_lvl(idraw),                           &
                                                                 irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                                                 rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 czero, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, -ctmp1*zra2b, nnz_lvr(idraw),                          &
                                                                 irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                                                 rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 cunity, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, ctmp2*zfb2a, nnz_lvl(idraw),                           &
                                                                 irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                                                 rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 czero, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                              call zlmat_zcoomv4(transb, nrho, -ctmp1*zrb2a, nnz_lvr(idraw),                          &
                                                                 irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                                                 rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                                                 cunity, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                                                 nlvrow, nlvcol)
                                          else             
                                              call zmat_zcoomm4('l', 'n', transb, ctmp2*zfa2b, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                                rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                                czero, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('r', 'n', transb, -ctmp1*zra2b, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                                rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                                cunity, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('l', 'n', transb, ctmp2*zfb2a, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                                rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                                czero, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                              call zmat_zcoomm4('r', 'n', transb, -ctmp1*zrb2a, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                                irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                                rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                                cunity, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                                cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                          end if
                                      end if
                                      !
                                      do mi=2,ndegen   ! mind the order here: closest acts first
                                         ntmp2 = ndegen - npermu(mi,lnr,ndegen) + 1
                                         ntmp3 = ndegen - npermu(mi,  1,ndegen) + 1
                                         idraw = nvec_tmp1_omp(ntmp2,ithread)
                                         ctmp2 = dcmplx(dble(nvec_tmp3_omp(ntmp3,ithread)), 0.d0)  ! gl(ntmp3)
                                         ctmp1 = dcmplx(dble(nvec_tmp4_omp(ntmp3,ithread)), 0.d0)  ! gr(ntmp3)
                                         !
                                         if (lcop_ad) then
                                             call zlmat_zcoomv3('n', nrho, ctmp2, nnz_lvl(idraw),                                     &
                                                                irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),              &
                                                                cmtmp7_omp(1,1,ithread), nrho, czero, cmtmp8_omp(1,1,ithread), nrho,  &
                                                                nlvrow, nlvcol)
                                             call zlmat_zcoomv3('n', nrho, -ctmp1, nnz_lvr(idraw),                                    &
                                                                irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),              &
                                                                cmtmp7_omp(1,1,ithread), nrho, cunity, cmtmp8_omp(1,1,ithread), nrho, &
                                                                nlvrow, nlvcol)
                                             cmtmp7_omp(1:nrho,1:nrho,ithread) = cmtmp8_omp(1:nrho,1:nrho,ithread) 
                                             call zlmat_zcoomv3('n', nrho, ctmp2, nnz_lvl(idraw),                                     &
                                                                irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),              &
                                                                cmtmp6_omp(1,1,ithread), nrho, czero, cmtmp8_omp(1,1,ithread), nrho,  &
                                                                nlvrow, nlvcol)
                                             call zlmat_zcoomv3('n', nrho, -ctmp1, nnz_lvr(idraw),                                    &
                                                                irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),              &
                                                                cmtmp6_omp(1,1,ithread), nrho, cunity, cmtmp8_omp(1,1,ithread), nrho, &
                                                                nlvrow, nlvcol)
                                             cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp8_omp(1:nrho,1:nrho,ithread) 
                                         else
                                             call zmat_zcoomm('l', 'n', 'n', nrho, ctmp2, cval_zgfl(1,idraw),                &
                                                              irow_zgf(1,idraw), icol_zgf(1,idraw), nnz_zgf(idraw),          &
                                                              cmtmp7_omp(1,1,ithread), nrho, czero,                          &
                                                              cmtmp8_omp(1,1,ithread), nrho, cval_tmp1_omp(1,ithread),       &
                                                              irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread),            &
                                                              cmtmpa_omp(1,1,ithread), nrho, cmtmpb_omp(1,1,ithread), nrho)
                                             call zmat_zcoomm('r', 'n', 'n', nrho, -ctmp1, cval_zgfr(1,idraw),               &
                                                              irow_zgf(1,idraw), icol_zgf(1,idraw), nnz_zgf(idraw),          &
                                                              cmtmp7_omp(1,1,ithread), nrho, cunity,                         &
                                                              cmtmp8_omp(1,1,ithread), nrho, cval_tmp1_omp(1,ithread),       &
                                                              irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread),            &
                                                              cmtmpa_omp(1,1,ithread), nrho, cmtmpb_omp(1,1,ithread), nrho)
                                             cmtmp7_omp(1:nrho,1:nrho,ithread) = cmtmp8_omp(1:nrho,1:nrho,ithread) 
                                             call zmat_zcoomm('l', 'n', 'n', nrho, ctmp2, cval_zgfl(1,idraw),                &
                                                              irow_zgf(1,idraw), icol_zgf(1,idraw), nnz_zgf(idraw),          &
                                                              cmtmp6_omp(1,1,ithread), nrho, czero,                          &
                                                              cmtmp8_omp(1,1,ithread), nrho, cval_tmp1_omp(1,ithread),       &
                                                              irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread),            &
                                                              cmtmpa_omp(1,1,ithread), nrho, cmtmpb_omp(1,1,ithread), nrho)
                                             call zmat_zcoomm('r', 'n', 'n', nrho, -ctmp1, cval_zgfr(1,idraw),               &
                                                              irow_zgf(1,idraw), icol_zgf(1,idraw), nnz_zgf(idraw),          &
                                                              cmtmp6_omp(1,1,ithread), nrho, cunity,                         &
                                                              cmtmp8_omp(1,1,ithread), nrho, cval_tmp1_omp(1,ithread),       &
                                                              irow_tmp1_omp(1,ithread), icol_tmp1_omp(1,ithread),            &
                                                              cmtmpa_omp(1,1,ithread), nrho, cmtmpb_omp(1,1,ithread), nrho)
                                             cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp8_omp(1:nrho,1:nrho,ithread) 
                                         end if
                                      end do
                                      dtmp1 = dble(ifactfront * ifactf)         ! fl * fact
                                      dtmp2 = dtmp1 * dble((-1)**(itier-1))     ! fr * fact
                                      dtmp1 = dtmp1 * dtmp3 * dtmp4
                                      dtmp2 = dtmp2 * dtmp3 * dtmp4
                                      ctmp1 = dtmp1 * zfb2a
                                      ctmp2 = dtmp2 * zrb2a
                                      ctmp3 = dtmp1 * zfa2b
                                      ctmp4 = dtmp2 * zrb2a
                                      cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) +      &
                                                                          cmtmp7_omp(1:nrho,1:nrho,ithread) * ctmp1
                                      cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) +      &
                                                                          cmtmp7_omp(1:nrho,1:nrho,ithread) * ctmp2
                                      cmtmp1_omp(1:nrho,1:nrho,ithread) = cmtmp1_omp(1:nrho,1:nrho,ithread) +      &
                                                                          cmtmp6_omp(1:nrho,1:nrho,ithread) * ctmp3
                                      cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) +      &
                                                                          cmtmp6_omp(1:nrho,1:nrho,ithread) * ctmp4
                                   end do
                                   !
                                   cycle loop_draw0
                               end if
                               !
                               dtmp1 = dble(ifactfront)
                               dtmp2 = dble(ifactrear )
                               ctmp1 = dcmplx(dtmp1*dtmp2, 0.d0)      ! fl*gr = fr*gl
                               dtmp1 = dble( (-1)**(itier-1) )
                               ctmp2 = ctmp1 * dtmp1                  ! fl*gl = fr*gr 
                               !
                               ! fl*gl*zgefnl*rcgtmp - fl*gr*rcgtmp*zgefnr
                               transb = 'n'
                               if (itype .ne. 0) transb = 'c'
                               if (lcop_ad) then
                                   !call zlmat_zcoomv4(transb, nrho, ctmp2*zfa2b, nnz_lvl(idraw),                           &
                                   !                   irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                   !                   rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                   !                   czero, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                   !                   nlvrow, nlvcol)
                                   !call zlmat_zcoomv4(transb, nrho, -ctmp1*zra2b, nnz_lvr(idraw),                          &
                                   !                   irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                   !                   rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                   !                   cunity, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                   !                   nlvrow, nlvcol)
                                   !call zlmat_zcoomv4(transb, nrho, ctmp2*zfb2a, nnz_lvl(idraw),                           &
                                   !                   irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                   !                   rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                   !                   czero, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                   !                   nlvrow, nlvcol)
                                   !call zlmat_zcoomv4(transb, nrho, -ctmp1*zrb2a, nnz_lvr(idraw),                          &
                                   !                   irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                   !                   rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                   !                   cunity, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                   !                   nlvrow, nlvcol)
                                   call zlmat_zcoomv5(transb, nrho, ctmp2*zfa2b, nnz_lvl(idraw), irow_lvl(1,idraw), icol_lvl(1,idraw), &
                                                      cval_lvl(1,idraw), rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                      czero, cmtmp7_omp(1,1,ithread), nrho,                                            &
                                                      nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                      &
                                                      ind_lvl(1,idraw), nnc_lvl(1,idraw))
                                   call zlmat_zcoomv5(transb, nrho, -ctmp1*zra2b, nnz_lvr(idraw), irow_lvr(1,idraw), icol_lvr(1,idraw), &
                                                      cval_lvr(1,idraw), rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),             &
                                                      cunity, cmtmp7_omp(1,1,ithread), nrho,                                            &
                                                      nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                       & 
                                                      ind_lvr(1,idraw), nnc_lvr(1,idraw))
                                   call zlmat_zcoomv5(transb, nrho, ctmp2*zfb2a, nnz_lvl(idraw), irow_lvl(1,idraw), icol_lvl(1,idraw), &
                                                      cval_lvl(1,idraw), rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                      czero, cmtmp6_omp(1,1,ithread), nrho,                                            &
                                                      nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                      &
                                                      ind_lvl(1,idraw), nnc_lvl(1,idraw))
                                   call zlmat_zcoomv5(transb, nrho, -ctmp1*zrb2a, nnz_lvr(idraw), irow_lvr(1,idraw), icol_lvr(1,idraw), &
                                                      cval_lvr(1,idraw), rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),             &
                                                      cunity, cmtmp6_omp(1,1,ithread), nrho,                                            &
                                                      nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                       & 
                                                      ind_lvr(1,idraw), nnc_lvr(1,idraw))
                               else
                                   call zmat_zcoomm4('l', 'n', transb, ctmp2*zfa2b, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                     irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                     rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                     czero, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                     cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                   call zmat_zcoomm4('r', 'n', transb, -ctmp1*zra2b, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                     irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                     rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                     cunity, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                     cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                   call zmat_zcoomm4('l', 'n', transb, ctmp2*zfb2a, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                     irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                     rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                     czero, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                     cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                   call zmat_zcoomm4('r', 'n', transb, -ctmp1*zrb2a, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                     irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                     rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                     cunity, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                     cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                               end if
!
                               ctmp1 = degen * zfb2a
                               ctmp2 = degen * dtmp1 * zrb2a
                               ctmp3 = degen * zfa2b
                               ctmp4 = degen * dtmp1 * zra2b
                               do nj=1,nrho
                                  do ni=1,nrho
                                     cmtmp4_omp(ni,nj,ithread) = cmtmp4_omp(ni,nj,ithread) +  &
                                                                 cmtmp7_omp(ni,nj,ithread) * ctmp1
                                     cmtmp5_omp(ni,nj,ithread) = cmtmp5_omp(ni,nj,ithread) +  &
                                                                 cmtmp7_omp(ni,nj,ithread) * ctmp2
                                     cmtmp1_omp(ni,nj,ithread) = cmtmp1_omp(ni,nj,ithread) +  &
                                                                 cmtmp6_omp(ni,nj,ithread) * ctmp3
                                     cmtmp2_omp(ni,nj,ithread) = cmtmp2_omp(ni,nj,ithread) +  &
                                                                 cmtmp6_omp(ni,nj,ithread) * ctmp4
                                  end do
                               end do
                               !
                               if (idegen_fast .eq. 1) then
                                   do mi=2,ndegen,1
                                      lnj        = indexcoef(lnl)
                                      itype      = itypecoef(lnl)
                                      ifactfront = ifactfcoef(lnl)
                                      ifactrear  = ifactrcoef(lnl)
                                      lnl        = lnl + 1
                                      lnq        = indcf_spa(lnj) 
                                      nnz2       = nnzcf_spa(lnj)
                                      if (nnz2 .eq. 0) cycle
                                      !
                                      ctmp1 = dcmplx(dble(ifactfront*ifactrear), 0.d0)      ! fl*gr = fr*gl
                                      dtmp1 = dble( (-1)**(itier-1) )
                                      ctmp2 = ctmp1 * dtmp1                                 ! fl*gl = fr*gr 
                                      idraw = nvec_tmp1_omp(mi,ithread)
                                      !
                                      transb = 'n'
                                      if (itype .ne. 0) transb = 'c'
                                      if (lcop_ad) then
                                          !call zlmat_zcoomv4(transb, nrho, ctmp2*zfa2b, nnz_lvl(idraw),                           &
                                          !                   irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                          !                   rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                          !                   czero, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                          !                   nlvrow, nlvcol)
                                          !call zlmat_zcoomv4(transb, nrho, -ctmp1*zra2b, nnz_lvr(idraw),                          &
                                          !                   irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                          !                   rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                          !                   cunity, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                          !                   nlvrow, nlvcol)
                                          !call zlmat_zcoomv4(transb, nrho, ctmp2*zfb2a, nnz_lvl(idraw),                           &
                                          !                   irow_lvl(1,idraw), icol_lvl(1,idraw), cval_lvl(1,idraw),             &
                                          !                   rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                          !                   czero, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                          !                   nlvrow, nlvcol)
                                          !call zlmat_zcoomv4(transb, nrho, -ctmp1*zrb2a, nnz_lvr(idraw),                          &
                                          !                   irow_lvr(1,idraw), icol_lvr(1,idraw), cval_lvr(1,idraw),             &
                                          !                   rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                   &
                                          !                   cunity, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho, &
                                          !                   nlvrow, nlvcol)
                                          call zlmat_zcoomv5(transb, nrho, ctmp2*zfa2b, nnz_lvl(idraw), irow_lvl(1,idraw), icol_lvl(1,idraw), &
                                                             cval_lvl(1,idraw), rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                             czero, cmtmp7_omp(1,1,ithread), nrho,                                            &
                                                             nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                      &
                                                             ind_lvl(1,idraw), nnc_lvl(1,idraw))
                                          call zlmat_zcoomv5(transb, nrho, -ctmp1*zra2b, nnz_lvr(idraw), irow_lvr(1,idraw), icol_lvr(1,idraw), &
                                                             cval_lvr(1,idraw), rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),             &
                                                             cunity, cmtmp7_omp(1,1,ithread), nrho,                                            &
                                                             nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                       & 
                                                             ind_lvr(1,idraw), nnc_lvr(1,idraw))
                                          call zlmat_zcoomv5(transb, nrho, ctmp2*zfb2a, nnz_lvl(idraw), irow_lvl(1,idraw), icol_lvl(1,idraw), &
                                                             cval_lvl(1,idraw), rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),            &
                                                             czero, cmtmp6_omp(1,1,ithread), nrho,                                            &
                                                             nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                      &
                                                             ind_lvl(1,idraw), nnc_lvl(1,idraw))
                                          call zlmat_zcoomv5(transb, nrho, -ctmp1*zrb2a, nnz_lvr(idraw), irow_lvr(1,idraw), icol_lvr(1,idraw), &
                                                             cval_lvr(1,idraw), rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),             &
                                                             cunity, cmtmp6_omp(1,1,ithread), nrho,                                            &
                                                             nlvrow, nlvcol, nnrc2v, nrho, ncrc2v, nrho,                                       & 
                                                             ind_lvr(1,idraw), nnc_lvr(1,idraw))
                                      else
                                          call zmat_zcoomm4('l', 'n', transb, ctmp2*zfa2b, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                            irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                            rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                            czero, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                            cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                          call zmat_zcoomm4('r', 'n', transb, -ctmp1*zra2b, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                            irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                            rcga(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                            cunity, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                            cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                          call zmat_zcoomm4('l', 'n', transb, ctmp2*zfb2a, cval_zgfl(1,idraw), nnz_zgf(idraw),  &
                                                            irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                            rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                            czero, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                   &
                                                            cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                          call zmat_zcoomm4('r', 'n', transb, -ctmp1*zrb2a, cval_zgfr(1,idraw), nnz_zgf(idraw), &
                                                            irow_zgf(1,idraw), icol_zgf(1,idraw),                               &
                                                            rcgb(lnq), nnz2, irowcf_spa(lnq), icolcf_spa(lnq),                  &
                                                            cunity, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,                  &
                                                            cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                                      end if
!
                                      ctmp1 = degen * zfb2a
                                      ctmp2 = degen * dtmp1 * zrb2a
                                      ctmp3 = degen * zfa2b
                                      ctmp4 = degen * dtmp1 * zra2b
                                      do nj=1,nrho
                                         do ni=1,nrho
                                            cmtmp4_omp(ni,nj,ithread) = cmtmp4_omp(ni,nj,ithread) +  &
                                                                        cmtmp7_omp(ni,nj,ithread) * ctmp1
                                            cmtmp5_omp(ni,nj,ithread) = cmtmp5_omp(ni,nj,ithread) +  &
                                                                        cmtmp7_omp(ni,nj,ithread) * ctmp2
                                            cmtmp1_omp(ni,nj,ithread) = cmtmp1_omp(ni,nj,ithread) +  &
                                                                        cmtmp6_omp(ni,nj,ithread) * ctmp3
                                            cmtmp2_omp(ni,nj,ithread) = cmtmp2_omp(ni,nj,ithread) +  &
                                                                        cmtmp6_omp(ni,nj,ithread) * ctmp4
                                         end do
                                      end do
                                   end do
                               end if
                           end if
!
                        end do loop_draw0
                     end do
                     if (lprams) then
                         call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                                               rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin,   eye , cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                               rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp1_omp(1,1,ithread), nrho, cunity,  &
                                               rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin,   eye , cmtmp2_omp(1,1,ithread), nrho, cunity,  &
                                               rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                     else
                         call zamsmm2_coo('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('r', transa, 'n', iorbs1, ispin,   eye , cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp1_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('r', transa, 'n', iorbs1, ispin,   eye , cmtmp2_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                     end if
                  end do
               end do
            end do
!
        else if (itier .gt. ntier .and. itier .lt. ntier_ad) then
            if (lad_fast) then
                ntmp1 = ncor_fast
            else
                ntmp1 = ncor_slow
            end if
            do isgn=1,nsgn
               transa = 'n'
               if (isgn .ne. 1) transa = 'c'
               do ispin=1,nspin
                  do iorbs1=1,norbs
                     call look4sysopr(isgn, iorbs1, ispin, isopr)
                     cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
                     cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
                     ! slow modes (or fast modes)
                     do ialf=1,nalf
                        iopr = (isopr - 1) * nalf + ialf
                        if (lscreen .and. jomit(iopr) .eq. 1) cycle
                        !
                        do imats=1,ntmp1
                           lnj        = indexcoef(lnl)
                           itype      = itypecoef(lnl)
                           ifactfront = ifactfcoef(lnl)
                           ifactrear  = ifactrcoef(lnl)
                           lnl        = lnl + 1
                           lnq  = indcf_spa(lnj)                   
                           nnz2 = nnzcf_spa(lnj)
                           if (nnz2 .eq. 0) cycle
                           !
                           dtmp1 = dble(ifactfront)
                           dtmp2 = dble(ifactrear )
                           ctmp1 = dcmplx(dtmp1, 0.d0)
                           ctmp2 = dcmplx(dtmp2, 0.d0)
                           transb = 'n'
                           if (itype .ne. 0) transb = 'c'
                           call zmat_accu_coo2dns_pair(transb, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcgb(lnq), &
                                                       ctmp1*zfb2a, cmtmp4_omp(1,1,ithread), nrho,                &
                                                       ctmp2*zrb2a, cmtmp5_omp(1,1,ithread), nrho)
                           call zmat_accu_coo2dns_pair(transb, nnz2, irowcf_spa(lnq), icolcf_spa(lnq), rcga(lnq), &
                                                       ctmp1*zfa2b, cmtmp1_omp(1,1,ithread), nrho,                &
                                                       ctmp2*zra2b, cmtmp2_omp(1,1,ithread), nrho)
                        end do 
                     end do
                     !if (.not. lad_fast) then
                         ! fast modes
                     !if (nnz_sop(isopr) .gt. 0) then
                     !if (nnz_lvlsop(isopr) .gt. 0 .or. nnz_lvrsop(isopr) .gt. 0) then
                         ! fast modes or slow modes, depending on lad_fast
                         ! zgefnl_sop*rcgtmp - rcgtmp*zgefnr_sop
                     if (lcop_ad) then
                         call zlmat_zcoomv4('n', nrho, zfa2b, nnz_lvlsop(isopr), irow_lvlsop(1,isopr), icol_lvlsop(1,isopr),   &
                                            cval_lvlsop(1,isopr), rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                            czero, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho,               &
                                            nlvrow, nlvcol)
                         call zlmat_zcoomv4('n', nrho, -zra2b, nnz_lvrsop(isopr), irow_lvrsop(1,isopr), icol_lvrsop(1,isopr),  &
                                            cval_lvrsop(1,isopr), rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                            cunity, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho,              &
                                            nlvrow, nlvcol)
                         call zlmat_zcoomv4('n', nrho, zfb2a, nnz_lvlsop(isopr), irow_lvlsop(1,isopr), icol_lvlsop(1,isopr),   &
                                            cval_lvlsop(1,isopr), rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                            czero, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho,               &
                                            nlvrow, nlvcol)
                         call zlmat_zcoomv4('n', nrho, -zrb2a, nnz_lvrsop(isopr), irow_lvrsop(1,isopr), icol_lvrsop(1,isopr),  &
                                            cval_lvrsop(1,isopr), rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                            cunity, cmtmp6_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho,              &
                                            nlvrow, nlvcol)
                     else
                         call zmat_zcoomm4('l', 'n', 'n', zfa2b, cval_sopl(1,isopr), nnz_sop(isopr),    &
                                           irow_sop(1,isopr), icol_sop(1,isopr),                        &
                                           rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                           czero, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,            &
                                           cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                         call zmat_zcoomm4('r', 'n', 'n', -zra2b, cval_sopr(1,isopr), nnz_sop(isopr),   &
                                           irow_sop(1,isopr), icol_sop(1,isopr),                        &
                                           rcga(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                           cunity, cmtmp7_omp(1,1,ithread), nrho, nrho, nrho,           &
                                           cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                         call zmat_zcoomm4('l', 'n', 'n', zfb2a, cval_sopl(1,isopr), nnz_sop(isopr),    &
                                           irow_sop(1,isopr), icol_sop(1,isopr),                        &
                                           rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                           czero, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,            &
                                           cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                         call zmat_zcoomm4('r', 'n', 'n', -zrb2a, cval_sopr(1,isopr), nnz_sop(isopr),   &
                                           irow_sop(1,isopr), icol_sop(1,isopr),                        &
                                           rcgb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp),            &
                                           cunity, cmtmp6_omp(1,1,ithread), nrho, nrho, nrho,           &
                                           cval_tmp1_omp(1,ithread), cval_tmp2_omp(1,ithread))  
                     end if
                     cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) +  &
                                                         cmtmp7_omp(1:nrho,1:nrho,ithread) * zfb2a
                     cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) +  &
                                                         cmtmp7_omp(1:nrho,1:nrho,ithread) * zrb2a
                     cmtmp1_omp(1:nrho,1:nrho,ithread) = cmtmp1_omp(1:nrho,1:nrho,ithread) +  &
                                                         cmtmp6_omp(1:nrho,1:nrho,ithread) * zfa2b
                     cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) +  &
                                                         cmtmp6_omp(1:nrho,1:nrho,ithread) * zra2b
!
                     if (lprams) then
                         call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                                               rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin,   eye , cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                               rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp1_omp(1,1,ithread), nrho, cunity,  &
                                               rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo_proj('r', transa, 'n', iorbs1, ispin,   eye , cmtmp2_omp(1,1,ithread), nrho, cunity,  &
                                               rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                     else
                         call zamsmm2_coo('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('r', transa, 'n', iorbs1, ispin,   eye , cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('l', transa, 'n', iorbs1, ispin,  -eye , cmtmp1_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                         call zamsmm2_coo('r', transa, 'n', iorbs1, ispin,   eye , cmtmp2_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                     end if
                  end do
               end do
            end do
!
        end if
    end if
!
! Truncation for the derivative of terminal-tier ADOs
!
    if (ltrun_der .and. itier .eq. ntier) then  ! now only for freq=0 case
       do isgn=1,nsgn
          transc = 'n'
          if (isgn .ne. 1) transc = 'c'
!
          do ispin=1,nspin
             do iorbs1=1,norbs
                cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
                cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
                cmtmpa_omp(1:nrho,1:nrho,ithread) = czero
                cmtmpb_omp(1:nrho,1:nrho,ithread) = czero
                do ialf=1,nalf
                   call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                   if (lscreen .and. jomit(iopr) .eq. 1) cycle
                   !
                   do imats=1,ncor
                      ifactf = ifactfcoef(lnl)
                      ifactr = ifactrcoef(lnl)
                      idraw  = indexcoef(lnl)
                      lnl = lnl + 1
                      cgamma_nj = cgamma_nk + cgamma(idraw)
                      if (igroundsteady .eq. 1) then
                          cgamma_nj = cgamma_nj + eye * dipm(idraw) * eleadinfty(ilead(idraw), mspin(idraw))
                      end if
                      call add2list(index_nk(1,ithread), mdim, ntier-1, index_nj(1,ithread), mdim, idraw, mj)
                      dtmp1 = dble(ifactf)
                      dtmp2 = dble(ifactr)
                      if (lscale) then
                         dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                         dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                      end if     
!
! (i*L - gamma_j)^(-1) = X * (diagonal part) * X^(-1) 
! Note: Use consecutive matrix-vector multiplication rather than matrix-matrix multiplication to save
!       computational time
!
                      cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
                      cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
                      do mi=1,ntier
                         lnj        = indexcoef(lnl)
                         itype      = itypecoef(lnl)
                         ifactfront = ifactfcoef(lnl)
                         ifactrear  = ifactrcoef(lnl)
                         lnl = lnl + 1
                         nn  = index_nj(mi,ithread)
                         lnq  = indcf_spa(lnj)
                         nnz2 = nnzcf_spa(lnj)
                         transa = 'n'
                         if (mpm(nn) .eq. 1) transa = 'c'
                         transb = 'n'
                         if (itype .ne. 0) transb = 'c'
                         if (.not. offcor) then
                            ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) *   eye 
                            ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * (-eye)
                            if (lprams) then
                                call zamsmm_coo_proj('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcga(lnq), nnz2,        &
                                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                     cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo_proj('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcga(lnq), nnz2,        &
                                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                     cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo_proj('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcgb(lnq), nnz2,       &
                                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                     cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo_proj('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcgb(lnq), nnz2,       &
                                                     irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                     cmtmp7_omp(1,1,ithread), nrho)
                            else
                                call zamsmm_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcga(lnq), nnz2,        &
                                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcga(lnq), nnz2,        &
                                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcgb(lnq), nnz2,       &
                                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                cmtmp7_omp(1,1,ithread), nrho)
                                call zamsmm_coo('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcgb(lnq), nnz2,       &
                                                irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                cmtmp7_omp(1,1,ithread), nrho)
                            end if
                         else
                            do iorbs2=1,norbs
                               ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
                               ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) *   eye
                               ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * (-eye)
                               if (lprams) then
                                   call zamsmm_coo_proj('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcga(lnq), nnz2,           &
                                                        irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                        cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo_proj('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcga(lnq), nnz2,           &
                                                        irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                        cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo_proj('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcgb(lnq), nnz2,          &
                                                        irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                        cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo_proj('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcgb(lnq), nnz2,          &
                                                        irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                        cmtmp7_omp(1,1,ithread), nrho)
                               else
                                   call zamsmm_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcga(lnq), nnz2,           &
                                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                   cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcga(lnq), nnz2,           &
                                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp2_omp(1,1,ithread), nrho,  &
                                                   cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcgb(lnq), nnz2,          &
                                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                   cmtmp7_omp(1,1,ithread), nrho)
                                   call zamsmm_coo('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcgb(lnq), nnz2,          &
                                                   irowcf_spa(lnq), icolcf_spa(lnq), cunity, cmtmp1_omp(1,1,ithread), nrho,  &
                                                   cmtmp7_omp(1,1,ithread), nrho)
                               end if
                            end do
                         end if
                      end do
!
                      call zgemv('c', nrho2, nrho2, -eye, zmatder_cf(1,1,isgnw), nrho2d, cmtmp2_omp(1,1,ithread), 1, &
                                 czero, cmtmp7_omp(1,1,ithread), 1) 
                      call zgemv('c', nrho2, nrho2, -eye, zmatder_cf(nrho2+1,1,isgnw), nrho2d, cmtmp1_omp(1,1,ithread), 1, &
                                 cunity, cmtmp7_omp(1,1,ithread), 1)
                      call zgemv('c', nrho2, nrho2, -eye, zmatder_cf(1,nrho2+1,isgnw), nrho2d, cmtmp2_omp(1,1,ithread), 1, &
                                 czero, cmtmp8_omp(1,1,ithread), 1) 
                      call zgemv('c', nrho2, nrho2, -eye, zmatder_cf(nrho2+1,nrho2+1,isgnw), nrho2d, cmtmp1_omp(1,1,ithread), 1, &
                                 cunity, cmtmp8_omp(1,1,ithread), 1)
                      mi = 1
                      do j=1,nrho
                         do i=1,nrho
                            cmtmp7_omp(j,i,ithread) = cmtmp7_omp(j,i,ithread) / (dvalder_cf(mi,isgnw) + eye * cgamma_nj)
                            cmtmp8_omp(j,i,ithread) = cmtmp8_omp(j,i,ithread) / (dvalder_cf(mi+nrho2,isgnw) + eye * cgamma_nj)
                            mi = mi + 1
                         end do
                      end do
                      call zgemv('n', nrho2, nrho2, cunity, zmatder_cf(1,1,isgnw), nrho2d, cmtmp7_omp(1,1,ithread), 1, &
                                 czero, cmtmp2_omp(1,1,ithread), 1)
                      call zgemv('n', nrho2, nrho2, cunity, zmatder_cf(1,nrho2+1,isgnw), nrho2d, cmtmp8_omp(1,1,ithread), 1, &
                                 cunity, cmtmp2_omp(1,1,ithread), 1)
                      call zgemv('n', nrho2, nrho2, cunity, zmatder_cf(nrho2+1,1,isgnw), nrho2d, cmtmp7_omp(1,1,ithread), 1, &
                                 czero, cmtmp1_omp(1,1,ithread), 1)
                      call zgemv('n', nrho2, nrho2, cunity, zmatder_cf(nrho2+1,nrho2+1,isgnw), nrho2d, cmtmp8_omp(1,1,ithread), 1, &
                                 cunity, cmtmp1_omp(1,1,ithread), 1)
                      cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp1
                      cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp6_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp2
                      cmtmpa_omp(1:nrho,1:nrho,ithread) = cmtmpa_omp(1:nrho,1:nrho,ithread) + cmtmp1_omp(1:nrho,1:nrho,ithread) * dtmp1
                      cmtmpb_omp(1:nrho,1:nrho,ithread) = cmtmpb_omp(1:nrho,1:nrho,ithread) + cmtmp1_omp(1:nrho,1:nrho,ithread) * dtmp2
                   end do
                end do
                if (lprams) then
                    call zamsmm2_coo_proj('l', transc, 'n', iorbs1, ispin, -cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo_proj('r', transc, 'n', iorbs1, ispin, -cunity, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                                          rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo_proj('l', transc, 'n', iorbs1, ispin,  cunity, cmtmpa_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo_proj('r', transc, 'n', iorbs1, ispin,  cunity, cmtmpb_omp(1,1,ithread), nrho, cunity,  &
                                          rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                else
                    call zamsmm2_coo('l', transc, 'n', iorbs1, ispin, -cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                     rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo('r', transc, 'n', iorbs1, ispin, -cunity, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                                     rhsa(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo('l', transc, 'n', iorbs1, ispin,  cunity, cmtmpa_omp(1,1,ithread), nrho, cunity,  &
                                     rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                    call zamsmm2_coo('r', transc, 'n', iorbs1, ispin,  cunity, cmtmpb_omp(1,1,ithread), nrho, cunity,  &
                                     rhsb(lnp), nnz, irowcf_spa(lnp), icolcf_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
                end if
             end do
          end do
       end do
    end if

!
  end do
!$omp end do nowait
!$omp end parallel
!
  call cpu_time(cpu2)
!  write(6,100)itier, nlast(itier), cpu2 - cpu1
!  call flush(6)
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
end do 
!
606 continue
!
deallocate(irow_tmp1, icol_tmp1, cval_tmp1, STAT=istat)
deallocate(irow_heff, icol_heff, cval_heff, STAT=istat)
deallocate(irow_tmp1_omp, icol_tmp1_omp, STAT=istat)
deallocate(irow_tmp2_omp, icol_tmp2_omp, STAT=istat)
deallocate(cval_tmp1_omp, cval_tmp2_omp, STAT=istat)
!
call cpu_time(cpu4)
!write(6,*)'calc_ax_cf_spa_omp: cpu time = ', cpu4 - cpu3
!call flush(6)
!
return
end subroutine calc_ax_cf_spa_omp
