subroutine calc_ax_jacobi_spa_omp(ldim, rcgtmp, rcgrhs)
!
! purpose : 
!           upate rho with Jacobi iteration process
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
complex*16, intent(in)    :: rcgtmp(*)
complex*16, intent(out)   :: rcgrhs(*)
integer                   :: itier, itype, dim1, ifactfront, ifactrear, istat
integer                   :: index_nk(MAXTIER, maxprocs)
integer*8                 :: lni, lnj, lnk, lnl, lnm, lnp, lnq
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs
integer                   :: iopr
integer                   :: nrho2, mi, mj
real*8                    :: cpu1, cpu2
real*8                    :: dtmp1, dtmp2
complex*16                :: cgamma_nk
complex*16                :: ctmp1, ctmp2
character*1               :: transa, transb, transc
integer                   :: nnz, nnz2
complex*16, allocatable   :: cval_tmp1(:)
complex*16, allocatable   :: cval_tmp1_omp(:,:)
integer                   :: ithread, icore
!
integer    :: itimes
data itimes /0/
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'calc_ax_jacobi_spa_omp: error entry, lsparse = ', lsparse
   stop
end if
if (ldim .ne. lunk_spa) then
   write(6,*)
   write(6,*)'calc_ax_jacobi_spa_omp: error dimension, ldim != lunk_spa ', ldim, lunk_spa
   stop
end if
if (itimes .eq. 0) then
   write(6,*)
   write(6,*)'calc_ax_jacobi_spa_omp: first entry'
   write(6,*)'calc_ax_jacobi_spa_omp: nnzjac = ', nnzjac
   write(6,*)'calc_ax_jacobi_spa_omp:   ejac = ', ejac*hbar
   call flush(6)
end if
!
dim1   = nrho**2
nrho2  = nrho**2
itimes = itimes + 1
!
allocate(cval_tmp1(dim1), STAT=istat)
allocate(cval_tmp1_omp(dim1,nthreads_omp), STAT=istat)
cval_tmp1(1:dim1) = czero
cval_tmp1_omp(1:dim1,1:nthreads_omp) = czero
!
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
lni = ind_spa(1)
nnz = nnz_spa(1)
call zmat_coo2dns(nnz, irow_spa(lni), icol_spa(lni), rcgtmp(lni), cmtmp1, nrho, nrho, nrho)
cmtmp1 = ejac * cmtmp1
cmtmp2 = czero
!
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4(1:nrho,1:nrho) = czero
    do ialf=1,nalf
       call look4operator(1, iorbs1, ispin, ialf, iopr)
       if (lscreen .and. jomit(iopr) .eq. 1) cycle
!
!$omp parallel default(shared) private(imats,dtmp1,nn,ithread,lni,lnj,nnz)
      ithread = omp_get_thread_num() + 1
      cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
!$omp do schedule(dynamic) 
      do imats=1,ncor
        if (lscale) then
           dtmp1 = dbsqrt(iorbs1,ispin,imats,ialf,1)
        else
           dtmp1 = 1.d0
        end if
        call look4drawer(1, ialf, iorbs1, ispin, imats, nn)
        lni = nfirst(2) - 1 + indexdraw(nn)
        lnj = ind_spa(lni)
        nnz = nnz_spa(lni)
        call zmat_coo2dns(nnz, irow_spa(lnj), icol_spa(lnj), rcgtmp(lnj),  &
                          cmtmp1_omp(1,1,ithread), nrho, nrho, nrho)
        cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) + dtmp1 * &
                                            cmtmp1_omp(1:nrho,1:nrho,ithread)
      end do
!$omp end do
!$omp end parallel
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
      end do
    end do
    call zamsmm_spa_omp('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp2, nrho, &
                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
    call zamsmm_spa_omp('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho, &
                        dmtmp1, nrho, dmtmp2, nrho, dmtmp3, nrho, dmtmp4, nrho)
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      cmtmp1(ni,nj) = cmtmp1(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
   end do
end do
cmtmp3 = cmtmp1
!
call zgemv('c', nrho2, nnzjac, cunity, zmatjac, nrho2, cmtmp1, 1, czero, cmtmp2, 1)
ntmp3 = 0
loop0: do mj=1,nrho
   do mi=1,nrho
      ntmp3 = ntmp3 + 1
      if (ntmp3 .gt. nnzjac) exit loop0
      ctmp1 = cunity / (eye * dvaljac(ntmp3) + ejac) - cunity / ejac
      cmtmp2(mi,mj) = cmtmp2(mi,mj) * ctmp1
   end do
end do loop0
call zgemv('n', nrho2, nnzjac, cunity, zmatjac, nrho2, cmtmp2, 1, czero, cmtmp1, 1)
cmtmp1 = cmtmp1 + cmtmp3 / ejac
call zmat_extract_coo(nnz, irow_spa(lni), icol_spa(lni), cval_tmp1, cmtmp1, nrho, nrho, nrho)
call zmat_add_coo('n', nnz, cunity, cval_tmp1, cunity, rcgrhs(lni))
!
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn) &
!$omp private(ntmp3,isgn,ispin,iorbs1,iorbs2,ialf,imats,dtmp1,dtmp2,i,j,ithread) &
!$omp private(lnp,lnq,nnz,nnz2,transa,transb,transc,ctmp1,ctmp2,mi,mj) &
!$omp private(iopr,lnm)
  ithread = omp_get_thread_num() + 1
!$omp do schedule(dynamic) 
  do lni=nfirst(itier), nlast(itier)
     if (nnz_spa(lni) .eq. 0) cycle 
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
     lnl = index_coef_ref(lni)
     if (lsimple) lnm = ioprref(lni)
!
     cgamma_nk = czero
     do ni=1,itier-1 
        index_nk(ni,ithread) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        cgamma_nk = cgamma_nk + cgamma(index_nk(ni,ithread))
        if (igroundsteady .ne. 0) then
           cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                         &
                       eleadinfty(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
        end if
     end do
!
! same itier contribution (N)
!
     lnp = ind_spa(lni)
     nnz = nnz_spa(lni)
     call zmat_add_coo('n', nnz, dcmplx(ejac,0.d0), rcgtmp(lnp), czero, cval_tmp1_omp(1,ithread))
!
! nearest lower tier contribution (N - 1)
!
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
           lnq = ind_spa(lnj)
           transa = 'n'
           if (itype .ne. 0) then
               write(6,*)
               write(6,*)'calc_ax_jacobi_spa_omp: unexpected itype ', itype
               stop
           end if
           call zmat_add_coo(transa, nnz, ctmp1, rcgtmp(lnq), cunity, cval_tmp1_omp(1,ithread))
           cycle
       end if
!------
!cycle ! -
!------
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
         call zamsmm1_coo('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgtmp(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                          irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
         call zamsmm1_coo('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgtmp(lnq), nnz2,  &
                          irow_spa(lnq), icol_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                          irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
      else
        do iorbs2=1,norbs
           ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
           ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
           ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
           call zamsmm1_coo('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgtmp(lnq), nnz2,     &
                            irow_spa(lnq), icol_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                            irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)                          
           call zamsmm1_coo('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgtmp(lnq), nnz2,     &
                            irow_spa(lnq), icol_spa(lnq), cunity, cval_tmp1_omp(1,ithread), nnz,  &
                            irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
!          
        end do
      end if
    end do
    call zmat_coo2dns(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1_omp(1,ithread),   &
                      cmtmp1_omp(1,1,ithread), nrho, nrho, nrho)
    cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) +  &
                                        cmtmp1_omp(1:nrho,1:nrho,ithread)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      do isgn=1,nsgn
        transa = 'n'
        if (isgn .ne. 1) transa = 'c'
        do ispin=1,nspin
!
          orbs0: do iorbs1=1,norbs
             if (lsimple .and. .not. lwalf) then
                iopr = ioprindex(lnm)
                lnm = lnm + 1
                if (iopr .lt. 0) cycle orbs0
             end if
!
             cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
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
                   lnj        = indexcoef(lnl)
                   itype      = itypecoef(lnl)
                   ifactfront = ifactfcoef(lnl)
                   ifactrear  = ifactrcoef(lnl)
                   lnl        = lnl + 1
!------
!cycle ! +
!------
! -eye * dtmp1 * op(rhotmp) at right, eye * dtmp2 * op(rhotmp) at left                   
                   dtmp1 = dble(ifactfront)
                   dtmp2 = dble(ifactrear )
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
                                          nrho, cmtmp4_omp(1,1,ithread), nrho)
                   call zmat_accu_coo2dns(transb, ctmp2, nnz2, irow_spa(lnq), icol_spa(lnq), rcgtmp(lnq),  &
                                          nrho, cmtmp5_omp(1,1,ithread), nrho)
                end do
             end do alf0
             call zamsmm2_coo('l', transa, 'n', iorbs1, ispin, cunity, cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                              cval_tmp1_omp(1,ithread), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm2_coo('r', transa, 'n', iorbs1, ispin, cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                              cval_tmp1_omp(1,ithread), nnz, irow_spa(lnp), icol_spa(lnp), cmtmp8_omp(1,1,ithread), nrho)
          end do orbs0
        end do
      end do
    end if
!
    call zmat_coo2dns(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1_omp(1,ithread),  &
                      cmtmp2_omp(1,1,ithread), nrho, nrho, nrho)
    cmtmp3_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread)
    call zgemv('c', nrho2, nnzjac, cunity, zmatjac, nrho2, cmtmp2_omp(1,1,ithread), 1,  &
                                                    czero, cmtmp1_omp(1,1,ithread), 1)
    ntmp3 = 0
    loop1: do mj=1,nrho
       do mi=1,nrho
          ntmp3 = ntmp3 + 1
          if (ntmp3 .gt. nnzjac) exit loop1
          if (cdabs(ejac - cgamma_nk) .gt. dpico) then
              ctmp1 = cunity / (eye * dvaljac(ntmp3) + ejac - cgamma_nk) -  &
                      cunity / (ejac - cgamma_nk)
          else
              ctmp1 = cunity / (eye * dvaljac(ntmp3))
          end if
          cmtmp1_omp(mi,mj,ithread) = cmtmp1_omp(mi,mj,ithread) * ctmp1
       end do
    end do loop1
    call zgemv('n', nrho2, nnzjac, cunity, zmatjac, nrho2, cmtmp1_omp(1,1,ithread), 1,  &
                                                    czero, cmtmp2_omp(1,1,ithread), 1)
    if (cdabs(ejac - cgamma_nk) .gt. dpico) then
        cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) +         &
                                            cmtmp3_omp(1:nrho,1:nrho,ithread) / (ejac - cgamma_nk)
    end if
    call zmat_extract_coo(nnz, irow_spa(lnp), icol_spa(lnp), cval_tmp1_omp(1,ithread),  &
                          cmtmp2_omp(1,1,ithread), nrho, nrho, nrho)
    call zmat_add_coo('n', nnz, cunity, cval_tmp1_omp(1,ithread), cunity, rcgrhs(lnp))
!
  end do
!$omp end do nowait
!$omp end parallel
!
  call cpu_time(cpu2)
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
end do 
!
300 continue
!
deallocate(cval_tmp1, STAT=istat)
deallocate(cval_tmp1_omp, STAT=istat)
!
return
end subroutine calc_ax_jacobi_spa_omp
