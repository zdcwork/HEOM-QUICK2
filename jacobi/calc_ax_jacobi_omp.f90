subroutine calc_ax_jacobi_omp(dim0, rcgtmp, rcgrhs)
!
! purpose : 
!           upate rho with Jacobi iteration process
!    
use matmod
use tmpmatmod
use omp_lib
use matmod_omp
use tmpmatmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer,    intent(in)    :: dim0
complex*16, intent(in)    :: rcgtmp(dim0,dim0,*)
complex*16, intent(out)   :: rcgrhs(dim0,dim0,*)  
integer                   :: itier, itype, dim1, ifactfront, ifactrear
integer                   :: index_nk(MAXTIER, maxprocs)
integer*8                 :: lni, lnj, lnk, lnl, lnm
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, iopr
integer                   :: nrho2, nwork, info, mi, mj, istat
real*8                    :: cpu1, cpu2
real*8                    :: dtmp1, dtmp2
complex*16                :: ctmp1, ctmp2
character*1               :: transa, transb
complex*16                :: cgamma_nk
integer                   :: ithread, icore
!real*8,     allocatable   :: rwork(:), tmpval(:)
!complex*16, allocatable   :: zwork(:)
integer                   :: itimes
data itimes /0/
!
if (itimes .eq. 0) then
   write(6,*)'calc_ax_jacobi_omp: first entry'
   write(6,*)'calc_ax_jacobi_omp: nnzjac = ', nnzjac
   write(6,*)'calc_ax_jacobi_omp:   ejac = ', ejac*hbar
   call flush(6)
end if
!
dim1   = nrho**2
nrho2  = nrho**2
itimes = itimes + 1
!
if (dim0 .ne. nrho) then
  write(6,*)
  write(6,*)' error! dim0 != nrho in calc_ax_jacobi_omp ', dim0, nrho
  stop
end if
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
rcgrhs(1:nrho,1:nrho,1:nunk) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
!
cmtmp1(1:nrho,1:nrho) = ejac * rcgtmp(1:nrho,1:nrho,1)
cmtmp2(1:nrho,1:nrho) = czero
!
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4(1:nrho,1:nrho) = czero
    do ialf=1,nalf
       call look4operator(1, iorbs1, ispin, ialf, iopr)
       if (lscreen .and. jomit(iopr) .eq. 1) cycle
!
!$omp parallel default(shared) private(imats,dtmp1,nn,lni,ithread)
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
        cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) + &
                                            rcgtmp(1:nrho,1:nrho,lni) * dtmp1
      end do
!$omp end do
!$omp end parallel
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
      end do
    end do
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      cmtmp1(ni,nj) = cmtmp1(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
   end do
end do
cmtmp3(1:nrho,1:nrho) = cmtmp1(1:nrho,1:nrho)
call zgemv('c', nrho2, nnzjac, cunity, zmatjac, nrho2, cmtmp1, 1, czero, cmtmp2, 1)
!
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
rcgrhs(1:nrho,1:nrho,1) = rcgrhs(1:nrho,1:nrho,1) + cmtmp1(1:nrho,1:nrho) + &
                          cmtmp3(1:nrho,1:nrho) / ejac
!
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn) &
!$omp private(ntmp3,isgn,ispin,iorbs1,iorbs2,ialf,imats,dtmp1,dtmp2,i,j,ithread) &
!$omp private(transa,transb,ctmp1,ctmp2,mi,mj) &
!$omp private(iopr,lnm)
  ithread = omp_get_thread_num() + 1
!$omp do schedule(dynamic) 
  do lni=nfirst(itier), nlast(itier)
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
          cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                             &
                      eleadinfty(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
       end if
    end do
!
! same itier contribution (N)
!
    cmtmp2_omp(1:nrho,1:nrho,ithread) = ejac * rcgtmp(1:nrho,1:nrho,lni)
!
! nearest lower tier contribution (N - 1)
!
    do ni=1,itier-1
      lnj        = indexcoef(lnl)
      itype      = itypecoef(lnl)
      ifactfront = ifactfcoef(lnl)
      ifactrear  = ifactrcoef(lnl)
      lnl        = lnl + 1
!------
!cycle ! -
!------
      nn = index_nk(ni,ithread)
!
      if (psdfff .and. mmfff(nn) .gt. 1) then
          ctmp1 = dble(ifactfront) * ctheta(mspin(nn),mcor(nn),ilead(nn),mpm(nn))
          if (itype .ne. 0) then
              write(6,*)
              write(6,*)'calc_ax_jacobi_omp: unexpected itype ', itype
              stop
          end if
          cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) + ctmp1 * rcgtmp(1:nrho,1:nrho,lnj)
          cycle
      end if
!
      transa = 'n'
      if (mpm(nn) .eq. 1) transa = 'c'
      transb = 'n'
      if (itype .ne. 0) transb = 'c'
!
      if (.not. offcor) then
         ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
         ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
         call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgtmp(1,1,lnj), nrho, cunity,             &
                          cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgtmp(1,1,lnj), nrho, cunity,             &
                          cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
      else
         do iorbs2=1,norbs
            ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
            ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
            ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
            call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgtmp(1,1,lnj), nrho, cunity,               &
                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
            call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgtmp(1,1,lnj), nrho, cunity,               &
                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         end do
      end if
    end do
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
!
             cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp7_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp8_omp(1:nrho,1:nrho,ithread) = czero
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
!------
!cycle ! +
!------
                 dtmp1 = dble(ifactfront)
                 dtmp2 = dble(ifactrear )
                 if (lscale) then
                    dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if     
                 if (itype .eq. 0) then
                   cmtmp5_omp(1:nrho, 1:nrho, ithread) = cmtmp5_omp(1:nrho, 1:nrho, ithread) +   &
                                                         rcgtmp(1:nrho, 1:nrho, lnj) * dtmp1
                   cmtmp6_omp(1:nrho, 1:nrho, ithread) = cmtmp6_omp(1:nrho, 1:nrho, ithread) +   &
                                                         rcgtmp(1:nrho, 1:nrho, lnj) * dtmp2
                 else 
                   cmtmp7_omp(1:nrho, 1:nrho, ithread) = cmtmp7_omp(1:nrho, 1:nrho, ithread) +   &
                                                         rcgtmp(1:nrho, 1:nrho, lnj) * dtmp1
                   cmtmp8_omp(1:nrho, 1:nrho, ithread) = cmtmp8_omp(1:nrho, 1:nrho, ithread) +   &
                                                         rcgtmp(1:nrho, 1:nrho, lnj) * dtmp2
                 end if
               end do
             end do alf0
             call zamsmm1_omp('l', transa, 'n', iorbs1, ispin, -eye, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                 &
                              cmtmp3_omp(1,1,ithread), nrho, cmtmp4_omp(1,1,ithread), nrho)
             call zamsmm1_omp('r', transa, 'n', iorbs1, ispin,  eye, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                 &
                              cmtmp3_omp(1,1,ithread), nrho, cmtmp4_omp(1,1,ithread), nrho)
             call zamsmm1_omp('l', transa, 'c', iorbs1, ispin, -eye, cmtmp7_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                 &
                              cmtmp3_omp(1,1,ithread), nrho, cmtmp4_omp(1,1,ithread), nrho)
             call zamsmm1_omp('r', transa, 'c', iorbs1, ispin,  eye, cmtmp8_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                 &
                              cmtmp3_omp(1,1,ithread), nrho, cmtmp4_omp(1,1,ithread), nrho)
          end do orbs0
        end do
      end do
    end if
!
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
    rcgrhs(1:nrho,1:nrho,lni) = rcgrhs(1:nrho,1:nrho,lni) + cmtmp2_omp(1:nrho,1:nrho,ithread)
    if (cdabs(ejac - cgamma_nk) .gt. dpico) then
        rcgrhs(1:nrho,1:nrho,lni) = rcgrhs(1:nrho,1:nrho,lni) +  &
                                    cmtmp3_omp(1:nrho,1:nrho,ithread) / (ejac - cgamma_nk)
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
300 continue
!
Return
end subroutine calc_ax_jacobi_omp
