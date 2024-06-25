subroutine calcderiv_omp(tt)
!
! purpose : calculate the RHS of EOM for rho
!    
! input   : rhotmp (not changed)
! output  : rhorhs (refreshed)
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
real*8,  intent(in) :: tt
!
integer    :: itier, itype, ifactfront, ifactrear, ifactf, ifactr
integer    :: index_nk(MAXTIER, maxprocs), index_nj(MAXTIER, maxprocs)
integer*8  :: lni, lnj, lnk, lnl, lnm
integer    :: ni, nj, nk, nl, nm, nn, ntmp3
integer    :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, dim1, idraw, iopr
integer    :: imode
integer    :: nrho2, nwork, info, mi, mj, istat, mdim
real*8     :: cpu1, cpu2, dtmp1, dtmp2, dtmp3, dtmp4
real*8     :: eshift(maxalf,maxspin)
complex*16 :: cgamma_nk, cgamma_nj, ctmp1, ctmp2, ctmp3, ctmp4
character*1  :: transa, transb, transc
integer    :: itimes
data itimes /0/
!
integer    :: ithread, icore
!
if (itimes .eq. 0) then
   write(6,*)'calcderiv_omp: first entry '
end if
!
dim1   = nrho**2
itimes = itimes + 1
!
if (ltrun_der) then
   mdim  = MAXTIER 
   nrho2 = nrho * nrho  
   call refresh_trun(tt, rhotmp(1,1,1))
end if
!
if (lhf) then
   call calchs(rhotmp(1,1,1))
end if
call calcdhs(tt, rhotmp(1,1,1))
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
rhorhs = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rhotmp(1,1,1), nrho, &
           czero, cmtmp2, nrho)
call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rhotmp(1,1,1), nrho, &
           cunity, cmtmp2, nrho)
rhorhs(1:nrho, 1:nrho, 1) = -eye * cmtmp2(1:nrho, 1:nrho)
!
cmtmp1 = czero
cmtmp2 = czero
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4 = czero
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
        cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) +         &
                                            rhotmp(1:nrho,1:nrho,lni) * dtmp1
      end do
!$omp end do
!$omp end parallel      
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
      end do
    end do
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp1, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
  end do
end do
do nj=1,nrho
  do ni=1,nrho
    ctmp1 = cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni)) + &
            cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
    rhorhs(ni,nj,1) = rhorhs(ni,nj,1) + ctmp1
  end do
end do
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,itype,ifactfront,ifactrear,nn) &
!$omp private(isgn,ispin,iorbs1,iorbs2,ialf,imats,ntmp3,ctmp1,ctmp2,i,j,ithread,imode) &
!$omp private(transa,transb,transc,idraw,cgamma_nj,mi,mj,ifactf,ifactr) &
!$omp private(lnm,iopr)
  ithread = omp_get_thread_num() + 1
!$omp do schedule(dynamic) 
  do lni=nfirst(itier),nlast(itier)
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
       cgamma_nk            = cgamma_nk + cgamma(index_nk(ni,ithread)) + eye * dipm(index_nk(ni,ithread)) *  &
                              eshift(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
    end do
!
! same itier contribution (N)
!
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rhotmp(1,1,lni),   &
               nrho,  czero, cmtmp2_omp(1,1,ithread), nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rhotmp(1,1,lni),   &
               nrho, cunity, cmtmp2_omp(1,1,ithread), nrho)
    rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) - eye * cmtmp2_omp(1:nrho, 1:nrho, ithread) + &
                                  cgamma_nk * rhotmp(1:nrho, 1:nrho, lni) 
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
    do ni=1,itier-1
      lnj        = indexcoef(lnl)
      itype      = itypecoef(lnl)
      ifactfront = ifactfcoef(lnl)
      ifactrear  = ifactrcoef(lnl)
      lnl        = lnl + 1
!
      nn = index_nk(ni,ithread)
!
      if (psdfff .and. mmfff(nn) .gt. 1) then
          ctmp1 = dble(ifactfront) * ctheta(mspin(nn),mcor(nn),ilead(nn),mpm(nn))
          if (itype .ne. 0) then
              write(6,*)
              write(6,*)'calcderiv_omp: unexpected itype ', itype
              stop
          end if
          rhorhs(1:nrho,1:nrho,lni) = rhorhs(1:nrho,1:nrho,lni) + ctmp1 * rhotmp(1:nrho,1:nrho,lnj)
          cycle
      end if
!      
      if (.not. offcor) then
        if (itype .eq. 0) then
          if (mpm(nn) .eq. 2) then
            call zamsmm1_omp('l', 'n', 'n', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
            call zamsmm1_omp('r', 'n', 'n', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                             
          else
            call zamsmm1_omp('l', 'c', 'n', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                             
            call zamsmm1_omp('r', 'c', 'n', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                             
          end if
        else
          if (mpm(nn) .eq. 2) then
            call zamsmm1_omp('l', 'n', 'c', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                             
            call zamsmm1_omp('r', 'n', 'c', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
          else
            call zamsmm1_omp('l', 'c', 'c', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                             
            call zamsmm1_omp('r', 'c', 'c', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                             dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,          &
                             cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                                                          
          end if
        end if
      else  ! offcor
        do iorbs2=1,norbs
          ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
          if (itype .eq. 0) then
            if (mpm(nn) .eq. 2) then
              call zamsmm1_omp('l', 'n', 'n', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                                                            
              call zamsmm1_omp('r', 'n', 'n', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
            else
              call zamsmm1_omp('l', 'c', 'n', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
              call zamsmm1_omp('r', 'c', 'n', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
            end if
          else
            if (mpm(nn) .eq. 2) then
              call zamsmm1_omp('l', 'n', 'c', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
              call zamsmm1_omp('r', 'n', 'c', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
            else
              call zamsmm1_omp('l', 'c', 'c', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
              call zamsmm1_omp('r', 'c', 'c', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                               dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho,   &
                               cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)                               
            end if
          end if
        end do
      end if
    end do
    rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) + cmtmp2_omp(1:nrho, 1:nrho, ithread)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
      do isgn=1,nsgn
        do ispin=1,nspin
          orbs0: do iorbs1=1,norbs
             if (lsimple .and. .not. lwalf) then
                iopr = ioprindex(lnm)
                lnm = lnm + 1
                if (iopr .lt. 0) cycle orbs0
             end if
             cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
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
                 lnl   = lnl + 1
!
                 ctmp1 = dble(ifactfront)
                 ctmp2 = dble(ifactrear )
                 if (lscale) then
                    ctmp1 = ctmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    ctmp2 = ctmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if
                 if (itype .eq. 0) then
                   cmtmp5_omp(1:nrho, 1:nrho, ithread) = cmtmp5_omp(1:nrho, 1:nrho, ithread) + rhotmp(1:nrho, 1:nrho, lnj) * ctmp1
                   cmtmp6_omp(1:nrho, 1:nrho, ithread) = cmtmp6_omp(1:nrho, 1:nrho, ithread) + rhotmp(1:nrho, 1:nrho, lnj) * ctmp2
                 else
                   do j=1,nrho
                     do i=1,nrho
                       cmtmp5_omp(i,j, ithread) = cmtmp5_omp(i,j, ithread) + dconjg(rhotmp(j,i,lnj)) * ctmp1
                       cmtmp6_omp(i,j, ithread) = cmtmp6_omp(i,j, ithread) + dconjg(rhotmp(j,i,lnj)) * ctmp2
                     end do
                   end do
                 end if
               end do
             end do alf0
             if (isgn .eq. 1) then
               call zamsmm1_omp('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp5_omp(1,1,ithread), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho, &
                                cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
               call zamsmm1_omp('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp6_omp(1,1,ithread), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho, &
                                cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)               
             else
               call zamsmm1_omp('l', 'c', 'n', iorbs1, ispin, -eye, cmtmp5_omp(1,1,ithread), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho, &
                                cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)               
               call zamsmm1_omp('r', 'c', 'n', iorbs1, ispin,  eye, cmtmp6_omp(1,1,ithread), nrho, cunity, cmtmp2_omp(1,1,ithread), nrho, &
                                cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)               
             end if
!
          end do orbs0
        end do
      end do
      rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) + cmtmp2_omp(1:nrho, 1:nrho, ithread)
    end if
!
! Truncation for the derivatives of terminal-tier ADOs
!
    if (ltrun_der .and. itier .eq. ntier) then
       cmtmp3_omp(1:nrho,1:nrho,ithread) = czero
       do isgn=1,nsgn
          transc = 'n'
          if (isgn .ne. 1) transc = 'c'
!
          do ispin=1,nspin
             do iorbs1=1,norbs
                cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
                cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
                do ialf=1,nalf
                   call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                   if (lscreen .and. jomit(iopr) .eq. 1) cycle
                   !
                   do imats=1,ncor
                      ifactf = ifactfcoef(lnl)
                      ifactr = ifactrcoef(lnl)
                      idraw  = indexcoef(lnl)
                      lnl = lnl + 1
                      cgamma_nj = cgamma_nk + cgamma(idraw) + eye * dipm(idraw) * eshift(ilead(idraw), mspin(idraw))
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
                      cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
                      do mi=1,ntier
                         lnj        = indexcoef(lnl)
                         itype      = itypecoef(lnl)
                         ifactfront = ifactfcoef(lnl)
                         ifactrear  = ifactrcoef(lnl)
                         lnl = lnl + 1
                         nn  = index_nj(mi,ithread)
                         transa = 'n'
                         if (mpm(nn) .eq. 1) transa = 'c'
                         transb = 'n'
                         if (itype .ne. 0) transb = 'c'
                         if (.not. offcor) then
                            ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
                            ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
                            call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rhotmp(1,1,lnj), nrho, cunity,               &
                                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rhotmp(1,1,lnj), nrho, cunity,               &
                                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                         else
                            do iorbs2=1,norbs
                               ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
                               ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
                               ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
                               call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), ctmp1, rhotmp(1,1,lnj), nrho, cunity,               &
                                                cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                               call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), ctmp2, rhotmp(1,1,lnj), nrho, cunity,               &
                                                cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            end do
                         end if
                      end do                   
                      if (cdabs(cgamma_nj) .gt. dpico) then
                         cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) -                   &
                                                             cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp1 / cgamma_nj
                         cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp6_omp(1:nrho,1:nrho,ithread) -                   &
                                                             cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp2 / cgamma_nj
                      end if
                      call zgemv('c', nrho2, nnzeig, cunity, zmatder, nrho2,                     &
                                 cmtmp2_omp(1,1,ithread), 1, czero, cmtmp1_omp(1,1,ithread), 1)
                      ntmp3 = 0
                      scale0: do mj=1,nrho
                         do mi=1,nrho
                            ntmp3 = ntmp3 + 1
                            if (ntmp3 .gt. nnzeig) exit scale0
                            if (cdabs(cgamma_nj) .gt. dpico) then
                               ctmp1 = cunity / (eye * dvalder(ntmp3) - cgamma_nj) + cunity / cgamma_nj
                            else
                               ctmp1 = -eye / dvalder(ntmp3)
                            end if
                            cmtmp1_omp(mi,mj,ithread) = cmtmp1_omp(mi,mj,ithread) * ctmp1
                         end do
                      end do scale0
                      call zgemv('n', nrho2, nnzeig, cunity, zmatder, nrho2, cmtmp1_omp(1,1,ithread), 1,  &
                                 czero, cmtmp2_omp(1,1,ithread), 1)
                      cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp1
                      cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp6_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp2
                   end do
                end do
                call zamsmm1_omp('l', transc, 'n', iorbs1, ispin, -eye, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp3_omp(1,1,ithread), nrho,                                                 &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                call zamsmm1_omp('r', transc, 'n', iorbs1, ispin,  eye, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp3_omp(1,1,ithread), nrho,                                                 &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
             end do
          end do
       end do
       rhorhs(1:nrho,1:nrho,lni) = rhorhs(1:nrho,1:nrho,lni) + cmtmp3_omp(1:nrho,1:nrho,ithread)
!       
    end if
!
  end do
!$omp end do
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
end subroutine calcderiv_omp
