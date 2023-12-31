subroutine calcderiv(tt)
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
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,  intent(in) :: tt
!
integer    :: itier, itype, ifactfront, ifactrear
integer    :: index_nk(MAXTIER)
integer*8  :: lni, lnj, lnk, lnl, lnm
integer    :: ni, nj, nk, nl, nm, nn, ntmp3
integer    :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, dim1
real*8     :: cpu1, cpu2, dtmp1, dtmp2, dtmp3, dtmp4
real*8     :: eshift(maxalf,maxspin)
complex*16 :: cgamma_nk, ctmp1, ctmp2, ctmp3, ctmp4
integer    :: itimes
data itimes /0/
!
if (lpara_omp) then
   call calcderiv_omp(tt)
   return
end if
!
dim1   = nrho**2
itimes = itimes + 1
!
if (lhf) then
   call calchs(rhotmp(1,1,1))
end if
call calcdhs(tt, rhotmp(1,1,1))
cmtmp3 = hs + dhs
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
      do imats=1,ncor
        if (lscale) then
           dtmp1 = dbsqrt(iorbs1,ispin,imats,ialf,1)
        else
           dtmp1 = 1.d0
        end if     
        call look4drawer(1, ialf, iorbs1, ispin, imats, nn)
        lni = nfirst(2) - 1 + nn
        cmtmp4(1:nrho, 1:nrho) = cmtmp4(1:nrho, 1:nrho) + rhotmp(1:nrho, 1:nrho, lni) * dtmp1
      end do
    end do
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp1, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
  end do
end do
do nj=1,nrho
   do ni=1,nrho
      ctmp1 = cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni)) +  &
              cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
      rhorhs(ni,nj,1) = rhorhs(ni,nj,1) + ctmp1
   end do
end do
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)
  do lni=nfirst(itier),nlast(itier)
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
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rhotmp(1,1,lni),   &
               nrho, czero, cmtmp2, nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rhotmp(1,1,lni),   &
               nrho, cunity, cmtmp2, nrho)
    rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) - eye * cmtmp2(1:nrho, 1:nrho) + &
                                  cgamma_nk * rhotmp(1:nrho, 1:nrho, lni) 
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2 = czero
    do ni=1,itier-1
      lnj        = indexcoef(lnl)
      itype      = itypecoef(lnl)
      ifactfront = ifactfcoef(lnl)
      ifactrear  = ifactrcoef(lnl)
      lnl        = lnl + 1
!
      nn = index_nk(ni)
      if (.not. offcor) then
        if (itype .eq. 0) then
          if (mpm(nn) .eq. 2) then
            call zamsmm1('l', 'n', 'n', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            call zamsmm1('r', 'n', 'n', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
          else
            call zamsmm1('l', 'c', 'n', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            call zamsmm1('r', 'c', 'n', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
          end if
        else
          if (mpm(nn) .eq. 2) then
            call zamsmm1('l', 'n', 'c', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            call zamsmm1('r', 'n', 'c', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
          else
            call zamsmm1('l', 'c', 'c', morbs(nn), mspin(nn), -eye * cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            call zamsmm1('r', 'c', 'c', morbs(nn), mspin(nn),  eye * cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                        dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
          end if
        end if
      else  ! offcor
        do iorbs2=1,norbs
          ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
          if (itype .eq. 0) then
            if (mpm(nn) .eq. 2) then
              call zamsmm1('l', 'n', 'n', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
              call zamsmm1('r', 'n', 'n', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            else
              call zamsmm1('l', 'c', 'n', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
              call zamsmm1('r', 'c', 'n', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            end if
          else
            if (mpm(nn) .eq. 2) then
              call zamsmm1('l', 'n', 'c', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
              call zamsmm1('r', 'n', 'c', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            else
              call zamsmm1('l', 'c', 'c', iorbs2, mspin(nn), -eye * cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactfront), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
              call zamsmm1('r', 'c', 'c', iorbs2, mspin(nn),  eye * cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * &
                          dble(ifactrear ), rhotmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            end if
          end if
        end do
      end if
    end do
    rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) + cmtmp2(1:nrho, 1:nrho)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2 = czero
      do isgn=1,nsgn
        do ispin=1,nspin
          do iorbs1=1,norbs
             cmtmp5 = czero
             cmtmp6 = czero
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
                 ctmp1 = dble(ifactfront)
                 ctmp2 = dble(ifactrear )
                 if (lscale) then
                    ctmp1 = ctmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    ctmp2 = ctmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if
                 if (itype .eq. 0) then
                   cmtmp5(1:nrho, 1:nrho) = cmtmp5(1:nrho, 1:nrho) + rhotmp(1:nrho, 1:nrho, lnj) * ctmp1
                   cmtmp6(1:nrho, 1:nrho) = cmtmp6(1:nrho, 1:nrho) + rhotmp(1:nrho, 1:nrho, lnj) * ctmp2
                 else
                   do j=1,nrho
                     do i=1,nrho
                       cmtmp5(i,j) = cmtmp5(i,j) + dconjg(rhotmp(j,i,lnj)) * ctmp1
                       cmtmp6(i,j) = cmtmp6(i,j) + dconjg(rhotmp(j,i,lnj)) * ctmp2
                     end do
                   end do
                 end if
               end do
             end do
             if (isgn .eq. 1) then
               call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp5, nrho, cunity, cmtmp2, nrho)
               call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp6, nrho, cunity, cmtmp2, nrho)
             else
               call zamsmm1('l', 'c', 'n', iorbs1, ispin, -eye, cmtmp5, nrho, cunity, cmtmp2, nrho)
               call zamsmm1('r', 'c', 'n', iorbs1, ispin,  eye, cmtmp6, nrho, cunity, cmtmp2, nrho)
             end if
!
          end do
        end do
      end do
      rhorhs(1:nrho, 1:nrho, lni) = rhorhs(1:nrho, 1:nrho, lni) + cmtmp2(1:nrho, 1:nrho)
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
end subroutine calcderiv
