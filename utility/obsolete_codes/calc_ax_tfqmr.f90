subroutine calc_ax_tfqmr(dim0, rcgtmp, rcgrhs)
!
! purpose : calculate the RHS of EOM for rho
!    
! input   : rcgtmp
! output  : rcgrhs
!         
use matmod_omp
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer,    intent(in)    :: dim0
complex*16, intent(in)    :: rcgtmp(dim0,dim0,*)
complex*16, intent(out)   :: rcgrhs(dim0,dim0,*)
integer                   :: itier, itype, dim1, ifactfront, ifactrear
integer                   :: index_nk(MAXTIER)
integer*8                 :: lni, lnj, lnk, lnl, lnm
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, jdraw
real*8                    :: cpu1, cpu2
real*8                    :: dtmp1, dtmp2
complex*16                :: ctmp1, ctmp2
character*1               :: transa, transb
complex*16                :: cgamma_nk
integer                   :: itimes
data itimes /0/
!
dim1   = nrho**2
itimes = itimes + 1
!
if (nthreads_omp > 1) then
   call calc_ax_tfqmr_omp(dim0, rcgtmp, rcgrhs)
   return 
end if
!
if (dim0 .ne. nrho) then
  write(6,*)
  write(6,*)' error! dim0 != nrho in calc_ax_tfqmr ', dim0, nrho
  stop
end if
!
if (lhf) then
   call calchs(rcgtmp(1,1,1))
end if
if (igroundsteady .eq. 0) then
   cmtmp3 = hs
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d20, rcgtmp(1,1,1))
   cmtmp3 = hs + dhs
else
  write(6,*)
  write(6,*)' error! unknown igroundsteady in calc_ax_tfqmr ', igroundsteady
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
rcgrhs(1:nrho, 1:nrho, 1:nunk) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcgtmp(1,1,1), nrho, &
           czero, cmtmp2, nrho)
call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcgtmp(1,1,1), nrho, &
           cunity, cmtmp2, nrho)
rcgrhs(1:nrho, 1:nrho, 1) = -eye * cmtmp2(1:nrho, 1:nrho)
!
!goto 333
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
        cmtmp4(1:nrho, 1:nrho) = cmtmp4(1:nrho, 1:nrho) + rcgtmp(1:nrho, 1:nrho, lni) * dtmp1
      end do
    end do
!    if (lcf) then
!       call zamsmm1('l', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
!    else 
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
!    end if
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  eye, cmtmp4, nrho, cunity, cmtmp2, nrho)
  end do
end do
!
do nj=1,nrho
  do ni=1,nrho
    rcgrhs(ni,nj,1) = rcgrhs(ni,nj,1) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
  end do
end do
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)
!
  do lni=nfirst(itier), nlast(itier)
    cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
    do ni=1,itier-1 
     index_nk(ni) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
     cgamma_nk    = cgamma_nk + cgamma(index_nk(ni))
!
     if (igroundsteady .ne. 0) then
        cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni)) * eleadinfty(ilead(index_nk(ni)), mspin(index_nk(ni)))
     end if
    end do
!
! same itier contribution (N)
!
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcgtmp(1,1,lni), nrho,  czero, cmtmp2, nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcgtmp(1,1,lni), nrho, cunity, cmtmp2, nrho)
    rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) - eye * cmtmp2(1:nrho,1:nrho) +  &
                                  cgamma_nk * rcgtmp(1:nrho, 1:nrho, lni)
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2 = czero
    lnl    = index_coef_ref(lni)
    do ni=1,itier-1
      lnj        = indexcoef(lnl)
      itype      = itypecoef(lnl)
      ifactfront = ifactfcoef(lnl)
      ifactrear  = ifactrcoef(lnl)
!      if (lcf) then
!         ifactfront = -ifactfront
!      end if
      lnl        = lnl + 1
!------
!cycle ! -
!------
      nn = index_nk(ni)
      transa = 'n'
      if (mpm(nn) .eq. 1) transa = 'c'
      transb = 'n'
      if (itype .ne. 0) transb = 'c'

      if (.not. offcor) then
         ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
         ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
         call zamsmm1('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgtmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
         call zamsmm1('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgtmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
      else
         do iorbs2=1,norbs
            ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)  !  linewidth(morbs(nn),iorbs2)
            ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront))
            ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) )
            call zamsmm1('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgtmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
            call zamsmm1('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgtmp(1,1,lnj), nrho, cunity, cmtmp2, nrho)
         end do
      end if
    end do
    rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) + cmtmp2(1:nrho, 1:nrho)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2 = czero
      do isgn=1,nsgn
        transa = 'n'
        if (isgn .ne. 1) transa = 'c'
        do ispin=1,nspin
          do iorbs1=1,norbs
             cmtmp5(1:nrho,1:nrho) = czero
             cmtmp6(1:nrho,1:nrho) = czero
             cmtmp7(1:nrho,1:nrho) = czero
             cmtmp8(1:nrho,1:nrho) = czero
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
!                 if (lcf) then
!                    ifactfront = -ifactfront
!                 end if
                 lnl        = lnl + 1
!
                 dtmp1 = dble(ifactfront)
                 dtmp2 = dble(ifactrear )
                 if (lscale) then
                    dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if     
                 if (itype .eq. 0) then
                    cmtmp5(1:nrho, 1:nrho) = cmtmp5(1:nrho, 1:nrho) + rcgtmp(1:nrho, 1:nrho, lnj) * dtmp1
                    cmtmp6(1:nrho, 1:nrho) = cmtmp6(1:nrho, 1:nrho) + rcgtmp(1:nrho, 1:nrho, lnj) * dtmp2
                 else 
                    cmtmp7(1:nrho, 1:nrho) = cmtmp7(1:nrho, 1:nrho) + rcgtmp(1:nrho, 1:nrho, lnj) * dtmp1
                    cmtmp8(1:nrho, 1:nrho) = cmtmp8(1:nrho, 1:nrho) + rcgtmp(1:nrho, 1:nrho, lnj) * dtmp2
                 end if
               end do
             end do
             call zamsmm1('l', transa, 'n', iorbs1, ispin, -eye, cmtmp5, nrho, cunity, cmtmp2, nrho)
             call zamsmm1('r', transa, 'n', iorbs1, ispin,  eye, cmtmp6, nrho, cunity, cmtmp2, nrho)
             call zamsmm1('l', transa, 'c', iorbs1, ispin, -eye, cmtmp7, nrho, cunity, cmtmp2, nrho)
             call zamsmm1('r', transa, 'c', iorbs1, ispin,  eye, cmtmp8, nrho, cunity, cmtmp2, nrho)
          end do
        end do
      end do
      rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) + cmtmp2(1:nrho, 1:nrho)
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
!if (.not. lcf) then
rcgrhs(1,1,1) = czero
do ni=1,nrho
  rcgrhs(1,1,1) = rcgrhs(1,1,1) + rcgtmp(ni,ni,1)
end do
!end if
!
300 continue
!
end subroutine calc_ax_tfqmr
