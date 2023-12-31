subroutine calc_ax_tfqmr_omp(dim0, rcgtmp, rcgrhs)
!
! purpose : calculate the RHS of EOM for rho
!    
use matmod
use hbmod
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
                           ! for lhb=.true., the last dimension is nunk*(nbath_hb+1)
integer                   :: itier, itype, dim1, ifactfront, ifactrear, ifactf, ifactr
integer                   :: index_nk(MAXTIER, maxprocs), index_nj(MAXTIER, maxprocs)
integer*8                 :: lni, lnj, lnk, lnl, lnm
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: imode
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, iopr
integer                   :: nrho2, nwork, info, mi, mj, istat, mdim
real*8                    :: cpu1, cpu2
real*8                    :: dtmp1, dtmp2
complex*16                :: ctmp1, ctmp2
character*1               :: transa, transb, transc
complex*16                :: cgamma_nk, cgamma_nj
integer                   :: ithread, icore
!real*8,     allocatable   :: rwork(:), tmpval(:)
!complex*16, allocatable   :: zwork(:)
integer                   :: itimes
data itimes /0/
!
if (itimes .eq. 0) then
   write(6,*)'calc_ax_tfqmr_omp: first entry'
end if
!
dim1   = nrho**2
itimes = itimes + 1
!
if (dim0 .ne. nrho) then
  write(6,*)
  write(6,*)' error! dim0 != nrho in calc_ax_tfqmr_omp ', dim0, nrho
  stop
end if
if (lhf) then 
   call calchs(rcgtmp(1,1,1))
end if
if (itimes .eq. 1 .and. ltrun_der) then
   call prelude_trun(rcgtmp(1,1,1))
end if
!
! Bug cleared here
! subroutine 'prelude_trun' should be called BEFORE assigning values to cmtmp3
! as cmtmp3 is reused in 'calchs' and 'calcdhs'
!
if (igroundsteady .eq. 0) then
   cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d20, rcgtmp(1,1,1))
   cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
else
   write(6,*)'calc_ax_tfqmr_omp: error! unknown igroundsteady ', igroundsteady
   stop
end if
!
! diagonalize system Liouvillian [Hamiltonian, *]
!
if (ltrun_der) then
   mdim  = MAXTIER 
   nrho2 = nrho * nrho  
!   call zhil2liou('l', 'n', nrho, cmtmp3, nrho, zlmtmp1, nrho2)
!   call zhil2liou('r', 'n', nrho, cmtmp3, nrho, zlmtmp2, nrho2)
!   zlmtmp1(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
!   nwork = nrho2**2
!   allocate(zwork(nwork), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
!   call zheev('V', 'u', nrho2, zlmtmp1, nrho2, tmpval, zwork, nwork, rwork, info)
!   if (info .ne. 0) then
!      write(6,*)'calc_ax_tfqmr_omp: error! diagonalization failed 1 ', info
!      stop
!   end if
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
rcgrhs(1:nrho,1:nrho,1:nunk) = czero
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
cmtmp2(1:nrho,1:nrho) = czero
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
        cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) +         &
                                            rcgtmp(1:nrho,1:nrho,lni) * dtmp1
      end do
!$omp end do
!$omp end parallel
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
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
! heat bath
!
if (lhb) then
   lni = 1
   cmtmp1(1:nrho,1:nrho) = czero
   do imode=1,nmode_hb
      call calc_tildeq_hb(lni, imode, rcgtmp(1,1,1), cmtmp2)
      call zgemm('n', 'n', nrho, nrho, nrho,  cunity, qa_hb(1,1,imode), nrho, cmtmp2, nrho, &
                 cunity, cmtmp1, nrho)      
      call zgemm('n', 'n', nrho, nrho, nrho, -cunity, cmtmp2, nrho, qa_hb(1,1,imode), nrho, &
                 cunity, cmtmp1, nrho)      
   end do
   rcgrhs(1:nrho,1:nrho,1) = rcgrhs(1:nrho,1:nrho,1) - cmtmp1(1:nrho,1:nrho)
end if
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn) &
!$omp private(ntmp3,isgn,imode,ispin,iorbs1,iorbs2,ialf,imats,dtmp1,dtmp2,i,j,ithread) &
!$omp private(transa,transb,transc,ctmp1,ctmp2,idraw,cgamma_nj,mi,mj,ifactf,ifactr) &
!$omp private(iopr,lnm)
  ithread = omp_get_thread_num() + 1
!$omp do schedule(dynamic) 
  do lni=nfirst(itier), nlast(itier)
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
       if (igroundsteady .ne. 0) then
          cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                             &
                      eleadinfty(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
       end if
    end do
!
! same itier contribution (N)
!
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcgtmp(1,1,lni), nrho,  czero,     &
               cmtmp2_omp(1,1,ithread), nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcgtmp(1,1,lni), nrho, cunity,     &
               cmtmp2_omp(1,1,ithread), nrho)
    rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) - eye * cmtmp2_omp(1:nrho,1:nrho,ithread) +  &
                                  cgamma_nk * rcgtmp(1:nrho, 1:nrho, lni)
!
! heat bath
! 
    if (lhb) then
       cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
       do imode=1,nmode_hb
          call calc_tildeq_hb(lni, imode, rcgtmp(1,1,1), cmtmp2_omp(1,1,ithread))
          call zgemm('n', 'n', nrho, nrho, nrho,  cunity, qa_hb(1,1,imode), nrho, cmtmp2_omp(1,1,ithread), nrho, &
                     cunity, cmtmp1_omp(1,1,ithread), nrho)
          call zgemm('n', 'n', nrho, nrho, nrho, -cunity, cmtmp2_omp(1,1,ithread), nrho, qa_hb(1,1,imode), nrho, &
                     cunity, cmtmp1_omp(1,1,ithread), nrho)
       end do
       rcgrhs(1:nrho,1:nrho,lni) = rcgrhs(1:nrho,1:nrho,lni) - cmtmp1_omp(1:nrho,1:nrho,ithread)
    end if
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
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
      nn = index_nk(ni,ithread)
!
      if (psdfff .and. mmfff(nn) .gt. 1) then
          ctmp1 = dble(ifactfront) * ctheta(mspin(nn),mcor(nn),ilead(nn),mpm(nn))
          if (itype .ne. 0) then
              write(6,*)
              write(6,*)'calc_ax_tfqmr_omp: unexpected itype ', itype
              stop
          end if
          rcgrhs(1:nrho,1:nrho,lni) = rcgrhs(1:nrho,1:nrho,lni) + ctmp1 * rcgtmp(1:nrho,1:nrho,lnj)
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
    rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) + cmtmp2_omp(1:nrho, 1:nrho,ithread)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
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
!                 if (lcf) then 
!                    ifactfront = -ifactfront
!                 end if
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
      rcgrhs(1:nrho, 1:nrho, lni) = rcgrhs(1:nrho, 1:nrho, lni) + cmtmp2_omp(1:nrho, 1:nrho,ithread)
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
                      cgamma_nj = cgamma_nk + cgamma(idraw)
                      if (igroundsteady .ne. 0) then
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
!                      if (cdabs(cgamma_nj) .gt. dpico) then
!                         do mj=1,nnzeig
!                            ctmp1 = cunity / (eye * dvalder(mj) - cgamma_nj) + cunity / cgamma_nj
!                            zlmtmp1_omp(1:nrho2,mj,ithread) = zmatder(1:nrho2,mj) * ctmp1
!                         end do
!                      else
!                         do mj=1,nnzeig
!                            ctmp1 = -eye / dvalder(mj)
!                            zlmtmp1_omp(1:nrho2,mj,ithread) = zmatder(1:nrho2,mj) * ctmp1
!                         end do
!                      end if
!
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
                            call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgtmp(1,1,lnj), nrho, cunity,               &
                                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgtmp(1,1,lnj), nrho, cunity,               &
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
       rcgrhs(1:nrho,1:nrho,lni) = rcgrhs(1:nrho,1:nrho,lni) + cmtmp3_omp(1:nrho,1:nrho,ithread)
!       
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
!if (.not. lcf) then
!if (ltrend) then
rcgrhs(nrho,nrho,nunk) = czero
do ni=1,nrho
   rcgrhs(nrho,nrho,nunk) = rcgrhs(nrho,nrho,nunk) + rcgtmp(ni,ni,1)
end do
!else 
!    rcgrhs(1,1,1) = czero
!    do ni=1,nrho
!       rcgrhs(1,1,1) = rcgrhs(1,1,1) + rcgtmp(ni,ni,1)
!    end do
!end if
!
300 continue
!
!if (ltrun_der) then
!   deallocate(zwork, rwork, tmpval, STAT=istat)
!end if
!
return
end subroutine calc_ax_tfqmr_omp
