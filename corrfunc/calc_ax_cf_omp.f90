subroutine calc_ax_cf_omp(dim0, rcga, rcgb, rhsa, rhsb)
!
! purpose : calculate the RHS of EOM for rho
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
complex*16, intent(in)    :: rcga(dim0,dim0,*), rcgb(dim0,dim0,*)
complex*16, intent(out)   :: rhsa(dim0,dim0,*), rhsb(dim0,dim0,*)
integer                   :: itier, itype, dim1, ifactfront, ifactrear, ifactf, ifactr
integer                   :: index_nk(MAXTIER, maxprocs), index_nj(MAXTIER, maxprocs)
integer*8                 :: lni, lnj, lnk, lnl, lnm
integer                   :: ni, nj, nk, nl, nm, nn, ntmp3
integer                   :: imode
integer                   :: i, j, iorbs1, iorbs2, ispin, imats, ialf, isgn, iorbs, idraw, jdraw
integer                   :: iopr
real*8                    :: cpu1, cpu2
real*8                    :: dfacta, dfactb
real*8                    :: dtmp1, dtmp2, dtmp3
complex*16                :: ctmp1, ctmp2
integer                   :: nrho2, mi, mj, istat, mdim, isgnw, nrho2d
complex*16                :: cgamma_nk, cgamma_nj
character*1               :: transa, transb, transc
complex*16                :: zfront, zrear, zf, zr
integer                   :: itimes
data itimes /0/
!
integer    :: ithread, icore
!
if (itimes .eq. 0) then
   write(6,*)'calc_ax_cf_omp: first entry'
end if
if (dim0 .ne. nrho) then
  write(6,*)
  write(6,*)' error! dim0 != nrho in calc_ax_cf_omp ', dim0, nrho
  stop
end if
!
dim1   = nrho**2
itimes = itimes + 1
cmtmp3(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!
!lnl = 1
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
if (ltrun_der) then
   mdim  = MAXTIER 
   nrho2 = nrho * nrho  
   nrho2d = nrho2 * 2
end if
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
rhsa(1:nrho,1:nrho,1:nunk) = czero
rhsb(1:nrho,1:nrho,1:nunk) = czero
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcga(1,1,1), nrho, &
           czero, cmtmp2, nrho)
call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcga(1,1,1), nrho, &
           cunity, cmtmp2, nrho)
rhsa(1:nrho,1:nrho,1) = -eye * cmtmp2(1:nrho,1:nrho)
call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcgb(1,1,1), nrho, &
           czero, cmtmp2, nrho)
call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcgb(1,1,1), nrho, &
           cunity, cmtmp2, nrho)
rhsb(1:nrho,1:nrho,1) = -eye * cmtmp2(1:nrho,1:nrho)
!
if (lfreq) then
   rhsa(1:nrho,1:nrho,1) = rhsa(1:nrho,1:nrho,1) + dfacta * dfreq * rcgb(1:nrho,1:nrho,1)
   rhsb(1:nrho,1:nrho,1) = rhsb(1:nrho,1:nrho,1) + dfactb * dfreq * rcga(1:nrho,1:nrho,1)
end if
!
!goto 333
cmtmp2 = czero
cmtmp1 = czero
do ispin=1,nspin
  do iorbs1=1,norbs
    cmtmp4(1:nrho,1:nrho) = czero
    cmtmp5(1:nrho,1:nrho) = czero
    do ialf=1,nalf
       call look4operator(1, iorbs1, ispin, ialf, iopr)
       if (lscreen .and. jomit(iopr) .eq. 1) cycle
!
!$omp parallel default(shared) private(imats,dtmp1,nn,lni,ithread)
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
        cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) +         &
                                            rcgb(1:nrho,1:nrho,lni) * dtmp1
        cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) +         &
                                            rcga(1:nrho,1:nrho,lni) * dtmp1
      end do
!$omp end do
!$omp end parallel
      do icore=1,nthreads_omp
         cmtmp4(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) + cmtmp4_omp(1:nrho,1:nrho,icore)
         cmtmp5(1:nrho,1:nrho) = cmtmp5(1:nrho,1:nrho) + cmtmp5_omp(1:nrho,1:nrho,icore)
      end do
    end do
    call zamsmm1('l', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin, -cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
    call zamsmm1('l', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
    call zamsmm1('r', 'n', 'n', iorbs1, ispin,  cunity, cmtmp5, nrho, cunity, cmtmp1, nrho)
  end do
end do
!
do nj=1,nrho
   do ni=1,nrho
      rhsa(ni,nj,1) = rhsa(ni,nj,1) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
      rhsb(ni,nj,1) = rhsb(ni,nj,1) + cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni))
   end do
end do
333 continue
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)

!$omp parallel default(shared) private(lni,lnj,lnk,lnl,cgamma_nk,ni,nj,itype,ifactfront,ifactrear,nn)  &
!$omp private(ntmp3,isgn,ispin,iorbs1,iorbs2,ialf,imats,dtmp1,dtmp2,i,j,ithread,imode)                 &
!$omp private(transa,transb,transc,ctmp1,ctmp2,idraw,cgamma_nj,mi,mj,ifactf,ifactr,zf,zr,zfront,zrear) &
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
!
     if (igroundsteady .eq. 1) then
        cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                             &
                    eleadinfty(ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
!
!     else if (igroundsteady .eq. 2) then
!        cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni,ithread)) *                             &
!                    eleadtime (ilead(index_nk(ni,ithread)), mspin(index_nk(ni,ithread)))
!
     end if
    end do
!
! same itier contribution (N)
!
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcga(1,1,lni), nrho,  czero,     &
               cmtmp2_omp(1,1,ithread), nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcga(1,1,lni), nrho, cunity,     &
               cmtmp2_omp(1,1,ithread), nrho)
    rhsa(1:nrho,1:nrho,lni) = rhsa(1:nrho,1:nrho,lni) - eye * cmtmp2_omp(1:nrho,1:nrho,ithread) +  &
                              cgamma_nk * rcga(1:nrho,1:nrho,lni)
    call zhemm('l', 'u', nrho, nrho,  cunity, cmtmp3, nrho, rcgb(1,1,lni), nrho,  czero,     &
               cmtmp2_omp(1,1,ithread), nrho)
    call zhemm('r', 'u', nrho, nrho, -cunity, cmtmp3, nrho, rcgb(1,1,lni), nrho, cunity,     &
               cmtmp2_omp(1,1,ithread), nrho)
    rhsb(1:nrho,1:nrho,lni) = rhsb(1:nrho,1:nrho,lni) - eye * cmtmp2_omp(1:nrho,1:nrho,ithread) +  &
                              cgamma_nk * rcgb(1:nrho,1:nrho,lni)
    if (lfreq) then
       rhsa(1:nrho,1:nrho,lni) = rhsa(1:nrho,1:nrho,lni) + dfacta * dfreq * rcgb(1:nrho,1:nrho,lni)
       rhsb(1:nrho,1:nrho,lni) = rhsb(1:nrho,1:nrho,lni) + dfactb * dfreq * rcga(1:nrho,1:nrho,lni)
    end if
!
! nearest lower tier contribution (N - 1)
!
    cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
    cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
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
              write(6,*)'calc_ax_cf_omp: unexpected itype ', itype
              stop
          end if
          cmtmp2_omp(1:nrho,1:nrho,ithread) = cmtmp2_omp(1:nrho,1:nrho,ithread) + ctmp1 * rcga(1:nrho,1:nrho,lnj)
          cmtmp1_omp(1:nrho,1:nrho,ithread) = cmtmp1_omp(1:nrho,1:nrho,ithread) + ctmp1 * rcgb(1:nrho,1:nrho,lnj)
          cycle
      end if
!
      transa = 'n'
      if (mpm(nn) .eq. 1) transa = 'c'
      transb = 'n'
      if (itype .ne. 0) transb = 'c'
      if (.not. offcor) then
         ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * (-eye)
         ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) *   eye
         call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcgb(1,1,lnj), nrho, cunity,               &
                          cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcgb(1,1,lnj), nrho, cunity,               &
                          cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcga(1,1,lnj), nrho, cunity,              &
                          cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcga(1,1,lnj), nrho, cunity,              &
                          cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
      else
         do iorbs2=1,norbs
            ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
            ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) * (-eye)
            ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) *   eye
            call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcgb(1,1,lnj), nrho, cunity,               &
                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
            call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcgb(1,1,lnj), nrho, cunity,               &
                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
            call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcga(1,1,lnj), nrho, cunity,              &
                             cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
            call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcga(1,1,lnj), nrho, cunity,              &
                             cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
         end do
      end if
    end do
    rhsa(1:nrho,1:nrho,lni) = rhsa(1:nrho,1:nrho,lni) + cmtmp2_omp(1:nrho,1:nrho,ithread)
    rhsb(1:nrho,1:nrho,lni) = rhsb(1:nrho,1:nrho,lni) + cmtmp1_omp(1:nrho,1:nrho,ithread)
!
! nearest higher tier contribution (N + 1)
!
    if (itier .lt. ntier) then
      cmtmp2_omp(1:nrho,1:nrho,ithread) = czero
      cmtmp1_omp(1:nrho,1:nrho,ithread) = czero
      do isgn=1,nsgn
        transc = 'n'
        if (isgn .ne. 1) transc = 'c'
        do ispin=1,nspin
          orbs0: do iorbs1=1,norbs
             if (lsimple .and. .not. lwalf) then
                iopr = ioprindex(lnm)
                lnm = lnm + 1
                if (iopr .lt. 0) cycle orbs0
             end if
             cmtmp5_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp6_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp3_omp(1:nrho,1:nrho,ithread) = czero
             cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
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
                 zfront = -cunity
                 zrear  = -cunity 
                 dtmp1 = dble(ifactfront)
                 dtmp2 = dble(ifactrear )
                 if (lscale) then
                    dtmp1 = dtmp1 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                    dtmp2 = dtmp2 * dbsqrt(iorbs1,ispin,imats,ialf,isgn)
                 end if     
                 if (itype .eq. 0) then
                    cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) + rcgb(1:nrho,1:nrho,lnj) * dtmp1
                    cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp6_omp(1:nrho,1:nrho,ithread) + rcgb(1:nrho,1:nrho,lnj) * dtmp2
                    cmtmp3_omp(1:nrho,1:nrho,ithread) = cmtmp3_omp(1:nrho,1:nrho,ithread) + rcga(1:nrho,1:nrho,lnj) * dtmp1
                    cmtmp4_omp(1:nrho,1:nrho,ithread) = cmtmp4_omp(1:nrho,1:nrho,ithread) + rcga(1:nrho,1:nrho,lnj) * dtmp2
                 else 
                   do j=1,nrho
                     do i=1,nrho
                       cmtmp5_omp(i,j,ithread) = cmtmp5_omp(i,j,ithread) + dconjg(rcgb(j,i,lnj)) * dtmp1
                       cmtmp6_omp(i,j,ithread) = cmtmp6_omp(i,j,ithread) + dconjg(rcgb(j,i,lnj)) * dtmp2
                       cmtmp3_omp(i,j,ithread) = cmtmp3_omp(i,j,ithread) + dconjg(rcga(j,i,lnj)) * dtmp1
                       cmtmp4_omp(i,j,ithread) = cmtmp4_omp(i,j,ithread) + dconjg(rcga(j,i,lnj)) * dtmp2
                     end do
                   end do
                 end if
               end do
             end do alf0
             call zamsmm1_omp('l', transc, 'n', iorbs1, ispin, zfront, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                   &
                              cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_omp('r', transc, 'n', iorbs1, ispin,  zrear, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp2_omp(1,1,ithread), nrho,                                                   &
                              cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_omp('l', transc, 'n', iorbs1, ispin, -zfront, cmtmp3_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp1_omp(1,1,ithread), nrho,                                                    &
                              cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
             call zamsmm1_omp('r', transc, 'n', iorbs1, ispin,  -zrear, cmtmp4_omp(1,1,ithread), nrho, cunity,  &
                              cmtmp1_omp(1,1,ithread), nrho,                                                    &
                              cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
          end do orbs0
        end do
      end do
      rhsa(1:nrho,1:nrho,lni) = rhsa(1:nrho,1:nrho,lni) + cmtmp2_omp(1:nrho,1:nrho,ithread)
      rhsb(1:nrho,1:nrho,lni) = rhsb(1:nrho,1:nrho,lni) + cmtmp1_omp(1:nrho,1:nrho,ithread)
    end if
!
    if (ltrun_der .and. itier .eq. ntier) then  ! now only for freq=0 case
       cmtmp3_omp(1:nrho,1:nrho,ithread) = czero
       cmtmp4_omp(1:nrho,1:nrho,ithread) = czero
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
!                      do mj=1,nrho2
!                         ctmp1 = cunity / (eye * dvalder(mj) - cgamma_nj)
!                         zlmtmp1_omp(1:nrho2,mj,ithread) = zmatder(1:nrho2,mj) * ctmp1
!                      end do
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
                         transa = 'n'
                         if (mpm(nn) .eq. 1) transa = 'c'
                         transb = 'n'
                         if (itype .ne. 0) transb = 'c'
                         if (.not. offcor) then
! note that ctmp1 is associated with eye, not (-eye), because it is coefficient of 'rcga' 
                            ctmp1 = cb(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) *   eye  
                            ctmp2 = cd(morbs(nn),mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * (-eye)
                            call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), ctmp1, rcga(1,1,lnj), nrho, cunity,               &
                                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), ctmp2, rcga(1,1,lnj), nrho, cunity,               &
                                             cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            call zamsmm1_omp('l', transa, transb, morbs(nn), mspin(nn), -ctmp1, rcgb(1,1,lnj), nrho, cunity,              &
                                             cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                            call zamsmm1_omp('r', transa, transb, morbs(nn), mspin(nn), -ctmp2, rcgb(1,1,lnj), nrho, cunity,              &
                                             cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                         else
                            do iorbs2=1,norbs
                               ntmp3 = (iorbs2 - 1) * norbs + morbs(nn)
                               ctmp1 = cb(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0, -dble(ifactfront)) *   eye
                               ctmp2 = cd(ntmp3,mspin(nn),mcor(nn),ilead(nn),mpm(nn)) * dcmplx(0.d0,  dble(ifactrear) ) * (-eye)
                               call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), ctmp1, rcga(1,1,lnj), nrho, cunity,               &
                                                cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                               call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), ctmp2, rcga(1,1,lnj), nrho, cunity,               &
                                                cmtmp2_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                               call zamsmm1_omp('l', transa, transb, iorbs2, mspin(nn), -ctmp1, rcgb(1,1,lnj), nrho, cunity,              &
                                                cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                               call zamsmm1_omp('r', transa, transb, iorbs2, mspin(nn), -ctmp2, rcgb(1,1,lnj), nrho, cunity,              &
                                                cmtmp1_omp(1,1,ithread), nrho, cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
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
!                      call zgemv('c', nrho2, nrho2, cunity, zmatder, nrho2, cmtmp2_omp(1,1,ithread), 1, czero, cmtmp7_omp(1,1,ithread), 1)
!                      call zgemv('n', nrho2, nrho2, cunity, zlmtmp1_omp(1,1,ithread), nrho2,    &
!                                 cmtmp7_omp(1,1,ithread), 1, czero, cmtmp2_omp(1,1,ithread), 1)
!                      call zgemv('c', nrho2, nrho2, cunity, zmatder, nrho2, cmtmp1_omp(1,1,ithread), 1, czero, cmtmp7_omp(1,1,ithread), 1)
!                      call zgemv('n', nrho2, nrho2, cunity, zlmtmp1_omp(1,1,ithread), nrho2,    &
!                                 cmtmp7_omp(1,1,ithread), 1, czero, cmtmp1_omp(1,1,ithread), 1)
                      cmtmp5_omp(1:nrho,1:nrho,ithread) = cmtmp5_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp1
                      cmtmp6_omp(1:nrho,1:nrho,ithread) = cmtmp6_omp(1:nrho,1:nrho,ithread) + cmtmp2_omp(1:nrho,1:nrho,ithread) * dtmp2
                      cmtmpa_omp(1:nrho,1:nrho,ithread) = cmtmpa_omp(1:nrho,1:nrho,ithread) + cmtmp1_omp(1:nrho,1:nrho,ithread) * dtmp1
                      cmtmpb_omp(1:nrho,1:nrho,ithread) = cmtmpb_omp(1:nrho,1:nrho,ithread) + cmtmp1_omp(1:nrho,1:nrho,ithread) * dtmp2
                   end do
                end do
                call zamsmm1_omp('l', transc, 'n', iorbs1, ispin, -cunity, cmtmp5_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp3_omp(1,1,ithread), nrho,                                                    &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                call zamsmm1_omp('r', transc, 'n', iorbs1, ispin, -cunity, cmtmp6_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp3_omp(1,1,ithread), nrho,                                                    &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                call zamsmm1_omp('l', transc, 'n', iorbs1, ispin,  cunity, cmtmpa_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp4_omp(1,1,ithread), nrho,                                                    &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
                call zamsmm1_omp('r', transc, 'n', iorbs1, ispin,  cunity, cmtmpb_omp(1,1,ithread), nrho, cunity,  &
                                 cmtmp4_omp(1,1,ithread), nrho,                                                    &
                                 cmtmp7_omp(1,1,ithread), nrho, cmtmp8_omp(1,1,ithread), nrho)
             end do
          end do
       end do
       rhsa(1:nrho,1:nrho,lni) = rhsa(1:nrho,1:nrho,lni) + cmtmp3_omp(1:nrho,1:nrho,ithread)
       rhsb(1:nrho,1:nrho,lni) = rhsb(1:nrho,1:nrho,lni) + cmtmp4_omp(1:nrho,1:nrho,ithread)
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
end subroutine calc_ax_cf_omp
