subroutine prelude_ad
use matmod
use hbmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer :: ni, nj, nk, istat, ncount, fact
integer :: iorbs, ispin, ialf, isgn, icor, iorbs2, idraw, isopr, iopr
logical :: lswap, lsame
real*8  :: dcor
real*8  :: dtmp1, dtmp2
!
integer, allocatable :: nvec_tmp1(:,:), nvec_tmp2(:)
real*8,  allocatable :: dcb_drawer(:)
integer*8 :: ifrac, lni, lnj
external ifrac
!
if (.not. lad) then
    write(6,*)
    write(6,*)'prelude_ad: error entry in prelude_ad ', lad
    stop
end if
!
if (ltrun_der) then
    write(6,*)
    write(6,*)'prelude_ad: error! lad=T and ltrun_der=T ', lad, ltrun_der
    write(6,*)'prelude_ad: code unavailable yet '
    stop
end if
if (lsimple) then
    write(6,*)
    write(6,*)'prelude_ad: error! lad=T and lsimple=T ', lad, lsimple
    write(6,*)'prelude_ad: code unavailable yet '
    stop
end if
if (lhb) then
    write(6,*)
    write(6,*)'prelude_ad: error! lad=T and lhb=T ', lad, lhb
    write(6,*)'prelude_ad: code unavailable yet '
    stop
end if
!if (.not. psdfff .or. numfff .le. 0) then
!    write(6,*)
!    write(6,*)'prelude_ad: error! lad=T but psdfff=F ', lad, psdfff
!    write(6,*)'prelude_ad: numfff = ', numfff
!    stop
!end if
if (.not. lsparse) then
    write(6,*)
    write(6,*)'prelude_ad: error! lad=T but lsparse=F ', lad, lsparse
    write(6,*)' Adiabatic HEOM now only works in SPARSE mode '
    stop
end if
if (lscale) then
    write(6,*)
    write(6,*)'prelude_ad: error! lad=T and lscale=T ', lad, lscale
    write(6,*)'prelude_ad: code unavailable yet '
    stop
end if
!
! ntier_ad should not be lower than ntier
if (ntier_ad .lt. ntier) then 
    write(6,*)
!    write(6,*)'prelude_ad: error! ntier_ad < ntier found ', ntier_ad, ntier
!    stop
    write(6,*)'prelude_ad: warning! ntier_ad < ntier found ', ntier_ad, ntier
    write(6,*)'prelude_ad: continue by setting ntier_ad = ntier '
    ntier_ad = ntier
    call flush(6)
end if
ntier_slow = ntier_ad - ntier
!
if (dgama_slow .le. 0.d0) then
    write(6,*)
    write(6,*)'prelude_ad: error! dgama_slow is non-positive ', dgama_slow
    stop
end if
!
allocate(ns2f(nvar), STAT=istat)
allocate(dgama_drawer(nvar), STAT=istat)
allocate(dcb_drawer(nvar), STAT=istat)
!
write(6,*)
write(6,*)'prelude_ad: Adiabatic approximation is invoked for '
write(6,*)'            itier from ', ntier, ' to ', ntier_ad 
if (lad_fast) then
    write(6,*)
    write(6,*)'prelude_ad: lad_fast=T found '
    write(6,*)'prelude_ad: save ADOs with all fast modes in upper tier '
end if
write(6,*)
write(6,*)'prelude_ad: number of all drawers ', nvar
!
if (lscba_ad) then
    write(6,*)
    write(6,*)'prelude_ad: lscba_ad=T found '
    write(6,*)'prelude_ad: use self-consistent Born approximation '
end if
!
if (lcop_ad) then
    write(6,*)
    write(6,*)'prelude_ad: lcop_ad=T '
    write(6,*)'prelude_ad: use COP (time-nonlocal) form to separate slow and fast modes'
else  
    write(6,*)
    write(6,*)'prelude_ad: lcop_ad=F '
    write(6,*)'prelude_ad: use POP (time-local) form to separate slow and fast modes'
end if
call flush(6)
!
do nk=1,nvar
   isgn  = mpm(nk)
   ialf  = ilead(nk)
   iorbs = morbs(nk)
   ispin = mspin(nk)
   icor  = mcor(nk)
   if (offcor) then
       iorbs2 = (iorbs - 1) * norbs + iorbs
   else
       iorbs2 = iorbs
   end if
   dgama_drawer(nk) = -dble(cgama(iorbs2,ispin,icor,ialf,isgn))
   dcb_drawer(nk)   = cdabs(   cb(iorbs2,ispin,icor,ialf,isgn))
   if (dgama_drawer(nk) .lt. 0.d0) then
       write(6,*)'prelude_ad: error! negative dissipation rate found '
       write(6,*)' drawer, rate ', nk, dgama_drawer(nk) 
       stop
   end if
end do
!
ncount = 0
do nk=1,nvar
   if (dgama_drawer(nk) .le. dgama_slow) then
       ncount = ncount + 1
   end if
end do
ndrawer_slow = ncount
!
! sort all dgama, save the order (from slow to fast) in array ns2f
!
do nk=1,nvar
   ns2f(nk) = nk
end do
do ni=1,nvar-1
   do nj=ni+1,nvar
      lswap = .false.
      if ( dgama_drawer(ns2f(ni)) .gt. dgama_drawer(ns2f(nj)) ) lswap = .true.       ! smaller rate is slower 
      if ( dabs(dgama_drawer(ns2f(ni)) - dgama_drawer(ns2f(nj))) .lt. dpico ) then
          if ( mmfff(ns2f(ni)) .lt. mmfff(ns2f(nj)) ) lswap = .true.                 ! same rate, higher power in t is slower 
      end if
      if (lswap) then
          nk       = ns2f(ni)
          ns2f(ni) = ns2f(nj)
          ns2f(nj) = nk
      end if
   end do
end do
!
!if (lad_fast .and. lset_fast) then
if (lset_fast) then
    if (ncor_fast .lt. 0 .or. ncor_fast .gt. ncor) then
        write(6,*)
        write(6,*)'prelude_ad: error with lset_fast=T ', lset_fast
        write(6,*)'prelude_ad:            ncor_fast=  ', ncor_fast
        stop
    end if
    ncor_slow    = ncor - ncor_fast
    ndrawer_slow = ncor_slow * nopr
    if (ndrawer_slow .eq. nvar) then
        dgama_slow = dgama_drawer(ns2f(nvar)) * 10.d0
    else if (ndrawer_slow .eq. 0) then
        dgama_slow = dgama_drawer(ns2f(1)) / 10.d0
    else
        dgama_slow = 5.d-1 * (dgama_drawer(ns2f(ndrawer_slow)) + dgama_drawer(ns2f(ndrawer_slow+1)))
    end if
    write(6,*)
    write(6,*)'prelude_ad: lset_fast=T, set dgama_slow= ', dgama_slow*hbar
    call flush(6)
end if
!
dtmp1 = dble(ndrawer_slow) / dble(nsopr)
dcor  = dtmp1 / dble(nalf)
write(6,*)
write(6,*)'prelude_ad: threshold rate for the slow modes '
write(6,*)' dgama_slow =               ', dgama_slow*hbar
write(6,*)' number of drawers (modes)  ', nvar
write(6,*)' number of slow modes       ', ndrawer_slow
write(6,*)' number of system operators ', nsopr
write(6,*)' slow modes per operator    ', dtmp1
write(6,*)'            per reservoir   ', dcor 
write(6,*)
write(6,*)'prelude_ad: print rates of all modes from slow to fast '
write(6,*)' ni   idrawer   dgama   power   isopr   ialf           '
do ni=1,nvar
   idraw = ns2f(ni)
   isopr = msopr(idraw)
   ialf  = ilead(idraw)
   write(6,1001)ni, idraw, dgama_drawer(idraw)*hbar, mmfff(idraw)-1, isopr, ialf, dcb_drawer(idraw)
end do
call flush(6)
1001 format(I4, 1x, I4, 1x, e14.6e3, 1x, 3(I2, 1x), 1x, e14.6e3)
!
! The present code presumes every operator iopr=(isgn,iorbs,ispin,ialf) 
! has the same number of slow modes, such a number is set to ncor_slow.
! Otherwise the present code does not execute adiabatic HEOM calculation.
!
! check whether each iopr is associated with the same ncor_slow
!
lsame = .true.
check0: do iopr=1,nopr
   ncount = 0
   do nk=1,ndrawer_slow
      if (mopr(ns2f(nk)) .eq. iopr) then
          ncount = ncount + 1
      end if
   end do
   if ( dabs(dcor - dble(ncount)) .gt. dpico ) then
       lsame = .false.
       write(6,*)
       write(6,*)'prelude_ad: dcor != dble(ncount) ', dcor, ncount
       exit check0
   end if
end do check0
!
if (lsame) then
    ncor_slow = int(dcor)
else
    write(6,*)
    write(6,*)'prelude_ad: error! slow modes are not even among iopr '
    write(6,*)'prelude_ad: code unavailable yet '
    stop
end if
ncor_fast = ncor - ncor_slow
!
write(6,*)
write(6,*)'prelude_ad: number of slow modes (ncor_slow) = ', ncor_slow
write(6,*)'prelude_ad: number of fast modes (ncor_fast) = ', ncor_fast
call flush(6)
!if (ncor_slow .le. 0) then
if (ncor_slow .lt. 0) then
    write(6,*)
    !write(6,*)'prelude_ad: error! ncor_slow <= 0 found ', ncor_slow
    write(6,*)'prelude_ad: error! ncor_slow < 0 found ', ncor_slow
    stop
end if
!
! re-order all the slow modes (drawers) so that for the same
! dissipation rate, the lowest power in t appears first,
! this is to facilitate the subsequent construction of hierarchy
! (for generating functionals of polynormial-exponential form)
!
do ni=1,ndrawer_slow-1
   do nj=ni+1,ndrawer_slow
      lswap = .false.
      if ( dabs(dgama_drawer(ns2f(ni)) - dgama_drawer(ns2f(nj))) .lt. dpico ) then
          if ( mmfff(ns2f(ni)) .gt. mmfff(ns2f(nj)) ) lswap = .true.    ! same rate, lower power in t appears earlier
      end if
      if (lswap) then
          nk       = ns2f(ni)
          ns2f(ni) = ns2f(nj)
          ns2f(nj) = nk
      end if
   end do
end do
if (ndrawer_slow .gt. 0) then
    write(6,*)
    write(6,*)'prelude_ad: print rates of all slow modes (re-ordered) '
    write(6,*)' ni   idrawer   dgama   power   isopr   ialf           '
    do ni=1,ndrawer_slow
       idraw = ns2f(ni)
       isopr = msopr(idraw)
       ialf  = ilead(idraw)
       write(6,1001)ni, idraw, dgama_drawer(idraw)*hbar, mmfff(idraw)-1, isopr, ialf
    end do
    call flush(6)
end if
!
! tag all drawers 
! nj = ns2f(ni) : the ni-th slowest mode is nj-th mode (drawer)
! ni = ms2f(nj) : the inverse mapping of ns2f
!                 the nj-th mode (drawer) is the ni-th slowest among all
!
allocate(ms2f(nvar), mslow(nvar), STAT=istat)
do ni=1,nvar                         ! running through all modes from slowest to fastest
   ms2f(ns2f(ni)) = ni
   if (ni .le. ndrawer_slow) then    
       mslow(ns2f(ni)) = 1           ! mslow(idraw) = 1 means idraw-th mode is slow
   else
       mslow(ns2f(ni)) = 0           ! mslow(idraw) = 0 means idraw-th mode is fast
   end if
end do
!
! further check all operators have the same slow modes
!
if (.not. allocated(icor_slow) .and. ncor_slow .gt. 0) allocate(icor_slow(ncor_slow), STAT=istat)
if (.not. allocated(icor_fast) .and. ncor_fast .gt. 0) allocate(icor_fast(ncor_fast), STAT=istat)
!iopr   = 1
!ncount = 0
!do nk=1,ndrawer_slow
!   if (mopr(ns2f(nk)) .eq. iopr) then
!       ncount = ncount + 1
!       if (ncount .le. ncor_slow) then
!           icor_slow(ncount) = mcor(ns2f(nk))
!       end if
!   end if
!end do
ncount = 0
do idraw=1,ncor  ! take iopr=1
   if (mslow(idraw) .eq. 1) then
       ncount = ncount + 1
       icor_slow(ncount) = idraw
   end if
end do
if (ncount .ne. ncor_slow) then
    write(6,*)
    write(6,*)'prelude_ad: error! ncount, ncor_slow ', ncount, ncor_slow
    stop
end if
!
do nk=1,ndrawer_slow
   call check_inum_vec(mcor(ns2f(nk)), icor_slow, ncor_slow, ni, nj)
   if (nj .ne. 0) then
       write(6,*)
       write(6,*)'prelude_ad: unknown slow mode found '
       write(6,*)'prelude_ad: nk, mcor ', nk, mcor(ns2f(nk))
       stop
   end if
end do
!
ncount = 0
do idraw=1,ncor  ! take iopr=1
   if (mslow(idraw) .eq. 0) then
       ncount = ncount + 1
       icor_fast(ncount) = idraw
   end if
end do
if (ncount .ne. ncor_fast) then
    write(6,*)
    write(6,*)'prelude_ad: error! ncount, ncor_fast ', ncount, ncor_fast
    stop
end if
!
write(6,*)
write(6,*)'prelude_ad: print icor of slow modes '
write(6,*)' ni  icor '
do ni=1,ncor_slow
   write(6,*)ni, icor_slow(ni)
end do
call flush(6)
!
write(6,*)
write(6,*)'prelude_ad: fast modes with the rates between   '
write(6,*)'        max(dgama) and max(dgama)*dratio_fast   ' 
write(6,*)'        are considered to be degenerate,        '
write(6,*)'        with dratio_fast = ', dratio_fast
write(6,*)
if (idegen_fast .eq. 0) then
    write(6,*)'prelude_ad: omit contributions of degenerate fast modes '
else if (idegen_fast .eq. 1) then
    write(6,*)'prelude_ad: symmetrize contributions of degenerate fast modes '
else if (idegen_fast .eq. 2) then
    write(6,*)'prelude_ad: adopt adiabatic approx. for all degenerate modes '
    write(6,*)'prelude_ad: omit contributions if number of degenerate modes '
    write(6,*)'            exceeds ndegen_fast = ', ndegen_fast
    if (ndegen_fast .le. 0) then
        write(6,*)
        write(6,*)'prelude_ad: error! unphysical ndegen_fast = ', ndegen_fast
        stop
    end if
else 
    write(6,*)'prelude_ad: error! unknown idegen_fast = ', idegen_fast
    stop
end if
call flush(6)
!
!!!!!!!!!
!!! test
!!!!!!!!!
!ncount = ntier_ad - 1
!lni    = ifrac(ncount)
!write(6,*)
!write(6,*)'prelude_ad: print all permutations of length ', ncount
!write(6,*)'            number of permutations ', lni
!call flush(6)
!!
!allocate(nvec_tmp1(ncount,lni), nvec_tmp2(lni), STAT=istat)
!do ni=1,ncount
!   nvec_tmp1(ni,1) = ni
!end do
!nvec_tmp2(1) = 1
!lnj = 1
!write(6,1002)lnj, (nvec_tmp1(nj,lnj), nj=1,ncount), nvec_tmp2(lnj)
!do lnj=2,lni
!   call genpermute(ncount, nvec_tmp1(1,lnj-1), nvec_tmp1(1,lnj), nvec_tmp2(lnj))
!   write(6,1002)lnj, (nvec_tmp1(nj,lnj), nj=1,ncount), nvec_tmp2(lnj)
!end do
!call flush(6)
!1002 format(I6, ' | ', 10(1x, I2))
!deallocate(nvec_tmp1, nvec_tmp2, STAT=istat)
!
ncount = ntier
lni    = ifrac(ncount)
!
allocate(npermu(ncount,lni,ncount), npermu_fact(lni,ncount), STAT=istat)
do nk=2,ncount
   lni = ifrac(nk)
   do ni=1,nk
      npermu(ni,1,nk) = ni
   end do
   npermu_fact(1,nk) = 1
   do lnj=2,lni
      call genpermute(nk, npermu(1,lnj-1,nk), npermu(1,lnj,nk), fact)
      npermu_fact(lnj,nk) = npermu_fact(lnj-1,nk) * fact
   end do
end do
!
deallocate(dcb_drawer, STAT=istat)
!
return
!
end subroutine prelude_ad
