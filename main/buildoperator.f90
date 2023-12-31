subroutine buildoperator
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: nunity, istat
integer :: ni, nj, nk, ispin
integer :: irowtmp, icoltmp, isgntmp
integer, allocatable :: irow(:,:), icol(:,:), isgn(:,:)
real*8 :: cpu1, cpu2
!
!write(6,*)
!write(6,*)' entering buildoperator '
!call cpu_time(cpu1)
!
call calcams(norbs, nspin, 1, 1, nrho, ams)
amsori(1:nrho, 1:nrho, 1, 1) = ams(1:nrho, 1:nrho)
!
if (nspin * norbs .eq. 1) return
!
if (nspin .eq. 1) then ! spinless system
 nunity = 2**(norbs - 1)
 allocate(irow(nunity,nspin), icol(nunity,nspin), isgn(nunity,nspin), STAT=istat)
 nk = 0
 do ni=1,nrho
  do nj=1,nrho
   if (dabs(ams(ni,nj) - 1.d0) .le. dnano) then
    nk = nk + 1
    irow(nk,1) = ni
    icol(nk,1) = nj
    isgn(nk,1) = int(ams(ni,nj))
   end if 
  end do
 end do
 if (nk .ne. nunity) then
  write(6,*)
  write(6,*)' error! wrong c1 '
  call amatout(nrho, nrho, ams)
  stop
 end if
!
 do ni=2,norbs
  dmtmp1(1:nrho, 1:nrho) = 0.d0
  do nj=1,nunity
   call swapms(nrho, norbs, nspin, irow(nj,1), icol(nj,1), ni, 1, &
               irowtmp, icoltmp, nk)
   dmtmp1(irowtmp, icoltmp) = dble(nk) * dble(isgn(nj,1))
   irow(nj,1) = irowtmp
   icol(nj,1) = icoltmp
   isgn(nj,1) = int(dmtmp1(irowtmp, icoltmp))
  end do
  amsori(1:nrho, 1:nrho, ni, 1) = dmtmp1(1:nrho, 1:nrho)
 end do
 deallocate(irow, icol, isgn, STAT=istat)
!
else
!  
 nunity = 2 * 4**(norbs - 1) 
 allocate(irow(nunity,nspin), icol(nunity,nspin), isgn(nunity,nspin), STAT=istat)
 nk = 0
 do ni=1,nrho 
  do nj=1,nrho
   if (dabs(ams(ni,nj) - 1.d0) .le. dnano) then
    nk = nk + 1
    irow(nk,1) = ni
    icol(nk,1) = nj
    isgn(nk,1) = int(ams(ni,nj))
   end if
  end do
 end do
 if (nk .ne. nunity) then
  write(6,*)
  write(6,*)' error! wrong c1u '
  call amatout(nrho, nrho, ams)
  stop
 end if
!
 dmtmp1(1:nrho, 1:nrho) = 0.d0
 do nj=1,nunity
  call swapspin(nrho, norbs, irow(nj,1), icol(nj,1), 1, irowtmp, &
                icoltmp, isgntmp)
  irow(nj,2) = irowtmp
  icol(nj,2) = icoltmp
  isgn(nj,2) = isgntmp * isgn(nj,1)
  dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,2))
 end do
 amsori(1:nrho, 1:nrho, 1, 2) = dmtmp1(1:nrho, 1:nrho)
!
 ispin=2
 do ni=2,norbs
   dmtmp1(1:nrho, 1:nrho) = 0.d0
   do nj=1,nunity
     call swapms(nrho, norbs, nspin, irow(nj,ispin), icol(nj,ispin), ni,         &
                 ispin, irowtmp, icoltmp, isgntmp)
     dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,ispin)) * dble(isgntmp)
     irow(nj,ispin) = irowtmp
     icol(nj,ispin) = icoltmp 
     isgn(nj,ispin) = int(dmtmp1(irowtmp, icoltmp))
   end do
   amsori(1:nrho, 1:nrho, ni, ispin) = dmtmp1(1:nrho, 1:nrho) 
!
   ispin = 3 - ispin
!
   dmtmp1(1:nrho, 1:nrho) = 0.d0
   do nj=1,nunity 
     call swapspin(nrho, norbs, irow(nj,3-ispin), icol(nj,3-ispin), ni, irowtmp, &
                   icoltmp, isgntmp)
     dmtmp1(irowtmp, icoltmp) = dble(isgn(nj,3-ispin)) * dble(isgntmp)
     irow(nj,ispin) = irowtmp
     icol(nj,ispin) = icoltmp
     isgn(nj,ispin) = int(dmtmp1(irowtmp, icoltmp))
   end do
   amsori(1:nrho, 1:nrho, ni, ispin) = dmtmp1(1:nrho, 1:nrho) 
 end do

 deallocate(irow, icol, isgn, STAT=istat)
!
end if
!
call cpu_time(cpu2)
!write(6,*)
!write(6,*)' leaving buildoperator ', cpu2 - cpu1
!call flush(6)
end subroutine buildoperator
!
!------------------------------------------------------------------
subroutine printams(iorbs, ispin)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: iorbs, ispin
!
if (iorbs .le. 0 .or. iorbs .gt. norbs) then
 write(6,*)
 write(6,*)' error! wrong input iorbs for printams ', iorbs
 stop
end if
if (nspin .eq. 2 .and. ispin .ne. 1 .and. ispin .ne. 2) then
 write(6,*)
 write(6,*)' error! wrong input ispin for printams ', ispin, nspin
 stop
end if
dmtmp1(1:nrho, 1:nrho) = amsori(1:nrho, 1:nrho, iorbs, ispin)
write(6,*)
write(6,*)' iorbs = ', iorbs, ' ispin = ', ispin
call amatout(nrho, nrho, dmtmp1)
!
end subroutine printams
