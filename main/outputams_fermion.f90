subroutine outputams_fermion
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
!integer, intent(in)  :: norbs, nspin, nrho
integer :: nunity, istat
integer :: ni, nj, nk, ispin, orbs, spin
integer :: irowtmp, icoltmp, isgntmp
integer, allocatable :: irow(:,:), icol(:,:), isgn(:,:)
real*8 :: cpu1, cpu2
real*8, allocatable :: amstmp(:,:), ams_all(:,:,:,:)
!
allocate(amstmp(nrho,nrho), ams_all(nrho,nrho,norbs,nspin), STAT=istat)
!
call calcams(norbs, nspin, 1, 1, nrho, amstmp)
ams_all(1:nrho, 1:nrho, 1, 1) = amstmp(1:nrho, 1:nrho)
!
if (nspin * norbs .eq. 1) goto 100
!
if (nspin .eq. 1) then ! spinless system
 nunity = 2**(norbs - 1)
 allocate(irow(nunity,nspin), icol(nunity,nspin), isgn(nunity,nspin), STAT=istat)
 nk = 0
 do ni=1,nrho
  do nj=1,nrho
   if (dabs(amstmp(ni,nj) - 1.d0) .le. dnano) then
    nk = nk + 1
    irow(nk,1) = ni
    icol(nk,1) = nj
    isgn(nk,1) = int(amstmp(ni,nj))
   end if 
  end do
 end do
 if (nk .ne. nunity) then
  write(6,*)
  write(6,*)'outputams_fermion: error! wrong c1 ', nk, nunity
  call amatout(nrho, nrho, amstmp)
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
  ams_all(1:nrho, 1:nrho, ni, 1) = dmtmp1(1:nrho, 1:nrho)
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
   if (dabs(amstmp(ni,nj) - 1.d0) .le. dnano) then
    nk = nk + 1
    irow(nk,1) = ni
    icol(nk,1) = nj
    isgn(nk,1) = int(amstmp(ni,nj))
   end if
  end do
 end do
 if (nk .ne. nunity) then
  write(6,*)
  write(6,*)'outputams_fermion: error! wrong c1u ', nk, nunity
  call amatout(nrho, nrho, amstmp)
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
 ams_all(1:nrho, 1:nrho, 1, 2) = dmtmp1(1:nrho, 1:nrho)
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
   ams_all(1:nrho, 1:nrho, ni, ispin) = dmtmp1(1:nrho, 1:nrho) 
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
   ams_all(1:nrho, 1:nrho, ni, ispin) = dmtmp1(1:nrho, 1:nrho) 
 end do

 deallocate(irow, icol, isgn, amstmp, STAT=istat)
!
end if
!
100 open(unit=76, file='ams_fermion.data', status='unknown')
write(76,*) norbs
write(76,*) nspin
!
do orbs=1, norbs
  do spin=1, nspin
    write(76,*) orbs
    write(76,*) spin
    do ni=1, nrho
      do nj=1, nrho
        write(76,200) ams_all(ni,nj,orbs,spin)
      end do
    end do
  end do
end do
!
deallocate(ams_all, STAT=istat)
!
200 format(f8.5)
!
close(76)
!
end subroutine outputams_fermion
