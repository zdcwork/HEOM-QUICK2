subroutine assigndrawer(nballs, nmem, mtable, mlen)
use matmod
implicit none
!
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)  :: nmem, nballs
integer*8, intent(out) :: mlen
integer,   intent(out) :: mtable(*)          ! output mlen*nballs elements
!
integer :: ni, nj, nk, nl, istat
integer, allocatable :: listncor(:,:), ipoint(:)
integer*8 :: ilist, nlist, lni
!
nlist = nmem**nballs
!
allocate(ipoint(nballs), STAT=istat)
allocate(listncor(nballs,nlist), STAT=istat)
!
ilist = 0
ipoint(1:nballs-1) = 1
ipoint(nballs)     = 0
point1: do
  nj = nballs
  check1: do
    ipoint(nj) = ipoint(nj) + 1
    if ( ipoint(nj) .le. nmem ) then
      exit check1
    else
      nj = nj - 1
      if (nj >= 1) then
        cycle check1
      else
        exit point1
      end if
    end if
  end do check1
  ipoint(nj+1:nballs) = 1 
  ilist = ilist + 1
  listncor(1:nballs,ilist) = ipoint(1:nballs)
end do point1
!
do ilist=1,nlist
   lni = (ilist - 1) * nballs
   mtable(lni+1:lni+nballs) = listncor(1:nballs,ilist)
end do
mlen = nlist
!
!write(6,*)
!write(6,*)'assigndrawer: nballs = ', nballs, 'nmem = ', nmem
!write(6,*)'assigndrawer: mlen   = ', mlen
!do ilist=1,mlen
!   lni = (ilist - 1) * nballs
!   write(6,'(I8, 1x, 10(I3, 1x))') ilist, (mtable(lni+ni), ni=1,nballs)
!end do
!call flush(6)
!
deallocate(ipoint, listncor, STAT=istat)
!
return
end subroutine assigndrawer
!
!---------
subroutine buildcombination(nballs, npick, mtable, mlen)
implicit none
!
integer, intent(in)    :: nballs, npick
integer*8, intent(out) :: mlen
integer, intent(out)   :: mtable(*)
integer                :: istat, ni, nj
integer*8              :: ilist, nlist, ifrac, ifrac2
external               :: ifrac, ifrac2
integer, allocatable   :: ipoint(:)
!
! mlen = C(nballs, npick)
!
mlen = ifrac2(nballs, npick) / ifrac(nballs - npick)
!
allocate(ipoint(npick), STAT=istat)
!
nlist = 0
do ni=1,npick
   ipoint(ni) = ni
end do
ipoint(npick) = npick - 1
!
point1: do
  nj = npick
  check1: do
    ipoint(nj) = ipoint(nj) + 1
    if ( ipoint(nj) .le. nballs-npick+nj ) then
      exit check1
    else
      nj = nj - 1
      if (nj >= 1) then
        cycle check1
      else
        exit point1
      end if
    end if
  end do check1
!
  do ni=nj+1,npick
     ipoint(ni) = ipoint(ni-1) + 1 
  end do
  nlist = nlist + 1
  ilist = (nlist - 1) * npick
  mtable(ilist+1:ilist+npick) = ipoint(1:npick)
end do point1
!
if (nlist .ne. mlen) then
   write(6,*)'buildcombination: error! nlist = ', nlist, ' mlen = ', mlen
   stop
end if
!write(6,*)
!write(6,*)'buildcombination: nballs = ', nballs, 'npick = ', npick
!write(6,*)'buildcombination: mlen   = ', mlen
!do nlist=1,mlen
!   ilist = (nlist - 1) * npick
!   write(6,'(I8, 1x, 10(I3, 1x))') nlist, (mtable(ilist+ni), ni=1,npick)
!end do
!call flush(6)
!
deallocate(ipoint, STAT=istat)
!
return
end subroutine buildcombination

