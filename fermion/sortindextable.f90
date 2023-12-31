subroutine sortindextable
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, allocatable :: mindex(:), kindex(:), mwork(:)
integer :: ni, nj, nk, nballs, itier, istat, nblock
integer*8 :: lni, lnj
real*8 :: cpu1, cpu2
!
write(6,*)
write(6,*)' entering sortindextable...        '
write(6,*)' this could be very time-consuming '
call flush(6)
!
nk = (ntier - 1) * nvar0
allocate(mindex(nk), kindex(nk), mwork(nk), STAT=istat)
do itier=2,ntier
  call cpu_time(cpu1)
  nballs = itier - 1
  nblock = nballs * nvar0
  write(6,*)
  write(6,*)' nballs= ', nballs, ' nblock= ', nblock
  call flush(6)
  do lni=ifirst(itier), ilast(itier)-2*nblock+1, nblock
    mindex(1:nblock) = indextable(lni:lni+nblock-1)
    do lnj=lni+nblock, ilast(itier)-nblock+1, nblock
      if ( mindex(1) .ge. indextable(lnj) ) then
        kindex(1:nblock) = indextable(lnj:lnj+nblock-1)
        call compare4sort(nballs, nblock, mindex(1), kindex(1), mwork(1))
        indextable(lni:lni+nblock-1) = mindex(1:nblock)
        indextable(lnj:lnj+nblock-1) = kindex(1:nblock)
      end if
    end do
  end do
  call cpu_time(cpu2)
  write(6,*)' indextable on tier ', itier, ' is sorted by ', cpu2 - cpu1, ' secs'
  call flush(6)
end do
deallocate(mindex, kindex, mwork, STAT=istat)
!
end subroutine sortindextable
!
subroutine compare4sort(nballs, nblock, mindex, kindex, mwork)
implicit none
integer, intent(in) :: nballs, nblock
integer, intent(inout) :: mindex(*), kindex(*), mwork(*)
integer :: ni, nj
logical :: lswap
!
lswap = .false.
comp: do ni=1,nballs
  if (mindex(ni) .gt. kindex(ni)) then
    lswap = .true. 
    exit comp
  else if (mindex(ni) .lt. kindex(ni)) then
    exit comp
  else 
    cycle comp
  end if
end do comp
!
if (lswap) then
  do ni=1,nblock
    nj = mindex(ni) 
    mindex(ni) = kindex(ni)
    kindex(ni) = nj
  end do
!  mwork (1:nblock) = mindex(1:nblock)
!  mindex(1:nblock) = kindex(1:nblock)
!  kindex(1:nblock) = mwork (1:nblock)
end if
end subroutine compare4sort
