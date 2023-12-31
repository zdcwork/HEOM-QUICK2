subroutine sortindextable2
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
write(6,*)' entering sortindextable2 ...        '
write(6,*)' this could be very time-consuming   '
call flush(6)
!
nk = nvar0**(ntier - 1) * (ntier - 1)
allocate(mindex(nk), kindex(nk), mwork(ntier*nvar0), STAT=istat)
do itier=2,ntier
  call cpu_time(cpu1)
  nballs = itier - 1
  nblock = nvar0**nballs * nballs
  write(6,*)
  write(6,*)' nballs= ', nballs, ' nblock= ', nblock
  call flush(6)
  do lni=ifirst(itier), ilast(itier)-2*nblock+1, nblock
    mindex(1:nblock) = indextable(lni:lni+nblock-1)
    loopj: do lnj=lni+nblock, ilast(itier)-nblock+1, nblock
      if ( mopr(mindex(1)) .lt. mopr(indextable(lnj)) ) then
        cycle loopj
      else 
        kindex(1:nblock) = indextable(lnj:lnj+nblock-1)
!        call compare4sort2(nballs, nvar0**nballs, nblock, mindex(1), kindex(1), mwork(1))
!        call compare4sort3(nballs, nvar0**nballs, nblock, mindex(1), kindex(1), mwork(1))
        call compare4sort3(nballs*nvar0, nblock/nballs/nvar0, nblock, mindex(1), kindex(1), mwork(1))
        indextable(lni:lni+nblock-1) = mindex(1:nblock)
        indextable(lnj:lnj+nblock-1) = kindex(1:nblock)
      end if
    end do loopj
  end do
  call cpu_time(cpu2)
  write(6,*)' indextable on tier ', itier, ' is sorted by ', cpu2 - cpu1, ' secs'
  call flush(6)
end do
deallocate(mindex, kindex, mwork, STAT=istat)
end subroutine sortindextable2
!
subroutine compare4sort2(nballs, nlist, nblock, mindex, kindex, mwork)
implicit none
integer, intent(in)    :: nballs, nlist, nblock
integer, intent(inout) :: mindex(*), kindex(*), mwork(*)
integer :: ni, nj, nk, nl, np, nq, ilist1, ilist2, ilist3, ilist4, ntmp1
logical :: lswap
nk = 1
loopj: do nj=1,nlist
  ilist2 = (nj - 1) * nballs + 1
  loopi: do ni=nk,nlist
    ilist1 = (ni - 1) * nballs + 1
    call swapornot(nballs, mindex(ilist1), kindex(ilist2), lswap)
    if (lswap) then
      do nl=1,nballs
        mwork(nl)           = mindex(ilist1-1+nl)
        mindex(ilist1-1+nl) = kindex(ilist2-1+nl)
      end do
      nk = ni + 1
! adjust the 2nd block
      do nl=ilist2,nblock-2*nballs+1,nballs
        kindex(nl:nl-1+nballs) = kindex(nl+nballs:nl-1+2*nballs)
      end do 
      kindex(nblock-nballs+1:nblock) = mwork(1:nballs)
      loopp: do np=nlist,nj+1,-1
        ilist3 = (np - 1) * nballs + 1
        ilist4 = ilist3 - nballs
        call swapornot(nballs, kindex(ilist4), kindex(ilist3), lswap)
        if (lswap) then
          do nl=1,nballs
            ntmp1 = kindex(ilist4-1+nl)
            kindex(ilist4-1+nl) = kindex(ilist3-1+nl)
            kindex(ilist3-1+nl) = ntmp1
          end do
        else
          exit loopp 
        end if
      end do loopp
      cycle loopi
    end if
  end do loopi
end do loopj
end subroutine compare4sort2
!
subroutine swapornot(nballs, mindex, kindex, lswap)
implicit none
integer, intent(in)  :: nballs
integer, intent(in)  :: mindex(*), kindex(*)
logical, intent(out) :: lswap
integer              :: ni, nj
lswap = .false.
loopi: do ni=1,nballs
  if (mindex(ni) .gt. kindex(ni)) then
    lswap = .true.
    exit loopi
  else if (mindex(ni) .lt. kindex(ni)) then
    exit loopi
  else
    cycle loopi
  end if
end do loopi
return
end subroutine swapornot
!
subroutine compare4sort3(nballs, nlist, nblock, mindex, kindex, mwork)
implicit none
integer, intent(in)    :: nballs, nlist, nblock
integer, intent(inout) :: mindex(*), kindex(*), mwork(*)
integer :: ni, nj, nk, nl, np, nq, ilist1, ilist2, ilist3, ilist4, ntmp1
integer :: iout1, iout2
logical :: lswap
integer :: iter
data iter /0/
iter = iter + 1
nk = 1
loopj: do nj=1,nlist
  ilist2 = (nj - 1) * nballs + 1
  10 if (nk .gt. nlist) return
  call findpost(nballs, kindex(ilist2), mindex(1), nk, nlist, iout1)
  if (iout1 .le. nlist) then
    ilist1 = (iout1 - 1) * nballs + 1
    do nl=1,nballs
      mwork(nl)           = mindex(ilist1-1+nl)
      mindex(ilist1-1+nl) = kindex(ilist2-1+nl)
    end do
    nk = nk + 1 ! ok
    if (nj .lt. nlist) then
      call findpost(nballs, mwork(1), kindex(1), nj+1, nlist, iout2)
      do np=nj,iout2-2
        ilist3 = (np - 1) * nballs + 1
        do nl=1,nballs
          kindex(ilist3-1+nl) = kindex(ilist3-1+nl+nballs)
        end do
      end do
      ilist4 = (iout2 - 2) * nballs + 1
      do nl=1,nballs
        kindex(ilist4-1+nl) = mwork(nl)
      end do
    else  ! nj=nlist
      do nl=1,nballs
        kindex(ilist2-1+nl) = mwork(nl)
      end do
    end if
    goto 10
  end if
end do loopj
end subroutine compare4sort3
!
subroutine findpost(nballs, recrd1, mindex, nstart, nend, nout)
implicit none
integer,  intent(in)  :: nballs, nstart, nend, recrd1(*), mindex(*)
integer,  intent(out) :: nout
! find the proper position of recrd1 in the list of sorted records (mindex)
integer               :: low, high, ni, nj, ilisti, ilistj, low0, high0
logical               :: lswap
!
ilisti = (nstart - 1) * nballs + 1
call swapornot(nballs, mindex(ilisti), recrd1(1), lswap)
if (lswap) then
  nout = nstart
  return
end if
ilisti = (nend - 1) * nballs + 1
call swapornot(nballs, recrd1(1), mindex(ilisti), lswap)
if (lswap) then
  nout = nend + 1
  return
end if
low  = 1
high = nend - nstart + 1
ni   = (low + high) / 2
low0 = 0
high0 = 0
if (high .lt. low) then
  write(6,*)
  write(6,*)' error! error in findpost! high, low : ', nstart, nend, low, high
  stop
end if
find: do while (low .le. high)
  if (low .eq. low0 .and. high .eq. high0) then
    nout = (nstart - 1) + low + 1
    exit find 
  else
    low0 = low
    high0 = high
  end if
  nj     = nstart - 1 + ni
  ilistj = (nj - 1) * nballs + 1
  call swapornot(nballs, recrd1(1), mindex(ilistj), lswap) ! low < recrd1 <= high
  if (lswap) then ! recrd1 > nj
    low = max(low0, ni)
  else
    high = min(high0, ni)
  end if
  ni = (low + high) / 2
end do find
if (low .gt. high) then
  write(6,*)
  write(6,*)' error in findpost! low > high found ', low, high
  stop
end if
return
end subroutine findpost

