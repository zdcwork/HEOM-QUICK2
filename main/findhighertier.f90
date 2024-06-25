subroutine findhighertier(itier, kindex, nplus, lout, iop, ffact, rfact)
use matmod
!
!   purpose : goes from {nk} (N) -> {nk, nl+1} (N + 1)
!            
!      note : N = itier - 1, (itier = 2, ..., ntier0)
!    
!     input 
!    kindex : the multi-component index vector to search (sorted) 
!     nplus : the drawer to be added, i.e., nl
!   idirect : 0, if nplus is to be added in front of input kindex
!             1, if nplus is to be appended to the rear of kindex
!         
!    output
!      lout : the starting position of relevant higher-tier counterpart in indextable
!       iop : the type of the outcome of search
!             0, if the {nk, nl+1} is found 
!             1, if the complement of {nk, nl+1} is found , i.e., \bar{nk, nl+1}
!      fact : factor, -1 or 1
!  
! algorithm : binary search algorithm
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)  :: itier, nplus, kindex(*)
integer*1, intent(out) :: iop, ffact, rfact
integer*8, intent(out) :: lout
integer   :: ni, nj, nk, nballs
integer*1 :: fact
integer   :: mindex(MAXTIER), iwork(MAXTIER)
integer*8 :: ipost
!
iop    = -1
nballs = itier - 1
!
if (itier .ge. ntier0) then
    write(6,*)
    write(6,*)' error! input itier in findhighertier ', itier
    stop
end if
!
if (nplus .lt. 0 .or. nplus .gt. nvar) then
 write(6,*)
 write(6,*)' error nplus input for findhighertier ', nplus
 stop
end if
!
! front
!
mindex(1) = nplus
mindex(2:itier) = kindex(1:itier-1)
!
ffact = 1
do ni=1,itier-1
  do nj=ni+1,itier
    if (mindex(ni) .gt. mindex(nj)) then
      nk         = mindex(ni)
      mindex(ni) = mindex(nj)
      mindex(nj) = nk
      ffact      = ffact * (-1)
    endif
  enddo
enddo
!
! rear
!
mindex(1:itier-1) = kindex(1:itier-1)
mindex(itier) = nplus
!
rfact = 1
do ni=1,itier-1
  do nj=ni+1,itier
    if (mindex(ni) .gt. mindex(nj)) then
      nk         = mindex(ni)
      mindex(ni) = mindex(nj)
      mindex(nj) = nk
      rfact      = rfact * (-1)
    endif
  enddo
enddo
!
! try {nk, nl+1} first
!
!
call findindex(ifirst(itier + 1), ilast(itier + 1), nballs + 1, mindex(1), MAXTIER, &
               iwork, ipost, 0, fact)
if (ipost .gt. 0) then
 lout  = ipost
 iop   = 0
 ffact = ffact * fact
 rfact = rfact * fact
 return
end if
!
! search the complement of {nk, nl+1}
!
call findindex(ifirst(itier + 1), ilast(itier + 1), nballs + 1, mindex(1), MAXTIER, &
               iwork, ipost, 1, fact)
if (ipost .gt. 0) then
 lout  = ipost
 iop   = 1
 ffact = ffact * fact
 rfact = rfact * fact
 return
end if
!
! error report if neither {nk, nl+1} nor its complement was found
!
return
!
write(6,*)
write(6,*)' error! findhighertier failed! ', itier, nplus
write(6,*)(kindex(ni), ni=1,itier-1)
write(6,*)(mindex(ni), ni=1,itier)
write(6,*)ipost, iop, ffact, rfact
write(6,*)ifirst(itier+1), ilast(itier+1), nballs+1
write(6,*)mpm(kindex(1)), morbs(kindex(1)), mspin(kindex(1))
write(6,*)mpm(nplus), morbs(nplus), mspin(nplus)
stop
end subroutine findhighertier
