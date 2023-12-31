subroutine findoprhigher(itier, kindex, nplus, lout, iop)
use matmod
!
!   purpose : goes from {nk} (N) -> {nk, nl+1} (N + 1)
!            
!      note : N = itier - 1, (itier = 2, ..., ntier)
!    
!     input 
!    kindex : the multi-component index vector to search (sorted) 
!     nplus : the drawer to be added, i.e., nl
!         
!    output
!      lout : the starting position of relevant higher-tier counterpart in oprtable
!       iop : the type of the outcome of search
!              0, if the {nk, nl+1} is found 
!              1, if the complement of {nk, nl+1} is found , i.e., \bar{nk, nl+1}
!             -1, if neither 0 or 1 is satisfied. 
!  
! algorithm : binary search algorithm
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)  :: itier, nplus, kindex(*)
integer*1, intent(out) :: iop
integer*8, intent(out) :: lout
integer   :: ni, nj, nk, nballs
integer   :: mindex(MAXTIER), iwork(MAXTIER)
integer*8 :: ipost
!
iop    = -1
nballs = itier - 1
!
if (itier .ge. ntier) then
 write(6,*)
 write(6,*)' error! input itier in findoprhigher ', itier, ntier
 stop
end if
!
if (nplus .lt. 0 .or. nplus .gt. nopr) then
 write(6,*)
 write(6,*)' error nplus input for findoprhigher ', nplus, nopr
 stop
end if
!
mindex(1) = nplus
mindex(2:itier) = kindex(1:itier-1)
!
! sort
!
do ni=1,itier-1
  do nj=ni+1,itier
    if (mindex(ni) .gt. mindex(nj)) then
      nk         = mindex(ni)
      mindex(ni) = mindex(nj)
      mindex(nj) = nk
    endif
  enddo
enddo
!
!if (itier .eq. 2) then
!  write(6,*) ' mindex ', (mindex(ni), ni=1,itier)
!end if
!
! try {nk, nl+1} first
!
call findoprindex(ioprfirst(itier+1), ioprlast(itier+1), nballs+1, mindex(1), MAXTIER, &
                  iwork, ipost, 0)
if (ipost .gt. 0) then
 lout  = ipost
 iop   = 0
 return
end if
!
! search the complement of {nk, nl+1}
!
call findoprindex(ioprfirst(itier+1), ioprlast(itier+1), nballs+1, mindex(1), MAXTIER, &
                  iwork, ipost, 1)
if (ipost .gt. 0) then
 lout  = ipost
 iop   = 1
 return
end if
!
! not found
!
lout = -1
iop  = -1 
return
!
end subroutine findoprhigher
