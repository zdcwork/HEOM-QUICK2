subroutine findanytier(itier, kindex, lout, iop, fact)
use matmod
!
!   purpose : goes from {nk} (N) -> {nk'} (N)
!
!      note : N = itier - 1, (itier = 2, ..., ntier) (N is # of balls)
!    
!     input 
!    kindex : the multi-component index vector to search (sorted) 
!
!  WARNING! the code does not check whether kindex is sorted or not
!             
!    output
!      lout : the starting position of relevant lower-tier counterpart in indextable
!       iop : the type of the outcome of search
!             0, if the {nk'} is found 
!             1, if the complement of {nk'} is found , i.e., \bar{nk'}
!
! algorithm : binary search algorithm
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)  :: itier, kindex(*)
integer*1, intent(out) :: iop, fact
integer*8, intent(out) :: lout
integer   :: ni, nj, nballs
integer   :: mindex(MAXTIER), iwork(MAXTIER)
integer*8 :: ipost
!
iop    = -1
nballs = itier - 1
!
if (itier .lt. 1) then
    write(6,*)
    write(6,*)'findanytier: error input itier ', itier
    stop
end if
if (itier .eq. 1) then
    iop  = 0
    lout = 1
    fact = 1
    return
end if
!
mindex(1:nballs) = kindex(1:nballs)
!
! try {nk'} first
!
call findindex(ifirst(itier), ilast(itier), nballs, mindex(1), MAXTIER,  &
               iwork, ipost, 0, fact)
if (ipost .gt. 0) then
    lout = ipost
    iop  = 0
    return
end if
!
! search the complement of {nk'}
! 
call findindex(ifirst(itier), ilast(itier), nballs, mindex(1), MAXTIER,  &
               iwork, ipost, 1, fact)
if (ipost .gt. 0) then
    lout = ipost
    iop  = 1 
    return
end if
!
! error report if neither {nk'} nor its complement was found
!
!return
!
write(6,*)
write(6,*)'findanytier: error! search has failed! ', itier
write(6,*)(kindex(ni), ni=1,nballs)
stop
!
end subroutine findanytier
