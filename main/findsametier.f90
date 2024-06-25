subroutine findsametier(itier, kindex, nomit, ndraw, lout, iop, fact)
use matmod
!
!   purpose : goes from {nk} (N) -> {nk'} (N)
!             same tier, replace nomit with nadd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   WARNING : the code is only correct when ndraw < kindex(nomit)    
!             it is NOT universally correct!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      note : N = itier - 1, (itier = 2, ..., ntier) (N is # of balls)
!    
!     input 
!    kindex : the multi-component index vector to search (sorted) 
!     nomit : the ball to omit, whose drawer index is kindex(nomit) 
!     ndraw : the draw index to be added 
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
integer,   intent(in)  :: itier, nomit, ndraw, kindex(*)
integer*1, intent(out) :: iop, fact
integer*8, intent(out) :: lout
integer   :: ni, nj, nballs
integer*1 :: fact1, fact2
integer   :: mindex(MAXTIER), iwork(MAXTIER)
integer*8 :: ipost
!
iop    = -1
nballs = itier - 1
!
if (itier .le. 1) then
    write(6,*)
    write(6,*)'findsametier: error input itier ', itier
    stop
end if
if (nomit .lt. 1 .or. nomit .gt. nballs) then
    write(6,*)
    write(6,*)'findsametier: error nomit ', nomit
    stop
end if
if (ndraw .lt. 1 .or. ndraw .gt. nvar) then
    write(6,*)
    write(6,*)'findsametier: error ndraw ', ndraw
    stop
end if
!
mindex(1:nballs) = kindex(1:nballs)
!
fact1 = 1
if (nomit .eq. 1) then
    mindex(1) = ndraw
else
    do ni=nomit-1,1,-1
       mindex(ni+1) = ndraw 
       if (mindex(ni) .gt. ndraw) then
          mindex(ni+1) = mindex(ni) 
          fact1 = -fact1
       else
          exit
       end if
       if (ni .eq. 1) mindex(ni) = ndraw
    end do
end if
!
! try {nk'} first
!
call findindex(ifirst(itier), ilast(itier), nballs, mindex(1), MAXTIER,  &
               iwork, ipost, 0, fact2)
if (ipost .gt. 0) then
    lout = ipost
    iop  = 0
    fact = fact1 * fact2
    return
end if
!
! search the complement of {nk'}
! 
call findindex(ifirst(itier), ilast(itier), nballs, mindex(1), MAXTIER,  &
               iwork, ipost, 1, fact2)
if (ipost .gt. 0) then
    lout = ipost
    iop  = 1 
    fact = fact1 * fact2
    return
end if
!
! error report if neither {nk'} nor its complement was found
!
!return
!
write(6,*)
write(6,*)'findsametier: error! search has failed! ', itier, nomit, ndraw
write(6,*)(kindex(ni), ni=1,nballs)
write(6,*)(mindex(ni), ni=1,nballs)
stop
!
end subroutine findsametier
