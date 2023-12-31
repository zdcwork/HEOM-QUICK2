subroutine findlowertier(itier, kindex, nomit, lout, iop, fact)
use matmod
!
!   purpose : goes from {nk} (N) -> {nk, nl-1} (N - 1)
!
!      note : N = itier - 1, (itier = 2, ..., ntier0) (N is # of balls)
!    
!     input 
!    kindex : the multi-component index vector to search (sorted) 
!     nomit : the drawer to be omitted, i.e., nl
!         
!    output
!      lout : the starting position of relevant lower-tier counterpart in indextable
!       iop : the type of the outcome of search
!             0, if the {nk, nl-1} is found 
!             1, if the complement of {nk, nl-1} is found , i.e., \bar{nk, nl-1}
!
! algorithm : binary search algorithm
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)  :: itier, nomit, kindex(*)
integer*1, intent(out) :: iop, fact
integer*8, intent(out) :: lout
integer   :: ni, nj, nballs
integer   :: mindex(MAXTIER), iwork(MAXTIER)
integer*8 :: ipost
!
iop    = -1
nballs = itier - 1
!
if (itier .le. 1) then
 write(6,*)
 write(6,*)' error input itier in findlowertier ', itier
 stop
end if
if (nomit .lt. 1 .or. nomit .gt. nballs) then
 write(6,*)
 write(6,*)' error nomit input for findlowertier ', nomit
 stop
end if
!
! quick return if possible
!
if (itier .eq. 2) then  ! lower tier contains no ball
 lout = 1
 iop  = 0
 return
end if
!
nj = 0
do ni=1,itier-1
  if (ni .eq. nomit) cycle
  nj = nj + 1
  mindex(nj) = kindex(ni)
end do 
!
! try {nk, nl-1} first
!
call findindex(ifirst(itier - 1), ilast(itier - 1), nballs - 1, mindex(1), MAXTIER, &
               iwork, ipost, 0, fact)
if (ipost .gt. 0) then
 lout = ipost
 iop  = 0
 return
end if
!
! search the complement of {nk, nl-1}
! 
call findindex(ifirst(itier - 1), ilast(itier - 1), nballs - 1, mindex(1), MAXTIER, &
               iwork, ipost, 1, fact)
if (ipost .gt. 0) then
 lout = ipost
 iop  = 1 
 return
end if
!
! error report if neither {nk, nl-1} nor its complement was found
!
!return
!
write(6,*)
write(6,*)' error! findlowertier fail! ', itier, nomit
write(6,*)(kindex(ni), ni=1,itier-1)
write(6,*)(mindex(ni), ni=1,nballs-1)
!write(6,*)ifirst(itier-1), ilast(itier-1), nballs-1
!
!write(6,*)
!write(6,*)'------------------------'
!do ipost=ifirst(itier-1), ilast(itier-1)-(nballs-1), nballs-1
!  write(6,*)(indextable(ipost+ni-1), ni=1,nballs-1), ' | ', &
!            ( iredex(indextable(ipost+ni-1)), ni=1,nballs-1 )
!end do
stop
!
end subroutine findlowertier
