subroutine sortdrawer(mindex, nball, ifact)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
! on exit, the drawers in array "mindex" follows the standard sequence
! # of interchanges of drawers made ---> nswap
! if nswap is an odd integer, the rearrangement results in a (-1.d0) factor (ifact)
! 
integer, intent(in) :: nball
integer, intent(inout) :: mindex(*)
integer, intent(out) :: ifact
!
integer :: ni, nj, ntmp1, ntmp2
!
ifact = 1
!
do ni=1,nball-1
  do nj=ni+1,nball
    call comparedrawer(mindex(ni), mindex(nj), ntmp1)
    if (ntmp1 .lt. 0) then
      ntmp2      = mindex(ni)
      mindex(ni) = mindex(nj)
      mindex(nj) = ntmp2
      ifact      = ifact * (-1)
    else if (ntmp1 .eq. 0) then
      write(6,*)
      write(6,*)' error in <sortdrawer.f90> '
      write(6,*)' identical operators found, drawers are ', mindex(ni), mindex(nj)
      stop
    end if
  end do
end do
!
end subroutine sortdrawer
!
!----------------------------------------------
subroutine comparedrawer(ivar, jvar, iout)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
!             for ivar comes in front of jvar 
!
! iout > 0,   if by standard sequence, the operator associated with ivar appears earlier 
!             than that with jvar, i.e., no swap is needed.
!      = 0,   if ivar and jvar correspond to the same operator.
!      < 0,   if the operator associcated with ivar appears later than that with jvar, 
!             i.e., a swap is needed.
!       
integer, intent(in) :: ivar, jvar
integer, intent(out) :: iout
!
!if (mpm(ivar) .eq. mpm(jvar) .and. morbs(ivar) .eq. morbs(jvar) .and.  &
!  mspin(ivar) .eq. mspin(jvar)) then
!  iout = 0
!  return
!end if
!
if (mpm(ivar) .gt. mpm(jvar)) then
  iout = -1  
  return
else if (mpm(ivar) .lt. mpm(jvar)) then
  iout = 1
  return
end if
!
if (mspin(ivar) .gt. mpm(jvar)) then
  iout = -1
  return
else if (mspin(ivar) .lt. mspin(jvar)) then
  iout = 1
  return
end if
!
if (morbs(ivar) .gt. morbs(jvar)) then
  iout = -1
  return
else if (morbs(ivar) .lt. morbs(jvar)) then
  iout = 1
  return
end if
!
iout = 0
return
!
end subroutine comparedrawer
