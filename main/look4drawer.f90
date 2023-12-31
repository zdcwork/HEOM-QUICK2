subroutine look4drawer(isgn, ialf, iorbs, ispin, imats, iout)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(in)  :: isgn, ialf, iorbs, ispin, imats
integer, intent(out) :: iout
integer              :: iopr
!
if (imats .lt. 0) then
  write(6,*)
  write(6,*)' error input in look4drawer ', imats
  stop
end if
!
!iout = (isgn - 1) * nalf * ncor * nspin * norbs   &
!       + (ispin - 1) * norbs * nalf * ncor        &
!       + (iorbs - 1) * nalf * ncor                &
!       + (ialf - 1) * ncor                        &
!       + imats 
call look4operator(isgn, iorbs, ispin, ialf, iopr)
iout = (iopr - 1) * ncor + imats
!
return
end subroutine look4drawer
!
!========
!
subroutine look4operator(isgn, iorbs, ispin, ialf, iout)
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(in)  :: isgn, iorbs, ispin, ialf
integer, intent(out) :: iout
integer              :: isopr
!
!iout = (isgn - 1) * nspin * norbs * nalf    &
!       + (ispin - 1) * norbs * nalf         &
!       + (iorbs - 1) * nalf                 &
!       + ialf
!
call look4sysopr(isgn, iorbs, ispin, isopr)
iout = (isopr - 1) * nalf + ialf
!
return
end subroutine look4operator
!

subroutine look4sysopr(isgn, iorbs, ispin, iout)
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(in) :: isgn, iorbs, ispin
integer, intent(out) :: iout
!
iout = (isgn - 1) * nspin * norbs + (ispin - 1) * norbs + iorbs
!
return
end subroutine look4sysopr
