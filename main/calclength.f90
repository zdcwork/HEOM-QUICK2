subroutine calclength(itier, len0, len2) 
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: itier
integer*8, intent(out) :: len0, len2
integer :: nballs
!
! purpose : 
!     (itier - 1) indistinguishable balls fall into (nvar) different drawers
!     nballs = itier - 1   
! output 
!     len : total number of possible cases
!
integer*8 :: ifrac, ifrac2
!
nballs = itier - 1
if (itier .eq. 1) then
  len2 = 1
else
  len2 = ifrac2(nopr + nballs - 1, max(nballs, nopr - 1)) / ifrac(min(nballs, nopr - 1))
end if
if (itier .gt. ntier) then
    len0 = len2 * ncor_slow**nballs
    if (lad_fast) then
        len0 = len2 * ncor_fast**nballs
    end if
else
    len0 = len2 * ncor**nballs
end if
!
end subroutine calclength
!
!---------------------------------------------
integer*8 function ifrac(nin)
! 
! purpose : ifrac = (nin)!
!
implicit none
integer, intent(in) :: nin
integer :: ni
!
ifrac = 1
do ni=1,nin
 ifrac = ifrac * ni
end do
end function ifrac
!---------------------------------------------
integer*8 function ifrac2(n1, n2)
!
! purpose : ifrac2 = (n1)! / (n2)!
!
implicit none
integer, intent(in) :: n1, n2
integer :: ni
!
ifrac2 = 1
do ni=n2+1, n1
 ifrac2 = ifrac2 * ni
end do
end function ifrac2

