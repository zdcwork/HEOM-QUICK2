integer function imodms(inp)
implicit none
include '../include/sizes'
include '../include/common'
integer, intent(in) :: inp
!
if (inp .le. 0 .or. inp .gt. nvar) then
 write(6,*)
 write(6,*)' error! wrong input for imodms ', inp, nvar
 stop
end if
!
if (mod(inp, norbs * nspin) .ne. 0) then
 imodms = inp / (norbs * nspin) + 1 
 return
else 
 imodms = inp / (norbs * nspin)
 return
end if
end function imodms
!
subroutine jmod(linear, dim0, irow, icol)
implicit none
!
! purpose : linear = (icol - 1) * dim0 + irow
!           linear -> (irow, icol) 
!
integer, intent(in)  :: linear, dim0
integer, intent(out) :: irow, icol
integer              :: ntmp1
!
if (linear .le. 0 .or. dim0 .le. 0) then
  write(6,*)
  write(6,*)' error in jmod ', linear, dim0
  stop
end if
!
ntmp1 = mod(linear, dim0)
if (ntmp1 .eq. 0) then
    irow = dim0
    icol = linear / dim0 
else
    irow = ntmp1
    icol = linear / dim0 + 1
end if
!
end subroutine jmod
