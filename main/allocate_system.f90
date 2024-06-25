subroutine allocate_system
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
allocate(hs(nrho,nrho), STAT=istat)
allocate(dhs(nrho,nrho), STAT=istat) 
allocate(rhosys(nrho,nrho), STAT=istat)
!
return
end subroutine allocate_system
