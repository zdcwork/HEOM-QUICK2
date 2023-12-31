subroutine refresh_rhosys
use matmod
implicit none
!
! refresh rhosys from rho
!
include '../include/sizes' 
include '../include/common'
!
if (lsparse) then
   call refresh_rhosys_spa
   return
end if
!
rhosys(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
!
return
end subroutine refresh_rhosys
