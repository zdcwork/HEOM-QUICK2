subroutine refresh_rhosys_spa
use matmod
use sparsemod
implicit none
!
! refresh rhosys from rho_spa
!
include '../include/sizes' 
include '../include/common'
!
integer*8  :: lni
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'refresh_rhosys_spa: error entry ', lsparse
   stop
end if
!
lni = ind_spa(1)
call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rho_spa(lni), rhosys, nrho, nrho, nrho)
!
return
end subroutine refresh_rhosys_spa
