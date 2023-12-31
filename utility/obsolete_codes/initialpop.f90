subroutine initialpop(itrun, iresi)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: itrun, iresi
!
! Initialize upop and vpop
!
! [1] dot(upop) = -i [Hs(t), upop]
!       upop(0) = c_ms (or c^dag_ms)
!
upop(1:nrho, 1:nrho, 1:norbs, 1:nspin) = dcmplx(amsall(1:nrho, 1:nrho, 1:norbs, 1:nspin), 0.d0)
!
! [2] dot(vpop) = exp( (gamma_k + i*dipm(k)*engyshift(k)) * t ) * upop
!       vpop(0) = 0  
!
! here vpop and/or vb, vd should have already been given proper values, either via
! ground state calculation or read from input file.
!
!if (itrun .eq. 1) then
!  vpop = czero
!end if
!
!if (iresi .eq. 1) then
!  vb = czero
!  vd = czero
!end if
!
end subroutine initialpop
