subroutine calcc12p(dout, rhoinp)
use matmod
use tmpmatmod
implicit none
!
include '../include/sizes'
include '../include/common'
!
! Calculate concurrence of two dots (norbs=2, nspin=2) C12
! The calculation follows Eqs.(15)-(16) of the following paper
! Quantum Information and Computation, Vol. 1, pp 27-44 (2001)
! This part of code is contributed by Shikuan Wang
!
real*8, intent(out) :: dout
complex*16, intent(in) :: rhoinp(nrho,*)
integer :: ni, istat, ndim
complex*16, allocatable :: rhoct(:,:)
!
if (norbs .ne. 2 .or. nspin .ne. 2) then
   write(6,*)'calc12p: error entrance, norbs, nspin ', norbs, nspin
   stop
end if
!
ndim = 4
allocate(rhoct(ndim,ndim), STAT=istat)
!
rhoct(1:ndim,1:ndim) = czero
rhoct(1,1) = rhoinp(6  ,6)
rhoct(1,2) = rhoinp(6  ,7)
rhoct(1,3) = rhoinp(6  ,10)
rhoct(1,4) = rhoinp(6  ,11)
rhoct(2,1) = rhoinp(7  ,6)
rhoct(2,2) = rhoinp(7  ,7)
rhoct(2,3) = rhoinp(7  ,10)
rhoct(2,4) = rhoinp(7  ,11)
rhoct(3,1) = rhoinp(10 ,6)
rhoct(3,2) = rhoinp(10 ,7)
rhoct(3,3) = rhoinp(10 ,10)
rhoct(3,4) = rhoinp(10 ,11)
rhoct(4,1) = rhoinp(11 ,6)
rhoct(4,2) = rhoinp(11 ,7)
rhoct(4,3) = rhoinp(11 ,10)
rhoct(4,4) = rhoinp(11 ,11)
!
call calcconcurr(ndim, rhoct, dout)
!
deallocate(rhoct, STAT=istat)
!
return
end subroutine calcc12p
