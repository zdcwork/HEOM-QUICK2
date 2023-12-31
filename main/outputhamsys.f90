subroutine outputhamsys
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'

integer :: i,j
real*8, allocatable :: tmpmat(:,:)

i=0
j=0
allocate(tmpmat(nrho,nrho))
!
open(unit=78, file='ham_sys.data', status='unknown')
!
write(78,*) nrho 
!write(78,*) nrho**2
tmpmat(1:nrho,1:nrho) = dble(hs(1:nrho,1:nrho)) * hbar
do i=1, nrho
  do j=1, nrho
    write(78,*) tmpmat(i,j)
  end do
end do
!
close(78)
deallocate(tmpmat)

end subroutine outputhamsys
