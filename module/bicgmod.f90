module bicgmod
implicit none
!
! bicg memory
!
complex*16, save,  allocatable :: sk(:,:,:), rk(:,:,:), dk(:,:,:), fk(:,:,:)
complex*16, save,  allocatable :: rcgrhs(:,:,:), rcgtmp(:,:,:)
!
end module bicgmod
