module diismod
implicit none
!
! diis memory
!
complex*16,   allocatable :: xnew(:,:,:), xold(:,:,:)
complex*16,   allocatable :: rhs(:,:,:)
!
end module diismod
