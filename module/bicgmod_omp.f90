module bicgmod_omp
implicit none
!
! openmp auxiliary memory for bicg
!
complex*16, save, allocatable :: rcgrhs_omp(:,:,:,:)
!
end module bicgmod_omp
