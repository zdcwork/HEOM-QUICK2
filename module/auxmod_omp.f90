module auxmod_omp
implicit none
!
! sparse matrix (liouville space)
!
integer, save, allocatable :: lentmp1_omp(:)
integer, save, allocatable :: lrtmp1_omp(:,:), lctmp1_omp(:,:)
real*8,  save, allocatable :: lvtmp1_omp(:,:)
!
end module auxmod_omp
