subroutine allocate_tmpauxmat_omp
use omp_lib
use matmod_omp
use matmod
use auxmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
allocate(lentmp1_omp(nthreads_omp), STAT=istat)
allocate(lrtmp1_omp(namsmax, nthreads_omp), STAT=istat)
allocate(lctmp1_omp(namsmax, nthreads_omp), STAT=istat)
allocate(lvtmp1_omp(namsmax, nthreads_omp), STAT=istat)
!
end subroutine allocate_tmpauxmat_omp
!
subroutine free_tmpauxmat_omp
use auxmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(lentmp1_omp, lrtmp1_omp, lctmp1_omp, lvtmp1_omp, STAT=istat)
!
end subroutine free_tmpauxmat_omp
