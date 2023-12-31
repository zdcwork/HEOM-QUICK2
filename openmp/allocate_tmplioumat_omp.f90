subroutine allocate_tmplioumat_omp
use omp_lib
use matmod_omp
use matmod
use tmplioumat_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat, ndim
!
ndim = nrho * nrho
!
allocate(dlvec1_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec2_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec3_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec4_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec5_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec6_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec7_omp(ndim,nthreads_omp), STAT=istat)
allocate(dlvec8_omp(ndim,nthreads_omp), STAT=istat)
dlvec1_omp = 0.d0
dlvec2_omp = 0.d0
dlvec3_omp = 0.d0
dlvec4_omp = 0.d0
dlvec5_omp = 0.d0
dlvec6_omp = 0.d0
dlvec7_omp = 0.d0
dlvec8_omp = 0.d0
!
end subroutine allocate_tmplioumat_omp
!
subroutine free_tmplioumat_omp
use tmplioumat_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(dlvec1_omp, dlvec2_omp, dlvec3_omp, dlvec4_omp,  &
           dlvec5_omp, dlvec6_omp, dlvec7_omp, dlvec8_omp,  &
           STAT=istat)
!
end subroutine free_tmplioumat_omp
