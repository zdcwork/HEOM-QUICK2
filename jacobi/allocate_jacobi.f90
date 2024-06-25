subroutine allocate_jacobi(istat)
use jacobimod
use matmod
use matmod_omp
use sparsemod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: istat
integer              :: nrho2
real*8               :: dtmp1
integer*8            :: ltmp1
!
nrho2 = nrho**2
if (nrho2 .gt. ndimbig) then
    write(6,*)
    write(6,*)'allocate_jacobi: error! nrho2 > ndimbig is found ', nrho2, ndimbig
    stop
end if
if (.not. allocated(zlmtmp1)) allocate(zlmtmp1(nrho2,nrho2), STAT=istat)
if (.not. allocated(zlmtmp2)) allocate(zlmtmp2(nrho2,nrho2), STAT=istat)
!
if (lsparse) then
   ltmp1 = lunk_spa * 2        ! both real and imaginary parts
   ltmp1 = ltmp1 + 2
else
   ltmp1 = nunk * nrho**2 * 2  ! both real and imaginary parts
end if
!
dtmp1 = dble(ltmp1 * 8) / dble(1024**2)
memojacobi = 0.d0
!
allocate(jacvecs(ltmp1,2), STAT=istat)
memojacobi = memojacobi + dtmp1 * 2.d0
if (istat .ne. 0) then
   write(6,*)'allocate_jacobi: error! memory allocation fail 1 '
   write(6,*)'memory desired ', memojacobi, ' MB'
   stop
end if
!
call allocate_tmpmem_omp
!
write(6,*)
write(6,*)' jacobi memory allocated ', memojacobi,          ' MB '
write(6,*)' Total memory used       ', memojacobi + memory, ' MB '
call flush(6)
!
return
end subroutine allocate_jacobi
!----------------
!
subroutine free_jacobi
use jacobimod
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(jacvecs, STAT=istat)
!
call free_tmpmem_omp
!
end subroutine free_jacobi
