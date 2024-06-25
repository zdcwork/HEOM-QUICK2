subroutine allocate_tfqmr(istat)
use tfqmrmod
use matmod
use matmod_omp
use sparsemod
use hbmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: istat
real*8               :: dtmp1
integer*8            :: ltmp1
!
if (lsparse) then
   ltmp1 = lunk_spa * 2        ! both real and imaginary parts
   ltmp1 = ltmp1 + 2           ! tail element (complex) for normalization of rho0
else
   ltmp1 = nunk * nrho**2 * 2  ! both real and imaginary parts
   if (lhb) then
      ltmp1 = ltmp1 + nunk * nrho**2 * nbath_hb * 2
   end if
end if
!
dtmp1     = dble(ltmp1 * 8) / dble(1024**2)
memotfqmr = 0.d0
!
allocate(qmrvecs(ltmp1,9), STAT=istat)
memotfqmr = memotfqmr + dtmp1 * 9.d0
if (istat .ne. 0) then
   write(6,*)'allocate_tfqmr: error! memory allocation fail 1 '
   write(6,*)'memory desired ', memotfqmr, ' MB'
   stop
end if
!
call allocate_tmpmem_omp
!
write(6,*)
write(6,*)' TFQMR memory allocated ', memotfqmr,          ' MB '
write(6,*)' Total memory used      ', memotfqmr + memory, ' MB '
call flush(6)
!
return
end subroutine allocate_tfqmr
!----------------
!
subroutine free_tfqmr
use tfqmrmod
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(qmrvecs, STAT=istat)
!
call free_tmpmem_omp
!
end subroutine free_tfqmr
