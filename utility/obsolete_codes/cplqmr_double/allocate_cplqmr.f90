subroutine allocate_cplqmr(istat)
use cplqmrmod
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: istat
real*8               :: dtmp1
integer*8            :: ltmp1
!
dtmp1      = dble(nunk * nrho**2 * 16) / dble(1024**2)
memocplqmr = 0.d0
!
ltmp1 = nunk * nrho**2 * 2    ! both real and imaginary parts
!
! setting dimension
!
nlim  = maxit0 / 2
maxpq = 3
maxvw = maxpq
mvec  = max0(maxvw, maxpq) + 2
!dimm  = (maxpq + maxvw + 2) * 2
dimm  = maxpq + maxvw + 2
!
nrow_dwk  = dimm
ncol_dwk  = 8 * dimm + 18
nrow_idx  = 6
ncol_idx  = nlim + 2
nrow_iwk  = dimm
ncol_iwk  = 13
ncol_vecs = 5 * mvec + 3
!
allocate(qmrvecs(ltmp1,ncol_vecs), STAT=istat)
memocplqmr = memocplqmr + dtmp1 * dble(ncol_vecs)
if (istat .ne. 0) then
   write(6,*)'allocate_cplqmr: error! memory allocation fail 1 '
   write(6,*)'memory desired ', memocplqmr, ' MB'
   stop
end if
!
allocate(iwk(nrow_iwk,ncol_iwk), STAT=istat)
allocate(idx(nrow_idx,ncol_idx), STAT=istat)
allocate(dwk(nrow_dwk,ncol_dwk), STAT=istat)
!
call allocate_aux
!
if (nthreads_omp > 1) then
   call allocate_tmpmem_omp
end if
!
write(6,*)
write(6,*)' CPLQMR memory allocated ', memocplqmr,          ' MB '
write(6,*)'  Total memory used      ', memocplqmr + memory, ' MB '
call flush(6)
end subroutine allocate_cplqmr
!----------------
!
subroutine free_cplqmr
use cplqmrmod
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(qmrvecs, STAT=istat)
deallocate(iwk, idx, dwk, STAT=istat)
!
call free_aux
!
if (nthreads_omp > 1) then
   call free_tmpmem_omp
end if
!
end subroutine free_cplqmr
