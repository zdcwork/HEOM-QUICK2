subroutine allocate_tfqmr(istat)
use tfqmrmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: istat
real*8               :: dtmp1
integer*8            :: ltmp1
!
dtmp1     = dble(nunk * 16 * nrho**2) / dble(1024**2)
memotfqmr = 0.d0
!
allocate(xqmr(nrho,nrho,nunk), bqmr(nrho,nrho,nunk), STAT=istat)
memotfqmr = memotfqmr + dtmp1 * 2.d0
if (istat .ne. 0) then
   write(6,*)'allocate_tfqmr: error! memory allocation fail 1 '
   write(6,*)'memory desired ', memotfqmr, ' MB'
   stop
end if
!
ltmp1 = nunk * nrho**2 
allocate(qmrvecs(ltmp1,9), STAT=istat)
memotfqmr = memotfqmr + dtmp1 * 9.d0
if (istat .ne. 0) then
   write(6,*)'allocate_tfqmr: error! memory allocation fail 2 '
   write(6,*)'memory desired ', memotfqmr, ' MB'
   stop
end if
!
write(6,*)
write(6,*)' TFQMR memory allocated ', memotfqmr,          ' MB '
write(6,*)' Total memory used      ', memotfqmr + memory, ' MB '
call flush(6)
end subroutine allocate_tfqmr
!----------------
!
subroutine free_tfqmr
use tfqmrmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(xqmr, bqmr, STAT=istat)
deallocate(qmrvecs, STAT=istat)
!
end subroutine free_tfqmr
