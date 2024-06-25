subroutine allocate_bicg(ierror, memo)
use bicgmod
use bicgmod_omp
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: ierror, memo
integer              :: istat, dim1
real*8               :: dtmp1, dtmp2
integer*8            :: ltmp1
!
ierror = -1
dtmp1  = dble(nunk * 16 * nrho**2) / dble(1024**2)
memobicg = 0.d0
!
allocate(rcgrhs(nrho,nrho,nunk), STAT=istat)
ierror   = istat
memobicg = memobicg + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 1 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
write(6,*)
if (dtmp1 .le. memmax) then
  memo = 0
  write(6,*)'allocate_bicg: sk, rk, dk, fk stored in RAM '
else
  memo = 1
  write(6,*)'allocate_bicg: size of rho is too large to put in RAM, abort...'
  stop
end if
call flush(6)
!
ltmp1 = nunk * nrho**2 
!
if (memo .eq. 0) then
  allocate(sk(nrho,nrho,nunk), rk(nrho,nrho,nunk), dk(nrho,nrho,nunk), fk(nrho,nrho,nunk), STAT=istat)
  ierror = istat
  memobicg = memobicg + dble(ltmp1 * 16) * 4.d0 / dble(1024**2)
end if
!
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 2 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
call allocate_aux
!
call allocate_tmpmem_omp
call allocate_tmplioumat_omp
!   allocate(rcgrhs_omp(nrho,nrho,nunk,nthreads_omp), STAT=istat)
!   memobicg = memobicg + dble(ltmp1 * 16 * nthreads_omp) / dble(1024**2)
!
if (ierror .eq. 0) then
   write(6,*)
   write(6,*)' BICG memory allocated ', memobicg,          ' MB '
   write(6,*)' Total memory used     ', memobicg + memory, ' MB '
   call flush(6)
end if
!
return
end subroutine allocate_bicg
!----------------
!
subroutine free_bicg(memo)
use bicgmod
use bicgmod_omp
use matmod
use matmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: memo
integer :: istat
!
deallocate(rcgrhs, STAT=istat)
deallocate(sk, rk, dk, fk, STAT=istat)
!
call free_aux
!
call free_tmpmem_omp
call free_tmplioumat_omp
call free_tmpauxmat_omp
!   deallocate(rcgrhs_omp, STAT=istat)
!
end subroutine free_bicg
