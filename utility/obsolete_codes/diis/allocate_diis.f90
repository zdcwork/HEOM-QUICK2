subroutine allocate_diis(ierror)
use diismod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: ierror
integer :: istat
real*8  :: dtmp1
integer*8 :: ltmp1, ltmp2
!
ierror = -1
dtmp1  = dble(nunk * 16 * nrho**2) / dble(1024**2)
memodiis = 0.d0
!
allocate(rhs(nrho,nrho,nunk), STAT=istat)
ierror   = istat
memodiis = memodiis + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 1 '
 write(6,*)' memory desired ', memodiis, ' MB  '
 stop
end if
!
ltmp1 = nunk * nrho**2 
ltmp2 = 2 * ltmp1
allocate(xnew(nrho,nrho,nunk), xold(nrho,nrho,nunk), STAT=istat)
ierror = istat
memodiis = memodiis + dble(ltmp1 * 16) * 2.d0 / dble(1024**2)
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 2 '
 write(6,*)' memory desired ', memodiis, ' MB  '
 stop
end if
!
if (ierror .eq. 0) then
 write(6,*)
 write(6,*)' DIIS memory allocated ', memodiis,          ' MB '
 write(6,*)' Total memory used     ', memodiis + memory, ' MB '
 call flush(6)
end if
end subroutine allocate_diis
!----------------
!
subroutine free_diis
use diismod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(xnew, xold, rhs, STAT=istat)
!
end subroutine free_diis
