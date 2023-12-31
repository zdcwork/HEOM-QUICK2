subroutine allocate_bicg(ierror)
use bicgmod
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer, intent(out) :: ierror
integer :: istat
real*8  :: dtmp1
!
ierror = -1
dtmp1  = dble(nunk * 8) / dble(1024**2)
memobicg = 0.d0
!
allocate(rhsrbicg(nunk), rhsibicg(nunk), STAT=istat)
ierror   = istat
memobicg = memobicg + dtmp1 * 2.d0
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 1 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
allocate(sk(6*nunk), rk(6*nunk), skp(6*nunk), rkp(6*nunk),  &
         dk(6*nunk), fk(6*nunk), tmp1bicg(6*nunk), STAT=istat)
ierror = istat
memobicg = memobicg + dtmp1 * 6.d0 * 7.d0
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 2 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
allocate(zetarhsr(nunk), phirhsr(nunk), psirhsr(nunk),         &
         zetatmpr(nunk), phitmpr(nunk), psitmpr(nunk), STAT=istat)        
ierror = istat
memobicg = memobicg + dtmp1 * 6.d0
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 3 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
allocate(zetarhsi(nunk), phirhsi(nunk), psirhsi(nunk),         &
         zetatmpi(nunk), phitmpi(nunk), psitmpi(nunk), STAT=istat)        
ierror = istat
memobicg = memobicg + dtmp1 * 6.d0
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 4 '
 write(6,*)' memory desired ', memobicg, ' MB  '
 stop
end if
!
if (ierror .eq. 0) then
 write(6,*)
 write(6,*)' BICG memory allocated ', memobicg,          ' MB '
 write(6,*)' Total memory used     ', memobicg + memory, ' MB '
 call flush(6)
end if
!
end subroutine allocate_bicg
!
!----------------
!
subroutine free_bicg
use bicgmod
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer :: istat
!
deallocate(rhsrbicg, rhsibicg, STAT=istat)
deallocate(sk, rk, skp, rkp, dk, fk, tmp1bicg, STAT=istat)
deallocate(zetarhsr, phirhsr, psirhsr, zetatmpr, phitmpr, psitmpr, STAT=istat)
deallocate(zetarhsi, phirhsi, psirhsi, zetatmpi, phitmpi, psitmpi, STAT=istat)
!
end subroutine free_bicg

