subroutine allocate_rk4(ierror)
use matmod
implicit none
!
include './include/sizes'
include './include/common'
!
integer, intent(out) :: ierror
integer :: istat
real*8  :: dtmp1
!
ierror = -1
dtmp1  = dble(nunk * 3 * 16) / dble(1024**2)
!
memork4 = 0.d0
!
allocate(zeta0(nunk), phi0(nunk), psi0(nunk), STAT=istat)
ierror = istat
memork4 = memork4 + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 1 '
 write(6,*)' memory desired ', memork4, ' MB  '
 stop
end if
!
allocate(zetatmp(nunk), phitmp(nunk), psitmp(nunk), STAT=istat)
ierror = istat
memork4 = memork4 + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 2 '
 write(6,*)' memory desired ', memork4, ' MB  '
 stop
end if
!
allocate(zetarhs(nunk), phirhs(nunk), psirhs(nunk), STAT=istat)
ierror = istat
memork4 = memork4 + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 3 '
 write(6,*)' memory desired ', memork4, ' MB  '
 stop
end if
!
if (ierror .eq. 0) then
 write(6,*)
 write(6,*)' RK4 memory allocated ', memork4,          ' MB '
 write(6,*)' Total memory used    ', memork4 + memory, ' MB '
 call flush(6)
end if
!
end subroutine allocate_rk4
!
!----------------------------------
subroutine free_rk4
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer :: istat
!
deallocate(zeta0, phi0, psi0, STAT=istat)
deallocate(zetatmp, phitmp, psitmp, STAT=istat)
deallocate(zetarhs, phirhs, psirhs, STAT=istat)
!
end subroutine free_rk4
