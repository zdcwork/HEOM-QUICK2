subroutine allocate_rk4(ierror)
use matmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: ierror
integer              :: istat
integer*8            :: lnum
real*8               :: dtmp1
!
ierror = -1
if (lsparse) then
   lnum = lunk_spa
else
   lnum = nunk * nrho**2 
end if
dtmp1 = dble(lnum * 16) / dble(1024**2)
!
memork4 = 0.d0
!
if (lsparse) then
   allocate(rho0_spa(lunk_spa), STAT=istat)
else
   allocate(rho0(nrho,nrho,nunk), STAT=istat)
end if
ierror = istat
memork4 = memork4 + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 1 '
 write(6,*)' memory desired ', memork4, ' MB  '
 stop
end if
!
if (lsparse) then
   allocate(rhotmp_spa(lunk_spa), STAT=istat)
else
   allocate(rhotmp(nrho,nrho,nunk), STAT=istat)
end if
ierror = istat
memork4 = memork4 + dtmp1
if (ierror .ne. 0) then
 write(6,*)
 write(6,*)' error! memory allocation fail 2 '
 write(6,*)' memory desired ', memork4, ' MB  '
 stop
end if
!
if (lsparse) then
   allocate(rhorhs_spa(lunk_spa), STAT=istat)
else
   allocate(rhorhs(nrho,nrho,nunk), STAT=istat)
end if
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
return
end subroutine allocate_rk4
!
!----------------------------------
subroutine free_rk4
use matmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
if (lsparse) then
   deallocate(rho0_spa, rhotmp_spa, rhorhs_spa, STAT=istat)
   deallocate(irowtd_hsys, icoltd_hsys, cvaltd_hsys, STAT=istat)
else
   deallocate(rho0, rhotmp, rhorhs, STAT=istat)
end if
!
end subroutine free_rk4
