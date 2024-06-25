subroutine initialize2
use matmod
implicit none
!
include './include/sizes'
include './include/common'
!
integer :: istat
integer*8 :: ntmp1
!
allocate(zeta(nunk), phi(nunk), psi(nunk), STAT=istat)
!
write(6,*)
if (istat==0) then
 write(6,*)' OK! Memory allocation for zeta, phi, and psi successful '
 write(6,*)' memory allocated : ', dble(nunk * 3 * 16) / dble(1024**2), ' MB'
else
 write(6,*)' error! Memory allocation for zeta, phi, and psi failed  '
 write(6,*)' memory needed    : ', dble(nunk * 3 * 16) / dble(1024**2), ' MB'
 stop
end if
call flush(6)
!
open(unit=15, file='variables.data', form='unformatted', status='old')
rewind(15)
read(15)ntmp1
!
if (ntmp1 .ne. nunk) then
 write(6,*)
 write(6,*)' error! bad input <variables.data> ', ntmp1, nunk
 stop
end if
!
read(15)(zeta(ntmp1), ntmp1=1,nunk)
read(15)(phi(ntmp1), ntmp1=1,nunk)
read(15)(psi(ntmp1), ntmp1=1,nunk)
close(15)
!
! psi(1) is actually not an unknown, psi(1) = cunity
!
psi(1) = cunity
!
memory = memory + dble(nunk * 3 * 16) / dble(1024**2)
!
end subroutine initialize2
