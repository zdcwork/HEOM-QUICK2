subroutine allocate_hb(istat)
use matmod
use matmod_omp
use hbmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(out) :: istat
real*8               :: dtmp1
integer*8            :: ltmp1
!
if (.not. lhb) then
   write(6,*)'allocate_hb: error! wrong entry ', lhb
   stop
end if
!
if (lsparse) then
   ltmp1 = lunk_hb_spa * 2                ! both real and imaginary parts
   ! to complete: lunk_hb_spa = ? 
else
   ltmp1 = nunk * nrho**2 * nbath_hb * 2  ! both real and imaginary parts
end if
!
dtmp1  = dble(ltmp1 * 8) / dble(1024**2)
memohb = 0.d0
!
if (lsparse) then
   ! to complete 
else                
   ! this part to be included in  <allocate_rk4.f90> ?
   allocate(rho_hb   (nrho,nrho,nunk,nbath_hb), STAT=istat)
!   allocate(rho_hb   (nrho,nrho,nunk,nbath_hb),    &
!            rho0_hb  (nrho,nrho,nunk,nbath_hb),    &
!            rhotmp_hb(nrho,nrho,nunk,nbath_hb),    &
!            rhorhs_hb(nrho,nrho,nunk,nbath_hb), STAT=istat)
   rho_hb(1:nrho,1:nrho,1:nunk,1:nbath_hb) = czero
!   memohb = memohb + dtmp1 * 4.d0
   memohb = memohb + dtmp1 
   if (istat .ne. 0) then
      write(6,*)'allocate_hb: error! memory allocation fail 2 '
      write(6,*)'memory desired ', memohb, ' MB'
      stop
   end if
end if
!
write(6,*)
write(6,*)' Heat bath memory allocated ', memohb,          ' MB '
write(6,*)' Total memory used          ', memohb + memory, ' MB '
call flush(6)
!
return
end subroutine allocate_hb
