subroutine intstp(nstep, tt, esys, esb, etot, jleftu, jleftd, jrightu, jrightd)
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)  :: nstep
real*8,  intent(in)  :: tt
real*8,  intent(out) :: jleftu, jleftd, jrightu, jrightd, etot, esb, esys
!
integer   :: ni, nj, nk
integer*8 :: lni
real*8    :: dt2, dt3, dt6
!
dt2 = dt / 2.d0
dt3 = dt / 3.d0
dt6 = dt / 6.d0
!
if (lsparse) then
   rho0_spa  (1:lunk_spa) = rho_spa(1:lunk_spa)
   rhotmp_spa(1:lunk_spa) = rho_spa(1:lunk_spa)
else 
   rho0  (1:nrho,1:nrho,1:nunk) = rho(1:nrho,1:nrho,1:nunk)
   rhotmp(1:nrho,1:nrho,1:nunk) = rho(1:nrho,1:nrho,1:nunk)
end if
!
! [1st] tt
!
if (lsparse) then
   call calcderiv_spa_omp(tt)
else
   call calcderiv_omp(tt)
end if
!
! evaluate current
!
call calcj(jleftu, jrightu, jleftd, jrightd)
call calcenergy(esys, esb, etot)
!
! calculate information compressibility
!
if (lrun_loc) then
   call calc_compress(dinfcom_loc, energy_loc, entropy_loc)
end if
!
if (lsparse) then
   rhotmp_spa(1:lunk_spa) = rho0_spa(1:lunk_spa) + dt2 * rhorhs_spa(1:lunk_spa)
   rho_spa   (1:lunk_spa) = rho_spa (1:lunk_spa) + dt6 * rhorhs_spa(1:lunk_spa)
else
   rhotmp(1:nrho,1:nrho,1:nunk) = rho0(1:nrho,1:nrho,1:nunk) + dt2 * rhorhs(1:nrho,1:nrho,1:nunk)  
   rho   (1:nrho,1:nrho,1:nunk) = rho (1:nrho,1:nrho,1:nunk) + dt6 * rhorhs(1:nrho,1:nrho,1:nunk)
end if
!
! [2nd] tt + dt2
! 
if (lsparse) then
   call calcderiv_spa_omp(tt + dt2)
else
   call calcderiv_omp(tt + dt2)
end if
!
if (lsparse) then
   rhotmp_spa(1:lunk_spa) = rho0_spa(1:lunk_spa) + dt2 * rhorhs_spa(1:lunk_spa)
   rho_spa   (1:lunk_spa) = rho_spa (1:lunk_spa) + dt3 * rhorhs_spa(1:lunk_spa)
else
   rhotmp(1:nrho,1:nrho,1:nunk) = rho0(1:nrho,1:nrho,1:nunk) + dt2 * rhorhs(1:nrho,1:nrho,1:nunk) 
   rho   (1:nrho,1:nrho,1:nunk) = rho (1:nrho,1:nrho,1:nunk) + dt3 * rhorhs(1:nrho,1:nrho,1:nunk)
end if
!
! [3rd] tt + dt2
! 
if (lsparse) then
   call calcderiv_spa_omp(tt + dt2)
else
   call calcderiv_omp(tt + dt2)
end if
!
if (lsparse) then
   rhotmp_spa(1:lunk_spa) = rho0_spa(1:lunk_spa) + dt  * rhorhs_spa(1:lunk_spa)
   rho_spa   (1:lunk_spa) = rho_spa (1:lunk_spa) + dt3 * rhorhs_spa(1:lunk_spa)
else
   rhotmp(1:nrho,1:nrho,1:nunk) = rho0(1:nrho,1:nrho,1:nunk) + dt  * rhorhs(1:nrho,1:nrho,1:nunk)  
   rho   (1:nrho,1:nrho,1:nunk) = rho (1:nrho,1:nrho,1:nunk) + dt3 * rhorhs(1:nrho,1:nrho,1:nunk)
end if
!
! [4th] tt + dt
!
if (lsparse) then
   call calcderiv_spa_omp(tt + dt)
else
   call calcderiv_omp(tt + dt)
end if
!
if (lsparse) then
   rho_spa(1:lunk_spa) = rho_spa(1:lunk_spa) + dt6 * rhorhs_spa(1:lunk_spa)
else
   rho(1:nrho,1:nrho,1:nunk) = rho(1:nrho,1:nrho,1:nunk) + dt6 * rhorhs(1:nrho,1:nrho,1:nunk)
end if
!
end subroutine intstp
