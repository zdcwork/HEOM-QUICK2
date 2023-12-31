subroutine intstp(nstep, tt, jleft, jright)
use matmod
implicit none
!
include './include/sizes'
include './include/common'
!
integer, intent(in) :: nstep
real*8,  intent(in) :: tt
real*8,  intent(out) :: jleft, jright
!
integer :: ni, nj, nk
integer :: nhalf
real*8 :: dt2, dt3, dt6
!
dt2 = dt / 2.d0
dt3 = dt / 3.d0
dt6 = dt / 6.d0
!
zeta0(1:nunk) = zeta(1:nunk)
phi0(1:nunk)  = phi(1:nunk)
psi0(1:nunk)  = psi(1:nunk)
!
zetatmp(1:nunk) = zeta0(1:nunk)
phitmp(1:nunk)  = phi0(1:nunk)
psitmp(1:nunk)  = psi0(1:nunk)
!
! [1st] tt
!
call calcderiv(tt)
!write(6,*)'here1'
!
! evaluate current
!
jleft  = 0.d0
jright = 0.d0
!
! sum over phi(alpha, (-), n; itier=1)
! 
nhalf = nvar / 2
!
do ni=1,ncor
 jleft  = jleft  - 2.d0 * dble(phi(nfirst(2) - 1 + nhalf + ni))
 jright = jright - 2.d0 * dble(phi(nfirst(2) - 1 + nhalf + ncor + ni))
end do
jleft  = jleft  * jconst
jright = jright * jconst 
!
zetatmp(1:nunk) = zeta0(1:nunk) + dt2 * zetarhs(1:nunk)
phitmp(1:nunk)  = phi0(1:nunk)  + dt2 * phirhs(1:nunk)
psitmp(1:nunk)  = psi0(1:nunk)  + dt2 * psirhs(1:nunk)
!
zeta(1:nunk) = zeta(1:nunk) + dt6 * zetarhs(1:nunk)
phi(1:nunk)  = phi(1:nunk)  + dt6 * phirhs(1:nunk)
psi(1:nunk)  = psi(1:nunk)  + dt6 * psirhs(1:nunk)
!
! [2nd] tt + dt2
! 
call calcderiv(tt + dt2)
!write(6,*)'here2'
!
zetatmp(1:nunk) = zeta0(1:nunk) + dt2 * zetarhs(1:nunk)
phitmp(1:nunk)  = phi0(1:nunk)  + dt2 * phirhs(1:nunk)
psitmp(1:nunk)  = psi0(1:nunk)  + dt2 * psirhs(1:nunk)
!
zeta(1:nunk) = zeta(1:nunk) + dt3 * zetarhs(1:nunk)
phi(1:nunk)  = phi(1:nunk)  + dt3 * phirhs(1:nunk)
psi(1:nunk)  = psi(1:nunk)  + dt3 * psirhs(1:nunk)
!
! [3rd] tt + dt2
! 
call calcderiv(tt + dt2)
!write(6,*)'here3'
!
zetatmp(1:nunk) = zeta0(1:nunk) + dt * zetarhs(1:nunk)
phitmp(1:nunk)  = phi0(1:nunk)  + dt * phirhs(1:nunk)
psitmp(1:nunk)  = psi0(1:nunk)  + dt * psirhs(1:nunk)
!
zeta(1:nunk) = zeta(1:nunk) + dt3 * zetarhs(1:nunk)
phi(1:nunk)  = phi(1:nunk)  + dt3 * phirhs(1:nunk)
psi(1:nunk)  = psi(1:nunk)  + dt3 * psirhs(1:nunk)
!
! [4th] tt + dt
!
call calcderiv(tt + dt)
!write(6,*)'here4'
!
zeta(1:nunk) = zeta(1:nunk) + dt6 * zetarhs(1:nunk)
phi(1:nunk)  = phi(1:nunk)  + dt6 * phirhs(1:nunk)
psi(1:nunk)  = psi(1:nunk)  + dt6 * psirhs(1:nunk)
!
end subroutine intstp
