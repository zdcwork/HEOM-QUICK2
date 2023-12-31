subroutine solvebicg(init)
use bicgmod
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer, intent(in) :: init
!
integer*8 :: lni, lnj, lnk
integer   :: ni, nj, nk
integer   :: istat, iconvg, iter, nhalf, ierror
integer, parameter :: maxit0 = 100
real*8             :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8             :: cpu1, cpu2, cpu3, cpu4
real*8,  parameter :: dprev1 = 1.d-7, dprev2 = 1.d-5
!
real*8     :: ddotlong, dnormlong, jleft, jright
real*8     :: alphak, betak, res1, res2
external   :: ddotlong, dnormlong
!
call cpu_time(cpu3)
write(6,*)
write(6,*)' entering solvebicg '
call flush(6)
!
ierror = -1
iconvg = 0
iter   = 0
!
! initialize
!
call inibicg
psi(1) = czero
!
zetatmpr(1:nunk) = dble(zeta(1:nunk))
zetatmpi(1:nunk) = dimag(zeta(1:nunk))
phitmpr(1:nunk)  = dble(phi(1:nunk))
phitmpi(1:nunk)  = dimag(phi(1:nunk))
psitmpr(1:nunk)  = dble(psi(1:nunk))
psitmpi(1:nunk)  = dimag(psi(1:nunk))
!
call calc_ax_bicg
!
!------debug
!dtmp1 = 0.d0
!dtmp2 = 0.d0
!dtmp3 = 0.d0
!dtmp4 = 0.d0
!dtmp5 = 0.d0
!dtmp6 = 0.d0
!do lni=1,nunk
! dtmp1 = max(dtmp1, dabs(zetarhsr(lni)))
! dtmp2 = max(dtmp2, dabs(zetarhsi(lni)))
! dtmp3 = max(dtmp3, dabs(phirhsr(lni)))
! dtmp4 = max(dtmp4, dabs(phirhsi(lni)))
! dtmp5 = max(dtmp5, dabs(psirhsr(lni)))
! dtmp6 = max(dtmp6, dabs(psirhsi(lni)))
!end do
!write(6,*)' initialization for bicg, max(A * x) = '
!write(6,*)dtmp1
!write(6,*)dtmp2
!write(6,*)dtmp3
!write(6,*)dtmp4
!write(6,*)dtmp5
!write(6,*)dtmp6
!call flush(6)
!------end of debug
!
rk(1:nunk)          = 0.d0 - zetarhsr(1:nunk)
rk(nunk+1:2*nunk)   = 0.d0 - zetarhsi(1:nunk)
rk(2*nunk+1:3*nunk) = rhsrbicg(1:nunk) - phirhsr(1:nunk)
rk(3*nunk+1:4*nunk) = rhsibicg(1:nunk) - phirhsi(1:nunk)
rk(4*nunk+1:5*nunk) = 0.d0 - psirhsr(1:nunk)
rk(5*nunk+1:6*nunk) = 0.d0 - psirhsi(1:nunk)
!
!---debug AX
!dtmp1 = 0.d0
!do lni=1, 6*nunk
!! dtmp1 = max(dtmp1, dabs(rk(lni)))
! if (dtmp1 .le. dabs(rk(lni))) then
!  dtmp1 = dabs(rk(lni))
!  lnj   = lni
! end if
!end do
!write(6,*)
!write(6,*)' max(rk) = ', dtmp1, lnj
!write(6,*)nunk
!!
!write(6,*)' real'
!do lni=1,nunk
! write(6,*)phirhsr(lni), rhsrbicg(lni), rk(2*nunk + lni)
!enddo
!write(6,*)' imag'
!do lni=1,nunk
! write(6,*)phirhsi(lni), rhsibicg(lni), rk(3*nunk + lni)
!enddo
!write(6,*)' psi  '
!do lni=1,nunk
! write(6,*)psirhsr(lni), psirhsi(lni)
!enddo
!stop
!--- end of debug
!
sk(1:6*nunk) = rk(1:6*nunk)
dk(1:6*nunk) = rk(1:6*nunk)
fk(1:6*nunk) = sk(1:6*nunk)
!
write(6,*)
call cpu_time(cpu1)
100 continue
!
res1 = ddotlong(6*nunk, sk, rk)
!
zetatmpr(1:nunk) = dk(1:nunk)
zetatmpi(1:nunk) = dk(nunk+1:2*nunk)
phitmpr(1:nunk)  = dk(2*nunk+1:3*nunk)
phitmpi(1:nunk)  = dk(3*nunk+1:4*nunk)
psitmpr(1:nunk)  = dk(4*nunk+1:5*nunk)
psitmpi(1:nunk)  = dk(5*nunk+1:6*nunk)
!
call calc_ax_bicg
!
tmp1bicg(1:nunk)          = zetarhsr(1:nunk)
tmp1bicg(nunk+1:2*nunk)   = zetarhsi(1:nunk)
tmp1bicg(2*nunk+1:3*nunk) = phirhsr(1:nunk)
tmp1bicg(3*nunk+1:4*nunk) = phirhsi(1:nunk)
tmp1bicg(4*nunk+1:5*nunk) = psirhsr(1:nunk)
tmp1bicg(5*nunk+1:6*nunk) = psirhsi(1:nunk)
!
res2 = ddotlong(6*nunk, fk, tmp1bicg)
!
alphak = res1 / res2
!
!dtmp1 = 0.d0
!do lni=1,6*nunk
! dtmp1 = max(dtmp1, dabs(alphak * dk(lni)))
!end do
!
!dtmp2 = dnormlong(6*nunk, dk)
!dtmp2 = dtmp2 * dabs(alphak)**2
!
!tmp1bicg(1:nunk)          = dble(zeta(1:nunk))
!tmp1bicg(nunk+1:2*nunk)   = dimag(zeta(1:nunk))
!tmp1bicg(2*nunk+1:3*nunk) = dble(phi(1:nunk))
!tmp1bicg(3*nunk+1:4*nunk) = dimag(phi(1:nunk))
!tmp1bicg(4*nunk+1:5*nunk) = dble(psi(1:nunk))
!tmp1bicg(5*nunk+1:6*nunk) = dimag(psi(1:nunk))
!
!dtmp3 = dnormlong(6*nunk, tmp1bicg)
!
dtmp4 = 0.d0
do lni=1,6*nunk
! dtmp4 = max(dtmp4, dabs(rk(lni)))
 if (dtmp4 .le. dabs(rk(lni))) then
  lnj   = lni
  dtmp4 = dabs(rk(lni))
 end if
end do
!
!if (dtmp2 / dtmp3 .le. dprev) iconvg = 1
!if (dtmp2 .le. dprev) iconvg = 1
!if (dtmp4 .le. dprev) iconvg = 1
if (dabs(res1) .le. dprev1 .and. dabs(res2) .le. dprev2) iconvg = 1
!
call cpu_time(cpu2)
!write(6,518)iter, dtmp2 / dtmp3, dtmp1, cpu2 - cpu1
!write(6,518)iter, dtmp4, dble(zeta(1)), res1, res2
write(6,518)iter, res1, res2, dble(zeta(1)), cpu2 - cpu1
!write(6,*)lnj, nunk
call flush(6)
call cpu_time(cpu1)
518 format(I6, 2x, 5(e16.6e2, 2x))
!
if (iconvg .eq. 1) then
 write(6,*)
 write(6,*)' converged'
 call flush(6)
 ierror = 0
 goto 200
end if
!
zeta(1:nunk) = zeta(1:nunk) + alphak * dcmplx(dk(1:nunk), dk(nunk+1:2*nunk))
phi(1:nunk)  = phi(1:nunk)  + alphak * dcmplx(dk(2*nunk+1:3*nunk), dk(3*nunk+1:4*nunk))
psi(1:nunk)  = psi(1:nunk)  + alphak * dcmplx(dk(4*nunk+1:5*nunk), dk(5*nunk+1:6*nunk))
!
tmp1bicg(1:nunk)          = zetarhsr(1:nunk)
tmp1bicg(nunk+1:2*nunk)   = zetarhsi(1:nunk)
tmp1bicg(2*nunk+1:3*nunk) = phirhsr(1:nunk)
tmp1bicg(3*nunk+1:4*nunk) = phirhsi(1:nunk)
tmp1bicg(4*nunk+1:5*nunk) = psirhsr(1:nunk)
tmp1bicg(5*nunk+1:6*nunk) = psirhsi(1:nunk)
!
rkp(1:6*nunk) = rk(1:6*nunk) - alphak * tmp1bicg(1:6*nunk)
!
zetatmpr(1:nunk) = fk(1:nunk)
zetatmpi(1:nunk) = fk(nunk+1:2*nunk)
phitmpr(1:nunk)  = fk(2*nunk+1:3*nunk)
phitmpi(1:nunk)  = fk(3*nunk+1:4*nunk)
psitmpr(1:nunk)  = fk(4*nunk+1:5*nunk)
psitmpi(1:nunk)  = fk(5*nunk+1:6*nunk)
!
call calc_adagx_bicg
!
tmp1bicg(1:nunk)          = zetarhsr(1:nunk)
tmp1bicg(nunk+1:2*nunk)   = zetarhsi(1:nunk)
tmp1bicg(2*nunk+1:3*nunk) = phirhsr(1:nunk)
tmp1bicg(3*nunk+1:4*nunk) = phirhsi(1:nunk)
tmp1bicg(4*nunk+1:5*nunk) = psirhsr(1:nunk)
tmp1bicg(5*nunk+1:6*nunk) = psirhsi(1:nunk)
!
skp(1:6*nunk) = sk(1:6*nunk) - alphak * tmp1bicg(1:6*nunk)
!
res1 = ddotlong(6*nunk, skp, rkp)
res2 = ddotlong(6*nunk, sk, rk)
!
betak = res1 / res2
!
dk(1:6*nunk) = rkp(1:6*nunk) + betak * dk(1:6*nunk)
fk(1:6*nunk) = skp(1:6*nunk) + betak * fk(1:6*nunk)
!
! regular update k index
!
rk(1:6*nunk) = rkp(1:6*nunk)
sk(1:6*nunk) = skp(1:6*nunk)
!
iter = iter + 1
if (iter .gt. maxit0) then
 write(6,*)
 write(6,*)' error! maxit0 reached without convergency ', iter
! stop
 goto 200
end if
!
goto 100
!
200 continue
!
jleft  = 0.d0
jright = 0.d0
nhalf = nvar / 2
do ni=1,ncor
 jleft  = jleft  - 2.d0 * dble(phi(nfirst(2) - 1 + nhalf + ni))
 jright = jright - 2.d0 * dble(phi(nfirst(2) - 1 + nhalf + ncor + ni))
end do
jleft  = jleft  * jconst
jright = jright * jconst 
!
write(6,*)
write(6,*)' steady state solved by BICG method '
write(6,*)
write(6,*)' occupation ', dble(zeta(1)), dimag(zeta(1))
write(6,*)' psi(1)     ', psi(1)
write(6,*)' jleft   =  ', jleft , ' nA '
write(6,*)' jright  =  ', jright, ' nA '
call cpu_time(cpu4)
!
psi(1) = cunity
!
write(6,*)
write(6,*)' leaving solvebicg ', cpu4 - cpu3
call flush(6)
end subroutine solvebicg
