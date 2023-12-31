subroutine bicg0(dim1, a, x, b, ierror)
implicit none
!
! Biconjugate gradient method  (BICG)
!
!   A * X = B
!
integer, intent(in) :: dim1
complex*16, intent(in) :: a(dim1, *), b(*)
integer, intent(out) :: ierror
complex*16, intent(out) :: x(dim1)
!
! ierror =  0  : diis successful, normal output.
!          -1  : maxit0 reached, not converged yet
!
integer :: ni, nj, nk
integer :: istat, iconvg, iter
integer, parameter :: maxit0 = 1000
real*8            :: dtmp1, dtmp2, dtmp3, dznrm2
real*8, parameter :: dprev = 1.d-9
!
complex*16 :: alphak, betak, res1, res2
complex*16, allocatable :: sk(:), rk(:), skp(:), rkp(:)
complex*16, allocatable :: dk(:), fk(:), xk(:), yk(:)
complex*16, allocatable :: ctmp1(:)
!
complex*16, parameter :: czero = (0.d0, 0.d0)
complex*16, parameter ::   eye = (0.d0, 1.d0)
complex*16, parameter :: cunity= (1.d0, 0.d0)
!
complex*16 :: zdotc
external zdotc, dznrm2
!
allocate(sk(dim1), rk(dim1), skp(dim1), rkp(dim1), STAT=istat)
!allocate(dk(dim1), fk(dim1), xk(dim1), yk(dim1), STAT=istat)
allocate(dk(dim1), fk(dim1), STAT=istat)
allocate(ctmp1(dim1), STAT=istat)
!
ierror = -1
iconvg = 0
iter   = 0
!
! initialize 
!
x(1:dim1) = czero
!xk(1:dim1) = czero
!yk(1:dim1) = czero
!
rk(1:dim1) = b(1:dim1)
sk(1:dim1) = dconjg(b(1:dim1))
dk(1:dim1) = rk(1:dim1)
fk(1:dim1) = sk(1:dim1)
!
100 continue
!
res1 = zdotc(dim1, sk, 1, rk, 1)
call zgemv('n', dim1, dim1, cunity, a, dim1, dk, 1, czero, ctmp1, 1)
!
!stop
!
res2 = zdotc(dim1, fk, 1, ctmp1, 1)
alphak = res1 / res2           
!
dtmp1 = 0.d0
do ni=1,dim1
 dtmp1 = max(dtmp1, cdabs(alphak * dk(ni)))
end do
x(1:dim1) = x(1:dim1) + alphak * dk(1:dim1)
!xk(1:dim1) = xk(1:dim1) + alphak * dk(1:dim1)
!yk(1:dim1) = yk(1:dim1) + dconjg(alphak) * fk(1:dim1)
!
!dtmp2 = dznrm2(dim1, rk, 1)
dtmp2 = dznrm2(dim1, dk, 1)
dtmp2 = dtmp2 * cdabs(alphak)**2
dtmp3 = dznrm2(dim1, x, 1)
!dtmp3 = dznrm2(dim1, xk, 1)
!
if (dtmp2 / dtmp3 .le. dprev) iconvg = 1
!
write(6,*)iter, dtmp2 / dtmp3, dtmp1
call flush(6)
!
call zgemv('n', dim1, dim1, cunity, a, dim1, dk, 1, czero, ctmp1, 1)  ! can be commented?
rkp(1:dim1) = rk(1:dim1) - alphak * ctmp1(1:dim1)
!
call zgemv('c', dim1, dim1, cunity, a, dim1, fk, 1, czero, ctmp1, 1)
skp(1:dim1) = sk(1:dim1) - dconjg(alphak) * ctmp1(1:dim1)
!
res1 = zdotc(dim1, skp, 1, rkp, 1)
res2 = zdotc(dim1, sk,  1, rk , 1)
!
betak = res1 / res2
!
dk(1:dim1) = rkp(1:dim1) + betak * dk(1:dim1)
fk(1:dim1) = skp(1:dim1) + dconjg(betak) * fk(1:dim1)
!
! regular update k index
!
rk(1:dim1) = rkp(1:dim1) 
sk(1:dim1) = skp(1:dim1)
!
if (iconvg .eq. 1) then
 write(6,*)
 write(6,*)' converged'
 call flush(6)
 ierror = 0
 goto 200
end if
!
iter = iter + 1
if (iter .gt. maxit0) then
 write(6,*)
 write(6,*)' error! maxit0 reached without convergency ', iter
 stop
end if
!
goto 100
!
200 continue
!
!x(1:dim1) = xk(1:dim1)
!
deallocate(ctmp1, STAT=istat)
!deallocate(sk, rk, skp, rkp, dk, fk, xk, yk, STAT=istat)
deallocate(sk, rk, skp, rkp, dk, fk, STAT=istat)
!
end subroutine bicg0
