subroutine diis0(dimvec,ldmat,d,g,X,ierror)
implicit none
!
!     Direct Inversion in the Iterative Subspace (DIIS)
!
!     Solve the linear homogeneous equation d * X = g
!
!     X = d0^(-1) * (g - d1 * X)
! 
!     if d0 is identity matrix, Xnew = g - (d - I) * Xold
!                                    = g + Xold - d * Xold
!
! note : upon exit, the coefficient matrix D is destroyed
!
integer, intent(in)       :: dimvec, ldmat
integer, intent(out)      :: ierror
!
! ierror =  0  : diis successful, normal output.
!          -1  : maxit0 reached, not converged yet
!
complex*16, intent(in)    :: g(*)
complex*16, intent(inout) :: d(ldmat, *)
complex*16, intent(out)   :: X(ldmat)
!
integer :: ni,nj,nk
integer :: iter, niter, istat, lwork, info
integer :: nitmax, dimmax
integer, parameter :: maxit0 = 1000
real*8            :: dtmp1, dtmp2
real*8, parameter :: dprev = 1.d-7
!
complex*16, parameter :: czero = (0.d0, 0.d0)
complex*16, parameter ::   eye = (0.d0, 1.d0)
complex*16, parameter :: cunity= (1.d0, 0.d0)
!
integer, allocatable :: ipiv(:)
real*8, allocatable ::  b(:), y(:), lambda(:), coef(:,:), work(:)
complex*16, allocatable :: err(:,:), xtmp(:,:), xold(:), xnew(:)
real*8 :: cpu1, cpu2
!
call cpu_time(cpu1)
ierror = 1
dimmax = dimvec
nitmax = 100
iter   = 1
!
allocate(b(nitmax), y(nitmax), lambda(nitmax), coef(nitmax,nitmax), STAT=istat)
allocate(err(dimmax,nitmax), xtmp(dimmax,nitmax), xold(dimmax), xnew(dimmax), STAT=istat)
!
write(6,*)
write(6,*)'dimvec, ldmat, dimmax ', dimvec, ldmat, dimmax
write(6,*)'nitmax                ', nitmax
call flush(6)
!
niter = 1
!
b(1)        = 1.d0
b(2:nitmax) = 0.d0
!
xtmp(1:dimvec, niter) = czero
xold(1:dimvec)        = czero
!
9  write(6,*)
write(6,*)'   NITER         error '
call flush(6)
!
10   continue
call formnewx(niter, dimvec, ldmat, dimmax, d, g, xold, xnew)
!
niter = niter + 1
!
dtmp1 = 0.d0
do ni=1,dimvec
 dtmp1           = dtmp1 + cdabs(xnew(ni) - xold(ni))
 err(ni,niter-1) = xnew(ni) - xold(ni)
enddo
!
write(6,*)iter, dtmp1
call flush(6)
!
if (dtmp1 .le. dprev) then
 X(1:dimvec) = xold(1:dimvec)
 write(6,*)
 write(6,*)' DIIS converged after ', niter, ' cycles.'
 call flush(6)
 goto 100
endif
!
if (iter .gt. maxit0) then
 write(6,*)
 write(6,*)' DIIS does not converge after ', iter, ' loops'
 call flush(6)
 ierror = -1 
 return
end if
!
if (niter .gt. nitmax) then
 niter = nitmax
 do nj=1,dimvec
  do ni=1,nitmax-1
   xtmp(nj,ni) = xtmp(nj,ni + 1)
    err(nj,ni) =  err(nj,ni + 1)
  enddo
 enddo
endif
!
xtmp(1:dimvec, niter) = xnew(1:dimvec)
!
coef(1,1)        = 0.d0
coef(1, 2:niter) = 1.d0
!
do ni=2,niter
 do nj=ni,niter
  coef(ni,nj) = 0.d0
  do nk=1,dimvec
   coef(ni,nj) = coef(ni,nj) + dble( dconjg(err(nk,nj-1)) *  &
                 err(nk,ni-1) )
  enddo
 enddo
enddo
!
! Solve the linear equation coef * y = b (coef is real symmetric, upper triangle)
!
lwork = niter**2
allocate(work(lwork), ipiv(niter), STAT=istat)
!
y(1:niter) = b(1:niter)
call dsysv('u', niter, 1, coef, nitmax, ipiv, y, nitmax, work, lwork, info)
!
if (info .ne. 0) then
 write(6,*)
 write(6,*)' error solving linear equation when processing diis ', niter
 write(6,*)info, nitmax
 stop
endif
!
deallocate(work, ipiv, STAT=istat)
!
dtmp2 = 0.d0
do ni=2,niter 
 lambda(ni) = y(ni)
 dtmp2      = dtmp2 + lambda(ni)
enddo
!
do ni=1,dimvec
 xold(ni) = czero
 do nj=2,niter
  xold(ni) = xold(ni) + lambda(nj) * xtmp(ni,nj)
 enddo
enddo
!
iter = iter + 1
goto 10
!
100  continue
!
ierror = 0
!
deallocate(b, y, lambda, coef, STAT=istat)
deallocate(err, xtmp, xold, xnew, STAT=istat)
!
call cpu_time(cpu2)
write(6,*)
write(6,*)' cpu time for diis ', cpu2 - cpu1
call flush(6)
!
end subroutine diis0
!
!--------------------------------
subroutine formnewx(niter,dimvec,ldmat,dimmax,d,g,xold,xnew) 
implicit none
!
integer, intent(in) :: niter, dimvec, ldmat, dimmax
complex*16, intent(in) :: d(ldmat,*), g(ldmat), xold(dimmax)
complex*16, intent(out) :: xnew(dimmax)
!
integer :: ni,nj
!
do ni=1,dimvec
  xnew(ni) = g(ni) + xold(ni)
  do nj=1,dimvec
     xnew(ni) = xnew(ni) - d(ni,nj) * xold(nj) 
  enddo
enddo
!
end subroutine formnewx
