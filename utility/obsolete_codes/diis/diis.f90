subroutine diis(iter, ierror)
use diismod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
!     Direct Inversion in the Iterative Subspace (DIIS)
!
!     Solve the linear homogeneous equation d * X = g
!
!     X = d0^(-1) * (g - d1 * X)
! 
!     If d0 is identity matrix, Xnew = g - (d - I) * Xold
!                                    = g + Xold - d * Xold
!
integer, intent(out) :: iter, ierror
!
! ierror =  0  : diis successful, normal output.
!          -1  : maxit0 reached, not converged yet
!
logical :: iexist
integer :: ni, nj, nk, i, j, k
integer :: niter, istat, lwork, info, iover, ineglect
integer :: nitmax
integer*8 :: dimmax, dimvec, lni
real*8    :: dprev 
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, cdotlong
!
integer, allocatable :: ipiv(:), ipoint(:)
real*8,  allocatable :: b(:), y(:), lambda(:), coef(:,:), work(:)
real*8,  allocatable :: coef0(:,:) 
complex*16 :: ztmp1, zdotc
external   :: zdotc, cdotlong
real*8 :: cpu1, cpu2, cpu3, cpu4, cpu5, cpu6
!
call cpu_time(cpu1)
ierror = 1
dimvec = nunk * nrho**2
dimmax = dimvec
nitmax = 300
iter   = 1
iover  = 0
dprev  = crit
!
allocate(b(nitmax), y(nitmax), lambda(nitmax), coef(nitmax,nitmax), STAT=istat)
allocate(ipoint(nitmax), coef0(nitmax,nitmax), STAT=istat)
!
write(6,*)
write(6,*)'nitmax         ', nitmax
write(6,*)'dimvec, dimmax ', dimvec, dimmax
call flush(6)
!
do ni=1,nitmax
  ipoint(ni) = ni          ! from 1 to min(nitmax,niter)-1 : record # from earliest to latest
                           ! ipoint( min(nitmax,niter) )   : record # to write
end do
!
niter = 1
inquire(file='diis.xtmp', exist=iexist, err=998)
998 continue
if (iexist) then
  open(unit=20, file='diis.xtmp', status='unknown', form='unformatted')
  close(unit=20, status="delete")
end if
!
inquire(file='diis.err', exist=iexist, err=999)
999 continue
if (iexist) then
  open(unit=21, file='diis.err', status='unknown', form='unformatted')
  close(unit=21, status="delete")
end if
!
open(unit=20, file='diis.xtmp', status='unknown', access='direct', form='unformatted', recl=16*dimmax)
open(unit=21, file='diis.err',  status='unknown', access='direct', form='unformatted', recl=16*dimmax)
!
b(1)        = 1.d0
b(2:nitmax) = 0.d0
!
rhs(1:nrho, 1:nrho, 1:nunk) = czero
write(20,rec=ipoint(niter)) (((rhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
xold(1:nrho, 1:nrho, 1:nunk) = rho(1:nrho, 1:nrho, 1:nunk)
do ni=1,nrho
  do nj=ni+1,nrho
    xold(ni,nj,1) = czero
  end do
  xold(ni,ni,1) = dcmplx(dble(xold(ni,ni,1)), 0.d0)
end do
!
9  write(6,*)
write(6,*)'   NITER         error '
call flush(6)
!
call cpu_time(cpu3)
10   continue
!
call calcxnew
!
niter = niter + 1
!
dtmp1 = 0.d0
do lni=1,nunk
  do nj=1,nrho
    do ni=1,nrho
      rho(ni,nj,lni) = xnew(ni,nj,lni) - xold(ni,nj,lni)
      dtmp1          = dtmp1 + cdabs(rho(ni,nj,lni))
    end do
  end do
end do
write(21,rec=ipoint(niter-1)) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
call cpu_time(cpu4)
write(6,519)iter, dtmp1, ineglect, lambda(ineglect), cpu4-cpu3
call cpu_time(cpu3)
call flush(6)
518 format(I6, 2x, 10(e12.4e2, 2x))
519 format(I6, 2x, e12.4e2, 2x, I6, 2x, 10(e12.4e2, 2x))
!
if (dtmp1 .le. dprev) then
 rho(1:nrho, 1:nrho, 1:nunk) = xold(1:nrho, 1:nrho, 1:nunk)
 write(6,*)
 write(6,*)' DIIS converged after ', niter, ' cycles.'
 call flush(6)
 goto 100
endif
!
if (iter .gt. maxit0) then
 rho(1:nrho, 1:nrho, 1:nunk) = xold(1:nrho, 1:nrho, 1:nunk)
 write(6,*)
 write(6,*)' DIIS fail to converge after ', iter, ' cycles.'
 call flush(6)
 ierror = -1 
 deallocate(ipoint, b, y, lambda, coef, coef0, STAT=istat)
 close(unit=20, status="delete")
 close(unit=21, status="delete")
 return
end if
!
if (niter .gt. nitmax) then
  iover = 1
  niter = nitmax
!
! rearrange ipoint (remove the earliest)
!
  nk = ipoint(1)
  do ni=1,nitmax-1
    ipoint(ni) = ipoint(ni+1)
  end do
  ipoint(nitmax) = nk
!
! rearrange ipoint (remove ineglect-th record)
!
!  nk = ipoint(ineglect) 
!  do ni=ineglect,nitmax-1
!    ipoint(ni) = ipoint(ni+1)
!  end do
!  ipoint(nitmax) = nk
!
endif
!
write(20,rec=ipoint(niter)) (((xnew(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
coef(1,1)        = 0.d0
coef(1, 2:niter) = 1.d0
!
if (niter .eq. 2) then
  read(21,rec=ipoint(1)) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
  coef(2,2) = cdotlong(dimvec, rho, rho)
else
  if (iover .eq. 0) then
    do ni=2,niter-1
      do nj=ni,niter-1
        coef(ni,nj) = coef0(ni,nj)
      end do
    end do
  else if (iover .eq. 1) then  ! remove the second row and column from coef0
    coef(1:niter, 1:niter) = coef0(1:niter, 1:ntier)
    do ni=3,niter
      coef(ni-1,1:niter) = coef(ni,1:niter)
    end do
    do ni=3,niter
      coef(1:niter,ni-1) = coef(1:niter,ni) 
    end do
!
!  else if (iover .eq. 1) then ! remove the ineglect-th row and column from coef0
!    coef(1:niter, 1:niter) = coef0(1:niter, 1:ntier)
!    do ni=ineglect+1,niter
!      coef(ni-1,1:niter) = coef(ni,1:niter)
!    end do
!    do ni=ineglect+1,niter
!      coef(1:niter,ni-1) = coef(1:niter,ni) 
!    end do
!
  else
    write(6,*)
    write(6,*)' error! unknown iover in diis ', iover 
    stop
  end if
  read(21,rec=ipoint(niter-1)) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
  coef(niter, niter) = cdotlong(dimvec, rho, rho)
  do ni=2,niter-1
    read(21,rec=ipoint(ni-1)) (((rhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
    coef(ni,niter) = cdotlong(dimvec, rhs, rho)
  end do
end if
!
! save coef to coef0
!
coef0(1:nitmax, 1:nitmax) = 0.d0
coef0(1:niter, 1:niter)   = coef(1:niter, 1:niter)
!
! Solve the linear equation coef * y = b (coef is real symmetric, upper triangle)
!
lwork = niter**2
allocate(work(lwork), ipiv(niter), STAT=istat)
y(1:niter) = b(1:niter)
call dsysv('u', niter, 1, coef, nitmax, ipiv, y, nitmax, work, lwork, info)
if (info .ne. 0) then
 write(6,*)
 write(6,*)' error solving linear equation when processing diis ', niter
 write(6,*)info, nitmax
 stop
endif
deallocate(work, ipiv, STAT=istat)
!
dtmp2 = 0.d0
dtmp3 = 1.d20
do ni=2,niter 
 lambda(ni) = y(ni)
 dtmp2      = dtmp2 + lambda(ni)
 if (dtmp3 .ge. dabs(lambda(ni))) then
   dtmp3    = dabs(lambda(ni))
   ineglect = ni 
 end if
enddo
!
xold(1:nrho, 1:nrho, 1:nunk) = czero
do nj=2,niter
  if (dabs(lambda(nj)) .ge. dnano) then
    read(20,rec=ipoint(nj)) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
    xold(1:nrho, 1:nrho, 1:nunk) = xold(1:nrho, 1:nrho, 1:nunk) + lambda(nj) * &
                                    rho(1:nrho, 1:nrho, 1:nunk) 
  end if
end do
iter = iter + 1
goto 10
!
100  continue
!
ierror = 0
deallocate(ipoint, b, y, lambda, coef, coef0, STAT=istat)
close(unit=20, status="delete")
close(unit=21, status="delete")
!
call cpu_time(cpu2)
write(6,*)
write(6,*)' cpu time for diis ', cpu2 - cpu1
call flush(6)
!
end subroutine diis
!--------------------------------
subroutine calcxnew
use diismod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
xnew(1:nrho, 1:nrho, 1:nunk) = xold(1:nrho, 1:nrho, 1:nunk)
xnew(1,1,1) = xnew(1,1,1) + cunity
!
call calc_ax_bicg(nrho, xold, rhs)
!
xnew(1:nrho, 1:nrho, 1:nunk) = xnew(1:nrho, 1:nrho, 1:nunk) -  &
                                rhs(1:nrho, 1:nrho, 1:nunk)
end subroutine calcxnew
