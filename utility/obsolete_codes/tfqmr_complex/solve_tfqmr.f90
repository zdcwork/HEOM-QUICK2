subroutine solve_tfqmr
use tfqmrmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer   :: ni, nj, nk, nl
integer   :: istat, iter, ierr, maxiter
integer   :: info0(4) 
integer   :: revcom, colx, colb
integer*8 :: lni, lnj, lnk, lnum
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8    :: cpu1, cpu2, cpu3, cpu4
real*8    :: dprev0, derr0
complex*16, allocatable :: zvectmp0(:)
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_tfqmr: entering solve_tfqmr, nrho ', nrho
call flush(6)
!
iter = 0
lnum = nunk * nrho**2
write(6,*)'solve_tfqmr: length of TFQMR vector ', lnum 
!
allocate(zvectmp0(lnum), STAT=istat)
!
dprev0 = crit**2        ! TFQMR working-subroutine uses residue norm
                        ! instead of square of norm as its convergence criterion
maxiter = maxit0 / 2    ! TFQMR has two iteration numbers
!
! initialize
!
info0(1:4) = 0
info0(1)   = 600        ! write iterative procedure to file channel #6
!info0(1)   = 10600       
qmrvecs    = czero
xqmr       = czero
bqmr       = czero
!
!xqmr(1:nrho,1:nrho,1:nunk) = rho(1:nrho,1:nrho,1:nunk)
!do ni=1,nrho
!   do nj=ni+1,nrho
!      xqmr(ni,nj,1) = czero
!   end do
!end do
!do ni=1,nrho
!   xqmr(ni,ni,1) = dcmplx(dble(xqmr(ni,ni,1)), 0.d0)
!end do
!call calc_ax_tfqmr_comp(nrho, xqmr, bqmr)
!bqmr = -bqmr
bqmr(1,1,1) = bqmr(1,1,1) + cunity   ! first element of B vector is 1, 
                                     ! due to trace(rho) = 1
call zmatvec(bqmr, nrho, nunk, qmrvecs(1,2), 0)
!
! Transpose-Free Quasi-Minimal Residual Algorithm
!
do 
  call zutfx(lnum, lnum, maxiter, qmrvecs, dprev0, info0)
  revcom = info0(2)
  colx   = info0(3)
  colb   = info0(4)
!  
  if (revcom .ne. 1) exit
  iter = iter + 1
  call calc_ax_tfqmr(nrho, nunk, qmrvecs(1,colx), qmrvecs(1,colb))

  call calc_ax_tfqmr(nrho, nunk, qmrvecs(1,1), zvectmp0(1))
  derr0 = 0.d0
  do lni=1,lnum
     derr0 = derr0 + (cdabs(zvectmp0(lni) - qmrvecs(lni,2)))**2
  end do
  write(6,*)iter, derr0
end do
!
if (iter > maxiter) then
   write(6,*)'solve_tfqmr: Warning! '
   write(6,*)'solve_tfqmr: TFQMR procedures NOT converge after', iter, ' iterations'
end if
ierr = info0(1)
write(6,*)
write(6,*)'solve_tfqmr: ierr = ', ierr
if (ierr .eq. 0) then
   write(6,*)'solve_tfqmr: TFQMR converged succesfully'
else if (ierr .eq. 1) then
   write(6,*)'solve_tfqmr: error! invalid reverse communication call'
else if (ierr .eq. 2) then
   write(6,*)'solve_tfqmr: error! invalid input for TFQMR'
else if (ierr .eq. 4) then
   write(6,*)'solve_tfqmr: TFQMR failed to converge after ', iter, ' steps'
else if (ierr .eq. 8) then
   write(6,*)'solve_tfqmr: TFQMR ran into an exact breakdown'
else 
   write(6,*)'solve_tfqmr: unknown ierr flag, ierr = ', ierr
end if
call flush(6)
!
! recover rho
!
call zmatvec(rho(1,1,1), nrho, nunk, qmrvecs(1,1), 1)
! 
! important: recover the upper half of rho by Hermicity
!  
!do ni=1,nrho
!  do nj=ni+1,nrho
!     rho(ni,nj,1) = dconjg(rho(nj,ni,1))
!  end do
!end do
deallocate(zvectmp0, STAT=istat)
!
! Post-TFQMR analysis begins here
! 
call outsteady(iter)
!
call cpu_time(cpu4)
write(6,*)
write(6,*)'solve_tfqmr: leaving solve_tfqmr ', cpu4 - cpu3
call flush(6)
end subroutine solve_tfqmr
!
