subroutine solve_cplqmr
! Coupled two-term Look-ahead Quasi-Minimal Residue method for  
! solving linear problem A * x = b
use cplqmrmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, nl
integer                 :: istat, iter, ierr, maxiter, iter0
integer, parameter      :: iter_show = 50
integer                 :: info0(4) 
integer                 :: revcom, colx, colb
integer*8               :: lni, lnj, lnk, lnum
integer                 :: iop0
real*8                  :: dtmp1, dtmp2
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: dprev0, derr0
real*8                  :: dnrm2_long
external dnrm2_long
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_cplqmr: entering solve_cplqmr, nrho ', nrho
call flush(6)
!
iter = 0
lnum = nunk * nrho**2 * 2
write(6,*)'solve_cplqmr: length of CPLQMR vector ', lnum 
dtmp1 = dble(nunk * nrho**2 * 16) / dble(1024**2)
!
dprev0 = crit           ! CPLQMR working-subroutine uses residue norm
                        ! as its convergence criterion
dprev0 = dprev0 * 1.d-3 ! use more rigorous criterion inside 'dutfx'
!
iter0   = nlim          ! CPLQMR has two iteration numbers
maxiter = iter0
!
! initialize
!
info0(1:4) = 0
info0(1)   = 600        ! write iterative procedure to file channel #6
!info0(1)   = 606  
!info0(1)   = 10600       
qmrvecs    = 0.d0
!
! Here, the user needs to make sure the zeroth-tier of rho (the
! density matrix) satisfies Hermicity. rho is obtained either by
! initialization or by reading the old result from disk.
!
do ni=1,nrho
   do nj=ni+1,nrho
      rho(ni,nj,1) = dconjg(rho(nj,ni,1))
   end do
   rho(ni,ni,1) = dcmplx(dble(rho(ni,ni,1)), 0.d0)
end do
call calc_ax_cplqmr(nrho, rho(1,1,1), qmrvecs(1,2))
qmrvecs(1:lnum,2) = -qmrvecs(1:lnum,2)
qmrvecs(1,2) = 1.d0 + qmrvecs(1,2)    ! first element of B vector is 1, 
!                                     ! due to trace(rho) = 1
dtmp2 = dnrm2_long(lnum, qmrvecs(1,2), 1)
write(6,*)
write(6,*)'solve_cplqmr: initial norm ', dtmp2
!
! Coupled two-term Look-ahead Quasi-Minimal Residue Algorithm
!
call pre_atrx_bicg
!
iop0 = 0
do 
  call ducpl(lnum, lnum, maxiter, maxpq, maxvw, dimm, mvec, dnorms,  &
             dwk, idx, iwk, qmrvecs, dprev0, info0)
  revcom = info0(2)
  colx   = info0(3)
  colb   = info0(4)
!  
  if (mod(iter, iter_show) .eq. 0) then
     call calc_ax_cplqmr(nrho, qmrvecs(1,1), qmrvecs(1,3))
     qmrvecs(1:lnum,3) = qmrvecs(1:lnum,3) - qmrvecs(1:lnum,2)
     derr0 = dnrm2_long(lnum, qmrvecs(1,3), 1)
     write(6,*)'iter_cplqmr ', iter, derr0
     if (derr0 .le. crit) then
        iop0 = 1
        exit
     end if
  end if
!
  if (revcom .ne. 1 .and. revcom .ne. 2) exit
  iter = iter + 1
!dtmp1 = 0.d0
!dtmp2 = 0.d0
!do ni=1,nrho
!   do nj=ni+1,nrho
!      nk = ((nj - 1) * nrho + ni - 1) * 2 + 1
!      dtmp1 = dtmp1 + qmrvecs(nk,colx)**2 + qmrvecs(nk+1,colx)**2
!      dtmp2 = dtmp2 + qmrvecs(nk,colb)**2 + qmrvecs(nk+1,colb)**2
!   end do
!end do
!write(6,*)'solve_cplqmr: revcom x b', revcom, dtmp1, dtmp2
!call flush(6)
  if (revcom .eq. 1) then
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_cplqmr: before ax ', dtmp1, dtmp2
     call calc_ax_cplqmr(nrho, qmrvecs(1,colx), qmrvecs(1,colb))
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_cplqmr: after  ax ', dtmp1, dtmp2
  else if (revcom .eq. 2) then 
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_cplqmr: before atrx ', dtmp1, dtmp2
     call calc_atrx_cplqmr(nrho, qmrvecs(1,colx), qmrvecs(1,colb))
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_cplqmr: after  atrx ', dtmp1, dtmp2
  else
     write(6,*)'solve_cplqmr: error! unknown revcom ', revcom
     stop
  end if
  call flush(6)
end do
!
if (maxiter > iter0) then
   write(6,*)'solve_cplqmr: error! maxiter = ', maxiter
   write(6,*)'solve_cplqmr: CPLQMR procedures NOT converge after', iter, ' iterations'
else 
   write(6,*)'solve_cplqmr: CPLQMR converged after', iter, ' iterations'
end if
ierr = info0(1)
if (iop0 .eq. 0) then
   write(6,*)
   write(6,*)'solve_cplqmr: ierr = ', ierr
   if (ierr .eq. 0) then
      write(6,*)'solve_cplqmr: CPLQMR converged succesfully'
   else if (ierr .eq. 1) then
      write(6,*)'solve_cplqmr: error! invalid reverse communication call'
   else if (ierr .eq. 2) then
      write(6,*)'solve_cplqmr: error! invalid input for CPLQMR'
   else if (ierr .eq. 4) then
      write(6,*)'solve_cplqmr: CPLQMR failed to converge after ', iter, ' steps'
   else if (ierr .eq. 8) then
      write(6,*)'solve_cplqmr: last block of CPLQMR could not be closed'
   else if (ierr .eq. 16) then
      write(6,*)'solve_cplqmr: an A-invariant subspace has been found'
   else if (ierr .eq. 32) then
      write(6,*)'solve_cplqmr: an A^T-invariant subspace has been found'
   else if (ierr .eq. 48) then
      write(6,*)'solve_cplqmr: both A- and A^T-invariant subspaces have been found'
   else if (ierr .eq. 64) then
      write(6,*)'solve_cplqmr: insufficient memory for vecs (mvec < maxpq+maxvw) '
   else if (ierr .eq. 128) then
      write(6,*)'solve_cplqmr: could not rebuild a regular vector '
   else 
      write(6,*)'solve_cplqmr: unknown ierr flag, ierr = ', ierr
   end if
end if
call flush(6)
!
! recover rho
!
lnj = 1
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         rho(ni,nj,lni) = rho(ni,nj,lni) + dcmplx(qmrvecs(lnj,1), qmrvecs(lnj+1,1))
         lnj = lnj + 2
      end do
   end do
end do
!do ni=1,nrho
!   do nj=ni+1,nrho
!      rho(ni,nj,1) = dconjg(rho(nj,ni,1))
!   end do
!end do
!
call calc_ax_cplqmr(nrho, rho(1,1,1), qmrvecs(1,3))
qmrvecs(1,3) = qmrvecs(1,3) - 1.d0
dtmp2 = dnrm2_long(lnum, qmrvecs(1,3), 1)
write(6,*)'solve_cplqmr: final norm ', dtmp2
! 
! important: recover the upper half of rho by Hermicity
!  
! Post-CPLQMR analysis begins here
! 
call outsteady(iter)
!
call cpu_time(cpu4)
write(6,*)
write(6,*)'solve_cplqmr: leaving solve_cplqmr ', cpu4 - cpu3
call flush(6)
end subroutine solve_cplqmr
!
