subroutine solve_jacobi
! Jacobi iteration method; see Eq.(5) of J. Chem. Phys. 147, 044105 (2017)
! solving linear problem A * x = b
use omp_lib
use matmod_omp
use jacobimod
use hbmod
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, nl
integer                 :: istat, iter, ierr, ispin, iorbs
integer                 :: nnz
integer, parameter      :: iter_show = 50
integer*8               :: lni, lnj, lnk, lnm, lnum
integer                 :: iop0
real*8                  :: dtmp1, dtmp2
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: dprev0, derr0, derrmin
real*8                  :: dnrm2_long
logical                 :: fexist
complex*16              :: ctmp1, ctmp2
complex*16, allocatable :: cvec1(:), cmat1(:,:)
external dnrm2_long
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_jacobi: entering solve_jacobi, nrho ', nrho
!
if (lhb .or. lhf) then
    write(6,*)
    write(6,*)'solve_jacobi: lhb=T or lhf=T not available. Abort...'
    stop
end if
if (ltrun_der) then
    write(6,*)
    write(6,*)'solve_jacobi: ltrun_der=T not available. Abort... '
    stop
end if
!
allocate(cvec1(nrho**2), cmat1(nrho,nrho), STAT=istat)
!
iter = 0
if (.not. lsparse) then
   lnum = nunk * nrho**2 * 2   ! last factor 2 accounts for both real and imaginary parts
else
   lnum = (lunk_spa + 1) * 2
end if
dtmp1 = dble(lnum * 8) / dble(1024**2)
write(6,*)'solve_jacobi: length of JACOBI vector ', lnum 
call flush(6)
!
dprev0 = crit           ! JACOBI working-subroutine uses residue norm
                        ! as its convergence criterion
                        ! maximal iteration numbers: maxit0
!
! initialize
!
jacvecs    = 0.d0
!
! Here, the user needs to make sure the zeroth-tier of rho (the
! density matrix) satisfies Hermicity. 
! rho should also satisfy the normalization condition that 
! trace(rho) = 1
! rho is obtained either by
! initialization or by reading the old result from disk
!
if (.not. lsparse) then
   do ni=1,nrho
      do nj=ni+1,nrho
         rho(ni,nj,1) = dconjg(rho(nj,ni,1))
      end do
      rho(ni,ni,1) = dcmplx(dble(rho(ni,ni,1)), 0.d0)
   end do
   ctmp1 = czero
   do ni=1,nrho
      ctmp1 = ctmp1 + rho(ni,ni,1)
   end do
   if (cdabs(ctmp1 - cunity) .gt. dpico) then
       write(6,*)
       write(6,*)'solve_jacobi: error! rho is not normalized ', ctmp1
       stop
   end if
   lnj = 1
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            jacvecs(lnj,1)   = dble (rho(ni,nj,lni))
            jacvecs(lnj+1,1) = dimag(rho(ni,nj,lni))
            lnj = lnj + 2
         end do 
      end do
   end do
   call calc_ax_tfqmr_omp(nrho, jacvecs(1,1), jacvecs(1,2))
!write(6,*)'checkpoint 1'
!call cmatout(nrho, nrho, rho(1,1,1), dmtmp1)
!call flush(6)
else 
   nnz = nnz_spa(1)
   lni = ind_spa(1)
   call zmat_herm_coo('u', nnz, rho_spa(lni), irow_spa(lni), icol_spa(lni), cmtmp1, nrho, nrho)
   ctmp1 = czero
   do ni=1,nnz
      if (irow_spa(lni-1+ni) .eq. icol_spa(lni-1+ni)) then
          ctmp1 = ctmp1 + rho_spa(lni-1+ni)
      end if
   end do
   if (cdabs(ctmp1 - cunity) .gt. dpico) then
       write(6,*)
       write(6,*)'solve_jacobi: Warning! rho is not normalized ', ctmp1
       call flush(6)
   end if
   lnj = 1
   do lni=1,lunk_spa
      jacvecs(lnj,1)   = dble (rho_spa(lni))
      jacvecs(lnj+1,1) = dimag(rho_spa(lni))
      lnj = lnj + 2
   end do
   call calc_ax_tfqmr_spa_omp(lunk_spa, jacvecs(1,1), jacvecs(1,2))
end if
call refresh_rhosys
!
call prelude_jacobi(rhosys(1,1))
!
write(6,*)
write(6,*)'solve_jacobi: rhosys before jacobi cycles '
call flush(6)
call cmatout(nrho, nrho, rhosys, dmtmp1)
!
jacvecs(lnum-1,2) = jacvecs(lnum-1,2) - 1.d0
!jacvecs(1,2) = jacvecs(1,2) - 1.d0
dtmp2 = dnrm2_long(lnum, jacvecs(1,2), 1)
write(6,*)
write(6,*)'solve_jacobi: initial norm ', dtmp2
if (dtmp2 .le. dprev0) then
   write(6,*)'solve_jacobi: initial norm smaller than criterion ', dprev0
   write(6,*)'solve_jacobi: no need to go through JACOBI process '
   goto 100
end if
derrmin = dtmp2
!
! Jacobi Iteraction Algorithm
!
iop0 = 0
do 
  if (lsparse) then
      call calc_ax_jacobi_spa_omp(lunk_spa, jacvecs(1,1), jacvecs(1,2))
  else
      call calc_ax_jacobi_omp(nrho, jacvecs(1,1), jacvecs(1,2))
  end if
!
  dtmp1 = 0.d0
  do lni=1,lunk_spa*2
     dtmp1 = dtmp1 + (jacvecs(lni,1) - jacvecs(lni,2))**2
  end do
  dtmp2 = 0.d0
  lnj = 1
  if (lsparse) then
      nnz = nnz_spa(1)
      lni = ind_spa(1)
      do ni=1,nnz
         if (irow_spa(lni-1+ni) .eq. icol_spa(lni-1+ni)) then
             dtmp2 = dtmp2 + jacvecs(lnj,2)
         end if
         lnj = lnj + 2
      end do
  else
      do nj=1,nrho
         do ni=1,nrho
            if (ni .eq. nj) dtmp2 = dtmp2 + jacvecs(lnj,2)
            lnj = lnj + 2
         end do
      end do
  end if
  write(6,1000)iter, dtmp2, dsqrt(dtmp1)
  call flush(6)
!
  jacvecs(1:lnum,1) = jacvecs(1:lnum,2)
!  
  if (mod(iter, iter_show) .eq. 0) then
     if (lsparse) then
         call calc_ax_tfqmr_spa_omp(lunk_spa, jacvecs(1,1), jacvecs(1,2))
     else
         call calc_ax_tfqmr_omp(nrho, jacvecs(1,1), jacvecs(1,2))
     end if
     jacvecs(lnum-1,2) = jacvecs(lnum-1,2) - 1.d0
!     jacvecs(1,2) = jacvecs(1,2) - 1.d0
     derr0 = dnrm2_long(lnum, jacvecs(1,2), 1)
     write(6,*)'solve_jacobi: residue norm ', iter, derr0
     call flush(6)
!
     if (derr0 .lt. derrmin) then
         derrmin = derr0
         if (lsparse) then
             open(unit=69, file='rho_spa.jac', form='binary', status='unknown')
             rewind(69)
             write(69)norbs, nspin, ncor, ntier, nalf
             write(69)nunk
             lnj = 1
             do lni=1,nunk
                nnz = nnz_spa(lni)
                do ni=1,nnz
                   cvec1(ni) = dcmplx(jacvecs(lnj,1), jacvecs(lnj+1,1))
                   lnj = lnj + 2
                end do
                if (nnz .gt. 0) then
                   write(69)(cvec1(ni), ni=1,nnz)
                end if
             end do
             close(69)
         else
             open(unit=70, file='variables.jac', form='unformatted', status='unknown')
             rewind(70)
             write(70)norbs, nspin, ncor, ntier, nalf
             write(70)nunk
             lnj = 1
             do lni=1,nunk
                do nj=1,nrho
                   do ni=1,nrho
                      cmat1(ni,nj) = dcmplx(jacvecs(lnj,1), jacvecs(lnj+1,1))
                      lnj = lnj + 2
                   end do
                end do
                write(70)((cmat1(ni,nj), ni=1,nrho), nj=1,nrho)
             end do
             close(70)
         end if
     end if
!
     if (derr0 .le. crit) then
        iop0 = 1
        exit
     end if
     fexist = .false.
     inquire(file="stopst", exist=fexist)  ! exit if 
     if (fexist) then
        iop0 = 2
        write(6,*)'solve_jacobi: file stopst detected, exit from JACOBI solver'
        exit
     end if
  end if
!
  if (iter .ge. maxit0) exit
  iter = iter + 1
end do
!
if (iop0 .eq. 1) then
    write(6,*)'solve_jacobi: JACOBI converged after', iter, ' iterations'
else if (iop0 .eq. 2) then
    write(6,*)'solve_jacobi: JACOBI stopped after', iter, ' iterations'
else 
    write(6,*)'solve_jacobi: JACOBI failed to converge after', iter, 'iterations'
end if
call flush(6)
!
! recover rho
!
if (lsparse) then
   lnj = 1
   do lni=1,lunk_spa
      rho_spa(lni) = dcmplx(jacvecs(lnj,1), jacvecs(lnj+1,1))
      lnj = lnj + 2
   end do
   call calc_ax_tfqmr_spa_omp(lunk_spa, rho_spa, jacvecs(1,2))
else
   lnj = 1
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            rho(ni,nj,lni) = dcmplx(jacvecs(lnj,1), jacvecs(lnj+1,1))
            lnj = lnj + 2
         end do
      end do
   end do
   call calc_ax_tfqmr_omp(nrho, rho(1,1,1), jacvecs(1,2))
end if
jacvecs(lnum-1,2) = jacvecs(lnum-1,2) - 1.d0
!jacvecs(1,2) = jacvecs(1,2) - 1.d0
dtmp2 = dnrm2_long(lnum, jacvecs(1,2), 1)
write(6,*)'solve_jacobi: final norm ', dtmp2
call flush(6)
!
100 continue
! 
! important: recover the upper half of rho by Hermicity
!  
! Post-JACOBI analysis begins here
! 
deallocate(cvec1, cmat1, STAT=istat)
!
call outsteady(iter)
call cpu_time(cpu4)
write(6,*)
write(6,*)'solve_jacobi: leaving solve_jacobi ', cpu4 - cpu3
call flush(6)
!
1000 format(2x, I6, 1x, f10.6, 1x, e20.8e3)
end subroutine solve_jacobi
!
