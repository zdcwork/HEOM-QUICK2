subroutine solve_cf_tfqmr(rhoinp, isgnb, ispinb, iorbsb, isgnw, zouth, zouta)
! Transpose-Free Quasi-Minimal Residue method 
! for system correlation functions
! and system spectral functions by 
! solving linear problem A * x = b
use omp_lib
use matmod_omp
use corrfuncmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent (in) :: rhoinp(nrho,nrho,*)
complex*16, intent (out):: zouth(nrho,nrho), zouta(nrho,nrho)
integer, intent(in)     :: isgnb, ispinb, iorbsb, isgnw
! isgnw = 1: +i * omega in the linear problem
! isgnw = 2: -i * omega in the linear problem
!
integer                 :: ni, nj, nk, nl
integer                 :: istat, iter, ierr, maxiter, iter0, ispin, iorbs, ialf
integer, parameter      :: iter_show = 50
integer                 :: info0(4) 
integer                 :: revcom, colx, colb
integer*8               :: lni, lnj, lnk, lnum, lrho
integer                 :: iop0
real*8                  :: dtmp1, dtmp2
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: dprev0, derr0
real*8                  :: dlamch, cdotlong, dnrm2_long
real*8                  :: occspin(2)
external dlamch, cdotlong, dnrm2_long
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_cf_tfqmr: entering solve_cf_tfqmr, nrho ', nrho
call flush(6)
!
iter = 0
lrho = nunk * nrho**2 * 2   ! factor of 2 for two types of symmetries: brhoh and brhoa
lnum = lrho * 2             ! factor of 2 for real and imaginary parts 
write(6,*)'solve_cf_tfqmr: length of TFQMR vector ', lnum 
dtmp1 = dble(nunk * 2 * nrho**2 * 16) / dble(1024**2)
!
dprev0 = crit_dos       ! TFQMR working-subroutine uses residue norm
                        ! as its convergence criterion
dprev0 = dprev0 * 1.d-3 ! use more rigorous criterion inside 'dutfx'
!
iter0  = maxit_dos / 2  ! TFQMR has two iteration numbers
maxiter = iter0
!
! lplus_cf enters as a flag in <calc_ax_cf>
if (isgnw .eq. 1) then
   lplus_cf = .true.
else if (isgnw .eq. 2) then
   lplus_cf = .false.
else 
   write(6,*)'solve_cf_tfqmr: error! unknown isgnw = ', isgnw
   stop
end if
lplus = lplus_cf
!
! initialize
!
call calchs(rhoinp(1,1,1))
if (igroundsteady .eq. 0) then          ! ground state
   dhs = czero 
else if (igroundsteady .eq. 1) then     ! steady state
   call calcdhs(1.d10, rhoinp(1,1,1))
else if (igroundsteady .eq. 2) then     ! transient time
!   call calcdhs(tnow, rhoinp(1,1,1))
!   write(6,*)'solve_cf_tfqmr: lead energy shift at time, ialf, ispin, eshift'
!   do ialf=1,nalf
!      do ispin=1,nspin
!         call getengyshift(tnow, ialf, ispin, eleadtime(ialf,ispin))
!         write(6,*)'solve_cf_tfqmr: ', tnow, ialf, ispin, eleadtime(ialf,ispin)*hbar
!      end do
!   end do
!   call flush(6)
   continue ! -- to check, in principle the current code only deals with static cases
else 
   write(6,*)'solve_cf_tfqmr: error! unknown igroundsteady = ', igroundsteady
   stop
end if
!
if (ltrun_der) call prelude_cf_trun
!
info0(1:4) = 0
info0(1)   = 600        ! write iterative procedure to file channel #6
!info0(1)   = 10600       
qmrvecs_cf = 0.d0
!
! Here, the user needs to make sure the zeroth-tier of rho (the
! density matrix) satisfies Hermicity. rho is obtained either by
! initialization or by reading the old result from disk
!
if (lhf .and. .not. lfixhs) then  
   !
   ! for effective single-electron Hamiltonian, TFQMR should start with rho=zero
   ! but in fact the iteration would be very difficult to converge 
   ! (tested on HF for single-level Hubbard model)
   !
   write(6,*)'solve_cf_tfqmr: HF with TFQMR does not work. Abort...'
   stop
end if
!
brhoh(1:nrho,1:nrho,1:nunk) = czero
brhoa(1:nrho,1:nrho,1:nunk) = czero
!
do lni=1,nunk
   ! Bh^- = B * rho^- + rho^- * B^+ = (Bh^+)^dag
   ! Ba^- = (-i) * (B * rho^- - rho^- * B^+) = (Ba^+)^dag 
   if (isgnb .eq. 1) then
! Bh
      call zamsmm1('l', 'c', 'n', iorbsb, ispinb, cunity,     &
                   rhoinp(1,1,lni), nrho,  czero, brhoh(1,1,lni), nrho)
      call zamsmm1('r', 'n', 'n', iorbsb, ispinb, cunity,     &
                   rhoinp(1,1,lni), nrho, cunity, brhoh(1,1,lni), nrho)
! Ba     
      call zamsmm1('l', 'c', 'n', iorbsb, ispinb, -eye,       &
                   rhoinp(1,1,lni), nrho,  czero, brhoa(1,1,lni), nrho)
      call zamsmm1('r', 'n', 'n', iorbsb, ispinb,  eye,       &
                   rhoinp(1,1,lni), nrho, cunity, brhoa(1,1,lni), nrho)
   else if (isgnb .eq. 2) then
! Bh
      call zamsmm1('l', 'n', 'n', iorbsb, ispinb, cunity,     &
                   rhoinp(1,1,lni), nrho,  czero, brhoh(1,1,lni), nrho)
      call zamsmm1('r', 'c', 'n', iorbsb, ispinb, cunity,     &
                   rhoinp(1,1,lni), nrho, cunity, brhoh(1,1,lni), nrho)
! Ba
      call zamsmm1('l', 'n', 'n', iorbsb, ispinb, -eye,       &
                   rhoinp(1,1,lni), nrho,  czero, brhoa(1,1,lni), nrho)
      call zamsmm1('r', 'c', 'n', iorbsb, ispinb,  eye,       &
                   rhoinp(1,1,lni), nrho, cunity, brhoa(1,1,lni), nrho)
   else
      write(6,*)'solve_cf_tfqmr: error! unknown isgnb ', isgnb
      stop
   end if
end do
!
! A * x = b, where b is (-brhoh, -brhoa)
! 
lnj = 1
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         qmrvecs_cf(lnj,   2) = -dble (brhoh(ni,nj,lni))
         qmrvecs_cf(lnj+1, 2) = -dimag(brhoh(ni,nj,lni))
         lnj = lnj + 2
      end do
   end do
end do
lnj = lrho + 1
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         qmrvecs_cf(lnj,   2) = -dble (brhoa(ni,nj,lni))
         qmrvecs_cf(lnj+1, 2) = -dimag(brhoa(ni,nj,lni))
         lnj = lnj + 2
      end do
   end do
end do
!
dtmp2 = dnrm2_long(lnum, qmrvecs_cf(1,2), 1)
write(6,*)
write(6,*)'solve_cf_tfqmr: initial norm ', dtmp2
if (dtmp2 .le. crit_dos) then
   write(6,*)'solve_cf_tfqmr: initial norm smaller than criterion ', crit_dos
   write(6,*)'solve_cf_tfqmr: no need to go through TFQMR process '
   goto 100
end if
!
! Transpose-Free Quasi-Minimal Residue Algorithm
!
iop0 = 0
!
! reset dutfx 
!
info0(1:4) = -9999
call dutfx(lnum, lnum, maxiter, qmrvecs_cf, dprev0, info0)
info0(1:4) = 0
info0(1)   = 600        ! write iterative procedure to file channel #6
!
do 
  call dutfx(lnum, lnum, maxiter, qmrvecs_cf, dprev0, info0)
  revcom = info0(2)
  colx   = info0(3)
  colb   = info0(4)
!  
  if (mod(iter, iter_show) .eq. 0) then
     call calc_ax_cf_omp(nrho, qmrvecs_cf(1,1), qmrvecs_cf(lrho+1,1), qmrvecs_cf(1,9), qmrvecs_cf(lrho+1,9))
     qmrvecs_cf(1:lnum,9) = qmrvecs_cf(1:lnum,9) - qmrvecs_cf(1:lnum,2)
     derr0 = dnrm2_long(lnum, qmrvecs_cf(1,9), 1)
     write(6,*)'iter_tfqmr ', iter, derr0
     if (derr0 .le. crit_dos) then
        iop0 = 1
        exit
     end if
  end if
!
  if (revcom .ne. 1) exit
  iter = iter + 1
  call calc_ax_cf_omp(nrho, qmrvecs_cf(1,colx), qmrvecs_cf(lrho+1,colx), qmrvecs_cf(1,colb), qmrvecs_cf(lrho+1,colb))
  call flush(6)
end do
!
if (maxiter > iter0) then
   write(6,*)'solve_cf_tfqmr: error! maxiter = ', maxiter
   write(6,*)'solve_cf_tfqmr: TFQMR procedures NOT converge after', iter, ' iterations'
else 
   write(6,*)'solve_cf_tfqmr: TFQMR converged after', iter, ' iterations'
end if
ierr = info0(1)
if (iop0 .eq. 0) then
   write(6,*)
   write(6,*)'solve_cf_tfqmr: ierr = ', ierr
   if (ierr .eq. 0) then
      write(6,*)'solve_cf_tfqmr: TFQMR converged succesfully'
   else if (ierr .eq. 1) then
      write(6,*)'solve_cf_tfqmr: error! invalid reverse communication call'
   else if (ierr .eq. 2) then
      write(6,*)'solve_cf_tfqmr: error! invalid input for TFQMR'
   else if (ierr .eq. 4) then
      write(6,*)'solve_cf_tfqmr: TFQMR failed to converge after ', iter, ' steps'
   else if (ierr .eq. 8) then
      write(6,*)'solve_cf_tfqmr: TFQMR ran into an exact breakdown'
   else 
      write(6,*)'solve_cf_tfqmr: unknown ierr flag, ierr = ', ierr
   end if
end if
call flush(6)
!
! extract information and output
!
call calc_ax_cf_omp(nrho, qmrvecs_cf(1,1), qmrvecs_cf(lrho+1,1), qmrvecs_cf(1,9), qmrvecs_cf(lrho+1,9))
!
lnj = 1
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         qmrvecs_cf(lnj,9)   = qmrvecs_cf(lnj,9)   + dble (brhoh(ni,nj,lni))
         qmrvecs_cf(lnj+1,9) = qmrvecs_cf(lnj+1,9) + dimag(brhoh(ni,nj,lni))
         lnj = lnj + 2
      end do
   end do
end do
lnj = lrho + 1
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         qmrvecs_cf(lnj,9)   = qmrvecs_cf(lnj,9)   + dble (brhoa(ni,nj,lni))
         qmrvecs_cf(lnj+1,9) = qmrvecs_cf(lnj+1,9) + dimag(brhoa(ni,nj,lni))
         lnj = lnj + 2
      end do
   end do
end do
dtmp2 = dnrm2_long(lnum, qmrvecs_cf(1,9), 1)
write(6,*)'solve_cf_tfqmr: final norm ', dtmp2
!
lnj = 1
do nj=1,nrho
   do ni=1,nrho
      zouth(ni,nj) = dcmplx(qmrvecs_cf(lnj,1), qmrvecs_cf(lnj+1,1))
      lnj = lnj + 2     
   end do
end do
lnj = lrho + 1
do nj=1,nrho
   do ni=1,nrho
      zouta(ni,nj) = dcmplx(qmrvecs_cf(lnj,1), qmrvecs_cf(lnj+1,1))
      lnj = lnj + 2     
   end do
end do
!
100 continue
! 
! check "rho" (exploit sparsity of brho)
if (lchkdos) then
   lnj = 1
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            brhoh(ni,nj,lni) = dcmplx(qmrvecs_cf(lnj,1), qmrvecs_cf(lnj+1,1))
            lnj = lnj + 2
         end do
      end do
   end do
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            brhoa(ni,nj,lni) = dcmplx(qmrvecs_cf(lnj,1), qmrvecs_cf(lnj+1,1))
            lnj = lnj + 2
         end do
      end do
   end do
end if
!
call cpu_time(cpu4)
!
write(6,*)'solve_cf_tfqmr: leaving solve_cf_tfqmr ', cpu4 - cpu3
call flush(6)
end subroutine solve_cf_tfqmr
!
