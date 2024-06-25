subroutine solve_tfqmr
! Transpose-Free Quasi-Minimal Residue method for  
! solving linear problem A * x = b
use omp_lib
use matmod_omp
use tfqmrmod
use matmod
use tmpmatmod
use sparsemod
use hbmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, nl
integer                 :: istat, iter, ierr, maxiter, iter0, ispin, iorbs, ibath, imode
integer                 :: nnz
integer, parameter      :: iter_show = 50
integer                 :: info0(4) 
integer                 :: revcom, colx, colb
integer*8               :: lni, lnj, lnk, lnum, lunk
integer                 :: iop0
integer                 :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
real*8                  :: dtmp1, dtmp2
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: dprev0, derr0, derrmin
real*8                  :: dlamch, cdotlong, dnrm2_long
real*8                  :: occspin(2)
complex*16              :: ctmp1, ctmp2
complex*16, allocatable :: cvec1(:)
logical                 :: fexist
external dlamch, cdotlong, dnrm2_long
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_tfqmr: entering solve_tfqmr, nrho ', nrho
!
iter = 0
if (.not. lsparse) then
   lnum = nunk * nrho**2 * 2   ! last factor 2 accounts for both real and imaginary parts
   if (lhb) then
      lnum = lnum + nunk * nrho**2 * nbath_hb * 2
   end if
else
   lnum = (lunk_spa + 1) * 2
end if
dtmp1 = dble(lnum * 8) / dble(1024**2)
write(6,*)'solve_tfqmr: length of TFQMR vector ', lnum 
call flush(6)
!
dprev0 = crit           ! TFQMR working-subroutine uses residue norm
                        ! as its convergence criterion
dprev0 = dprev0 * 1.d-3 ! use more rigorous criterion inside 'dutfx'
!
iter0  = maxit0 / 2     ! TFQMR has two iteration numbers
maxiter = iter0
derrmin = 1.d10
!
! initialize
!
info0(1:4) = 0
info0(1)   = 600        ! write iterative procedure to file channel #6
!info0(1)   = 10600       
qmrvecs    = 0.d0
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
   if (lsparse) then
      rho_spa = czero
   else
      rho = czero
   end if
   write(6,*)'solve_tfqmr: HF with TFQMR does not work. Abort...'
   stop
end if
!
if (lssp .and. nstep_ssp .eq. 2) then
    lnj = 1
    do lni=1,lunk_spa+1
       qmrvecs(lnj,2)   = dble (rhs_spa(lni))
       qmrvecs(lnj+1,2) = dimag(rhs_spa(lni))
       lnj = lnj + 2
    end do
    qmrvecs(1:lnum,2) = -qmrvecs(1:lnum,2)
    dtmp2 = dnrm2_long(lnum, qmrvecs(1,2), 1)
    goto 102
end if
!
if (.not. lsparse) then
   do ni=1,nrho
      do nj=ni+1,nrho
         rho(ni,nj,1) = dconjg(rho(nj,ni,1))
      end do
      rho(ni,ni,1) = dcmplx(dble(rho(ni,ni,1)), 0.d0)
   end do
   lnj = 1
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            qmrvecs(lnj,1)   = dble (rho(ni,nj,lni))
            qmrvecs(lnj+1,1) = dimag(rho(ni,nj,lni))
            lnj = lnj + 2
         end do 
      end do
   end do
   if (lhb) then
      do ibath=1,nbath_hb
         do lni=1,nunk
            do nj=1,nrho
               do ni=1,nrho
                  qmrvecs(lnj,1)   = dble (rho_hb(ni,nj,lni,ibath))
                  qmrvecs(lnj+1,1) = dimag(rho_hb(ni,nj,lni,ibath))
                  lnj = lnj + 2
               end do
            end do
         end do
      end do
   end if
   call calc_ax_tfqmr_omp(nrho, qmrvecs(1,1), qmrvecs(1,2))
   if (lhb) then
      do ibath=1,nbath_hb
         call calc_ax_hb_omp(ibath, nrho, qmrvecs(1,1), qmrvecs(1,2))
      end do
   end if   
   qmrvecs(1:lnum,1) = 0.d0
else 
   nnz = nnz_spa(1)
   lni = ind_spa(1)
!   if (.not. lreadi) then
   call zmat_herm_coo('u', nnz, rho_spa(lni), irow_spa(lni), icol_spa(lni), cmtmp1, nrho, nrho)
!   end if
   call calc_ax_tfqmr_spa_omp(lunk_spa, rho_spa, qmrvecs(1,2))
end if
call refresh_rhosys
write(6,*)
write(6,*)'solve_tfqmr: rhosys before tfqmr cycles '
call flush(6)
call cmatout(nrho, nrho, rhosys, dmtmp1)
!
if (lssp .and. nstep_ssp .eq. 1) then
    qmrvecs(lnum-1,2) = qmrvecs(lnum-1,2) - 1.d0
    dtmp2 = dnrm2_long(lnum, qmrvecs(1,2), 1)
    write(6,*)
    write(6,*)'solve_tfqmr: initial norm ', dtmp2
    call flush(6)
    open(unit=74, file='rhs_spa.per', form='binary', status='unknown')
    rewind(74)
    write(74)norbs, nspin, ncor, ntier0, nalf
    write(74)nunk
    do lni=1,nunk
       nnz = nnz_spa(lni)
       lnj = ind_spa(lni) * 2 - 2
       if (nnz .gt. 0) then
           write(74) (dcmplx(qmrvecs(lnj+ni*2-1,2),qmrvecs(lnj+ni*2,2)), ni=1,nnz)
       end if
    end do
    write(74)dcmplx(qmrvecs(lnum-1,2), qmrvecs(lnum,2))
    close(74)
    write(6,*)'solve_tfqmr: <rhs_spa.per> updated '
    call flush(6)
    goto 101
end if
!
qmrvecs(1:lnum,2) = -qmrvecs(1:lnum,2)
!if (ltrend) then
!    qmrvecs(lnum-1,2) = 1.d0 + qmrvecs(lnum-1,2)   ! last element of B vector is trace(rho) = 1
!else
!    qmrvecs(1,2) = 1.d0 + qmrvecs(1,2)             ! first element of B vector is trace(rho) = 1
!end if
qmrvecs(lnum-1,2) = 1.d0 + qmrvecs(lnum-1,2)   ! last element of B vector is trace(rho) = 1
dtmp2 = dnrm2_long(lnum, qmrvecs(1,2), 1)
!
102 continue
!
write(6,*)
write(6,*)'solve_tfqmr: initial norm ', dtmp2
if (dtmp2 .le. crit) then
   write(6,*)'solve_tfqmr: initial norm smaller than criterion ', crit
   write(6,*)'solve_tfqmr: no need to go through TFQMR process '
   goto 100
end if
!
! Transpose-Free Quasi-Minimal Residue Algorithm
!
iop0 = 0
do 
  call dutfx(lnum, lnum, maxiter, qmrvecs, dprev0, info0)
  revcom = info0(2)
  colx   = info0(3)
  colb   = info0(4)
!  
  if (mod(iter, iter_show) .eq. 0) then
     if (lsparse) then
        call calc_ax_tfqmr_spa_omp(lunk_spa, qmrvecs(1,1), qmrvecs(1,9))
     else
        call calc_ax_tfqmr_omp(nrho, qmrvecs(1,1), qmrvecs(1,9))
        if (lhb) then
           do ibath=1,nbath_hb
              call calc_ax_hb_omp(ibath, nrho, qmrvecs(1,1), qmrvecs(1,9))
           end do
        end if   
     end if
     qmrvecs(1:lnum,9) = qmrvecs(1:lnum,9) - qmrvecs(1:lnum,2)
     derr0 = dnrm2_long(lnum, qmrvecs(1,9), 1)
     write(6,*)'iter_tfqmr ', iter, derr0
     if (derr0 .le. crit) then
        iop0 = 1
        exit
     end if
!
     if (lsparse .and. derr0 .lt. derrmin) then
         derrmin = derr0
         allocate(cvec1(nrho**2), STAT=istat)
         open(unit=71, file='rho_spa.chk', status='unknown', form='binary')
         rewind(71)
         write(71)norbs, nspin, ncor, ntier0, nalf
         write(71)nunk
         lnk = 1
         do lni=1,nunk
            nnz = nnz_spa(lni)
            lnj = ind_spa(lni)
            if (nnz .gt. 0) then
                do ni=1,nnz
                   cvec1(ni) = rho_spa(lnj-1+ni) + dcmplx(qmrvecs(lnk,1), qmrvecs(lnk+1,1))
                   lnk = lnk + 2
                end do
                write(71)(cvec1(ni), ni=1,nnz)
            end if
         end do
         close(71)
         deallocate(cvec1, STAT=istat)
         write(6,*)'solve_tfqmr: the file <rho_spa.chk> updated '
         call flush(6)
     end if
!
     !fexist = .false.
     !inquire(file="stopst", exist=fexist)  ! exit do 
     !if (fexist) then
     !   iop0 = 2
     !   write(6,*)'solve_tfqmr: file stopst detected, exit from TFQMR solver'
     !   exit
     !end if
  end if
!
  fexist = .false.
  inquire(file="stopst", exist=fexist)  ! exit do 
  if (fexist) then
     iop0 = 2
     write(6,*)'solve_tfqmr: file stopst detected, exit from TFQMR solver'
     exit
  end if
!
  if (revcom .ne. 1) exit
  iter = iter + 1
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_tfqmr: before ax', dtmp1, dtmp2
  if (lsparse) then 
     call calc_ax_tfqmr_spa_omp(lunk_spa, qmrvecs(1,colx), qmrvecs(1,colb))
  else
     call calc_ax_tfqmr_omp(nrho, qmrvecs(1,colx), qmrvecs(1,colb))
     if (lhb) then
        do ibath=1,nbath_hb
           call calc_ax_hb_omp(ibath, nrho, qmrvecs(1,colx), qmrvecs(1,colb))
        end do
     end if   
  end if
!dtmp1 = dnrm2_long(lnum, qmrvecs(1,colx), 1)
!dtmp2 = dnrm2_long(lnum, qmrvecs(1,colb), 1)
!write(6,*)'solve_tfqmr: after  ax', dtmp1, dtmp2
  call flush(6)
end do
!
if (maxiter > iter0) then
   write(6,*)'solve_tfqmr: error! maxiter = ', maxiter
   write(6,*)'solve_tfqmr: TFQMR procedures NOT converge after', iter, ' iterations'
else 
   if (iop0 .eq. 1) then
      write(6,*)'solve_tfqmr: TFQMR converged after', iter, ' iterations'
   else if (iop0 .eq. 2) then
      write(6,*)'solve_tfqmr: TFQMR stopped   after', iter, ' iterations'
   end if
end if
ierr = info0(1)
if (iop0 .eq. 0) then
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
end if
call flush(6)
!
! recover rho
!
if (lsparse) then
   lnj = 1
   do lni=1,lunk_spa
      rho_spa(lni) = rho_spa(lni) + dcmplx(qmrvecs(lnj,1), qmrvecs(lnj+1,1))
      lnj = lnj + 2
   end do
   call refresh_rhosys
   !call calc_ax_tfqmr_spa_omp(lunk_spa, rho_spa, qmrvecs(1,9))
   if (lssp) then
       call calc_ax_tfqmr_spa_omp(lunk_spa, qmrvecs(1,1), qmrvecs(1,9))
       qmrvecs(1:lnum,9) = qmrvecs(1:lnum,9) - qmrvecs(1:lnum,2)
   else
       call calc_ax_tfqmr_spa_omp(lunk_spa, rho_spa, qmrvecs(1,9))
   end if
else
   lnj = 1
   do lni=1,nunk
      do nj=1,nrho
         do ni=1,nrho
            rho(ni,nj,lni) = rho(ni,nj,lni) + dcmplx(qmrvecs(lnj,1), qmrvecs(lnj+1,1))
            qmrvecs(lnj,1) = dble(rho(ni,nj,lni))
            qmrvecs(lnj+1,1) = dimag(rho(ni,nj,lni))
            lnj = lnj + 2
         end do
      end do
   end do
   if (lhb) then
      do ibath=1,nbath_hb
         do lni=1,nunk
            do nj=1,nrho
               do ni=1,nrho
                  rho_hb(ni,nj,lni,ibath) = rho_hb(ni,nj,lni,ibath) + dcmplx(qmrvecs(lnj,1), qmrvecs(lnj+1,1))
                  qmrvecs(lnj,1) = dble(rho_hb(ni,nj,lni,ibath))
                  qmrvecs(lnj+1,1) = dimag(rho_hb(ni,nj,lni,ibath))
                  lnj = lnj + 2
               end do
            end do
         end do
      end do
   end if
   call calc_ax_tfqmr_omp(nrho, qmrvecs(1,1), qmrvecs(1,9))
   if (lhb) then
      do ibath=1,nbath_hb
         call calc_ax_hb_omp(ibath, nrho, qmrvecs(1,1), qmrvecs(1,9))
      end do
   end if   
end if
!
if (lssp) then
    !lnj = 1
    !do lni=1,lunk_spa+1
    !   qmrvecs(lnj,9)   = qmrvecs(lnj,9)   + dble (rhs_spa(lni)) 
    !   qmrvecs(lnj+1,9) = qmrvecs(lnj+1,9) + dimag(rhs_spa(lni))
    !   lnj = lnj + 2
    !end do
    !qmrvecs(lnum-1,9) = qmrvecs(lnum-1,9) - 1.d0
    open(unit=73, file='rho_spa.per', form='binary', status='unknown')
    rewind(73)
    write(73)norbs, nspin, ncor, ntier0, nalf
    write(73)nunk
    do lni=1,nunk
       nnz = nnz_spa(lni)
       lnj = ind_spa(lni) 
       if (nnz .gt. 0) then
           write(73) (rho_spa(lnj-1+ni), ni=1,nnz)
       end if
    end do
    close(73)
    write(6,*)'solve_tfqmr: <rho_spa.per> updated '
    call flush(6)
else
    qmrvecs(lnum-1,9) = qmrvecs(lnum-1,9) - 1.d0
!    if (ltrend) then
!        qmrvecs(lnum-1,9) = qmrvecs(lnum-1,9) - 1.d0
!    else
!        qmrvecs(1,9) = qmrvecs(1,9) - 1.d0
!    end if
end if
dtmp2 = dnrm2_long(lnum, qmrvecs(1,9), 1)
write(6,*)'solve_tfqmr: final norm ', dtmp2
!
100 continue
!
!if (lssp .and. nstep_ssp .eq. 2) then
!    lnj = 1
!    do lni=1,lunk_spa
!       rho_spa(lni) = rho_spa(lni) + dcmplx(qmrvecs(lnj,1), qmrvecs(lnj+1,1))
!       lnj = lnj + 2
!    end do
!    call refresh_rhosys
!end if
! 
! important: recover the upper half of rho by Hermicity
!  
! Post-TFQMR analysis begins here
! 
call outsteady(iter)
!
101 continue
!
call cpu_time(cpu4)
write(6,*)
write(6,*)'solve_tfqmr: leaving solve_tfqmr ', cpu4 - cpu3
call flush(6)
end subroutine solve_tfqmr
!
