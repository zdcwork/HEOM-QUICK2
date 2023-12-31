subroutine solve_cf_tfqmr_spa(rhoinp, isgnb, ispinb, iorbsb, isgnw, zouth, zouta)
! Transpose-Free Quasi-Minimal Residue method 
! for system correlation functions
! and system spectral functions by 
! solving linear problem A * x = b
use omp_lib
use matmod
use matmod_omp
use corrfuncmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent (in) :: rhoinp(*)
complex*16, intent (out):: zouth(nrho,*), zouta(nrho,*)
integer, intent(in)     :: isgnb, ispinb, iorbsb, isgnw
! isgnw = 1: +i * omega in the linear problem
! isgnw = 2: -i * omega in the linear problem
!
logical                 :: lexist
integer                 :: ni, nj, nk, nl
integer                 :: istat, iter, ierr, maxiter, iter0, ispin, iorbs, ialf
integer, parameter      :: iter_show = 50
integer                 :: info0(4) 
integer                 :: revcom, colx, colb
integer*8               :: lni, lnj, lnk, lnum, lrho
integer                 :: iop0, nnz, nnz2
real*8                  :: dtmp1, dtmp2
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: dprev0, derr0
real*8                  :: dlamch, cdotlong, dnrm2_long
real*8                  :: occspin(2)
character*1             :: transa, transb
complex*16, allocatable :: cvec1(:)
complex*16, allocatable :: cmat1(:,:), cmat2(:,:)
external dlamch, cdotlong, dnrm2_long
!
call cpu_time(cpu3)
write(6,*)
write(6,*)'solve_cf_tfqmr_spa: entering solve_cf_tfqmr_spa, nrho ', nrho
call flush(6)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'solve_cf_tfqmr_spa: error entry ', lsparse
   stop
end if
!
iter = 0
lrho = lunkcf_spa * 2       ! factor of 2 for two types of symmetries: brhoh and brhoa
lnum = lrho * 2             ! factor of 2 for real and imaginary parts 
write(6,*)'solve_cf_tfqmr_spa: length of TFQMR vector ', lnum 
dtmp1 = dble(lrho * 16) / dble(1024**2)
!
dprev0 = crit_dos       ! TFQMR working-subroutine uses residue norm
                        ! as its convergence criterion
dprev0 = dprev0 * 1.d-3 ! use more rigorous criterion inside 'dutfx'
!
iter0  = maxit_dos / 2  ! TFQMR has two iteration numbers
maxiter = iter0
!
! lplus_cf enters as a flag in <_cf>
if (isgnw .eq. 1) then
   lplus_cf = .true.
else if (isgnw .eq. 2) then
   lplus_cf = .false.
else 
   write(6,*)'solve_cf_tfqmr_spa: error! unknown isgnw = ', isgnw
   stop
end if
lplus = lplus_cf
!
! initialize
!
allocate(cvec1(nrho**2), STAT=istat)
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), STAT=istat)
cvec1(1:nrho**2)     = czero
cmat1(1:nrho,1:nrho) = czero
cmat2(1:nrho,1:nrho) = czero
!
lni = ind_spa(1)
nnz = nnz_spa(1)
call zmat_coo2dns(nnz, irow_spa(lni), icol_spa(lni), rhoinp(lni), cmat1, nrho, nrho, nrho)
call calchs(cmat1)
if (igroundsteady .eq. 0) then          ! ground state
    dhs(1:nrho,1:nrho) = czero 
else if (igroundsteady .eq. 1) then     ! steady state
    call calcdhs(1.d10, cmat1)
else if (igroundsteady .eq. 2) then     ! transient time
    call calcdhs(tnow, cmat1)
!   continue ! -- to check, in principle the current code only deals with static cases
!   if (lad) then
!       write(6,*)
!       write(6,*)'solve_cf_tfqmr_spa: lad=T and igroundsteady=2 found '
!       write(6,*)'solve_cf_tfqmr_spa: code unavailable yet '
!       stop
!   end if
else 
   write(6,*)'solve_cf_tfqmr_spa: error! unknown igroundsteady = ', igroundsteady
   stop
end if
!
if (ltrun_der) call prelude_cf_trun
!
if (lad) then
    if ( (.not. allocated(zgefnl)) .or. (.not. allocated(zgefnr)) .or.   &
         (.not. allocated(zgelvl)) .or. (.not. allocated(zgelvr)) ) then
        write(6,*)
        write(6,*)'solve_cf_tfqmr_spa: call prepare_ad '
        call flush(6)
        call prepare_ad
    end if
end if
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
   write(6,*)'solve_cf_tfqmr_spa: HF with TFQMR does not work. Abort...'
   stop
end if
!
! Bh^- = B * rho^- + rho^- * B^+ = (Bh^+)^dag
! Ba^- = (-i) * (B * rho^- - rho^- * B^+) = (Ba^+)^dag 
!
brhoh_spa(1:lunkcf_spa) = czero
brhoa_spa(1:lunkcf_spa) = czero
if (isgnb .eq. 1) then
   transa = 'c'
   transb = 'n'
else if (isgnb .eq. 2) then
   transa = 'n'
   transb = 'c'
else 
   write(6,*)
   write(6,*)'solve_cf_tfqmr_spa: error! unknown isgnb ', isgnb
   stop
end if
do lni=1,nunk
   if ( nnz_spa(lni) .eq. 0 .or. nnzcf_spa(lni) .eq. 0 ) cycle
   lnj  = ind_spa(lni) 
   nnz  = nnz_spa(lni)
   lnk  = indcf_spa(lni)
   nnz2 = nnzcf_spa(lni)
! Bh
   if (lproj) then
       call zamsmm_coo_proj('l', transa, 'n', iorbsb, ispinb, cunity, rhoinp(lnj), nnz,             &
                            irow_spa(lnj), icol_spa(lnj),  czero, cmat1, nrho, cmat2, nrho)
       call zamsmm_coo_proj('r', transb, 'n', iorbsb, ispinb, cunity, rhoinp(lnj), nnz,             &
                            irow_spa(lnj), icol_spa(lnj), cunity, cmat1, nrho, cmat2, nrho)
   else
       call zamsmm_coo('l', transa, 'n', iorbsb, ispinb, cunity, rhoinp(lnj), nnz,             &
                       irow_spa(lnj), icol_spa(lnj),  czero, cmat1, nrho, cmat2, nrho)
       call zamsmm_coo('r', transb, 'n', iorbsb, ispinb, cunity, rhoinp(lnj), nnz,             &
                       irow_spa(lnj), icol_spa(lnj), cunity, cmat1, nrho, cmat2, nrho)
   end if
   call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), brhoh_spa(lnk), cmat1, nrho, nrho, nrho)
! Ba
   if (lproj) then
       call zamsmm_coo_proj('l', transa, 'n', iorbsb, ispinb, -eye  , rhoinp(lnj), nnz,             &
                            irow_spa(lnj), icol_spa(lnj),  czero, cmat1, nrho, cmat2, nrho)  
       call zamsmm_coo_proj('r', transb, 'n', iorbsb, ispinb,  eye  , rhoinp(lnj), nnz,             &
                            irow_spa(lnj), icol_spa(lnj), cunity, cmat1, nrho, cmat2, nrho)  
   else
       call zamsmm_coo('l', transa, 'n', iorbsb, ispinb, -eye  , rhoinp(lnj), nnz,             &
                       irow_spa(lnj), icol_spa(lnj),  czero, cmat1, nrho, cmat2, nrho)  
       call zamsmm_coo('r', transb, 'n', iorbsb, ispinb,  eye  , rhoinp(lnj), nnz,             &
                       irow_spa(lnj), icol_spa(lnj), cunity, cmat1, nrho, cmat2, nrho)  
   end if
   call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), brhoa_spa(lnk), cmat1, nrho, nrho, nrho)
end do
!
! A * x = b, where b is (-brhoh, -brhoa)
! 
lnj = 1
do lni=1,lunkcf_spa
   qmrvecs_cf(lnj,   2) = -dble (brhoh_spa(lni))
   qmrvecs_cf(lnj+1, 2) = -dimag(brhoh_spa(lni))
   lnj = lnj + 2
end do
lnj = lrho + 1
do lni=1,lunkcf_spa
   qmrvecs_cf(lnj,   2) = -dble (brhoa_spa(lni))
   qmrvecs_cf(lnj+1, 2) = -dimag(brhoa_spa(lni))
   lnj = lnj + 2
end do
!
lexist = .false.
if (lread_dos) then
    if (isgnb .eq. 1) then
        write(6,*)'solve_cf_tfqmr_spa: inquiring file <rho1_spa_cf.sav> '
        inquire(file='rho1_spa_cf.sav', exist=lexist)
        if (lexist) then
            write(6,*)'solve_cf_tfqmr_spa: trying to read <rho1_spa_cf.sav> '
            open(unit=75, file='rho1_spa_cf.sav', form='binary', status='old')
        end if
    else if (isgnb .eq. 2) then
        write(6,*)'solve_cf_tfqmr_spa: inquiring file <rho2_spa_cf.sav> '
        inquire(file='rho2_spa_cf.sav', exist=lexist)
        if (lexist) then
            write(6,*)'solve_cf_tfqmr_spa: trying to read <rho2_spa_cf.sav> '
            open(unit=75, file='rho2_spa_cf.sav', form='binary', status='old')
        end if
    else
        write(6,*)'solve_cf_tfqmr_spa: error! unknown isgnb ', isgnb
        stop
    end if
    if (lexist) then
        rewind(75)
        read(75)lni
        if (lni .eq. lnum) then
            read(75)(qmrvecs_cf(lni,1), lni=1,lnum)
            write(6,*)'solve_cf_tfqmr_spa: file loaded succesfully '
        else
            lexist = .false.
            write(6,*)'solve_cf_tfqmr_spa: incompatible file, abort reading '
        end if
        call flush(6)
        close(75)
    else 
        write(6,*)'solve_cf_tfqmr_spa: file does not exist '
        call flush(6)
    end if
    if (lexist) then
        call calc_ax_cf_spa_omp(lunkcf_spa, qmrvecs_cf(1,1), qmrvecs_cf(lrho+1,1), qmrvecs_cf(1,9), qmrvecs_cf(lrho+1,9))
        qmrvecs_cf(1:lnum,2) = qmrvecs_cf(1:lnum,2) - qmrvecs_cf(1:lnum,9)
    end if
end if
!
dtmp2 = dnrm2_long(lnum, qmrvecs_cf(1,2), 1)
write(6,*)
write(6,*)'solve_cf_tfqmr_spa: initial norm ', dtmp2
if (dtmp2 .le. crit_dos) then
   write(6,*)'solve_cf_tfqmr_spa: initial norm smaller than criterion ', crit_dos
   write(6,*)'solve_cf_tfqmr_spa: no need to go through TFQMR process '
   goto 100
end if
!
if (lexist) then
    qmrvecs_cf(1:lnum,1) = 0.d0
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
     call calc_ax_cf_spa_omp(lunkcf_spa, qmrvecs_cf(1,1), qmrvecs_cf(lrho+1,1), qmrvecs_cf(1,9), qmrvecs_cf(lrho+1,9))
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
  call calc_ax_cf_spa_omp(lunkcf_spa, qmrvecs_cf(1,colx), qmrvecs_cf(lrho+1,colx), qmrvecs_cf(1,colb), qmrvecs_cf(lrho+1,colb))
  call flush(6)
end do
!
if (maxiter > iter0) then
   write(6,*)'solve_cf_tfqmr_spa: error! maxiter = ', maxiter
   write(6,*)'solve_cf_tfqmr_spa: TFQMR procedures NOT converge after', iter, ' iterations'
else 
   write(6,*)'solve_cf_tfqmr_spa: TFQMR converged after', iter, ' iterations'
end if
ierr = info0(1)
if (iop0 .eq. 0) then
   write(6,*)
   write(6,*)'solve_cf_tfqmr_spa: ierr = ', ierr
   if (ierr .eq. 0) then
      write(6,*)'solve_cf_tfqmr_spa: TFQMR converged succesfully'
   else if (ierr .eq. 1) then
      write(6,*)'solve_cf_tfqmr_spa: error! invalid reverse communication call'
   else if (ierr .eq. 2) then
      write(6,*)'solve_cf_tfqmr_spa: error! invalid input for TFQMR'
   else if (ierr .eq. 4) then
      write(6,*)'solve_cf_tfqmr_spa: TFQMR failed to converge after ', iter, ' steps'
   else if (ierr .eq. 8) then
      write(6,*)'solve_cf_tfqmr_spa: TFQMR ran into an exact breakdown'
   else 
      write(6,*)'solve_cf_tfqmr_spa: unknown ierr flag, ierr = ', ierr
   end if
end if
call flush(6)
!
! extract information and output
!
if (lexist) then
    if (isgnb .eq. 1) then
        open(unit=75, file='rho1_spa_cf.sav', form='binary', status='old')
    else if (isgnb .eq. 2) then
        open(unit=75, file='rho2_spa_cf.sav', form='binary', status='old')
    else
        write(6,*)'solve_cf_tfqmr_spa: error! unknown isgnb ', isgnb
        stop
    end if
    rewind(75)
    read(75)lni
    read(75)(qmrvecs_cf(lni,2), lni=1,lnum)
    close(75)
    qmrvecs_cf(1:lnum,1) = qmrvecs_cf(1:lnum,1) + qmrvecs_cf(1:lnum,2)
end if
!
call calc_ax_cf_spa_omp(lunkcf_spa, qmrvecs_cf(1,1), qmrvecs_cf(lrho+1,1), qmrvecs_cf(1,9), qmrvecs_cf(lrho+1,9))
!
lnj = 1
do lni=1,lunkcf_spa
   qmrvecs_cf(lnj,9)   = qmrvecs_cf(lnj,9)   + dble (brhoh_spa(lni))
   qmrvecs_cf(lnj+1,9) = qmrvecs_cf(lnj+1,9) + dimag(brhoh_spa(lni))
   lnj = lnj + 2
end do
lnj = lrho + 1
do lni=1,lunkcf_spa
   qmrvecs_cf(lnj,9)   = qmrvecs_cf(lnj,9)   + dble (brhoa_spa(lni))
   qmrvecs_cf(lnj+1,9) = qmrvecs_cf(lnj+1,9) + dimag(brhoa_spa(lni))
   lnj = lnj + 2
end do
dtmp2 = dnrm2_long(lnum, qmrvecs_cf(1,9), 1)
write(6,*)'solve_cf_tfqmr_spa: final norm ', dtmp2
call flush(6)
!
100 continue
!
nnz = nnzcf_spa(1)
lni = 1
do ni=1,nnz
   cvec1(ni) = dcmplx( qmrvecs_cf(lni,1), qmrvecs_cf(lni+1,1) )
   lni = lni + 2
end do
lnj = indcf_spa(1) 
call zmat_coo2dns(nnz, irowcf_spa(lnj), icolcf_spa(lnj), cvec1, zouth, nrho, nrho, nrho)
!
lni = lrho + 1
do ni=1,nnz
   cvec1(ni) = dcmplx( qmrvecs_cf(lni,1), qmrvecs_cf(lni+1,1) )
   lni = lni + 2
end do
lnj = indcf_spa(1) 
call zmat_coo2dns(nnz, irowcf_spa(lnj), icolcf_spa(lnj), cvec1, zouta, nrho, nrho, nrho)
!
! output to file
!
if (isgnb .eq. 1) then
    open(unit=75, file='rho1_spa_cf.sav', form='binary', status='unknown')
else if (isgnb .eq. 2) then
    open(unit=75, file='rho2_spa_cf.sav', form='binary', status='unknown')
else
    write(6,*)'solve_cf_tfqmr_spa: error! unknown isgnb ', isgnb
    stop
end if
rewind(75)
write(75)lnum
write(75)(qmrvecs_cf(lni,1), lni=1,lnum)
close(75)
!
deallocate(cvec1, STAT=istat)
deallocate(cmat1, cmat2, STAT=istat)
!
call cpu_time(cpu4)
!
write(6,*)'solve_cf_tfqmr_spa: leaving solve_cf_tfqmr_spa ', cpu4 - cpu3
call flush(6)
end subroutine solve_cf_tfqmr_spa
!
