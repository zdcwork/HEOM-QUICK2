subroutine initialize_ssp(length0)
use matmod
use sparsemod
implicit none
!
include '../include/sizes'
include '../include/common'
!
integer*8, intent(in) :: length0(*)
integer               :: istat
integer               :: ni, nj, nread
integer               :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
integer               :: nnz
integer*8             :: lni, lnj, lnk, lunk
logical               :: lexist
real*8, allocatable   :: dmat1(:,:)
real*8                :: dtmp1
!
if (.not. lssp) then
    write(6,*)
    write(6,*)'initialize_ssp: error! '
    write(6,*)'initialize_ssp: lssp=F found, wrong entry '
    stop
end if
if (.not. lsparse) then
    write(6,*)
    write(6,*)'initialize_ssp: lssp=T & lsparse=F found '
    write(6,*)'initialize_ssp: code not ready yet '
    stop
end if
!
inquire(file='rho_spa.per', exist=lexist, err=990)
990 continue
if (.not. lexist) then
    write(6,*)
    write(6,*)'initialize_ssp: error! '
    write(6,*)'initialize_ssp: file <rho_spa.per> does not exist! '
    stop
end if
if (nstep_ssp .eq. 2) then
    inquire(file='rhs_spa.per', exist=lexist, err=991)
    991 continue
    if (.not. lexist) then
        write(6,*)
        write(6,*)'initialize_ssp: error! '
        write(6,*)'initialize_ssp: file <rhs_spa.per> does not exist! '
        stop
    end if
    allocate(rhs_spa(lunk_spa+1), STAT=istat)
end if
!
allocate(rsdm(norbs,norbs,nspin), STAT=istat)
!
! Initialize density matrix (and auxiliary density matrix) here, and
! note that normalization constraint for rho^0 should be satisfied at this stage.
!
open(unit=73, file='rho_spa.per', form='binary', status='unknown')
rewind(73)
if (nstep_ssp .eq. 2) then
    open(unit=74, file='rhs_spa.per', form='binary', status='unknown')
    rewind(74)
end if
!
read(73)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5  ! corresponding to norbs, nspin, ncor, jtier, and nalf, respectively
read(73)lunk                               ! number of unknowns up to jtier (jtier == ntier0)
!
if (ntmp1 .ne. norbs .or. ntmp2 .ne. nspin .or. ntmp3 .ne. ncor .or. &
    ntmp4 .ne. ntier0 .or. ntmp5 .ne. nalf) then 
    write(6,*)'initialize_ssp: error input file <rho_spa.per> '
    write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
    write(6,*)norbs, nspin, ncor, ntier0, nalf
    stop
end if
!
if (iread_spa .eq. 1) then
    ! note added on 2021-11-26
    ! reference and present job should have exactly same sparsity pattern
    ! so report error if iread_spa != 2
    deallocate(ind_spa0, nnz_spa0, irow_spa0, icol_spa0, STAT=istat)
    write(6,*)
    write(6,*)'initialize_ssp: error! '
    write(6,*)'initialize_ssp: iread_spa = 1 found '
    write(6,*)'initialize_ssp: sparsity pattern must be identical '
    stop
end if
!
lnj = 0
do ni=1,ntmp4
   lnj = lnj + length0(ni)
end do
if (lnj .ne. nunk .or. lunk .ne. nunk) then
   write(6,*)'initialize_ssp: error! incompatible <rho_spa.per> '
   write(6,*)lnj, lunk, nunk
   stop
end if
do lni=1,nunk
   nnz = nnz_spa(lni)
   lnk = ind_spa(lni)
   if (nnz .gt. 0) then
      read(73) (rho_spa(lnk-1+ni), ni=1,nnz)
   end if
end do
close(73)
write(6,*)'initialize_ssp: <rho_spa.per> read successfully '
call flush(6)
call refresh_rhosys
write(6,*)
write(6,*)'initialize_ssp: rhosys after read'
call flush(6)
allocate(dmat1(nrho,nrho), STAT=istat)
call cmatout(nrho, nrho, rhosys, dmat1)
deallocate(dmat1, STAT=istat)
!
if (nstep_ssp .eq. 2) then
    read(74)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5  ! corresponding to norbs, nspin, ncor, jtier, and nalf, respectively
    read(74)lunk                               ! number of unknowns up to jtier (jtier == ntier0)
    if (ntmp1 .ne. norbs .or. ntmp2 .ne. nspin .or. ntmp3 .ne. ncor .or. &
        ntmp4 .ne. ntier0 .or. ntmp5 .ne. nalf) then 
        write(6,*)'initialize_ssp: error input file <rhs_spa.per> '
        write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
        write(6,*)norbs, nspin, ncor, ntier0, nalf
        stop
    end if
    if (lunk .ne. nunk) then
       write(6,*)'initialize_ssp: error! incompatible <rhs_spa.per> '
       write(6,*)lunk, nunk
       stop
    end if
    do lni=1,nunk
       nnz = nnz_spa(lni)
       lnk = ind_spa(lni)
       if (nnz .gt. 0) then
          read(74) (rhs_spa(lnk-1+ni), ni=1,nnz)
       end if
    end do
    read(74)rhs_spa(lunk_spa+1)
    close(74)
    write(6,*)'initialize_ssp: <rhs_spa.per> read successfully '
    dtmp1 = 0.d0
    do lni=1,lunk_spa+1
       dtmp1 = dtmp1 + dble(rhs_spa(lni))**2 + dimag(rhs_spa(lni))**2
    end do
    dtmp1 = dsqrt(dtmp1)
    write(6,*)'initialize_ssp: norm of rhs_spa: ', dtmp1
    call flush(6)
end if
!
! Initialize system Hamiltonian hs as zero
!
hs(1:nrho,1:nrho) = czero
!
return
end subroutine initialize_ssp
