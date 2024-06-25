subroutine initialize2_spa(length0)
use matmod
use sparsemod
implicit none
!
include '../include/sizes'
include '../include/common'
!
integer*8, intent(in) :: length0(*)
integer               :: istat
integer               :: ni, nj
integer               :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
integer               :: nnz
integer*8             :: lni, lnj, lnk, lunk
logical               :: lsmall 
complex*16, allocatable :: cvec1(:)
complex*16, allocatable :: cmat1(:,:)
!
allocate(rsdm(norbs,norbs,nspin), STAT=istat)
!
! Initialize density matrix (and auxiliary density matrix) here, and
! note that normalization constraint for rho^0 should be satisfied at this stage.
!
open(unit=50, file='rho_spa.sav', form='binary', status='unknown')
rewind(50)
read(50)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5  ! corresponding to norbs, nspin, ncor, jtier, and nalf, respectively
read(50)lunk                               ! number of unknowns up to jtier (jtier <= ntier0)
!
if (ntmp1 .ne. norbs .or. ntmp2 .ne. nspin .or. ntmp3 .ne. ncor .or. &
    ntmp5 .ne. nalf  .or. ntmp4 .le. 1) then 
    write(6,*)
    write(6,*)'initialize2_spa: error input file <rho_spa.sav> '
    stop
end if
!if (lunk .ne. nunk) then
!    write(6,*)'initialize2_spa: error input file <rho_spa.sav> '
!    write(6,*)' lunk, nunk: ', lunk, nunk
!    stop
!end if
!
lsmall = .false.
if (ntmp4 < ntier0) lsmall = .true.
!
if (iread_spa .eq. 1) then   ! here sparsity pattern of hs differs from saved 
    ! note added on Jan 24, 2020: (see also the note in <prelude_spa.f90>)
    ! the reference job should have same number of ADOs, and those ADOs
    ! should be more sparse than the present job, so we should have
    ! nunk == lunk .and. lunk_spa >= lunk_spa0
    if (lunk .ne. nunk .or. lunk_spa .lt. lunk_spa0) then
        write(6,*)'initialize2_spa: error! lunk != nunk found; or '
        write(6,*)'initialize2_spa:        lunk_spa < lunk_spa0   '
        write(6,*)lunk_spa, lunk_spa0
        stop
    end if
    rho_spa = czero
    allocate(cvec1(nrho*nrho), cmat1(nrho,nrho), STAT=istat)
    read0: do lni=1,nunk
       if (ind_spa0(lni)-1+nnz_spa0(lni) .gt. lunk_spa) then
           write(6,*)
           write(6,*)'initialize2_spa: error! <rho_spa.sav> too large '
           write(6,*)ind_spa0(lni)-1+nnz_spa0(lni), lunk_spa
           stop
       end if
       if (nnz_spa0(lni) .eq. 0) cycle
       cvec1 = czero
       read(50) (cvec1(ni), ni=1,nnz_spa0(lni))
       lnj = ind_spa0(lni)
       lnk = ind_spa(lni)
       call zmat_coo2dns(nnz_spa0(lni), irow_spa0(lnj), icol_spa0(lnj), cvec1,         &
                         cmat1, nrho, nrho, nrho)
       call zmat_extract_coo(nnz_spa(lni), irow_spa(lnk), icol_spa(lnk), rho_spa(lnk), &
                             cmat1, nrho, nrho, nrho)
    end do read0
    deallocate(cvec1, cmat1, STAT=istat)
    deallocate(ind_spa0, nnz_spa0, irow_spa0, icol_spa0, STAT=istat)
    goto 100
end if
!
lnj = 0
do ni=1,min0(ntmp4,ntier0)
   lnj = lnj + length0(ni)
end do
if (lsmall) then
    do lni=1,lunk
       nnz = nnz_spa(lni)
       lnk = ind_spa(lni)      
       if (nnz .gt. 0) then
          read(50) (rho_spa(lnk-1+ni), ni=1,nnz)
       end if
    end do
    rho_spa(ind_spa(lunk+1):lunk_spa) = czero
else 
   if (lnj .ne. nunk) then
      write(6,*)'initialize2_spa: error! incompatible <rho_spa.sav> ', ntmp4, ntier0, lnj, nunk
      stop
   end if
   do lni=1,nunk
      nnz = nnz_spa(lni)
      lnk = ind_spa(lni)
      if (nnz .gt. 0) then
         read(50) (rho_spa(lnk-1+ni), ni=1,nnz)
      end if
   end do
end if
!
100 continue
!
close(50)
!
! Initialize system Hamiltonian hs as zero
!
hs(1:nrho,1:nrho) = czero
!
return
end subroutine initialize2_spa
