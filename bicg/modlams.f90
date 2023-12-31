subroutine modlams(korbs, kspin, ksgn, klr, knt, iop)
!
! treat Liouville matrix elements corresponding to the upper-half of 
! rho in consideration of the Hermicity of reduced density matrix
!
use auxmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! iop = 0, real part (column transfer with sign preserved)
!     = 1, imag part (column transfer with sign inversed)
!     = 2, row removal
!
integer, intent(in) :: iop, korbs, kspin, ksgn, klr, knt
integer             :: ni, nj, nk, nl
!
call lcpams(korbs, kspin, ksgn, klr, knt)
!
do nj=1,lentmp1
  if (iop .eq. 0 .or. iop .eq. 1) then
    call jmod(lctmp1(nj), nrho, nk, nl)    
    if (nk .lt. nl) then
      if (iop .eq. 0) then                ! see 'modifylmat'
        lctmp1(nj) = (nk - 1) * nrho + nl
      else                                ! iop .eq. 1 
        lctmp1(nj) = (nk - 1) * nrho + nl
        lvtmp1(nj) = -lvtmp1(nj)
      end if
    end if
  elseif (iop .eq. 2) then
    call jmod(lrtmp1(nj), nrho, nk, nl)
    if (nk .lt. nl) then
      lvtmp1(nj) = 0.d0
    end if
  end if
end do
end subroutine modlams
!
subroutine modlams_omp(korbs, kspin, ksgn, klr, knt, iop, ithread)
!
! treat Liouville matrix elements corresponding to the upper-half of 
! rho in consideration of the Hermicity of reduced density matrix
!
use auxmod
use auxmod_omp
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! iop = 0, real part (column transfer with sign preserved)
!     = 1, imag part (column transfer with sign inversed)
!     = 2, row removal
!
integer, intent(in) :: iop, korbs, kspin, ksgn, klr, knt, ithread
integer             :: ni, nj, nk, nl
!
call lcpams_omp(korbs, kspin, ksgn, klr, knt, ithread)
!
do nj=1,lentmp1_omp(ithread)
  if (iop .eq. 0 .or. iop .eq. 1) then
    call jmod(lctmp1_omp(nj,ithread), nrho, nk, nl)    
    if (nk .lt. nl) then
      if (iop .eq. 0) then                ! see 'modifylmat'
        lctmp1_omp(nj,ithread) = (nk - 1) * nrho + nl
      else                                ! iop .eq. 1 
        lctmp1_omp(nj,ithread) = (nk - 1) * nrho + nl
        lvtmp1_omp(nj,ithread) = -lvtmp1_omp(nj,ithread)
      end if
    end if
  elseif (iop .eq. 2) then
    call jmod(lrtmp1_omp(nj,ithread), nrho, nk, nl)
    if (nk .lt. nl) then
      lvtmp1_omp(nj,ithread) = 0.d0
    end if
  end if
end do
end subroutine modlams_omp
!
subroutine lcpams(korbs, kspin, ksgn, klr, knt)
use auxmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn, klr, knt
integer             :: ni
ni           = lenams(korbs, kspin, ksgn, klr, knt)
lentmp1      = ni
lrtmp1(1:ni) = lrowams(1:ni, korbs, kspin, ksgn, klr, knt)
lctmp1(1:ni) = lcolams(1:ni, korbs, kspin, ksgn, klr, knt)
lvtmp1(1:ni) = lvalams(1:ni, korbs, kspin, ksgn, klr, knt)
end subroutine lcpams
!
subroutine lcpams_omp(korbs, kspin, ksgn, klr, knt, ithread)
use auxmod
use auxmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn, klr, knt, ithread
integer             :: ni
ni                        = lenams(korbs, kspin, ksgn, klr, knt)
lentmp1_omp(ithread)      = ni
lrtmp1_omp(1:ni, ithread) = lrowams(1:ni, korbs, kspin, ksgn, klr, knt)
lctmp1_omp(1:ni, ithread) = lcolams(1:ni, korbs, kspin, ksgn, klr, knt)
lvtmp1_omp(1:ni, ithread) = lvalams(1:ni, korbs, kspin, ksgn, klr, knt)
end subroutine lcpams_omp
!
