module cplqmrmod
implicit none
!
! cplqmr memory
!
! variables used by library subroutine DUCPL
!
integer :: maxpq, maxvw, dimm, mvec, nlim
integer :: ncol_dwk, ncol_idx, ncol_iwk, ncol_vecs
integer :: nrow_dwk, nrow_idx, nrow_iwk
real*8  :: dnorms(2)
!
! use double precision cplqmr method
!
real*8, allocatable :: qmrvecs(:,:)
!
! arrays
!
integer, allocatable :: iwk(:,:), idx(:,:)
real*8, allocatable :: dwk(:,:)
!
end module cplqmrmod
