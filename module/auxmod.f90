module auxmod
implicit none
!
! auxiliary memory
!
real*8,  save, allocatable :: lhsrnl(:,:), lhsrnr(:,:), lhsinl(:,:), lhsinr(:,:)
!
! sparse matrix (liouville space)
!
integer, save, allocatable :: lenams(:,:,:,:,:)     ! norbs,nspin,nsgn,nlr,nnt
integer, save, allocatable :: lrowams(:,:,:,:,:,:)
integer, save, allocatable :: lcolams(:,:,:,:,:,:)
real*8,  save, allocatable :: lvalams(:,:,:,:,:,:) 
!
integer                    :: lentmp1
integer, save, allocatable :: lrtmp1(:), lctmp1(:)
real*8,  save, allocatable :: lvtmp1(:)
!
end module auxmod
