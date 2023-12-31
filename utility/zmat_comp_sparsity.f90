subroutine zmat_comp_sparsity(a, ma, na, lda, b, mb, nb, ldb, isame, icomp, nnz_a, nnz_b, dsmall_spa)
implicit none
!
! isame = 1 : a and b have same sparsity topology
! icomp = 1 : a is a subset of b (every nonzero element in a has its correspondence in b)
!             (opposite relation not necessarily true)
!
integer,   intent(in)  :: ma, na, lda, mb, nb, ldb
integer,   intent(out) :: isame, icomp, nnz_a, nnz_b
complex*16, intent(in) :: a(lda,*), b(ldb,*)
real*8, intent(in)     :: dsmall_spa
integer                :: ni, nj
!real*8, parameter      :: dpico = 1.d-12, dsmall = 1.d-32
!
if (ma .le. 0 .or. na .le. 0 .or. mb .le. 0 .or. nb .le. 0) then
   write(6,*)
   write(6,*)'zmat_comp_sparsity: wrong size ', ma, na, mb, nb
   stop
end if
if (ma .ne. mb .or. na .ne. nb) then
   write(6,*)
   write(6,*)'zmat_comp_sparsity: wrong dimension ', ma, na, mb, nb
   stop
end if
!
nnz_a = 0
nnz_b = 0
do nj=1,na
   do ni=1,ma
      if ( cdabs(a(ni,nj)) .ge. dsmall_spa ) nnz_a = nnz_a + 1
      if ( cdabs(b(ni,nj)) .ge. dsmall_spa ) nnz_b = nnz_b + 1
   end do
end do
isame = 1
if (nnz_a .ne. nnz_b) isame = 0
icomp = 1
check1: do nj=1,na
   do ni=1,ma
      if ( cdabs(a(ni,nj)) .ge. dsmall_spa ) then
         if ( cdabs(b(ni,nj)) .lt. dsmall_spa ) then
            isame = 0
            icomp = 0
            exit check1
         end if
      end if
   end do
end do check1
!
return
end subroutine zmat_comp_sparsity

