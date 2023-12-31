subroutine zmat_add_coo(transa, nnz, alpha, cvala, beta, cvalb)
implicit none
!
! B = alpha * op(A) + beta * B
! op(A) can only be A or complex conjugate of A
! A and B are in COO format, and have same sparsity topology
!
character*1, intent(in)    :: transa
integer,     intent(in)    :: nnz
complex*16,  intent(in)    :: alpha, beta, cvala(*)
complex*16,  intent(inout) :: cvalb(*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni
!
if (nnz .eq. 0) return
if (nnz .lt. 0) then
   write(6,*)
   write(6,*)'zmat_add_coo: error! negative nnz found ', nnz
   stop
end if
!
if ( lsame(transa, 'N') ) then
   do ni=1,nnz
      cvalb(ni) = beta * cvalb(ni) + alpha * cvala(ni)
   end do
else if ( lsame(transa, 'R') ) then
   do ni=1,nnz
      cvalb(ni) = beta * cvalb(ni) + alpha * dconjg(cvala(ni))
   end do
else
   write(6,*)
   write(6,*)'zmat_add_coo: error! wrong transa ', transa
   stop
end if
!
return
end subroutine zmat_add_coo
