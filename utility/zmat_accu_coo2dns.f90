subroutine zmat_accu_coo2dns(transa, alpha, nnz, irowa, icola, cvala, ndim, bmat, ldb)
implicit none
!
! alpha * op(A) is accumulated to dense matrix B, A is in COO sparse form
! 
character*1, intent(in)    :: transa
integer,     intent(in)    :: nnz, irowa(*), icola(*), ndim, ldb
complex*16,  intent(in)    :: alpha, cvala(*)
complex*16,  intent(inout) :: bmat(ldb,*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni
!
if (nnz .eq. 0) return
if (nnz .lt. 0) then
    write(6,*)
    write(6,*)'zmat_accu_coo2dns: error! negative nnz found ', nnz
    stop
end if
!
if ( lsame(transa, 'N') ) then
   do ni=1,nnz
      bmat(irowa(ni),icola(ni)) = bmat(irowa(ni),icola(ni)) + alpha * cvala(ni)
   end do
else if ( lsame(transa, 'R') ) then
   do ni=1,nnz
      bmat(irowa(ni),icola(ni)) = bmat(irowa(ni),icola(ni)) + alpha * dconjg(cvala(ni))
   end do
else if ( lsame(transa, 'T') ) then
   do ni=1,nnz
      bmat(icola(ni),irowa(ni)) = bmat(icola(ni),irowa(ni)) + alpha * cvala(ni)
   end do
else if ( lsame(transa, 'C') ) then
   do ni=1,nnz
      bmat(icola(ni),irowa(ni)) = bmat(icola(ni),irowa(ni)) + alpha * dconjg(cvala(ni))
   end do
else
   write(6,*)'zmat_accu_coo2dns: error! unknown transa ', transa
   stop
end if
!
return
end subroutine zmat_accu_coo2dns
