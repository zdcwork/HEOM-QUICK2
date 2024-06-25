subroutine zmat_accu_coo2dns_pair(transa, nnz, irowa, icola, cvala, coefb, bmat, ldb, &
                                  coefc, cmat, ldc)
implicit none
!
! coefb * op(A) is accumulated to dense matrix B, A is in COO sparse form,
! and
! coefc * op(A) is accumulated to dense matrix C  
! 
character*1, intent(in)    :: transa
integer,     intent(in)    :: nnz, irowa(*), icola(*), ldb, ldc
complex*16,  intent(in)    :: coefb, coefc, cvala(*)
complex*16,  intent(inout) :: bmat(ldb,*), cmat(ldc,*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni
!
if (nnz .eq. 0) return
if (nnz .lt. 0) then
    write(6,*)
    write(6,*)'zmat_accu_coo2dns_pair: error! negative nnz found ', nnz
    stop
end if
!
if ( lsame(transa, 'N') ) then
   do ni=1,nnz
      bmat(irowa(ni),icola(ni)) = bmat(irowa(ni),icola(ni)) + coefb * cvala(ni)
      cmat(irowa(ni),icola(ni)) = cmat(irowa(ni),icola(ni)) + coefc * cvala(ni)
   end do
else if ( lsame(transa, 'R') ) then
   do ni=1,nnz
      bmat(irowa(ni),icola(ni)) = bmat(irowa(ni),icola(ni)) + coefb * dconjg(cvala(ni))
      cmat(irowa(ni),icola(ni)) = cmat(irowa(ni),icola(ni)) + coefc * dconjg(cvala(ni))
   end do
else if ( lsame(transa, 'T') ) then
   do ni=1,nnz
      bmat(icola(ni),irowa(ni)) = bmat(icola(ni),irowa(ni)) + coefb * cvala(ni)
      cmat(icola(ni),irowa(ni)) = cmat(icola(ni),irowa(ni)) + coefc * cvala(ni)
   end do
else if ( lsame(transa, 'C') ) then
   do ni=1,nnz
      bmat(icola(ni),irowa(ni)) = bmat(icola(ni),irowa(ni)) + coefb * dconjg(cvala(ni))
      cmat(icola(ni),irowa(ni)) = cmat(icola(ni),irowa(ni)) + coefc * dconjg(cvala(ni))
   end do
else
   write(6,*)'zmat_accu_coo2dns_pair: error! unknown transa ', transa
   stop
end if
!
return
end subroutine zmat_accu_coo2dns_pair
