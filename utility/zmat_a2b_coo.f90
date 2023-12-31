subroutine zmat_a2b_coo(transa, nnz, irowa, icola, zvala, irowb, icolb, zvalb)
implicit none
!
! B = op(A), op(A) can be A, A^dag, A^t, or A^*
! copy sparse matrix op(A) in COO format to matrix B in COO format
!
character*1, intent(in)  :: transa
integer,     intent(in)  :: nnz, irowa(*), icola(*)
integer,     intent(out) :: irowb(*), icolb(*)
complex*16,  intent(in)  :: zvala(*) 
complex*16,  intent(out) :: zvalb(*) 
logical                  :: lsame
external                 :: lsame
!
if (nnz .eq. 0) return
!
if (nnz .lt. 0) then
   write(6,*)
   write(6,*)'zmat_a2b_coo: error! nnz < 0 found ', nnz
   stop
end if
!
if ( lsame(transa, 'N') ) then
   irowb(1:nnz) = irowa(1:nnz)    
   icolb(1:nnz) = icola(1:nnz)    
   zvalb(1:nnz) = zvala(1:nnz)
!
else if ( lsame(transa, 'T') ) then
   irowb(1:nnz) = icola(1:nnz)
   icolb(1:nnz) = irowa(1:nnz)
   zvalb(1:nnz) = zvala(1:nnz)
!
else if ( lsame(transa, 'C') ) then
   irowb(1:nnz) = icola(1:nnz)
   icolb(1:nnz) = irowa(1:nnz)
   zvalb(1:nnz) = dconjg(zvala(1:nnz))
!
else if ( lsame(transa, 'R') ) then
   irowb(1:nnz) = irowa(1:nnz)    
   icolb(1:nnz) = icola(1:nnz)    
   zvalb(1:nnz) = dconjg(zvala(1:nnz))
!
else 
   write(6,*)
   write(6,*)'zmat_a2b_coo: error! unknown input transa ', transa
   stop
end if
!
return
end subroutine zmat_a2b_coo
