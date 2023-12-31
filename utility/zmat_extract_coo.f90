subroutine zmat_extract_coo(nnz, irowa, icola, cvala, a, m, n, lda)
implicit none
!
integer,    intent(in)  :: nnz, irowa(*), icola(*), m, n, lda
complex*16, intent(in)  :: a(lda,*)
complex*16, intent(out) :: cvala(*)
integer                 :: ni
!
if (nnz .eq. 0) return
if (nnz .lt. 0) then
   write(6,*)
   write(6,*)'zmat_extract_coo: error, negative nnz found ', nnz
   stop
end if
do ni=1,nnz
   cvala(ni) = a(irowa(ni), icola(ni))
end do
!
return
end subroutine zmat_extract_coo
