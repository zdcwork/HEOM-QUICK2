subroutine zmat_coo2dns(nnz, irowa, icola, cvala, a, m, n, lda)
implicit none
!
! convert a sparse matrix in COO format to dense matrix and output
!
integer,    intent(in)  :: nnz, m, n, lda, irowa(*), icola(*) 
complex*16, intent(in)  :: cvala(*)
complex*16, intent(out) :: a(lda,*)
integer                 :: ni
!
! allow `nnz = 0' case : return without action 
!
if (nnz .lt. 0 .or. nnz .gt. m*n) then
    write(6,*)
    write(6,*)'zmat_coo2dns: error size nnz ', nnz, m, n
    stop
end if
a(1:m,1:n) = dcmplx(0.d0, 0.d0)
do ni=1,nnz
   a(irowa(ni),icola(ni)) = cvala(ni)
end do
!
return
end subroutine zmat_coo2dns
!
subroutine zmat_dns2coo(nnz, irowa, icola, cvala, a, m, n, lda, dsmall_spa)
implicit none
!
! convert a dense matrix to a sparse matrix in COO format and output
! criterion for a significant value is : dsmall_spa 
!
integer,    intent(in)  :: m, n, lda
integer,    intent(out) :: nnz, irowa(*), icola(*) 
real*8,     intent(in)  :: dsmall_spa
complex*16, intent(out) :: cvala(*)
complex*16, intent(in)  :: a(lda,*)
integer                 :: ni, nj
!real*8, parameter       :: dpico = 1.d-12, dsmall = 1.d-32
!
nnz = 0
do nj=1,n ! column first
   do ni=1,m
!      if ( cdabs(a(ni,nj)) .ge. dsmall ) then
      if ( cdabs(a(ni,nj)) .ge. dsmall_spa ) then
         nnz = nnz + 1
         irowa(nnz) = ni
         icola(nnz) = nj
         cvala(nnz) = a(ni,nj)
      end if
   end do 
end do
!
return
end subroutine zmat_dns2coo
!
