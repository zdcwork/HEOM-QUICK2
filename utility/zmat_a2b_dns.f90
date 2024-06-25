subroutine zmat_a2b_dns(transa, n, a, lda, b, ldb)
implicit none
!
! copy dense square matrix op(A) to matrix B 
! matrix dimension : (n,n)
! for non-square matrix, use MKL subroutine mkl_zimatcopy
!
character*1, intent(in)  :: transa
integer,     intent(in)  :: n, lda, ldb
complex*16,  intent(in)  :: a(lda,*)
complex*16,  intent(out) :: b(ldb,*)
integer                  :: ni, nj
logical                  :: lsame
external                 :: lsame
!
if ( lsame(transa, 'N') ) then
   b(1:n,1:n) = a(1:n,1:n)    
!
else if ( lsame(transa, 'T') ) then   ! transpose
   do nj=1,n
      do ni=1,n
         b(ni,nj) = a(nj,ni)
      end do
   end do
!
else if ( lsame(transa, 'R') ) then   ! complex conjugate
   b(1:n,1:n) = dconjg( a(1:n,1:n) )
!
else if ( lsame(transa, 'C') ) then   ! tranpose conjugate
   do nj=1,n
      do ni=1,n
         b(ni,nj) = dconjg(a(nj,ni))
      end do
   end do
!  
else
   write(6,*)
   write(6,*)'zmat_a2b_dns: error! unknown input transa ', transa
   stop
end if
!
return
end subroutine zmat_a2b_dns
