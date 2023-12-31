subroutine zmat_zcoomm(lefta, transa, transb, n, alpha, cvala, irowa, icola, nnz, &
                       b, ldb, beta, c, ldc,                                      &
                       cvaltmp, irowtmp, icoltmp, btmp, ldbtmp, ctmp, ldctmp)
implicit none
!
! An interface subroutine to realize C = alpha * op(A) * op(B) + beta * C
! op(X) = X, X^t (transpose of X), or X^dag (complex conjugate of X)
! A is in COO sparse format, and B, C are dense matrices
! All matrices are square matrices, and have the same dimension of (n,n)
!
character*1, intent(in)    :: lefta, transa, transb
integer,     intent(in)    :: n, nnz, ldb, ldc, irowa(*), icola(*)
integer,     intent(in)    :: ldbtmp, ldctmp
complex*16,  intent(in)    :: alpha, beta, cvala(*), b(ldb,*)
complex*16,  intent(inout) :: c(ldc,*), cvaltmp(*), btmp(ldbtmp,*), ctmp(ldctmp,*) 
integer,     intent(inout) :: irowtmp(*), icoltmp(*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj
character*1                :: transbt
character                  :: matdescra(6)
!
if (nnz .eq. 0) then
   c(1:n,1:n) = beta * c(1:n,1:n)
   return
end if
!
if (nnz .lt. 0 .or. nnz .gt. n*n) then
   write(6,*)
   write(6,*)'zmat_zcoomm: error! wrong nnz found ', nnz, n
   stop
end if
!
matdescra(1) = 'g'
matdescra(4) = 'f'
!
if ( lsame(lefta, 'L') ) then
   call zmat_a2b_coo('n', nnz, irowa, icola, cvala, irowtmp, icoltmp, cvaltmp)
   call zmat_a2b_dns(transb, n, b, ldb, btmp, ldbtmp)
   call mkl_zcoomm(transa, n, n, n, alpha, matdescra, cvaltmp, irowtmp, icoltmp, &
                   nnz, btmp, ldbtmp, beta, c, ldc)
!
else if ( lsame(lefta, 'R') ) then
! transpose A and B
   if ( lsame(transb, 'N') ) then
      transbt = 'T'
   else if ( lsame(transb, 'T') ) then
      transbt = 'N'
   else if ( lsame(transb, 'C') ) then
      transbt = 'R'
   else if ( lsame(transb, 'R') ) then
      transbt = 'C'
   else 
      write(6,*)
      write(6,*)'zmat_zcoomm: error! unknown transb ', transb
      stop
   end if
!
   call zmat_a2b_coo('t', nnz, irowa, icola, cvala, irowtmp, icoltmp, cvaltmp)
   call zmat_a2b_dns(transbt, n, b, ldb, btmp, ldbtmp)
   call mkl_zcoomm(transa, n, n, n, alpha, matdescra, cvaltmp, irowtmp, icoltmp, &
                   nnz, btmp, ldbtmp, dcmplx(0.d0, 0.d0), ctmp, ldctmp)
   do nj=1,n
      do ni=1,n
         c(ni,nj) = ctmp(nj,ni) + beta * c(ni,nj)
      end do
   end do
!
else
   write(6,*)
   write(6,*)'zmat_zcoomm: error! unknown lefta ', lefta
   stop
end if
!
return
end subroutine zmat_zcoomm
