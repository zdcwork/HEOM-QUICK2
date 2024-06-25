subroutine zmat_zcoomm3(alpha, a, nnza, irowa, icola, b, nnzb, irowb, icolb,  &
                         beta, c, ndim, ldc)
implicit none
!
! alpha * A * B + beta * C ==> C
! A and B are all in sparse COO format, and C is a dense matrix
!
! IMPORTANT: 
!    For COO storage, nonzeros with smaller column index appear earlier in the array
!    For efficiency, do not check the following facts:
!        ldc should be large enough, and ndim should be appropriate 
!
! Note:
!    To maximize speed, remove flags such as lefta, transa, transb
!   (no operations on A, B, C except for matrix multiplication)
!
integer,     intent(in)    :: nnza, nnzb, irowa(*), icola(*), irowb(*), icolb(*)
integer,     intent(in)    :: ndim, ldc
complex*16,  intent(in)    :: a(*), b(*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
integer                    :: ni, nj
integer                    :: mi, mj, mk
complex*16                 :: ctmp1
!
if (nnza .lt. 0 .or. nnzb .lt. 0) then
    write(6,*)
    write(6,*)'zmat_zcoomm3: error! wrong nnz found ', nnza, nnzb
    stop
end if
!
c(1:ndim,1:ndim) = beta * c(1:ndim,1:ndim) 
if (nnza .eq. 0 .or. nnzb .eq. 0) return
!
do ni=1,nnzb
   mj = icolb(ni)
   mk = irowb(ni)
   ctmp1 = b(ni) * alpha
   loop1: do nj=1,nnza
!      if (icola(nj) .gt. mk) exit loop1
      if (icola(nj) .ne. mk) cycle
      mi = irowa(nj)
!      c(mi,mj) = c(mi,mj) + a(nj) * b(ni) * alpha
      c(mi,mj) = c(mi,mj) + a(nj) * ctmp1
   end do loop1
end do
!
return
end subroutine zmat_zcoomm3
