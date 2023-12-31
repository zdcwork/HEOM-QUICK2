subroutine zlmat_zcoomv1(transb, ndim, alpha, a, lda, b, nnzb, irowb, icolb, beta, c, ldc, cval2)
implicit none
!
! alpha * A * op(B) + beta * op(C) ==> C
! A is a full Liouville matrix, B is a Hilbert-space matrix (Liouville-space vector) in sparse COO format, and C is a full matrix
!
! IMPORTANT: 
!    For COO storage, nonzeros with smaller column index appear earlier in the array
!
! Note: 
!
character*1, intent(in)    :: transb
integer,     intent(in)    :: ndim, lda, nnzb, ldc, irowb(*), icolb(*)
complex*16,  intent(in)    :: a(lda,*), b(*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
complex*16,  intent(inout) :: cval2(*)
logical                    :: lcb
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj, nk, ndim2
integer                    :: ij, mn
complex*16                 :: ctmp1
real*8, parameter          :: dsmall = 1.d-32
!
if (lda .le. 0 .or. ndim .lt. 0 .or. nnzb .lt. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zlmat_zcoomv1: error! wrong para found ', ndim, lda, nnzb, ldc
    stop
end if
!
c(1:ndim,1:ndim) = c(1:ndim,1:ndim) * beta
if (ndim .eq. 0 .or. nnzb .eq. 0) return
if (cdabs(alpha) .lt. dsmall) return
ndim2 = ndim**2
!
if (lsame(transb, 'N') .or. lsame(transb, 'R')) then       ! no need to transpose
    lcb = .false.  
else if (lsame(transb, 'C') .or. lsame(transb, 'T')) then  ! need to transpose
    lcb = .true.
else
    write(6,*)
    write(6,*)'zlmat_zcoomv1: error! unknown transb ', transb
    stop
end if
!
if (lsame(transb, 'N') .or. lsame(transb, 'T')) then       ! no conjugation
    cval2(1:nnzb) = b(1:nnzb)
else                                                       ! conjugation
    cval2(1:nnzb) = dconjg(b(1:nnzb))
end if
!
if (lcb) then
    do nk=1,nnzb
       mn = (irowb(nk) - 1) * ndim + icolb(nk)
       ctmp1 = cval2(nk) * alpha
       ij = 0 
       do nj=1,ndim
          do ni=1,ndim
             ij = ij + 1
             c(ni,nj) = c(ni,nj) + a(ij,mn) * ctmp1
          end do
       end do 
    end do
else
    do nk=1,nnzb
       mn = (icolb(nk) - 1) * ndim + irowb(nk)
       ctmp1 = cval2(nk) * alpha
       ij = 0 
       do nj=1,ndim
          do ni=1,ndim
             ij = ij + 1
             c(ni,nj) = c(ni,nj) + a(ij,mn) * ctmp1
          end do
       end do 
    end do
end if
!
return
end subroutine zlmat_zcoomv1
