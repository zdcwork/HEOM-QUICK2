subroutine zlmat_zcoomv4(transb, ndim, alpha, nnza, irowa, icola, cvala,   &
                         b, nnzb, irowb, icolb, beta, c, ldc, cmat1, ldc1, nvrow, nvcol)
implicit none
!
! alpha * A * op(B) + beta * op(C) ==> C
! A is a Liouville matrix in sparse COO format,
! B is a Hilbert-space matrix (Liouville-space vector) in sparse COO format, and C is a full matrix
!
! IMPORTANT: 
!    For COO storage, nonzeros with smaller column index appear earlier in the array
!
! Note: 
!
character*1, intent(in)    :: transb
integer,     intent(in)    :: ndim, nnza, nnzb, ldc, irowb(*), icolb(*), irowa(*), icola(*), ldc1
integer,     intent(in)    :: nvrow(*), nvcol(*)
complex*16,  intent(in)    :: cvala(*), b(*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
complex*16,  intent(inout) :: cmat1(ldc1,*)
!
if (nnza .lt. 0 .or. ndim .lt. 0 .or. nnzb .lt. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zlmat_zcoomv4: error! wrong para found ', ndim, nnza, nnzb, ldc
    stop
end if
!
if (ndim .eq. 0) return
if (nnza .eq. 0 .or. nnzb .eq. 0) then
    c(1:ndim,1:ndim) = c(1:ndim,1:ndim) * beta
    return
end if
!
call zmat_coo2dns(nnzb, irowb, icolb, b, cmat1, ndim, ndim, ldc1)
!
call zlmat_zcoomv3(transb, ndim, alpha, nnza, irowa, icola, cvala, cmat1, ldc1, beta, c, ldc, nvrow, nvcol)
!
return
end subroutine zlmat_zcoomv4
