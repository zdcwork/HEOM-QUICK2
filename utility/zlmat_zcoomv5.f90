subroutine zlmat_zcoomv5(transb, ndim, alpha, nnza, irowa, icola, cvala,   &
                         b, nnzb, irowb, icolb, beta, c, ldc,              &
                         nlvrow, nlvcol, nnrc2v, ldc1, ncrc2v, ldc2, inda, nnca)
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
integer,     intent(in)    :: ndim, nnza, nnzb, ldc, ldc1, ldc2, irowb(*), icolb(*), irowa(*), icola(*)
integer,     intent(in)    :: inda(*), nnca(*)
integer,     intent(in)    :: nnrc2v(ldc1,*), ncrc2v(ldc2,*), nlvrow(*), nlvcol(*)
complex*16,  intent(in)    :: cvala(*), b(*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj, nk 
integer                    :: i, j, ii, jj
complex*16                 :: ctmp1
real*8, parameter          :: dsmall = 1.d-32
!
if (nnza .lt. 0 .or. ndim .lt. 0 .or. nnzb .lt. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zlmat_zcoomv5: error! wrong para found ', ndim, nnza, nnzb, ldc
    stop
end if
!
if (ndim .eq. 0) return
c(1:ndim,1:ndim) = c(1:ndim,1:ndim) * beta
!
if (nnza .eq. 0 .or. nnzb .eq. 0) return
if (cdabs(alpha) .lt. dsmall) return
!
if (lsame(transb, 'N')) then       ! no need to transpose
    do nk=1,nnzb
       j = nnrc2v(irowb(nk), icolb(nk))
       if (nnca(j) .eq. 0) cycle
       ctmp1 = b(nk) * alpha
       nj = inda(j)
       do ni=1,nnca(j)
          ii = nlvrow(irowa(nj))
          jj = nlvcol(irowa(nj))
!          if (j .ne. icola(nj)) then
!              write(6,*)'zlmat_zcoomv5: error! j, icola(nj) ', j, icola(nj)
!              write(6,*)transb
!              stop
!          end if
          c(ii,jj) = c(ii,jj) + ctmp1 * cvala(nj)
          nj = nj + 1
       end do
    end do
else if (lsame(transb, 'C')) then
    do nk=1,nnzb
       j = ncrc2v(irowb(nk), icolb(nk))
       if (nnca(j) .eq. 0) cycle
       ctmp1 = dconjg(b(nk)) * alpha
       nj = inda(j)
       do ni=1,nnca(j)
          ii = nlvrow(irowa(nj))
          jj = nlvcol(irowa(nj))
!          if (j .ne. icola(nj)) then
!              write(6,*)'zlmat_zcoomv5: error! j, icola(nj) ', j, icola(nj)
!              write(6,*)transb
!              stop
!          end if
          c(ii,jj) = c(ii,jj) + ctmp1 * cvala(nj)
          nj = nj + 1
       end do
    end do
else
    write(6,*)
    write(6,*)'zlmat_zcoomv5: error! unknown transb ', transb
    stop
end if
!
return
end subroutine zlmat_zcoomv5
