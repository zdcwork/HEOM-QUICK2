subroutine zlmat_zcoomv3(transb, ndim, alpha, nnza, irowa, icola, cvala, b, ldb, beta, c, ldc, nvrow, nvcol)
implicit none
!
! alpha * A * op(B) + beta * op(C) ==> C
! A is a sparse Louville matrix in the COO format, 
! B is a Hilbert-space matrix (Liouville-space vector), and C is a full matrix
!
! IMPORTANT: 
!
! Note: 
!
character*1, intent(in)    :: transb
integer,     intent(in)    :: ndim, nnza, ldb, ldc, irowa(*), icola(*), nvrow(*), nvcol(*)
complex*16,  intent(in)    :: cvala(*), b(ldb,*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
logical                    :: lcb, lconjg
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj, nk 
integer                    :: m, n, i, j
complex*16                 :: ctmp1, ctmp2
real*8, parameter          :: dsmall = 1.d-32
!
if (nnza .lt. 0 .or. ndim .lt. 0 .or. ldb .le. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zlmat_zcoomv3: error! wrong para found ', ndim, nnza, ldb, ldc
    stop
end if
!
c(1:ndim,1:ndim) = c(1:ndim,1:ndim) * beta
if (ndim .eq. 0 .or. nnza .eq. 0) return
if (cdabs(alpha) .lt. dsmall) return
!
if (lsame(transb, 'N') .or. lsame(transb, 'R')) then       ! no need to transpose
    lcb = .false.  
else if (lsame(transb, 'C') .or. lsame(transb, 'T')) then  ! need to transpose
    lcb = .true.
else
    write(6,*)
    write(6,*)'zlmat_zcoomv3: error! unknown transb ', transb
    stop
end if
!
if (lsame(transb, 'N') .or. lsame(transb, 'T')) then       ! no conjugation
    lconjg = .false.
else                                                       ! conjugation
    lconjg = .true.
end if
!
!do nk=1,nnza
!   call jmod(icola(nk), ndim, m, n)
!   if (lconjg .and. lcb) then
!       ctmp1 = dconjg(b(n,m))
!   else if ((.not. lconjg) .and. lcb) then
!       ctmp1 = b(n,m)
!   else if ((.not. lconjg) .and. (.not. lcb)) then
!       ctmp1 = b(m,n)
!   else 
!       ctmp1 = dconjg(b(m,n))
!   end if
!   call jmod(irowa(nk), ndim, i, j)
!   c(i,j) = c(i,j) + cvala(nk) * ctmp1 * alpha
!end do
!
if (lcb .and. lconjg) then
    do nk=1,nnza
       !call jmod(icola(nk), ndim, m, n)
       !call jmod(irowa(nk), ndim, i, j)
       m = nvrow(icola(nk))
       n = nvcol(icola(nk))
       i = nvrow(irowa(nk))
       j = nvcol(irowa(nk))
       ctmp1 = dconjg(b(n,m))
       c(i,j) = c(i,j) + cvala(nk) * ctmp1 * alpha
    end do
else if ((.not. lconjg) .and. lcb) then
    do nk=1,nnza
       !call jmod(icola(nk), ndim, m, n)
       !call jmod(irowa(nk), ndim, i, j)
       m = nvrow(icola(nk))
       n = nvcol(icola(nk))
       i = nvrow(irowa(nk))
       j = nvcol(irowa(nk))
       ctmp1 = b(n,m)
       c(i,j) = c(i,j) + cvala(nk) * ctmp1 * alpha
    end do
else if ((.not. lconjg) .and. (.not. lcb)) then
    do nk=1,nnza
       !call jmod(icola(nk), ndim, m, n)
       !call jmod(irowa(nk), ndim, i, j)
       m = nvrow(icola(nk))
       n = nvcol(icola(nk))
       i = nvrow(irowa(nk))
       j = nvcol(irowa(nk))
       ctmp1 = b(m,n)
       c(i,j) = c(i,j) + cvala(nk) * ctmp1 * alpha
    end do
else
    do nk=1,nnza
       !call jmod(icola(nk), ndim, m, n)
       !call jmod(irowa(nk), ndim, i, j)
       m = nvrow(icola(nk))
       n = nvcol(icola(nk))
       i = nvrow(irowa(nk))
       j = nvcol(irowa(nk))
       ctmp1 = dconjg(b(m,n))
       c(i,j) = c(i,j) + cvala(nk) * ctmp1 * alpha
    end do
end if
!
return
end subroutine zlmat_zcoomv3
