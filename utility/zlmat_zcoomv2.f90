subroutine zlmat_zcoomv2(transb, ndim, alpha, a, lda, b, ldb, beta, c, ldc)
implicit none
!
! alpha * A * op(B) + beta * op(C) ==> C
! A is a full Liouville matrix, B is a Hilbert-space matrix (Liouville-space vector), and C is a full matrix
!
! IMPORTANT: 
!
! Note: 
!
character*1, intent(in)    :: transb
integer,     intent(in)    :: ndim, lda, ldb, ldc
complex*16,  intent(in)    :: a(lda,*), b(ldb,*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*)
logical                    :: lcb, lconjg
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj, nk, ndim2
integer                    :: ij, mn, mi, mj
complex*16                 :: ctmp1
real*8, parameter          :: dsmall = 1.d-32
!
if (lda .le. 0 .or. ndim .lt. 0 .or. ldb .le. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zlmat_zcoomv2: error! wrong para found ', ndim, lda, ldb, ldc
    stop
end if
!
c(1:ndim,1:ndim) = c(1:ndim,1:ndim) * beta
if (ndim .eq. 0) return
if (cdabs(alpha) .lt. dsmall) return
ndim2 = ndim**2
!
if (lsame(transb, 'N') .or. lsame(transb, 'R')) then       ! no need to transpose
    lcb = .false.  
else if (lsame(transb, 'C') .or. lsame(transb, 'T')) then  ! need to transpose
    lcb = .true.
else
    write(6,*)
    write(6,*)'zlmat_zcoomv2: error! unknown transb ', transb
    stop
end if
!
if (lsame(transb, 'N') .or. lsame(transb, 'T')) then       ! no conjugation
    lconjg = .false.
else                                                       ! conjugation
    lconjg = .true.
end if
!
if (lcb) then
    do mj=1,ndim
       do mi=1,ndim
          mn = (mi - 1) * ndim + mj
          if (lconjg) then
              ctmp1 = dconjg(b(mj,mi)) * alpha
          else
              ctmp1 = b(mj,mi) * alpha
          end if
          ij = 0 
          do nj=1,ndim
             do ni=1,ndim
                ij = ij + 1
                c(ni,nj) = c(ni,nj) + a(ij,mn) * ctmp1
             end do
          end do 
       end do
    end do
else
    do mj=1,ndim
       do mi=1,ndim
          mn = (mj - 1) * ndim + mi
          if (lconjg) then
              ctmp1 = dconjg(b(mi,mj)) * alpha
          else
              ctmp1 = b(mi,mj) * alpha
          end if
          ij = 0 
          do nj=1,ndim
             do ni=1,ndim
                ij = ij + 1
                c(ni,nj) = c(ni,nj) + a(ij,mn) * ctmp1
             end do
          end do 
       end do
    end do
end if
!
return
end subroutine zlmat_zcoomv2
