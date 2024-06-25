subroutine zmat_herm_coo(uplo, nnz, a, irowa, icola, zmtmp, n, ldtmp)
implicit none
!
! use upper or lower half matrix to build the full hermitian matrix 
! the transpose counterpart of (i,j) nonzero element (j,i) should already
!     be included in the COO topology 
!
character*1, intent(in)    :: uplo
integer,     intent(in)    :: nnz, n, ldtmp, irowa(*), icola(*) 
complex*16,  intent(inout) :: a(*), zmtmp(ldtmp,*)
integer                    :: ni, nj
real*8                     :: dtmp1
logical                    :: lsame
external                   :: lsame
!
if ( .not. lsame(uplo, 'U') .and. .not. lsame(uplo, 'L') ) then
   write(6,*)
   write(6,*)'zmat_herm_coo: error, unknown uplo ', uplo
   stop
end if
!
zmtmp(1:n, 1:n) = dcmplx(0.d0, 0.d0)
!
if ( lsame(uplo, 'U') ) then
   do ni=1,nnz
      if ( icola(ni) .gt. irowa(ni) ) then
          zmtmp(irowa(ni), icola(ni)) = a(ni)
          zmtmp(icola(ni), irowa(ni)) = dconjg(a(ni)) 
      end if
      if ( icola(ni) .eq. irowa(ni) ) then
          dtmp1 = dble( a(ni) )
          zmtmp(irowa(ni), icola(ni)) = dcmplx(dtmp1, 0.d0)
      end if
   end do 
else 
   do ni=1,nnz
      if ( icola(ni) .lt. irowa(ni) ) then
          zmtmp(irowa(ni), icola(ni)) = a(ni)
          zmtmp(icola(ni), irowa(ni)) = dconjg(a(ni)) 
      end if
      if ( icola(ni) .eq. irowa(ni) ) then
          dtmp1 = dble( a(ni) )
          zmtmp(irowa(ni), icola(ni)) = dcmplx(dtmp1, 0.d0)
      end if
   end do 
end if
!
call zmat_extract_coo(nnz, irowa, icola, a, zmtmp, n, n, ldtmp)
!
return
end subroutine zmat_herm_coo
