subroutine h2lvec(dimh, hmat, lhmat, lvec)
implicit none
!
! dimh-by-dimh Hilbert space matrix to dimh*dimh Liouville space vector
!
integer, intent(in) :: dimh, lhmat
real*8,  intent(in) :: hmat(lhmat,*)
real*8,  intent(out) :: lvec(*)
integer :: ni, nj, nk
!
nk = 0
do nj=1,dimh
  do ni=1,dimh
    nk = nk + 1
    lvec(nk) = hmat(ni,nj)
  end do
end do
end subroutine h2lvec
!
subroutine h2lmln(dimh, hmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), l(left), n(normal)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from left) to
! dimh**2-by-dimh**2 Liouville space matrix LA
!
! Y = A * X  => LY = LA * LX
!
integer, intent(in) :: dimh, lhmat, llmat
real*8,  intent(in) :: hmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
!
do nj=1,dimh
  do ni=1,dimh
    do nm=1,dimh
      n1 = (nj - 1) * dimh + ni
      n2 = (nj - 1) * dimh + nm
      lmat(n1,n2) = hmat(ni,nm)
    end do
  end do
end do
end subroutine h2lmln
!
subroutine h2lmrn(dimh, hmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), r(right), n(normal)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from right) to
! dimh**2-by-dimh**2 Liouville space matrix LA
!
! Y = X * A  => LY = LA * LX
!
integer, intent(in) :: dimh, lhmat, llmat
real*8,  intent(in) :: hmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
!
do nj=1,dimh
  do ni=1,dimh
    do nn=1,dimh
      n1 = (nj - 1) * dimh + ni
      n2 = (nn - 1) * dimh + ni
      lmat(n1,n2) = hmat(nn,nj)
    end do
  end do
end do
end subroutine h2lmrn
!
subroutine h2lmlt(dimh, hmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), l(left), t(transpose)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from left) to
! dimh**2-by-dimh**2 Liouville space matrix LA
!
! Y = A * X^T  => LY = LA * LX
!
integer, intent(in) :: dimh, lhmat, llmat
real*8,  intent(in) :: hmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
!
do nj=1,dimh
  do ni=1,dimh
    do nn=1,dimh
      n1 = (nj - 1) * dimh + ni
      n2 = (nn - 1) * dimh + nj
      lmat(n1,n2) = hmat(ni,nn)
    end do
  end do
end do
end subroutine h2lmlt
!
subroutine h2lmrt(dimh, hmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), r(right), t(transpose)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from right) to
! dimh**2-by-dimh**2 Liouville space matrix LA
!
! Y = X^T * A => LY = LA * LX
!
integer, intent(in) :: dimh, lhmat, llmat
real*8,  intent(in) :: hmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
!
do nj=1,dimh 
  do ni=1,dimh
    do nm=1,dimh
      n1 = (nj - 1) * dimh + ni
      n2 = (ni - 1) * dimh + nm
      lmat(n1,n2) = hmat(nm,nj)
    end do
  end do
end do
end subroutine h2lmrt
!
subroutine modifylmat(dimh, lmat, llmat, iop)
implicit none
integer, intent(in) :: dimh, llmat, iop
real*8,  intent(inout) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, n1, n2, n3
!
dim1 = dimh**2
do ni=1,dimh
  do nj=ni+1,dimh
    n2 = (nj - 1) * dimh + ni
    n3 = (ni - 1) * dimh + nj
    do n1=1,dim1
      if (iop .eq. 0) then           ! col : real part
        lmat(n1, n3) = lmat(n1, n3) + lmat(n1, n2)
        lmat(n1, n2) = 0.d0
      else if (iop .eq. 1) then      ! col : imaginary part
        lmat(n1, n3) = lmat(n1, n3) - lmat(n1, n2)
        lmat(n1, n2) = 0.d0
      else                           ! row : removal
        lmat(n2, n1) = 0.d0
      end if
    end do
  end do
end do
end subroutine modifylmat
!
subroutine h2lmlrn(dimh, ahmat, bhmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), l(left), r(right), n(normal)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from left) and
! dimh-by-dimh Hilbert space matrix B (multiplied from right) to
! dimh**2-by-dimh**2 Liouville space matrix L_AB
!
! Y = A * X * B  => LY = L_AB * LX
!
integer, intent(in)  :: dimh, lhmat, llmat
real*8,  intent(in)  :: ahmat(lhmat,*), bhmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
do nj=1,dimh
  do ni=1,dimh
    n1 = (nj - 1) * dimh + ni
    do nn=1,dimh
      do nm=1,dimh
        n2 = (nn - 1) * dimh + nm
        lmat(n1,n2) = ahmat(ni,nm) * bhmat(nn,nj)
      end do
    end do
  end do
end do
end subroutine h2lmlrn
!
subroutine h2lmlrt(dimh, ahmat, bhmat, lhmat, lmat, llmat)
implicit none
!
! m(matrix), l(left), r(right), n(normal)
!
! dimh-by-dimh Hilbert space matrix A (multiplied from left) and
! dimh-by-dimh Hilbert space matrix B (multiplied from right) to
! dimh**2-by-dimh**2 Liouville space matrix L_AB
!
! Y = A * X^T * B  => LY = L_AB * LX
!
integer, intent(in)  :: dimh, lhmat, llmat
real*8,  intent(in)  :: ahmat(lhmat,*), bhmat(lhmat,*)
real*8,  intent(out) :: lmat(llmat,*)
integer :: dim1
integer :: ni, nj, nm, nn, n1, n2
!
dim1 = dimh**2
lmat(1:dim1, 1:dim1) = 0.d0
do nj=1,dimh
  do ni=1,dimh
    n1 = (nj - 1) * dimh + ni
    do nn=1,dimh
      do nm=1,dimh
        n2 = (nn - 1) * dimh + nm
        lmat(n1,n2) = ahmat(ni,nn) * bhmat(nm,nj)
      end do
    end do
  end do
end do
end subroutine h2lmlrt
!
