subroutine l2hvec(dimh, hmat, lhmat, lvec)
implicit none
!
! dimh-by-dimh Hilbert space matrix from dimh*dimh Liouville space vector
!
integer, intent(in) :: dimh, lhmat
real*8,  intent(in) :: lvec(*)
real*8,  intent(out) :: hmat(lhmat,*)
integer :: ni, nj, nk
!
do nj=1,dimh
  do ni=1,dimh
    hmat(ni,nj) = lvec( (nj-1)*dimh+ni )
  end do
end do
end subroutine l2hvec
!
