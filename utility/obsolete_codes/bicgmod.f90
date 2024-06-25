module bicgmod
implicit none
!
! bicg memory
!
real*8,       allocatable :: rhsrbicg(:), rhsibicg(:)
real*8,       allocatable :: sk(:), rk(:), skp(:), rkp(:), dk(:), fk(:), tmp1bicg(:)
!
real*8,       allocatable :: zetatmpr(:), phitmpr(:), psitmpr(:)
real*8,       allocatable :: zetatmpi(:), phitmpi(:), psitmpi(:)
real*8,       allocatable :: zetarhsr(:), phirhsr(:), psirhsr(:)
real*8,       allocatable :: zetarhsi(:), phirhsi(:), psirhsi(:)
!
end module bicgmod
