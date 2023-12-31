subroutine chebyocc(rhoinp, occu, occd, sigmau, sigmad)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in)  :: rhoinp(nrho,*)
real*8,     intent(out) :: occu, occd
complex*16, intent(out) :: sigmau(*), sigmad(*)
integer :: ni, nj, nk, istat
!
occu = 0.d0
occd = 0.d0
!
! population = tr( c^\dag_{ms} c_{ms} \rho )
!
do nk=1,nspin
  do ni=1,norbs
    cmtmp1(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, ni, nk), 0.d0)
    cmtmp2(1:nrho, 1:nrho) = cmtmp1(1:nrho, 1:nrho)
    call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho,       &
               czero, cmtmp3, nrho)
    call zhemm('r', 'u', nrho, nrho, cunity, rhoinp(1,1), nrho, cmtmp3, nrho, czero, &
               cmtmp2, nrho)
    if (nk .eq. 1) then
      sigmau(ni) = czero
      do nj=1,nrho
       sigmau(ni) = sigmau(ni) + cmtmp2(nj,nj)
      end do
      occu = occu + dble(sigmau(ni))
    else 
      sigmad(ni) = czero
      do nj=1,nrho
       sigmad(ni) = sigmad(ni) + cmtmp2(nj,nj)
      end do
      occd = occd + dble(sigmad(ni))  
    end if
  end do
end do
!
end subroutine chebyocc
