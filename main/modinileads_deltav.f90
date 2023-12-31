subroutine modinileads_deltav
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, itier
integer*8 :: lni, lnj
real*8 :: dtmp1
complex*16 :: ctmp1
!
rho0 = rho
!
do itier=2,ntier
  do lni=nfirst(itier),nlast(itier) 
    ctmp1 = czero
    do ni=1,itier-1
       nj = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
       ctmp1 = ctmp1 + dcmplx(0.d0,  dipm(nj) * deltav_amps(ilead(nj), mspin(nj))) /       &
                       dcmplx(2.d0, -dipm(nj) * deltav_amps(ilead(nj), mspin(nj)))
    end do
    rho0(1:nrho, 1:nrho, lni) = (cunity + ctmp1) * rho(1:nrho, 1:nrho, lni)
  end do
end do
!
cmtmp1(1:nrho, 1:nrho) = rho0(1:nrho, 1:nrho, 1) - rho(1:nrho, 1:nrho, 1) ! rho(0^+) - rho(0) should be zero
dtmp1 = 0.d0
do nj=1,nrho
  do ni=1,nrho
    dtmp1 = max(dtmp1, cdabs(cmtmp1(ni,nj)))
  end do
end do
!
write(6,*)
write(6,*)' max(rho(0^+) - rho(0)) = ', dtmp1
call flush(6)
!
rho = rho0
!
end subroutine modinileads_deltav

