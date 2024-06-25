subroutine modinileads_spa
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, itier, nnz
integer*8 :: lni, lnj
real*8 :: dtmp1
complex*16 :: ctmp1
!
rho0_spa = rho_spa
!
do itier=2,ntier
  do lni=nfirst(itier),nlast(itier) 
    ctmp1 = czero
    do ni=1,itier-1
       nj = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
!      if (ilead(nj) .eq. 1) then
!        ctmp1 = ctmp1 + dcmplx(0.d0, dipm(nj) * ampL) / dcmplx(2, -dipm(nj) * ampL)
!      else if (ilead(nj) .eq. 2) then
!        ctmp1 = ctmp1 + dcmplx(0.d0, dipm(nj) * ampR) / dcmplx(2, -dipm(nj) * ampR)
!      else
!        write(6,*)
!        write(6,*)' error in modinileads 1', ilead(nj)
!        stop
!      end if
       ctmp1 = ctmp1 + dcmplx(0.d0,  dipm(nj) * amps(ilead(nj), mspin(nj))) /       &
                       dcmplx(2.d0, -dipm(nj) * amps(ilead(nj), mspin(nj)))
    end do
    nnz = nnz_spa(lni)
    lnj = ind_spa(lni)
    rho0_spa(lnj:lnj-1+nnz) = (cunity + ctmp1) * rho_spa(lnj:lnj-1+nnz)
  end do
end do
!
nnz = nnz_spa(1)
lnj = ind_spa(1)
dtmp1 = 0.d0
do ni=1,nnz
   dtmp1 = max(dtmp1, cdabs(rho0_spa(lnj-1+ni) - rho_spa(lnj-1+ni)))
end do
!
write(6,*)
write(6,*)' max(rho(0^+) - rho(0)) = ', dtmp1
call flush(6)
!
rho_spa = rho0_spa
!
end subroutine modinileads_spa

