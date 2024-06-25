subroutine zamsmm(sidea, transa, transb, korbs, kspin, zalpha, b, ldb, zbeta, c, ldc)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent(in) :: sidea, transa, transb
integer, intent(in)    :: korbs, kspin, ldb, ldc
complex*16, intent(in) :: zalpha, zbeta, b(ldb,*)
complex*16, intent(inout) :: c(ldc,*)
!
integer :: ni, nj
logical :: lsame
external :: lsame
logical :: nota, notb, conja, conjb, lefta, righta
character*1 :: ctra
!
lefta  = lsame(sidea, 'L')
righta = lsame(sidea, 'R')
nota   = lsame(transa, 'N')
notb   = lsame(transb, 'N')
conja  = lsame(transa, 'C')
conjb  = lsame(transb, 'C')
!
if ( (.not. nota) .and. (.not. conja) ) then
  write(6,*)
  write(6,*)' error input for transa in zamsmm', transa
  stop
end if
!
if ( (.not. notb) .and. (.not. conjb) ) then
  write(6,*)
  write(6,*)' error input for transb in zamsmm', transb
  stop
end if
!
if ( (.not. lefta) .and. (.not. righta) ) then
  write(6,*)
  write(6,*)' error input for sidea in zamsmm', sidea
  stop
end if
!
if (lefta .and. nota .and. notb) then
  ctra  = 'N'
  denmbr(1:nrho, 1:nrho) =  dble(zalpha) *  dble(b(1:nrho, 1:nrho)) - &
                           dimag(zalpha) * dimag(b(1:nrho, 1:nrho))
  denmbi(1:nrho, 1:nrho) =  dble(zalpha) * dimag(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) *  dble(b(1:nrho, 1:nrho))
else if (lefta .and. nota .and. conjb) then
  ctra = 'N'
  do nj=1,nrho
    do ni=1,nrho
      denmbr(ni,nj) =  dble(zalpha) *  dble(b(nj,ni)) + dimag(zalpha) * dimag(b(nj,ni))     
      denmbi(ni,nj) = -dble(zalpha) * dimag(b(nj,ni)) + dimag(zalpha) *  dble(b(nj,ni))
    end do
  end do
else if (lefta .and. conja .and. notb) then
  ctra = 'C'
  denmbr(1:nrho, 1:nrho) =  dble(zalpha) *  dble(b(1:nrho, 1:nrho)) - &
                           dimag(zalpha) * dimag(b(1:nrho, 1:nrho))
  denmbi(1:nrho, 1:nrho) =  dble(zalpha) * dimag(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) *  dble(b(1:nrho, 1:nrho))
else if (lefta .and. conja .and. conjb) then
  ctra = 'C'
  do nj=1,nrho
    do ni=1,nrho
      denmbr(ni,nj) =  dble(zalpha) *  dble(b(nj,ni)) + dimag(zalpha) * dimag(b(nj,ni))     
      denmbi(ni,nj) = -dble(zalpha) * dimag(b(nj,ni)) + dimag(zalpha) *  dble(b(nj,ni))
    end do
  end do
else if (righta .and. nota .and. notb) then
  ctra = 'C'
  do nj=1,nrho
    do ni=1,nrho
      denmbr(ni,nj) = dble(zalpha) *  dble(b(nj,ni)) - dimag(zalpha) * dimag(b(nj,ni))     
      denmbi(ni,nj) = dble(zalpha) * dimag(b(nj,ni)) + dimag(zalpha) *  dble(b(nj,ni))
    end do
  end do
else if (righta .and. nota .and. conjb) then
  ctra = 'C'
  denmbr(1:nrho, 1:nrho) =  dble(zalpha) *  dble(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) * dimag(b(1:nrho, 1:nrho)) 
  denmbi(1:nrho, 1:nrho) = -dble(zalpha) * dimag(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) *  dble(b(1:nrho, 1:nrho))
else if (righta .and. conja .and. notb) then
  ctra = 'N'
  do nj=1,nrho
    do ni=1,nrho
      denmbr(ni,nj) = dble(zalpha) *  dble(b(nj,ni)) - dimag(zalpha) * dimag(b(nj,ni))     
      denmbi(ni,nj) = dble(zalpha) * dimag(b(nj,ni)) + dimag(zalpha) *  dble(b(nj,ni))
    end do
  end do
else if (righta .and. conja .and. conjb) then
  ctra = 'N'
  denmbr(1:nrho, 1:nrho) =  dble(zalpha) *  dble(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) * dimag(b(1:nrho, 1:nrho)) 
  denmbi(1:nrho, 1:nrho) = -dble(zalpha) * dimag(b(1:nrho, 1:nrho)) + &
                           dimag(zalpha) *  dble(b(1:nrho, 1:nrho))
else
  write(6,*)
  write(6,*)' error! unsolved situation encountered in zamsmm '
  write(6,*)lefta, righta, nota, conja, notb, conjb
  stop
end if
!
call mkl_dcoomm(ctra, nrho, nrho, nrho, 1.d0, 'gogogo', valams(1,korbs,kspin), rowams(1,korbs,kspin), &
                colams(1,korbs,kspin), nams, denmbr, nrho, 0.d0, denmcr, nrho)
call mkl_dcoomm(ctra, nrho, nrho, nrho, 1.d0, 'gogogo', valams(1,korbs,kspin), rowams(1,korbs,kspin), &
                colams(1,korbs,kspin), nams, denmbi, nrho, 0.d0, denmci, nrho)
!
if (lefta) then
  c(1:nrho, 1:nrho) = dcmplx(denmcr(1:nrho, 1:nrho), denmci(1:nrho, 1:nrho)) + zbeta * c(1:nrho, 1:nrho)
else
  do nj=1,nrho
    do ni=1,nrho
      c(ni,nj) = zbeta * c(ni,nj) + dcmplx(denmcr(nj,ni), denmci(nj,ni))
    end do
  end do
end if
!
end subroutine zamsmm
