subroutine calchs
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni
real*8 :: edot
real*8 :: eup, edown, uu, ac
namelist / para1 / eup, edown, uu, ac, capa
namelist / para2 / edot
!
! single-site, spin-splitted
! 
! H_S = eup * c^\dag_up * c_up + edown * c^\dag_down * c_down +
!        uu * c^\dag_up * c_up * c^\dag_down * c_down
!
hs(1:nrho, 1:nrho) = czero
!
!--------------------- for (norbs=1, nspin=1) ----------------
if (nspin .eq. 1 .and. norbs .eq. 1) then
  edot = 0.d0
  rewind(5)
  read(5, para2, end=91)
  91 continue
  epara1 = edot
  write(6,*)
  write(6,*)' parameters for Hamiltonian '
  write(6,*)' edot = ', edot
  edot = edot / hbar
!   
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0)
  call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(edot, 0.d0), cmtmp1, nrho, cmtmp1, nrho, &
              czero, cmtmp2, nrho)
  hs = hs + cmtmp2
end if
!---------------------end of area for (norbs=1, nspin=1) ----------------
!
if (nspin .eq. 2 .and. norbs .eq. 1) then
  eup   = 0.d0
  edown = 0.d0
  uu    = 0.d0
  ac    = 0.d0
  capa  = 0.d0
  rewind(5)
  read(5, para1, end=101)
  101 continue
rewind(5)
read(5,para2, end=109)
109 continue
  epara1 = eup
  epara2 = edown
  epara3 = uu
  epara4 = ac
  epara6 = edot

  if (capa .le. -dnano) then
    write(6,*)
    write(6,*)' error! negative capacitance read from input ', capa
    stop
  end if
!
  write(6,*)
  write(6,*)' parameters for Hamiltonian '
!
  if (dabs(ac) .ge. dnano) write(6,112) ac
  if (capa .ge. dnano) then
    ec     = qe0 * 5.d-1 / capa
    epara5 = ec
    uu     = 5.d-1 * ec         ! Ec = 2 * U
    epara3 = uu
    write(6,113)capa, ec
    write(6,111)eup, edown, uu
  else
    write(6,111)eup, edown, uu
  end if
  111 format('  eup = ', f14.8, 2x, ' edown = ', f14.8, 2x, ' U = ', f14.8, 2x)
  112 format('   ac = ', f14.8, 2x)
  113 format(' capa = ', f14.8, 2x, ' ec = ', f14.8, 2x)
  call flush(6)
!
  eup   = eup   / hbar
  edown = edown / hbar
  uu    = uu    / hbar
  ec    = ec    / hbar
  ac    = ac    / hbar / 2.d0
!
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0)  ! c_up
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 2), 0.d0)  ! c_down
!
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp3, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)
  hs = hs + eup * cmtmp3 + edown * cmtmp4
!
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, &
             czero, cmtmp5, nrho)
  hs = hs + uu * cmtmp5
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho, &
             czero, cmtmp3, nrho)
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp1, nrho, &
             czero, cmtmp4, nrho)
 hs = hs + ac *( cmtmp3 +  cmtmp4)

!
!  if (capa .ge. dnano) then
!    cmtmp1 = cmtmp3 + cmtmp4  ! Ndot
!    call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
!               czero, cmtmp2, nrho) 
!    hs = hs + ec * cmtmp2
!  end if
end if
!
dmtmp1 = dble(hs) * hbar
write(6,*)
write(6,*)' real part of equilibrium Hamiltonian '
call amatout(nrho, nrho, dmtmp1)
call flush(6)
!
end subroutine calchs
