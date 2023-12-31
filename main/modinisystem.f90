subroutine modinisystem
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer   :: ni, nj, nk, lmax, nkmax
integer*8 :: lni, lnj
real*8    :: cpu1, cpu2
real*8    :: dtmp1, dtmp2, dtmp3, derr
!
namelist / inifield2 / vdot
namelist / inifield3 / v01, v02
namelist / inifield1 / vupdn
!
if (lsparse) then
   write(6,*)'modinisystem: error! not yet available in sparse mode '
   stop
end if
!
call cpu_time(cpu1)
!write(6,*)
!write(6,*)' entering modinisystem '
!call flush(6)
!
rho0   = rho
rhotmp = rho
cmtmp2 = czero
!
if (nspin .eq. 1 .and. norbs .eq. 1) then
  vdot = 0.d0
  rewind(5)
  read(5, inifield2, end=100) 
  100 continue
  vdot = vdot / hbar
  write(6,*)
  write(6,*)' initial external field strength vdot = ', vdot*hbar
  call flush(6)
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, 1, 1), 0.d0)
  call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(vdot, 0.d0), cmtmp1, nrho, &
             cmtmp1, nrho, czero, cmtmp2, nrho)
end if
!
if (nspin .eq. 2 .and. norbs .eq. 1) then
  vupdn = 0.d0
  rewind(5)
  read(5, inifield1, end=102)
  102 continue
  vupdn = vupdn / hbar
  write(6,*)
  write(6,*)' initial external field strength vupdn = ', vupdn*hbar    ! spin-flip term
  call flush(6)
  cmtmp5(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, 1, 1), 0.d0)  ! c_up
  cmtmp6(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, 1, 2), 0.d0)  ! c_dn
!
! real field       | cos(wt) 
!
  call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(vupdn, 0.d0), cmtmp5, nrho, &
             cmtmp6, nrho, czero, cmtmp1, nrho)
  do nj=1,nrho
    do ni=1,nrho
      cmtmp2(ni,nj) = cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni))  ! Hermitianize
    end do
  end do
!
! imaginary field  | sin(wt)
!
!  call zgemm('c', 'n', nrho, nrho, nrho,  cunity, cmtmp5, nrho, cmtmp6, nrho,  czero, cmtmp1, nrho)
!  call zgemm('c', 'n', nrho, nrho, nrho, -cunity, cmtmp6, nrho, cmtmp5, nrho, cunity, cmtmp1, nrho)
!  cmtmp2 = cmtmp1 * dcmplx(0.d0, vupdn)
end if
!
if (nspin .eq. 1 .and. norbs .eq. 2) then
  v01 = 0.d0
  v02 = 0.d0
  rewind(5)
  read(5, inifield3, end=101)
  101 continue
  v01 = v01 / hbar
  v02 = v02 / hbar
  write(6,*)
  write(6,*)' initial external field strength '
  write(6,*)' v01 = ', v01*hbar, ' v02 = ', v02*hbar
  call flush(6)
  cmtmp5(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, 1, 1), 0.d0) ! c_1
  cmtmp6(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, 2, 1), 0.d0) ! c_2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp5, nrho, cmtmp5, nrho, &
             czero, cmtmp1, nrho)
  cmtmp2 = cmtmp2 + v01 * cmtmp1
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp6, nrho, cmtmp6, nrho, &
             czero, cmtmp1, nrho)
  cmtmp2 = cmtmp2 + v02 * cmtmp1
end if
!
! V ===> cmtmp2            instantaneous change of system Hamiltonian
cmtmp2 = 5.d-1 * cmtmp2  ! the 1/2 factor comes from Heaviside step function 
                         ! Theta(t - t0) 
                         ! DO NOT CHANGE THIS MATRIX FROM NOW ON!!
!
cmtmp6 = cmtmp2 * hbar
write(6,*)
write(6,*)' Dirac-Delta-type excitation field for the reduced system (t0 --> t0+) '
call cmatout(nrho, nrho, cmtmp6, dmtmp1)
!
! V matrix (cmtmp2) should be Hermitian, check hermicity
!
derr = 0.d0
do ni=1,nrho
  do nj=1,ni
    derr = derr + cdabs(cmtmp2(ni,nj) - dconjg(cmtmp2(nj,ni)))
  end do
end do
if (derr .ge. dnano) then
  write(6,*)
  write(6,*)' error in modinisystem '
  write(6,*)' Perturbation system Hamiltonian is found not Hermitian ', derr
  stop
end if
!
nkmax  = 0
!
do lni=1,nunk
! A_1 = (-i)[V, P(t0)]
  call zhemm('l', 'u', nrho, nrho, -eye, cmtmp2, nrho, rhotmp(1,1,lni), nrho, &
              czero, cmtmp1, nrho)
  call zhemm('r', 'u', nrho, nrho,  eye, cmtmp2, nrho, rhotmp(1,1,lni), nrho, &
             cunity, cmtmp1, nrho)   
  dtmp1 = 0.d0
  dtmp2 = 0.d0
  do nj=1,nrho
    do ni=1,nrho
      dtmp2 = dtmp2 + cdabs(cmtmp1(ni,nj))
    end do
  end do
  cmtmp5 = cmtmp1
! A_{k+1} = (-i)[V, A_k]
! sum A_k over k from 1 to k_max ===> cmtmp1
!                          A_k   ===> cmtmp5
!                          A_k+1 ===> cmtmp6
  lmax = 1000
  nk   = 0
  loop1: do
    nk    = nk + 1
    dtmp1 = dtmp1 + dtmp2
    if (dabs(dtmp2/dtmp1) .le. dpico) then
      if (nkmax .le. nk) nkmax = nk
      exit loop1
    else if (nk .ge. lmax) then
      write(6,*)
      write(6,*)' error in modinisystem '
      write(6,*)' maximal terms reached, ', lmax, ' not converged '
      stop
    end if
    call zhemm('l', 'u', nrho, nrho, -eye, cmtmp2, nrho, cmtmp5, nrho,  czero, cmtmp6, nrho)
    call zhemm('r', 'u', nrho, nrho,  eye, cmtmp2, nrho, cmtmp5, nrho, cunity, cmtmp6, nrho)
    cmtmp5 = cmtmp6                !   A_k => A_k+1
    cmtmp1 = cmtmp1 + cmtmp6       !   sum = sum + A_k 
    dtmp2  = 0.d0
    do nj=1,nrho
      do ni=1,nrho
        dtmp2 = dtmp2 + cdabs(cmtmp6(ni,nj))
      end do
    end do
  end do loop1
!
  rho0(1:nrho, 1:nrho, lni) = rho0(1:nrho, 1:nrho, lni) + cmtmp1(1:nrho, 1:nrho)
end do
!
cmtmp6(1:nrho, 1:nrho) = rho0(1:nrho, 1:nrho, 1) - rho(1:nrho, 1:nrho, 1)  ! rho(0^+) - rho(0) should be zero
rho = rho0
!
dtmp1 = 0.d0
do nj=1,nrho
  do ni=1,nrho
    dtmp1 = max(dtmp1, cdabs(cmtmp6(ni,nj)))
  end do
end do
!
write(6,*)
write(6,*)' expansion for time-ordered propagator converged '
write(6,*)' in modinisystem after ', nkmax, ' cycles        '
write(6,*)' max(rho(0^+) - rho(0)) = ', dtmp1
call flush(6)
!
call cpu_time(cpu2)
!write(6,*)
!write(6,*)' leaving modinisystem ', cpu2 - cpu1
!call flush(6)
!
end subroutine modinisystem
