subroutine solvediis
use diismod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8 :: lni, lnj, lnk, ltmp1, ltmp2
integer   :: ni, nj, nk, nl
integer   :: istat, iconvg, iter, ierror
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8    :: dsav1
real*8    :: dnorm1, dnorm2, ratio12
real*8    :: cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cpu7, cpu8
real*8    :: ddotlong, cdotlong
real*8    :: jleftu, jrightu, jleftd, jrightd, jleft, jright
real*8    :: occu, occd, occ
external  :: ddotlong, cdotlong
complex*16, allocatable :: sigmau(:), sigmad(:)
!
call cpu_time(cpu3)
write(6,*)
write(6,*)' entering solvediis, nrho = ', nrho
call flush(6)
!
ierror = -1
iconvg = 0
iter   = 0
ltmp1  = nunk * nrho**2
ltmp2  = ltmp1 * 2
!
write(6,*)
write(6,*)' length of diis vector ', ltmp2
!
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
!
call diis(iter, ierror)
!
if (ierror .ne. 0) then
  write(6,*)
  write(6,*)' error! DIIS not converged! '
end if
!
call calcj(jleftu, jrightu, jleftd, jrightd)
jleft  = jleftu  + jleftd
if (nalf .eq. 2) then
  jright = jrightu + jrightd
end if
!
call calcocc(occu, occd, sigmau, sigmad)
occ = occu + occd
!
write(6,*)
write(6,*)' steady state solved by DIIS method '
write(6,*)
write(6,*)' occupation, spin up = ', occu
if (nspin .ne. 1) then
  write(6,*)'                down = ', occd
end if
if (funits .eq. 0) then
  write(6,*)'              jleft  = ', jleft , ' nA '
  if (nalf .eq. 2) then
    write(6,*)'             jright  = ', jright, ' nA '
  end if
else if (funits .eq. 1) then
  write(6,*)'              jleft  = ', jleft , ' pA '
  if (nalf .eq. 2) then
    write(6,*)'             jright  = ', jright, ' pA '
  end if
end if
!
if (nspin .ne. 1) then
  write(6,*)
  write(6,*)' spin current '
  if (funits .eq. 0) then
    write(6,*)'              jleftu = ', jleftu , ' nA '
    write(6,*)'              jleftd = ', jleftd , ' nA '
    if (nalf .eq. 2) then
      write(6,*)'             jrightu = ', jrightu, ' nA '
      write(6,*)'             jrightd = ', jrightd, ' nA '
    end if
  else if (funits .eq. 1) then
    write(6,*)'              jleftu = ', jleftu, ' pA '
    write(6,*)'              jleftd = ', jleftd, ' pA '
    if (nalf .eq. 2) then
      write(6,*)'             jrightu = ', jrightu, ' pA '
      write(6,*)'             jrightd = ', jrightd, ' pA '
    end if
  end if
end if
call flush(6)
!
! output for plotting iv curve
!
if (norbs .eq. 1 .and. nspin .eq. 2) then
  write(6,*)
  if (nalf .eq. 1) then
    write(6,222)iter, engyshift(1,1)*hbar, occu, occd, jleft
  else if (nalf .eq. 2) then
    write(6,222)iter, engyshift(1,1)*hbar, engyshift(2,1)*hbar, occu, occd, jleft, jright
  end if
  call flush(6)
end if
222 format(' IVOUTPUT ', I4, 2x, 6(e14.6e2, 2x))
!
if (norbs .eq. 1 .and. nspin .eq. 1) then
  write(6,*)
  if (nalf .eq. 1) then
    write(6,222)iter, epara1, occ, jleft
  else if (nalf .eq. 2) then
    write(6,222)iter, epara1, occ, jleft, jright
  end if
  call flush(6)
end if
!
call cpu_time(cpu4)
deallocate(sigmau, sigmad, STAT=istat)
!
write(6,*)
write(6,*)' leaving solvediis ', cpu4 - cpu3
call flush(6)
end subroutine solvediis
!
