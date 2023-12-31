subroutine solvebicg2
use bicgmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8 :: lni, lnj, lnk, ltmp1, ltmp2, dimvec
integer   :: ni, nj, nk, nl, i, j
integer   :: istat, iconvg, iter, ierror
logical   :: iexist
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8    :: dsav1, errini, erriter
real*8    :: dnorm1, dnorm2, ratio12
real*8    :: cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cpu7, cpu8
real*8    :: dprev1, dprev2, dprev3
real*8    :: ddotlong, cdotlong
real*8    :: alphak, betak, res1, res2, occu, occd, occ
external  :: ddotlong, cdotlong
complex*16, allocatable :: sigmau(:), sigmad(:)
!
call cpu_time(cpu3)
write(6,*)
write(6,*)' entering solvebicg2, nrho = ', nrho
call flush(6)
!
inquire(file='bicg.xtmp', exist=iexist, err=998)
998 continue
if (iexist) then
  open(unit=24, file='bicg.xtmp', status='unknown', form='unformatted')
  close(unit=24, status="delete")
end if
!
dimvec = nunk * nrho**2
open(unit=24, file='bicg.xtmp', status='unknown', access='direct', form='unformatted', recl=16*dimvec)
!
! 1st record : sk
! 2nd record : rk
! 3rd record : dk
! 4th record : fk
! 5th record : rho
! 6th record : rcgrhs
!
ierror = -1
iconvg = 0
iter   = 0
ltmp1  = nunk * nrho**2
ltmp2  = ltmp1 * 2
dprev3 = crit
!
write(6,*)
write(6,*)' length of bicg vector ', ltmp2
!
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
!---------------------debug
! initialize
!
call pre_atrx_bicg
!
write(24,rec=5)(((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
do ni=1,nrho
  do nj=ni+1,nrho
    rho(ni,nj,1) = czero
  end do
end do
do ni=1,nrho
  rho(ni,ni,1) = dcmplx(dble(rho(ni,ni,1)), 0.d0)
end do
!
rcgtmp = rho
call calc_ax_bicg(nrho, rcgtmp, rcgrhs)
!
do lni=1,nunk
  do ni=1,nrho
    do nj=1,nrho
      if (lni .eq. 1 .and. ni .eq. 1 .and. nj .eq. 1) then
        rho(1,1,1)     = cunity - rcgrhs(1,1,1)
      else
        rho(ni,nj,lni) = -rcgrhs(ni,nj,lni)
      end if
    end do
  end do
end do
dnorm1 = cdotlong(ltmp1, rho(1,1,1), rho(1,1,1))
!
dtmp1 = 0.d0
do lni=1,nunk
  do ni=1,nrho
    do nj=1,nrho
      dtmp1 = max(dtmp1, dabs(dble(rho(ni,nj,lni))))
      dtmp1 = max(dtmp1, dabs(dimag(rho(ni,nj,lni))))
    end do
  end do
end do
write(6,*)' max(rk) = ', dtmp1
errini = dtmp1
!
do ni=1,4
  write(24, rec=ni) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
end do
!
write(6,*)
call cpu_time(cpu1)
!
100 continue
!
read(24, rec=3) (((rcgtmp(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
call cpu_time(cpu5)
call calc_ax_bicg(nrho, rcgtmp, rcgrhs)
call cpu_time(cpu6)
write(24, rec=6) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
read(24, rec=4) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
res2 = cdotlong(ltmp1, rho(1,1,1), rcgrhs(1,1,1))
read(24, rec=1) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
read(24, rec=2) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
res1 = cdotlong(ltmp1, rho(1,1,1), rcgrhs(1,1,1))
alphak = res1 / res2
dsav1  = res1
!
read(24, rec=3) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
read(24, rec=5) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
rho = rho + rcgrhs * alphak
do ni=1,nrho
  do nj=ni+1,nrho
    rho(ni,nj,1) = dconjg(rho(nj,ni,1))
  end do
end do
write(24, rec=5) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
read(24, rec=6) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
read(24, rec=2) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
rho = rho - rcgrhs * alphak 
write(24, rec=2) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
dnorm2 = cdotlong(ltmp1, rho(1,1,1), rho(1,1,1))
!
read(24, rec=4) (((rcgtmp(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
call cpu_time(cpu7)
call calc_atrx_bicg(nrho, rcgtmp, rcgrhs)
call cpu_time(cpu8)
!
read(24, rec=1) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
rho = rho - rcgrhs * alphak
write(24, rec=1) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
ratio12 = dnorm2 / dnorm1
!
!if (ratio12 .le. dprev3) iconvg = 1
erriter = dsqrt(dnorm2)
if (erriter .le. dprev3) iconvg = 1
!
read(24, rec=5) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
call calcocc(occu, occd, sigmau, sigmad)
!
call cpu_time(cpu2)
if (nspin .eq. 1) then
!  write(6,518)iter, res1, res2, ratio12, occu, cpu2 - cpu1
  write(6,518)iter, res1, res2, erriter, occu, cpu2 - cpu1
else 
!  write(6,518)iter, res1, res2, ratio12, occu, occd, cpu2 - cpu1, cpu6 - cpu5, cpu8 - cpu7
  write(6,518)iter, res1, res2, erriter, occu, occd, cpu2 - cpu1, cpu6 - cpu5, cpu8 - cpu7
end if
call flush(6)
call cpu_time(cpu1)
518 format(I6, 2x, 10(e12.4e2, 2x))
!
if (iconvg .eq. 1) then
 write(6,*)
 write(6,*)' converged'
 call flush(6)
 ierror = 0
 goto 200
end if
!
read(24, rec=1) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
read(24, rec=2) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
res1 = cdotlong(ltmp1, rho(1,1,1), rcgrhs(1,1,1))
res2 = dsav1
!
betak = res1 / res2
!
read(24, rec=3) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
rho = rcgrhs + rho * betak
write(24, rec=3) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
read(24, rec=4) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
read(24, rec=1) (((rcgrhs(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
rho = rcgrhs + rho * betak
write(24, rec=4) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
!
111 format(I4, 2x, ' check ', 2x, 6(e15.6e3,2x))
call flush(6)
!
iter = iter + 1
if (iter .gt. maxit0) then
 write(6,*)
 write(6,*)' error! maxit0 reached without convergency ', iter
! stop
 goto 200
end if
!
goto 100
!
200 continue
!
read(24, rec=5) (((rho(i,j,lni), i=1,nrho), j=1,nrho), lni=1,nunk)
close(unit=24, status="delete")
!
!
write(6,*)
write(6,*)'solvebicg2: steady state solved by BICG method '
!
call outsteady(iter)
!
call cpu_time(cpu4)
deallocate(sigmau, sigmad, STAT=istat)
!
write(6,*)
write(6,*)' leaving solvebicg2 ', cpu4 - cpu3
call flush(6)
end subroutine solvebicg2
!
