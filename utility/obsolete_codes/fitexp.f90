subroutine fitexp(iorbs, ispin, ialf)
use matmod
use tmpmatmod
use fitmod
use random
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(in) :: iorbs, ispin, ialf
integer             :: ni, nj, nk, istat, imats, prin
integer*8           :: nmax1, nmax0
real*8              :: dpre0, dtmp1, dtmp2, dtmp3
real*8              :: cpu1, cpu2, cpu3, cpu4
real*8              :: tmpmatsu, tmpb, tmpgama, dintm0, dintmn, rerror
real*8, external    :: ferror, dlamch, praxis
real*8              :: t0, machep, h0, fmin, final
!
!write(6,*)
!write(6,*)' entering fitexp '
!call flush(6)
call cpu_time(cpu1)
!
dpre0 = dnano
!
! Evaluate # of Matsubara terms for "exact" C_{\mu s}^{\alpha \sigma}(t).
! Consider band-center = equilibrium fermi energy = 0, where 
!  "cb"s are pure imaginary and "cgama" are real numbers.
!
dtmp1 = 0.d0
imats = 1
do
  tmpmatsu = (2.d0 * dble(imats) - 1.d0) * pi * dinvbeta(ialf) 
  tmpb     =  -2.d0 * dxip(iorbs,ispin,ialf) * dinvbeta(ialf) /      &
              (gammap(iorbs,ispin,ialf)**2 - tmpmatsu**2) / hbar**2
  tmpgama  = -tmpmatsu / hbar
  dtmp2    = dabs(tmpb / tmpgama)
  dtmp1    = dtmp1 + dtmp2
  dtmp2    = dtmp2 / dtmp1
  if (dtmp2 .le. dpre0) then
    nmax0 = imats 
    exit
  else
    imats = imats + 1
  end if
end do
!
if (detail) then
  write(6,*)
  write(6,*)' # of Matsubara terms for exact correlation function: ', nmax0
  write(6,*)' # of Matsubara terms considered explicilty for heom: ', nmats0
end if
if (nmax0 .le. nmats0) then
  write(6,*)
  write(6,*)' error! nmax0 <= nmats0 found, no need for fitexp ', nmax0, nmats0
  stop
end if
call flush(6)
nmax1 = nmax0 - nmats0
!
if (nfitexp .ge. nmax1) then
  write(6,*)
  write(6,*)' error! wrong nfitexp in fitexp ', nfitexp, nmax1
  write(6,*)' iorbs, ispin, ialf : ', iorbs, ispin, ialf
  stop
end if
noriexp = nmax1
allocate(amatsu(nmax0), ab(nmax0), agama(nmax0), STAT=istat)
!
do imats=nmats0+1,nmax0
  amatsu(imats - nmats0) = (2.d0 * dble(imats) - 1.d0) * pi * dinvbeta(ialf)
  ab(imats - nmats0)     = -2.d0 * dxip(iorbs,ispin,ialf) * dinvbeta(ialf) /         &
                           (gammap(iorbs,ispin,ialf)**2 - amatsu(imats - nmats0)**2) &
                           / hbar**2
end do
agama(1:nmax1) = -amatsu(1:nmax1) / hbar
!
! initialize fitb and fitgama
!
fitb(1:nfitexp)    = ab(1:nfitexp)  
fitgama(1:nfitexp) = agama(1:nfitexp)
!
!  Be careful on the usage of random numbers (reproducible?)
!
!do ni=1,nfitexp
!  fitb(ni)    = random_normal()
!  fitgama(ni) = -abs(random_normal()) * 1.d2
!end do
!
fitx(1:nfitexp)           = fitb(1:nfitexp)
fitx(nfitexp+1:nfitexp*2) = fitgama(1:nfitexp)
!
t0     = fitdpre
machep = dlamch('e')
h0     = 0.d0        ! h0 is assumed mainly due to fitgama
do ni=1,nfitexp
  h0 = h0 + dabs(fitgama(nfitexp-ni+1) - agama(noriexp-ni+1))
end do
if (detail) then
  prin = 1
else 
  prin = 0
end if
!
call cpu_time(cpu3)
final = praxis(t0, machep, h0, nfitexp*2, prin, fitx, ferror, fmin)
call cpu_time(cpu4)
fitb(1:nfitexp)    = fitx(1:nfitexp)
fitgama(1:nfitexp) = fitx(nfitexp+1:nfitexp*2)
if (detail) then
  write(6,*)
  write(6,*)' mininum found by praxis : ', final 
  do imats=1,nfitexp
    write(6,1000)imats, fitb(imats), fitgama(imats)
  end do
  call flush(6)
end if
1000 format(1x, I4, 2x, 2(f14.8, 1x))
!
dintm0 = 0.d0
do ni=1,noriexp
  dintm0 = dintm0 + 5.d-1 * ab(ni)**2 / agama(ni)
  do nj=1,ni-1
    dintm0 = dintm0 + 2.d0 * ab(ni) * ab(nj) / (agama(ni) + agama(nj))
  end do
end do
dtmp1 = 0.d0
do ni=1,nfitexp
  dtmp1 = dtmp1 + 5.d-1 * fitb(ni)**2 / fitgama(ni)
  do nj=1,ni-1
    dtmp1 = dtmp1 + 2.d0 * fitb(ni) * fitb(nj) / (fitgama(ni) + fitgama(nj))
  end do
end do
dtmp2 = 0.d0
do ni=1,nfitexp
  do nj=1,noriexp
    dtmp2 = dtmp2 + fitb(ni) * ab(nj) / (fitgama(ni) + agama(nj))
  end do
end do
dintmn = dintm0 + dtmp1 - 2.d0 * dtmp2
rerror = dintmn / dintm0
if (detail) then
  write(6,*)
  write(6,*)' fitting exponential functions result '
  write(6,*)' M = ', noriexp, ' N = ', nfitexp
  write(6,*)' original residue, I(M,0) = ', dintm0
  write(6,*)' current  residue, I(M,N) = ', dintmn
  write(6,*)' relative ratio ', rerror
  call flush(6)
else
  write(6,*)
  write(6,*)' ispin, iorbs, ialf : ', ispin, iorbs, ialf
  write(6,1001)noriexp, nfitexp, final
  write(6,1002)dintm0, dintmn, rerror
  do imats=1,nfitexp
    write(6,1000)imats, fitb(imats), fitgama(imats)
  end do
  call flush(6)
end if
1001 format(' M = ', I8, 1x, ' N = ', I3, 1x, ' praxis = ', e14.6e3) 
1002 format(' I(M,0) = ', e14.6e2, 1x, ' I(M,N) = ', e14.6e2, 1x, ' ratio = ', e14.6e2)
!
deallocate(amatsu, ab, agama, STAT=istat)
call cpu_time(cpu2)
!write(6,*)
!write(6,*)' leaving fitexp  ', cpu2 - cpu1
!write(6,*)' time for praxis ', cpu4 - cpu3
!call flush(6)
end subroutine fitexp
!
!-----------------------------------
real*8 function ferror(x,n) 
use fitmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(in) :: n
real*8,  intent(in) :: x(n)
integer             :: ni, nj
integer*8           :: lnm
real*8              :: dtmp1
!
ferror             = 0.d0
fitb(1:nfitexp)    = x(1:nfitexp)
fitgama(1:nfitexp) = x(nfitexp+1:nfitexp*2)
!
do ni=1,nfitexp
  dtmp1 = 0.d0
  do lnm=1,noriexp
    dtmp1 = dtmp1 + ab(lnm) / (agama(lnm) + fitgama(ni))
  end do 
  do nj=1,nfitexp
    dtmp1 = dtmp1 - fitb(nj) / (fitgama(nj) + fitgama(ni))
  end do
  ferror = ferror + dtmp1**2
!  
  dtmp1 = 0.d0
  do lnm=1,noriexp
    dtmp1 = dtmp1 + ab(lnm) / (agama(lnm) + fitgama(ni))**2
  end do
  do nj=1,nfitexp
    dtmp1 = dtmp1 - fitb(nj) / (fitgama(nj) + fitgama(ni))**2
  end do
  ferror = ferror + dtmp1**2
end do
!
return
end function ferror
