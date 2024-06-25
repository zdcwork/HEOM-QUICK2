subroutine plotcorr
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
integer :: ni, nj, nk
integer :: isgn, ialf, ispin, iorbs, iorbsp
real*8  :: xmin, xmax, ymin, ymax, xx, yy
real*8  :: dtmp1, dtmp2, dtmp3, dtmp4
complex*16 :: zw, ctmp1, ctmp2, ctmp3
!
if (offcor) then
   write(6,*)
   write(6,*)'plotcorr: error! plotcorr for offcor not available '
   stop
end if
if (tour .le. 0.d0) then
   write(6,*)
   write(6,*)'plotcorr: error! tour should be larger than zero ', tour
   stop
end if
if (lband) then
   write(6,*)
   write(6,*)'plotcorr: error! plotcorr for lband not available '
   stop
end if
isgn   = 1
ialf   = 1
ispin  = 1
iorbs  = 1
if (isgn .eq. 1) then
  nypmin = max(0, nypmin)
else if (isgn .eq. 2) then
  nypmax = min(nypmax, 0)
end if
if (nxpmin .ge. nxpmax .or. nypmin .ge. nypmax) then
  write(6,*)
  write(6,*)' error! wrong plotting range!   '
  write(6,*)' nxpmin, nxpmax, nypmin, nypmax '
  write(6,*)nxpmin, nxpmax, nypmin, nypmax
  stop
end if
xmin = dble(nxpmin) * dpgrid
xmax = dble(nxpmax) * dpgrid
ymin = dble(nypmin) * dpgrid
ymax = dble(nypmax) * dpgrid
!
write(6,*)
write(6,*)' plotting J(z)*f(z)*e^{izt} in complex plane '
write(6,*)' (1) isgn, ialf, ispin, iorbs                '
write(6,*)isgn, ialf, ispin, iorbs
write(6,*)' (2) plotting range : '
write(6,*)' tour         = ', tour
write(6,*)' (xmin, xmax) = ', xmin, xmax
write(6,*)' (ymin, ymax) = ', ymin, ymax
if (mfdjob) then
  write(6,*)' (3) yshift(ialf) = ', yshift(ialf)
end if
call flush(6)
!
open(unit=35, file='cwreal.data', status='unknown')
rewind(35)
open(unit=36, file='cwimag.data', status='unknown')
rewind(36)
!
if (ialf .eq. 1) then
  if (ispectral .eq. 1) then
    dtmp1 = egaul
  else if (ispectral .eq. 2) then
    dtmp1 = esincl
  end if
else if (ialf .eq. 2) then
  if (ispectral .eq. 1) then
    dtmp1 = egaur
  else if (ispectral .eq. 2) then
    dtmp1 = esincr
  end if
end if
do ni=nxpmin,nxpmax,1
  xx = dble(ni) * dpgrid
  do nj=nypmin,nypmax,1
    yy = dble(nj) * dpgrid
    zw = dcmplx(xx, yy)
    if (ispectral .eq. 0) then
      ctmp1 = 1.d0 / ((zw / bandwidth(ialf))**2 + 1.d0)
    else if (ispectral .eq. 1) then
      ctmp2 = -(zw / dtmp1)**2
      ctmp1 = cdexp(ctmp2)
    else if (ispectral .eq. 2) then
      if (cdabs(zw / dtmp1) .le. dnano) then
        ctmp1 = cunity
      else
        ctmp1 = ( cdsin(zw / dtmp1) / (zw / dtmp1) )**2
      end if
    end if
    ctmp1 = ctmp1 * linewidth(ialf,iorbs) / hbar
    call foccfunc(isgn, ialf, zw, ctmp2)
    ctmp3 = cdexp(dpm(isgn) * eye * zw * tour / hbar) * ctmp1 * ctmp2 ! dimensionless
    dtmp3 = dble(ctmp3)
    dtmp4 = dimag(ctmp3)
    if (dabs(dtmp3) .le. dpico) dtmp3 = 0.d0
    if (dabs(dtmp4) .le. dpico) dtmp4 = 0.d0
    write(35,518)ni, nj, dtmp3
    write(36,518)ni, nj, dtmp4
  end do
end do
close(35)
close(36)
518 format(I4, 2x, I4, 2x, e21.6e3, 2x)
return
end subroutine plotcorr
