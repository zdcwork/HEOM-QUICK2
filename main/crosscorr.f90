subroutine crosscorr(cout, ialf, iorbs, cin)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
integer, intent(in) :: ialf, iorbs
complex*16, intent(in)  :: cin
complex*16, intent(out) :: cout
real*8     :: aa, phiphi
complex*16 :: ctmp1, ztmp1, ztmp2, ztmp3
!
! this subroutine works for "megaflux"
! feature : maintain the Hermicity of output with any complex cin
!           therefore, when cin returns to be real, ai and phi_i should be zero.
!
if (.not. megaflux) then
  write(6,*)
  write(6,*)' error! wrong entry for crosscorr! '
  stop
end if
if (ialf .eq. 1) then
  aa     = aoffL
  phiphi = phioffL
else
  aa     = aoffR
  phiphi = phioffR
end if
if (lafreq) then  ! dcmplx(ar,ai)
  if (natype .eq. 1) then
    cout = aa * cdexp(- (cin - aomega)**2 / anu**2 )
  else
    write(6,*)
    write(6,*)' error! unknown natype used ', natype
    stop
  end if
else
  cout = dcmplx(aa, 0.d0)
end if
if (lphifreq) then ! dcmplx(phi_r,phi_i)
  if (nphitype .eq. 1) then
    ctmp1 = (cin - phiomega) / phinu
  else if (nphitype .eq. 2) then
    if (mfdjob .and. dabs(yshift(ialf)/phinu) .ge. 1.d0-1.d-3) then
      write(6,*)' error! in <crosscorr.f90> '
      write(6,*)' additional poles should be accounted for ', phinu, yshift(ialf)
      stop
    end if
    ctmp1 = phiomega / ((cin / phinu)**2 + 1.d0)**2
  else 
    write(6,*)
    write(6,*)' error! unknown nphitype used ', nphitype
    stop
  end if
else
  ctmp1 = dcmplx(phiphi, 0.d0)
end if
if (iorbs .eq. 2) then          
  cout = cout * cdexp(-eye * ctmp1)
else if (iorbs .eq. 3) then                   ! a*exp(i*phi) with a=dcmplx(ar,ai), phi=dcmplx(phi_r,phi_i)
  cout = cout * cdexp(eye * ctmp1)
else if (iorbs .eq. 1 .or. iorbs .eq. 4) then
  cout = cunity
else
  write(6,*)
  write(6,*)' error! wrong iorbs in crosscorr.f90 ', iorbs
  stop
end if
!write(99,*)' line1 ', iorbs, ctmp1, cout
!write(99,*)' line2 ', cdexp(-eye * ctmp1), cdexp(eye * ctmp1)
return
end subroutine crosscorr
