subroutine specdensfunc(ialf, ispin, iorbs, iorbs2, zfreq, zout)
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer,    intent(in)  :: ialf, ispin, iorbs, iorbs2
complex*16, intent(in)  :: zfreq
complex*16, intent(out) :: zout
integer    :: nk
real*8     :: dtmp1, dtmp2
complex*16 :: ctmp1, ctmp2
!
zout = czero
!
if (lband) then
   ctmp1 = zfreq / hbar
   call calc_spectral_band(ispin, ialf, ctmp1, ctmp2)
   zout = zlwmat(iorbs,iorbs2,ispin,ialf) * ctmp2
   return 
end if
!
if (lspcor) then
    dtmp1 = linewidth(ialf,iorbs)
    dtmp2 = dlw_sp(iorbs,ispin,ialf)
    if (lequileads .and. ialf .eq. 1) then
        do nk=2,nleads
           dtmp1 = dtmp1 + linewidth(nk,iorbs)
           dtmp2 = dtmp2 + dlw_sp(iorbs,ispin,nk)
        end do
    end if
end if
!
if (ispectral .eq. 0) then
  zout = zlwmat(iorbs,iorbs2,ispin,ialf) / ( ((zfreq - bcenter(ialf,ispin)) / bandwidth(ialf))**2 + 1.d0 )
  if (megaflux .and. (lafreq .or. lphifreq)) then
     nk = (iorbs2 - 1) * norbs + iorbs
     call crosscorr(ctmp1, ialf, nk, zfreq)
     zout = ctmp1 * linewidth(ialf,1) / ( ((zfreq - bcenter(ialf,ispin)) / bandwidth(ialf))**2 + 1.d0 )
  end if
  if (lspcor) then
      zout = zout / dtmp1 * dtmp2
  end if
else if (ispectral .eq. 1) then
  if (ialf .eq. 1) then
    zout = zlwmat(iorbs,iorbs2,ispin,ialf) * cdexp(-(zfreq / egaul)**2)
  else if (ialf .eq. 2) then
    zout = zlwmat(iorbs,iorbs2,ispin,ialf) * cdexp(-(zfreq / egaur)**2)
  end if
  if (lspcor) then
      zout = zout / dtmp1 * dtmp2
  end if
else if (ispectral .eq. 2) then
  if (ialf .eq. 1) then
    if (cdabs(zfreq / esincl) .le. dnano) then
      zout = zlwmat(iorbs,iorbs2,ispin,ialf)
    else 
      zout = zlwmat(iorbs,iorbs2,ispin,ialf) * ( cdsin(zfreq / esincl) / (zfreq / esincl) )**2 
    end if
  else if (ialf .eq. 2) then
    if (cdabs(zfreq / esincr) .le. dnano) then
      zout = zlwmat(iorbs,iorbs2,ispin,ialf)
    else 
      zout = zlwmat(iorbs,iorbs2,ispin,ialf) * ( cdsin(zfreq / esincr) / (zfreq / esincr) )**2 
    end if
  end if
  if (lspcor) then
      zout = zout / dtmp1 * dtmp2
  end if
end if
!
return
end subroutine specdensfunc
!

subroutine foccfunc(isgn, ialf, zfreq, zout)
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer,    intent(in)  :: isgn, ialf
complex*16, intent(in)  :: zfreq
complex*16, intent(out) :: zout
complex*16              :: ctmp1, ctmp2
!
zout  = czero
ctmp2 = zfreq / dinvbeta(ialf)
if (dble(ctmp2) .lt. 15.d0) then
  ctmp1 = 1.d0 / (1.d0 + cdexp(ctmp2))
else
  ctmp1 = cdexp(-ctmp2) / (1.d0 + cdexp(-ctmp2))
end if
if (isgn .eq. 1) then
  zout = ctmp1
  return
else
  zout = cunity - ctmp1
  return
end if
end subroutine foccfunc
