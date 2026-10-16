subroutine fermifunc(isgn, ialf, cin, cout)
use matmod
use tmpmatmod
implicit none
!
! cout is actually independent of isgn, as long as cppole(ni,1) and cppole(ni,2) are 
! complex conjugate to each other 
!
include "../include/sizes"
include "../include/common"
!
! Matsubara (MSD), Polynomial (PFD) or Pade (PSD) expansion for 
! Fermi funtion 1/(e^x + 1)
!
integer    :: ni, nj, nk, ialf, isgn
complex*16 :: cin, cout
!
cout = dcmplx(5.d-1, 0.d0)
if (ldiff_spe_dec .or. ldiff_spe_dec_sp) then
  select case (ncount_scheme(ialf))
  case (1,2)
    do ni=1,nmats_alf(ialf)
       cout = cout + cpcoef_alf(ni,ialf) / (cin + cppole_alf(ni,isgn,ialf)) + &
                     cpcoef_alf(ni,ialf) / (cin - cppole_alf(ni,isgn,ialf))
    end do
  case (5,6)
    do ni=1,nmats
       cout = cout + cpcoef(ni) / (cin + cppole(ni,isgn)) + &
                     cpcoef(ni) / (cin - cppole(ni,isgn))
    end do
  case (8)
    do ni=1,nmats_alf(ialf)
       cout = cout - 1.d0 / (cin + dcmplx(0.d0, dble(2 * ni - 1) * pi)) - &
                     1.d0 / (cin - dcmplx(0.d0, dble(2 * ni - 1) * pi))
    end do
  case default
    write(6,*)'fermifunc: unsupported decomposition for lead ', ialf, ncount_scheme(ialf)
    stop
  end select
else if (pfdjob .or. psdjob .or. psdlor .or. psdfff) then
  do ni=1,nmats
     cout = cout + cpcoef(ni) / (cin + cppole(ni,isgn)) + &
                   cpcoef(ni) / (cin - cppole(ni,isgn))
  end do
else
  do ni=1,nmats
     cout = cout - 1.d0 / (cin + dcmplx(0.d0, dble(2 * ni - 1) * pi)) - &
                   1.d0 / (cin - dcmplx(0.d0, dble(2 * ni - 1) * pi))
  end do
end if
!
end subroutine fermifunc
