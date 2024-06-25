subroutine calc_spectral_band(ispin, ialf, cin, cout)
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
! evaluate reservoir spectral function of multiple Lorentzian functions
! enter only for lband = .true.
!
integer    :: ispin, ialf, ni
complex*16 :: cin, cout
!
cout = czero
do ni=1,nband
   cout = cout + band_coef(ni,ispin,ialf) * band_width(ni,ispin,ialf)**2 /     &
          ( (cin - band_cent(ni,ispin,ialf))**2 + band_width(ni,ispin,ialf)**2 )
end do
!
end subroutine calc_spectral_band

!---------------

subroutine zxip_band(malf)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer :: norbs2, malf, istat
integer :: ialf, iorbs, iorbs2, ispin, ni, nj
!
namelist / coupling / readcp, readmat, lrcmplx
!
norbs2 = norbs**2
if (offcor) then
    zxip(1:norbs2,1:nspin,1:malf) = czero
    do ialf=1,malf
       do iorbs=1,norbs
          ni = (iorbs - 1) * norbs + iorbs
          zxip(ni,1:nspin,ialf) =  linewidth(ialf,iorbs) * 5.d-1 
       end do
    end do
    readcp  = .false.
    readmat = .false.
    lrcmplx = .false.
    rewind(5)
    read(5, coupling, end=112)
    112 continue
    if (readcp) then
      do ialf=1,malf
        dmtmp2(1:norbs,1:norbs) = 0.d0
        if (readmat) then
          do ni=1,norbs
            if (lrcmplx) then
               read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs), (dmtmp2(ni,nj), nj=1,norbs)
            else
               read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs)
            end if
            if (istat .ne. 0) then
               write(6,*)'zxip_band: error reading linewidth, ialf = ', ialf, ' iorbs = ', ni
               stop
            end if
          end do
        else
          if (lrcmplx) then
             read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,norbs), nj=1,norbs), ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs)
          else 
             read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,norbs), nj=1,norbs)
          end if
          if (istat .ne. 0) then
             write(6,*)'zxip_band: error reading linewidth matrix for ialf = ', ialf
             stop
          end if
        end if
        do nj=1,norbs
           do ni=1,norbs
              zxip((nj-1)*norbs+ni,1:nspin,ialf) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj)) * 5.d-1
           end do
        end do
      end do
    end if
!
else
    do ialf=1,malf
       do iorbs=1,norbs
          zxip(iorbs,1:nspin,ialf) = linewidth(ialf,iorbs) * 5.d-1
       end do
    end do
end if
!
end subroutine zxip_band

