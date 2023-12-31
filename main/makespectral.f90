subroutine makespectral
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer             :: ialf, ispin, mu, nu, imats, isgn, istat, iorbs, iorbs2
integer             :: ni, nj, nk
real*8, allocatable :: etmp1(:), etmp2(:)
real*8              :: dtmp1, dtmp2
complex*16          :: ctmp1, ctmp2, ctmp3
namelist / coupling / readcp, readmat, lrcmplx
!
! Build zlwmat for ispectral != 0 cases,           or
!                  ispectral .eq. 0 .and. megaflux
!
if (ispectral .eq. 0) then
  if (megaflux .and. (lafreq .or. lphifreq)) then
    do ialf=1,nalf
       zlwmat(1:norbs,1:norbs,1:nspin,ialf) = dcmplx(linewidth(ialf,1), 0.d0)
    end do
  end if
else 
  zlwmat(1:norbs,1:norbs,1:nspin,1:nalf) = czero
  if (offcor) then
    do ialf=1,nalf
      do iorbs=1,norbs
        zlwmat(iorbs,iorbs,1:nspin,ialf) = dcmplx(linewidth(ialf,iorbs), 0.d0)
      end do
    end do
    readcp  = .false.
    readmat = .false.
    lrcmplx = .false.
    rewind(5)
    read(5, coupling, end=112)
    112 continue
    if (readcp) then
      do ialf=1,nalf
        dmtmp2(1:norbs,1:norbs) = 0.d0
        if (readmat) then
          do ni=1,norbs
            if (lrcmplx) then
               read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs), (dmtmp2(ni,nj), nj=1,norbs)
            else
               read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs)
            end if
            if (istat .ne. 0) then
               write(6,*)'makespectral: error reading linewidth matrix, ialf, iorbs ', ialf, ni
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
             write(6,*)'makespectral: error reading linewidth matrix for ialf ', ialf
             stop
          end if
        end if
        do nj=1,norbs
          do ni=1,norbs
            zlwmat(ni,nj,1:nspin,ialf) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
          end do
        end do
      end do
    end if
  else
    do ialf=1,nalf
      do iorbs=1,norbs
        zlwmat(iorbs,iorbs,1:nspin,ialf) = dcmplx(linewidth(ialf,iorbs), 0.d0)
      end do
    end do
  end if
  if (doubledot) then  ! iorbs=1 only coupled to ialf=1, norbs only to ialf=2
    zlwmat = czero
    zlwmat(1,1,1:nspin,1)         = dcmplx(linewidth(1,1), 0.d0)
    zlwmat(norbs,norbs,1:nspin,2) = dcmplx(linewidth(2,norbs), 0.d0)
  end if
end if
!
! Important : make sure spec. dens. for each freq. is non-nega.-definite.
!
if (nfreq .gt. 0) then
  spectral = czero
  allocate(etmp1(nalf), etmp2(nalf), STAT=istat)
!
  write(6,*)
  if (ispectral .eq. 0) then
    write(6,*)' Lorentzian func. is used for lead spec. dens. func. '
  else if (ispectral .eq. 1) then
    write(6,*)'   Gaussian func. is used for lead spec. dens. func. '
    write(6,*)'   egaul = ', egaul
    etmp1(1) = egaul
    if (nalf .ne. 1) then
      write(6,*)'   egaur = ', egaur 
      etmp1(2) = egaur
    end if
  else if (ispectral .eq. 2) then
    write(6,*)'Square-sinc func. is used for lead spec. dens. func. '
    etmp1(1) = esincl
    if (nalf .ne. 1) then
      etmp1(2) = esincr
    end if
  else
    write(6,*)' error in <makespectral>! unknown ispectral ', ispectral
    stop
  end if
  call flush(6)
!
! # of all lead states : 
!      integration of dos(energy) within energy interval under consideration
!
  etmp2(1:nalf) = 0.d0
  do imats=1,nfreq
    if (ispectral .eq. 0) then
      do ialf=1,nalf
        etmp2(ialf) = etmp2(ialf) + dkln(imats) / ((wkln(imats) / bandwidth(ialf))**2 + 1.d0)
      end do
    else if (ispectral .eq. 1) then
      do ialf=1,nalf
        etmp2(ialf) = etmp2(ialf) + dkln(imats) * dexp(-(wkln(imats) / etmp1(ialf))**2)
      end do
    else if (ispectral .eq. 2) then
      do ialf=1,nalf
        if (dabs(wkln(imats) / etmp1(ialf)) .le. dnano) then
          dtmp1 = 1.d0
        else
          dtmp1 = dsin(wkln(imats) / etmp1(ialf)) / (wkln(imats) / etmp1(ialf))
        end if
        etmp2(ialf) = etmp2(ialf) + dkln(imats) * dtmp1**2
      end do
    else
      write(6,*)' error! unknown ispectral ', ispectral
      stop
    end if
  end do
  write(6,*)
  write(6,*)' int[dos(energy),energy] = ', (etmp2(ialf), ialf=1,nalf)
  call flush(6)
!
  do isgn=1,nsgn
    do ialf=1,nalf
      do ispin=1,nspin
        do iorbs=1,norbs
          do iorbs2=1,norbs
            do imats=1,nfreq
              ctmp1 = dcmplx(wkln(imats), dpm(isgn) * yshift(ialf))
              if (ispectral .eq. 0) then
                ctmp2 = cunity / ((ctmp1 / bandwidth(ialf))**2 + cunity)
              else if (ispectral .eq. 1) then
                ctmp3 = -(ctmp1 / etmp1(ialf))**2
                ctmp2 = cdexp(ctmp3)
              else if (ispectral .eq. 2) then
                if (cdabs(ctmp1 / etmp1(ialf)) .le. dnano) then
                  ctmp2 = cunity
                else
                  ctmp2 = ( cdsin(ctmp1 / etmp1(ialf)) / (ctmp1 / etmp1(ialf)) )**2
                end if
              end if
              if (isgn .eq. 1) then
                spectral(isgn,ialf,ispin,iorbs,iorbs2,imats) = zlwmat(iorbs2,iorbs,ispin,ialf) * ctmp2
              else
                spectral(isgn,ialf,ispin,iorbs,iorbs2,imats) = zlwmat(iorbs,iorbs2,ispin,ialf) * ctmp2
              end if
              if (ispectral .eq. 0 .and. megaflux .and. (lafreq .or. lphifreq)) then
                if (isgn .eq. 1) then
                  nk = (iorbs - 1) * norbs + iorbs2   ! transpose for plus
                else
                  nk = (iorbs2 - 1) * norbs + iorbs
                end if
                call crosscorr(ctmp2, ialf, nk, ctmp1)
                spectral(isgn,ialf,ispin,iorbs,iorbs2,imats) = spectral(isgn,ialf,ispin,iorbs,iorbs2,imats) * ctmp2
              end if
            end do
          end do
        end do
      end do
    end do
  end do
  deallocate(etmp1, etmp2, STAT=istat)
end if
end subroutine makespectral
