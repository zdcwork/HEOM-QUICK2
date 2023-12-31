subroutine prechebyshev
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer          :: ni, nj, nk, dim1, i, j, k, itier, ispin, ialf
integer*8        :: lni, lnj
real*8           :: tt, drmax, drmin, dimax, dimin, dpre, tdelta0, dt0
real*8           :: dtmp1, dtmp2, dtmp3
real*8, external :: bessj
real*8           :: eshift(maxalf,maxspin)
complex*16       :: cgamma_nk
namelist / chebyshev / jcheby, cheby_dpre
!
jcheby     = 10
cheby_dpre = 1.d-21
rewind(5)
read(5, chebyshev, end=101) 
101 continue
if (jcheby .le. 0) then
  write(6,*)
  write(6,*)' error! jcheby <= 0 ! ', jcheby
  stop
end if
!
dim1    = nrho**2
tt      = dnano
tdelta0 = 20.d0
dpre    = cheby_dpre
!
if (lhf) then
   call calchs(rho(1,1,1))
end if
if (igroundsteady .eq. 0) then
  cmtmp3 = hs
else 
  call calcdhs(tt, rho(1,1,1))  ! initial time is t=0^+
  cmtmp3 = hs + dhs
  !do ispin=1,lspin
  !  call getengyshift(tt, 1, ispin, eshift(1,ispin))
  !end do
  !if (nalf .eq. 2) then
  !  do ispin=1,lspin
  !    call getengyshift(tt, 2, ispin, eshift(2,ispin))
  !  end do
  !end if
  do ialf=1,nalf
     do ispin=1,nspin
        call getengyshift(tt, ialf, ispin, eshift(ialf,ispin))
     end do
  end do
end if
!
drmax = -1.d20
drmin =  1.d20
dimax = -1.d20
dimin =  1.d20
!
dmtmp1 =  dble( -eye * cmtmp3)
dmtmp2 = dimag( -eye * cmtmp3)
call h2lmln(nrho, dmtmp1, nrho, dlmat1, dim1)
call h2lmrn(nrho, dmtmp1, nrho, dlmat2, dim1)
call h2lmln(nrho, dmtmp2, nrho, dlmat3, dim1)
call h2lmrn(nrho, dmtmp2, nrho, dlmat4, dim1)
dlmat1 = dlmat1 - dlmat2  ! real part (diagonal in Liouville)
dlmat3 = dlmat3 - dlmat4  
dlmat2 = dlmat3           ! imag part (off-diagonal)
!
do itier=1,ntier
  do lni=nfirst(itier), nlast(itier)
    cgamma_nk = czero
    do ni=1,itier-1
      nk        = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
      cgamma_nk = cgamma_nk + cgamma(nk)
      if (igroundsteady .ne. 0) then
         !if (.not. leadspin) then
         !  cgamma_nk = cgamma_nk + eye * dipm(nk) * eshift(ilead(nk), 1) 
         !else
           cgamma_nk = cgamma_nk + eye * dipm(ni) * eshift(ilead(nk), mspin(nk))
         !end if
      end if
    end do
    do j=1,nrho
      do i=1,nrho
        if (i .eq. j) then
          cmtmp5(i,j) = cgamma_nk
        else
          cmtmp5(i,j) = czero
        end if
        dmtmp3(i,j) =  dble(cmtmp5(i,j))
        dmtmp4(i,j) = dimag(cmtmp5(i,j))
      end do
    end do
    call h2lmln(nrho, dmtmp3, nrho, dlmat3, dim1)  
    call h2lmln(nrho, dmtmp4, nrho, dlmat4, dim1)
    dlmat3 = dlmat3 + dlmat1  ! diagonal 
    dlmat4 = dlmat4 + dlmat2  ! off-diagonal
    do k=1,dim1
      drmax = max(drmax, dlmat3(k,k))
      drmin = min(drmin, dlmat3(k,k))
      dimax = max(dimax, dlmat4(k,k))
      dimin = min(dimin, dlmat4(k,k))
    end do
  end do
end do
!
cheby_delta = drmax - drmin
tdelta      = cheby_delta * dt
dt0         = tdelta0 / cheby_delta
write(6,*)
write(6,*)' drmax       = ', drmax
write(6,*)' drmin       = ', drmin
write(6,*)' cheby_delta = ', cheby_delta
write(6,*)' tdelta      = ', tdelta,  ' dt  = ', dt
write(6,*)
write(6,*)' recommended tdelta and time-step: '
write(6,*)' tdelta0     = ', tdelta0, ' dt0 = ', dt0
call flush(6)
!
dtmp1 = 0.d0
dtmp3 = 1.d0
lni   = 0
write(6,*)
do while (dtmp3 .ge. dpre) 
  dtmp2 = bessj(lni, tdelta)
  dtmp1 = dtmp1 + dabs(dtmp2)
  dtmp3 = dabs(dtmp2) / dtmp1
  write(6,1000)lni, dtmp2
  call flush(6)
  lni = lni + 1
end do
ncheby = lni
1000 format(1x, I12, 2x, e20.8e3)
!
write(6,*)' # of Bessel functions : ', ncheby
write(6,*)' for precision         : ', dpre
call flush(6)
!
end subroutine prechebyshev
