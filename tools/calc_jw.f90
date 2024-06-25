program calc_jw
implicit none
!
integer, parameter :: nmax1 = 10000
real*8  :: dtmp1, dtmp2, dtmp3, ww
real*8  :: dtmp4, dtmp5, dtmp6
real*8  :: pi, hbar
integer :: istat
integer :: ni, nj, nk
complex*16 :: jw(nmax1)
real*8  ::  omega(nmax1)
!
pi   = 3.1415926535d0
hbar = 0.65822d0
!
open(unit=10, file='jwa.dat', status='unknown') ! half-fourier transform of Re[J(t)]
rewind(10)
open(unit=11, file='jwb.dat', status='unknown') ! half-fourier transform of Im[J(t)]
rewind(11)
!
read(10, *, iostat=istat) ww, dtmp1, dtmp2, dtmp3
if (istat .ne. 0) then
   write(6,*)' error reading <jwa.dat> 1 '
   stop
end if
read(11, *, iostat=istat) ww, dtmp4, dtmp5, dtmp6
if (istat .ne. 0) then
   write(6,*)' error reading <jwb.dat> 1 '
   stop
end if
!
nk = 1
!
jw(nk) = dcmplx(dtmp1, dtmp5) / pi / hbar
omega(nk) = ww
!
loop0: do
   read(10, *, iostat=istat) ww, dtmp1, dtmp2, dtmp3
   if (istat .ne. 0) then
      exit loop0
   end if 
   read(11, *, iostat=istat) ww, dtmp4, dtmp5, dtmp6
   if (istat .ne. 0) then
      write(6,*)' error reading <jwb.dat> ', nk + 1
      stop
   end if
   nk = nk + 1   
   jw(nk) = dcmplx(dtmp1, dtmp5) / pi / hbar
   omega(nk) = ww
end do loop0
!
close(10)
close(11)
!
open(unit=12, file='jw_sys.dat', status='unknown')
rewind(12)
do ni=nk,1,-1
   write(12,518) -omega(ni), dble(jw(ni))-dimag(jw(ni))
end do
ni = 1
if (dabs(omega(ni)) .gt. 1.d-9) then
   write(12,518) omega(ni), dble(jw(ni))+dimag(jw(ni))
end if
do ni=2,nk
   write(12,518) omega(ni), dble(jw(ni))+dimag(jw(ni))
end do
close(12)
!
518 format(5(2x, e18.6e3))
!
end program calc_jw
