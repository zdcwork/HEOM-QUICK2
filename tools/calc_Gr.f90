program calc_Gr
implicit none
!
integer, parameter :: nmax1 = 10000
real*8  :: dtmp1, dtmp2, dtmp3, ww
real*8  :: dtmp4, dtmp5, dtmp6
real*8  :: pi, hbar
integer :: istat
integer :: ni, nj, nk
complex*16 :: jwr(nmax1), jwi(nmax1)
complex*16 :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6
real*8  ::  omega(nmax1)
real*8  :: bandwidth, linewidth, edot
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
jwr(nk) = dcmplx(dtmp1, dtmp2) / hbar
jwi(nk) = dcmplx(dtmp4, dtmp5) / hbar
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
   jwr(nk) = dcmplx(dtmp1, dtmp2) / hbar
   jwi(nk) = dcmplx(dtmp4, dtmp5) / hbar
   omega(nk) = ww
end do loop0
!
close(10)
close(11)
!
write(6,*)'bandwidth, linewidth, edot '
read(*,*)bandwidth, linewidth, edot
!
open(unit=12, file='Gr_sys.dat', status='unknown')
rewind(12)
!
open(unit=13, file='self_energy_ee.dat', status='unknown')
rewind(13)
!
!open(unit=14, file='Gr_nonint.dat', status='unknown')
!rewind(14)
!
do ni=nk,1,-1
   ctmp1 = dcmplx( -dimag(jwr(ni))-dble(jwi(ni)), -dble(jwr(ni))+dimag(jwi(ni)) )
   ctmp5 = dcmplx(1.d0,0.d0) / ctmp1
   ctmp3 = linewidth * 5.d-1 * bandwidth / dcmplx( (-omega(ni)), bandwidth )
   write(12,518) -omega(ni), dble(ctmp1), dimag(ctmp1), dble(ctmp5), dimag(ctmp5)
   ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (-omega(ni)) - edot - ctmp3
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (-omega(ni)) - edot 
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (-omega(ni)) 
   write(13,518) -omega(ni), dble(ctmp2), dimag(ctmp2)
   ctmp4 = dcmplx(1.d0, 0.d0) / ((-omega(ni)) - edot - ctmp3)
   ctmp6 = -omega(ni) - edot - ctmp3
   !write(14,518) -omega(ni), dble(ctmp4), dimag(ctmp4), dble(ctmp6), dimag(ctmp6)
end do
ni = 1
if (dabs(omega(ni)) .gt. 1.d-9) then
   ctmp1 = dcmplx( dimag(jwr(ni))-dble(jwi(ni)), -dble(jwr(ni))-dimag(jwi(ni)) )
   ctmp5 = dcmplx(1.d0,0.d0) / ctmp1
   ctmp3 = linewidth * 5.d-1 * bandwidth / dcmplx( omega(ni), bandwidth )
   write(12,518) omega(ni), dble(ctmp1), dimag(ctmp1), dble(ctmp5), dimag(ctmp5)
   ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) - edot - ctmp3
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) - edot 
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) 
   write(13,518) omega(ni), dble(ctmp2), dimag(ctmp2)
   ctmp4 = dcmplx(1.d0, 0.d0) / (omega(ni) - edot - ctmp3)
   ctmp6 = omega(ni) - edot - ctmp3
   !write(14,518) omega(ni), dble(ctmp4), dimag(ctmp4), dble(ctmp6), dimag(ctmp6)
end if
do ni=2,nk
   ctmp1 = dcmplx( dimag(jwr(ni))-dble(jwi(ni)), -dble(jwr(ni))-dimag(jwi(ni)) )
   ctmp5 = dcmplx(1.d0,0.d0) / ctmp1
   ctmp3 = linewidth * 5.d-1 * bandwidth / dcmplx( omega(ni), bandwidth )
   write(12,518) omega(ni), dble(ctmp1), dimag(ctmp1), dble(ctmp5), dimag(ctmp5)
   ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) - edot - ctmp3
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) - edot 
   !ctmp2 = -dcmplx(1.d0, 0.d0) / ctmp1 + (omega(ni)) 
   write(13,518) omega(ni), dble(ctmp2), dimag(ctmp2)
   !ctmp4 = 1.d0 / (omega(ni) - edot - ctmp3)
   ctmp4 = dcmplx(1.d0, 0.d0) / (omega(ni) - edot - ctmp3)
   ctmp6 = omega(ni) - edot - ctmp3
   !write(14,518) omega(ni), dble(ctmp4), dimag(ctmp4), dble(ctmp6), dimag(ctmp6)
end do
close(12)
close(13)
!close(14)
!
518 format(5(2x, e18.6e3))
!
end program calc_Gr
