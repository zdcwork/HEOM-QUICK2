program plotspectral
implicit none
!
integer :: i, j, k, ne
integer, parameter :: nmax1 = 20000
real*8  :: dtmp1, dtmp2, dtmp3
real*8  :: w0, w1, w2
real*8  :: energy(nmax1), dos(nmax1,3)
real*8  :: emin, emax, de, ee
!
w0 = 5.d0
w1 = w0 / dsqrt(dlog(2.d0))
w2 = w0 / 1.39156d0
!
write(6,*)' w0 = ', w0
write(6,*)' w1 = ', w1
write(6,*)' w2 = ', w2
!
emin = -5.d1
emax = -emin
de   = 1.d-2
ee   = emin
k    = 0
!
do while (ee .le. emax) 
  k = k + 1
  energy(k) = ee
  dos(k,1)  = 1.d0 / (1.d0 + (ee / w0)**2) 
  dos(k,2)  = dexp(-(ee / w1)**2)
  if (dabs(ee / w2) .le. 1.d-9) then
    dos(k,3) = 1.d0
  else
    dos(k,3) = (dsin(ee / w2) / (ee / w2))**2
  end if
  ee = ee + de
end do
!
open(unit=10, file='spec.dat', status='unknown')
rewind(10)
do j=1,k
  write(10,518)energy(j), (dos(j,i),i=1,3)
end do
518 format(2x, 4(e15.6e3, 2x))
close(10)
!
end program plotspectral
