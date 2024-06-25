program testmemory
implicit none
!
integer    :: ni, nj
integer*1  :: nni
integer*8  :: lni
real*8     :: dtmp1
complex*16 :: ctmp1
!
intrinsic dcmplx
!

lni = 24
ni = lni
write(6,*)lni, ni
stop



nj = 2
ni=-1
lni=-1
nni=(-1)**(nj-1)
!
dtmp1 = 1.d10
ctmp1 = dcmplx(dtmp1,dtmp1)
!


write(6,*)'ni  ', ni, kind(ni), sizeof(ni)
write(6,*)'nni ', nni, kind(nni), sizeof(nni)
write(6,*)'lni ', lni, kind(lni), sizeof(lni)
write(6,*)'dtmp1 ', dtmp1, kind(dtmp1), sizeof(dtmp1)
write(6,*)'ctmp1 ', ctmp1, kind(ctmp1), sizeof(ctmp1)
!
open(unit=10, file='integer_1.unf', form='unformatted', status='unknown')
open(unit=20, file='integer_1.bin', form='binary', status='unknown')
rewind(10)
rewind(20)
do nj=1,1024
   write(10)nni
   write(20)nni
end do
close(10)
close(20)
!
open(unit=11, file='integer_4.unf', form='unformatted', status='unknown')
open(unit=21, file='integer_4.bin', form='binary', status='unknown')
rewind(11)
rewind(21)
do nj=1,1024
   write(11)ni
   write(21)ni
end do
close(11)
close(21)
!
open(unit=12, file='integer_8.unf', form='unformatted', status='unknown')
open(unit=22, file='integer_8.bin', form='binary', status='unknown')
rewind(12)
rewind(22)
do nj=1,1024
   write(12)lni
   write(22)lni
end do
close(12)
close(22)
!
open(unit=13, file='real_8.unf', form='unformatted', status='unknown')
open(unit=23, file='real_8.bin', form='binary', status='unknown')
rewind(13)
rewind(23)
do nj=1,1024
   write(13)dtmp1
   write(23)dtmp1
end do
close(13)
close(23)
!
open(unit=14, file='complex_16.unf', form='unformatted', status='unknown')
open(unit=24, file='complex_16.bin', form='binary', status='unknown')
rewind(14)
rewind(24)
do nj=1,1024
   write(14)ctmp1
   write(24)ctmp1
end do
close(14)
close(24)
!
end program testmemory
