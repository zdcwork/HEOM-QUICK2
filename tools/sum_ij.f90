program sum_ij
implicit none
!
integer :: ncount, nalf, n1, n2, istat
real*8  :: t1, t2, dt, tt
real*8, allocatable :: di(:), dj(:), dje(:), dtmp(:)


!
open(7, file='Iw.dat', status='unknown')
open(8, file='JHw.dat', status='unknown')
open(9, file='JEw.dat', status='unknown')
!
rewind(7)
rewind(8)
rewind(9)
!
read(7,*)ncount, nalf
!
if (ncount .le. 1 .or. nalf .le. 0) then
   write(*,*)' error! ncount, nalf ', ncount, nalf
   stop
end if
!
allocate(di(nalf), dj(nalf), dje(nalf), dtmp(nalf), STAT=istat)
di = 0.d0
dj = 0.d0
dje = 0.d0
dtmp = 0.d0
!
do n1=1,ncount
   read(7,*) tt, (dtmp(n2), n2=1,nalf)
   if (n1 .eq. 1) t1 = tt
   if (n1 .eq. 2) t2 = tt 
   di(1:nalf) = di(1:nalf) + dtmp(1:nalf)
end do
dt = t2 - t1
!
do n1=1,ncount
   read(8,*) tt, (dtmp(n2), n2=1,nalf)
   dj(1:nalf) = dj(1:nalf) + dtmp(1:nalf)
end do
!
do n1=1,ncount
   read(9,*) tt, (dtmp(n2), n2=1,nalf)
   dje(1:nalf) = dje(1:nalf) + dtmp(1:nalf)
end do
!
di(1:nalf) = di(1:nalf) * dt
dj(1:nalf) = dj(1:nalf) * dt
dje(1:nalf) = dje(1:nalf) * dt
!
close(7)
close(8)
close(9)
!
write(*,*)
write(*,*)' electronic current (nA) '
write(*,1000) (di(n1), n1=1,nalf)
write(*,*)
write(*,*)' heat current (pW) '
write(*,1000) (dj(n1), n1=1,nalf)
write(*,*)
write(*,*)' energy current (pW) '
write(*,1000) (dje(n1), n1=1,nalf)
!
1000 format(6(2x, e15.6e3))

!
deallocate(di, dj, dje, dtmp,STAT=istat)

end program sum_ij
