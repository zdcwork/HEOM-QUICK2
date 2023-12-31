subroutine bosefunc0(zin, zout)
implicit none
!
! 1/(1-e^(-z))
! 
! zin and zout are dimensionless numbers
!
complex*16, intent(in)  :: zin
complex*16, intent(out) :: zout
real*8                  :: dtmp0
!
if (cdabs(zin) .le. 1.d-15) then
   write(6,*)
   write(6,*)'bosefunc0: error! input energy cannot be zero ', zin
   stop
end if
!
dtmp0 = dble(zin)
!
if (dtmp0 < -200.d0) then
   zout = dcmplx(0.d0, 0.d0)
else if (dtmp0 > 200.d0) then
   zout = dcmplx(1.d0, 0.d0)
else
   zout = 1.d0 / (1.d0 - cdexp(-zin))
end if
!
return
!
end subroutine bosefunc0
