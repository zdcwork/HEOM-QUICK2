subroutine fermi(dengy, mu, temp, dout)
implicit none
real*8, intent(in)  :: dengy, mu, temp  ! all in same (energy) unit
real*8, intent(out) :: dout
real*8              :: x
!
x = (dengy - mu) / temp
if (x > 200.0) then 
   dout = 0.d0
else if (x < -200.0) then
   dout = 1.d0
else
   dout = 1.d0 / (1.d0 + dexp(x))
end if
return
end subroutine fermi
