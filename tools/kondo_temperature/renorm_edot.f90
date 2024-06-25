program renorm
implicit none
!
integer :: ni, nj, nk
real*8  :: edot, edotp, gamma, U, crit, etmp
real*8, parameter :: pi = 3.1415926d0
!
edot = -0.2d0
U    =  2.d0
gamma = 0.1d0
!
crit = 1.d-3
!
nk = 1000
!
edotp = edot
do ni=1,nk
   etmp = edotp
   edotp = edot - gamma/ pi * log(-U / edotp)
   if (dabs(etmp - edotp) .le. crit) then
      write(6,*)'converged, iter, crit, edot, edotp '
      write(6,*)ni, crit, edot, edotp
      call flush(6)
      exit
   end if
end do
if (dabs(etmp - edotp) .gt. crit) then
   write(6,*)'failed to converge after ', nk, ' cycles'
   write(6,*)'uncertainty is ', dabs(etmp - edotp)
end if


end program renorm
