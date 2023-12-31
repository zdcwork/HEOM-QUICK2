subroutine zmat_trace(cin, dim0, ndim, zout)
implicit none
!
integer, intent(in)     :: dim0, ndim
complex*16, intent(in)  :: cin(dim0,*)
complex*16, intent(out) :: zout
!
integer :: ni
!
zout = dcmplx(0.d0, 0.d0)
do ni=1,ndim
   zout = zout + cin(ni,ni)
end do
!
return
end subroutine zmat_trace
