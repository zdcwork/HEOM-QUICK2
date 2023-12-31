subroutine checksparsity(cin1, cin2, nmat, ndim, dcrit, lspa, ndiff)
implicit none
!
integer,    intent(in)   :: nmat, ndim
complex*16, intent(in)   :: cin1(ndim,*), cin2(ndim,*)
logical,    intent(out)  :: lspa
real*8,     intent(in)   :: dcrit
integer,    intent(out)  :: ndiff
integer                  :: ni, nj
!
if (ndim < nmat) then
   write(6,*)'checksparsity: error! ndim < nmat ', ndim, nmat
   stop
end if
!
ndiff = 0
lspa  = .false.
!
do nj=1,nmat
   do ni=1,nmat
      if ( cdabs(cin1(ni,nj)) .gt. dcrit .and. cdabs(cin2(ni,nj)) .le. dcrit .or.   &
           cdabs(cin2(ni,nj)) .gt. dcrit .and. cdabs(cin1(ni,nj)) .le. dcrit ) then
          ndiff = ndiff + 1
      end if
   end do
end do
!
if (ndiff .eq. 0) lspa = .true.
!
return
end subroutine checksparsity
