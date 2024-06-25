subroutine calc_exp_hmat(ndim, zmatin, dim0, dalpha, kmax, zmatout)
implicit none
!
! exp(alpha*zmatin) ~= 1 + \sum_{k=1}^{kmax} (dalpha*zmatin)^k / (k!)
!
integer, intent(in) :: ndim, dim0, kmax
real*8,  intent(in) :: dalpha
complex*16, intent(in)  :: zmatin(dim0,*)
complex*16, intent(out) :: zmatout(dim0,*)
!
integer :: ni, istat
real*8  :: dtmp1
complex*16, allocatable :: cmat1(:,:), cmat2(:,:)
!
if (kmax .le. 0) then
   write(6,*)'calc_exp_hmat: error, kmax <= 0 found ', kmax
   stop
end if
!
allocate(cmat1(ndim,ndim), cmat2(ndim,ndim), STAT=istat)
!
zmatout(1:ndim,1:ndim) = dcmplx(0.d0, 0.d0)
do ni=1,ndim
   zmatout(ni,ni) = dcmplx(1.d0, 0.d0)
end do
cmat1(1:ndim,1:ndim) = zmatout(1:ndim,1:ndim) ! unit matrix
!
dtmp1 = 1.d0
!
do ni=1,kmax
   dtmp1 = dtmp1 / dble(ni)
   call zgemm('n', 'n', ndim, ndim, ndim, dcmplx(dalpha,0.d0), cmat1, ndim, zmatin, dim0, dcmplx(0.d0,0.d0), cmat2, ndim)
   zmatout(1:ndim,1:ndim) = zmatout(1:ndim,1:ndim) + dtmp1 * cmat2(1:ndim,1:ndim)
   cmat1(1:ndim,1:ndim) = cmat2(1:ndim,1:ndim)
end do
!
deallocate(cmat1, cmat2, STAT=istat)
return
end subroutine calc_exp_hmat
