subroutine calc_exp_hmat_diag(ndim, zmatin, dim0, dalpha, zmatout, n1, n2)
implicit none
!
! calculate normalized exponential 
!   rho = exp(-dalpha * zmatin) / trace(exp(...))
! by diagonalizing Hermitian matrix zmatin
!
! dalpha is actually inverse temperature, do not use it
! instead, use a reference temperature 1/Delta (Delta is the energy range of zmatin)
! to allow excitation within the whole energy spectrum
!
integer, intent(in) :: ndim, dim0, n1, n2
real*8,  intent(in) :: dalpha
complex*16, intent(in)  :: zmatin(dim0,*)
complex*16, intent(out) :: zmatout(dim0,*)
!
integer :: ni, nj, nk, istat, info, nwork
logical :: lherm
real*8  :: dtmp1, dtmp2, emin, dinvt
complex*16 :: ctmp1
real*8,     allocatable :: rwork(:), tmpval(:)
complex*16, allocatable :: zwork(:,:)
complex*16, allocatable :: cmat1(:,:), cmat2(:,:)
!
if (dalpha .lt. 0.d0) then
    write(6,*)
    write(6,*)'calc_exp_hmat_diag: error! negative dalpha ', dalpha
    stop
end if
!
allocate(cmat1(ndim,ndim), cmat2(ndim,ndim), STAT=istat)
!
call checkhermicity(zmatin, dim0, cmat1, ndim, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)
    write(6,*)'calc_exp_hmat_diag: zmatin is not Hermitian ', dtmp1
    stop
end if
!
allocate(zwork(ndim,ndim), rwork(3*ndim-2), tmpval(ndim), STAT=istat)
nwork = ndim**2
cmat1(1:ndim,1:ndim) = zmatin(1:ndim,1:ndim)
call zheev('V', 'u', ndim, cmat1, ndim, tmpval, zwork, nwork, rwork, info)
if (info .ne. 0) then
    write(6,*)'calc_exp_hmat_diag: error! diagonalization failed 1 ', info
    stop
end if
!
!do nk=1,ndim
!   write(6,*)nk, tmpval(nk)
!end do
!stop
!
zwork(1:ndim,1:ndim) = cmat1(1:ndim,1:ndim)
emin  = tmpval(n1)
dinvt = 1.d0 / (tmpval(n2) - emin + 1.d-6)
cmat2(1:ndim,1:ndim) = dcmplx(0.d0, 0.d0)
!do nk=1,ndim
do nk=n1,n2,1
   !dtmp1 = dexp(-dalpha * (tmpval(nk) - emin))
   dtmp1 = dexp(-dinvt * (tmpval(nk) - emin))
   do nj=1,ndim
      ctmp1 = dconjg(zwork(nj,nk)) * dtmp1
      do ni=1,ndim
         !cmat2(ni,nj) = cmat2(ni,nj) + cdabs(zwork(ni,nk) * ctmp1)
         cmat2(ni,nj) = cmat2(ni,nj) + zwork(ni,nk) * ctmp1
      end do
   end do
end do
dtmp2 = 0.d0
do ni=1,ndim
   dtmp2 = dtmp2 + dble(cmat2(ni,ni))
end do
cmat2(1:ndim,1:ndim) = cmat2(1:ndim,1:ndim) / dtmp2
!
zmatout(1:ndim,1:ndim) = cmat2(1:ndim,1:ndim)
!
deallocate(zwork, rwork, tmpval, STAT=istat)
deallocate(cmat1, cmat2, STAT=istat)
return
end subroutine calc_exp_hmat_diag
