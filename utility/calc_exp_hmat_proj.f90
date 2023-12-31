subroutine calc_exp_hmat_proj(ndim, zmatin, dim0, dalpha, zmatout, nstate, istate)
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
integer, intent(in) :: ndim, dim0, nstate, istate(*)
real*8,  intent(in) :: dalpha
complex*16, intent(in)  :: zmatin(dim0,*)
complex*16, intent(out) :: zmatout(dim0,*)
!
integer :: ni, nj, nk, istat, info, nwork
integer :: n1, n2, nm
logical :: lherm
real*8  :: dtmp1, dtmp2, emin, dinvt
complex*16 :: ctmp1
real*8,     allocatable :: rwork(:), tmpval(:)
complex*16, allocatable :: zwork(:,:)
complex*16, allocatable :: cmat1(:,:), cmat2(:,:)
!
if (dalpha .lt. 0.d0) then
    write(6,*)
    write(6,*)'calc_exp_hmat_proj: error! negative dalpha ', dalpha
    stop
end if
!
allocate(cmat1(ndim,ndim), cmat2(ndim,ndim), STAT=istat)
!
call checkhermicity(zmatin, dim0, cmat1, ndim, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)
    write(6,*)'calc_exp_hmat_proj: zmatin is not Hermitian ', dtmp1
    stop
end if
!
allocate(zwork(ndim,ndim), rwork(3*ndim-2), tmpval(ndim), STAT=istat)
nwork = ndim**2
cmat1(1:ndim,1:ndim) = zmatin(1:ndim,1:ndim)
call zheev('V', 'u', ndim, cmat1, ndim, tmpval, zwork, nwork, rwork, info)
if (info .ne. 0) then
    write(6,*)'calc_exp_hmat_proj: error! diagonalization failed 1 ', info
    stop
end if
!
zwork(1:ndim,1:ndim) = cmat1(1:ndim,1:ndim)
n1 = istate(1)
n2 = istate(nstate)
emin  = tmpval(n1)
dinvt = 1.d0 / (tmpval(n2) - emin + 1.d-6)
cmat2(1:ndim,1:ndim) = dcmplx(0.d0, 0.d0)
do nm=1,nstate
   nk = istate(nm)
   dtmp1 = dexp(-dinvt * (tmpval(nk) - emin))
   do nj=1,ndim
      ctmp1 = dconjg(zwork(nj,nk)) * dtmp1
      do ni=1,ndim
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
end subroutine calc_exp_hmat_proj
