subroutine zpout_psd(mpole, msgn, zeig, zcoef, zpole)
!
! purpose : find poles and coefficients for PSD scheme
! 
! f(z) = 0.5 + sum_i zcoef * 2 * z / (z^2 - zpole^2)
! with z = beta * w
!
implicit none
!
integer, intent (in)     :: mpole, msgn
complex*16, intent (out)  :: zeig(*), zcoef(*), zpole(mpole,*)
integer              :: istat, info, ni, nj, nk
integer              :: dim0, dim1, n0, n1
integer              :: lwork
real*8               :: dtmp0
real*8, allocatable     :: delta(:,:), dmat0(:,:), dmat1(:,:)
real*8, allocatable     :: dval0(:), dval1(:), dwork(:)
real*8, allocatable     :: eval0(:), eval1(:), mattmp1(:,:)
real*8, allocatable     :: tmpa(:), dcoefeta(:), mattmp0(:,:)
complex*16            :: cunity
!
n0   = mpole
n1   = mpole - 1
dim0 = 2 * mpole
dim1 = dim0 - 1
cunity = (1.d0, 0.d0)
!
lwork = 5 * dim0
allocate(delta(dim0,dim0), dmat0(dim0,dim0), dmat1(dim1,dim1), mattmp0(dim0,dim0), mattmp1(dim1,dim1), STAT=istat)
allocate(dval0(dim0), dval1(dim1), dwork(lwork), STAT=istat)
!
! Kronecker delta function
!
delta = 0.d0
do ni=1,dim0
   delta(ni,ni) = 1.d0
end do
!
! construct real symmetric matrices
!
dmat0 = 0.d0
do nj=1,dim0
   do ni=1,dim0
      if (ni .eq. nj + 1 .or. ni .eq. nj - 1) then
         dtmp0 = (2.d0 * dble(ni) - 1.d0) * (2.d0 * dble(nj) - 1.d0)
         dmat0(ni,nj) = dmat0(ni,nj) + 1.d0 / dsqrt(dtmp0)
      end if
   end do
end do
!
dmat1 = 0.d0
do nj=1,dim1
   do ni=1,dim1
      if (ni .eq. nj + 1 .or. ni .eq. nj - 1) then
         dtmp0 = (2.d0 * dble(ni) + 1.d0) * (2.d0 * dble(nj) + 1.d0)
         dmat1(ni,nj) = dmat1(ni,nj) + 1.d0 / dsqrt(dtmp0)
      end if
   end do
end do
!write(6,*)'lambda'
!do ni = 1, dim0
!      write(6,10)dmat0(ni,1),dmat0(ni,2),dmat0(ni,3),dmat0(ni,4)
!end do
!write(6,*)
!write(6,*)'tilde lambda'
!do ni = 1, dim1
!      write(6,11)dmat1(ni,1),dmat1(ni,2),dmat1(ni,3)
!end do
!write(6,*)
!10 format(4(f9.5))
!11 format(3(f9.5))
!
! eigenvalues 
!
call dsyev('N', 'U', dim0, dmat0, dim0, dval0, dwork, lwork, info)
if (info .ne. 0) then
   write(6,*)'zpout_psd: error! diagonalization failed for dmat0', info
   stop
end if
call dsyev('N', 'U', dim1, dmat1, dim1, dval1, dwork, lwork, info)
if (info .ne. 0) then
   write(6,*)'zpout_psd: error! diagonalization failed for dmat1', info
   stop
end if
!
allocate(eval0(n0), eval1(n1), STAT=istat)
allocate(tmpa(n0), dcoefeta(n0), STAT=istat)
!
eval0(1:n0) = dval0(1:n0)
eval1(1:n1) = dval1(1:n1)
do ni=1,n0
   eval0(ni) = 4.d0 / eval0(ni)**2  ! eigXi
end do
do ni=1,n1
   eval1(ni) = 4.d0 / eval1(ni)**2  ! eigZeta
end do
!
zeig(1:n0) = dcmplx(eval0(1:n0), 0.d0)
!
tmpa = 1.d0
do nj=1,n0
   do nk=1,nj-1
      tmpa(nj) = tmpa(nj) * (eval1(nk) - eval0(nj)) / (eval0(nk) - eval0(nj)) 
   end do
   do nk=nj+1,n0
      tmpa(nj) = tmpa(nj) * (eval1(nk - 1) - eval0(nj)) / (eval0(nk) - eval0(nj))
   end do
end do

!
!write(6,*)'zpout_psd: show dcoefeta'
do ni=1,n0
   dcoefeta(ni) = 5.d-1 * dble(n0 * (2 * n0 + 1)) * tmpa(ni)
   !write(6,100) ni, dcoefeta(ni)
   if (isnan(dcoefeta(ni))) then
      write(6,*)'zpout_psd: error! NAN found for dcoefeta(', ni, ')'
      stop
   end if
end do
100 format('zpout_psd: dcoefeta(', I4, ') = ', e14.6e3)
!
zcoef(1:n0) = -cunity * dcoefeta(1:n0)
!
do ni=1,n0
   dtmp0 = dsqrt(eval0(ni))
   !zpole(ni) = dcmplx(0.d0, dabs(dtmp0))
   zpole(ni,1) = dcmplx(0.d0, dabs(dtmp0))
   zpole(ni,2) = dconjg(zpole(ni,1))
end do

!
deallocate(delta, dmat0, dmat1, STAT=istat)
deallocate(dval0, dval1, dwork, STAT=istat)
deallocate(eval0, eval1, STAT=istat)
deallocate(tmpa, dcoefeta, STAT=istat)
!
end subroutine
