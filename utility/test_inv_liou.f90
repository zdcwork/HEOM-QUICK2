subroutine test_inv_liou(hinp)
use random
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in)  :: hinp(nrho,*)
integer                 :: ni, nj, nk, istat, nrho2, info, nwork
logical                 :: lherm
real*8                  :: dtmp1
complex*16              :: gamma_j, ctmp1
real*8,     allocatable :: tmpval(:), rwork(:)
complex*16, allocatable :: work(:)
complex*16, allocatable :: chmat1(:,:), chmat2(:,:), chmat3(:,:)
complex*16, allocatable :: clmat1(:,:), clmat2(:,:), clmat3(:,:)
!
nrho2 = nrho**2
allocate(chmat1(nrho,nrho), chmat2(nrho,nrho), chmat3(nrho,nrho), STAT=istat)
allocate(clmat1(nrho2,nrho2), clmat2(nrho2,nrho2), clmat3(nrho2,nrho2), STAT=istat)
chmat1(1:nrho,1:nrho)   = czero
chmat2(1:nrho,1:nrho)   = czero
chmat3(1:nrho,1:nrho)   = czero
clmat1(1:nrho2,1:nrho2) = czero
clmat2(1:nrho2,1:nrho2) = czero
clmat3(1:nrho2,1:nrho2) = czero
!
call checkhermicity(hinp, nrho, chmat1, nrho, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)
   write(6,*)'test_inv_liou: error! hinp is not Hermitian '
   stop
end if
!
gamma_j = dcmplx(random_normal(), random_normal())
do nj=1,nrho
   do ni=1,nrho
      chmat1(ni,nj) = dcmplx(random_normal(), random_normal())
   end do
end do
!
! [hinp, *] --> clmat1
!
call zhil2liou('l', 'n', nrho, hinp, nrho, clmat1, nrho2)
call zhil2liou('r', 'n', nrho, hinp, nrho, clmat2, nrho2)
clmat1(1:nrho2,1:nrho2) =  clmat1(1:nrho2,1:nrho2) - clmat2(1:nrho2,1:nrho2)
call checkhermicity(clmat1, nrho2, clmat2, nrho2, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)
   write(6,*)'test_inv_liou: error! Liouvillian of [hinp,*] is not Hermitian '
   stop
end if
!
! diagonalize clmat1, and then obtain 
! (-i*clmat1 + gamma_j)^(-1) = i * (clmat1 + i * gamma_j)^(-1)
!
nwork = nrho2**2
allocate(work(nwork), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
call zheev('V', 'u', nrho2, clmat1, nrho2, tmpval, work, nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)
   write(6,*)'test_inv_liou: error! diagonalization failed 1 ', info
   stop
end if
clmat2(1:nrho2,1:nrho2) = clmat1(1:nrho2,1:nrho2)
do nj=1,nrho2
   do ni=1,nrho2
      clmat1(ni,nj) = clmat1(ni,nj) / (tmpval(nj) + eye * gamma_j)
   end do
end do
call zgemm('n', 'c', nrho2, nrho2, nrho2, eye, clmat1, nrho2, clmat2, nrho2, czero, clmat3, nrho2)
deallocate(tmpval, work, rwork, STAT=istat)
!
call zgemv('n', nrho2, nrho2, cunity, clmat3, nrho2, chmat1(1,1), 1, czero, chmat2(1,1), 1)
!
! verify -i[hinp, chmat2] + gamma_j * chmat2 = chmat1 
! 
chmat3(1:nrho,1:nrho) = gamma_j * chmat2(1:nrho,1:nrho) - chmat1(1:nrho,1:nrho)
call zhemm('l', 'u', nrho, nrho, -eye, hinp, nrho, chmat2, nrho, cunity, chmat3, nrho)
call zgemm('n', 'n', nrho, nrho, nrho, eye, chmat2, nrho, hinp, nrho, cunity, chmat3, nrho)
call cmaxmat(nrho, nrho, chmat3, nrho, dtmp1)
write(6,*)
if ( dabs(dtmp1) .lt. dpico ) then
   write(6,*)'test_inv_liou: test passed '
else
   write(6,*)'test_inv_liou: test failed, max deviation ', dtmp1
end if
call flush(6)
!
deallocate(chmat1, chmat2, chmat3, STAT=istat)
deallocate(clmat1, clmat2, clmat3, STAT=istat)
!
return
end subroutine test_inv_liou
