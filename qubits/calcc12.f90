subroutine calcc12(dout, rhoinp)
use matmod
use tmpmatmod
implicit none
!
include '../include/sizes'
include '../include/common'
!
! Calculate concurrence of two dots (norbs=2, nspin=2) C12
! The calculation follows Eq.(S2) of SM of RPL 111, 246807 (2013)
!
real*8, intent(out) :: dout
complex*16, intent(in) :: rhoinp(nrho,*)
integer :: ni
real*8  :: s12, dtmp1
!
if (.not. lspin3d .or. norbs .ne. 2 .or. nspin .ne. 2) then
   write(6,*)'calcc12: error entrance, lspin3d, norbs, nspin ', &
             lspin3d, norbs, nspin
   stop
end if
!
call calcspinproduct(1, 2, s12, rhoinp)
!
! N1 = n1u + n1d - 2 n1u*n1d
!
cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,1,1), 0.d0) ! a1u
call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
           czero, cmtmp2, nrho) ! n1u
cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,1,2), 0.d0) ! a1d
call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
           czero, cmtmp3, nrho) ! n1d
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp3, nrho, &
           czero, cmtmp1, nrho) ! n1u * n1d 
cmtmp4(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho) + cmtmp3(1:nrho,1:nrho) - &
                        2.d0 * cmtmp1(1:nrho,1:nrho)
!
! N2 = n2u + n2d - 2 n2u*n2d
!
cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,1), 0.d0) ! a2u
call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
           czero, cmtmp2, nrho) ! n2u
cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,2), 0.d0) ! a2d
call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
           czero, cmtmp3, nrho) ! n2d
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp3, nrho, &
           czero, cmtmp1, nrho) ! n2u * n2d 
cmtmp5(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho) + cmtmp3(1:nrho,1:nrho) - &
                        2.d0 * cmtmp1(1:nrho,1:nrho)
! N1 * N2
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp5, nrho, &
           czero, cmtmp1, nrho)
! <N1 * N2>
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, rhoinp, nrho, &
           czero, cmtmp2, nrho)
dtmp1 = 0.d0
do ni=1,nrho
   dtmp1 = dtmp1 + dble(cmtmp2(ni,ni))
end do
!
dout = max(0.d0, -5.d-1 - 2.d0 * s12 / dtmp1)
!
return
end subroutine calcc12
