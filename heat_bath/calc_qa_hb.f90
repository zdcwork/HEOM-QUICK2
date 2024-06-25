subroutine calc_qa_hb
use matmod
use tmpmatmod
use hbmod
!
implicit none
include '../include/sizes'
include '../include/common'
!
logical :: lherm
integer :: imode, iorbs, ispin
real*8  :: dtmp1
!
if (.not. lhb) then
   write(6,*)'calc_qa_hb: error! wrong entry ', lhb
   stop
end if
!
! Q_a = sum_s a_s^+ a_s  
!
! For the present case, the i-th mode corresponds to electron occupation operator for i-th system level.
! In future, Q_a may extend to inclue hopping term between i-th and j-th system levels.
! For the latter case, nmodes_nb != norbs (rather, nmodes_nb = norbs**2)
!
do imode=1,nmode_hb
   iorbs = imode
   cmtmp1(1:nrho,1:nrho) = czero
   do ispin=1,nspin
      cmtmp2(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,iorbs,ispin) * cunity
      call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, cunity, cmtmp1, nrho)
   end do
   qa_hb(1:nrho,1:nrho,imode) = cmtmp1(1:nrho,1:nrho)
!
   call checkhermicity(qa_hb(1,1,imode), nrho, cmtmp1, nrho, lherm, dtmp1)
   if (.not. lherm) then
      write(6,*)'calc_qa_hb: error! qa_hb not hermitian for mode ', imode
      stop
   end if
end do

end subroutine calc_qa_hb
