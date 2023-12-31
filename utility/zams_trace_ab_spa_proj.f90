subroutine zams_trace_ab_spa_proj(transb, korbs, kspin, nnz, b, irowb, icolb, ztrace)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! calculate and output trace(ams * op(b)), b in sparse COO format
! op(b) can be b, b^t (transpose of b)
!
character*1, intent(in)   :: transb
integer, intent(in)       :: korbs, kspin, nnz, irowb(*), icolb(*)
complex*16, intent(in)    :: b(*)
complex*16, intent(out)   :: ztrace
integer                   :: mi, mk, nk, ni
logical                   :: lsame
external                  :: lsame
complex*16                :: ctmp1
!
ztrace = czero
!
if ( lsame(transb, 'N') ) then
   do nk=1,nnz_ams(korbs,kspin)
      mi = row_ams(nk,korbs,kspin)
      mk = col_ams(nk,korbs,kspin)
      ctmp1 = cval_ams(nk,korbs,kspin)
      do ni=1,nnz
         if (irowb(ni) .ne. mk .or. icolb(ni) .ne. mi) cycle
         ztrace = ztrace + ctmp1 * b(ni)
      end do
   end do
else if ( lsame(transb, 'T') ) then
   do nk=1,nnz_ams(korbs,kspin)
      mi = row_ams(nk,korbs,kspin)
      mk = col_ams(nk,korbs,kspin)
      ctmp1 = cval_ams(nk,korbs,kspin)
      do ni=1,nnz
         if (irowb(ni) .ne. mi .or. icolb(ni) .ne. mk) cycle
         ztrace = ztrace + ctmp1 * b(ni)
      end do
   end do
else
   write(6,*)
   write(6,*)'zams_trace_ab_spa_proj: error! only transb = N or T allowed ', transb
   stop
end if
!
return
end subroutine zams_trace_ab_spa_proj
