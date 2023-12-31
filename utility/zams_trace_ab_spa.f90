subroutine zams_trace_ab_spa(transb, korbs, kspin, nnz, b, irowb, icolb, ztrace)
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
!
ztrace = czero
!
if ( lsame(transb, 'N') ) then
   do nk=1,nams
      mi = rowams(nk,korbs,kspin)
      mk = colams(nk,korbs,kspin)
      do ni=1,nnz
         if (irowb(ni) .ne. mk .or. icolb(ni) .ne. mi) cycle
         if (sgnams(nk,korbs,kspin) .eq. 0) then
            ztrace = ztrace + b(ni)
         else
            ztrace = ztrace - b(ni)
         end if
      end do
   end do
else if ( lsame(transb, 'T') ) then
   do nk=1,nams
      mi = rowams(nk,korbs,kspin)
      mk = colams(nk,korbs,kspin)
      do ni=1,nnz
         if (irowb(ni) .ne. mi .or. icolb(ni) .ne. mk) cycle
         if (sgnams(nk,korbs,kspin) .eq. 0) then
            ztrace = ztrace + b(ni)
         else
            ztrace = ztrace - b(ni)
         end if
      end do
   end do
else
   write(6,*)
   write(6,*)'zams_trace_ab_spa: error! only transb = N or T allowed ', transb
   stop
end if
!
return
end subroutine zams_trace_ab_spa
