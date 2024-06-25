subroutine zams_trace_ab(transb, korbs, kspin, b, ldb, ztrace)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! calculate and output trace(ams * op(b))
! op(b) can be b, b^t (transpose of b)
!
character*1, intent(in)   :: transb
integer, intent(in)       :: korbs, kspin, ldb
complex*16, intent(in)    :: b(ldb,*)
complex*16, intent(out)   :: ztrace
integer                   :: mi, mk, nk
logical                   :: lsame
external                  :: lsame
!
ztrace = czero
!
if ( lsame(transb, 'N') ) then
   do nk=1,nams
      mi = rowams(nk,korbs,kspin)
      mk = colams(nk,korbs,kspin)
      if (sgnams(nk,korbs,kspin) .eq. 0) then
         ztrace = ztrace + b(mk,mi)
      else
         ztrace = ztrace - b(mk,mi)
      end if
   end do
else if ( lsame(transb, 'T') ) then
   do nk=1,nams
      mi = rowams(nk,korbs,kspin)
      mk = colams(nk,korbs,kspin)
      if (sgnams(nk,korbs,kspin) .eq. 0) then
         ztrace = ztrace + b(mi,mk)
      else
         ztrace = ztrace - b(mi,mk)
      end if
   end do
else
   write(6,*)
   write(6,*)'zams_trace_ab: error! only transb = N or T allowed ', transb
   stop
end if
!
return
end subroutine zams_trace_ab
