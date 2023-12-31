subroutine zams_trace_ab_proj(transb, korbs, kspin, b, ldb, ztrace)
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
   do nk=1,nnz_ams(korbs,kspin)
      mi = row_ams(nk,korbs,kspin)
      mk = col_ams(nk,korbs,kspin)
      ztrace = ztrace + cval_ams(nk,korbs,kspin) * b(mk,mi)
   end do
else if ( lsame(transb, 'T') ) then
   do nk=1,nnz_ams(korbs,kspin)
      mi = row_ams(nk,korbs,kspin)
      mk = col_ams(nk,korbs,kspin)
      ztrace = ztrace + cval_ams(nk,korbs,kspin) * b(mi,mk)
   end do
else
   write(6,*)
   write(6,*)'zams_trace_ab_proj: error! only transb = N or T allowed ', transb
   stop
end if
!
return
end subroutine zams_trace_ab_proj
