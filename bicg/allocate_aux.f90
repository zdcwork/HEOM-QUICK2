subroutine allocate_aux
use auxmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer              :: istat, dim1
!
!-----------------------------
! Auxiliary memory for calculation of transpose(A) * X  - obsolete
!
dim1 = nrho**2
!
! ***name rule***  [labcde]
! (1) l : liouville, (2-3) ab : hs (Hs), an (a_ms), ac (a^dag_ms), bn (related to cb), dn (related to cd)         
! (4) c : r for real,   i for imaginary
! (5) d : n for normal, t for transpose
! (6) e : l for left,   r for right 
! 
allocate(lhsrnl(dim1,dim1), lhsrnr(dim1,dim1), lhsinl(dim1,dim1), lhsinr(dim1,dim1), STAT=istat)
!----------------------------
!
end subroutine allocate_aux
!----------------
!
subroutine free_aux
use auxmod
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(lhsrnl, lhsrnr, lhsinl, lhsinr, STAT=istat)
deallocate(lenams, lrowams, lcolams, lvalams, STAT=istat)
deallocate(lrtmp1, lctmp1, lvtmp1, STAT=istat)
!
end subroutine free_aux
