subroutine calc_ax_tfqmr(nrho, nunk, zvecin, zvecout)
use tfqmrmod
implicit none
!
integer, intent (in)     :: nrho
integer*8, intent (in)   :: nunk
complex*16, intent (in)  :: zvecin(*)
complex*16, intent (out) :: zvecout(*)
!
call zmatvec(xqmr(1,1,1), nrho, nunk, zvecin(1), 1)
!
call calc_ax_tfqmr_comp(nrho, xqmr, bqmr)
!
call zmatvec(bqmr(1,1,1), nrho, nunk, zvecout(1), 0)
!
end subroutine calc_ax_tfqmr
