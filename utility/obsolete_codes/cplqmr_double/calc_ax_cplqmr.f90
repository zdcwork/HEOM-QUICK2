subroutine calc_ax_cplqmr(nrho, xqmr, bqmr)
!
! purpose : calculate the RHS of EOM for rho
!    
! input/output : rcgtmp, rcgrhs (refreshed)
!         
implicit none
!
integer,    intent(in)  :: nrho
complex*16, intent(in)  :: xqmr(nrho,nrho,*)
complex*16, intent(out) :: bqmr(nrho,nrho,*)
!
call calc_ax_tfqmr(nrho, xqmr, bqmr)
!
return
end subroutine calc_ax_cplqmr
