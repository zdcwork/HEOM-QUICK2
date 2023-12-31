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
integer                 :: ni, nj
!
call calc_ax_tfqmr(nrho, xqmr, bqmr)
!
do ni=1,nrho
   do nj=ni+1,nrho
      bqmr(ni,nj,1) = dcmplx(0.d0, 0.d0)
   end do
   bqmr(ni,ni,1) = dcmplx(dble(bqmr(ni,ni,1)), 0.d0)
end do
!
end subroutine calc_ax_cplqmr
