subroutine checkrho_spa
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs, nnz
logical                 :: lcheck
integer                 :: isgn, iorbs, ispin, ialf, idraw
integer*8               :: lni, lnp
real*8                  :: dtmp1
real*8,     allocatable :: tmpval(:), rwork(:)
complex*16, allocatable :: work(:)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'checkrho_spa: error entry '
   stop
end if
!
lcheck = .true.
write(6,*)
write(6,*)'checkrho_spa: check reduced density matrix '
!
!dtmp1 = 0.d0
!do nj=1,nrho
!   do ni=1,nrho
!      dtmp1 = max(dtmp1, cdabs(rhosys(ni,nj)))
!   end do
!end do
!write(6,*)'checkrho_spa: max element in rho '
!write(6,1211)0, dtmp1
!call flush(6)
!
do ni=1,nrho
  if (dble(rhosys(ni,ni)) .lt. 0.d0) then
    write(6,*)'checkrho_spa: Warning, negative diag. ', ni, dble(rhosys(ni,ni))
  end if
end do
call flush(6)
!
! rho0 should be a positive definite matrix
!
cmtmp1(1:nrho,1:nrho) = rhosys(1:nrho,1:nrho)
!
lwork = nrho**2
allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
call zheev('N', 'u', nrho, cmtmp1, nrho, tmpval, work, lwork, rwork, info)
if (info .ne. 0) then
   write(6,*)
   write(6,*)'checkrho_spa: error! diagnalization failed in checkrho_spa ', info
   stop
end if
!
do ni=1,nrho
   if (tmpval(ni) .lt. 0.d0) then
      write(6,*)'checkrho_spa: Warning, negative eigen.', ni, tmpval(ni)
      lcheck = .false.
   end if
end do
call flush(6)
deallocate(work ,rwork, tmpval, STAT=istat)
!
write(6,*)
write(6,*)'checkrho_spa: positivite-definiteness check complete! '
if (lcheck) then
  write(6,*)'checkrho_spa: rho is found OK! '
else
  write(6,*)'checkrho_spa: Warning! rho is not positive definite! '
end if
call flush(6)
!
write(6,*)
write(6,*)'checkrho_spa: max element in rho and ADOs '
do itier=1,ntier0
   dtmp1 = 0.d0
   do lni=nfirst(itier), nlast(itier)
      lnp = ind_spa(lni)
      nnz = nnz_spa(lni)
      do ni=1,nnz
         dtmp1 = max(dtmp1, cdabs(rho_spa(lnp+ni-1)))
      end do
   end do
   write(6,1211)itier, dtmp1
end do
call flush(6)
1211 format(' MAX ',I2,1x,e20.12e3)
!
write(6,*)
write(6,*)'checkrho_spa: sum of abs_val of rho and ADOs '
do itier=1,ntier0
   dtmp1 = 0.d0
   do lni=nfirst(itier), nlast(itier)
      lnp = ind_spa(lni)
      nnz = nnz_spa(lni)
      do ni=1,nnz
         dtmp1 = dtmp1 + cdabs(rho_spa(lnp+ni-1))
      end do
   end do
   write(6,1212)itier, dtmp1
end do
call flush(6)
1212 format(' SUM ',I2,1x,e20.12e3)
!
return
end subroutine checkrho_spa
