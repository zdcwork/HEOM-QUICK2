subroutine liouzhemm(ialr, irhont, zalpha, nrho, zamat, lamat, zrho, lrho, &
                     const, zout, lout, zwork, lwork)
implicit none
!
! purpose : calculation in the Liouville space associated with
!           Y=A*rho or Y=rho*A or Y=A*rho^\dag or Y=rho^\dag*A. However, the 
!           Liouville coefficient matrix L_A should be transposed. 
!           A and rho are complex.
!
!           zout = const * zout + L^t_zamat * zrho * alpha
!
!           where zamat is a Hermitian matrix in Hilbert space (A^dag = A)
!
!           ialr   = 1 : A to the left of rho
!                  = 2 : A to the right of rho
!          
!           irhont = 0 : rho is in normal form
!                 != 0 : rho is in its Hermitian conjugated form
!
integer, intent(in)       :: ialr, irhont, lamat, nrho, lrho, lout, lwork
complex*16, intent(in)    :: const, zalpha
complex*16, intent(in)    :: zamat(lamat,*), zrho(lrho,*)
complex*16, intent(inout) :: zout(lout,*), zwork(lwork,*)
integer                   :: ni, nj
complex*16                :: cunity
!
cunity = dcmplx(1.d0, 0.d0)
!
if (lwork .lt. lrho) then
   write(6,*)'liouzhemm: error! lwork < lrho ', lwork, lrho
   stop
end if
!
if (ialr .eq. 1 .and. irhont .eq. 0) then  ! LY = LA * Lrho
! Y = A^dag * rho
   call zhemm('l', 'l', nrho, nrho, zalpha, zamat, lamat, zrho, lrho,        &
              const, zout, lout)
   return
!
else if (ialr .eq. 1 .and. irhont .ne. 0) then  ! LY = LA * Lrho^dag
! Y = rho^dag * A
   do nj=1,nrho
      do ni=1,nrho
         zwork(ni,nj) = dconjg(zrho(nj,ni))
      end do
   end do
   call zhemm('r', 'l', nrho, nrho, zalpha, zamat, lamat, zwork, lwork,      &
              const, zout, lout)
   return
!
else if (ialr .eq. 2 .and. irhont .eq. 0) then  ! LY = Lrho * LA
! Y = rho * A^dag
   call zhemm('r', 'l', nrho, nrho, zalpha, zamat, lamat, zrho, lrho,        &
              const, zout, lout)
   return
!
else if (ialr .eq. 2 .and. irhont .ne. 0) then  ! LY = Lrho^dag * LA
! Y = A * rho^dag
   do nj=1,nrho
      do ni=1,nrho
         zwork(ni,nj) = dconjg(zrho(nj,ni))
      end do
   end do
   call zhemm('l', 'l', nrho, nrho, zalpha, zamat, lamat, zwork, lwork,      &
              const, zout, lout)
   return
!
else
   write(6,*)
   write(6,*)' error in liouzhemm, unknown ialr= ', ialr
   stop
end if
!
end subroutine liouzhemm
