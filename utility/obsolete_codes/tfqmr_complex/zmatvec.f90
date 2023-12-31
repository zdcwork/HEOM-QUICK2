subroutine zmatvec(zmat, nrho, nunk, zvec, iop)
implicit none
!
! transfer between matrix zmat of size (nrho,nrho,nunk) and 
!                  vector zvec of size (nunk*nrho**2)
! iop = 0 : zmat ==> zvec
!     = 1 : zvec ==> zmat
!
integer,   intent(in)     :: nrho, iop
integer*8, intent(in)     :: nunk
complex*16, intent(inout) :: zmat(nrho,nrho,*)
complex*16, intent(inout) :: zvec(*)
!
integer                   :: ni, nj, nk
integer*8                 :: lni, lnj
!
if (iop .eq. 0) then
   do lni=1,nunk
      lnj = (lni - 1) * nrho**2
      nk = 0
      do nj=1,nrho
         do ni=1,nrho
            nk = nk + 1
            zvec(lnj + nk) = zmat(ni,nj,lni)
         end do
      end do
   end do
else if (iop .eq. 1) then
   do lni=1,nunk
      lnj = (lni - 1) * nrho**2
      nk = 0
      do nj=1,nrho
         do ni=1,nrho
            nk = nk + 1
            zmat(ni,nj,lni) = zvec(lnj + nk)
         end do
      end do
   end do
else 
   write(6,*)'zmatvec: error! unknown iop = ', iop
   stop
end if
return
end subroutine zmatvec
