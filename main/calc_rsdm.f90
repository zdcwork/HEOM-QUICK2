subroutine calc_rsdm(rhoinp, nrho0, norbs0, nspin0, rsdmout)
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer,    intent(in)  :: nrho0, norbs0, nspin0
complex*16, intent(in)  :: rhoinp(nrho0,*)
complex*16, intent(out) :: rsdmout(norbs0,norbs0,*)
integer                 :: iorbs, ispin, iorbs2, ni, nj
complex*16              :: ctmp1
!
rsdmout(1:norbs0,1:norbs0,1:nspin0) = czero
!
! P_ij = tr(a_j^dag a_i rho)
!
do ispin=1,nspin0
   do iorbs2=1,norbs0
      call zamsmm1('l', 'n', 'n', iorbs2, ispin, cunity, rhoinp(1,1), nrho0, czero, cmtmp1, nrho0)
      do iorbs=1,norbs0
         call zamsmm1('l', 'c', 'n', iorbs, ispin, cunity, cmtmp1, nrho0, czero, cmtmp2, nrho0)
         ctmp1 = czero
         do ni=1,nrho0 
            ctmp1 = ctmp1 + cmtmp2(ni,ni)
         end do
         rsdmout(iorbs2,iorbs,ispin) = ctmp1  
      end do
   end do
end do
!
! hermitianize
!
do ispin=1,nspin
   do ni=1,norbs
      do nj=1,ni
         ctmp1 = (rsdmout(ni,nj,ispin) + dconjg(rsdmout(nj,ni,ispin))) * 5.d-1
         rsdmout(ni,nj,ispin) = ctmp1
         rsdmout(nj,ni,ispin) = dconjg(ctmp1)
      end do
   end do
end do
!
end subroutine calc_rsdm
