subroutine calcocc(occu, occd, sigmau, sigmad)
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,     intent(out) :: occu, occd
complex*16, intent(out) :: sigmau(*), sigmad(*)
integer                 :: ni, nj, nk, istat, nnz
integer*8               :: lni
complex*16              :: ctmp1
!
occu = 0.d0
occd = 0.d0
!
! population = tr( c^\dag_{ms} c_{ms} \rho )
!
if (.not. lsparse) then
   do nk=1,nspin
      do ni=1,norbs
         call zamsmm1('l', 'n', 'n', ni, nk, cunity, rho(1,1,1), nrho, czero, cmtmp1, nrho) 
         call zams_trace_ab('t', ni, nk, cmtmp1, nrho, ctmp1)
         if (nk .eq. 1) then
            sigmau(ni) = ctmp1
            occu = occu + dble(ctmp1)
         else
            sigmad(ni) = ctmp1
            occd = occd + dble(ctmp1)
         end if
      end do
   end do
else 
   nnz = nnz_spa(1)
   lni = ind_spa(1)
   call zmat_coo2dns(nnz, irow_spa(lni), icol_spa(lni), rho_spa(lni), cmtmp1, nrho, nrho, nrho)
   do nk=1,nspin
      do ni=1,norbs
         !if (lprams) then
         !    call zamsmm_coo_proj('l', 'n', 'n', ni, nk, cunity, rho_spa(lni), nnz, irow_spa(lni), icol_spa(lni), &
         !                         czero, cmtmp1, nrho, cmtmp2, nrho)                
         !    call zams_trace_ab_proj('t', ni, nk, cmtmp1, nrho, ctmp1)
         !else
         !    call zamsmm_coo('l', 'n', 'n', ni, nk, cunity, rho_spa(lni), nnz, irow_spa(lni), icol_spa(lni), &
         !                    czero, cmtmp1, nrho, cmtmp2, nrho)                
         !    call zams_trace_ab('t', ni, nk, cmtmp1, nrho, ctmp1)
         !end if
         call zhemm('l', 'u', nrho, nrho, cunity, znms(1,1,ni,nk), nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
         ctmp1 = czero
         do nj=1,nrho
            ctmp1 = ctmp1 + cmtmp2(nj,nj)
         end do
         if (nk .eq. 1) then
            sigmau(ni) = ctmp1
            occu = occu + dble(ctmp1)
         else
            sigmad(ni) = ctmp1
            occd = occd + dble(ctmp1)
         end if
      end do
   end do
end if
!
return
end subroutine calcocc
