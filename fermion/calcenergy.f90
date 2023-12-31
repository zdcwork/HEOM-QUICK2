subroutine calcenergy(dengy_sys, dengy_sb, dengy_tot)
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,  intent(out) :: dengy_sys, dengy_sb, dengy_tot
integer*8            :: lni, lnj
integer              :: ni, nj, nk, iorbs, ispin, ialf, imats, isgn, nnz
real*8               :: dtmp1
complex*16           :: ctmp1, ctmp2
!
! E_sys = tr(rho * H_sys)
!
dengy_sys = 0.d0
!
if (.not. lsparse) then
   cmtmp1(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
else 
   lni = ind_spa(1)
   call zmat_coo2dns(nnz_spa(1), irow_spa(lni), icol_spa(lni), rho_spa(lni), cmtmp1, nrho, nrho, nrho)
end if
call zgemm('n', 'n', nrho, nrho, nrho, cunity, hs, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
if (igroundsteady .gt. 1) then
    call zgemm('n', 'n', nrho, nrho, nrho, cunity, dhs, nrho, cmtmp1, nrho, cunity, cmtmp2, nrho)
end if
do ni=1,nrho
   dengy_sys = dengy_sys + dble(cmtmp2(ni,ni))
end do
!
! E_sb = \sum_{n,m,s} tr( c_{ms}\rho^{+,1}_{n,m,s,L} + \rho^{-,1}_{n,m,s,L} c^\dag_{ms})
!      = 2 Re \sum_{m,s,n} tr ( c_{ms}\rho^{+,1}_{n,m,s,L} )
!
dengy_sb = 0.d0
!
if (.not. lsparse) then
   do ispin=1,nspin
      do iorbs=1,norbs
         do ialf=1,nalf
            cmtmp2(1:nrho,1:nrho) = czero
            do imats=1,ncor
               if (lscale) then
                  dtmp1 = dbsqrt(iorbs,ispin,imats,ialf,1)
               else 
                  dtmp1 = 1.d0
               end if
               call look4drawer(1, ialf, iorbs, ispin, imats, ni)
               lni = nfirst(2) - 1 + indexdraw(ni)
               cmtmp2(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho) + rho(1:nrho,1:nrho,lni) * dtmp1
            end do
            call zams_trace_ab('n', iorbs, ispin, cmtmp2, nrho, ctmp1)
            dengy_sb = dengy_sb + 2.d0 * dble(ctmp1)
         end do
      end do
   end do
else
   do ispin=1,nspin
      do iorbs=1,norbs
         do ialf=1,nalf
            ctmp2 = czero
            do imats=1,ncor
               if (lscale) then
                  dtmp1 = dbsqrt(iorbs,ispin,imats,ialf,1)
               else 
                  dtmp1 = 1.d0
               end if
               call look4drawer(1, ialf, iorbs, ispin, imats, ni)
               lni = nfirst(2) - 1 + indexdraw(ni)
               nnz = nnz_spa(lni)
               lnj = ind_spa(lni)
               if (lprams) then
                   call zams_trace_ab_spa_proj('n', iorbs, ispin, nnz, rho_spa(lnj), irow_spa(lnj), icol_spa(lnj), ctmp1)
               else
                   call zams_trace_ab_spa('n', iorbs, ispin, nnz, rho_spa(lnj), irow_spa(lnj), icol_spa(lnj), ctmp1)
               end if
               ctmp2 = ctmp2 + ctmp1 * dtmp1
            end do
            dengy_sb = dengy_sb + 2.d0 * dble(ctmp2)
         end do
      end do
   end do
end if
606 continue
!
dengy_tot = dengy_sys + dengy_sb 
!
return
end subroutine calcenergy
