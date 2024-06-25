subroutine calcj(jleftu, jrightu, jleftd, jrightd)
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,  intent(out) :: jleftu, jrightu, jleftd, jrightd
integer*8            :: lni, lnj
integer              :: ni, nj, nk, iorbs, ispin, ialf, imats, isgn, nnz
real*8               :: dtmp1
complex*16           :: ctmp1, ctmp2
!
! J_L(t) = i \sum_{n,m,s} tr( c_{ms}\rho^{+,1}_{n,m,s,L} - \rho^{-,1}_{n,m,s,L} c^\dag_{ms})
!        = -2 Im \sum_{m,s,n} tr ( c_{ms}\rho^{+,1}_{n,m,s,L} )
!
jt = 0.d0
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
            jt(ialf,ispin) = jt(ialf,ispin) - 2.d0 * dimag(ctmp1)
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
            jt(ialf,ispin) = jt(ialf,ispin) - 2.d0 * dimag(ctmp2)
         end do
      end do
   end do
end if
606 continue
!
jt = jt * jconst
!
jleads = 0.d0
do ialf=1,nalf
   do ispin=1,nspin
      jleads(ialf) = jleads(ialf) + jt(ialf,ispin)
   end do
end do
!
jleftu  = jt(1,1)
if (nalf .ne. 1) jrightu = jt(2,1)
if (nspin .eq. 1) then
  jleftd  = 0.d0
  if (nalf .ne. 1) jrightd = 0.d0
else
  jleftd  = jt(1,2)
  if (nalf .ne. 1) jrightd = jt(2,2)
end if
!
end subroutine calcj
