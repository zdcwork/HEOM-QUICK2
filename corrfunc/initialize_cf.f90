subroutine initialize_cf(tt)
use matmod
use corrfuncmod
use sparsemod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8, intent(out) :: tt
integer             :: ni, nj, ialf, ispin
character*1         :: transa, transb
integer*8           :: lni, lnj, lnk
integer             :: nnz, nnz2
!
tt = 0.d0
!
! system Hamiltonian
!
call refresh_rhosys
call calchs(rhosys)
if (igroundsteady .eq. 0) then
   dhs(1:nrho,1:nrho) = czero
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d10, rhosys)
else if (igroundsteady .eq. 2) then
   call calcdhs(tnow, rhosys)
   do ialf=1,nalf
      do ispin=1,nspin
         call getengyshift(tnow, ialf, ispin, eleadtime(ialf,ispin))
         write(6,*)'initialize_cf: lead energy shift at time, ialf, ispin, eshift'
         write(6,*)'initialize_cf: ', tnow, ialf, ispin, eleadtime(ialf,ispin)*hbar
      end do
   end do
   call flush(6)
else
   write(6,*)'initialize_cf: error! unknown igroundsteady = ', igroundsteady
   stop
end if
!
! Bh^- = B * rho^- + rho^- * B^+ = (Bh^+)^dag
! Ba^- = (-i) * (B * rho^- - rho^- * B^+) = (Ba^+)^dag 
!
if (lsparse) then
   brhoh_spa(1:lunkcf_spa) = czero
   brhoa_spa(1:lunkcf_spa) = czero
   transa = 'c'
   transb = 'n'
   if (b_ams_isgn .eq. 2) then
      transa = 'n'
      transb = 'c'
   end if
   do lni=1,nunk
      if ( nnzcf_spa(lni) .eq. 0 .or. nnz_spa(lni) .eq. 0 ) cycle
!
      lnj  = ind_spa(lni)
      nnz  = nnz_spa(lni)  
      lnk  = indcf_spa(lni)
      nnz2 = nnzcf_spa(lni)
! Bh
      if (lproj) then
          call zamsmm_coo_proj('l', transa, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo_proj('r', transb, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      else
          call zamsmm_coo('l', transa, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo('r', transb, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      end if
      call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), brhoh_spa(lnk), cmtmp1, nrho, nrho, nrho)
! Ba
      if (lproj) then
          call zamsmm_coo_proj('l', transa, 'n', b_ams_iorbs, b_ams_ispin, -eye  , rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo_proj('r', transb, 'n', b_ams_iorbs, b_ams_ispin,  eye  , rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      else
          call zamsmm_coo('l', transa, 'n', b_ams_iorbs, b_ams_ispin, -eye  , rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo('r', transb, 'n', b_ams_iorbs, b_ams_ispin,  eye  , rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      end if
      call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), brhoa_spa(lnk), cmtmp1, nrho, nrho, nrho)
   end do
else
   brhoh(1:nrho,1:nrho,1:nunk) = czero
   brhoa(1:nrho,1:nrho,1:nunk) = czero
   do lni=1,nunk 
      if (b_ams_isgn .eq. 1) then   
   ! Bh
         call zamsmm1('l', 'c', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho,  czero, brhoh(1,1,lni), nrho)
         call zamsmm1('r', 'n', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho, cunity, brhoh(1,1,lni), nrho)
   ! Ba     
         call zamsmm1('l', 'c', 'n', b_ams_iorbs, b_ams_ispin, -eye,       &
                      rho(1,1,lni), nrho,  czero, brhoa(1,1,lni), nrho)
         call zamsmm1('r', 'n', 'n', b_ams_iorbs, b_ams_ispin,  eye,       &
                      rho(1,1,lni), nrho, cunity, brhoa(1,1,lni), nrho)
      else if (b_ams_isgn .eq. 2) then
   ! Bh
         call zamsmm1('l', 'n', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho,  czero, brhoh(1,1,lni), nrho) 
         call zamsmm1('r', 'c', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho, cunity, brhoh(1,1,lni), nrho)
   ! Ba
         call zamsmm1('l', 'n', 'n', b_ams_iorbs, b_ams_ispin, -eye,       &
                      rho(1,1,lni), nrho,  czero, brhoa(1,1,lni), nrho)
         call zamsmm1('r', 'c', 'n', b_ams_iorbs, b_ams_ispin,  eye,       &
                      rho(1,1,lni), nrho, cunity, brhoa(1,1,lni), nrho)
      else 
         write(6,*)'initialize_cf: error! unknown b_ams_isgn ', b_ams_isgn
         stop
      end if
   end do
end if
!
if (.not. ljw_dos) goto 101
!
! Bdh^- = B^+ * rho^- + rho^- * B = (Bdh^+)^dag
! Bda^- = (-i) * (B^+ * rho^- - rho^- * B) = (Bda^+)^dag 
!
if (lsparse) then
   bdrhoh_spa(1:lunkcf_spa) = czero
   bdrhoa_spa(1:lunkcf_spa) = czero
   transa = 'n'
   transb = 'c'
   if (b_ams_isgn .eq. 2) then
      transa = 'c'
      transb = 'n'
   end if
   do lni=1,nunk
      if ( nnz_spa(lni) .eq. 0 .or. nnzcf_spa(lni) .eq. 0 ) cycle
      lnj  = ind_spa(lni)
      nnz  = nnz_spa(lni)  
      lnk  = indcf_spa(lni)
      nnz2 = nnzcf_spa(lni)
! Bdh
      if (lproj) then
          call zamsmm_coo_proj('l', transa, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo_proj('r', transb, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      else
          call zamsmm_coo('l', transa, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo('r', transb, 'n', b_ams_iorbs, b_ams_ispin, cunity, rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      end if
      call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), bdrhoh_spa(lnk),          &
                            cmtmp1, nrho, nrho, nrho)
! Bda
      if (lproj) then
          call zamsmm_coo_proj('l', transa, 'n', b_ams_iorbs, b_ams_ispin, -eye  , rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo_proj('r', transb, 'n', b_ams_iorbs, b_ams_ispin,  eye  , rho_spa(lnj), nnz,  &
                               irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      else
          call zamsmm_coo('l', transa, 'n', b_ams_iorbs, b_ams_ispin, -eye  , rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj),  czero, cmtmp1, nrho, cmtmp2, nrho)
          call zamsmm_coo('r', transb, 'n', b_ams_iorbs, b_ams_ispin,  eye  , rho_spa(lnj), nnz,  &
                          irow_spa(lnj), icol_spa(lnj), cunity, cmtmp1, nrho, cmtmp2, nrho)
      end if
      call zmat_extract_coo(nnz2, irowcf_spa(lnk), icolcf_spa(lnk), bdrhoa_spa(lnk),          &
                            cmtmp1, nrho, nrho, nrho)
   end do
else 
   bdrhoh(1:nrho,1:nrho,1:nunk) = czero
   bdrhoa(1:nrho,1:nrho,1:nunk) = czero
   do lni=1,nunk
      if (b_ams_isgn .eq. 1) then
   ! Bdh
         call zamsmm1('l', 'n', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho,  czero, bdrhoh(1,1,lni), nrho)
         call zamsmm1('r', 'c', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho, cunity, bdrhoh(1,1,lni), nrho)
   ! Bda     
         call zamsmm1('l', 'n', 'n', b_ams_iorbs, b_ams_ispin, -eye,       &
                      rho(1,1,lni), nrho,  czero, bdrhoa(1,1,lni), nrho)
         call zamsmm1('r', 'c', 'n', b_ams_iorbs, b_ams_ispin,  eye,       &
                      rho(1,1,lni), nrho, cunity, bdrhoa(1,1,lni), nrho)
      else if (b_ams_isgn .eq. 2) then
   ! Bdh
         call zamsmm1('l', 'c', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho,  czero, bdrhoh(1,1,lni), nrho) 
         call zamsmm1('r', 'n', 'n', b_ams_iorbs, b_ams_ispin, cunity,     &
                      rho(1,1,lni), nrho, cunity, bdrhoh(1,1,lni), nrho)
   ! Bda
         call zamsmm1('l', 'c', 'n', b_ams_iorbs, b_ams_ispin, -eye,       &
                      rho(1,1,lni), nrho,  czero, bdrhoa(1,1,lni), nrho)
         call zamsmm1('r', 'n', 'n', b_ams_iorbs, b_ams_ispin,  eye,       &
                      rho(1,1,lni), nrho, cunity, bdrhoa(1,1,lni), nrho)
      else 
         write(6,*)'initialize_cf: error! unknown b_ams_isgn ', b_ams_isgn
         stop
      end if
   end do
end if
!
101 continue
!
! check if resume from a previous broken job
!
if (lresume_cf .and. icont_cf .eq. 1) then
   write(6,*)'initialize_cf: resume from a previous td job for dos '
   call resume_cf(tt, 1)
end if
!
if (a_ams_isgn .eq. 1) then
   do nj=1,nrho
      do ni=1,nrho
         a_ams(ni,nj) = dconjg(zamsall(nj,ni,a_ams_iorbs,a_ams_ispin))
!         a_ams(ni,nj) = dconjg(amsori(nj,ni,a_ams_iorbs,a_ams_ispin) * cunity)
      end do
   end do
else if (a_ams_isgn .eq. 2) then
   a_ams(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,a_ams_iorbs,a_ams_ispin)
!   a_ams(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,a_ams_iorbs,a_ams_ispin) * cunity
else 
   write(6,*)'initialize_cf: error! unknown a_ams_isgn ', a_ams_isgn
   stop
end if
!
if (b_ams_isgn .eq. 1) then
   do nj=1,nrho
      do ni=1,nrho
         b_ams(ni,nj) = dconjg(zamsall(nj,ni,b_ams_iorbs,b_ams_ispin))
!         b_ams(ni,nj) = dconjg(amsori(nj,ni,b_ams_iorbs,b_ams_ispin) * cunity)
      end do
   end do
   b_ams_dag(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,b_ams_iorbs,b_ams_ispin)
!   b_ams_dag(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,b_ams_iorbs,b_ams_ispin) * cunity
else if (b_ams_isgn .eq. 2) then
   b_ams(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,b_ams_iorbs,b_ams_ispin)
!   b_ams(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,b_ams_iorbs,b_ams_ispin) * cunity
   do nj=1,nrho
      do ni=1,nrho
         b_ams_dag(ni,nj) = dconjg(zamsall(nj,ni,b_ams_iorbs,b_ams_ispin))
!         b_ams_dag(ni,nj) = dconjg(amsori(nj,ni,b_ams_iorbs,b_ams_ispin) * cunity)
      end do
   end do
else 
   write(6,*)'initialize_cf: error! unknown b_ams_isgn ', b_ams_isgn
   stop
end if
!
end subroutine initialize_cf
