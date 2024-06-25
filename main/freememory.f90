subroutine allocatememory
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat, dim1
!
allocate(ams(nrho,nrho), STAT=istat)
allocate(amsall(nrho,nrho,norbs,nspin), STAT=istat)
allocate(amsori(nrho,nrho,norbs,nspin), STAT=istat)
allocate(zamsall(nrho,nrho,norbs,nspin), STAT=istat)
allocate(znms(nrho,nrho,norbs,nspin), STAT=istat)
allocate(denmbr(nrho,nrho), denmbi(nrho,nrho), STAT=istat)
allocate(denmcr(nrho,nrho), denmci(nrho,nrho), STAT=istat)
allocate(jt(nalf,nspin), pocc(norbs,nspin), jleads(nalf), STAT=istat)
!
dim1 = nrho**2
!
! allocate auxiliary matrices (Liouville space)
!
allocate(dlmat1(dim1,dim1), dlmat2(dim1,dim1), dlmat3(dim1,dim1), dlmat4(dim1,dim1), &
         dlmat5(dim1,dim1), dlmat6(dim1,dim1),                                       &
         dlvec1(dim1), dlvec2(dim1), dlvec3(dim1), dlvec4(dim1),                     &
         dlvec5(dim1), dlvec6(dim1), dlvec7(dim1), dlvec8(dim1),                     &
         STAT=istat) 
!
end subroutine allocatememory
!
!------------
subroutine freememory
use matmod
use tmpmatmod
use sparsemod
use hbmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(indexcoef, itypecoef, ifactfcoef, ifactrcoef, STAT=istat)
deallocate(index_coef_ref, STAT=istat)
!
if (.not. lsparse) then
   deallocate(rho, STAT=istat)
end if
if (allocated(rhodiag)) deallocate(rhodiag, STAT=istat)
!
deallocate(hs, dhs, rsdm, rhosys, STAT=istat)
deallocate(indextable, STAT=istat)
deallocate(filtertable, STAT=istat)
deallocate(oprtable, STAT=istat)
deallocate(ifirst, ilast, STAT=istat)
deallocate(nfirst, nlast, STAT=istat)
deallocate(ioprfirst, ioprlast, STAT=istat)
deallocate(noprfirst, noprlast, STAT=istat)
deallocate(indexdraw, STAT=istat)
!
deallocate(eleadinfty, eleadtime, STAT=istat)
deallocate(amps, tchar, omegas, amps_dc, STAT=istat)
deallocate(cgamma, cgama, cb, cd, STAT=istat)
if (lscale) deallocate(dbsqrt, dbinvsqrt, STAT=istat)
deallocate(dmatsu, dbeta, dwp, gammap, dlwidth, STAT=istat)
if (allocated(dlw_sp)) deallocate(dlw_sp, STAT=istat)
!
if (nmats .gt. 0) then
  deallocate(dimpf, STAT=istat)
end if
deallocate(dxip, zxip, zlwmat, STAT=istat)
!
deallocate(iredex, dipm, ilead, dpm, STAT=istat)
deallocate(morbs, mspin, mmats, mpm, msopr, STAT=istat)
deallocate(mdrude, mfreq, mcor, mopr, mnfff, mmfff, STAT=istat)
deallocate(jorbs, jspin, jpm, jalf, jredex, jomit, STAT=istat)
deallocate(kmats, kpmats, dmemexp, STAT=istat)
!
deallocate(ams, STAT=istat)
deallocate(amsall, amsori, zamsall, znms, STAT=istat)
deallocate(rowams, colams, sgnams, valams, STAT=istat)
deallocate(denmbr, denmbi, denmcr, denmci, STAT=istat)
deallocate(jt, jleads, pocc, STAT=istat)
!
deallocate(iatmp1, iatmp2, STAT=istat)
deallocate(dmtmp1, dmtmp2, dmtmp3, dmtmp4, dmtmp5, dmtmp6, STAT=istat)
deallocate(cmtmp1, cmtmp2, cmtmp3, cmtmp4, cmtmp5, cmtmp6, STAT=istat) 
deallocate(cmtmp7, cmtmp8, cmtmp9, STAT=istat)
deallocate(cmtmp10, cmtmp11, cmtmp12, cmtmp13, cmtmp14, STAT=istat)
deallocate(cmtmp15, cmtmp16, cmtmp17, cmtmp18, STAT=istat)
deallocate(cmtmpa, cmtmpb, STAT=istat)
deallocate(zmtmp1, zmtmp2, STAT=istat)
if (allocated(zlmtmp1)) deallocate(zlmtmp1, STAT=istat)
if (allocated(zlmtmp2)) deallocate(zlmtmp2, STAT=istat)
if (allocated(zlmtmp3)) deallocate(zlmtmp3, STAT=istat)
if (allocated(zlmtmp4)) deallocate(zlmtmp4, STAT=istat)
deallocate(dmtmp1_sp, dmtmp2_sp, STAT=istat)
!
deallocate(dlmat1, dlmat2, dlmat3, dlmat4, dlmat5, dlmat6, STAT=istat)
deallocate(dlvec1, dlvec2, dlvec3, dlvec4, STAT=istat)
deallocate(dlvec5, dlvec6, dlvec7, dlvec8, STAT=istat)
!
! MFD
! 
if (mfdjob) then
  deallocate(wl0, wl1, wkln, dkln, STAT=istat)
  deallocate(spectral, focc, STAT=istat)
end if
!
! PFD and PSD
!
if (pfdjob .or. psdjob .or. psdlor .or. psdfff .or. psfjob) then
   deallocate(cpeigv, cpcoef, cppole, STAT=istat)
end if
!
! PSD_LOR
!
if (psdlor) then
   deallocate(dlor_coef, dlor_cent, dlor_width, STAT=istat)
end if
!
! PSD_FFF
!
if (psdfff) then
   deallocate(mpfff, nfff1d, mfff1d, afff, bfff, temfff, gfff, cfff, dfff, thefff, STAT=istat)
   deallocate(ctheta, STAT=istat)
end if
!
if (lband) then
   deallocate(band_coef, band_cent, band_width, STAT=istat)
end if
!
if (lspin3d) then
   deallocate(pauli, sopr, sdot, STAT=istat)
   if (lbfield3d) deallocate(dbf3d, STAT=istat)
end if
!
if (lsparse) then
   deallocate(ams_spa, ia_spa, ja_spa, STAT=istat)
   deallocate(irow_hsys, icol_hsys, cval_hsys, STAT=istat)
   deallocate(ind_spa, nnz_spa, STAT=istat)
   deallocate(irow_spa, icol_spa, rho_spa, STAT=istat)
end if
!
if (lssp) then
    deallocate(rhs_spa, STAT=istat)
end if
!
if (ltrun_der) then
   deallocate(dvalder, zmatder, STAT=istat)
end if
!
if (allocated(dvaljac)) deallocate(dvaljac, STAT=istat)
if (allocated(zmatjac)) deallocate(zmatjac, STAT=istat)
!
if (lhb) then
   deallocate(dwidth_hb, dcenter_hb, dcouple_hb, zpeta_hb, zppole_hb, zpeigv_hb, STAT=istat)
   deallocate(cb_hb, cd_hb, cgamma_hb, STAT=istat)
   deallocate(qa_hb, STAT=istat)
   deallocate(rho_hb, STAT=istat)
end if
!
if (lsimple) then
   deallocate(ioprindex, ioprref, STAT=istat)
end if
!
if (lad) then
    deallocate(ns2f, ms2f, mslow, icor_slow, dgama_drawer, STAT=istat)
    deallocate(zgefnl, zgefnr, zgelvl, zgelvr, STAT=istat)
    deallocate(zgefnl_sop, zgefnr_sop, zgelvl_sop, zgelvr_sop, STAT=istat)
    deallocate(nnz_zgf, irow_zgf, icol_zgf, cval_zgfl, cval_zgfr, STAT=istat)
    deallocate(nnz_sop, irow_sop, icol_sop, cval_sopl, cval_sopr, STAT=istat)
    deallocate(nnz_lvl, irow_lvl, icol_lvl, cval_lvl, STAT=istat)
    deallocate(nnz_lvr, irow_lvr, icol_lvr, cval_lvr, STAT=istat)
    deallocate(ind_lvl, ind_lvr, nnc_lvl, nnc_lvr, STAT=istat)
    deallocate(nnz_lvlsop, irow_lvlsop, icol_lvlsop, cval_lvlsop, STAT=istat)
    deallocate(nnz_lvrsop, irow_lvrsop, icol_lvrsop, cval_lvrsop, STAT=istat)
    deallocate(npermu, npermu_fact, STAT=istat)
    deallocate(icor_fast, STAT=istat)
    deallocate(nlvrow, nlvcol, nnrc2v, ncrc2v, STAT=istat)
end if
!
if (lproj) then
    deallocate(cproj, cproj_bar, STAT=istat)
    deallocate(nnz_ams, cval_ams, row_ams, col_ams, STAT=istat)
    deallocate(istate_proj, STAT=istat)
end if
!
end subroutine freememory
