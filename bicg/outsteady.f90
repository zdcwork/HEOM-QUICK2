subroutine outsteady(iter)
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent (in)    :: iter
integer                 :: ni, nj, nk, nl, istat, ispin, ialf, iorbs, iorbs2, nnz
integer*8               :: lni, lnj
real*8                  :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8                  :: jleftu, jrightu, jleftd, jrightd, jleft, jright
real*8                  :: occu, occd, occ
real*8                  :: dunitt, dunitj
real*8                  :: dspin2(3)
complex*16, allocatable :: sigmau(:), sigmad(:)
integer,    allocatable :: irow_tmp1(:), icol_tmp1(:)
complex*16, allocatable :: cval_tmp1(:)
complex*16, allocatable :: cmat1(:,:), cmat2(:,:), cmat3(:,:), cmat4(:,:)
complex*16, allocatable :: cmat5(:,:), cmat6(:,:)
namelist / thermo / thermopower
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), STAT=istat)
allocate(cmat3(nrho,nrho), cmat4(nrho,nrho), STAT=istat)
allocate(cmat5(nrho,nrho), cmat6(nrho,nrho), STAT=istat)
!
thermopower = .false.
rewind(5)
read(5, thermo, end=101)
101 continue
!
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
if (lsparse) then
   nk = nrho**2
   allocate(irow_tmp1(nk), icol_tmp1(nk), cval_tmp1(nk), STAT=istat)
end if
!
call calcj(jleftu, jrightu, jleftd, jrightd)
jleft  = jleftu  + jleftd
if (nalf .gt. 1) then
   jright = jrightu + jrightd
end if
call calcocc(occu, occd, sigmau, sigmad)
occ = occu + occd
!
write(6,*)
write(6,*)'outsteady: output steady state properties'
write(6,*)
write(6,*)' occupation, spin up = ', occu
if (nspin .ne. 1) then
  write(6,*)'                down = ', occd
end if
!
call calcocc_detail
write(6,*)
write(6,*)' occupancy on all orbitals '
do iorbs=1,norbs
   do ispin=1,nspin
      write(6,518)iorbs, ispin, pocc(iorbs,ispin)
   end do
end do
call flush(6)
518 format(' iorb= ', I2, ' ispin= ', I2, ' nele= ', f8.3)
!
write(6,*)
write(6,*)'outsteady: steady electric current: '
if (funits .eq. 0) then
   if (nalf .le. 2) then
      write(6,*)'              jleft  = ', jleft , ' nA '
      if (nalf .gt. 1) then
         write(6,*)'             jright  = ', jright, ' nA '
      end if
   else ! multi-leads
      do ialf=1,nalf
         write(6,*)' j(ialf = ', ialf, ' )= ', jleads(ialf), ' nA '
      end do
   end if
else if (funits .eq. 1) then
   if (nalf .le. 2) then
      write(6,*)'              jleft  = ', jleft , ' pA '
      if (nalf .gt. 1) then
         write(6,*)'             jright  = ', jright, ' pA '
      end if
   else 
      do ialf=1,nalf
         write(6,*)' j(ialf = ', ialf, ' )= ', jleads(ialf), ' pA '
      end do
   end if
else if (funits .eq. 2) then
   if (nalf .le. 2) then
      write(6,*)'              jleft  = ', jleft , ' a.u. A '
      if (nalf .gt. 1) then
         write(6,*)'             jright  = ', jright, ' a.u. A '
      end if
   else
      do ialf=1,nalf
         write(6,*)' j(ialf = ', ialf, ' )= ', jleads(ialf), ' a.u. A '
      end do
   end if
end if
!
if (nspin .ne. 1) then
  write(6,*)
  write(6,*)' spin current '
  if (funits .eq. 0) then
     if (nalf .le. 2) then
        write(6,*)'              jleftu = ', jleftu , ' nA '
        write(6,*)'              jleftd = ', jleftd , ' nA '
        if (nalf .gt. 1) then
           write(6,*)'             jrightu = ', jrightu, ' nA '
           write(6,*)'             jrightd = ', jrightd, ' nA '
        end if
     else 
        do ialf=1,nalf
           write(6,*)' j_up  (ialf = ', ialf, ' )= ', jt(ialf,1), ' nA '
           write(6,*)' j_down(ialf = ', ialf, ' )= ', jt(ialf,2), ' nA '
        end do
     end if
  else if (funits .eq. 1) then
     if (nalf .le. 2) then
        write(6,*)'              jleftu = ', jleftu, ' pA '
        write(6,*)'              jleftd = ', jleftd, ' pA '
        if (nalf .gt. 1) then
           write(6,*)'             jrightu = ', jrightu, ' pA '
           write(6,*)'             jrightd = ', jrightd, ' pA '
        end if
     else  
        do ialf=1,nalf
           write(6,*)' j_up  (ialf = ', ialf, ' )= ', jt(ialf,1), ' pA '
           write(6,*)' j_down(ialf = ', ialf, ' )= ', jt(ialf,2), ' pA '
        end do
     end if
  end if
  else if (funits .eq. 2) then
   if (nalf .le. 2) then
      write(6,*)'              jleftu = ', jleftu, ' a.u. A '
      write(6,*)'              jleftd = ', jleftd, ' a.u. A'
      if (nalf .gt. 1) then
         write(6,*)'             jrightu = ', jrightu, ' a.u. A'
         write(6,*)'             jrightd = ', jrightd, ' a.u. A'
      end if
   else
      do ialf=1,nalf
         write(6,*)' j_up  (ialf = ', ialf, ' )= ', jt(ialf,1), ' a.u. A'
         write(6,*)' j_down(ialf = ', ialf, ' )= ', jt(ialf,2), ' a.u. A'
      end do
   end if
end if
call flush(6)
!
! output for plotting iv curve
!
!if (norbs .eq. 1 .and. nspin .eq. 2) then
if (nspin .eq. 2) then
    write(6,*)
    if (nalf .eq. 1) then
      write(6,222)iter, eleadinfty(1,1)*hbar, occu, occd, jleft
    else if (nalf .eq. 2) then
      write(6,222)iter, eleadinfty(1,1)*hbar, eleadinfty(2,1)*hbar, occu, occd, jleft, jright
    else ! multi-leads
      write(6,222)iter, (eleadinfty(ialf,1)*hbar, ialf=1,nalf), occu, occd, (jleads(ialf), ialf=1,nalf)
    end if
    call flush(6)
end if
222 format(' IVOUTPUT ', I5, 2x, 30(e15.6e3, 2x))
!
if (norbs .eq. 1 .and. nspin .eq. 1) then
    write(6,*)
    if (nalf .eq. 1) then
      write(6,222)iter, epara1, occ, jleft
    else if (nalf .eq. 2) then
      write(6,222)iter, epara1, occ, jleft, jright
    else ! multi-leads
      write(6,222)iter, epara1, occ, (jleads(ialf), ialf=1,nalf) 
    end if
    call flush(6)
end if
!
! energy
!
call calcenergy(dtmp1, dtmp2, dtmp3)
write(6,*)
write(6,*)'outsteady: Esys = ', dtmp1*hbar
write(6,*)'           Esb  = ', dtmp2*hbar
write(6,*)'           Etot = ', dtmp3*hbar
call flush(6)
!
! impurity susceptibility (2d spin)
!
if (nspin .gt. 1 .and. .not. lspin3d) then
   cmat1 = czero
   cmat2 = czero
   do iorbs=1,norbs
      cmat3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 1), 0.d0)
      call zamsmm1('l', 'c', 'n', iorbs, 1, cunity, cmat3, nrho, cunity, cmat1, nrho)
      cmat4(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 2), 0.d0)
      call zamsmm1('l', 'c', 'n', iorbs, 2, cunity, cmat4, nrho, cunity, cmat2, nrho)
   end do
   cmat1 = cmat1 - cmat2          ! n_up - n_down
   cmat3 = cmat1
   if (lproj) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmat1, nrho, czero, cmat2, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmat2, nrho, czero, cmat1, nrho)
   end if
   if (lsparse) then
      nnz = nnz_spa(1)
      lni = ind_spa(1)
      call zmat_zcoomm('r', 'n', 'n', nrho, dcmplx(0.5d0, 0.d0), rho_spa(lni), irow_spa(lni),    &
                       icol_spa(lni), nnz, cmat1, nrho, czero, cmat2, nrho,                      &
                       cval_tmp1, irow_tmp1, icol_tmp1, cmat5, nrho, cmat6, nrho)
   else
      call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.5d0, 0.d0), cmat1, nrho, rho(1,1,1), nrho, &
                 czero, cmat2, nrho)   
   end if
   dtmp1 = 0.d0
   do ni=1,nrho
      dtmp1 = dtmp1 + dble(cmat2(ni,ni))
   end do
   write(6,*)
   write(6,*)'outsteady: < Sz > = ', dtmp1
   dtmp2 = dtmp1
!
   cmat1 = cmat3
   call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.25d0, 0.d0), cmat1, nrho, cmat1, nrho,  &
              czero, cmat2, nrho)   ! (n_up - n_down)^2 / 4
   if (lproj) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmat2, nrho, czero, cmat1, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmat1, nrho, czero, cmat2, nrho)
   end if
   if (lsparse) then
      nnz = nnz_spa(1)
      lni = ind_spa(1)
      call zmat_zcoomm('r', 'n', 'n', nrho, cunity, rho_spa(lni), irow_spa(lni), icol_spa(lni), &
                       nnz, cmat2, nrho, czero, cmat1, nrho,                                    &
                       cval_tmp1, irow_tmp1, icol_tmp1, cmat5, nrho, cmat6, nrho)
   else
      call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmat2, nrho, rho(1,1,1), nrho,            &
                 czero, cmat1, nrho)
   end if
   dtmp1 = 0.d0
   do ni=1,nrho
      dtmp1 = dtmp1 + dble(cmat1(ni,ni))
   end do
   write(6,*)'outsteady: <Sz^2> = ', dtmp1
   call flush(6)
!
   if (lbfield) then
      write(6,*)
      write(6,*)'outsteady: calculation of local magnetic susceptibility '
      if (dabs(dbfield) .le. dpico) then
         write(6,*)'outsteady: dbfield too small, abort ...'
      else
!         write(6,*)'outsteady: chi_B = <Sz> / dbfield = ', dtmp2 / (dbfield * hbar)
         write(6,*)'outsteady: chi_B * T(1)           = ', dtmp2 / (dbfield * hbar) * dinvbeta(1)
      end if
      call flush(6)
   end if  
end if
!
! 3d spin 
!
if (lspin3d) then
   write(6,*)
   write(6,*)'outsteady: spin moment of each orbital     '
   write(6,*)'outsteady: iorbs, <Sx_i>, <Sy_i>, <Sz_i>   '
   if (lsparse) then
      nnz = nnz_spa(1)
      lni = ind_spa(1)
      call zmat_coo2dns(nnz, irow_spa(lni), icol_spa(lni), rho_spa(lni), cmat1, nrho, nrho, nrho)
   else 
      cmat1(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
   end if
   do iorbs=1,norbs
      call calcspinmoment(iorbs, sdot(1,iorbs), cmat1)
      !write(6,'(I3, 1x, 3(e15.6e3, 1x))') iorbs, (sdot(ni,iorbs), ni=1,3)
      write(6,228)iorbs, (sdot(ni,iorbs), ni=1,3)
   end do
   write(6,*)
   write(6,*)'outsteady: total spin-square moment       '
   write(6,*)'outsteady: <Sx^2>, <Sy^2>, <Sz^2>, <S^2>  '
   call calctotalspin2(dspin2, cmat1)
   !write(6,'(1x, 4(e15.6e3, 1x))') (dspin2(ni), ni=1,3), dspin2(1)+dspin2(2)+dspin2(3)
   write(6,229) (dspin2(ni), ni=1,3), dspin2(1)+dspin2(2)+dspin2(3)
   write(6,*)
   write(6,*)'outsteady: iorbs, jorbs, <S_i * S_j>       '
   do iorbs=1,norbs-1
      do iorbs2=iorbs+1,norbs
         call calcspinproduct(iorbs, iorbs2, dtmp1, cmat1)
         !write(6,'(I3, 1x, I3, 1x, e15.6e3, 1x)') iorbs, iorbs2, dtmp1
         write(6,230)iorbs, iorbs2, dtmp1
      end do
   end do
   if (norbs .eq. 2 .and. nspin .eq. 2) then
      call calcc12(dtmp1, cmat1)
      write(6,*)
      write(6,*)'outsteady: concurrence C12 (via spin-corr) = ', dtmp1
      call calcc12p(dtmp1, cmat1)
      write(6,*)
      write(6,*)'outsteady: concurrence C12 (via spin-flip) = ', dtmp1
   end if
   call flush(6)
end if
!
! output for thermoelectrical analysis
!
if (thermopower) then
  write(6,*)
  if (norbs .eq. 1 .and. nspin .eq. 2) then
    if (nalf .eq. 1) then
      write(6,227)iter, dinvbeta(1), eleadinfty(1,1)*hbar, epara1, epara2, occu, occd, jleft
    else if (nalf .eq. 2) then
      write(6,227)iter, dinvbeta(1), dinvbeta(2), eleadinfty(1,1)*hbar, eleadinfty(2,1)*hbar,  &
                  epara1, epara2, occu, occd, jleft, jright
    else ! multi-leads
      write(6,227)iter, (dinvbeta(ialf), ialf=1,nalf), (eleadinfty(ialf,1)*hbar, ialf=1,nalf), &
                  epara1, epara2, occu, occd, (jleads(ialf), ialf=1,nalf)
    end if
  else if (norbs .eq. 1 .and. nspin .eq. 1) then
    if (nalf .eq. 1) then
      write(6,227)iter, dinvbeta(1), eleadinfty(1,1)*hbar, epara1, occ, jleft
    else if (nalf .eq. 2) then
      write(6,227)iter, dinvbeta(1), dinvbeta(2), eleadinfty(1,1)*hbar, eleadinfty(2,1)*hbar,  &
                  epara1, occ, jleft, jright
    else ! multi-leads
      write(6,227)iter, (dinvbeta(ialf), ialf=1,nalf), (eleadinfty(ialf,1)*hbar, ialf=1,nalf), &
                  epara1, occ, (jleads(ialf), ialf=1,nalf)
    end if
  end if
  call flush(6)
end if
223 format(' THERMO ', I4, 1x, e15.6e2, 1x, e10.3e2, 1x, 2(f8.4, 1x), 2(e15.6e2, 1x))
224 format(' THERMO ', I4, 1x, 2(e15.6e2, 1x), 2(e10.3e2, 1x), 2(f8.4, 1x), 4(e15.6e2, 1x))
225 format(' THERMO ', I4, 1x, e15.6e2, 1x, e10.3e2, 1x, f8.4, 1x, 2(e15.6e2, 1x))
226 format(' THERMO ', I4, 1x, 2(e15.6e2, 1x), 2(e10.3e2, 1x), f8.4, 1x, 4(e15.6e2, 1x))
227 format(' THERMO ', I4, 1x, 30(e15.6e3, 1x))
228 format(' SM  ', I3, 1x, 3(e15.6e3, 1x))
229 format(' S2  ', 1x, 4(e15.6e3, 1x))
230 format(' S12 ', I3, 1x, I3, 1x, e15.6e3, 1x)
!
dunitt = 1.d0
dunitj = 1.d0
if (runits .eq. 1) then
  dtmp1 = 0.d0
  do ni=1,nalf
    dtmp2 = 0.d0
    do nj=1,norbs
      dtmp2 = max(dtmp2, dabs(linewidth(ni,nj)))
    end do
    dtmp1 = dtmp1 + dtmp2 / hbar
  end do
  dunitt = dtmp1                                          ! hbar / gamma      : unit of time    (inverse)
  dunitj = 1.d0 / jconst / dtmp1                          ! qe * gamma / hbar : unit of current (inverse)
!
  write(6,*)
  write(6,*)' *********** outputs in relative units **********'
  write(6,*)'         current in unit of  e*Gamma/hbar        '
  write(6,*)
  if (nalf .le. 2) then
     write(6,*)' jleft  = ', jleft * dunitj
     if (nalf .ne. 1) write(6,*)' jright = ', jright * dunitj
  else
     do ialf=1,nalf
        write(6,*)' jleads(ialf = ', ialf, ' )= ', jleads(ialf) * dunitj
     end do
  end if
  if (nspin .ne. 1) then
    write(6,*)
    write(6,*)' spin current '
    if (nalf .le. 2) then
       write(6,*)' jleftu  = ', jleftu * dunitj
       write(6,*)' jleftd  = ', jleftd * dunitj
       if (nalf .eq. 2) then
         write(6,*)' jrightu = ', jrightu * dunitj
         write(6,*)' jrightd = ', jrightd * dunitj
       end if
    else 
       do ialf=1,nalf
          write(6,*)' j_up  (ialf = ', ialf, ' )= ', jt(ialf,1) * dunitj
          write(6,*)' j_down(ialf = ', ialf, ' )= ', jt(ialf,2) * dunitj
       end do
    end if
  end if
  write(6,*)
  write(6,*)' ************************************************'
  call flush(6)
end if
!
if (norbs .eq. 2 .and. nspin .eq. 1) then
  write(6,*)
  if (nalf .eq. 1) then
    write(6,222)iter, eleadinfty(1,1)*hbar, occ, jleft*dunitj
  else if (nalf .eq. 2) then
    write(6,222)iter, eleadinfty(1,1)*hbar, eleadinfty(2,1)*hbar, occ, jleft*dunitj, jright*dunitj
  else ! multi-leads
    write(6,222)iter, (eleadinfty(ialf,1)*hbar, ialf=1,nalf), occ, (jleads(ialf)*dunitj, ialf=1,nalf)
  end if
  call flush(6)
end if
!
deallocate(sigmau, sigmad, STAT=istat)
if (lsparse) then
   deallocate(irow_tmp1, icol_tmp1, cval_tmp1, STAT=istat)
end if
deallocate(cmat1, cmat2, cmat3, cmat4, cmat5, cmat6, STAT=istat)
!
return
end subroutine outsteady
!
