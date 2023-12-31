subroutine corrfunc(iorbs0, ispin0, dt0, tmax0, ljw0, lfreq0, lread0, freq0, maxit00, crit00)
!
! Calculate 
!                 C_{AB}(t) = << A | G_s(t) | B rho >>
! if ljw0 = TRUE, calculate 
!                 C_{AB}(t) + [C_{A^+ B^+}(t)]^*
!
use matmod
use corrfuncmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
logical, intent (in) :: ljw0, lfreq0, lread0
integer, intent (in) :: iorbs0, ispin0, maxit00
real*8,  intent (in) :: dt0, tmax0, freq0, crit00
integer              :: ni, nj, nstep, ialf, isgnb, istat
integer              :: nnz
integer*8            :: lni
logical              :: fexist
real*8               :: tt
real*8               :: dt1, dt2, dt3, dt6
real*8               :: dbrhor, dbrhoi
real*8               :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6, dtmp7
real*8, allocatable  :: dvec1(:), dvec2(:), dvec3(:)
complex*16           :: ctmp1, ctmp2, ctmp3, ctmp4
complex*16, allocatable :: cmat1(:,:), cmat2(:,:), cmat3(:,:), cmat4(:,:)
complex*16, allocatable :: cmat5(:,:), cmat6(:,:), cmat7(:,:), cmat8(:,:)
complex*16, allocatable :: cvec1(:)
!
namelist / resumedos / icont_cf, lresume_cf, nresume_cf
!
iorbs_dos = iorbs0
ispin_dos = ispin0
dt_dos    = dt0
tmax_dos  = tmax0
ljw_dos   = ljw0
lfreq_dos = lfreq0
lread_dos = lread0
freq_dos  = freq0
maxit_dos = maxit00
crit_dos  = crit00
!
! check input
!
write(6,*)
write(6,*)'corrfunc: entering subroutine'
call flush(6)
if (iorbs_dos .le. 0 .or. iorbs_dos .gt. nrho) then
   write(6,*)'corrfunc: error! invalid value of iorbs_dos ', iorbs_dos
   stop
end if
!
if (ispin_dos .le. 0 .or. ispin_dos .gt. nspin) then
   write(6,*)'corrfunc: error! invalid value of ispin_dos ', ispin_dos
   stop
end if
!
write(6,*)'corrfunc: parameters used '
write(6,*)'corrfunc: iorbs_dos = ', iorbs_dos
write(6,*)'corrfunc: ispin_dos = ', ispin_dos
write(6,*)'corrfunc: ljw_dos   = ', ljw_dos
write(6,*)'corrfunc: lfreq_dos = ', lfreq_dos
!
lfreq = lfreq_dos
if (lfreq_dos) then
   freq_dos = freq_dos / hbar
   dfreq    = freq_dos
!
   write(6,*)'corrfunc: freq_dos  = ', freq_dos * hbar, ' eV'
   write(6,*)'corrfunc: maxit_dos = ', maxit_dos
   write(6,*)'corrfunc: crit_dos  = ', crit_dos
else 
   write(6,*)'corrfunc: dt_dos    = ', dt_dos
   write(6,*)'corrfunc: tmax_dos  = ', tmax_dos
end if
call flush(6)
!
if (lad) then
    if (.not. lfreq_dos) then
        write(6,*)
        write(6,*)'corrfunc: lad=T and lfreq_dos=F found '
        write(6,*)'corrfunc: code unavailable yet        '
        stop
    end if
end if
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), STAT=istat)
allocate(cmat3(nrho,nrho), cmat4(nrho,nrho), STAT=istat)
allocate(cmat5(nrho,nrho), cmat6(nrho,nrho), STAT=istat)
allocate(cmat7(nrho,nrho), cmat8(nrho,nrho), STAT=istat)
allocate(cvec1(nalf), STAT=istat)
allocate(dvec1(nalf), dvec2(nalf), dvec3(nalf), STAT=istat)
!cmat1 = czero
!cmat2 = czero
!cmat3 = czero
!cmat4 = czero
!cmat5 = czero
!cmat6 = czero
!cmat7 = czero
!cmat8 = czero
!cvec1 = czero
!dvec1 = 0.d0
!dvec2 = 0.d0
!dvec3 = 0.d0
!
if (.not. lfreq_dos) then
    lresume_cf = .true.
    icont_cf   = 0
    nresume_cf = 50
    rewind(5)
    read(5, resumedos, end=101)
    101 continue
end if
!
if (lsparse) then
    call prelude_cf_spa
end if
!
call allocate_cf
!
! initialization
!
a_ams_isgn  = 1
a_ams_iorbs = iorbs_dos
a_ams_ispin = ispin_dos
b_ams_isgn  = 2
b_ams_iorbs = iorbs_dos
b_ams_ispin = ispin_dos
!
isgnb = 2
if (ljw_dos) then
   isgnb = 1
end if
!
! steady state: carefully set bias voltage to zero for equilibrium state
!
!igroundsteady = 1
!
if (lfreq_dos) then
   write(6,*)
   write(6,*)'corrfunc: frequency-domain approach for system dynamic properties'
   call flush(6)
!
   if (ljw_dos) then
      !
      ! Re[ int_0^\infty C_{B^+ B} e^{-iwt} dt ] / pi / hbar
      ! 
      if (lsparse) then
         call solve_cf_tfqmr_spa(rho_spa, 2, b_ams_ispin, b_ams_iorbs, 2, cmat3, cmat4)
      else
         call solve_cf_tfqmr(rho, 2, b_ams_ispin, b_ams_iorbs, 2, cmat3, cmat4)
      end if
      cmat5 = dcmplx(0.5d0, 0.d0) * cmat3 + dcmplx(0.d0, 0.5d0) * cmat4
      cmat1(1:nrho,1:nrho) = cmat5(1:nrho,1:nrho)
      if (lproj) then
          call zamsmm1_proj('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat5, nrho, czero, cmat3, nrho)
      else
          call zamsmm1('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat5, nrho, czero, cmat3, nrho)
      end if
      dtmp1 = 0.d0
      do ni=1,nrho
         dtmp1 = dtmp1 + dble(cmat3(ni,ni))
      end do          
      dtmp1 = dtmp1 / pi / hbar
      !
      ! Re[ int_0^\infty C_{B B^+} e^{iwt} dt ] / pi / hbar
      ! 
      if (lsparse) then
         call solve_cf_tfqmr_spa(rho_spa, 1, b_ams_ispin, b_ams_iorbs, 1, cmat3, cmat4)
      else
         call solve_cf_tfqmr(rho, 1, b_ams_ispin, b_ams_iorbs, 1, cmat3, cmat4)
      end if
      cmat5 = dcmplx(0.5d0, 0.d0) * cmat3 + dcmplx(0.d0, 0.5d0) * cmat4
      cmat2(1:nrho,1:nrho) = cmat5(1:nrho,1:nrho)
      if (lproj) then
          call zamsmm1_proj('l', 'n', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                            cmat5, nrho, czero, cmat3, nrho)
      else
          call zamsmm1('l', 'n', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat5, nrho, czero, cmat3, nrho)
      end if
      dtmp2 = 0.d0
      do ni=1,nrho
         dtmp2 = dtmp2 + dble(cmat3(ni,ni))
      end do          
      dtmp2 = dtmp2 / pi / hbar
      write(6,*)
      write(6,1001)freq_dos * hbar, dtmp1 + dtmp2, dtmp1, dtmp2
!
! the dimension of A(w) is inverse of energy, and in unit of hbar^(-1), so need to be
! divided by hbar to output value in SI unit
!
! retarded Green's function
! G^r(\w) = -i \int_0^\infty [C_{B B^+}(t) e^{iwt} + dconjg( C_{B^+ B} e^{-iwt} )] dt
!
      if (lproj) then
          call zamsmm1_proj('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                            cmat1, nrho, czero, cmat3, nrho)
          call zamsmm1_proj('l', 'n', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                            cmat2, nrho, czero, cmat4, nrho)
      else
          call zamsmm1('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat1, nrho, czero, cmat3, nrho)
          call zamsmm1('l', 'n', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat2, nrho, czero, cmat4, nrho)
      end if
      ctmp1 = czero
      do ni=1,nrho
         ctmp1 = ctmp1 + dconjg(cmat3(ni,ni)) + cmat4(ni,ni)
      end do    
      ctmp1 = -eye * ctmp1 / hbar   ! G^r(w) in SI unit (inverse of energy) 
      ctmp2 = cunity / ctmp1        ! [G^r(w)]^{-1} in SI unit (energy)
      write(6,*)
      write(6,1003)freq_dos * hbar, dble(ctmp1), dimag(ctmp1), dble(ctmp2), dimag(ctmp2)
!
! self-energy due to electron-electron interactions
! 
      ctmp2 = ctmp2 / hbar          ! [G^r(w)]^(-1) in atomic unit 
!
      cvec1(1:nalf) = czero
      if (lband) then
         do ialf=1,nalf
            dtmp1 = 0.d0
            if (igroundsteady .eq. 1) then
               dtmp1 = eleadinfty(ialf,ispin_dos)
            end if
            do ni=1,nband
               if (lequileads) then
                  cvec1(ialf) = cvec1(ialf) + dlwidth(iorbs_dos,ispin_dos) * band_coef(ni,ispin_dos,ialf) *    &
                                band_width(ni,ispin_dos,ialf) * 5.d-1 / hbar /                                 &
                                dcmplx(freq_dos - dtmp1 - band_cent(ni,ispin_dos,ialf), band_width(ni,ispin_dos,ialf))
               else
                  cvec1(ialf) = cvec1(ialf) + linewidth(ialf,iorbs_dos) * band_coef(ni,ispin_dos,ialf) *       &
                                band_width(ni,ispin_dos,ialf) * 5.d-1 / hbar /                                 &
                                dcmplx(freq_dos - dtmp1 - band_cent(ni,ispin_dos,ialf), band_width(ni,ispin_dos,ialf))
               end if
            end do
         end do
      else 
         do ialf=1,nalf
            dtmp1 = 0.d0
            if (igroundsteady .eq. 1) then
               dtmp1 = eleadinfty(ialf,ispin_dos)
            end if
            if (lequileads) then
               cvec1(ialf) = dlwidth(iorbs_dos,ispin_dos) * bandwidth(ialf) * 5.d-1 / hbar**2 /                &
                             dcmplx(freq_dos - dtmp1 - bcenter(ialf,ispin_dos) / hbar , bandwidth(ialf) / hbar)
            else 
               cvec1(ialf) = linewidth(ialf,iorbs_dos) * bandwidth(ialf) * 5.d-1 / hbar**2 /                   &
                             dcmplx(freq_dos - dtmp1 - bcenter(ialf,ispin_dos) / hbar , bandwidth(ialf) / hbar)
            end if
         end do
      end if
      ctmp3 = czero
      do ialf=1,nalf
         ctmp3 = ctmp3 + cvec1(ialf)
      end do
!
      ctmp4 = -ctmp2 + freq_dos - ctmp3  ! Self_ee(w) - engy_level (engy_level is the 1e-orbital energy)
      write(6,*)
      write(6,1004) freq_dos*hbar, dble(ctmp4) * hbar, dimag(ctmp4) * hbar,  &
                                   dble(ctmp3) * hbar, dimag(ctmp3) * hbar
!
! lesser Green function
! G^<(\w) = 2i Re[ \int_0^\infty C_{B^+ B} e^{-iwt} dt]
!
      ctmp1 = czero
      do ni=1,nrho
         ctmp1 = ctmp1 + cmat3(ni,ni)
      end do    
      ctmp2 = eye * 2.d0 * dble(ctmp1) / hbar   ! G^<(w) in SI unit (inverse of energy) 
      write(6,*)
      write(6,1005)freq_dos*hbar, dble(ctmp2), dimag(ctmp2)
!
! electronic current 
! I_alpha = e/h \int d\ep i * 2*\Delta(\ep) * ( G^<(ep) + 2*i*f_alpha(\ep)*Im[G^r(\ep)] )
! heat current
! J_alpha = 1/h \int d\ep (\ep - \mu) * {i * 2*\Delta(\ep) * ( G^<(ep) + 2*i*f_alpha(\ep)*Im[G^r(\ep)] )}
!
      ctmp3 = czero
      do ni=1,nrho
         ctmp3 = ctmp3 + dconjg(cmat3(ni,ni)) + cmat4(ni,ni)
      end do    
      ctmp2 = -eye * ctmp3   ! G^r(w) 
!
      do ialf=1,nalf
         dtmp1 = 0.d0
         if (igroundsteady .eq. 1) then
            dtmp1 = eleadinfty(ialf,ispin_dos)
         end if
         call fermi(freq_dos, dtmp1, dinvbeta(ialf)/hbar, dtmp2)
         dtmp5 = -dimag(cvec1(ialf)) * 2.d0   ! 2*Delta(ep)
         dtmp3 = -2.d0 * dble(ctmp1)          ! i*G^<(ep)
         dtmp4 = -2.d0 * dtmp2 * dimag(ctmp2) ! i*2*i*f_alpha(ep)*Im[G^r(ep)]  
         dtmp2 = dtmp5 * (dtmp3 + dtmp4)      ! integrand for I_alpha -- dimensionless (incoming electric current)
                                              !                          (-dtmp2 is the incoming particle current)
         dtmp6 = -(freq_dos - dtmp1) * dtmp2  ! integrand for JH_alpha -- energy       (incoming heat current)
         dtmp7 = -freq_dos * dtmp2            ! integrand for JE_alpha -- energy       (incoming energy current)
                                              ! the minus sign on the rhs of above eqs for JH and JE is due to 
                                              ! the negative charge carried by an electron: e = - q_e
         dvec1(ialf) = dtmp2
         dvec2(ialf) = dtmp6 * hbar           ! in SI unit (meV)
         dvec3(ialf) = dtmp7 * hbar           ! in SI unit (meV)
      end do
      write(6,*)
      write(6,1006)freq_dos*hbar, (dvec1(ialf)*dijconst, ialf=1,nalf) ! output * incr_freq (in meV) and sum over energy --> nA
      write(6,1007)freq_dos*hbar, (dvec2(ialf)*dijconst, ialf=1,nalf) ! output * incr_freq (in meV) and sum over energy --> pW (pico-Watt)
      write(6,1008)freq_dos*hbar, (dvec3(ialf)*dijconst, ialf=1,nalf) ! output * incr_freq (in meV) and sum over energy --> pW (pico-Watt)
      call flush(6)
!
   else 
      !
      ! 2.0 * Re[ int_0^\infty C_{B^+ B} e^{iwt} dt ] / pi / hbar
      ! 
      if (lsparse) then
         call solve_cf_tfqmr_spa(rho_spa, 2, b_ams_ispin, b_ams_iorbs, 1, cmat3, cmat4)
      else
         call solve_cf_tfqmr(rho, 2, b_ams_ispin, b_ams_iorbs, 1, cmat3, cmat4)
      end if
      cmat5 = dcmplx(0.5d0, 0.d0) * cmat3 + dcmplx(0.d0, 0.5d0) * cmat4
      if (lproj) then
          call zamsmm1_proj('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                            cmat5, nrho, czero, cmat3, nrho)
      else
          call zamsmm1('l', 'c', 'n', a_ams_iorbs, a_ams_ispin, cunity,     &
                       cmat5, nrho, czero, cmat3, nrho)
      end if
      dtmp1 = 0.d0
      do ni=1,nrho
         dtmp1 = dtmp1 + dble(cmat3(ni,ni))
      end do          
      dtmp1 = dtmp1 / pi / hbar * 2.d0
      write(6,*)
      write(6,1002)freq_dos * hbar, dtmp1
   end if
   call flush(6)
!
   if ( (.not. lsparse) .and. lchkdos ) then
      write(6,*)
      write(6,*)'corrfunc: check sparsity for final brhoh '
      call checkdos(isgnb, b_ams_ispin, b_ams_iorbs, brhoh)
      write(6,*)
      write(6,*)'corrfunc: check sparsity for final brhoa '
      call checkdos(isgnb, b_ams_ispin, b_ams_iorbs, brhoa)
      call flush(6)
   end if
!
   goto 912
end if
!
1001 format('CORRFUNC_ AW   ', 2x, e15.6e3, 4(2x, e15.6e3))
1002 format('CORRFUNC_ CW   ', 2x, e15.6e3, 4(2x, e15.6e3))
1003 format('CORRFUNC_ GRW  ', 2x, e15.6e3, 4(2x, e15.6e3))
1004 format('CORRFUNC_ SELFW', 2x, e15.6e3, 4(2x, e15.6e3))
1005 format('CORRFUNC_ GLW  ', 2x, e15.6e3, 4(2x, e15.6e3))
1006 format('CORRFUNC_ IW   ', 2x, e15.6e3, 6(2x, e15.6e3))
1007 format('CORRFUNC_ JHW  ', 2x, e15.6e3, 6(2x, e15.6e3))
1008 format('CORRFUNC_ JEW  ', 2x, e15.6e3, 6(2x, e15.6e3))
!
call initialize_cf(tt)
!
! time propagation
! 
open(unit=38, file='dos.data', status='unknown')
rewind(38)
!
dt1 = dt_dos 
dt2 = dt_dos / 2.d0
dt3 = dt_dos / 3.d0
dt6 = dt_dos / 6.d0
!
100 continue
!
! Bh^- = B * rho^- + rho^- * B^+ = (Bh^+)^dag
! Ba^- = (-i) * (B * rho^- - rho^- * B^+) = (Ba^+)^dag 
! Brho = (Bh + i * Ba) / 2
!
! C_{AB}(t) = tr(A * Brho(t))
!
if (lsparse) then
   nnz = nnzcf_spa(1)
   lni = indcf_spa(1)
   call zmat_coo2dns(nnz, irowcf_spa(lni), icolcf_spa(lni), brhoh_spa(lni), cmat7, nrho, nrho, nrho)
   call zmat_coo2dns(nnz, irowcf_spa(lni), icolcf_spa(lni), brhoa_spa(lni), cmat8, nrho, nrho, nrho)
   cmat3(1:nrho,1:nrho) = 0.5d0 * cmat7(1:nrho,1:nrho) + 0.5d0 * eye * cmat8(1:nrho,1:nrho)
else
   cmat3(1:nrho,1:nrho) = 0.5d0 * brhoh(1:nrho,1:nrho,1) + 0.5d0 * eye * brhoa(1:nrho,1:nrho,1)
end if
call zgemm('n', 'n', nrho, nrho, nrho, cunity, a_ams, nrho, cmat3, nrho, czero, cmat4, nrho)
ctmp1 = czero
do ni=1,nrho
   ctmp1 = ctmp1 + cmat4(ni,ni)
end do
!
ctmp2 = czero
if (ljw_dos) then
! Bdrho = (Bdh + i * Bda) / 2
!
! [C_{A^+ B^+}(t)]^* = dconjg( tr(A^+ * Bdrho(t)) )
   if (lsparse) then
      nnz = nnzcf_spa(1)
      lni = indcf_spa(1)
      call zmat_coo2dns(nnz, irowcf_spa(lni), icolcf_spa(lni), bdrhoh_spa(lni), cmat7, nrho, nrho, nrho)
      call zmat_coo2dns(nnz, irowcf_spa(lni), icolcf_spa(lni), bdrhoa_spa(lni), cmat8, nrho, nrho, nrho)
      cmat3(1:nrho,1:nrho) = 0.5d0 * cmat7(1:nrho,1:nrho) + 0.5d0 * eye * cmat8(1:nrho,1:nrho)
   else
      cmat3(1:nrho,1:nrho) = 0.5d0 * bdrhoh(1:nrho,1:nrho,1) + 0.5d0 * eye * bdrhoa(1:nrho,1:nrho,1)
   end if
   call zgemm('c', 'n', nrho, nrho, nrho, cunity, a_ams, nrho, cmat3, nrho, czero, cmat4, nrho) 
   do ni=1,nrho
      ctmp2 = ctmp2 + dconjg(cmat4(ni,ni))
   end do
end if
!
dbrhor = dble (ctmp1 + ctmp2)
dbrhoi = dimag(ctmp1 + ctmp2)
if (ljw_dos) then
   write(6, 518)tt, dbrhor, dbrhoi, dble(ctmp1), dimag(ctmp1), dble(ctmp2), dimag(ctmp2)
   write(38,518)tt, dbrhor, dbrhoi, dble(ctmp1), dimag(ctmp1), dble(ctmp2), dimag(ctmp2)
else 
   write(6, 518)tt, dbrhor, dbrhoi
   write(38,518)tt, dbrhor, dbrhoi
end if
!
! here, use rho as brhotmp to save memory 
! important: original rho is no longer in memory!
!
! solve for Bh
if (lsparse) then
   brho0_spa(1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa)
   brho1_spa(1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa)
   rhoq_spa (1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa)
   rhop_spa (1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa)
else
   brho0(1:nrho,1:nrho,1:nunk) = brhoh(1:nrho,1:nrho,1:nunk)
   brho1(1:nrho,1:nrho,1:nunk) = brhoa(1:nrho,1:nrho,1:nunk)
   rho  (1:nrho,1:nrho,1:nunk) = brhoh(1:nrho,1:nrho,1:nunk)
   rhop (1:nrho,1:nrho,1:nunk) = brhoa(1:nrho,1:nrho,1:nunk)
end if
!
! 1st
!
if (lsparse) then
   call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
   rhoq_spa (1:lunkcf_spa) = brho0_spa(1:lunkcf_spa) + dt2 * brhorhs_spa(1:lunkcf_spa)
   rhop_spa (1:lunkcf_spa) = brho1_spa(1:lunkcf_spa) + dt2 * rhoprhs_spa(1:lunkcf_spa)
   brhoh_spa(1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa) + dt6 * brhorhs_spa(1:lunkcf_spa)
   brhoa_spa(1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa) + dt6 * rhoprhs_spa(1:lunkcf_spa)
else
   call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
   rho  (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   rhop (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
   brhoh(1:nrho,1:nrho,1:nunk) = brhoh  (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   brhoa(1:nrho,1:nrho,1:nunk) = brhoa  (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
end if
!
! 2nd
!
if (lsparse) then
   call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
   rhoq_spa (1:lunkcf_spa) = brho0_spa(1:lunkcf_spa) + dt2 * brhorhs_spa(1:lunkcf_spa)
   rhop_spa (1:lunkcf_spa) = brho1_spa(1:lunkcf_spa) + dt2 * rhoprhs_spa(1:lunkcf_spa)
   brhoh_spa(1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa) + dt3 * brhorhs_spa(1:lunkcf_spa)
   brhoa_spa(1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa) + dt3 * rhoprhs_spa(1:lunkcf_spa)
else
   call calc_ax_cf_omp(nrho, rho, rhop, brhorhs,rhoprhs)
   rho  (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   rhop (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
   brhoh(1:nrho,1:nrho,1:nunk) = brhoh  (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   brhoa(1:nrho,1:nrho,1:nunk) = brhoa  (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
end if
!
! 3rd
!
if (lsparse) then
   call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
   rhoq_spa (1:lunkcf_spa) = brho0_spa(1:lunkcf_spa) + dt1 * brhorhs_spa(1:lunkcf_spa)
   rhop_spa (1:lunkcf_spa) = brho1_spa(1:lunkcf_spa) + dt1 * rhoprhs_spa(1:lunkcf_spa)
   brhoh_spa(1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa) + dt3 * brhorhs_spa(1:lunkcf_spa)
   brhoa_spa(1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa) + dt3 * rhoprhs_spa(1:lunkcf_spa)
else 
   call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
   rho  (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt1 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   rhop (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt1 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
   brhoh(1:nrho,1:nrho,1:nunk) = brhoh  (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   brhoa(1:nrho,1:nrho,1:nunk) = brhoa  (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
end if
!
! 4th
!
if (lsparse) then
   call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
   brhoh_spa(1:lunkcf_spa) = brhoh_spa(1:lunkcf_spa) + dt6 * brhorhs_spa(1:lunkcf_spa)
   brhoa_spa(1:lunkcf_spa) = brhoa_spa(1:lunkcf_spa) + dt6 * rhoprhs_spa(1:lunkcf_spa)
else 
   call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
   brhoh(1:nrho,1:nrho,1:nunk) = brhoh  (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                 brhorhs(1:nrho,1:nrho,1:nunk)
   brhoa(1:nrho,1:nrho,1:nunk) = brhoa  (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                 rhoprhs(1:nrho,1:nrho,1:nunk)
end if
!
if (ljw_dos) then
! solve for Bdh and Bda
   if (lsparse) then
      brho0_spa(1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa)
      brho1_spa(1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa)
      rhoq_spa (1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa)
      rhop_spa (1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa)
   else
      brho0(1:nrho,1:nrho,1:nunk) = bdrhoh(1:nrho,1:nrho,1:nunk)
      brho1(1:nrho,1:nrho,1:nunk) = bdrhoa(1:nrho,1:nrho,1:nunk)
      rho  (1:nrho,1:nrho,1:nunk) = bdrhoh(1:nrho,1:nrho,1:nunk)
      rhop (1:nrho,1:nrho,1:nunk) = bdrhoa(1:nrho,1:nrho,1:nunk)
   end if
!
! 1st
!
   if (lsparse) then
      call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
      rhoq_spa  (1:lunkcf_spa) = brho0_spa (1:lunkcf_spa) + dt2 * brhorhs_spa(1:lunkcf_spa)
      rhop_spa  (1:lunkcf_spa) = brho1_spa (1:lunkcf_spa) + dt2 * rhoprhs_spa(1:lunkcf_spa)
      bdrhoh_spa(1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa) + dt6 * brhorhs_spa(1:lunkcf_spa)
      bdrhoa_spa(1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa) + dt6 * rhoprhs_spa(1:lunkcf_spa)
   else
      call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
      rho   (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      rhop  (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
      bdrhoh(1:nrho,1:nrho,1:nunk) = bdrhoh (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      bdrhoa(1:nrho,1:nrho,1:nunk) = bdrhoa (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
   end if
!
! 2nd
!
   if (lsparse) then
      call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
      rhoq_spa  (1:lunkcf_spa) = brho0_spa (1:lunkcf_spa) + dt2 * brhorhs_spa(1:lunkcf_spa)
      rhop_spa  (1:lunkcf_spa) = brho1_spa (1:lunkcf_spa) + dt2 * rhoprhs_spa(1:lunkcf_spa)
      bdrhoh_spa(1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa) + dt3 * brhorhs_spa(1:lunkcf_spa)
      bdrhoa_spa(1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa) + dt3 * rhoprhs_spa(1:lunkcf_spa)
   else
      call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
      rho   (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      rhop  (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt2 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
      bdrhoh(1:nrho,1:nrho,1:nunk) = bdrhoh (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      bdrhoa(1:nrho,1:nrho,1:nunk) = bdrhoa (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
   end if
!
! 3rd
!
   if (lsparse) then
      call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
      rhoq_spa  (1:lunkcf_spa) = brho0_spa (1:lunkcf_spa) + dt1 * brhorhs_spa(1:lunkcf_spa)
      rhop_spa  (1:lunkcf_spa) = brho1_spa (1:lunkcf_spa) + dt1 * rhoprhs_spa(1:lunkcf_spa)
      bdrhoh_spa(1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa) + dt3 * brhorhs_spa(1:lunkcf_spa)
      bdrhoa_spa(1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa) + dt3 * rhoprhs_spa(1:lunkcf_spa)
   else
      call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
      rho   (1:nrho,1:nrho,1:nunk) = brho0  (1:nrho,1:nrho,1:nunk) + dt1 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      rhop  (1:nrho,1:nrho,1:nunk) = brho1  (1:nrho,1:nrho,1:nunk) + dt1 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
      bdrhoh(1:nrho,1:nrho,1:nunk) = bdrhoh (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      bdrhoa(1:nrho,1:nrho,1:nunk) = bdrhoa (1:nrho,1:nrho,1:nunk) + dt3 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
   end if
!
! 4th
!
   if (lsparse) then
      call calc_ax_cf_spa_omp(lunkcf_spa, rhoq_spa, rhop_spa, brhorhs_spa, rhoprhs_spa)
      bdrhoh_spa(1:lunkcf_spa) = bdrhoh_spa(1:lunkcf_spa) + dt6 * brhorhs_spa(1:lunkcf_spa)
      bdrhoa_spa(1:lunkcf_spa) = bdrhoa_spa(1:lunkcf_spa) + dt6 * rhoprhs_spa(1:lunkcf_spa)
   else
      call calc_ax_cf_omp(nrho, rho, rhop, brhorhs, rhoprhs)
      bdrhoh(1:nrho,1:nrho,1:nunk) = bdrhoh (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                     brhorhs(1:nrho,1:nrho,1:nunk)
      bdrhoa(1:nrho,1:nrho,1:nunk) = bdrhoa (1:nrho,1:nrho,1:nunk) + dt6 *   &
                                     rhoprhs(1:nrho,1:nrho,1:nunk)
   end if
end if
!
! save intermediate result 
! 
tt    = tt + dt_dos
nstep = nstep + 1
!
if (lresume_cf) then
   if (nstep .gt. 0 .and. mod(nstep, nresume_cf) .eq. 0) then
      call resume_cf(tt, 2)
   end if
end if
! 
! stop if needed
!
inquire(file="stopcf", exist=fexist)
if (fexist) then
   call resume_cf(tt, 2)
   write(6,*)
   write(6,*)'corrfunc: STOPCF detected. calculation terminates '
   call flush(6) 
   goto 911
end if
!
if (tt .lt. tmax_dos) goto 100
!
911 continue
!
close(38)
518 format(f12.6, 2x, 10(e15.6e3, 2x))
!
if ( (.not. lsparse) .and. lchkdos ) then
   write(6,*)
   write(6,*)'corrfunc: check sparsity for final brhoh '
   call checkdos(2, b_ams_ispin, b_ams_iorbs, brhoh)
   write(6,*)
   write(6,*)'corrfunc: check sparsity for final brhoa '
   call checkdos(2, b_ams_ispin, b_ams_iorbs, brhoa)
   if (ljw_dos) then
      write(6,*)
      write(6,*)'corrfunc: check sparsity for final bdrhoh '
      call checkdos(1, b_ams_ispin, b_ams_iorbs, bdrhoh)
      write(6,*)
      write(6,*)'corrfunc: check sparsity for final bdrhoa '
      call checkdos(1, b_ams_ispin, b_ams_iorbs, bdrhoa)
   end if
   call flush(6)
end if
!
912 continue
!
call free_cf
deallocate(cmat1, cmat2, STAT=istat)
deallocate(cmat3, cmat4, STAT=istat)
deallocate(cmat5, cmat6, STAT=istat)
deallocate(cmat7, cmat8, STAT=istat)
deallocate(cvec1, STAT=istat)
deallocate(dvec1, dvec2, dvec3, STAT=istat)
!
write(6,*)'corrfunc: leaving subroutine'
call flush(6)
!
end subroutine corrfunc
