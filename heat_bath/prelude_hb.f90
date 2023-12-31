subroutine prelude_hb
use matmod
use tmpmatmod
use hbmod
!
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, nl, istat, imode, icor, idrude, ifreq, iorbs, ispin
integer                 :: ntmp1, ntmp2
integer                 :: nmatmp, nw0tmp, nwatmp
real*8                  :: wmin0, wmax0, yshft0
real*8                  :: dtmp1, dtmp2, dtmp3
real*8, allocatable     :: w0tmp(:), w1tmp(:), wkltmp(:), dkltmp(:)
real*8                  :: ttmin, ttmax, dtt, tt
complex*16              :: ctmp1, ctmp2, ctmp3
complex*16, allocatable :: zw0tmp(:)
!
namelist  / heat_bath / lhb, lscba, dinvbeta_hb, npade_hb, ndrude_hb, itype_hb
!
lhb         = .false.
lscba       = .false.
dinvbeta_hb = 1.d-1
npade_hb    = 8 
ndrude_hb   = 1
nmode_hb    = norbs   ! Number of coupling modes = Number of system orbitals
itype_hb    = 1       ! 1 for Drude model, 2 for super-Drude model, 3 for Einstein model
! 
! The system-bath interaction Hamiltonian assumes the form of
! Hsb = \sum_i Q_i F_i, with Q_i = a^+_i a_i (occupation number operator)                      
!
rewind(5)
read(5, heat_bath, end=101)
101 continue
!
if (.not. lhb) return
if (lscale) then
   Write(6,*)'prelude_hb: error! lscale=.true. unavailable for heat bath yet '
   stop
end if
write(6,*)
write(6,*)'prelude_hb: Consider a heat bath coupled to the impurity '
call flush(6)
!
if (ndrude_hb .lt. 1 .or. npade_hb .lt. 1) then
   write(6,*)'prelude_hb: error input! ndrude_hb, npade_hb ', ndrude_hb, npade_hb
   stop
end if
if (itype_hb .ne. 1 .and. itype_hb .ne. 2 .and. itype_hb .ne. 3) then
   write(6,*)'prelude_hb: error! unknown itype_hb ', itype_hb
   stop
end if
!
allocate(dwidth_hb(ndrude_hb), dcenter_hb(ndrude_hb), dcouple_hb(ndrude_hb,nmode_hb), STAT=istat)
allocate(zpeta_hb(npade_hb), zppole_hb(npade_hb), zpeigv_hb(npade_hb), STAT=istat)
allocate(qa_hb(nrho,nrho,nmode_hb), STAT=istat)
!
read(5,*)(dwidth_hb(ni), ni=1,ndrude_hb)
if (itype_hb .eq. 3) then
   read(5,*)(dcenter_hb(ni), ni=1,ndrude_hb)
else 
   dcenter_hb(1:ndrude_hb) = 0.d0
end if
do ni=1,nmode_hb
   read(5,*)(dcouple_hb(nj,ni), nj=1,ndrude_hb)
end do
!
dinvbeta_hb = dinvbeta_hb / hbar
dbeta_hb = 1.d0 / dinvbeta_hb
dwidth_hb(1:ndrude_hb) = dwidth_hb(1:ndrude_hb) / hbar
dcenter_hb(1:ndrude_hb) = dcenter_hb(1:ndrude_hb) / hbar
dcouple_hb(1:ndrude_hb,1:nmode_hb) = dcouple_hb(1:ndrude_hb,1:nmode_hb) / hbar
!
! output
!
write(6,*)'prelude_hb: temperature of heat bath: ', dinvbeta_hb*hbar
if (itype_hb .eq. 1) then
   write(6,*)'prelude_hb: the Drude model is adopted for heat bath              '
   write(6,*)'prelude_hb: J(w) = sum_i Gamma_i/Omega_i * w/[1+(w/Omega_i)**2]   '
else if (itype_hb .eq. 2) then
   write(6,*)'prelude_hb: the super-Drude model is adopted for heat bath        '
   write(6,*)'prelude_hb: J(w) = sum_i Gamma_i/Omega_i * w/[1+(w/Omega_i)**2]**2'
else if (itype_hb .eq. 3) then
   write(6,*)'prelude_hb: the Einstein phonon model is adopted for heat bath    '
   write(6,*)'prelude_hb: J(w) = sum_i Gamma_i/Omega_i * w *                    '
   write(6,*)'prelude_hb:        { 1 / [1 + ((w - w_i)/Omega_i)**2] +           '
   write(6,*)'prelude_hb:          1 / [1 + ((w + w_i)/Omega_i)**2] }           ' 
else
   write(6,*)'prelude_hb: error! unknown itype_hb ', itype_hb
   stop
end if
write(6,*)'prelude_hb: [N-1/N] Pade approximation is used '
write(6,*)'prelude_hb: number of  Pade poles    ', npade_hb
write(6,*)'prelude_hb: number of Drude poles    ', ndrude_hb
write(6,*)
write(6,*)'prelude_hb: system coupling Hamiltonian is n_i = sum_s a_{is}^+ a_{is}'
write(6,*)'prelude_hb: number of coupling modes ', nmode_hb
if (itype_hb .eq. 1 .or. itype_hb .eq. 2) then
   write(6,*)'prelude_hb: {Omega_i} and {Gamma_i} for bath spectral function '
   write(6,*)'prelude_hb: imode, Omega_i, Gamma_i '
   do idrude=1,ndrude_hb
      write(6,201)idrude, dwidth_hb(idrude)*hbar, (dcouple_hb(idrude,imode)*hbar, imode=1,nmode_hb)
   end do
else if (itype_hb .eq. 3) then 
   write(6,*)'prelude_hb: {Omega_i}, {w_i}, {Gamma_i} for bath spectral fnc '
   write(6,*)'prelude_hb: imode, Omega_i, w_i, Gamma_i '
   do idrude=1,ndrude_hb
      write(6,201)idrude, dwidth_hb(idrude)*hbar, dcenter_hb(idrude)*hbar, &
                  (dcouple_hb(idrude,imode)*hbar, imode=1,nmode_hb)
   end do
else 
   write(6,*)'prelude_hb: error! unknown itype_hb ', itype_hb
   stop
end if
201 format(I3, 2x, 16(1x, f14.6))
call flush(6)
!
! Pade decomposition of Bose function 1/(1-e^(-x))
!
call zpout_psd_bose(npade_hb, zpeigv_hb, zpeta_hb, zppole_hb)
!
! debugging region
!ctmp1 = dcmplx(1.d1, 1.d1)
!call bosefunc(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp2)
!call bosefunc0(ctmp1, ctmp3)
!write(6,*)
!write(6,*)'prelude_hb: debug result: '
!write(6,*)'prelude_hb: Pade approximate = ', ctmp2
!write(6,*)'prelude_hb: exact value      = ', ctmp3
!call flush(6)
!stop
! end of debugging region
!
! Only consider diagonal correlation function
!
ncor_hb = npade_hb + ndrude_hb
if (itype_hb .eq. 3) then
   ncor_hb = npade_hb + ndrude_hb * 2
end if
!
nn_hb = npade_hb
nk_hb = 0
np_hb = 0
if (itype_hb .eq. 1) then
   nn_hb = nn_hb + ndrude_hb
else if (itype_hb .eq. 2) then
   nn_hb = nn_hb + ndrude_hb
   np_hb = np_hb + ndrude_hb 
else if (itype_hb .eq. 3) then
   nk_hb = nk_hb + ndrude_hb 
else
   write(6,*)'prelude_hb: error! unknown itype_hb ', itype_hb
end if
nkp_hb = nn_hb + 2 * nk_hb + np_hb
nbath_hb = 2 * nmode_hb * nkp_hb
!
! In the present version, consider only the diagonal correlation function for the heat bath
!
allocate(cb_hb(ncor_hb,nmode_hb), cd_hb(ndrude_hb,nmode_hb), cgamma_hb(ncor_hb,nmode_hb), STAT=istat)
!
loop0: do imode=1,nmode_hb
   do icor=1,npade_hb
      ctmp1 = zppole_hb(icor) * dinvbeta_hb
      call calc_spectral_hb(itype_hb, imode, ctmp1, ctmp2)
      cb_hb(icor,imode) = -2.d0 * eye * dinvbeta_hb * zpeta_hb(icor) * ctmp2
      cgamma_hb(icor,imode) = -eye * zppole_hb(icor) * dinvbeta_hb
   end do
!
   if (itype_hb .eq. 3) then
      icor = npade_hb
      do idrude=1,ndrude_hb
         icor  = icor + 1
         ctmp1 = (dcenter_hb(idrude) - eye * dwidth_hb(idrude)) * dbeta_hb
         call bosefunc(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp2) 
         cb_hb(icor,imode) = dcouple_hb(idrude,imode) *                               &
                             dcmplx(dcenter_hb(idrude), -dwidth_hb(idrude)) * ctmp2
         cgamma_hb(icor,imode) = -eye * dcmplx(dcenter_hb(idrude), -dwidth_hb(idrude))
         icor  = icor + 1
         ctmp1 = -(dcenter_hb(idrude) + eye * dwidth_hb(idrude)) * dbeta_hb
         call bosefunc(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp2) 
         cb_hb(icor,imode) = -dcouple_hb(idrude,imode) *                              &
                              dcmplx(dcenter_hb(idrude), dwidth_hb(idrude)) * ctmp2
         cgamma_hb(icor,imode) = eye * dcmplx(dcenter_hb(idrude), dwidth_hb(idrude))
      end do
      cycle loop0
   end if
!
   do idrude=1,ndrude_hb
      icor = idrude + npade_hb
      ctmp1 = -eye * dwidth_hb(idrude) * dbeta_hb
      call bosefunc(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp2) 
      if (itype_hb .eq. 1) then
         cb_hb(icor,imode) = -eye * dcouple_hb(idrude,imode) * dwidth_hb(idrude) * ctmp2
      else if (itype_hb .eq. 2) then
         call bosefunc_der(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp3) 
         cb_hb(icor,imode) = 5.d-1 * dbeta_hb * dcouple_hb(idrude,imode) * dwidth_hb(idrude)**2 * ctmp3
      else 
         write(6,*)'prelude_hb: error! unknown itype_hb ', itype_hb
         stop
      end if
      cgamma_hb(icor,imode) = dcmplx(-dwidth_hb(idrude), 0.d0)
   end do
end do loop0
!
if (itype_hb .eq. 2) then
    do imode=1,nmode_hb
       do idrude=1,ndrude_hb
          ctmp1 = -eye * dwidth_hb(idrude) * dbeta_hb
          call bosefunc(npade_hb, zpeta_hb, zppole_hb, ctmp1, ctmp2) 
          cd_hb(idrude,imode) = -5.d-1 * eye * dcouple_hb(idrude,imode) * dwidth_hb(idrude)**2 * ctmp2
       end do
    end do
end if
!
! check C(t) (imode=1)
!
write(6,*)
write(6,*)'prelude_hb: icor, cb, gamma for imode = 1'
do icor=1,ncor_hb
   write(6,521)icor, dble(cb_hb(icor,1))*hbar**2, dimag(cb_hb(icor,1))*hbar**2, &
               dble(cgamma_hb(icor,1))*hbar, dimag(cgamma_hb(icor,1))*hbar
end do
if (itype_hb .eq. 2) then
   write(6,*)'prelude_hb: idrude, cd for imode = 1'
   do icor=1,ndrude_hb
      write(6,521)icor, dble(cd_hb(icor,1))*hbar**3, dimag(cd_hb(icor,1))*hbar**3
   end do
end if
call flush(6)
521 format(I4, 1x, 4(2x, f14.8))
!
imode = 1
!nw0tmp = 4000
nw0tmp = 16000
wmin0  = -4.d3
wmax0  =  4.d3
yshft0 = -5.d-1 * pi * dinvbeta_hb  ! in unit of hbar
!
do idrude=1,ndrude_hb
   yshft0 = max(yshft0, -dwidth_hb(idrude) * 5.d-1)
end do
do icor=1,npade_hb
   yshft0 = max(yshft0, dimag(zppole_hb(icor)) * dinvbeta_hb * 5.d-1)
end do
!
nwatmp = nw0tmp * ngl8
dtmp1  = (wmax0 - wmin0) / dble(nw0tmp)
allocate(w0tmp(nw0tmp), w1tmp(nw0tmp), STAT=istat)
allocate(wkltmp(nwatmp), dkltmp(nwatmp), zw0tmp(nwatmp), STAT=istat)
w0tmp(1) = wmin0
do ni=1,nw0tmp-1
  w1tmp(ni)   = w0tmp(ni) + dtmp1
  w0tmp(ni+1) = w1tmp(ni)
end do
w1tmp(nw0tmp) = wmax0
do nl=1,nw0tmp
  do nk=1,ngl8
    ntmp1 = (nl - 1) * ngl8 + nk
    wkltmp(ntmp1) = (w1tmp(nl) + w0tmp(nl) + (w1tmp(nl) - w0tmp(nl)) * dzgl(nk)) * 5.d-1
    dkltmp(ntmp1) = (w1tmp(nl) - w0tmp(nl)) * dwgl(nk) * 5.d-1
  end do  
end do
wkltmp(1:nwatmp) = wkltmp(1:nwatmp) / hbar  ! scale to be in unit of hbar
dkltmp(1:nwatmp) = dkltmp(1:nwatmp) / hbar
!
do ifreq=1,nwatmp
  ctmp1 = dcmplx(wkltmp(ifreq), yshft0)
  call calc_spectral_hb(itype_hb, imode, ctmp1, ctmp2)
  call bosefunc0(ctmp1 * dbeta_hb, ctmp3)
  zw0tmp(ifreq) = dkltmp(ifreq) / pi * ctmp2 * ctmp3 
end do
!
ttmin  = 0.d0
ttmax  = 2.d0
dtt    = 1.d-3
tt     = ttmin 
!
dtmp1 = 0.d0
dtmp2 = 0.d0
open(unit=64, file='corr_hb.data', status='unknown')
rewind(64)
15 continue
!
ctmp1 = czero
ctmp2 = czero
do icor=1,ncor_hb
   ctmp1 = ctmp1 + cb_hb(icor,imode) * cdexp(cgamma_hb(icor,imode) * tt)
end do
if (itype_hb .eq. 2) then
   do idrude=1,ndrude_hb
      icor = idrude + npade_hb
      ctmp1 = ctmp1 + cd_hb(idrude,imode) * tt * cdexp(cgamma_hb(icor,imode) * tt)
   end do
end if
do ifreq=1,nwatmp
   ctmp2 = ctmp2 + zw0tmp(ifreq) * cdexp(-eye * dcmplx(wkltmp(ifreq), yshft0) * tt)
end do
dtmp1 = dtmp1 + cdabs(ctmp2)
if (tt .gt. dnano) then
   dtmp2 = dtmp2 + cdabs(ctmp1 - ctmp2)
end if
write(64,520)tt, dble(ctmp1), dimag(ctmp1), dble(ctmp2), dimag(ctmp2)
call flush(64)
520 format(f12.6, 2x, 10(e13.4e3, 2x))
!
tt = tt + dtt
if (tt .le. ttmax) goto 15
!
close(64)
!
write(6,*)
write(6,"(' prelude_hb: relative discrepancy to exact corrfunc', 2x, f8.3, 2x, ' %')") dtmp2/dtmp1*1.d2
write(6,*)'prelude_hb: exact,  exact - approximate '
write(6,*)dtmp1, dtmp2
call flush(6)
!
deallocate(w0tmp, w1tmp, wkltmp, dkltmp, zw0tmp, STAT=istat)
!
!write(6,*)'prelude_hb: leaving subroutine'
!call flush(6)
!stop
!
return
!
end subroutine prelude_hb
