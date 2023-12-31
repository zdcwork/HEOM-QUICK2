subroutine evaluatepara(init)
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent(out) :: init
!
integer                 :: ni, nj, nk, nl, nm, nn, np, nq, n0, n2, nr, n2t
integer                 :: malf, nrho2
integer                 :: nhalf, linear, noprhalf, norbs2, iorbsa, iorbsb, isopr
integer                 :: istat, ialf, isgn, iorbs, iorbsp, ispin, imats, iorbs2
integer                 :: ifreq, iopr, jsgn, ilor, ifff, idraw, icor
integer                 :: ntmp1, ntmp2
logical                 :: ltmp1
real*8                  :: ttmin, ttmax, dtt, tt, dsum, cpu1, cpu2
real*8                  :: dtmp1, dtmp2, dtmp3, dtmp4
real*8                  :: wmin0, wmax0, yshft0
real*8                  :: tpsf
complex*16              :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, ctmp7, ctmp8
integer                 :: nmatmp, nw0tmp, nwatmp
real*8,     allocatable :: dmatmp(:), w0tmp(:), w1tmp(:), wkltmp(:), dkltmp(:)
complex*16, allocatable :: zbmtmp(:), zw0tmp(:)
logical                 :: lread_bcenter
!
namelist / field    / fieldtype, lreadomega, t_off, leshift_dc
namelist / ddots    / doubledot
namelist / dfield2  / dfieldtype, aD, wD, dedot
namelist / dfield1  / dfieldtype, tflip, wflip
namelist / dfield3  / dfieldtype, dedot1, dedot2, wdot1, wdot2, forcesteady, tdot12, jdot12
namelist / dfield4  / dfieldtype, dedot1, dedot2, wdot1, wdot2, forcesteady, tson, tsoff, pha
namelist / dfield5  / egate, ton_gate, toff_gate
namelist / vgate    / gatetype, gpara1, gpara2, gpara3, gpara4
namelist / dfield_hubbard / ldfield_hub, dfieldtype, dedot1, dedot2, tson, tsoff, gamadot
namelist / bfield   / lbfield, dbfield, tbfield    ! B-field only applies to transient cases (turn-on)
namelist / leadpara / lread_bcenter
namelist / coupling / readcp, readmat, lrcmplx
namelist / flux0    / megaflux, aoffL, aoffR, phioffL, phioffR, lafreq, lphifreq,        &
                      natype, nphitype,                                                  &
                      anu, aomega, atheta, phinu, phiomega
namelist / wgrids   / nwfermi, nwocc, nwvir, ifmval, irange,                             &
                      wmin, wmax, wfermid, wfermiu, dfmval, wrange, ebreak
namelist / specfunc / ispectral, egaul, egaur, esincl, esincr, chkcd, skipcorr,          &
                      chkcorr, chksgn, chkalf, chkorbs, chkorbsp, chkspin
namelist / pcontour / lpcontour, nxpmin, nxpmax, nypmin, nypmax, tour, dpgrid
namelist / band     / lband, nband, lband_same
namelist / spin3d   / lspin3d, lbfield3d, lbsite
namelist / zfs      / lzfs, d_xx, d_yy, d_zz, d_xx_1, d_yy_1, d_zz_1, d_xx_2, d_yy_2, d_zz_2  ! zero-field splitting
namelist / bathcorr / offcor, lspcor
!
! experimental feature: renormalize all equivalent leads for equilibrium properties
!              (such as occupation, reduced density matrix, and system spectral function)
! IMPORTANT:
! present limitations at lequileads=.true.
! 1. does not work for 'megaflux=.true. cases
! 2. checkspectral only works for ialf=1
!
namelist / combineleads / lequileads
!
nsgn     = 2                      ! 1 for plus, 2 for minus
nalf     = 2                      ! ialf=1 -> left lead;  ialf=2 -> right lead
nfreq    = 0
nlor     = 0
numfff   = 0
onelead  = .false.
!
!open(5, file='input.in', status='unknown')
!open(5, status='unknown')                            ! can be any file name
rewind(5)
read(5,*)init                                        ! job type    (steady|transient|more)
read(5,*)ntier                                       ! N_trun + 1  (anchor tier)
read(5,*)nmats                                       ! # of Matsubara terms (or N-pfd)
read(5,*)norbs                                       ! # of system levels
read(5,*)nspin                                       ! # of spin directions considered for system (1|2)
read(5,*)nalf                                        ! # of leads
read(5,*)(bandwidth(ni), ni=1,nalf)                  ! Band width
read(5,*)((linewidth(ni,nj), nj=1,norbs), ni=1,nalf) ! system-lead coupling
                                                     ! (data for each lead are grouped together)
read(5,*)(dinvbeta(ni), ni=1,nalf)                   ! temperature
read(5,*)((engyshift(ni,nj), nj=1,nspin), ni=1,nalf) ! variation in chemical potential for ith lead for jth spin
                                                     ! (data for each lead are grouped together)
read(5,*)tmax                                        ! length of time propagation
read(5,*)dt                                          ! time step for propagation
write(6,*)
write(6,*)'evaluatepara: input file read successfully! '
call flush(6)
!
offcor = .false.
lspcor = .false.
rewind(5)
read(5, bathcorr, end=130)
130 continue
if (offcor .and. lspcor) then
    write(6,*)
    write(6,*)'evaluatepara: error! offcor=T and lspcor=T found '
    write(6,*)'evaluatepara: code not available yet             '
    stop
end if
if (offcor) then
  write(6,*)
  write(6,*)'evaluatepara: offcor=T is invoked              '
  write(6,*)'evaluatepara: off-diagonal coupling is present '
  call flush(6)
end if
if (lspcor) then
    allocate(dlw_sp(norbs,nspin,nalf), STAT=istat)
    do ialf=1,nalf
       read(5,*)((dlw_sp(iorbs,ispin,ialf), iorbs=1,norbs), ispin=1,nspin)
    end do
    write(6,*)
    write(6,*)'evaluatepara: lspcor=T is invoked              '
    write(6,*)'evaluatepara: spin-polarized sys-bath coupling '
!
! sparsity is determined by linewidth, so linewidth cannot be zero
!
    if (lsparse) then
        do ialf=1,nalf
           do ispin=1,nspin
              do iorbs=1,norbs
                 if ( dabs(dlw_sp(iorbs,ispin,ialf)) .gt. dpico  .and.  &
                      dabs(linewidth(ialf,iorbs)) .lt. dpico ) then
                      write(6,*)
                      write(6,*)'evaluparate: wrong sparsity, check linewidth '
                      write(6,*)iorbs, ispin, ialf
                      write(6,*)dlw_sp(iorbs,ispin,ialf), linewidth(ialf,iorbs)
                      stop
                 end if
              end do
           end do
        end do
    end if
    call flush(6)
end if
!
! read specified reservoir spectral function if required
!
lband = .false.
lband_same = .true.
nband = 1
rewind(5)
read(5, band, end=122)
122 continue
if (lband) then
   if (nband .le. 0) then
      write(6,*)'evaluatepara: error! nband <= 0 found ', nband
      stop
   end if
   if (.not. psdjob) then
      write(6,*)'evaluatepara: error! presently band only works with psdjob '
      stop
   end if
!
   allocate(band_coef(nband,nspin,nalf), band_cent(nband,nspin,nalf),  &
            band_width(nband,nspin,nalf), STAT=istat)
!
! read right after the line of $band ... $end
!
   if (lband_same) then
      do nk=1,nband
         read(5,*,iostat=istat)band_coef(nk,1,1), band_width(nk,1,1), band_cent(nk,1,1)
         if (istat .ne. 0) then
            write(6,*)'evaluatapara: error reading band data, lband_same = ', lband_same
            stop
         end if
      end do
      do ni=1,nalf
         do nj=1,nspin
            if (.not. (ni .eq. 1 .and. nj .eq. 1)) then
               band_coef (1:nband,ni,nj) = band_coef (1:nband,1,1)
               band_width(1:nband,ni,nj) = band_width(1:nband,1,1)
               band_cent (1:nband,ni,nj) = band_cent (1:nband,1,1)
            end if
         end do
      end do
   else
      do ni=1,nalf
         do nj=1,nspin
            do nk=1,nband
               read(5,*,iostat=istat)band_coef(nk,nj,ni), band_width(nk,nj,ni), band_cent(nk,nj,ni)
               if (istat .ne. 0) then
                  write(6,*)'evaluatapara: error reading band data, lband_same = ', lband_same
                  stop
               end if
            end do
         end do
      end do
   end if
!
! Note: bandwidth and linewidth are NOT scaled by hbar, but band_width and band_cent are scaled by hbar
!
   band_width(1:nband,1:nspin,1:nalf) = band_width(1:nband,1:nspin,1:nalf) / hbar
   band_cent (1:nband,1:nspin,1:nalf) = band_cent (1:nband,1:nspin,1:nalf) / hbar
   write(6,*)
   write(6,*)'evaluatepara: reservoir spectral function specified by '
   write(6,*)'evaluatepara: combination of Lorentzian functions      '
   write(6,*)'evaluatepara: ialf, ispin, iband, coef, center, width  '
   do ni=1,nalf
      do nj=1,nspin
         do nk=1,nband
            write(6,519)ni, nj, nk, band_coef(nk,nj,ni), band_cent(nk,nj,ni)*hbar, band_width(nk,nj,ni)*hbar
         end do
      end do
   end do
   call flush(6)
!
end if
519 format(I2, 1x, I2, 1x, I3, 1x, 3(e14.6e3, 1x))
!
! save nalf to malf
!
malf   = nalf
nleads = malf
!
allocate(dlwidth(norbs,nspin), STAT=istat)
dlwidth(1:norbs,1:nspin) = 0.d0
if (lspcor) then
    do ialf=1,nalf
       dlwidth(1:norbs,1:nspin) = dlwidth(1:norbs,1:nspin) + dlw_sp(1:norbs,1:nspin,ialf)
    end do
else
    do ialf=1,nalf
       do iorbs=1,norbs
          dlwidth(iorbs,1:nspin) = dlwidth(iorbs,1:nspin) + linewidth(ialf,iorbs)
       end do
    end do
end if
!
! read band center if needed
!
lread_bcenter = .false.
rewind(5)
read(5, leadpara, end=106)
106 continue
bcenter(1:malf,1:nspin) = 0.d0
if (lread_bcenter) then
   read(5,*) ((bcenter(ni,nj), nj=1,nspin), ni=1,malf)
end if
write(6,*)'evaluatepara: bcenter(nalf,nspin) '
write(6,*) ((bcenter(ni,nj), nj=1,nspin), ni=1,malf)
call flush(6)
!
! bcenter must be zero for psdfff=T; see J. Chem. Phys. 152, 064107 (2020), page 3
if (psdfff) then
    do ni=1,malf
       do nj=1,nspin
          if (dabs(bcenter(ni,nj)) .gt. dpico) then
              write(6,*)'evaluatepara: error! bcenter != 0 found for psdfff '
              write(6,*)' ialf= ', ni, ' ispin= ', nj, ' bcenter= ', bcenter(ni,nj)
              stop
          end if
       end do
    end do
end if
!
! set bcenter=0 for psfjob=T (for certain itype_psf)
if (psfjob) then
    if (itype_psf .eq. 1) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=1                '
        nmats = npsf1
        write(6,*)'evaluatepara: set nmats = ', npsf1
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf1
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf1 '
        write(6,*)'evaluatepara:     rpsf1 = ', rpsf1
        tpsf = tpsf1
        write(6,*)'evaluatepara:     tpsf1 = ', tpsf1
    else if (itype_psf .eq. 2) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=2                '
        nmats = npsf2
        write(6,*)'evaluatepara: set nmats = ', npsf2
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf2
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf2 '
        write(6,*)'evaluatepara:     rpsf2 = ', rpsf2
        tpsf = tpsf2
        write(6,*)'evaluatepara:     tpsf2 = ', tpsf2
    else if (itype_psf .eq. 3) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=3                '
        nmats = npsf3
        write(6,*)'evaluatepara: set nmats = ', npsf3
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf3
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf3 '
        write(6,*)'evaluatepara:     rpsf3 = ', rpsf3
        tpsf = tpsf3
        write(6,*)'evaluatepara:     tpsf3 = ', tpsf3
    else if (itype_psf .eq. 4) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=4                '
        nmats = npsf4
        write(6,*)'evaluatepara: set nmats = ', npsf4
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf4
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf4 '
        write(6,*)'evaluatepara:     rpsf4 = ', rpsf4
        tpsf = tpsf4
        write(6,*)'evaluatepara:     tpsf4 = ', tpsf4
    else if (itype_psf .eq. 5) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=5                '
        nmats = npsf5
        write(6,*)'evaluatepara: set nmats = ', npsf5
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf5
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf5 '
        write(6,*)'evaluatepara:     rpsf5 = ', rpsf5
        tpsf = tpsf5
        write(6,*)'evaluatepara:     tpsf5 = ', tpsf5
    else if (itype_psf .eq. 6) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=6                '
        nmats = npsf6
        write(6,*)'evaluatepara: set nmats = ', npsf6
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf6
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf6 '
        write(6,*)'evaluatepara:     rpsf6 = ', rpsf6
        tpsf = tpsf6
        write(6,*)'evaluatepara:     tpsf6 = ', tpsf6
    else if (itype_psf .eq. 7) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=7                '
        nmats = npsf7
        write(6,*)'evaluatepara: set nmats = ', npsf7
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf7
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf7 '
        write(6,*)'evaluatepara:     rpsf7 = ', rpsf7
        tpsf = tpsf7
        write(6,*)'evaluatepara:     tpsf7 = ', tpsf7
    else if (itype_psf .eq. 8) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=8                '
        nmats = npsf8
        write(6,*)'evaluatepara: set nmats = ', npsf8
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf8
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf8 '
        write(6,*)'evaluatepara:     rpsf8 = ', rpsf8
        tpsf = tpsf8
        write(6,*)'evaluatepara:     tpsf8 = ', tpsf8
    else if (itype_psf .eq. 9) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=9                '
        nmats = npsf9
        write(6,*)'evaluatepara: set nmats = ', npsf9
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf9
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf9 '
        write(6,*)'evaluatepara:     rpsf9 = ', rpsf9
        tpsf = tpsf9
        write(6,*)'evaluatepara:     tpsf9 = ', tpsf9
    else if (itype_psf .eq. 10) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=10                '
        nmats = npsf10
        write(6,*)'evaluatepara: set nmats = ', npsf10
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf10
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf10 '
        write(6,*)'evaluatepara:     rpsf10 = ', rpsf10
        tpsf = tpsf10
        write(6,*)'evaluatepara:     tpsf10 = ', tpsf10
    else if (itype_psf .eq. 11) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=11                '
        nmats = npsf11
        write(6,*)'evaluatepara: set nmats = ', npsf11
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf11
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf11 '
        write(6,*)'evaluatepara:     rpsf11 = ', rpsf11
        tpsf = tpsf11
        write(6,*)'evaluatepara:     tpsf11 = ', tpsf11
    else if (itype_psf .eq. 12) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=12                '
        nmats = npsf12
        write(6,*)'evaluatepara: set nmats = ', npsf12
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf12
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf12 '
        write(6,*)'evaluatepara:     rpsf12 = ', rpsf12
        tpsf = tpsf12
        write(6,*)'evaluatepara:     tpsf12 = ', tpsf12
    else if (itype_psf .eq. 13) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=13                '
        nmats = npsf13
        write(6,*)'evaluatepara: set nmats = ', npsf13
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf13
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf13 '
        write(6,*)'evaluatepara:     rpsf13 = ', rpsf13
        tpsf = tpsf13
        write(6,*)'evaluatepara:     tpsf13 = ', tpsf13
    else if (itype_psf .eq. 14) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=14                '
        nmats = npsf14
        write(6,*)'evaluatepara: set nmats = ', npsf14
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf14
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf14 '
        write(6,*)'evaluatepara:     rpsf14 = ', rpsf14
        tpsf = tpsf14
        write(6,*)'evaluatepara:     tpsf14 = ', tpsf14
    else if (itype_psf .eq. 15) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=15                '
        nmats = npsf15
        write(6,*)'evaluatepara: set nmats = ', npsf15
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf15
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf15 '
        write(6,*)'evaluatepara:     rpsf15 = ', rpsf15
        tpsf = tpsf15
        write(6,*)'evaluatepara:     tpsf15 = ', tpsf15
    else if (itype_psf .eq. 16) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=16                '
        nmats = npsf16
        write(6,*)'evaluatepara: set nmats = ', npsf16
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf16
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf16 '
        write(6,*)'evaluatepara:     rpsf16 = ', rpsf16
        tpsf = tpsf16
        write(6,*)'evaluatepara:     tpsf16 = ', tpsf16
    else if (itype_psf .eq. 17) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=17                '
        nmats = npsf17
        write(6,*)'evaluatepara: set nmats = ', npsf17
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf17
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf17 '
        write(6,*)'evaluatepara:     rpsf17 = ', rpsf17
        tpsf = tpsf17
        write(6,*)'evaluatepara:     tpsf17 = ', tpsf17
    else if (itype_psf .eq. 18) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=18                '
        nmats = npsf18
        write(6,*)'evaluatepara: set nmats = ', npsf18
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf18
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf18 '
        write(6,*)'evaluatepara:     rpsf18 = ', rpsf18
        tpsf = tpsf18
        write(6,*)'evaluatepara:     tpsf18 = ', tpsf18
    else if (itype_psf .eq. 19) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=19                '
        nmats = npsf19
        write(6,*)'evaluatepara: set nmats = ', npsf19
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf19
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf19 '
        write(6,*)'evaluatepara:     rpsf19 = ', rpsf19
        tpsf = tpsf19
        write(6,*)'evaluatepara:     tpsf19 = ', tpsf19
   else if (itype_psf .eq. 20) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=20                '
        nmats = npsf20
        write(6,*)'evaluatepara: set nmats = ', npsf20
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf20
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf20 '
        write(6,*)'evaluatepara:     rpsf20 = ', rpsf20
        tpsf = tpsf20
        write(6,*)'evaluatepara:     tpsf20 = ', tpsf20
    else if (itype_psf .eq. 21) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=21                '
        nmats = npsf21
        write(6,*)'evaluatepara: set nmats = ', npsf21
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf21
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf21 '
        write(6,*)'evaluatepara:     rpsf21 = ', rpsf21
        tpsf = tpsf21
        write(6,*)'evaluatepara:     tpsf21 = ', tpsf21
    else if (itype_psf .eq. 22) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=22                '
        nmats = npsf22
        write(6,*)'evaluatepara: set nmats = ', npsf22
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf22
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf22 '
        write(6,*)'evaluatepara:     rpsf22 = ', rpsf22
        tpsf = tpsf22
        write(6,*)'evaluatepara:     tpsf22 = ', tpsf22
    else if (itype_psf .eq. 23) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=23                '
        nmats = npsf23
        write(6,*)'evaluatepara: set nmats = ', npsf23
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf23
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf23 '
        write(6,*)'evaluatepara:     rpsf23 = ', rpsf23
        tpsf = tpsf23
        write(6,*)'evaluatepara:     tpsf23 = ', tpsf23
    else if (itype_psf .eq. 24) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=24                '
        nmats = npsf24
        write(6,*)'evaluatepara: set nmats = ', npsf24
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf24
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf24 '
        write(6,*)'evaluatepara:     rpsf24 = ', rpsf14
        tpsf = tpsf24
        write(6,*)'evaluatepara:     tpsf24 = ', tpsf24
    else if (itype_psf .eq. 25) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=25                '
        nmats = npsf25
        write(6,*)'evaluatepara: set nmats = ', npsf25
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf25
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf25 '
        write(6,*)'evaluatepara:     rpsf25 = ', rpsf25
        tpsf = tpsf25
        write(6,*)'evaluatepara:     tpsf25 = ', tpsf25
    else if (itype_psf .eq. 26) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=26                '
        nmats = npsf26
        write(6,*)'evaluatepara: set nmats = ', npsf26
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf26
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf26 '
        write(6,*)'evaluatepara:     rpsf26 = ', rpsf26
        tpsf = tpsf26
        write(6,*)'evaluatepara:     tpsf26 = ', tpsf26
    else if (itype_psf .eq. 27) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=27                '
        nmats = npsf27
        write(6,*)'evaluatepara: set nmats = ', npsf27
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf27
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf27 '
        write(6,*)'evaluatepara:     rpsf27 = ', rpsf27
        tpsf = tpsf27
        write(6,*)'evaluatepara:     tpsf27 = ', tpsf27
    else if (itype_psf .eq. 28) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=28                '
        nmats = npsf28
        write(6,*)'evaluatepara: set nmats = ', npsf28
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf28
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf28 '
        write(6,*)'evaluatepara:     rpsf28 = ', rpsf28
        tpsf = tpsf28
        write(6,*)'evaluatepara:     tpsf28 = ', tpsf28
    else if (itype_psf .eq. 29) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=29                '
        nmats = npsf29
        write(6,*)'evaluatepara: set nmats = ', npsf29
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf29
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf29 '
        write(6,*)'evaluatepara:     rpsf29 = ', rpsf29
        tpsf = tpsf29
        write(6,*)'evaluatepara:     tpsf29 = ', tpsf29
    else if (itype_psf .eq. 30) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=30                '
        nmats = npsf30
        write(6,*)'evaluatepara: set nmats = ', npsf30
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf30
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf30 '
        write(6,*)'evaluatepara:     rpsf30 = ', rpsf30
        tpsf = tpsf30
        write(6,*)'evaluatepara:     tpsf30 = ', tpsf30
    else if (itype_psf .eq. 31) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=31                '
        nmats = npsf31
        write(6,*)'evaluatepara: set nmats = ', npsf31
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf31
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf31 '
        write(6,*)'evaluatepara:     rpsf31 = ', rpsf31
        tpsf = tpsf31
        write(6,*)'evaluatepara:     tpsf31 = ', tpsf31
    else if (itype_psf .eq. 32) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=32                '
        nmats = npsf32
        write(6,*)'evaluatepara: set nmats = ', npsf32
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf32
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf32 '
        write(6,*)'evaluatepara:     rpsf32 = ', rpsf32
        tpsf = tpsf32
        write(6,*)'evaluatepara:     tpsf32 = ', tpsf32
    else if (itype_psf .eq. 33) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=33                '
        nmats = npsf33
        write(6,*)'evaluatepara: set nmats = ', npsf33
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf33
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf33 '
        write(6,*)'evaluatepara:     rpsf33 = ', rpsf33
        tpsf = tpsf33
        write(6,*)'evaluatepara:     tpsf33 = ', tpsf33
    else if (itype_psf .eq. 34) then
        bcenter(1:malf,1:nspin) = 0.d0
        write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
        write(6,*)'evaluatepara:     itype_psf=34                '
        nmats = npsf34
        write(6,*)'evaluatepara: set nmats = ', npsf34
        bandwidth(1:malf) = dinvbeta(1:malf) * rpsf34
        write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf34 '
        write(6,*)'evaluatepara:     rpsf34 = ', rpsf34
        tpsf = tpsf34
        write(6,*)'evaluatepara:     tpsf34 = ', tpsf34
     else if (itype_psf .eq. 35) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=35                '
         nmats = npsf35
         write(6,*)'evaluatepara: set nmats = ', npsf35
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf35
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf35 '
         write(6,*)'evaluatepara:     rpsf35 = ', rpsf35
         tpsf = tpsf35
         write(6,*)'evaluatepara:     tpsf35 = ', tpsf35
     else if (itype_psf .eq. 36) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=36                '
         nmats = npsf36
         write(6,*)'evaluatepara: set nmats = ', npsf36
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf36
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf36 '
         write(6,*)'evaluatepara:     rpsf36 = ', rpsf36
         tpsf = tpsf36
         write(6,*)'evaluatepara:     tpsf36 = ', tpsf36
     else if (itype_psf .eq. 37) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=37                '
         nmats = npsf37
         write(6,*)'evaluatepara: set nmats = ', npsf37
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf37
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf37 '
         write(6,*)'evaluatepara:     rpsf37 = ', rpsf37
         tpsf = tpsf37
         write(6,*)'evaluatepara:     tpsf37 = ', tpsf37
     else if (itype_psf .eq. 38) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=38                '
         nmats = npsf38
         write(6,*)'evaluatepara: set nmats = ', npsf38
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf38
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf38 '
         write(6,*)'evaluatepara:     rpsf38 = ', rpsf38
         tpsf = tpsf38
         write(6,*)'evaluatepara:     tpsf38 = ', tpsf38
     else if (itype_psf .eq. 39) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=39                '
         nmats = npsf39
         write(6,*)'evaluatepara: set nmats = ', npsf39
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf39
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf39 '
         write(6,*)'evaluatepara:     rpsf39 = ', rpsf39
         tpsf = tpsf39
         write(6,*)'evaluatepara:     tpsf39 = ', tpsf39
     else if (itype_psf .eq. 40) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=40                '
         nmats = npsf40
         write(6,*)'evaluatepara: set nmats = ', npsf40
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf40
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf40 '
         write(6,*)'evaluatepara:     rpsf40 = ', rpsf40
         tpsf = tpsf40
         write(6,*)'evaluatepara:     tpsf40 = ', tpsf40
     else if (itype_psf .eq. 41) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=41                '
         nmats = npsf41
         write(6,*)'evaluatepara: set nmats = ', npsf41
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf41
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf41 '
         write(6,*)'evaluatepara:     rpsf41 = ', rpsf41
         tpsf = tpsf41
         write(6,*)'evaluatepara:     tpsf41 = ', tpsf41
     else if (itype_psf .eq. 42) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=42                '
         nmats = npsf42
         write(6,*)'evaluatepara: set nmats = ', npsf42
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf42
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf42 '
         write(6,*)'evaluatepara:     rpsf42 = ', rpsf42
         tpsf = tpsf42
         write(6,*)'evaluatepara:     tpsf42 = ', tpsf42
     else if (itype_psf .eq. 43) then
         bcenter(1:malf,1:nspin) = 0.d0
         write(6,*)'evaluatepara: set bcenter=0 for psfjob=T and '
         write(6,*)'evaluatepara:     itype_psf=43                '
         nmats = npsf43
         write(6,*)'evaluatepara: set nmats = ', npsf43
         bandwidth(1:malf) = dinvbeta(1:malf) * rpsf43
         write(6,*)'evaluatepara: set bandwidth = dinvbeta * rpsf43 '
         write(6,*)'evaluatepara:     rpsf43 = ', rpsf43
         tpsf = tpsf43
         write(6,*)'evaluatepara:     tpsf43 = ', tpsf43
   else
        write(6,*)'evaluatepara: error! unknown itype_psf= ', itype_psf
        stop
    end if
    call flush(6)
end if
!
! check if whether leads are equivalent
!
lequileads = .false.
rewind(5)
read(5, combineleads, end=119)
119 continue
if (lequileads) then
   !if (init .ne. 2 .and. init .ne. 3) then
   if (init .ne. 2 .and. init .ne. 3 .and. init .ne. 4) then
       !write(6,*)'evaluatepara: error! lequileads = T invoked falsely, only init=2 or 3 allowed'
       write(6,*)'evaluatepara: error! lequileads = T invoked falsely, only init=2,3,4 allowed'
       stop
   end if
end if
if (lequileads) then
   if (lband) then
      do ni=1,malf-1
         do nj=1,nspin
            do nk=1,nband
               if ( dabs(band_coef(nk,nj,ni) - band_coef(nk,nj,ni+1)) .gt. dpico .or.  &
                    dabs(band_cent(nk,nj,ni) - band_cent(nk,nj,ni+1)) .gt. dpico .or.  &
                    dabs(band_width(nk,nj,ni) - band_width(nk,nj,ni+1)) .gt. dpico ) then
                  lequileads = .false.
                  write(6,*)'evaluatepara: error! lequileads = T invoked falsely, check band input'
                  stop
               end if
            end do
         end do
         if ( dabs(dinvbeta(ni) - dinvbeta(ni+1)) .gt. dpico ) then
             lequileads = .false.
             write(6,*)'evaluatepara: error! lequileads = T invoked falsely, check dinvbeta'
             stop
         end if
      end do
   else
      do ni=1,malf-1
         if ( dabs(bandwidth(ni) - bandwidth(ni+1)) .gt. dpico .or.    &
              dabs( dinvbeta(ni) -  dinvbeta(ni+1)) .gt. dpico ) then
             lequileads = .false.
             write(6,*)'evaluatepara: error! lequileads = T invoked falsely, check bandwidth and dinvbeta'
             stop
         end if
      end do
   end if
end if
if (lequileads) then
   do ni=1,malf
      do nj=1,nspin
         if (dabs(engyshift(ni,nj)) .gt. dpico) then
             lequileads = .false.
             write(6,*)'evaluatepara: error! lequileads = T invoked falsely, engyshift=0 is required'
             stop
         end if
      end do
   end do
end if
if (lequileads) then
   do ni=1,malf-1
      do nj=1,nspin
         if (dabs(bcenter(ni,nj) - bcenter(ni+1,nj)) .gt. dpico) then
            write(6,*)'evaluatapara: error! lequileads = T invoked falsely, check bcenter '
            stop
         end if
      end do
   end do
end if
if (lequileads) then
   write(6,*)'evaluatepara: identical leads are to be combined as one lead'
   call flush(6)
end if
!
! check validity of input files
!
if (norbs .le. 0 .or. norbs .gt. maxorb) then
   write(6,*)'evaluatepara: error! wrong norbs from input ', norbs, maxorb
   stop
end if
if (nspin .ne. 1 .and. nspin .ne. 2) then
   write(6,*)'evaluatepara: error! wrong nspin from input ', nspin
   write(6,*)'evaluatepara: currently only nspin=1 or 2 is allowed '
   stop
end if
!
if (nalf .eq. 1) then
   onelead = .true.
end if
if (nalf .le. 0) then
   write(6,*)'evaluatepara: error nalf = ', nalf
   stop
end if
if (nalf .gt. maxalf) then
   write(6,*)'evaluatepara: error! nalf > maxalf '
   write(6,*)'evaluatepara: nalf, maxalf ', nalf, maxalf
   stop
end if
if (nalf .gt. 2) then
   write(6,*)'evaluatepara: multi-leads case, nalf = ', nalf
end if
!
if (nmats .lt. 0 .or. nmats .gt. maxmats) then
   write(6,*)
   write(6,*)'evaluatepara: error! nmats, maxmats ', nmats, maxmats
   stop
end if
!
if (lequileads) then
   nalf = 1
end if
!
if (psdlor) then
    call prelude_psdlor
end if
!
if (psdfff) then
    call prelude_psdfff
end if
!
if (psfjob) then
    continue
!    call prelude_psfjob
end if
!
if (lsimple) then
   if (lwalf) then
      ntmp1 = nsgn * norbs * nspin * nalf
   else
      ntmp1 = nsgn * norbs * nspin
   end if
   if (ntier-1 .gt. ntmp1) then
      ntier = ntmp1 + 1
      write(6,*)
      write(6,*)'evaluatepara: number of tiers exceeds limit '
      write(6,*)'maximal ntier allowed is ', ntier, ' with lsimple=T '
      stop
   end if
end if
!
! output job parameters
!
write(6,*)
write(6,*)'evaluatepara: system parameters used '
write(6,*)'evaluatepara: init  = ', init
write(6,*)'evaluatepara: ntier = ', ntier
write(6,*)'evaluatepara: nmats = ', nmats
if (psdlor) then
   write(6,*)'evaluatepara: nlor  = ', nlor
end if
if (psdfff) then
    write(6,*)'evaluatepara: numfff= ', numfff
end if
write(6,*)'evaluatepara: norbs = ', norbs
write(6,*)'evaluatepara: nspin = ', nspin
write(6,*)'evaluatepara: nalf  = ', nalf
if (lequileads) then
   write(6,*)'evaluatepara: malf  = ', malf
end if
!
if (.not. lband) then
   write(6,*)'evaluatepara: bandwidth(nalf) '
   write(6,1027) (bandwidth(ni), ni=1,malf)
   write(6,*)'evaluatepara: bcenter(nalf,nspin) '
   write(6,1027) ((bcenter(ni,nj), nj=1,nspin), ni=1,malf)
end if
if (lspcor) then
    write(6,*)'evaluatepara: ialf, dlw_sp(norbs,nspin,ialf) '
    do ialf=1,malf
       write(6,1028)ialf, ((dlw_sp(iorbs,ispin,ialf), iorbs=1,norbs), ispin=1,nspin)
    end do
else
    write(6,*)'evaluatepara: linewidth(nalf,norbs) '
    write(6,1027) ((linewidth(ni,nj), nj=1,norbs), ni=1,malf)
end if
write(6,*)'evaluatepara: dinvbeta(nalf) '
write(6,1027) (dinvbeta(ni), ni=1,malf)
write(6,*)'evaluatepara: engyshift(nalf,nspin) '
write(6,1027) ((engyshift(ni,nj), nj=1,nspin), ni=1,malf)
!
if (psdlor) then
    write(6,*)'evaluatepara: PSD_LOR coef, width, center '
    do ialf=1,nalf
       write(6,*)'evaluatepara: ialf = ', ialf
       do ni=1,nlor
          write(6,1027)dlor_coef(ni,ialf), dlor_width(ni,ialf), dlor_cent(ni,ialf)
       end do
    end do
end if
1027 format('evaluatepara: ', 12(e14.6e3, 2x))
1028 format('evaluatepara: ', I2, 1x, 12(e14.6e3, 2x))
call flush(6)
!
! ispectral   : type of lead spectral density function, its general form is
!               linewidth * dos(energy), with dos(energy) involving the density of
!               states at a particular energy but is dimensionless.
! dos(energy) :
!    ispectral = 0 : default, Lorentzian function : 1/((energy/bandwidth)**2+1)
!              = 1 :            Gaussian function : exp(-(energy/egau)**2)
!              = 2 :      Square of Sinc function : (sin(energy/a)/(energy/a))**2
!        otherwise : not defined
!
ispectral = 0
egaul     = 5.d0
egaur     = 5.d0
esincl    = 5.d0
esincr    = 5.d0
chkcorr   = .false.
chkcd     = .false.
skipcorr  = .false.
chksgn    = 1
chkalf    = 1
chkspin   = 1
chkorbs   = 1
chkorbsp  = 1
rewind(5)
read(5, specfunc, end=118)
118 continue
if ( .not. skipcorr .and.                                                                                         &
   ( chksgn .gt. nsgn .or. chksgn .le. 0 .or. chkalf .gt. nalf .or. chkalf .le. 0 .or. chkspin .gt. nspin .or.    &
     chkspin .le. 0 .or. chkorbs .gt. norbs .or. chkorbs .le. 0 .or. chkorbsp .gt. norbs .or. chkorbsp .le. 0 ) ) then
  write(6,*)
  write(6,*)' error! error index when checking correlation function! '
  write(6,*)' chksgn, chkalf, chkspin, chkorbs, chkorbsp, nalf, nspin, norbs '
  write(6,*)chksgn, chkalf, chkspin, chkorbs, chkorbsp, nalf, nspin, norbs
  stop
end if
if (ispectral .ne. 0) then
   if (malf .ne. 1 .and. malf .ne. 2) then
      write(6,*)'evaluatepara: multi-lead case unavailable for this case '
      write(6,*)'evaluatepara: ispectral = ', ispectral, ' nalf = ', malf
      stop
   end if
end if
if (ispectral .ne. 0 .and. (psdlor .or. psdfff)) then
   write(6,*)'evaluatepara: error! psdlor/psdfff works only with ispectral=0', ispectral
   stop
end if
if (ispectral .ne. 0 .and. lband) then
   write(6,*)'evaluatepara: error! lband works only with ispectral=0', ispectral
   stop
end if
!
yshift(1:nalf) = 1.d10
if (mfdjob) then
  nwfermi = 1           ! # of grids within fermi function range (wfermid, wfermiu)
  nwocc   = 1           ! # of grids for occupied conduction bands
  nwvir   = 1           ! # of grids for virtual  conduction bands
  ifmval  = 0           ! ifmval = 0 : wfermid and wfermiu are determined by fmval
                        !       != 0 : wfermid and wfermiu can override the limits from fmval
                        !              (in cases where wfermid and wfermiu are manually assigned)
  irange  = 1           ! irange = 0 : wmin and wmax are determined by bandwidth and wrange
                        !       != 0 : wmin and wmax can override the limits from wrange
                        !              (in cases where wmin and wmax are manually assigned)
  ebreak  = 0.d0        ! small energy shift breaking the possible symmetry of wfermid/u
  yshift(1:nalf) = 0.d0 ! shift of integration contour from real axis
  dtmp3 = 0.d0
  do ni=1,nalf
     dtmp3 = max(dtmp3, bandwidth(ni))
  end do
!
! 1/((wmin/W)^2 + 1) = wrange ==> wmin = -dsqrt(1/wrange - 1) * W
!
  wrange  = 1.d-1
  dtmp2   = dsqrt(1.d0 / wrange - 1.d0)
  wmin    = -dtmp2 * dtmp3
  wmax    = -wmin
  dfmval  = 0.999d0      ! fermi function range defined by f(wfermid) = dfmval
                         !                               & f(wfermiu) = 1 - dfmval
  dtmp1 = -1.d10
  dtmp4 = dlog(1.d0 / dfmval - 1.d0)
  do ni=1,nalf
     dtmp1 = max(dtmp1, dinvbeta(ni) * dtmp4)
  end do
!
  wfermid =  dtmp1
  wfermiu = -dtmp1
!
  rewind(5)
  read(5, wgrids, end=114)
  114 continue
  if (nwfermi .le. 0 .or. nwocc .lt. 0 .or. nwvir .lt. 0) then
    write(6,*)
    write(6,*)' error input! nwfermi, nwocc, nwvir : ', nwfermi, nwocc, nwvir
    stop
  end if
  nwl   = nwfermi + nwocc + nwvir
  nfreq = nwl * ngl8
  if (nfreq .lt. 0 .or. nfreq .gt. maxfreq) then
    write(6,*)
    write(6,*)' error for Frequency dispersion! nfreq, maxfreq ', nfreq, maxfreq
    stop
  end if
!
  if (nwocc .eq. 0) then
    wmin = -1.d10
  else if (irange .eq. 0) then
    dtmp2 = dsqrt(1.d0 / wrange - 1.d0)
    wmin  = -dtmp2 * dtmp3
  end if
  if (nwvir .eq. 0) then
    wmax = 1.d10
  else if (irange .eq. 0) then
    dtmp2 = dsqrt(1.d0 / wrange - 1.d0)
    wmax  = dtmp2 * dtmp3
  end if
  if (ifmval .eq. 0) then
     dtmp4 = dlog(1.d0 / dfmval - 1.d0)
     dtmp1 = -1.d10
     do ni=1,nalf
        dtmp1 = max(dtmp1, dinvbeta(ni) * dtmp4)
     end do
     wfermid =  dtmp1
     wfermiu = -dtmp1
  end if
  if (abs(ebreak) .ge. 1.d-2) then
    write(6,*)
    write(6,*)' warning! ebreak possibly too large '
    call flush(6)
  end if
  wfermid = wfermid + ebreak
  wfermiu = wfermiu + ebreak
  if ( .not. (wmin .lt. wfermid .and. wfermid .lt. 0.d0 .and.   &
              wfermiu .gt. 0.d0 .and. wfermiu .lt. wmax) ) then
    write(6,*)
    write(6,*)' error input! check wmin < wfermid < 0 < wfermiu < wmax '
    write(6,*)' wmin, wfermid, wfermiu, wmax : ', wmin, wfermid, wfermiu, wmax
    stop
  end if
!
  allocate(wl0(nwl), wl1(nwl), STAT=istat)
  allocate(wkln(nfreq), dkln(nfreq), STAT=istat)
  if (nwocc .gt. 0) then
    dtmp1 = wmin
    dtmp2 = wfermid
    dtmp3 = (wfermid - wmin) / dble(nwocc)
    wl0(1) = wmin
    do ni=1,nwocc-1
      wl1(ni)   = wl0(ni) + dtmp3
      wl0(ni+1) = wl1(ni)
    end do
    wl1(nwocc) = wfermid
  end if
  dtmp1 = wfermid
  dtmp2 = wfermiu
  dtmp3 = (wfermiu - wfermid) / dble(nwfermi)
  wl0(nwocc+1) = wfermid
  do ni=nwocc+1,nwocc+nwfermi-1
    wl1(ni)   = wl0(ni) + dtmp3
    wl0(ni+1) = wl1(ni)
  end do
  wl1(nwocc+nwfermi) = wfermiu
  if (nwvir .gt. 0) then
    dtmp1 = wfermiu
    dtmp2 = wmax
    dtmp3 = (wmax - wfermiu) / dble(nwvir)
    wl0(nwocc+nwfermi+1) = wfermiu
    do ni=nwocc+nwfermi+1,nwl-1
      wl1(ni)   = wl0(ni) + dtmp3
      wl0(ni+1) = wl1(ni)
    end do
    wl1(nwl) = wmax
  end if
!
! yshift can be any value between (2*nmats-1)*pi*T and (2*nmats+1)*pi*T
! choose it to be 2*nmats*pi*T
!
  if (nmats .gt. 0) then
     do ni=1,nalf
        yshift(ni) = dble(2 * nmats) * pi * dinvbeta(ni)
     end do
  else
     do ni=1,nalf
        yshift(ni) = pi * dinvbeta(ni) * 5.d-1
     end do
  end if
  call flush(6)
!
  write(6,*)
  write(6,*)' parameters determined as follows: '
  write(6,*)' ispectral = ', ispectral
  write(6,*)' nwocc     = ', nwocc
  write(6,*)' nwfermi   = ', nwfermi
  write(6,*)' nwvir     = ', nwvir
  write(6,*)' ifmval    = ', ifmval
  write(6,*)' irange    = ', irange
  write(6,*)' dfmval    = ', dfmval
  write(6,*)' wrange    = ', wrange
  write(6,*)' wmin      = ', wmin
  write(6,*)' wfermid   = ', wfermid
  write(6,*)' wfermiu   = ', wfermiu
  write(6,*)' wmax      = ', wmax
  write(6,*)' ebreak    = ', ebreak
  write(6,*)' yshift    = ', (yshift(ialf), ialf=1,nalf)
  write(6,*)
  write(6,*)' partition of frequency domain: '
  do nl=1,nwl
    write(6,618)nl, wl0(nl), wl1(nl)
  end do
  call flush(6)
!
  do nl=1,nwl
    do nk=1,ngl8
      ntmp1 = (nl - 1) * ngl8 + nk
      wkln(ntmp1) = (wl1(nl) + wl0(nl) + (wl1(nl) - wl0(nl)) * dzgl(nk)) * 5.d-1
      dkln(ntmp1) = (wl1(nl) - wl0(nl)) * dwgl(nk) * 5.d-1
    end do
  end do
  write(6,*)
  write(6,*)' grid points used in frequency domain '
  do ni=1,nfreq
    write(6,618)ni, wkln(ni), dkln(ni)
  end do
  call flush(6)
  618 format(1x, I5, 2x, f12.6, 2x, f12.6, 2x)
!
  allocate(spectral(nsgn,nalf,nspin,norbs,norbs,nfreq), STAT=istat)
  allocate(focc(nsgn,nalf,nspin,nfreq), STAT=istat)
end if
!
! # of drawers
!
ndrude = 0
if (ispectral .eq. 0) then
!
!  if the y-position of drude pole is between 0 and yshift, then consider it,
!   otherwise do not consider it
!  some possible problematic cases with two leads : T1 != T2, and the
!   drude pole is between yshift(1) and yshift(2), in this case, the
!   program terminates, and input should be adjusted to avoid this
!
  if (bandwidth(1) .le. yshift(1)) then
    ntmp1 = 1
  else
    ntmp1 = 0
  end if
  nj = 1
  do ni=2,nalf
     ntmp2 = 1
     if (bandwidth(ni) .gt. yshift(ni)) then
        ntmp2 = 0
     end if
     if (ntmp2 .ne. ntmp1) then
        write(6,*)'evaluatepara: error! number of drude poles has to be same for all leads! '
        write(6,*)'evaluatepara: ialf = ', nj, ' bandwidth, yshift ', bandwidth(nj), yshift(nj)
        write(6,*)'evaluatepara: ialf = ', ni, ' bandwidth, yshift ', bandwidth(ni), yshift(ni)
        write(6,*)'evaluatepara: change your input to avoid this pathological case '
        stop
     end if
  end do
  ndrude = ntmp1
!
  if (psfjob) then
      ndrude = 0
  end if
else if (ispectral .eq. 1) then
  ndrude = 0                      ! for Gaussian, no pole
else if (ispectral .eq. 2) then
  ndrude = 0                      ! for Square-sinc, no pole
else
  write(6,*)
  write(6,*)' error! unknown ispectral ', ispectral
  stop
end if
!
if (lband) then
   ndrude = nband
end if
!
ncor  = ndrude + nmats + nfreq + nlor + numfff
nopr  = nsgn * norbs * nspin * nalf
nsopr = nsgn * norbs * nspin
write(6,*)
write(6,*)' total basis funcs to unravel bath correlation func ', ncor
write(6,*)' # of Drude poles      ', ndrude
if (psdjob .or. psdlor .or. psdfff) then
    write(6,*)' # of Pade       terms ', nmats
else if (psfjob) then
    write(6,*)' # of Prony      terms ', nmats
else
    write(6,*)' # of Matsubara  terms ', nmats
end if
if (nfreq .gt. 0) then
    write(6,*)' # of Frequency  terms ', nfreq
end if
if (psdlor) write(6,*)' # of Lorentzian terms ', nlor
if (psdfff) write(6,*)' # of Fano       terms ', numfff
call flush(6)
if (nmats+nfreq .le. 0) then
  write(6,*)
  write(6,*)' error! no Matsubara and Frequency terms! ', nmats, nfreq
  stop
end if
!
nvar0  = ncor
nvar1  = ncor
nvar   = nopr * nvar0        ! # of drawers
nrho   = (nspin * 2)**norbs  ! nrho >= norbs
norbs2 = norbs**2
nrho2  = nrho**2
!
allocate(dmtmp1(nrho,nrho), dmtmp2(nrho,nrho), STAT=istat)
allocate(dmtmp3(nrho,nrho), dmtmp4(nrho,nrho), STAT=istat)
allocate(dmtmp5(nrho,nrho), dmtmp6(nrho,nrho), STAT=istat)
allocate(cmtmp1(nrho,nrho), cmtmp2(nrho,nrho), STAT=istat)
allocate(cmtmp3(nrho,nrho), cmtmp4(nrho,nrho), STAT=istat)
allocate(cmtmp5(nrho,nrho), cmtmp6(nrho,nrho), STAT=istat)
allocate(cmtmp7(nrho,nrho), cmtmp8(nrho,nrho), cmtmp9(nrho,nrho), STAT=istat)
allocate(cmtmpa(nrho,nrho), cmtmpb(nrho,nrho), STAT=istat)
allocate(zmtmp1(nrho,nrho), zmtmp2(nrho,nrho), STAT=istat)
allocate(cmtmp10(nrho,nrho), cmtmp11(nrho,nrho), cmtmp12(nrho,nrho), STAT=istat)
allocate(cmtmp13(nrho,nrho), cmtmp14(nrho,nrho), cmtmp15(nrho,nrho), STAT=istat)
allocate(cmtmp16(nrho,nrho), cmtmp17(nrho,nrho), cmtmp18(nrho,nrho), STAT=istat)
!
if (ltrun_der .or. lad) then
    if (nrho2 .gt. ndimbig) then
        write(6,*)
        write(6,*)'evaluatepara: error! nrho2 > ndimbig is found ', nrho2, ndimbig
        stop
    end if
    allocate(zlmtmp1(nrho2,nrho2), zlmtmp2(nrho2,nrho2), STAT=istat)
    allocate(zlmtmp3(nrho2,nrho2), zlmtmp4(nrho2,nrho2), STAT=istat)
end if
!
allocate(dmtmp1_sp(norbs,norbs), dmtmp2_sp(norbs,norbs), STAT=istat)
!
do ni=1,nalf
 tkel(ni) = dinvbeta(ni) * tconst
end do
if (funits .eq. 1) then
  tkel(1:nalf) = tkel(1:nalf) * 1.d-3
end if
write(6,*)
write(6,*)'evaluatepara: Temperature of leads (Kelvin) '
write(6,*) (tkel(ni), ni=1,nalf)
call flush(6)
!
! megaflux: inter-level correlation considered for ispectral=0, and offcor=.true.
!           consider only the case of (norbs=2, nspin=1, nalf<=2)
!           In this case,
!           assume linewidth(ialf,1) = linewidth(ialf,2),
!           and in practice, use only linewidth(ialf,1) to construct nondiagonal
!           bath spectral density matrix
!
megaflux = .false.
aoffL    = 1.d0
aoffR    = 1.d0
phioffL  = 0.d0
phioffR  = 0.d0
lafreq   = .false.
lphifreq = .false.
natype   = 1
nphitype = 1
anu      = 1.d10
aomega   = 0.d0
atheta   = 0.d0
phinu    = 1.d10
phiomega = 0.d0
rewind(5)
read(5, flux0, end=113)
113 continue
if (megaflux .and. lequileads) then
   write(6,*)'evaluatepara: error! lequileads not ready for megaflux calculations'
   stop
end if
if (megaflux) then
  write(6,*)
  write(6,*)' Inter-level correlation due to magnetic flux is considered '
!  if (.not. (ispectral .eq. 0 .and. offcor) .or. .not. (norbs .eq. 2 .and. nspin .eq. 1) .or. &
  if (.not. (ispectral .eq. 0 .and. offcor) .or. .not. (norbs .eq. 2) .or. &
      .not. (nalf .le. 2) ) then
    megaflux = .false.
    write(6,*)
    write(6,*)' error in <evaluatepara.f90>! megaflux invoked unexpectedly '
    write(6,*)' offcor, norbs, nspin, nalf ', offcor, norbs, nspin, nalf
    stop
  end if
  write(6,*)
  write(6,*)' linewidth(ialf,iorbs) set to linewidth(ialf,1) for iorbs=1,2 '
  if (megaflux .and. (lafreq .or. lphifreq)) then
    write(6,*)
    write(6,*)' cross-correlation is considered to be frequency-dependent '
    write(6,*)'   lafreq = ', lafreq,   '   natype = ', natype
    write(6,*)' lphifreq = ', lphifreq, ' nphitype = ', nphitype
  end if
  call flush(6)
end if
!
! (serial) double quantum dot case (norbs=2, nalf=2, offcor=false (diagonal linewidth matrix only))
!
doubledot = .false.
rewind(5)
read(5, ddots, end=199)
199 continue
if (doubledot) then
   if (norbs .ne. 2) then
      write(6,*)'evaluatepara: error! doubledot is invoked improperly, (norbs = 2) ', norbs
      stop
   end if
   if (malf .ne. 2) then
      write(6,*)'evaluatepara: error! doubledot is invoked improperly, (nalf = 2)  ', malf
      stop
   end if
   if (offcor) then
      write(6,*)'evaluatepara: error! doubledot is invoked improperly, (offcor=F)  ', offcor
      stop
   end if
!
   linewidth(1,2) = 0.d0  ! first  dot does not couple to right lead
   linewidth(2,1) = 0.d0  ! second dot does not couple to left  lead
end if
!
! allocate necessary memory for storage of constants
!
allocate(dpm(nsgn), STAT=istat)
allocate(dmatsu(nmats,nalf), dbeta(nalf), STAT=istat)
allocate(cgamma(nvar), dipm(nvar), ilead(nvar), iredex(nvar), STAT=istat)
allocate(eleadinfty(nalf,nspin), eleadtime(nalf,nspin), STAT=istat)
allocate(amps(malf,nspin), tchar(malf,nspin), omegas(malf,nspin), STAT=istat)
allocate(amps_dc(malf,nspin), STAT=istat)
!
!  IMPORTANT : nvar = # of drawers,
!  For Spectral  terms, cgamma(idraw) depends explicitly on system-levels
!  For Matsubara terms, cgamma(idraw) does not depend on system-levels
!  For Frequency terms, cgamma(idraw) does not depend on system-levels
!
allocate(morbs(nvar), mspin(nvar), mmats(nvar), mpm(nvar), msopr(nvar), STAT=istat)
allocate(mdrude(nvar), mfreq(nvar), mcor(nvar), mopr(nvar), STAT=istat)
allocate(mnfff(nvar), mmfff(nvar), STAT=istat)
allocate(jpm(nopr), jspin(nopr), jorbs(nopr), jalf(nopr), jredex(nopr), jomit(nopr), STAT=istat)
allocate(kmats(ncor), kpmats(ncor), dmemexp(ncor), STAT=istat)
!
if (offcor) then
  nvar2 = norbs2
else
  nvar2 = norbs
end if
!
allocate(  dwp(nvar2,nspin,malf), gammap(nvar2,nspin,malf), STAT=istat )
allocate(   cb(nvar2,nspin,ncor,nalf,nsgn), cd(nvar2,nspin,ncor,nalf,nsgn), STAT=istat )
allocate(cgama(nvar2,nspin,ncor,nalf,nsgn), STAT=istat)
if (psdfff) allocate(ctheta(nspin,ncor,nalf,nsgn), STAT=istat)  ! note: ctheta has no dependence on (iorbs,iorbsp)
if (nmats .gt. 0) then
  allocate(dimpf(nvar2,nspin,nmats,nalf,nsgn), STAT=istat)
end if
allocate(dxip(nvar2,nspin,malf), STAT=istat)
allocate(zxip(nvar2,nspin,malf), STAT=istat)
allocate(zlwmat(norbs,norbs,nspin,malf), STAT=istat)
!
! calculate Matsubara poles / PFD poles
!
dpm(1) =  1.d0
dpm(2) = -1.d0
do nj=1,nalf
  do ni=1,nmats
    dmatsu(ni,nj) = (2.d0 * dble(ni) - 1.d0) * pi * dinvbeta(nj)   ! Warning, only for cases where miu0=0.
  enddo
  dbeta(nj) = 1.d0 / dinvbeta(nj)
end do
if (pfdjob .or. psdjob .or. psdlor .or. psdfff .or. psfjob) then
   allocate(cpeigv(nmats), cpcoef(nmats), cppole(nmats, nsgn), STAT=istat)
end if
if (pfdjob) then
   call zpout_pfd(nmats, nsgn, cpeigv, cpcoef, cppole)
end if
if (psdjob .or. psdlor .or. psdfff) then
   call zpout_psd(nmats, nsgn, cpeigv, cpcoef, cppole)
end if
if (psfjob) then
    call zpout_psf(nmats, nsgn, cpcoef, cppole)
end if
!
if (ispectral .eq. 0) then  ! Detailed options for Lorentzian (Drude) type
  if (offcor) then
    if (megaflux) then
!     J_L(engy) = linewidth_L * MatrixA_L * DOS(engy) * bandwidth_L * 5.d-1
!                 MatrixA_L = (1.d0            aoffL * exp(eye*phioffL))
!                             (aoffL * exp(-eye*phioffL)           1.d0)
!                 Same for R
      zxip(1:norbs2,1:nspin,1:nalf) = czero
      zxip(1,1:nspin,1) = cunity
      zxip(4,1:nspin,1) = cunity
      zxip(3,1:nspin,1) = aoffL * cdexp(eye * phioffL)
      zxip(2,1:nspin,1) = dconjg(zxip(3,1:nspin,1))
      do iorbs=1,norbs
        do iorbs2=1,norbs
           zlwmat(iorbs,iorbs2,1:nspin,1) = zxip((iorbs2-1)*norbs+iorbs,1:nspin,1) * linewidth(1,1)
        end do
      end do
      zxip(1:norbs2,1:nspin,1) = zxip(1:norbs2,1:nspin,1) * linewidth(1,1) * bandwidth(1)**2 * 5.d-1
!
      if (nalf .ne. 1) then
         zxip(1,1:nspin,2) = cunity
         zxip(4,1:nspin,2) = cunity
         zxip(3,1:nspin,2) = aoffR * cdexp(eye * phioffR)
         zxip(2,1:nspin,2) = dconjg(zxip(3,1:nspin,2))
         do iorbs=1,norbs
           do iorbs2=1,norbs
              zlwmat(iorbs,iorbs2,1:nspin,2) = zxip((iorbs2-1)*norbs+iorbs,1:nspin,2) * linewidth(2,1)
           end do
         end do
         zxip(1:norbs2,1:nspin,2) = zxip(1:norbs2,1:nspin,2) * linewidth(2,1) * &
                                    bandwidth(2)**2 * 5.d-1
      end if
    else
      zxip(1:norbs2,1:nspin,1:malf) = czero
      zlwmat(1:norbs,1:norbs,1:nspin,1:malf) = czero
      do ialf=1,malf
        do iorbs=1,norbs
          zxip((iorbs-1)*norbs+iorbs,1:nspin,ialf) = linewidth(ialf,iorbs) *   &
                                                     bandwidth(ialf)**2 * 5.d-1
          zlwmat(iorbs,iorbs,1:nspin,ialf) = dcmplx(linewidth(ialf,iorbs), 0.d0)
        end do
      end do
      readcp  = .false.
      readmat = .false.
      lrcmplx = .false.
      rewind(5)
      read(5, coupling, end=112)
      112 continue
      if (readcp) then
        do ialf=1,malf
           dmtmp2(1:norbs,1:norbs) = 0.d0
           if (readmat) then
              do ni=1,norbs
                 if (lrcmplx) then
                    read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs), (dmtmp2(ni,nj), nj=1,norbs)
                 else
                    read(5,*,iostat=istat) (dmtmp1(ni,nj), nj=1,norbs)
                 end if
                 if (istat .ne. 0) then
                    write(6,*)'evaluatepara: error reading linewidth, ialf = ', ialf, ' iorbs = ', ni
                    stop
                 end if
              end do
           else
              if (lrcmplx) then
                 read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,norbs), nj=1,norbs), ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs)
              else
                 read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,norbs), nj=1,norbs)
              end if
              if (istat .ne. 0) then
                 write(6,*)'evaluatepara: error reading linewidth matrix for ialf = ', ialf
                 stop
              end if
           end if
           do nj=1,norbs
              do ni=1,norbs
                 zxip((nj-1)*norbs+ni,1:nspin,ialf) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj)) * bandwidth(ialf)**2 * 5.d-1
                 zlwmat(ni,nj,1:nspin,ialf) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
              end do
           end do
        end do
      end if
    end if
  else
    zlwmat(1:norbs,1:norbs,1:nspin,1:malf) = czero
    do ialf=1,malf
      do iorbs=1,norbs
        zxip(iorbs,1:nspin,ialf) = linewidth(ialf,iorbs) * bandwidth(ialf)**2 * 5.d-1
        zlwmat(iorbs,iorbs,1:nspin,ialf) = dcmplx(linewidth(ialf,iorbs), 0.d0)
      end do
    end do
  end if
!
  do nj=1,malf
    do nk=1,nvar2
      do nm=1,nspin
        dwp(nk,nm,nj)    = bcenter(nj,nm)
        gammap(nk,nm,nj) = bandwidth(nj)
      end do
    end do
  end do
!
!if (offcor) then
!    write(6,*)
!    write(6,*)'evaluatepara: check zxip '
!    write(6,*)'evaluatepara: ialf, ispin, nk, zxip(nk,ispin,ialf) '
!    do ialf=1,malf
!       do ispin=1,nspin
!          do nk=1,nvar2
!             write(6,*)ialf, ispin, nk, zxip(nk,ispin,ialf)
!          end do
!       end do
!    end do
!    call flush(6)
!end if
!
  if (.not. megaflux) then
    call checkspectral ! check spec. dens. func. (symm. and nonnega.-definite)
  else
    write(6,*)
    write(6,*)' Check spectral density of leads '
    if (dabs(aoffL) .gt. 1.d0) then
      write(6,*)' error! abs(aoffL) > 1 ', aoffL
      stop
    end if
    do ni=1,norbs
      do nj=1,norbs
        cmtmp1(ni,nj) = zxip((nj-1)*norbs+ni,1,1) / (bandwidth(1)**2 * 5.d-1)
      end do
    end do
    write(6,*)' coupling matrix for system-lead_L '
    call cmatout(norbs, nrho, cmtmp1, dmtmp1)
    if (nalf .ne. 1) then
      if (dabs(aoffR) .gt. 1.d0) then
        write(6,*)' error! abs(aoffR) > 1 ', aoffR
        stop
      end if
      do ni=1,norbs
        do nj=1,norbs
          cmtmp1(ni,nj) = zxip((nj-1)*norbs+ni,1,2) / (bandwidth(2)**2 * 5.d-1)
        end do
      end do
      write(6,*)' coupling matrix for system-lead_R '
      call cmatout(norbs, nrho, cmtmp1, dmtmp1)
    end if
    call flush(6)
  end if
!
  if (doubledot) then  ! iorbs=1 only coupled to ialf=1, norbs only to ialf=2
     zxip   = czero
     zlwmat = czero
     zxip(1,1:nspin,1)             = dcmplx(linewidth(1,1) * bandwidth(1)**2 * 5.d-1, 0.d0)
     zlwmat(1,1,1:nspin,1)         = dcmplx(linewidth(1,1), 0.d0)
     zxip(norbs,1:nspin,2)         = dcmplx(linewidth(2,norbs) * bandwidth(2)**2 * 5.d-1, 0.d0)
     zlwmat(norbs,norbs,1:nspin,2) = dcmplx(linewidth(2,norbs), 0.d0)
  end if
end if
!
if (lspcor) then   ! currently lspcor only works for offcor=F (2020-03-25)
    do ni=1,malf
       zxip(1:norbs,1:nspin,ni) = dlw_sp(1:norbs,1:nspin,ni) * bandwidth(ni)**2 * 5.d-1
    end do
end if
!
if (lequileads) then
   do ni=2,malf
      zxip(1:nvar2,1:nspin,1)           = zxip(1:nvar2,1:nspin,1) +             &
                                          zxip(1:nvar2,1:nspin,ni)
      zlwmat(1:norbs,1:norbs,1:nspin,1) = zlwmat(1:norbs,1:norbs,1:nspin,1) +   &
                                          zlwmat(1:norbs,1:norbs,1:nspin,ni)
   end do
!   if (lspcor) then
!       do ni=2,malf
!          dlw_sp(1:norbs,1:nspin,1) = dlw_sp(1:norbs,1:nspin,1) + &
!                                      dlw_sp(1:norbs,1:nspin,ni)
!       end do
!   end if
end if
!
call makespectral
!
! sort all (ncor) memory exponents in decreasing value
!
do ni=1,ncor    ! using 1st reservoir information to sort
   kmats(ni) = ni
   if (ni .le. ndrude) then
      if (lband) then
         dmemexp(ni) = band_width(ni,1,1) * hbar
      else
         if (ndrude .gt. 1) then
            write(6,*)'evaluatepara: error creating dmemexp'
            stop
         end if
         dmemexp(ni) = gammap(1,1,1)
      end if
!
   else if (ni .le. ndrude+nmats .and. nmats .gt. 0) then
      imats = ni - ndrude
      if (pfdjob .or. psdjob .or. psdlor) then
         dmemexp(ni) = dinvbeta(1) * dabs(dimag(cppole(imats,1)))
      else if (psdfff) then
         dmemexp(ni) = temfff(1) * dabs(dimag(cppole(imats,1)))
      else if (psfjob) then
         dtmp1 = dinvbeta(1) / tpsf
         dmemexp(ni) = cdabs(cppole(imats,1)) * dtmp1 * hbar
      else
         dmemexp(ni) = dmatsu(imats,1)
      end if
!
   else if (ni .gt. ndrude+nmats .and. nfreq .gt. 0) then
      dmemexp(ni) = dabs(yshift(1))
!
   else if (ni .gt. ndrude+nmats .and. nlor .gt. 0) then
      ilor = ni - ndrude - nmats
      dmemexp(ni) = dlor_width(ilor,1)
!
   else if (ni .gt. ndrude+nmats .and. numfff .gt. 0) then
       ifff = ni - ndrude - nmats
       dmemexp(ni) = gfff(ifff,1) * hbar
   end if
end do
! bubble sort
! j = kmats(i), ith longest memory corresponds to jth element in overall ncor memory components
! i --> j ==  rank --> element
do ni=1,ncor-1
   do nj=ni+1,ncor
      if (dmemexp(ni) .lt. dmemexp(nj)) then
         dtmp1 = dmemexp(ni)
         dmemexp(ni) = dmemexp(nj)
         dmemexp(nj) = dtmp1
         np = kmats(ni)
         kmats(ni) = kmats(nj)
         kmats(nj) = np
      end if
   end do
end do
!
! j = kpmats(i), ith memory component ranks the jth longest in memory
! i -- > j == element --> rank
!
do ni=1,ncor  ! ni-th element
   do nj=1,ncor  ! nj-th rank
      if (kmats(nj) .eq. ni) exit
   end do
   kpmats(ni) = nj
end do
!
write(6,*)
write(6,*)'evaluatepara: sorting memory components in decreasing decay rate '
write(6,*)'evaluatepara: rank, kmats, dmemexp                               '
do ni=1,ncor
   write(6,'(I3, 1x, I5, 1x, e15.6e3)')ni, kmats(ni), dmemexp(ni)
end do
write(6,*)'evaluatepara: imats, kpmats                                      '
do ni=1,ncor
   write(6,'(I3, 1x, I5, 1x, e15.6e3)')ni, kpmats(ni)
end do
call flush(6)
!
! Sequence of operators
!    +             +             -             -    |  nsgn  (highest priority)
!  1 -> nspin                                       |  nspin (high    priority)
!  1 -> norbs                                       |  norbs (low     priority)
!  L -> R (nalf)                                    |  norbs (lowest  priority)
!
noprhalf = nopr / 2
do nl=1,nopr
  if (nl .le. noprhalf) then
    jredex(nl) = nl + noprhalf
  else
    jredex(nl) = nl - noprhalf
  end if
end do
!
nq = 0
do isgn=1,nsgn
  do ispin=1,nspin
    do iorbs=1,norbs
      do ialf=1,nalf
        nq = nq + 1
        jpm(nq)   = isgn
        jspin(nq) = ispin
        jorbs(nq) = iorbs
        jalf(nq)  = ialf
      end do
    end do
  end do
end do
if (nq .ne. nopr) then
  write(6,*)
  write(6,*)'evaluatepara: error counting operators ', nq, nopr
  stop
end if
!
!  VERY IMPORTANT  !!!
!  Sequence of drawers :  (priority decreases from upper to lower)
!    +                 +                 -                 -          |  nsgn  (highest priority)
!  1 -> nspin                                                         |  nspin (high    priority)
!  1 -> norbs                                                         |  norbs (medium  priority)
!    L        =====>   R        =====>   L        =====>   R          |  nalf  (low     priority)
!  1 -> ndrude, ndrude+1 -> ndrude+nmats, ndrude+nmats+1 -> ncor      |  ncor  (lowest  priority)
!
nq = 0
do ni=1,nsgn
  do ispin=1,nspin
    do iorbs=1,norbs
      call look4sysopr(ni, iorbs, ispin, isopr)
      do nj=1,nalf
        call look4operator(ni, iorbs, ispin, nj, iopr)
        do nk=1,ncor
          nq = nq + 1
          mopr(nq) = iopr
          msopr(nq) = isopr
          if (ni .eq. 1) then
            dipm(nq) = 1.d0
          else
            dipm(nq) = -1.d0
          end if
          mpm(nq)    = ni
          ilead(nq)  = nj
          morbs(nq)  = iorbs
          mspin(nq)  = ispin
          mcor(nq)   = nk
          if (ndrude .gt. 0 .and. nk .le. ndrude) then
            mdrude(nq) = nk
          else
            mdrude(nq) = 0
          end if
          if (nmats .gt. 0 .and. nk .gt. ndrude .and. nk .le. ndrude+nmats) then
            mmats(nq) = nk - ndrude
          else
            mmats(nq) = 0
          end if
          if (nfreq .gt. 0 .and. nk .gt. ndrude+nmats) then
            mfreq(nq) = nk - ndrude - nmats
          else
            mfreq(nq) = 0
          end if
          if (numfff .gt. 0 .and. nk .gt. ndrude+nmats) then
              mnfff(nq) = nfff1d(nk - ndrude - nmats)
              mmfff(nq) = mfff1d(nk - ndrude - nmats)
          else
              mnfff(nq) = 0
              mmfff(nq) = 1
          end if
        end do
      end do
    end do
  end do
end do
if (nq .ne. nvar) then
  write(6,*)
  write(6,*)' error! something wrong when counting drawers ', nq, nvar
  stop
end if
!
! Find correspondence between index k and bar{k} where
! k = (alpha, diamond, n, m, s), while bar{k} = (alpha, -diamond, n, m, s)
!
nhalf = nvar / nsgn
do nl=1,nvar
 if (nl .le. nhalf) then
   iredex(nl) = nl + nhalf
 else
   iredex(nl) = nl - nhalf
 end if
end do
!
if (nfreq .gt. 0) then
  do isgn=1,nsgn
    do ialf=1,nalf
      do ispin=1,nspin
        do imats=1,nfreq
          ctmp2 = dcmplx(wkln(imats), dpm(isgn) * yshift(ialf)) / dinvbeta(ialf)
          if (dble(ctmp2) .lt. 15.d0) then
            ctmp1 = 1.d0 / (1.d0 + cdexp(ctmp2))
          else
            ctmp1 = cdexp(-ctmp2) / (1.d0 + cdexp(-ctmp2))
          end if
          if (isgn .eq. 1) then
            focc(isgn,ialf,ispin,imats) = ctmp1
          else
            focc(isgn,ialf,ispin,imats) = cunity - ctmp1
          end if
        end do
      end do
    end do
  end do
  if (pfdjob .or. psdjob .or. psdlor .or. psdfff .or. psfjob) then
     write(6,*)'evaluatepara: error! nfreq > 0 is not allowed with PFD or PSD scheme'
     write(6,*)'evaluatepara: nfreq  = ', nfreq,  ' pfdjob = ', pfdjob
     write(6,*)'evaluatepara: psdjob = ', psdjob, ' psdlor = ', psdlor
     write(6,*)'evaluatepara: psdfff = ', psdfff, ' psfjob = ', psfjob
     stop
  end if
end if
!
call cpu_time(cpu1)
!
! For cgama (decaying constant), it is essential to have its rhs expression not depend on system level,
! i.e., cgama does not involve the index "iorbs, iorbs2" explicitly on the rhs.
!
! In other words, the current code architecture cannot handle cases where cgama varies upon different
! (iorbs,iorbs2) (distinct off-diagonal dissipation rate). To cope with this, the drawers need to be
! redefined according to all distinct cgama. This requires the current (norbs,ncor) index of drawers
! to be replaced by (norbs,norbs2,ncor) with possible multiple-counting.
!                                 ----------------- obsolete comments?
!
! NOTE: special treatment for lband=.true.
!
if (lband) then
   call zxip_band(malf)
   if (lequileads) then
      do ni=2,malf
         zxip(1:nvar2,1:nspin,1) = zxip(1:nvar2,1:nspin,1) + zxip(1:nvar2,1:nspin,ni)
      end do
   end if
   do ni=1,nsgn
      do nm=1,nspin
         do nk=1,nvar2
            if (nvar2 .eq. norbs2) then
               call jmod(nk, norbs, iorbs, iorbs2)
               nq = (iorbs - 1) * norbs + iorbs2   ! transpose of nk
               if (ni .eq. 1) then
                 n2 = nq                           ! take the transpose for sigma=+
                 n2t = nk
               else
                 n2 = nk
                 n2t = nq
               end if
            else
               iorbs  = nk
               iorbs2 = nk
               nq     = nk
               n2     = nk
               n2t    = nk
            end if
            do nj=1,nalf
               do np=1,ncor
                  if (np .le. ndrude) then
                     cgama(nk,nm,np,nj,ni) = dcmplx( -band_width(np,nm,nj), dpm(ni) * band_cent(np,nm,nj) )
                     if (lfermi_exact) then
                        cb(nk,nm,np,nj,ni) = zxip(n2,nm,nj) * band_width(np,nm,nj) * band_coef(np,nm,nj) /             &
                                             ( 1.d0 + cdexp(dbeta(nj) * hbar *                                         &
                                             dcmplx(dpm(ni) * band_cent(np,nm,nj), band_width(np,nm,nj))) ) / hbar
                        cd(nk,nm,np,nj,ni) = dconjg( zxip(n2t,nm,nj) * band_width(np,nm,nj) * band_coef(np,nm,nj) ) /  &
                                             dconjg( 1.d0 + cdexp(dbeta(nj) * hbar *                                   &
                                             dcmplx(dpm(nsgn+1-ni) * band_cent(np,nm,nj), band_width(np,nm,nj))) ) / hbar
                     else
                        ctmp1 = dbeta(nj) * hbar * dcmplx(dpm(ni) * band_cent(np,nm,nj), band_width(np,nm,nj))
                        call fermifunc(ni, nj, ctmp1, ctmp2)
                        cb(nk,nm,np,nj,ni) = zxip(n2,nm,nj) * band_width(np,nm,nj) * band_coef(np,nm,nj) * ctmp2 / hbar
                        ctmp1 = dbeta(nj) * hbar * dcmplx(dpm(nsgn+1-ni) * band_cent(np,nm,nj), band_width(np,nm,nj))
                        call fermifunc(nsgn+1-ni, nj, ctmp1, ctmp2)
                        cd(nk,nm,np,nj,ni) = dconjg(zxip(n2t,nm,nj) * band_width(np,nm,nj) * band_coef(np,nm,nj) * ctmp2) / hbar
                     end if
!
                  else if (np .le. ndrude+nmats .and. nmats .gt. 0) then
! PSD scheme
                     imats = np - ndrude
                     ctmp1 = eye * 2.d0 * cpcoef(imats) * zxip(n2,nm,nj) * dinvbeta(nj) / hbar**2
                     ctmp4 = eye * 2.d0 * cpcoef(imats) * zxip(n2t,nm,nj) * dinvbeta(nj) / hbar**2
                     cgama(nk,nm,np,nj,ni) = dpm(ni) * eye * dinvbeta(nj) * cppole(imats,ni) / hbar
                     cb(nk,nm,np,nj,ni)    = czero
                     cd(nk,nm,np,nj,ni)    = czero
                     do nr=1,nband
                        ctmp2 = (cppole(imats,ni) * dinvbeta(nj) / hbar - band_cent(nr,nm,nj))**2 +           &
                                band_width(nr,nm,nj)**2
                        ctmp3 = (cppole(imats,nsgn+1-ni) * dinvbeta(nj) / hbar - band_cent(nr,nm,nj))**2 +    &
                                band_width(nr,nm,nj)**2
                        cb(nk,nm,np,nj,ni) = cb(nk,nm,np,nj,ni) +                                             &
                                             ctmp1 * band_width(nr,nm,nj)**2 * band_coef(nr,nm,nj) / ctmp2
                        cd(nk,nm,np,nj,ni) = cd(nk,nm,np,nj,ni) + dconjg(                                     &
                                             ctmp4 * band_width(nr,nm,nj)**2 * band_coef(nr,nm,nj) / ctmp3 )
                     end do
                  else
! no other possibilities
                     write(6,*)'evaluatepara: error! ncor > ndrude + nmats found for lband=.true.'
                     write(6,*)'evaluatepara: ncor, ndrude, nmats', ncor, ndrude, nmats
                     stop
                  end if
               end do  ! ncor
            end do     ! nalf
         end do        ! nvar2
      end do           ! nspin
   end do              ! nsgn
!
   goto 999
end if
!
if (ispectral .eq. 0) then
  do ni=1,nsgn
    do nm=1,nspin
      do nk=1,nvar2
        if (nvar2 .eq. norbs2) then
          call jmod(nk, norbs, iorbs, iorbs2)
          nq = (iorbs - 1) * norbs + iorbs2   ! transpose of nk
          if (ni .eq. 1) then
            n2 = nq                           ! take the transpose for sigma=+
            n2t = nk                          ! n2t is transpose index of n2 (for opposite sign of sigma)
          else
            n2 = nk
            n2t = nq
          end if
        else
          iorbs  = nk
          iorbs2 = nk
          nq = nk
          n2 = nk
          n2t = nk
        end if
        do nj=1,nalf
          do np=1,ncor
            if (np .le. ndrude) then
              cgama(nk,nm,np,nj,ni) = dcmplx( -gammap(nk,nm,nj), dpm(ni) * dwp(nk,nm,nj) ) / hbar
              if (megaflux .and. (lafreq .or. lphifreq)) then   ! to modify, note taken on 2013-10-22
                 ctmp1 = dcmplx(dwp(nk,nm,nj), dpm(ni) * gammap(nk,nm,nj))
                 call crosscorr(ctmp2, nj, nk, ctmp1)
                 zxip(nk,nm,nj) = ctmp2 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                 call crosscorr(ctmp3, nj, nq, ctmp1)
                 zxip(nq,nm,nj) = ctmp3 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
              end if
              if (lfermi_exact) then
                 cb(nk,nm,np,nj,ni) = zxip(n2,nm,nj) / gammap(n2,nm,nj) / ( 1.d0 + cdexp(dbeta(nj) *         &
                                      dcmplx(dpm(ni) * dwp(n2,nm,nj), gammap(n2,nm,nj))) ) / hbar**2
! modified on 09 Oct, 2012 (ver=-0.882), modified again on 2013-10-21
! see JCP 128, 234703 (2008) Eq. (3.6)
! note: Here nk should better be replaced by n2 or n2t
                 cd(nk,nm,np,nj,ni) = dconjg( zxip(n2t,nm,nj) / gammap(n2t,nm,nj) ) /                        &
                                      dconjg( 1.d0 + cdexp(dbeta(nj) *                                       &
                                      dcmplx(dpm(nsgn+1-ni) * dwp(n2t,nm,nj), gammap(n2t,nm,nj))) ) /        &
                                      hbar**2
              else
                 ctmp1 = dbeta(nj) * dcmplx(dpm(ni) * dwp(n2,nm,nj), gammap(n2,nm,nj))
                 if (psdfff) then
                     ctmp1 = 1.d0 / temfff(nj) * dcmplx(dpm(ni) * dwp(n2,nm,nj), gammap(n2,nm,nj))
                 end if
                 call fermifunc(ni, nj, ctmp1, ctmp2)
                 cb(nk,nm,np,nj,ni) = zxip(n2,nm,nj) / gammap(n2,nm,nj) * ctmp2 / hbar**2
!
                 ctmp1 = dbeta(nj) * dcmplx(dpm(nsgn+1-ni) * dwp(n2t,nm,nj), gammap(n2t,nm,nj))
                 if (psdfff) then
                     ctmp1 = 1.d0 / temfff(nj) * dcmplx(dpm(nsgn+1-ni) * dwp(n2t,nm,nj), gammap(n2t,nm,nj))
                 end if
                 call fermifunc(nsgn+1-ni, nj, ctmp1, ctmp2)
                 cd(nk,nm,np,nj,ni) = dconjg(zxip(n2t,nm,nj) / gammap(n2t,nm,nj) * ctmp2) / hbar**2
              end if
              if (psdfff) then
                  do ifff=1,nfff
                    cb(nk,nm,np,nj,ni) = cb(nk,nm,np,nj,ni) +                          &
                                         eye * bfff(ifff) * zxip(n2,nm,nj) *           &
                                         dfff(ifff,nm,nj,ni) / hbar**3
                    cd(nk,nm,np,nj,ni) = cd(nk,nm,np,nj,ni) + dconjg(                  &
                                         eye * bfff(ifff) * zxip(n2t,nm,nj) *          &
                                         dfff(ifff,nm,nj,nsgn+1-ni) ) / hbar**3
                  end do
              end if
            else if (np .le. ndrude+nmats .and. nmats .gt. 0) then
              imats = np - ndrude
              if (pfdjob .or. psdjob .or. psdlor) then
! PFD/PSD scheme
                  cgama(nk,nm,np,nj,ni) = dpm(ni) * eye * dinvbeta(nj) * cppole(imats,ni) / hbar
                  if (megaflux .and. (lafreq .or. lphifreq)) then
                     ! this part needs to be checked, noted on 2013-10-21
                     ctmp1 = cppole(imats,ni) * dinvbeta(nj)
                     call crosscorr(ctmp2, nj, nk, ctmp1)
                     zxip(nk,nm,nj) = ctmp2 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                     call crosscorr(ctmp3, nj, nq, ctmp1)
                     zxip(nq,nm,nj) = ctmp3 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                  end if
                  cb(nk,nm,np,nj,ni) = eye * 2.d0 * cpcoef(imats) * zxip(n2,nm,nj) * dinvbeta(nj) /           &
                                       ( (cppole(imats,ni) * dinvbeta(nj) - dwp(n2,nm,nj))**2 +               &
                                         gammap(n2,nm,nj)**2 ) / hbar**2
                  cd(nk,nm,np,nj,ni) = dconjg(eye * 2.d0 * cpcoef(imats) * zxip(n2t,nm,nj) * dinvbeta(nj)) /  &
                                       dconjg( ((cppole(imats,nsgn+1-ni) * dinvbeta(nj) -                     &
                                               dwp(n2t,nm,nj))**2 + gammap(n2t,nm,nj)**2) ) / hbar**2
              else if (psdfff) then
! PSD-FFF scheme (PSD done at the reference temperature)
                  cgama(nk,nm,np,nj,ni) = dpm(ni) * eye * temfff(nj) * cppole(imats,ni) / hbar
                  cb(nk,nm,np,nj,ni) = eye * 2.d0 * cpcoef(imats) * zxip(n2,nm,nj) * temfff(nj) /             &
                                       ( (cppole(imats,ni) * temfff(nj) - dwp(n2,nm,nj))**2 +                 &
                                         gammap(n2,nm,nj)**2 ) / hbar**2
                  cd(nk,nm,np,nj,ni) = dconjg(eye * 2.d0 * cpcoef(imats) * zxip(n2t,nm,nj) * temfff(nj)) /    &
                                       dconjg( ((cppole(imats,nsgn+1-ni) * temfff(nj) -                       &
                                               dwp(n2t,nm,nj))**2 + gammap(n2t,nm,nj)**2) ) / hbar**2
!
              else if (psfjob) then
! PSF scheme
                  dtmp1 = dinvbeta(nj) / tpsf
                  cgama(nk,nm,np,nj,ni) = dtmp1 * cppole(imats,ni)
                  ctmp1 = zxip(n2,nm,nj) / bandwidth(nj)**2 * 2.d0
                  cb(nk,nm,np,nj,ni) = dtmp1 * cpcoef(imats) * ctmp1
                  ctmp1 = zxip(n2t,nm,nj) / bandwidth(nj)**2 * 2.d0
                  cd(nk,nm,np,nj,ni) = dconjg(dtmp1 * cpcoef(imats) * ctmp1)
!
              else
                cgama(nk,nm,np,nj,ni) = dcmplx(-dmatsu(imats,nj), 0.d0) / hbar
                if (megaflux .and. (lafreq .or. lphifreq)) then
                   ctmp1 = dcmplx(0.d0, dpm(ni) * dmatsu(imats,nj))
                   call crosscorr(ctmp2, nj, nk, ctmp1)
                   zxip(nk,nm,nj) = ctmp2 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                   call crosscorr(ctmp3, nj, nq, ctmp1)
                   zxip(nq,nm,nj) = ctmp3 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                end if
                cb(nk,nm,np,nj,ni)    = -2.d0 * eye * zxip(n2,nm,nj) * dinvbeta(nj) /                          &
                                         ( dcmplx(dwp(n2,nm,nj), -dpm(ni) * dmatsu(imats,nj))**2 +             &
                                           gammap(n2,nm,nj)**2 ) / hbar**2
! modified on 09 Oct, 2012 (ver=-0.882), modified again on 2013-10-21
                cd(nk,nm,np,nj,ni)    =  dconjg( -2.d0 * eye * zxip(n2t,nm,nj) * dinvbeta(nj) ) /              &
                                         dconjg(                                                               &
                                         ( dcmplx(dwp(n2t,nm,nj), -dpm(nsgn+1-ni) * dmatsu(imats,nj))**2 +     &
                                           gammap(n2t,nm,nj)**2 ) ) / hbar**2
              end if
              dimpf(nk,nm,imats,nj,ni) = dsqrt(cdabs(cb(nk,nm,np,nj,ni))) / cdabs(cgama(nk,nm,np,nj,ni))
            else if (np .gt. ndrude+nmats .and. nfreq .gt. 0) then
! MFD scheme
! offcor case need to be checked
              ifreq = np - ndrude - nmats
              cgama(nk,nm,np,nj,ni) = dcmplx( -yshift(nj), dpm(ni) * wkln(ifreq) ) / hbar
              cb(nk,nm,np,nj,ni)    = dkln(ifreq) * spectral(ni,nj,nm,iorbs,iorbs2,ifreq) *                  &
                                      focc(ni,nj,nm,ifreq) / (2.d0 * pi) / hbar**2
              cd(nk,nm,np,nj,ni)    = dkln(ifreq) * spectral(ni,nj,nm,iorbs,iorbs2,ifreq) *                  &
                                      (1.d0 - focc(ni,nj,nm,ifreq)) / (2.d0 * pi) / hbar**2
!
            else if (np .gt. ndrude+nmats .and. nlor .gt. 0) then
! PSD_LOR scheme
              ilor = np - ndrude - nmats
              cgama(nk,nm,np,nj,ni) = dcmplx( -dlor_width(ilor,nj), dpm(ni) * dlor_cent(ilor,nj) ) / hbar
              if (megaflux .and. (lafreq .or. lphifreq)) then
                 ctmp1 = dcmplx(dlor_cent(ilor,nj), dpm(ni) * dlor_width(ilor,nj))
                 call crosscorr(ctmp2, nj, nk, ctmp1)
                 zxip(nk,nm,nj) = ctmp2 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
                 call crosscorr(ctmp3, nj, nq, ctmp1)
                 zxip(nq,nm,nj) = ctmp3 * linewidth(nj,1) * bandwidth(nj)**2 * 5.d-1
              end if
              ctmp1 = zxip(n2,nm,nj) / dlor_width(ilor,nj) * dlor_coef(ilor,nj)
! modified on 2013-10-22
! offcor case needs to be checked
! mind the minus sign between ni = 1 and ni = 2 -- definition of dlor_coef
              ctmp2 = zxip(n2t,nm,nj) / dlor_width(ilor,nj) * dlor_coef(ilor,nj)
              if (ni .eq. 2) then
                 cb(nk,nm,np,nj,ni) = -ctmp1 / hbar**2
                 cd(nk,nm,np,nj,ni) =  dcmplx(ctmp2) / hbar**2  ! dconjg of cb for ni==1, also transpose for offcor case
              else
                 cb(nk,nm,np,nj,ni) =  ctmp1 / hbar**2
                 cd(nk,nm,np,nj,ni) = -dcmplx(ctmp2) / hbar**2
              end if
!
            else if (np .gt. ndrude+nmats .and. numfff .gt. 0) then
! PSD-FFF scheme
              ifff = np - ndrude - nmats
              cgama(nk,nm,np,nj,ni) = dcmplx( -gfff(ifff,nj), 0.d0 )
              cb(nk,nm,np,nj,ni) = eye * bfff(nfff1d(ifff)) * zxip(n2,nm,nj) *                          &
                                   cfff(mfff1d(ifff),nfff1d(ifff),nm,nj,ni) / hbar**3
              cd(nk,nm,np,nj,ni) = dconjg( eye * bfff(nfff1d(ifff)) * zxip(n2t,nm,nj) *                 &
                                           cfff(mfff1d(ifff),nfff1d(ifff),nm,nj,nsgn+1-ni) ) / hbar**3

            end if
          end do
        end do
      end do
    end do
  end do
!
!  if (lspcor) then
!      do ni=1,nsgn
!         do nm=1,nspin
!            do nk=1,nvar2
!               do nj=1,nalf
!                  do np=1,ncor
!                     cb(nk,nm,np,nj,ni) = cb(nk,nm,np,nj,ni) / linewidth(nj,nk) * dlw_sp(nk,nm,nj)
!                     cd(nk,nm,np,nj,ni) = cd(nk,nm,np,nj,ni) / linewidth(nj,nk) * dlw_sp(nk,nm,nj)
!                  end do
!               end do
!            end do
!         end do
!      end do
!  end if
!
  if (psdfff) then
      ctheta = czero
      dtmp1 = 0.d0
      dtmp2 = 1.d20
      write(6,*)
      write(6,*)'evaluatepara: print cb and gama '
      nq = 0
      do ni=1,nsgn
         do nm=1,nspin
            do nj=1,nalf
               do np=1,ncor
                  ctmp1 = czero
                  do nk=1,nvar2
                     nq = nq + 1
                     write(6,1101)nq, dble(cb(nk,nm,np,nj,ni)), dimag(cb(nk,nm,np,nj,ni)),  &
                                  dble(cgama(nk,nm,np,nj,ni)), dimag(cgama(nk,nm,np,nj,ni))
                     ctmp2 = cb(nk,nm,np  ,nj,ni)
                     ctmp3 = cb(nk,nm,np-1,nj,ni)
                     if ( cdabs(ctmp1) .lt. dpico .and. cdabs(ctmp3) .gt. dpico ) then
                         ctmp1 = ctmp2 / ctmp3
                     end if
                  end do
                  if (np .le. ndrude+nmats) cycle
                  ifff = np - ndrude - nmats
                  if (mfff1d(ifff) .eq. 1) cycle
                  ctheta(nm,np,nj,ni) = dble(mfff1d(ifff) - 1) * ctmp1
                  dtmp1 = max(dtmp1, cdabs(ctheta(nm,np,nj,ni)))
                  dtmp2 = min(dtmp2, cdabs(ctheta(nm,np,nj,ni)))
               end do
            end do
         end do
      end do
      write(6,*)
      write(6,*)'evaluatepara: max(ctheta) = ', dtmp1
      write(6,*)'evaluatepara: min(ctheta) = ', dtmp2
      call flush(6)
  end if
!
else
! offcor case need to be checked
  do ni=1,nsgn
    do nm=1,nspin
      do nk=1,nvar2
        if (nvar2 .eq. norbs2) then
          call jmod(nk, norbs, iorbs, iorbs2)
          nq = (iorbs - 1) * norbs + iorbs2   ! transpose of nk
          if (ni .eq. 1) then
            n2 = nq                           ! take the transpose for sigma=+
            iorbsa = iorbs2
            iorbsb = iorbs
          else
            n2 = nk
            iorbsa = iorbs
            iorbsb = iorbs2
          end if
        else
          iorbs  = nk
          iorbs2 = nk
          iorbsa = nk
          iorbsb = nk
          nq = nk
          n2 = nk
        end if
        do nj=1,nalf
          do np=1,ncor
            if (np .le. ndrude) then
              write(6,*)
              write(6,*)' error! unexpected pole for spec. func. ', np, ndrude, ispectral
              stop
            else if (np .le. ndrude+nmats .and. nmats .gt. 0) then
              imats = np - ndrude
              if (pfdjob) then
                ! add this part later
                write(6,*)'evaluatepara: error! pfdjob for ispectral != 0 not available yet'
                stop
              else if (psdjob) then
                ! add this part later
                write(6,*)'evaluatepara: error! psdjob for ispectral != 0 not available yet'
                stop
              else
                cgama(nk,nm,np,nj,ni) = dcmplx(-dmatsu(imats,nj), 0.d0) / hbar
                ctmp1 = dcmplx(0.d0, dpm(ni) * dmatsu(imats,nj))
                call specdensfunc(nj, nm, iorbsa, iorbsb, ctmp1, ctmp2)
                cb(nk,nm,np,nj,ni) = dcmplx(0.d0, -dinvbeta(nj)) * ctmp2 / hbar**2
!
                ctmp1 = dcmplx(0.d0, dpm(nsgn-ni+1) * dmatsu(imats,nj))
                call specdensfunc(nj, nm, iorbsb, iorbsa, ctmp1, ctmp2)
                cd(nk,nm,np,nj,ni) = dcmplx(0.d0,  dinvbeta(nj)) * dconjg(ctmp2) / hbar**2
              end if
              dimpf(nk,nm,imats,nj,ni) = dsqrt(cdabs(cb(nk,nm,np,nj,ni))) / cdabs(cgama(nk,nm,np,nj,ni))
            else if (np .gt. ndrude+nmats .and. nfreq .gt. 0) then
              ifreq = np - ndrude - nmats
              cb(nk,nm,np,nj,ni)    = dkln(ifreq) * spectral(ni,nj,nm,iorbs,iorbs2,ifreq) *     &
                                      focc(ni,nj,nm,ifreq) / (2.d0 * pi) / hbar**2
              cd(nk,nm,np,nj,ni)    = dkln(ifreq) * spectral(ni,nj,nm,iorbs,iorbs2,ifreq) *     &
                                      (1.d0 - focc(ni,nj,nm,ifreq)) / (2.d0 * pi) / hbar**2
              cgama(nk,nm,np,nj,ni) = dcmplx( -yshift(nj), dpm(ni) * wkln(ifreq) ) / hbar
            end if
          end do
        end do
      end do
    end do
  end do
end if
!
!write(6,*)' output cb parameters '
!write(6,*)' isgn, ialf, icor, ispin, iorbs2, cb '
!do ni=1,nsgn
!   do nj=1,nalf
!      do nk=1,ncor
!         do nm=1,nspin
!            do nn=1,nvar2
!               write(6,1221)ni, nj, nk, nm, nn, dble(cb(nn,nm,nk,nj,ni)), dimag(cb(nn,nm,nk,nj,ni))
!            end do
!         end do
!      end do
!   end do
!end do
!call flush(6)
!
999 continue
!
call cpu_time(cpu2)
!
nq = 0
do isgn=1,nsgn
  do ispin=1,nspin
    do iorbs=1,norbs
      if (offcor) then
         iorbs2 = (iorbs - 1) * norbs + iorbs
      else
         iorbs2 = iorbs
      end if
      do ialf=1,nalf
        do imats=1,ncor
          nq = nq + 1
          cgamma(nq) = cgama(iorbs2,ispin,imats,ialf,isgn) ! cgama is constant, irrespective of any iorbsp in (iorbs,iorbsp) -- obsolete comment?
        end do
      end do
    end do
  end do
end do
if (nq .ne. nvar) then
  write(6,*)
  write(6,*)' error when making cgamma ', nq, nvar
  stop
end if
if (.not. lfermi_exact) then
   write(6,*)
   if (pfdjob) then
      write(6,*)'evaluatepara: fermi function is approximated by PFD expansion '
   else if (psdjob) then
      write(6,*)'evaluatepara: fermi function is approximated by PSD expansion '
   else if (psdlor) then
      write(6,*)'evaluatepara: fermi function is approximated by PSD_LOR expansion '
   else if (psdfff) then
      write(6,*)'evaluatepara: fermi function is approximated by PSD_FFF expansion '
   else
      write(6,*)'evaluatepara: fermi function is approximated by Matsubara expansion '
   end if
   call flush(6)
end if
!
! print parameters
!
write(6,*)
write(6,*)' Matsubara/Pade/Prony parameters '
call flush(6)
isgn   = chksgn
ialf   = chkalf
ispin  = chkspin
if (offcor) then
  iorbs = (chkorbsp - 1) * norbs + chkorbs
else
  iorbs = chkorbs
end if
if (nmats .gt. 0) then
  dsum = 0.d0
  do imats=1,nmats
    dsum = dsum + dimpf(iorbs,ispin,imats,ialf,isgn)
  end do
  write(6,*)
  write(6,201)isgn, ialf, iorbs, ispin
  do np=1,nmats
    nq = np + ndrude
    write(6,203)np, dble(cb(iorbs,ispin,nq,ialf,isgn)), dimag(cb(iorbs,ispin,nq,ialf,isgn)),  &
                dble (cgama(iorbs,ispin,nq,ialf,isgn)) * hbar,                                &
                dimag(cgama(iorbs,ispin,nq,ialf,isgn)) * hbar,                                &
                dimpf(iorbs,ispin,np,ialf,isgn) / dsum * 1.d2
  end do
end if
201 format(' isgn=', I2, 1x, ' nalf=', I2, 1x, ' iorbs=', I2, 1x, ' ispin=', I2)
202 format(' cb', 2(2x, f14.8), ' cgamma', 2(1x, f14.8), ' impact', 1x, f13.6, ' %', f6.2)
203 format(' mats=', I4, 2x, ' cb', 2(1x, f10.4), ' cgamma', 2(1x, f10.4), ' impact %', f8.3)
!
if (nfreq .gt. 0) then
  write(6,*)
  write(6,*)' Frequency parameters: cb, cd, cgamma '
  do imats=1,nfreq
    nq = imats + ndrude + nmats
    write(6,1004)imats, dble(cb(iorbs,ispin,nq,ialf,isgn)), dimag(cb(iorbs,ispin,nq,ialf,isgn)), &
                        dble(cd(iorbs,ispin,nq,ialf,isgn)), dimag(cd(iorbs,ispin,nq,ialf,isgn)), &
                        dble (cgama(iorbs,ispin,nq,ialf,isgn)) * hbar,                           &
                        dimag(cgama(iorbs,ispin,nq,ialf,isgn)) * hbar
  end do
end if
call flush(6)
1004 format(I4, 1x, ' cb', 2(2x, f14.8), ' cd', 2(2x, f14.8), ' cgamma', 2(2x, f14.8))
!
! Time-dependent external field applied on the leads
!
engyshift = engyshift / hbar
!
fieldtype  = 0
!lreadomega = .false.
lreadomega = .true.
leshift_dc = .false.
t_off      = 1.d99   ! careful! to calculate steady current under fixed voltage
                     ! t_off should be larger than the below value used to calculate eleadinfty
rewind(5)
read(5, field, end=103)
103 continue
!
if (malf .gt. 2) then
   if (fieldtype .ne. 0 .and. fieldtype .ne. 1 .and. fieldtype .ne. 2 .and. fieldtype .ne. -1) then
      write(6,*)'evaluatepara: multi-lead unavailable for the following fieldtype'
      write(6,*)'evaluatepara: fieldtype = ', fieldtype, ' nalf = ', malf
      stop
   end if
end if
if (t_off .lt. 1.d10) then
   write(6,*)'evaluatepara: t_off = ', t_off
end if
!
! read tchar
! for fieldtype = 0, tchar is the characteristic decaying time
! for fieldtype = 1, tchar is the period of sinusoidal voltage
! for fieldtype = 2, tchar is the characteristic decaying time
!
if (fieldtype .eq. 0) then
   read(5,*,iostat=istat) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
   if (istat .ne. 0) then
      write(6,*)'evaluatepara: error reading tchar from input '
      stop
   end if
else if (fieldtype .eq. 1) then
   if (lreadomega) then ! read in as eV or meV (energy), and converted to frequency (divided by hbar)
      read(5,*,iostat=istat) ((omegas(ni,nj), nj=1,nspin), ni=1,malf)
      if (istat .ne. 0) then
          write(6,*)'evaluatepara: error when reading omegas '
          stop
      end if
      omegas(1:malf,1:nspin) = omegas(1:malf,1:nspin) / hbar
      do ni=1,malf
         do nj=1,nspin
            if (dabs(omegas(ni,nj)) .le. dpico) then
               tchar(ni,nj) = 1.d12
            else
               tchar(ni,nj) = 2.d0 * pi / omegas(ni,nj)
            end if
         end do
      end do
   else
      read(5,*,iostat=istat) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
      if (istat .ne. 0) then
          write(6,*)'evaluatepara: error when reading tchar '
          stop
      end if
      do ni=1,malf
         do nj=1,nspin
            if (dabs(tchar(ni,nj)) .le. dpico) then
               omegas(ni,nj) = 1.d12
            else
               omegas(ni,nj) = 2.d0 * pi / tchar(ni,nj)
            end if
         end do
      end do
   end if
   write(6,*)'evaluatepara: period of sinusoidal function is T_period = '
   write(6,*) ((2.d0 * pi / omegas(ni,nj), nj=1,nspin), ni=1,malf)
   if (leshift_dc) then
       read(5,*,iostat=istat) ((amps_dc(ni,nj), nj=1,nspin), ni=1,malf)
       if (istat .ne. 0) then
           write(6,*)'evaluatepara: error when reading amps_dc '
           stop
       end if
       amps_dc(1:malf,1:nspin) = amps_dc(1:malf,1:nspin) / hbar
       write(6,*)'evaluatepara: dc amplitudes '
       write(6,*) ((amps_dc(ni,nj)*hbar, nj=1,nspin), ni=1,malf)
   end if
   call flush(6)
else if (fieldtype .eq. 2) then
   read(5,*,iostat=istat) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
   if (istat .ne. 0) then
      write(6,*)'evaluatepara: error reading tchar from input '
      stop
   end if
end if
!
write(6,*)
if (fieldtype .eq. 0) then
   amps(1:malf,1:nspin) = engyshift(1:malf,1:nspin)
   write(6,*)'evaluatepara: exponential voltages are applied, parameters are '
   write(6,*)'evaluatepara: amplitude, amps(nspin, malf) '
   write(6,*) ((amps(ni,nj) * hbar, nj=1,nspin), ni=1,malf)
   write(6,*)'evaluatepara: characteristic decay time, tchar(nspin, malf) '
   write(6,*) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
   call flush(6)
!
else if (fieldtype .eq. 1) then
   amps(1:malf,1:nspin) = engyshift(1:malf,1:nspin)
   write(6,*)'evaluatepara: sine voltages are applied, parameters are '
   write(6,*)'evaluatepara: amplitude, amps(nspin, malf) '
   write(6,*) ((amps(ni,nj) * hbar, nj=1,nspin), ni=1,malf)
   write(6,*)'evaluatepara: period time of sinusoidal voltage, tchar(nspin, malf) '
   write(6,*) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
   write(6,*)'evaluatapara: frequency of sinusoidal voltage, omegas(nspin, malf) '
   write(6,*) ((omegas(ni,nj) * hbar, nj=1,nspin), ni=1,malf)
   call flush(6)
!
else if (fieldtype .eq. 2) then
   amps(1:malf,1:nspin) = engyshift(1:malf,1:nspin)
   write(6,*)'evaluatepara: pulsed energy shifts are used on top of    '
   write(6,*)'evaluatepara: exponential voltages with parameters being '
   write(6,*)'evaluatepara: amplitude, amps(nspin, malf) '
   write(6,*) ((amps(ni,nj) * hbar, nj=1,nspin), ni=1,malf)
   write(6,*)'evaluatepara: characteristic decay time, tchar(nspin, malf) '
   write(6,*) ((tchar(ni,nj), nj=1,nspin), ni=1,malf)
   call flush(6)
!
else if (fieldtype .eq. -1) then
   amps(1:malf,1:nspin) = engyshift(1:malf,1:nspin)
   write(6,*)'evaluatepara: delta-type voltages are applied, parameters are '
   write(6,*)'evaluatepara: amplitude, amps(nspin, malf) '
   write(6,*) ((amps(ni,nj) * hbar, nj=1,nspin), ni=1,malf)
   call flush(6)
else
   write(6,*)' error! unknown external voltage type ', fieldtype
   stop
end if
!
! calculate lead energy shift as time goes to infinity
!
do ialf=1,nalf
   do ispin=1,nspin
      call getengyshift(1.d10, ialf, ispin, eleadinfty(ialf,ispin))
      write(6,*)'evaluatepara: lead energy shift at t=infty: ialf, ispin, eshift'
      write(6,*)'evaluatepara: ', ialf, ispin, eleadinfty(ialf,ispin) * hbar
   end do
end do
call flush(6)
!
lspin3d     = .true.
lbfield3d   = .false.
lbsite      = .false.
rewind(5)
read(5, spin3d, end=213)
213 continue
if (.not. lspin3d) then
    write(6,*)
    write(6,*)'evaluatepara: Warning! lspin3d=F used, spin-operators are absent '
    call flush(6)
end if
!
lzfs = .false.
d_xx = 0.d0
d_yy = 0.d0
d_zz = 0.d0
d_xx_1 = 0.d0
d_yy_1 = 0.d0
d_zz_1 = 0.d0
d_xx_2 = 0.d0
d_yy_2 = 0.d0
d_zz_2 = 0.d0
rewind(5)
read(5, zfs, end=214)
214 continue
!
! Time-dependent external field applied directly to the system
!
ldfield_hub = .false.
dfieldtype  = 0
tson        = 0.d0
tsoff       = 1.d10
dedot1      = 0.d0
dedot2      = 0.d0
gamadot     = 0.d0
rewind(5)
read(5, dfield_hubbard, end=412)
412 continue
if (ldfield_hub) then
   if (dfieldtype .eq. 1) then
      write(6,*)'evaluatepara: parameters for field exerted on Hubbard system'
      write(6,*)'evaluatepara: dfieldtype = ', dfieldtype
      write(6,*)'evalautepara: dedot1  = ', dedot1
      if (norbs > 1) then
         write(6,*)'evalautepara: dedot2  = ', dedot2
      end if
      write(6,*)'evaluatepara: tson    = ', tson
      write(6,*)'evaluatepara: tsoff   = ', tsoff
      write(6,*)'evaluatepara: gamadot = ', gamadot
      call flush(6)
      dedot1  = dedot1 / hbar
      dedot2  = dedot2 / hbar
      gamadot = gamadot / hbar
      if (tson .gt. tsoff) then
         write(6,*)'evaluatepara: error! tson > tsoff found ', tson, tsoff
         stop
      end if
   else
      write(6,*)'evaluatepara: error! unknown dfieldtype = ', dfieldtype
      write(6,*)'evaluatepara:        for Hubbard system   '
      stop
   end if
   goto 413
end if
!
gatetype = 0
gpara1 = 0.d0
gpara2 = 0.d0
gpara3 = 0.d0
gpara4 = 0.d0
rewind(5)
read(5, vgate, end=110)
110 continue
write(6,*)
write(6,*)'evaluatepara: gatetype = ', gatetype
write(6,*)'evaluatepara: gpara1 = ', gpara1
write(6,*)'evaluatepara: gpara2 = ', gpara2
write(6,*)'evaluatepara: gpara3 = ', gpara3
write(6,*)'evaluatepara: gpara4 = ', gpara4
call flush(6)
if (gatetype .eq. 1) then
   gpara1 = gpara1 / hbar
   gpara2 = gpara2 / hbar
   gpara3 = gpara3 / hbar
end if
!
dfieldtype  = 0
forcesteady = .false.
if (nspin .eq. 1 .and. norbs .eq. 1) then
  dedot      = 0.d0
  dfieldtype = 0
  aD         = 1.d-20
  wD         = 1.d0
  rewind(5)
  read(5, dfield2, end=104)
  104 continue
  dedot = dedot / hbar
  wD    = wD / hbar
  if (dabs(dedot) .ge. dnano) then
    write(6,*)
    if (dfieldtype .eq. 0) then
      write(6,*)' exponential dot level shift '
      write(6,*)' dedot = ', dedot*hbar, ' aD = ', aD
    else if (dfieldtype .eq. 1) then
      write(6,*)' sinusoidal  dot level shift '
      write(6,*)' dedot = ', dedot*hbar, ' wD = ', wD*hbar
    else
      write(6,*)' error! unknown dfieldtype read ', dfieldtype
      stop
    end if
    call flush(6)
  end if
else if (nspin .eq. 2 .and. norbs .eq. 1) then
   egate = 0.d0
   ton_gate = 0.1d0
   toff_gate = 0.2d0
   rewind(5)
   read(5, dfield5, end=301)
   301 continue
   egate = egate / hbar
   write(6,*)
   write(6,*)'evaluatepara: egate = ', egate*hbar
   write(6,*)'evaluatepara: ton_gate = ', ton_gate
   write(6,*)'evaluatepara: toff_gate = ', toff_gate
   call flush(6)
!
   if (lspin3d) then
      dfieldtype  = 0
      dedot1      = 0.d0
      dedot2      = 0.d0
      wdot1       = 1.d5
      wdot2       = 1.d5
      tdot12      = 0.d0
      jdot12      = 0.d0
      rewind(5)
      read(5, dfield3, end=216)
216 continue
      dedot1 = dedot1 / hbar
! dfieldtype = 1, wdot used as t_off on the dot (unit in time)
! otherwise, wdot is scaled by hbar
      if ( .not. dfieldtype .eq. 1) then
         wdot1 = wdot1 / hbar
      end if
      write(6,*)
      if (dfieldtype .eq. 0) then
         write(6,*)'no external field is applied to the dot'
      else if (dfieldtype .eq. 1) then
         write(6,*)'step type external field is applied to the dot'
         write(6,*)'dedot =', dedot1 * hbar, 'turn-off time =', wdot1
      else if (dfieldtype .eq. 2) then
         write(6,*)'cos(wt)-type filed is applied to the dot'
         write(6,*)' dedot1 =', dedot1 * hbar, ' wdot1 =', wdot1 * hbar
      else
         write(6,*)'error! unknown dfieldtype read', dfieldtype
         stop
      end if
   else
!
      dfieldtype = 0
      tflip      = 0.d0
      wflip      = 1.d0
      rewind(5)
      read(5, dfield1, end=105)
105 continue
      tflip = tflip / hbar
      wflip = wflip / hbar
      write(6,*)
      if (dfieldtype .eq. 0) then
         write(6,*)' no external field is applied to the dot '
      else if (dfieldtype .eq. 1 .and. dabs(tflip) .ge. dnano) then
         write(6,*)' sin(wt)-type spin-flip field is applied to the dot '
         write(6,*)' tflip = ', tflip*hbar, ' wflip = ', wflip*hbar
      else if (dfieldtype .eq. 2) then
         write(6,*)' cos(wt)-type spin-flip field is applied to the dot '
         write(6,*)' tflip = tupdn (see input file) ', ' wflip = ', wflip*hbar
      else if (dfieldtype .eq. 3 .and. dabs(tflip) .ge. dnano) then
         write(6,*)' exp(iwt)-type spin-flip field is applied to the dot '
         write(6,*)' tflip = ', tflip*hbar, ' wflip = ', wflip*hbar
      else
         write(6,*)' error! unknown dfieldtype read ', dfieldtype
         stop
      end if
   end if
  call flush(6)
else if (nspin .eq. 2 .and. norbs .eq. 2) then
  dfieldtype = 0
  dedot1     = 0.d0
  dedot2     = 0.d0
  wdot1      = 1.d5
  wdot2      = 1.d5
  tdot12     = 0.d0
  jdot12     = 0.d0
  rewind(5)
  read(5, dfield3, end=116)
  116 continue
  dedot1 = dedot1 / hbar
  dedot2 = dedot2 / hbar
! wdot used as turn-off time for dfieldtype == 1
  if (.not. dfieldtype .eq. 1) then
     wdot1  = wdot1 / hbar
     wdot2  = wdot2 / hbar
  end if
  tdot12 = tdot12 / hbar
  jdot12 = jdot12 / hbar
  write(6,*)
  if (dfieldtype .eq. 0) then
    write(6,*)' no field is applied on the dot '
  else if (dfieldtype .eq. 1) then
    write(6,*)' step-type external field is applied on the dot:'
    write(6,*)' dedot1 =', dedot1 * hbar, 'turn-off time =', wdot1
    write(6,*)' dedot2 =', dedot2 * hbar, 'turn-off time =', wdot2
    write(6,*)' interdot coupling tdot12 is turned on and kept as constant '
    write(6,*)' coupling strength tdot12 = ', tdot12*hbar
    write(6,*)' interdot exchange jdot12 is turned on and kept as constant '
    write(6,*)' coupling strength jdot12 = ', jdot12*hbar
  else if (dfieldtype .eq. 2) then
    write(6,*)' cos(wt)-type field is applied to the two dots '
    write(6,*)' dedot1 = ', dedot1*hbar, ' wdot1 = ', wdot1*hbar
    write(6,*)' dedot2 = ', dedot2*hbar, ' wdot2 = ', wdot2*hbar
    write(6,*)' interdot coupling tdot12 is turned on and kept as constant '
    write(6,*)' coupling strength tdot12 = ', tdot12*hbar
    write(6,*)' interdot exchange jdot12 is turned on and kept as constant '
    write(6,*)' coupling strength jdot12 = ', jdot12*hbar
  else
    write(6,*)' error! unknown dfieldtype read ', dfieldtype
    stop
  end if
  call flush(6)
else if (nspin .eq. 1 .and. norbs .eq. 2) then
  dfieldtype = 0
  dedot1     = 0.d0
  dedot2     = 0.d0
  wdot1      = 1.d0
  wdot2      = 1.d0
  tson       = 0.d0
  tsoff      = 1.d0
  pha        = 0.d0 ! phase coefficent for levels' oscillator
  rewind(5)
  read(5, dfield4, end=117)
  117 continue
  dedot1 = dedot1 / hbar
  dedot2 = dedot2 / hbar
  wdot1  = wdot1 / hbar
  wdot2  = wdot2 / hbar
  pha    = pha / hbar
  write(6,*)
  if (dfieldtype .eq. 0) then
    write(6,*)' no external field is applied to the dot '
  else if (dfieldtype .eq. 2) then
    write(6,*)' cos(wt)-type field is applied to the two dots '
    write(6,*)' dedot1 = ', dedot1*hbar, ' wdot1 = ', wdot1*hbar
    write(6,*)' dedot2 = ', dedot2*hbar, ' wdot2 = ', wdot2*hbar
  else if (dfieldtype .eq. 4) then
    write(6,*)' step on_off field is applied to the two dots '
    write(6,*)' dedot1 = ', dedot1*hbar, ' dedot2 = ', dedot2*hbar
    write(6,*)' t_on   = ', tson,        ' t_off  = ', tsoff
    write(6,*)' pha    = ', pha*hbar
  else
    write(6,*)' error! unknown dfieldtype read ', dfieldtype
    stop
  end if
  call flush(6)
end if
!
413 continue
!
! Step-like magnetic field for transient dynamics
!
lbfield = .false.
dbfield = 0.d0
tbfield = 1.d-20
rewind(5)
read(5, bfield, end=121)
121 continue
!
dbfield = dbfield / hbar
if (lbfield) then
   write(6,*)
   write(6,*)'evaluatepara: effective B-field for transient dynamics '
   write(6,*)'evaluatepara: H'' = - dbfield * Sz                     '
   write(6,*)'evaluatepara: dbfield = ', dbfield * hbar
   write(6,*)'evaluatepara: tbfield = ', tbfield
   call flush(6)
end if
!
! calc approx corr. func. C_{iorbs,ispin,ialf,isgn}(t), and compare to
!      exact values, and then output
!
if (skipcorr) then
  write(6,*)
  write(6,*)' checking of correlation function is skipped '
  call flush(6)
  goto 518
end if
ltmp1 = .false.
if (ispectral .eq. 0 .and. ltmp1) then
  nmatmp = 50000  ! Warning! may not be enough for very low temperature!
  if (nmatmp .le. nmats) then
    write(6,*)
    write(6,*)' error! nmatmp <= nmats found! ', nmatmp, nmats
    stop
  end if
  allocate(dmatmp(nmatmp), zbmtmp(nmatmp), STAT=istat)
else
  nw0tmp = 1500
  wmin0  = -4.d2
  wmax0  =  4.d2
  yshft0 = 1.d10
  do ni=1,nalf
     yshft0 = min(yshft0, 5.d-1 * dinvbeta(ni) * pi)   ! yshft0 should be adjusted so that all other possible
                                                       !  poles are excluded between y=0 and y=yshft0
  end do
!
  if (lband) then
     do ni=1,nalf
        do nj=1,nspin
           do nk=1,nband
              yshft0 = min(yshft0, band_width(nk,nj,ni) * 5.d-1 * hbar)
           end do
        end do
     end do
  else
     if (ispectral .eq. 0) then
        do ni=1,nalf
           yshft0 = min(yshft0, bandwidth(ni) * 5.d-1)
        end do
     end if
  end if
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
end if
!
ttmin  = 0.d0
ttmax  = 2.d0
dtt    = 1.d-3
tt     = ttmin
isgn   = chksgn
ialf   = chkalf
ispin  = chkspin
iorbs  = chkorbs
iorbsp = chkorbsp
if (offcor) then
  iorbs2 = (iorbsp - 1) * norbs + iorbs
  if (isgn .eq. 1) then
    iorbsa = iorbsp
    iorbsb = iorbs
        n2 = (iorbs - 1) * norbs + iorbsp
  else
    iorbsa = iorbs
    iorbsb = iorbsp
        n2 = iorbs2
  end if
else
  iorbs2 = iorbs
  iorbsa = iorbs
  iorbsb = iorbs
      n2 = iorbs
end if
write(6,*)
write(6,*)' C_{iorbs,iorbsp,ispin,ialf,isgn}(t) written to <corr.data> '
if (chkcd) then
  write(6,*)' checking C^{\bar{\sigma}}(t) '
else
  write(6,*)' checking C^\sigma(t)         '
end if
write(6,204)iorbs, iorbsp, ispin, ialf, isgn
call flush(6)
204 format(' iorbs=', I2, 2x, ' iorbsp=', I2, 2x, ' ispin=', I2, 2x, ' ialf=', I2, 2x, ' isgn=', I2)
!
if (ispectral .eq. 0 .and. ltmp1) then
  do imats=1,nmatmp
    dmatmp(imats) = (2.d0 * dble(imats) - 1.d0) * pi * dinvbeta(ialf)
    call specdensfunc(ialf, ispin, iorbsa, iorbsb, dcmplx(0.d0,dpm(isgn)*dmatmp(imats)), ctmp2)
    zbmtmp(imats) = dcmplx(0.d0, -dinvbeta(ialf)) * ctmp2 / hbar**2
  end do
  ctmp3 = dcmplx( -gammap(iorbs2,ispin,ialf), dpm(isgn) * dwp(iorbs2,ispin,ialf) ) / hbar
  ctmp4 = zxip(n2,ispin,ialf) / gammap(iorbs2,ispin,ialf) / (1.d0 + cdexp(dbeta(ialf) *     &
          dcmplx(dpm(isgn) * dwp(iorbs2,ispin,ialf), gammap(iorbs2,ispin,ialf))) ) / hbar**2
else
  do ifreq=1,nwatmp
    ctmp1 = dcmplx(wkltmp(ifreq), dpm(isgn) * yshft0)
    call specdensfunc(ialf, ispin, iorbsa, iorbsb, ctmp1, ctmp2)
    call foccfunc(isgn, ialf, ctmp1, ctmp3)
    zw0tmp(ifreq) = dkltmp(ifreq) / (2.d0 * pi) * ctmp2 * ctmp3 / hbar**2
  end do
end if
!
dtmp1 = 0.d0
dtmp2 = 0.d0
open(unit=22, file='corr.data', status='unknown')
rewind(22)
15 continue
ctmp1 = czero
ctmp5 = czero
jsgn  = nsgn - isgn + 1
if (psdfff) then
    do imats=1,ncor
       if (imats .le. ndrude+nmats) then
          ctmp1 = ctmp1 + cb(iorbs2,ispin,imats,ialf,isgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,isgn) * tt)
          ctmp5 = ctmp5 + cd(iorbs2,ispin,imats,ialf,jsgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,jsgn) * tt)
       else
          ifff = mfff1d(imats - ndrude - nmats) - 1
          dtmp3 = 1.d0
          if (ifff .gt. 0) dtmp3 = tt**ifff
          ctmp1 = ctmp1 + cb(iorbs2,ispin,imats,ialf,isgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,isgn) * tt) * dtmp3
          ctmp5 = ctmp5 + cd(iorbs2,ispin,imats,ialf,jsgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,jsgn) * tt) * dtmp3
       end if
    end do
else
    do imats=1,ncor
      ctmp1 = ctmp1 + cb(iorbs2,ispin,imats,ialf,isgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,isgn) * tt)
      ctmp5 = ctmp5 + cd(iorbs2,ispin,imats,ialf,jsgn) * cdexp(cgama(iorbs2,ispin,imats,ialf,jsgn) * tt)
    end do
end if
ctmp2 = czero
if (ispectral .eq. 0 .and. ltmp1) then
  do imats=1,nmatmp
    ctmp2 = ctmp2 + zbmtmp(imats) * dexp(-dmatmp(imats) * tt / hbar)
  end do
  ctmp2 = ctmp2 + ctmp4 * cdexp(ctmp3 * tt)
else
  do ifreq=1,nwatmp
     ctmp2 = ctmp2 + zw0tmp(ifreq) * cdexp( dcmplx(0.d0, dpm(isgn)) * dcmplx(wkltmp(ifreq),  &
             dpm(isgn) * yshft0) * tt / hbar )
  end do
  ctmp6 = dconjg(ctmp2)
end if
if (.not. chkcd) then
!  dtmp1 = dtmp1 + cdabs(ctmp1)
  dtmp1 = dtmp1 + cdabs(ctmp2)
  dtmp2 = dtmp2 + cdabs(ctmp1 - ctmp2)
  write(22,520)tt, dble(ctmp1), dimag(ctmp1), dble(ctmp2), dimag(ctmp2)
else
!  dtmp1 = dtmp1 + cdabs(ctmp5)
  dtmp1 = dtmp1 + cdabs(ctmp6)
  dtmp2 = dtmp2 + cdabs(ctmp5 - ctmp6)
  write(22,520)tt, dble(ctmp5), dimag(ctmp5), dble(ctmp6), dimag(ctmp6)
  end if
call flush(22)
520 format(f12.6, 2x, 10(e13.4e3, 2x))
tt = tt + dtt
if (tt .le. ttmax) goto 15
close(22)
write(6,*)
write(6,"('relative discrepancy to exact corrfunc ', 2x, f8.3, 2x, ' %')") dtmp2/dtmp1*1.d2
write(6,*)'exact,  exact - approximate ', dtmp1, dtmp2
call flush(6)
if (ispectral .eq. 0) then
  deallocate(dmatmp, zbmtmp, STAT=istat)
else
  deallocate(w0tmp, w1tmp, wkltmp, dkltmp, zw0tmp, STAT=istat)
end if
!
lpcontour = .false.
nxpmin    = -200
nxpmax    =  200
nypmin    = -50
nypmax    =  50
dpgrid    =  5.d-1
tour      =  1.d0
rewind(5)
read(5, pcontour, end=120)
120 continue
if (lpcontour) then
   call plotcorr
end if
518 continue
if (chkcorr) then
  write(6,*)
  write(6,*)' chkcorr=.true. found, program terminates after checking '
  write(6,*)' correlation functional expansion.                       '
  stop
end if
!
! Scale all ADOs to dimensionless matrices
!       rho_new = Times_k cdabs(cb_k)**(-0.5) rho_old
!
if (lscale) then
   allocate(dbsqrt(norbs,nspin,ncor,nalf,nsgn), dbinvsqrt(norbs,nspin,ncor,nalf,nsgn), STAT=istat)
   do isgn=1,nsgn
      do ialf=1,nalf
         do imats=1,ncor
            do ispin=1,nspin
               do iorbs=1,norbs
                  if (offcor) then
                     iorbs2 = (iorbs - 1) * norbs + iorbs
                  else
                     iorbs2 = iorbs
                  end if
                  dtmp1 = cdabs(cb(iorbs2,ispin,imats,ialf,isgn))
                  !if (dtmp1 .le. dpico) then
                  !   write(6,*)
                  !   write(6,*)' error! problem with scaling coefficients '
                  !   write(6,*)' abs(cb) = ', dtmp1, ' with               '
                  !   write(6,*)' isgn, ialf, imats, ispin, iorbs2 = ', isgn, ialf, imats, ispin, iorbs2
                  !   stop
                  !end if
!
!                  if (psdfff .and. imats .gt. ndrude+nmats) then
!                      ntmp1 = mfff1d(imats - ndrude - nmats) - 1
!                      dtmp1 = cdabs(cb(iorbs2,ispin,imats-ntmp1,ialf,isgn))
!                  end if
                  dbsqrt(iorbs,ispin,imats,ialf,isgn) = dsqrt(dtmp1)
                  if (dtmp1 .gt. dpico) then
                     dbinvsqrt(iorbs,ispin,imats,ialf,isgn) = 1.d0 / dbsqrt(iorbs,ispin,imats,ialf,isgn)
                  else
                     dbinvsqrt(iorbs,ispin,imats,ialf,isgn) = 0.d0
                  end if
               end do
!               do nk=1,nvar2
!                  if (offcor) then
!                     call jmod(nk, norbs, iorbs, iorbs2)
!                  else
!                     iorbs = nk
!                  end if
!                  dtmp1 = dbinvsqrt(iorbs,ispin,imats,ialf,isgn)
!                  cb(nk,ispin,imats,ialf,isgn) = cb(nk,ispin,imats,ialf,isgn) * dtmp1
!                  cd(nk,ispin,imats,ialf,isgn) = cd(nk,ispin,imats,ialf,isgn) * dtmp1
!               end do
            end do
         end do
      end do
   end do
!
   do isgn=1,nsgn
      do ialf=1,nalf
         do imats=1,ncor
            do ispin=1,nspin
               do nk=1,nvar2
                  if (offcor) then
                     call jmod(nk, norbs, iorbs, iorbs2)
                  else
                     iorbs = nk
                  end if
                  dtmp1 = dbinvsqrt(iorbs,ispin,imats,ialf,isgn)
                  cb(nk,ispin,imats,ialf,isgn) = cb(nk,ispin,imats,ialf,isgn) * dtmp1
                  cd(nk,ispin,imats,ialf,isgn) = cd(nk,ispin,imats,ialf,isgn) * dtmp1
               end do
            end do
         end do
      end do
   end do
!
   if (psdfff) then
       ctheta = czero
       dtmp1 = 0.d0
       dtmp2 = 1.d20
       write(6,*)
       write(6,*)'evaluatepara: print cb and gama (lscale = T) '
       nq = 0
       do ni=1,nsgn
          do nm=1,nspin
             do nj=1,nalf
                do np=1,ncor
                   ctmp1 = czero
                   do nk=1,nvar2
                      nq = nq + 1
                      write(6,1101)nq, dble(cb(nk,nm,np,nj,ni)), dimag(cb(nk,nm,np,nj,ni)),  &
                                   dble(cgama(nk,nm,np,nj,ni)), dimag(cgama(nk,nm,np,nj,ni))
                      ctmp2 = cb(nk,nm,np  ,nj,ni)
                      ctmp3 = cb(nk,nm,np-1,nj,ni)
                      if ( cdabs(ctmp1) .lt. dpico .and. cdabs(ctmp3) .gt. dpico ) then
                          ctmp1 = ctmp2 / ctmp3
                      end if
                   end do
                   if (np .le. ndrude+nmats) cycle
                   ifff = np - ndrude - nmats
                   if (mfff1d(ifff) .eq. 1) cycle
                   ctheta(nm,np,nj,ni) = dble(mfff1d(ifff) - 1) * ctmp1
                   dtmp1 = max(dtmp1, cdabs(ctheta(nm,np,nj,ni)))
                   dtmp2 = min(dtmp2, cdabs(ctheta(nm,np,nj,ni)))
                end do
             end do
          end do
       end do
       write(6,*)
       write(6,*)'evaluatepara: max(ctheta) = ', dtmp1
       write(6,*)'evaluatepara: min(ctheta) = ', dtmp2
       call flush(6)
   end if
end if
!
write(6,*)' output cb parameters '
write(6,*)' isgn, ialf, icor, ispin, iorbs2, cb '
do ni=1,nsgn
   do nj=1,nalf
      do nk=1,ncor
         do nm=1,nspin
            do nn=1,nvar2
               write(6,1221)ni, nj, nk, nm, nn, dble(cb(nn,nm,nk,nj,ni)), dimag(cb(nn,nm,nk,nj,ni))
            end do
         end do
      end do
   end do
end do
call flush(6)
!
1101 format(I4, 1x, 4(f14.6, 1x))
1221 format(5(I3, 1x), 2(f14.8, 2x))
!
! analyze possible pre-screening of dissipative operators
!
allocate(indexdraw(nvar), STAT=istat)
jomit(1:nopr) = 0
nq = 0
np = 0
do isgn=1,nsgn
   do ispin=1,nspin
      do iorbs=1,norbs
         do ialf=1,nalf
            nq = nq + 1
            if (offcor) then
                dtmp1 = 0.d0
                do iorbsp=1,norbs
                   n2 = (iorbsp - 1) * norbs + iorbs
                   dtmp1 = max( dtmp1, cdabs(zxip(n2,ispin,ialf)) )
                end do
            else
                dtmp1 = cdabs(zxip(iorbs,ispin,ialf))
            end if
            if (dtmp1 .lt. dsmall) then
                jomit(nq) = 1
                np = np + 1
            end if
         end do
      end do
   end do
end do
write(6,*)
write(6,*)'evaluatepara: number of decoupled dissipative operators ', np
call flush(6)
!
np = 0
do isgn=1,nsgn
   do ispin=1,nspin
      do iorbs=1,norbs
         do ialf=1,nalf
            call look4operator(isgn, iorbs, ispin, ialf, iopr)
            do icor=1,ncor
               idraw = (iopr - 1) * ncor + icor
               if (lscreen .and. jomit(iopr) .eq. 1) then
                   indexdraw(idraw) = 0
               else
                   np = np + 1
                   indexdraw(idraw) = np
               end if
            end do
         end do
      end do
   end do
end do
write(6,*)'evaluatepara: number of decoupled drawers               ', nvar-np
call flush(6)
!
dcbmax = 0.d0
do isgn=1,nsgn
   do ialf=1,nalf
      do imats=1,ncor
         do ispin=1,nspin
            do nk=1,nvar2
               dcbmax = max(dcbmax, cdabs(cb(nk,ispin,imats,ialf,isgn)))
            end do
         end do
      end do
   end do
end do
!
end subroutine evaluatepara
