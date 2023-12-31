program heom
use matmod
use matmod_omp
use tmpmatmod
use sparsemod
use hbmod
use omp_lib
use sj_timer
implicit none
include '../include/sizes'
include '../include/common'
!
!------------------------------------------------------------------------------------
!   File Channels
! 
!  Unit        Filename             Format      Purpose
!   5          user-defined          yes        input file 
!   6          user-defined          yes        general output information
!  12          indextable.tmp        bin        scratch file for indextable
!  13          curr.data             yes        current versus time, j(t)
!  14          indextable.sav        bin        backuped indextable
!  15          variables(_td).sav     no        backuped unknowns rho
!  16          popu.data             yes        level population versus time, occ(t)
!  16(2)       coefindex.tmp         bin        scratch file for building coef-index
!  17          coefindex.data        bin        backuped coef-index
!  18          coefindex.data        bin        backuped coef-index to read              
!  20          diis.xtmp              no        scratch file for diis purpose 
!  21          diis.err               no        scratch file for diis purpose
!  22          corr.data             yes        correlation functions versus time, C(t)
!                                               (used by residue correction)
!  23          TAPE.resume           bin        backup file for resuming a td job
!  24          bicg.xtmp              no        scratch file for bicg purpose 
!  25          variables.gr(.st)      no        backuped unknowns rho for ground (steady) state
!  26          oprdiff.tmp            no        see <buildcoefindex.f90> (obsolete)
!  27          oprdiff.data           no        see <buildcoefindex.f90> (obsolete)
!  28          rhodiag.data          yes        diagonal element of rho(:,:,1) (real part)
!  29          rhotime.data          yes        rho(t) (complex)
!  30          ado.data              yes        first ADO rho(:,:,nfirst(itierchk)+inumchk) 
!  31          rho0plus.data          no        rho(0^+)
!  32          stoptd                yes        flag file to stop and save break point
!  33          poccdetail.data       yes        pocc(1:norbs,1:nspin) versus time
!  35          cwreal.data           yes        dble(corr(z))
!  36          cwimag.data           yes        imag(corr(z))
!  37          auxindex.data         bin        auxiliary index file
!  38          dos.data              yes        correlation function for a system level
!  39          TAPE.resume_cf         no        backup file for resuming a time-propagation job
!                                               for the evaluation of system correlation functions
!  40          hamil_sys.data         no        system Hamiltonian file
!  41          rsdm_sys.data          no        system reduced single-eletron density matrix
!  42          psd_lor.data          yes        fitting parameters for reservoir spectral function
!  43          chi_loc.data          yes        local impurity <Sz> and <Sz^2> with time
!  44          ado_big.data          yes        index information for ADOs larger than some value
!  45          filterindex.data       no        coef-index for filtered variables 
!  46          filterindex.tmp        no        scratch file for <filterindex.data>
!  47          local_temp            yes        evolution of local temperature 
!  48          rhocoo.tmp            bin        (row,col) positions of nonzero elements of rho 
!  49          rhoval.tmp            bin        temporary values of nonzero elements of rho
!  50          rho_spa(_td).sav      bin        values of all unknowns in sparse format
!  51          rho_spa.gr(.st)       bin        backup ground (steady) state unknowns in sparse format
!  52          sparse_index.data     bin        number of nonzeros and index of nonzeros 
!  53          sparse_info.data      bin        (row,col) position of nonzeros
!  54          sparse_index_cf.data  bin        number of nonzeros and index of nonzeros 
!  55          sparse_info_cf.data   bin        (row,col) position of nonzeros
!  56          sparse_index_td.data  bin        number of nonzeros and index of nonzeros 
!  57          sparse_info_td.data   bin        (row,col) position of nonzeros
!  61          sdot1.data             yes       3d spin moment on dot 1
!  62          sdot2.data             yes       3d spin moment on dot 2
!  63          spin_ddot.data         yes       spin properties versus time (t, S12, <S2>, C12, C12p) 
!  64          corr_hb.data           yes       correlation function for heat bath 
!  65          rho_hb.data           bin        auxiliary density matrices for coupled heat bath 
!  66          oprindex.tmp          bin        information only for lsimple=T
!  67          oprindex.data         bin        information only for lsimple=T
!  68          rhoini.data            no        initial random rho for save 
!  69          rho_spa.jac           bin        optimal rho_spa generated by Jacobi iteration
!  70          variables.jac          no        optimal rho     generated by Jacobi iteration
!  71          rho_spa.chk           bin        values of all unknowns in sparse format
!  72          engy.data             yes        engy_sys, engy_sb, engy_tot versus time
!  73          rho_spa.per           bin        rho0 in sparse form in steady-state perturbation (SSP)
!  74          rhs_spa.per           bin        (L0+dL)*rho0 in sparse form in SSP method
!  75          rho1(2)_spa_cf.sav    bin        saved values of system correlation functions
!!!!!  ./tools/read_output_para/ directory has the python programs to read the following output files !!!!!
!  76          ams_fermion.data      yes        saved the annihilation fermionic operators for the system defined in input file
!  77          ams_boson.data        yes        saved the annihilation operators for bosonic system with same DOFs of spin and orbital
!  78          ham_sys.data          yes        saved the system Hamiltonlian
!  79          res_corr.data         yes        saved the results of spectral decomposition schemes for time-domain bath correlation functions.
!  80          RDO_and_ADO.data      yes        saved the nonzero elements of reduced system density operator and 1st-ADOs for stationary calculation
!                                               if you need all ADOs, change 'do itier=1,2' to 'do itier=1,ntier0'
!!!!!  ./tools/read_output_para/ directory has the python programs to read the above output files     !!!!!
!------------------------------------------------------------------------------------
!
character*10            :: date, time
integer                 :: itier, istat, init, nstep, memo, ispin, ialf, iorbs, ibath
integer                 :: ni, nj, nk, nbtmp, ndtmp, nballs, nnz
integer*8               :: length0(maxtier), length1(maxtier), length2(maxtier), length3(maxtier)
integer*8               :: len0, len1, len2, ntmp1, ntmp2, ntmp3, ntmp4, ifrac, ifrac2
integer*8               :: lni, lnj
integer                 :: iorbs_dos, ispin_dos, maxit_dos
real*8                  :: dt_dos, tmax_dos, freq_dos, crit_dos
logical                 :: iexist, ldos, ljw_dos, lfreq_dos, lread_dos
logical                 :: fexist
real*8                  :: cpu1, cpu2, cpu3, cpu4, cpu5, cpu6
real*8                  :: dtmp1, dtmp2, dtmp3, dtmp4, dunitt, dunitj
real*8                  :: dspin2(3)
real*8                  :: tt, jleft, jright, occ, eshift(maxalf, maxspin)
real*8                  :: jleftu, jleftd, jrightu, jrightd, occu, occd, sz, sz2
real*8                  :: esys, esb, etot
real*8,     allocatable :: tmpmat(:,:)   
external                :: ifrac, ifrac2
complex*16, allocatable :: sigmau(:), sigmad(:)
complex*16              :: ctmp1
complex*16, allocatable :: cmat1(:,:)
!
integer   :: methodss                ! steady state method option : 
                                     ! methodss = 0, use Bicg  method (BiConjugate Gradient)
                                     ! methodss = 1, use DIIS  method (Direct Inversion in Iterative Subspace), now obsolete
                                     ! methodss = 2, use TFQMR method (Transpose-Free Quasi-Minimal-Residue), now default
                                     ! methodss = 3, use CPL-QMR method (Coupled two-term Look-ahead QMR), still has problem
                                     ! methodss = 4, use JACOBI method (Jacobi iteration method)
namelist / method  / methodss
namelist / units   / funits, runits  ! funits : unit option for physical quantities : 
                                     !      0 : energy-eV,  time-fs, freq-10^15 Hz, current-nA (         for molecular devices) 
                                     !          heat_current-nWatt
                                     !      1 : energy-meV, time-ps, freq-10^12 Hz, current-pA (default, for quantum dots     ) 
                                     !          heat_current-fWatt
                                     ! runits : relative units for current and time (when plotting j(t)) 
                                     !      0 for no, 1 for yes
namelist / resume   / icont, lresume, nresume
namelist / converge / maxit0, crit, ltrend
namelist / guess0   / grandom, lhseig, lhsdeg, lbzman, iguess, nhseig, toldeg
                                     ! grandom=true, use random numbers for the initial guess of rho (default depends on methodss)
                                     ! iguess = 0 (default), diagonal elements are all the same, and off-diagonal elements are all zero
                                     ! for norbs=2 & nspin=2, different values of iguess correspond to different spin states
namelist / inifield / lfld
namelist / debug    / lchkado, itierchk, inumchk, lchkmm, idiffchk, lchksame, lchkbig, dchkbig, nchkcount, nchklong, lchksparse, lchkdos, lchkbigado, loutput
namelist / tdjob    / tdmethod
namelist / jobinfo  / mfdjob, pfdjob, psdjob, psdlor, psdfff, psfjob, lscale, lsymlor, lsdraw, lsparse, &
                      ltrun_der, lsimple, lwalf, lscreen, lanahs, lspdhs, itype_fff, itype_psf
namelist / filter   / lfilter, nfilter, nfilter_count, nfilter_long
namelist / dos      / ldos, iorbs_dos, ispin_dos, dt_dos, tmax_dos, ljw_dos, lfreq_dos, lread_dos, freq_dos, maxit_dos, crit_dos
namelist / td4ground / ltd4gr, jconvg, tref
namelist / localtemp / linfcom, nstep_loc
namelist / deltav   / ldeltav
namelist / adiabatic / lad, lad_fast, lset_fast, lscba_ad, lcop_ad, ltd2st, ntier_ad, &
                       idegen_fast, ndegen_fast, dgama_slow, lcheck_ad, dratio_fast, ncor_fast
namelist / perturb   / lssp, nstep_ssp                   
namelist / project   / lproj, lprocc, nproji, nprojf, lprams, nprocc
namelist / omitado   / lomit, dfact_omit
!
call cpu_time(cpu1)
!
call genoutput
!
call sj_timer_start(1)
!
! set some global variables
!
tnow   = 0.d0
!lcf    = .false.
mfdjob    = .false.
pfdjob    = .false.
psdjob    = .true.
psdlor    = .false.
lscale    = .false.
lsymlor   = .false.
lsdraw    = .true.
lsparse   = .false.
ltrun_der = .false.
lfermi_exact = .true.
lsimple   = .false.
lwalf     = .true.
psdfff    = .false.
psfjob    = .false.
lscreen   = .false.
lanahs    = .true.  
lspdhs    = .false.
itype_fff = 1
itype_psf = 1
rewind(5)
read(5, jobinfo, end=111)
111 continue
write(6,*)
!
call get_omp_env
!
! truncation scheme
!
write(6,*)
if (ltrun_der) then
   write(6,*)'main: Truncation of hierarchy is for derivatives of ADOS '
else
   write(6,*)'main: COP truncation is adopted for the hierarchy '
end if
!
! make sure only one of MFD, PFD, PSD, PSD_LOR, PSD_FFF, PSF schemes is invoked
!
if (psfjob) then
    psdfff = .false.
    psdlor = .false.
    psdjob = .false. 
    pfdjob = .false.
    mfdjob = .false.
else if (psdfff) then
    psdlor = .false.
    psdjob = .false. 
    pfdjob = .false.
    mfdjob = .false.
else if (psdlor) then
    psdjob = .false. 
    pfdjob = .false.
    mfdjob = .false.
else if (psdjob) then
    pfdjob = .false.
    mfdjob = .false. 
else if (pfdjob) then
    mfdjob = .false. 
end if
!
istat = 0
if (mfdjob) istat = istat + 1
if (pfdjob) istat = istat + 1
if (psdjob) istat = istat + 1
if (psdlor) istat = istat + 1
if (psdfff) istat = istat + 1
if (psfjob) istat = istat + 1
if (istat .gt. 1) then
    write(6,*)'heom: error! spectral decomposition scheme conflict!'
    write(6,*)'heom: mfdjob = ', mfdjob
    write(6,*)'heom: pfdjob = ', pfdjob
    write(6,*)'heom: psdjob = ', psdjob
    write(6,*)'heom: psdlor = ', psdlor
    write(6,*)'heom: psdfff = ', psdfff
    write(6,*)'heom: psfjob = ', psfjob
    stop
end if
write(6,*)
if (mfdjob) then
    write(6,*)' Multi-Frequency-Dispersed (MFD) Scheme is used '
else if (pfdjob) then 
    write(6,*)' Polynomial-Function-Decomposition (PFD) Scheme is used '
else if (psdjob) then 
    write(6,*)' Pade-Spectral-Decomposition (PSD) Scheme is used '
else if (psdlor) then
    write(6,*)' Pade with Lorentzian-Residue-Fitting Scheme is used '
else if (psdfff) then
    write(6,*)' Pade with Fano-Function-Fitting Scheme is used '
else if (psfjob) then
    write(6,*)' Prony-Spectral-Fitting (PSF) Scheme is used '
    write(6,*)' itype_psf = ', itype_psf
else
    write(6,*)' Matsubara-Spectral-Decomposition (MSD) Scheme is used '
end if
if (lscale) then
    write(6,*)' Scale all ADOs to be dimensionless density operators '
end if
call flush(6)
if (pfdjob .or. psdjob .or. psdlor .or. psdfff .or. psfjob) then
    lfermi_exact = .false.
end if
!
if (lsimple) then
   write(6,*)
   write(6,*)'main: lsimple = .true. '
   write(6,*)
   if (lwalf) then
      write(6,*)'main: lwalf = T: operator includes ialf index '
   else
      write(6,*)'main: lwalf = F: operator excludes ialf index '
   end if
   if (ltrun_der) then
      write(6,*)
      write(6,*)'main: error! Derivative mode unavailable for lsimple=T '
      stop
   end if
   call flush(6)
end if
if (psdfff) then
    if (ltrun_der) then
        write(6,*)
        write(6,*)'main: error! psdfff=T and ltrun_der=T found '
        write(6,*)'main: code unavailable yet '
        stop
    end if
end if
!
if (lscreen) then
    write(6,*)
    write(6,*)'main: pre-screen all decoupled dissipative bath modes '
    call flush(6)
end if
!
memory        = 0.d0
igroundsteady = 0
itdskip       = 1
!
funits = 1
runits = 0
rewind(5)
read(5, units, end=102)
102 continue
write(6,*)
if (funits .eq. 0) then
    write(6,*)'units for molecular devices are employed                '
    write(6,*)'energy~ eV, time~fs, freq.~10^15Hz, j~nA, jheat~nW, capa.~aF '
    write(6,*)
    write(6,*)'if you want the input and output values in the atomic unit'
    write(6,*)'please add the following statement in the input file'
    write(6,*)'"','$units funits = 2 $end','"'
else if (funits .eq. 1) then
    write(6,*)'units for quantum dots are employed                      '
    write(6,*)'energy~meV, time~ps, freq.~10^12Hz, j~pA, jheat~fW, capa.~fF '
    write(6,*)'if you want the input and output values in the atomic unit'
    write(6,*)'please add the following statement in the input file'
    write(6,*)'"','$units funits = 2 $end','"'
else if (funits .eq. 2) then
    write(6,*)'input and output in the atomic unit                    '
    write(6,*)'energy~1 Hartree, time~1 a.u. sec'
    write(6,*)'1 Hartree  = 2.721138625E+001 eV'
    write(6,*)'1 a.u. sec = 2.418884327E-002 fs'
    write(6,*)'hbar =',hbar_AU,' Hartree * a.u. sec'
    write(6,*)
    write(6,*)'if you want to convert the SI unit in input file into the atomic unit'
    write(6,*)'divide energy parameters by a factor of 27.21138625'
    write(6,*)'multiply time parameters by a factor of 27.21138625'
    write(6,*)
    write(6,*)'if you want to convert the atomic unit in output files into the SI unit'
    write(6,*)'multiply energy and current by a factor of 27.21138625'
    write(6,*)'divide time parameters by a factor of 27.21138625'
else
    write(6,*)'main: error! unknown funits= ', funits
    stop
end if
!
if (runits .eq. 1) then
  write(6,*)
  write(6,*)' relative units for output j(t) is invoked '
end if
call flush(6)
!
lresume = .true.
icont   = 0
nresume = 500
rewind(5)
read(5, resume, end=104)
104 continue
!
!onelead = .false.
!rewind(5)
!read(5, lead, end=105) 
!105 continue
!
lad         = .false. 
lcheck_ad   = .false.
lad_fast    = .true. 
lset_fast   = .true. 
lscba_ad    = .true. 
lcop_ad     = .true. 
ltd2st      = .false.
ntier_ad    = 1          ! it is essential that when lad=F, ntier_ad < ntier (note added on Aug 2, 2020: no need)
idegen_fast = 1      
!ndegen_fast = maxtier
ndegen_fast = 2
dgama_slow  = 0.2d0    
dratio_fast = 0.9d0
ncor_fast   = 0
rewind(5)
read(5, adiabatic, end=120)
120 continue
dgama_slow = dgama_slow / hbar
!
lssp       = .false.
nstep_ssp  = 0
rewind(5)
read(5, perturb, end=121)
121 continue
if (lssp) then
    write(6,*)
    write(6,*)'main: lssp=T found '
    write(6,*)'main: steady-state perturbation method is used '
    write(6,*)'main: nstep_ssp = ', nstep_ssp
    if (nstep_ssp .ne. 1 .and. nstep_ssp .ne. 2) then
        write(6,*)
        write(6,*)'main: error! '
        write(6,*)'main: unknown nstep_ssp = ', nstep_ssp
        stop
    end if
    call flush(6)
end if
!
lomit = .false.
dfact_omit = 4.d-6  ! to modify --- adapt to different spectral decomposition schemes
rewind(5)
read(5, omitado, end=122) 
122 continue
if (lomit) then
    write(6,*)
    write(6,*)'main: lomit=T found '
    write(6,*)'main: omit ADOs which are presumed to be small '
    write(6,*)'main: dfact_omit = ', dfact_omit
    call flush(6)
    if (lad) then
        write(6,*)'main: error! lad=T found '
        write(6,*)'main: code not ready yet for lomit=T and lad=T '
        stop
    end if
    if (lscale) then
        write(6,*)'main: error! lscale=T found '
        write(6,*)'main: code not ready yet for lomit=T and lscale=T '
        stop
    end if
end if
! 
! read basic information from input file
! and calculate reservoir correlation functions
!
call evaluatepara(init)
!
if (lssp) then
    if (init .ne. 7) then
        init = 7
        write(6,*)'main: set init = 7 '
    end if
    call flush(6)
end if
!
call prelude_hb
!                      
write(6,*)
write(6,*)' highest tier to achieve         : ', ntier - 1
write(6,*)' number of leads                 : ', nalf
write(6,*)' number of signs                 : ', nsgn
write(6,*)' number of Drude terms           : ', ndrude
write(6,*)' number of Mastubara terms       : ', nmats
if (mfdjob) then
  write(6,*)' number of frequency terms/grids : ', nfreq
end if
write(6,*)' number of exponential functions : ', ncor
write(6,*)' number of energy levels         : ', norbs
write(6,*)' number of spin considered       : ', nspin
!write(6,*)' number of spin for leads        : ', lspin
write(6,*)
write(6,*)' number of balls available       : ', ntier - 1
write(6,*)' number of operators             : ', nopr
write(6,*)' number of drawers               : ', nvar
write(6,*)' dimension of Hamiltonian        : ', nrho
call flush(6)
!
if (ntier .le. 1) then
 write(6,*)
 write(6,*)' error! no auxiliary term is allowed ', ntier
 write(6,*)' namely, no balls available.         '
 stop
end if
!
ndrawer_slow = 0 
ncor_slow = 0
if (lad) then
    call prelude_ad
    if (lcheck_ad) then
        write(6,*)
        write(6,*)'main: lcheck_ad=T found, stop the code after '
        write(6,*)'      checking <prelude_ad.f90>              '
        call flush(6)
        goto 302
    end if
end if
! 
! calculate the size of hierarchy
!
ntier0 = max(ntier, ntier_ad)
length0(1:ntier0) = 0
length1(1:ntier0) = 0
length2(1:ntier0) = 0
length3(1:ntier0) = 0
!
ntmp1 = 0
ntmp2 = 0
ntmp3 = 0
ntmp4 = 0
write(6,*)
write(6,*)' itier, # of cases, index memory '
call flush(6)
do itier=1,ntier0
   nballs = itier - 1
   call calclength(itier, len0, len2)   ! len0 : # of drawer-sets; len2 : # of operator-sets
   if (itier .eq. 1) then
       length0(itier) = len0
       length1(itier) = length0(itier)
       length2(itier) = len2
       length3(itier) = length2(itier)
   else
       if ( mod(nballs, 2) .ne. 0 ) then
           length0(itier) = len0 / 2
           length2(itier) = len2 / 2
       else
           ndtmp = nopr / 2            ! # of pairs
           nbtmp = nballs / 2              
           len1  = ifrac2(ndtmp + nbtmp - 1, max(nbtmp, ndtmp - 1)) / ifrac(min(nbtmp, ndtmp - 1))
           length2(itier) = (len1 + len2) / 2
           if (itier .gt. ntier) then
               if (lad_fast) then
                   length0(itier) = length2(itier) * ncor_fast**nballs
               else
                   length0(itier) = length2(itier) * ncor_slow**nballs
               end if
           else
               length0(itier) = length2(itier) * ncor**nballs
           end if
       end if
       length1(itier) = nballs * length0(itier)
       length3(itier) = nballs * length2(itier)
   end if
   ntmp1 = ntmp1 + length0(itier)  ! length0 and length1 here are upper bounds 
   ntmp2 = ntmp2 + length1(itier)  ! (symmetry not fully considered yet)
   ntmp3 = ntmp3 + length2(itier)
   ntmp4 = ntmp4 + length3(itier)
!  
   write(6,728)itier, length2(itier), length3(itier), length0(itier), length1(itier)
   call flush(6)
enddo
728 format(I4, 2x, 4(I12, 2x))
!
! \bar{k} is the complement of multi-component index k (with opposite sign)
!
nunk      = ntmp1   ! symmetry between {k} and {\bar{k}} has been taken into account.
ntable    = ntmp2   ! indextable is compressed by a ratio of 1/2 (approximately) because of symmetry.
noprunk   = ntmp3
noprtable = ntmp4
write(6,*)
write(6,*)'  total number of matrices is ', nunk 
write(6,*)'  total number of unknowns is ', nunk * nrho**2
write(6,*)'    est. length of indextable ', ntable
write(6,*)'    number of combinations is ', noprunk
write(6,*)'     length of operator-table ', noprtable
write(6,*)'      highest tier considered ', ntier
call flush(6)
!
allocate(indextable(ntable),  STAT=istat)
allocate(oprtable(noprtable), STAT=istat)
allocate(ifirst(ntier0), ilast(ntier0), STAT=istat)
allocate(nfirst(ntier0), nlast(ntier0), STAT=istat)
allocate(ioprfirst(ntier0), ioprlast(ntier0), STAT=istat)
allocate(noprfirst(ntier0), noprlast(ntier0), STAT=istat)
allocate(iatmp1(ntier0), iatmp2(ntier0), STAT=istat)
!
write(6,*) 
if (istat .eq. 0) then
  write(6,*)' OK! memory allocation for indextable successful.'
  write(6,*)' memory for indextable : ', dble(ntable * 4) / dble(1024**2), ' MB'
else
  write(6,*)' error! memory allocation for indextable failed. '
  write(6,*)' memory for indextable : ', dble(ntable * 4) / dble(1024**2), ' MB'
  stop
endif
call flush(6)
!
memory = memory + dble(ntable * 4) / dble(1024**2)
!
lfilter = .false.
nfilter = 10          ! currently no use
nfilter_count = 2     !  = 0: no filter, turn off filtering, if > (ntier-1), the whole ntier is screened, error
nfilter_long  = 3
rewind(5)
read(5, filter, end=112) 
112 continue
if (lfilter .and. lad) then
    write(6,*)
    write(6,*)'main: error! lfilter=T and lad=T found ', lfilter, lad
    write(6,*)'main: code unavailable yet '
    stop
end if
if (nfilter_count .lt. 0) then
   write(6,*)'main: error! nfilter_count < 0 found ', nfilter_count
   stop
end if
if (nfilter_long .lt. 0) then
   write(6,*)'main: error! nfilter_long < 0 found ', nfilter_long
   stop
end if
if (nfilter_count .eq. 0) then  ! turn off filter 
   lfilter = .false.
   write(6,*)'main: nfilter_count == 0 found, turn off filter '
   call flush(6)
end if
!
if (.not. lsdraw) then
    write(6,*)
    write(6,*)'main: lsdraw = .false. found '
    write(6,*)'main: discard all ADOs involving identical drawers '
    write(6,*)'main: code unavailable yet   '
    stop
end if
!
do itier=1,ntier0
   if (lfilter .and. itier .eq. ntier) then
      call ball2drawer_filter(itier, length0(1), length2(1))
      cycle
   end if
   call ball2drawer(itier, length0(1), length2(1))
end do
call sortindextable3
call examindextable(length0(1))
!
ntable = 1
do itier=2,ntier0
  ntable = ntable + length0(itier) * (itier - 1)
end do
call clearindextable(length0(1))
!
if (lfilter) then
   call buildfilter
end if
!
write(6,*)
write(6,*)' itier,  ifirst,  ilast '
do itier=1,ntier0
 write(6,*)itier, ifirst(itier), ilast(itier)
end do
write(6,*)
write(6,*)' itier,  nfirst,  nlast '
do itier=1,ntier0
 write(6,*)itier, nfirst(itier), nlast(itier)
enddo
write(6,*)
write(6,*)' itier,  ioprfirst,  ioprlast '
do itier=1,ntier0
  write(6,*)itier, ioprfirst(itier), ioprlast(itier)
end do
write(6,*)
write(6,*)' itier,  noprfirst,  noprlast '
do itier=1,ntier0
  write(6,*)itier, noprfirst(itier), noprlast(itier)
end do
call flush(6)
!
ntable = ilast(ntier0)
nunk   = nlast(ntier0)  ! symmetry between {k} and {\bar{k}} has been taken into account.
write(6,*)
write(6,*)'  actual length of indextable ', ntable
write(6,*)'  total number of unknowns is ', nunk * nrho**2
!
allocate(tmpmat(nrho,nrho), STAT=istat)
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
allocate(cmat1(nrho, nrho), STAT=istat)
!
! read options
! 
methodss = 2
rewind(5)
read(5, method, end=201)
201 continue
if (lfilter .and. methodss .ne. 2) then
   write(6,*)'main: error! for lfilter=T, only methodss=2 available now'
   write(6,*)'main: methodss = ', methodss
   stop
end if
!
grandom = .true.
iguess  = 0
lhseig  = .false.
lhsdeg  = .true.
!lbzman  = .false.
lbzman  = .true.  ! use Boltzmann as default initial guess (2022-05-06)
nhseig  = 1
toldeg  = 1.d-3
lreadi  = .false.
if (methodss .eq. 2 .or. methodss .eq. 4) then  ! tfqmr uses zero as initial guess
   grandom = .false.
end if
rewind(5)
read(5, guess0, end=203)
203 continue
if (lhseig) then
    if (lbzman) then
        write(6,*)'main: lhseig=T found, override lbzman '
        write(6,*)'main: build initial rhosys from a system eigenvector '
        call flush(6)
        lbzman = .false.
    end if
end if
!
call allocatememory
call buildoperator
call checkanticomm
call arrangeoperator
!
if (lhb) then
!   call allocate_hb(istat)
   call calc_qa_hb
end if
!
! 3d spin operator
!
! namelist spin3d already read in <evaluatepara.f90>
if (lspin3d) then
   call buildspin3d
end if
!
! Find connections with upper- and lower-tier ADOs.
!
call buildcoefindex
!
select case (init)
case (0)
   igroundsteady = 0
case (1)
   igroundsteady = 2
case (2)
   igroundsteady = 1
case (3)
   igroundsteady = 1
case (4)
   igroundsteady = 2
case (5)
   igroundsteady = 0
case (6)
   igroundsteady = 2
case (7)
   igroundsteady = 1
case default
   write(6,*)'main: error! unknown init = ', init
   stop
end select
!
call allocate_system
!
! construct system projector
!
lproj  = .false.
lprocc = .false. 
nprocc = 1
nproji = 1
nprojf = nrho
lprgen = .false.
ldopro = .false.
lprams = .true.  ! note added on 2022-04-26: by default use projected a_ms, if lproj=T
rewind(5)
read(5,project,end=115)
115 continue
if (lproj) then
    if (.not. lsparse) then
        write(6,*)'main: error! lsparse = ', lsparse
        write(6,*)'main: lproj=T now only works in SPARSE mode '
        stop
    end if
    if (lprocc) then 
! projector based on occupation number    
        if (nprocc .le. 0) then
            write(6,*)'main: error! '
            write(6,*)'main: unexpected nprocc = ', nprocc
            stop
        end if
    else 
! projector based on energy
        if (.not. (1 .le. nproji .and. nproji .le. nprojf .and. nprojf .le. nrho)) then
            write(6,*)'main: error! '
            write(6,*)'main: erroneous nproji, nprojf: ', nproji, nprojf
            stop
        end if
    end if
    call buildprojector
    lanahs = .true.  ! force lanahs=T to initialize rho_spa in <analyze_hs.f90>
    if (.not. lhseig) then
        lbzman = .true.  ! force lbzman=T to initialize rho_spa in <analyze_hs.f90>
    end if
else
    lprams = .false.
end if
! 
! SPARSE mode : sparsity of each ADO is predetermined
!
if (lsparse) then
   write(6,*)
   write(6,*)'main:   SPARSE MODE is turned on  '
   call flush(6)
   call prelude_spa
end if
!
! Initialize unknowns (specify job types)
!
! init = 0  :  igroundsteady = 0, itdskip = 1 (solve linear problem for ground state, 
!              and use the ground state info as initial to carry out td propagation)
!
! init = 1  :  igroundsteady = 2, itdskip = 1 (read input file to obtain initial values
!              of unknowns and then carry out td propagation, does not solve any linear problem)
! 
! init = 2  :  igroundsteady = 1, itdskip = 0 (solve linear problem for steady state, no
!              td job is carried out, initial guess for steady state rho is generated by program)
!
! init = 3  :  igroundsteady = 1, itdskip = 0 (solve linear problem for steady state, no
!              td job is carried out, initial guess for steady state rho is read from file if available)
!
! init = 4  :  igroundsteady = 2, itdskip = 1 (try to solve steady state or ground state via
!              td propagation approach)
!
! init = 5  :  igroundsteady = 0, itdskip = 1 (read input file to obtain initial values
!              of unknowns, solve linear problem for ground state, and use the ground state info
!              as initial to carry out td propagation)
! 
! init = 6  :  igroundsteady = 2, itdskip = 0 (read break-point information at a certain time
!              during time evolution, and calculate the 'adiabatic' transient system spectral 
!              function (obsolete??)
!              Note (2021-11-26): option reactivated with the development of adiabatic terminator
!
! init = 7  :  igroundsteady = 1, itdskip = 0 (solve linear problem for steady state 
!              by using perturbation method, no td job is carried out)
!
write(6,*)
if (init .eq. 0) then
   call initialize1
   write(6,*)' initialization done '
   write(6,*)' solve linear problem for ground state, and use it '
   write(6,*)' as initial state for subsequent td propagation.   '
!   igroundsteady = 0
   itdskip       = 1
!
elseif (init .eq. 1) then       
   call initialize2(length0(1))      
   write(6,*)' initialization done                               '
   if (lsparse) then
      write(6,*)' <rho_spa.sav> read successfully                '
   end if
   write(6,*)' do not solve linear problem, td propagation is    '
   write(6,*)' carried out straightforwardly.                    '
!   igroundsteady = 2 
   itdskip       = 1
   lreadi = .true.
!
elseif (init .eq. 2) then
   call initialize1
   write(6,*)' initialization done                               '
   write(6,*)' solve linear problem for steady state, no td job! '
!   igroundsteady = 1
   itdskip       = 0
!
   if (dfieldtype .ne. 0 .and. (.not. forcesteady)) then
     write(6,*)
     write(6,*)' error! dfieldtype!=0 used! ', dfieldtype
     write(6,*)' no steady state exists under this field '
     stop
   end if
!
elseif (init .eq. 3) then
   if (lsparse) then
      inquire(file='rho_spa.sav', exist=iexist, err=999)
   else
      inquire(file='variables.sav', exist=iexist, err=999)
   end if
   999 continue
   if (iexist) then
       call initialize2(length0(1))
       if (lsparse) then
          write(6,*)' <rho_spa.sav> read successfully              '
       end if
       lreadi = .true.
   else
       call initialize1
       lreadi = .false. 
   end if
   write(6,*)' initialization done                               '
!   igroundsteady = 1
   itdskip       = 0
!
   if (dfieldtype .ne. 0 .and. (.not. forcesteady)) then
      write(6,*)
      write(6,*)' error! dfieldtype!=0 used! ', dfieldtype
      write(6,*)' no steady state exists under this field '
      stop
   end if
!
else if (init .eq. 4) then
   call initialize1
   write(6,*)' initialization done                               '
   write(6,*)' try to solve steady state via td propagation.     '
!   igroundsteady = 2
   itdskip       = 1
!
else if (init .eq. 5) then
   call initialize2(length0(1))      
   write(6,*)' initialization done                               '
   if (lsparse) then
      write(6,*)' <rho_spa.sav> read successfully                '
   end if
   write(6,*)' solve linear problem for ground state, and use it '
   write(6,*)' as initial state for subsequent td propagation.   '
!   igroundsteady = 0
   itdskip       = 1
   lreadi = .true.
!
else if (init .eq. 6) then
   if (lhb) then
      write(6,*)
      write(6,*)'main: error! lhb=T not avaiable for init=6 '
      stop
   end if
   call initialize1
   inquire(file='TAPE.resume', exist=iexist, err=998)
   998 continue
   if (iexist) then
     call resumejob(tt, 1)
     tnow = tt
     write(6,*)'main: <TAPE.resume> read successfully              '
     write(6,*)'main: tnow = ', tnow
   else
     write(6,*)'main: error! cannot find <TAPE.resume>, abort ...  '
     stop
   end if
!   igroundsteady = 2
   itdskip       = 0
   lreadi = .true.
!
else if (init .eq. 7) then
   itdskip = 0
   call initialize_ssp(length0(1))
   lreadi = .true.
!
else 
   write(6,*)
   write(6,*)' error! unknown initialization flag ', init
   stop
end if
call flush(6)
!
! make ground state Hamiltonian
!
tt = 0.d0
call refresh_rhosys
lprths = .true.
ldopro = .false.  ! turn off projector for <analyze_hs>
call calchs(rhosys)
if (lanahs) then
    call analyze_hs(hs)  
!    call refresh_rhosys
end if

!write(6,*)
!write(6,*)'main: rhosys after analyze_hs '
!call cmatout(nrho, nrho, rhosys, dmtmp1)
!stop

if (lproj .and. lprgen) then ! construct projected system Hamiltonian
    ldopro = .true.
    call calchs(rhosys) 
end if
!
lchksame = .false.
lchkmm   = .false.
idiffchk = 5
lchkado  = .false.
lchkbig  = .false.
lchksparse = .false.
lchkdos  = .false.
lchkbigado = .false.
dchkbig  = 1.d-5
nchkcount = 2
nchklong  = 2
itierchk = 2
inumchk  = 0
loutput = .true.
rewind(5)
read(5, debug, end=106)
106 continue
!
if (lsparse) then
   lchkado = .false.
   write(6,*)
   write(6,*)'main: lchkado is turned off in sparse mode '
end if
!
!
! solve steady state (or ground state) 
!
if (igroundsteady .ge. 2) goto 90
if (init .eq. 0 .and. icont .eq. 1) goto 90
!
maxit0  = 10000
crit    = 1.d-8
!ltrend  = .true.  ! note added on 2022-01-01: only works for methodss==2 (TFQMR)
ltrend  = .false. ! note added on 2022-04-24: change the default back to FALSE
                  ! note (2022-05-16): ltrend becomes obsolete
rewind(5)
read(5, converge, end=202)
202 continue
write(6,*)
write(6,*)' steady state parameters  '
write(6,*)' maximal cycles allowed : ', maxit0 
! criterion is the square root of residue norm
write(6,*)' convergence criteria   : ', crit 
!if (lproj) then
!    ltrend = .true.
!    write(6,*)' lproj=T found, set ltrend=T '
!    call flush(6)
!end if
!if (ltrend) then
!    write(6,*)' put trace(rho) = 1 at the end of vector rho '
!else
!    write(6,*)' put trace(rho) = 1 at the beginning of rho '
!end if
!call flush(6)
!!
if (methodss .eq. 0) then
   if (lsparse) then
       write(6,*)
       write(6,*)'main: error! BICG solver in SPARSE mode not available '
       stop
   end if
   if (lad) then
       write(6,*)
       write(6,*)'main: error! BICG solver with lad=T is not available '
       stop
   end if
   call allocate_bicg(istat, memo)
   if (memo .eq. 0) then
      call solvebicg
   else
      write(6,*)
      write(6,*)' error! nozero memo found : ', memo
      write(6,*)' memory request too demanding, abort... '
      stop
   end if
   call free_bicg(memo)
!
!else if (methodss .eq. 1) then   ! DIIS is very inefficient, and is no longer used and supported
!  write(6,*)
!  write(6,*)' Warning! '
!  write(6,*)' DIIS is inefficient and possibly erroreous. Be cautious with results! '
!  call flush(6)
!  call allocate_diis(istat)
!  call solvediis
!  call free_diis
!
else if (methodss .eq. 2) then
! Note that TFQMR uses residue norm as criterion
  cpu3 = omp_get_wtime()
  call allocate_tfqmr(istat)
  call solve_tfqmr
  call free_tfqmr
  cpu4 = omp_get_wtime()
  write(6,*)
  write(6,*)'main: CPU time for TFQMR solver ', cpu4-cpu3
  call flush(6)
!
!else if (methodss .eq. 3) then
!  call allocate_cplqmr(istat)
!  call solve_cplqmr
!  call free_cplqmr
!
else if (methodss .eq. 4) then
! Note that JACOBI uses residue norm as criterion
   if (lad) then
       write(6,*)
       write(6,*)'main: error! JACOBI solver with lad=T is not available '
       stop
   end if
   if (lssp) then
       write(6,*)
       write(6,*)'main: error! JACOBI solver with lssp=T is not available '
       stop
   end if
   call allocate_jacobi(istat)
   call solve_jacobi
   call free_jacobi
!
else 
  write(6,*)
  write(6,*)'heom: error! unknown steady-state option, methodss = ', methodss
  stop
end if
!
! check steady state density matrix
!
write(6,*)
if (igroundsteady .eq. 0) then
  write(6,*)'  equilibrium reduce density matrix rho '
else if (igroundsteady .eq. 1) then
  write(6,*)' steady-state reduce density matrix rho '
end if
call refresh_rhosys
call cmatout(nrho, nrho, rhosys, dmtmp1)
call checkrho
!
if ( (.not. lsparse) .and. lchksparse ) then
   call checkado
end if
!
! output ground state density matrix rho 
!
if (lsparse) then
   if (igroundsteady .eq. 0) then
      open(unit=51, file='rho_spa.gr', form='binary', status='unknown')
   else if (igroundsteady .eq. 1) then
      open(unit=51, file='rho_spa.st', form='binary', status='unknown')
   end if
   if (igroundsteady .eq. 0 .or. igroundsteady .eq. 1) then
      rewind(51)
      write(51)norbs, nspin, ncor, ntier0, nalf
      write(51)nunk
      do lni=1,nunk
         nnz = nnz_spa(lni)
         lnj = ind_spa(lni)
         if (nnz .gt. 0) then
            write(51)(rho_spa(lnj-1+ni), ni=1,nnz)
         end if
      end do
      close(51)
   end if
else
   if (igroundsteady .eq. 0) then
      open(unit=25, file='variables.gr', form='unformatted', status='unknown')
   else if (igroundsteady .eq. 1) then
      open(unit=25, file='variables.st', form='unformatted', status='unknown')
   end if
   if (igroundsteady .eq. 0 .or. igroundsteady .eq. 1) then
      rewind(25)
      write(25)norbs, nspin, ncor, ntier, nalf
      write(25)nunk
      do lni=1,nunk
         write(25)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
      end do
      close(25)
   end if
end if
!
90 continue
!
if (lhf) then
   call calc_rsdm(rho, nrho, norbs, nspin, rsdm)
   write(6,*)
   write(6,*)'heom: rsdm before time evolution '
   do ispin=1,nspin
      write(6,*)'heom: ispin = ', ispin
      call cmatout(norbs, norbs, rsdm(1,1,ispin), dmtmp1_sp)
   end do
end if
!
if (loutput) then
   call outputhamsys
   call outputcorr
   call outputams_fermion
   if (itdskip .eq. 0) then
      call outputnnz
   end if
   write(6,*)
   write(6,*)'main: output to files finished'
   write(6,*)
!   stop
end if
!
if (itdskip .eq. 0) goto 200   ! skip td if required
!if (lsparse) then
!   write(6,*)
!   write(6,*)'main: time evolution in SPARSE mode not available yet ' 
!   call flush(6)
!   goto 200
!end if
!
! start time propagation 
!
call refresh_rhosys
if (lsparse) then
   call prelude_td_spa
end if
!
call allocate_rk4(istat)
nstep = 0
!
ldeltav = .false.
rewind(5)
read(5, deltav, end=310)
310 continue
if (lequileads .and. ldeltav) then
   write(6,*)
   write(6,*)'main: error! deltav is not ready for lequileads = TRUE ' 
   stop
end if
if (ldeltav) then
   allocate(deltav_amps(nalf, nspin), STAT=istat)
   read(5,*) ((deltav_amps(ni,nj), nj=1,nspin), ni=1,nalf)
   deltav_amps = deltav_amps / hbar
   write(6,*)
   write(6,*)'main: delta voltage is applied as a perturbation '
   write(6,*)'main: amplitude, amps(nspin, malf) '
   write(6,*) ((deltav_amps(ni,nj)*hbar, nj=1,nspin), ni=1,nalf)
   call flush(6)
end if
!
write(6,*)
open(unit=13, file='curr.data', status='unknown')
rewind(13)
open(unit=16, file='popu.data', status='unknown')
rewind(16)
open(unit=72, file='engy.data', status='unknown')
rewind(72)
open(unit=28, file='rhodiag.data', status='unknown')
rewind(28)
open(unit=29, file='rhotime.data', form='formatted', status='unknown')
rewind(29)
open(unit=33, file='poccdetail.data', status='unknown')
rewind(33)
if (lbfield .and. nspin .gt. 1) then
   open(unit=43, file='chi_loc.data', status='unknown')
   rewind(43)
end if
if (lspin3d) then
   open(unit=61, file='sdot1.data', status='unknown')
   rewind(61)
   if (norbs .eq. 2 .and. nspin .eq. 2) then
      open(unit=62, file='sdot2.data', status='unknown')
      rewind(62)
      open(unit=63, file='spin_ddot.data', status='unknown')
      rewind(63)
   end if
end if
! 
if (lchkado) then
  write(6,*)
  write(6,*)' A check for ADO is required for tier: ', itierchk
  if (itierchk .lt. 1 .or. itierchk .gt. ntier) then
    write(6,*)' error itierchk ', itierchk, ntier
    stop
  end if
  write(6,*)' #ADO : ', nfirst(itierchk) + inumchk
  if (nfirst(itierchk) + inumchk .gt. nlast(itierchk)) then
    write(6,*)' error inumchk ', nfirst(itierchk)+inumchk, nlast(itierchk)
    stop
  end if
  call flush(6)
  open(unit=30, file='ado.data', status='unknown')
  rewind(30)
end if
!
if (runits .eq. 1) then                                   ! gamma / hbar      : gamma is the maximal total linewidth 
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
else
  dunitt = 1.d0
  dunitj = 1.d0
end if
!
if (lresume .and. icont .eq. 1) then                      ! resume from previous broken point 
  write(6,*)
  write(6,*)' resume a previous time-dependent job '
  call resumejob(tt, 1)
  tresume = tt
else                                                      ! start from t_0 (very beginning)
  lfld = .false.
  rewind(5)
  read(5, inifield, end=107)
  107 continue
  if (lfld) then
    write(6,*)
    write(6,*)' a delta-function-type field is imposed on the QD at t=t0 '
    call flush(6)
    call modinisystem
  end if
!
  if (fieldtype .eq. -1) then
     write(6,*)
     write(6,*)' a delta-function-type voltage pulse is applied to the lead '
     call flush(6)
     if (lsparse) then
        call modinileads_spa
     else
        call modinileads
     end if
  end if
!
  if (ldeltav) then
     if (fieldtype .eq. -1) then
        write(6,*)
        write(6,*)'main: error! fieldtype cannnot be -1 if ldeltav = TRUE '
        stop
     end if
     if (lsparse) then
        call modinileads_deltav_spa
     else
        call modinileads_deltav
     end if
     deallocate(deltav_amps, STAT=istat)
  end if
!
!  open(unit=31, file='rho0plus.data', form='unformatted', status='unknown')
!  rewind(31)
!  write(31)nrho
!  write(31)nunk
!  do lni=1,nunk
!    write(31)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
!  end do
!  close(31)
!
end if
!
write(6,*)
write(6,*)' initial rho ', tt
call cmatout(nrho, nrho, rhosys, dmtmp1)
!
write(6,*)
write(6,*)' TIME DEPENDENT CURRENT '
write(6,*)
call flush(6)
!
call cpu_time(cpu5)
!
tdmethod = 0
rewind(5)
read(5, tdjob, end=108)
108 continue
!
ltd4gr = .false.
jconvg = 1.d-2
tref   = 5.d0
rewind(5)
read(5, td4ground, end=412)
412 continue
if (ltd4gr) then
   write(6,*)'heom: use time propagation to solve ground state '
   write(6,*)'heom: jconvg = ', jconvg
   write(6,*)'heom: tref   = ', tref
   call flush(6)
end if
!
linfcom = .false.
nstep_loc = 10
rewind(5)
read(5, localtemp, end=114)
114 continue
if (linfcom) then
   write(6,*)
   write(6,*)'main: information compressibility is to be calculated '
   open(unit=47, file='local_temp.data', status='unknown')
   rewind(47)
end if
call flush(6)
!
write(6,*)
if (tdmethod .eq. 0) then
  write(6,*)' 4th-order Runge-Kutta algorithm is employed for time evolution '
else if (tdmethod .eq. 1) then
  write(6,*)' Chebyshev expansion algorithm is employed for time evolution '
  write(6,*)' (for step and delta function voltages) '
  if (lsparse) then
     write(6,*)
     write(6,*)'main: Chebyshev propagator in sparse mode not available yet '
     write(6,*)'main: Abort ... '
     stop
  end if
  if ( .not. (fieldtype .eq. 0 .or. fieldtype .eq. -1) ) then
     write(6,*)' error! applied voltage not compatible with Chebyshev method ', fieldtype
     stop
  end if
  if (fieldtype .eq. 0) then
     do ni=1,nalf
        do nj=1,nspin
           if (tchar(ni,nj) .gt. dpico) then
              write(6,*)'heom: tchar ~ 0 is required for fieldtype=0 and Chebyshev method'
              write(6,*)'heom: ialf = ', ni, ' ispin = ', nj, ' tchar = ', tchar(ni,nj)
              stop
           end if
        end do
     end do
  end if
  call prechebyshev
else 
  write(6,*)' error! unknown tdmethod', tdmethod
  stop
end if
write(6,*)
call flush(6)
!
call allocate_tmpmem_omp
!
if (tdmethod .eq. 0) then
  100 continue
!
  if (lresume) then
    if (nstep .ge. 0 .and. mod(nstep, nresume) .eq. 0) then
      call resumejob(tt, 2)
      write(29,518)tt, ((dble (rhosys(ni,nj)), ni=1,nj), nj=1,nrho),   &
                       ((dimag(rhosys(ni,nj)), ni=1,nj), nj=1,nrho)        
    end if
  end if
!
  inquire(file="stoptd", exist=fexist)
  if (fexist) then
    call resumejob(tt, 2)
    write(6,*)
    write(6,*)'main: file <stoptd> detected, calculation terminates '
    call flush(6)
    goto 911
  end if
!
  call calcocc(occu, occd, sigmau, sigmad)
  occ = occu + occd
!call calcj(jleftu, jrightu, jleftd, jrightd)
!write (6,518)tt, jleftu, jleftd, jleftu+jleftd, occu, occd, 0.d0
!call flush(6)
!stop
!
  if (nspin .gt. 1) then
     cmtmp1 = czero
     cmtmp2 = czero
     do iorbs=1,norbs
        cmtmp3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 1), 0.d0)
        call zamsmm1('l', 'c', 'n', iorbs, 1, cunity, cmtmp3, nrho, cunity, cmtmp1, nrho)
        cmtmp4(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 2), 0.d0)
        call zamsmm1('l', 'c', 'n', iorbs, 2, cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
     end do
     cmtmp1 = cmtmp1 - cmtmp2          ! n_up - n_down
!
     call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.5d0, 0.d0), cmtmp1, nrho, rhosys, nrho,  &
                czero, cmtmp2, nrho)   
     sz = 0.d0
     do ni=1,nrho
        sz = sz + dble(cmtmp2(ni,ni))
     end do
!
     call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.25d0, 0.d0), cmtmp1, nrho, cmtmp1, nrho,  &
                czero, cmtmp2, nrho)   ! (n_up - n_down)^2 / 4
     call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, rhosys, nrho,                &
                czero, cmtmp1, nrho)
     sz2 = 0.d0
     do ni=1,nrho
        sz2 = sz2 + dble(cmtmp1(ni,ni))
     end do
  end if
!
  lrun_loc = .false.
  if (linfcom .and. mod(nstep, nstep_loc) .eq. 0) then
     lrun_loc = .true.
  end if
!
  call cpu_time(cpu3)
  call intstp(nstep, tt, esys, esb, etot, jleftu, jleftd, jrightu, jrightd)
  call refresh_rhosys
  call cpu_time(cpu4)
!
  jleft  = jleftu  + jleftd
  if (nalf .gt. 1) jright = jrightu + jrightd 
!
  do ialf=1,nalf
     do ispin=1,nspin
        call getengyshift(tt, ialf, ispin, eshift(ialf,ispin))
     end do
  end do
  eshift = eshift * hbar
!
  if (lrun_loc) then
     write(47,518)tt, 1.d0/dinfcom_loc*hbar, dinfcom_loc/hbar, energy_loc*hbar, entropy_loc
     call flush(47)
  end if
!
  if (nspin .eq. 1) then
    if (nalf .eq. 1) then
       write (6,518)tt, jleftu, occu, cpu4-cpu3
       write(13,518)tt*dunitt, eshift(1,1), jleft*dunitj, occu
       write(16,518)tt, eshift(1,1), occu
       write(28,518)tt, (dble(rhosys(ni,ni)), ni=1,nrho)
       if (lchkado) then
         write(30,518)tt, ((dble(rho(ni,nj,nfirst(itierchk)+inumchk)), ni=1,nrho), nj=1,nrho), &
                         ((dimag(rho(ni,nj,nfirst(itierchk)+inumchk)), ni=1,nrho), nj=1,nrho)
       end if
    else if (nalf .eq. 2) then
       write (6,518)tt, jleftu, jrightu, jleftu+jrightu, occu, cpu4-cpu3
       write(13,518)tt*dunitt, eshift(1,1), eshift(2,1), jleft*dunitj, jright*dunitj, &
                    (jleft+jright)*dunitj, occu
       write(16,518)tt, eshift(1,1), eshift(2,1), occu
    else  ! multi-leads 
       write (6,518)tt, (jleads(ialf), ialf=1,nalf), occu, cpu4-cpu3
       write(13,518)tt*dunitt, (eshift(ialf,1), ialf=1,nalf), (jleads(ialf)*dunitj, ialf=1,nalf), occu
       write(16,518)tt, (eshift(ialf,1), ialf=1,nalf), occu
    end if
  else
    if (nalf .eq. 1) then
       write (6,518)tt, jleftu, jleftd, jleftu+jleftd, occu, occd, cpu4-cpu3
       write(13,518)tt*dunitt, eshift(1,1), jleftu*dunitj, jleftd*dunitj,              &
                    (jleftu+jleftd)*dunitj, (jleftu-jleftd)*dunitj, occ
       write(16,518)tt, eshift(1,1), occu, occd, occ
       write(28,518)tt, (dble(rhosys(ni,ni)), ni=1,nrho)
       if (lchkado) then
         write(30,518)tt, ((dble(rho(ni,nj,nfirst(itierchk)+inumchk)), ni=1,nrho), nj=1,nrho), &
                         ((dimag(rho(ni,nj,nfirst(itierchk)+inumchk)), ni=1,nrho), nj=1,nrho)
       end if
    else if (nalf .eq. 2) then
       write (6,518)tt, jleft, jright, occu, occd, cpu4-cpu3
       write(13,518)tt*dunitt, eshift(1,1), eshift(2,1), jleft*dunitj, jright*dunitj, &
                    (jleft+jright)*dunitj, occu, occd, occ
       write(16,518)tt, eshift(1,1), eshift(2,1), occu, occd, occ
    else  ! multi-leads
       write (6,518)tt, (jt(ialf,1), jt(ialf,2), ialf=1,nalf), occu, occd, cpu4-cpu3
       write(13,518)tt*dunitt, (eshift(ialf,1), eshift(ialf,2), ialf=1,nalf),                  &
                    (jt(ialf,1)*dunitj, jt(ialf,2)*dunitj, ialf=1,nalf), occu, occd, occ
       write(16,518)tt, (eshift(ialf,1), eshift(ialf,2), ialf=1,nalf), occu, occd, occ
    end if
  end if
  write(72,518)tt, esys*hbar, esb*hbar, etot*hbar
!
  call calcocc_detail
  write(33,518)tt, ((pocc(ni,nj), ni=1,norbs), nj=1,nspin)
!
  if (lbfield .and. nspin .gt. 1) then
     write(43,518)tt, dbfield * (1.d0 - dexp(-tt / tbfield)), sz, sz2  ! time, bfield, <Sz>, <Sz^2>
     call flush(43)
  end if
! output spin moment on each dot if lspin3d=true
  if (lspin3d) then
      if (lsparse) then
         nnz = nnz_spa(1)
         lni = ind_spa(1)
         call zmat_coo2dns(nnz, irow_spa(lni), icol_spa(lni), rho_spa(lni), cmat1,  &
                           nrho, nrho, nrho)
      else
         cmat1(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
      end if
      do iorbs=1,norbs
         call calcspinmoment(iorbs, sdot(1,iorbs), cmat1)
      end do
      write(61,518) tt, (sdot(ni,1), ni=1,3)
      if (norbs .eq. 2 .and. nspin .eq. 2) then
         write(62,518) tt, (sdot(ni,2), ni=1,3)
         call calctotalspin2(dspin2, cmat1)
         dtmp2 = dspin2(1) + dspin2(2) + dspin2(3)
         call calcspinproduct(1, 2, dtmp1, cmat1)
         call calcc12(dtmp3, cmat1)
         call calcc12p(dtmp4, cmat1)
         write(63,518) tt, dtmp1, (dspin2(ni), ni=1,3), dtmp2, dtmp3, dtmp4
      end if
  end if
!
! output rho(time)
!
!  write(29)tt, ((rhosys(ni,nj), ni=1,nrho), nj=1,nrho)
  519 format(' go ', f12.6, 2x, 10(e15.6e3, 2x))
  call flush(6)
  call flush(13)
  call flush(16)
  call flush(72)
  call flush(28)
  call flush(33)
  if (lchkado) call flush(30)
  if (lspin3d) then
     call flush(61)
     if (norbs .eq. 2 .and. nspin .eq. 2) then
        call flush(62)
        call flush(63)
     end if
  end if
!
! check whether converged in td-for-ground case 
!
  if (ltd4gr .and. nstep > 1 .and. tt > tref) then
     dtmp1 = 0.d0
     do ialf=1,nalf  
        do ispin=1,nspin
           dtmp1 = max(dtmp1, dabs(jt(ialf,ispin)))
        end do
     end do   
     write(6,*)'td4gr', tt, dtmp1, jconvg
     if (dtmp1 < jconvg) then
        write(6,*)
        write(6,*)'heom: convergence criterion for current achieved for td4ground'
        write(6,*)'heom: jconvg = ', jconvg, ' max_jt = ', dtmp1
        call flush(6)
!
        goto 911
     end if
  end if
!
  tt    = tt + dt
  nstep = nstep + 1
  if (tt .lt. tmax) goto 100
!
else if (tdmethod .eq. 1) then
  if (lresume) then
    nresume = min(nresume, 5) 
  end if
  call chebyshev
  tt = tmax
end if
911 continue
!
open(unit=40, file='hamil_sys.data', form='unformatted', status='unknown')
rewind(40)
write(40) ((hs(ni,nj), ni=1,nrho), nj=1,nrho)
close(40)
write(6,*)'heom: system Hamiltonian written to <hamil_sys.data> '
call flush(6)
!
call cpu_time(cpu6)
call free_rk4
call free_tmpmem_omp
!
deallocate(tmpmat, STAT=istat)
deallocate(sigmau, sigmad, STAT=istat)
deallocate(cmat1, STAT=istat)
!
close(13)
close(16)
close(72)
close(28)
close(29)
close(33)
if (lbfield .and. nspin .gt. 1) close(43)
if (lchkado) close(30)
if (linfcom) close(47)
if (lspin3d) then
   close(61)
   if (norbs .eq. 2 .and. nspin .eq. 2) then
      close(62)
      close(63)
   end if
end if
!
! save final results to file
!
do ialf=1,nalf
   do ispin=1,nspin
      call getengyshift(tt, ialf, ispin, eshift(ialf,ispin))
   end do
end do
eshift = eshift * hbar
write(6,*)
write(6,1010)' TDIV  ', tt, ((eshift(ialf,ispin), ispin=1,nspin), ialf=1,nalf), &
                                ((jt(ialf,ispin), ispin=1,nspin), ialf=1,nalf)
1010 format(A7, 2x, 32(e15.6e3, 1x))
!
write(6,*)
write(6,*)' final rho ', tt
call cmatout(nrho, nrho, rhosys, dmtmp1)
call flush(6)
!
if (lhf) then
   call calc_rsdm(rho, nrho, norbs, nspin, rsdm)
   write(6,*)
   write(6,*)'heom: final rsdm '
   do ispin=1,nspin
      write(6,*)'heom: ispin = ', ispin
      call cmatout(norbs, norbs, rsdm(1,1,ispin), dmtmp1_sp)
   end do
end if
!
! check final density matrix and save to files
!
call checkrho
200 continue
!
!write(6,*)
!write(6,1010)' TDRHO ', tt, ((dble(rhosys(ni,nj)), ni=1,nj), nj=1,nrho),         &
!                           ((dimag(rhosys(ni,nj)), ni=1,nj), nj=1,nrho)           
!
open(unit=14, file='indextable.sav', form='binary', status='unknown')
rewind(14)
write(14)ntable
write(14)(indextable(ntmp2), ntmp2=1,ntable)
close(14)
!
if (lsparse) then
   if (itdskip .eq. 0) then
      open(unit=50, file='rho_spa.sav', form='binary', status='unknown')
   else
      open(unit=50, file='rho_spa_td.sav', form='binary', status='unknown')
   end if
   rewind(50)
   write(50)norbs, nspin, ncor, ntier0, nalf
   write(50)nunk
   do lni=1,nunk
      nnz = nnz_spa(lni)
      lnj = ind_spa(lni)
      if (nnz .gt. 0) then
         write(50)(rho_spa(lnj-1+ni), ni=1,nnz)
      end if
   end do
   close(50)
else 
   if (itdskip .eq. 0) then
      open(unit=15, file='variables.sav', form='unformatted', status='unknown')
   else
      open(unit=15, file='variables_td.sav', form='unformatted', status='unknown')
   end if
   rewind(15)
   write(15)norbs, nspin, ncor, ntier, nalf
   write(15)nunk
   do lni=1,nunk
     write(15)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
   end do
   close(15)
!
   if (lhb) then
      open(unit=65, file='rho_hb.data', form='binary', status='unknown')
      rewind(65)
      write(65)norbs, nspin, ncor, ntier, nalf
      write(65)nmode_hb, ncor_hb, nbath_hb
      write(65)nunk
      do ibath=1,nbath_hb
         do lni=1,nunk
            write(65)((rho_hb(ni,nj,lni,ibath), ni=1,nrho), nj=1,nrho)
         end do
      end do
      close(65)
   end if
end if
!
! correlation function calculation 
!
ldos      = .false.
ljw_dos   = .true.
iorbs_dos = 1
ispin_dos = 1
!lfreq_dos = .false.
lfreq_dos = .true.  ! change default to true (note on 2022-05-15)
lread_dos = .true. 
freq_dos  = 0.d0
dt_dos    = 0.02d0
tmax_dos  = 100.d0
maxit_dos = 10000
crit_dos  = 1.d-6
rewind(5)
read(5, dos, end=301)
301 continue
if (ldos) then
   call corrfunc(iorbs_dos, ispin_dos, dt_dos, tmax_dos, ljw_dos, lfreq_dos, lread_dos, freq_dos, maxit_dos, crit_dos)
end if
302 continue
!
! Free memory
!
call freememory
!
! check IO status
!
!call checkio
!
call cpu_time(cpu2)
write(6,*)
write(6,*)' calculation done!             '
write(6,*)' cpu time for time evolution   ', cpu6 - cpu5
write(6,*)' cpu time (overall)            ', cpu2 - cpu1
!
call sj_timer_stop(1)
!
write(6,*)
write(6,*)'   general memory cost estimated ', memory,     ' MB '
write(6,*)'       rk4 memory cost estimated ', memork4,    ' MB '
write(6,*)'      bicg memory cost estimated ', memobicg,   ' MB '
write(6,*)'     tfqmr memory cost estimated ', memotfqmr,  ' MB '
write(6,*)'    jacobi memory cost estimated ', memojacobi, ' MB '
write(6,*)'   corrfun memory cost estimated ', memocf,     ' MB '
if (lhb) then
   write(6,*)' heat bath memory cost estimated ', memohb, ' MB '
end if
!
call date_and_time(date, time)
write(6,*)
write(6,1002)date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)
1002 format(/,' Job finished on ', /,                                   &
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2,/)
!
call flush(6)
close(5)
518 format(f12.4, 2x, 32(e15.6e3, 1x))
!
end program heom
