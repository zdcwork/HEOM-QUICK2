module matmod
implicit none
!
integer,    save, allocatable :: indextable(:)
integer,    save, allocatable :: filtertable(:)
integer,    save, allocatable :: oprtable(:)
integer*8,  save, allocatable :: indexcoef(:), index_diff_ref(:), index_coef_ref(:)
integer*1,  save, allocatable :: itypecoef(:)
integer*1,  save, allocatable :: ifactfcoef(:), ifactrcoef(:)
integer*8,  save, allocatable :: ioprindex(:), ioprref(:)
logical,    save, allocatable :: idiffcoef(:)
integer*8,  save, allocatable :: ifirst(:), ilast(:)     ! first and last position in indextable for each tier
integer*8,  save, allocatable :: nfirst(:), nlast(:)     ! first and last position in array rho
integer*8,  save, allocatable :: ioprfirst(:), ioprlast(:)  
integer*8,  save, allocatable :: noprfirst(:), noprlast(:)  
integer,    save, allocatable :: indexdraw(:)
!
logical                       :: lprths                  ! T: print out info when executing 'calchs'; F: silent mode 
complex*16, save, allocatable :: hs(:,:)                 ! Hamiltonian with dimension (nrho, nrho) (should be Hermitian)
complex*16, save, allocatable :: dhs(:,:)                ! delta Hamiltonian with dimension (nrho, nrho) (should be Hermitian)
complex*16, save, allocatable :: rho(:,:,:)              ! dimension (nrho, nrho, all_possible_cases_for_balls_and_drawers)
complex*16, save, allocatable :: rhosys(:,:)             ! reduced density matrix 
!
real*8,     save, allocatable :: ams(:,:)                ! a_{ms} and a^\dag_{ms}
real*8,     save, allocatable :: amsall(:,:,:,:)
real*8,     save, allocatable :: amsori(:,:,:,:)
integer,    save, allocatable :: nnz_ams(:,:)
integer,    save, allocatable :: row_ams(:,:,:), col_ams(:,:,:)
integer,    save, allocatable :: rowams(:,:,:), colams(:,:,:)
integer,    save, allocatable :: sgnams(:,:,:)           ! 0 for plus sign, and others for minus sign
real*8,     save, allocatable :: valams(:,:,:)
complex*16, save, allocatable :: cval_ams(:,:,:), zamsall(:,:,:,:)
!
logical                       :: lreadini, lreadmat, lreadfile
real*8,     save, allocatable :: rhodiag(:)
complex*16, save, allocatable :: znms(:,:,:,:)           ! occupancy operators
!
! rk4 memory
!
complex*16, save, allocatable :: rho0(:,:,:), rhotmp(:,:,:), rhorhs(:,:,:)
!
integer,    save, allocatable :: iredex(:), morbs(:), mspin(:), mmats(:), ilead(:), mpm(:), msopr(:)
integer,    save, allocatable :: mdrude(:), mfreq(:), mcor(:), mopr(:)
integer,    save, allocatable :: mnfff(:), mmfff(:)
integer,    save, allocatable :: jpm(:), jspin(:), jorbs(:), jalf(:), jredex(:), jomit(:)
integer,    save, allocatable :: kmats(:), kpmats(:)
real*8,     save, allocatable :: dmemexp(:)              ! memory components 
!
real*8,     save, allocatable :: eleadinfty(:,:), eleadtime(:,:)
real*8,     save, allocatable :: dimpf(:,:,:,:,:)
real*8,     save, allocatable :: dbsqrt(:,:,:,:,:), dbinvsqrt(:,:,:,:,:)
complex*16, save, allocatable :: cgama(:,:,:,:,:), cb(:,:,:,:,:), cd(:,:,:,:,:), cgamma(:)
!
real*8,     save, allocatable :: dmatsu(:,:), dbeta(:), dlwidth(:,:)    ! dlwidth(norbs,nspin) overall linewidth for an orbital
real*8,     save, allocatable :: dwp(:,:,:), dbwp(:,:,:), gammap(:,:,:), dxip(:,:,:)
real*8,     save, allocatable :: dlw_sp(:,:,:)
complex*16, save, allocatable :: zxip(:,:,:)
!
complex*16, save, allocatable :: cpeigv(:), cpcoef(:), cppole(:,:)
!
real*8,     save, allocatable :: dpm(:)                  ! dpm(1)=1.d0, dpm(2)=-1.d0
real*8,     save, allocatable :: dipm(:)                 ! +1 for diamond=+, -1 for diamond=-
real*8,     save, allocatable :: jt(:,:), jleads(:)      ! jt(nalf,nspin), jleads(nalf)
real*8,     save, allocatable :: pocc(:,:)
!
! Chebyshev memory
!
integer                       :: jcheby                  ! # of current points output per dt 
integer*8                     :: ncheby
real*8                        :: cheby_delta, tdelta
real*8                        :: cheby_dpre
!
! POP related memory
!
complex*16,       allocatable :: cpop(:,:,:,:,:,:,:)
complex*16,       allocatable :: vech(:,:), vinvh(:,:)
real*8,           allocatable :: valh(:)
real*8,           allocatable :: lpopr(:,:), lpopi(:,:)
real*8,           allocatable :: lbpopr(:,:), lbpopi(:,:), ldpopr(:,:), ldpopi(:,:)
! Resume
real*8                        :: tresume
!
! time-dependent energy shift (applied voltage) (see getenergyshift.f90 for details)
!
integer                       :: fieldtype               ! 0 for exponential, -1 for delta-pulse
logical                       :: lreadomega              ! .true.  frequency is read for voltage applied to electrodes
                                                         ! .false. charactersitic time is read
logical                       :: leshift_dc              ! .true.  add dc component of energy shift  
real*8                        :: t_off                   ! turn off all applied voltage at t_off
!
! time-dependent energy shift for QD (system) 
!
logical                       :: forcesteady             ! .true. : override all limits for steady jobs
logical                       :: ldfield_hub              
integer                       :: gatetype
real*8                        :: gpara1, gpara2, gpara3, gpara4
integer                       :: dfieldtype              ! 0 for exponential, 1 for sinusoidal
! for norbs=1 .and. nspin=1
real*8                        :: dedot, aD, wD           ! time-dependent external field for the system
real*8                        :: vdot                    ! initial field (Dirac delta function)
! for norbs=1 .and. nspin=2
real*8                        :: tflip, wflip            ! spin-filp amplitude and frquency
real*8                        :: egate, ton_gate, toff_gate
! for norbs=2 .and. nspin=2  
real*8                        :: dedot1, dedot2, wdot1, wdot2, tdot12, jdot12
! for norbs=2 .and. npsin=1
real*8                        :: tson, tsoff, pha
! for Hubbard model
real*8                        :: gamadot
!
! system Hamiltonian
logical                       :: lreadhs, lfixhs, lgenhs ! lgenhs = .true.: construct hamiltonian using a general input form
! for Hubbard model
logical                       :: lhubbard, lhf, lspinsymm, lorbsymm, lspindiff
real*8                        :: vcoup
!
! reduced single-electron density matrix
complex*16, save, allocatable :: rsdm(:,:,:)
!
! LEAD property
!
logical,    save              :: onelead
logical,    save              :: lequileads              ! leads are equivalent: temperature, spectral density are
                                                         !                       identical and no temperature
integer,    save              :: ispectral
integer,    save              :: nleads                  ! number of leads designated in the original input file
integer                       :: chksgn, chkalf, chkspin, chkorbs, chkorbsp
logical                       :: chkcd, skipcorr
logical                       :: chkcorr                 ! check correlation function expansion
real*8                        :: egaul, egaur            ! for Gaussian-type dos(energy)
real*8                        :: esincl, esincr          ! for Square-sinc-type dos(energy)
real*8,     save, allocatable :: amps(:,:), tchar(:,:)   ! amplitude and characteristic decay time for each lead
real*8,     save, allocatable :: omegas(:,:)             ! spherical frequency of driving sinusoidal voltage
real*8,     save, allocatable :: amps_dc(:,:)            ! dc components
!
! Time-evolution related variables
!
integer                       :: tdmethod                ! 0 for 4th-order Runge-Kutta, 1 for Chebyshev
logical                       :: ltd4gr                  ! true, use time-propagation to solve ground state
real*8                        :: jconvg                  ! when ltd4gr=T, jconvg is the convergence criterion for current
real*8                        :: tref                    ! check jconvg criterion after t > tref
real*8                        :: tnow                    ! for transient system spectral function
!
! Frequency domain discretization info
!
real*8,     save, allocatable :: wl0(:), wl1(:)
real*8,     save, allocatable :: wkln(:), dkln(:) 
complex*16, save, allocatable :: focc(:,:,:,:)           ! f^\sigma_{\alpha,s}(\omega_m)
complex*16, save, allocatable :: spectral(:,:,:,:,:,:)   ! \Gamma_{\sigma,\alpha,s,\mu,\nu}(\omega_m)
complex*16, save, allocatable :: zlwmat(:,:,:,:)         ! complex linewidth matrix (see <makespectral>)
!
! Frequency domain approach to system dynamic properties
!
logical,    save              :: lfreq, lplus
real*8,     save              :: dfreq
!
! Lorentzian fitting for residual reservoir spectral function
!
real*8,     save, allocatable :: dlor_coef(:,:), dlor_cent(:,:), dlor_width(:,:)
!
! Magnetic field
!
logical                       :: lbfield                 ! whether a turn-on B-field is applied
real*8                        :: dbfield                 ! effective B-field strength, H = - B-field * Sz
real*8                        :: tbfield                 ! turn on time for B-field
!
! Lorentzian fitting for reservoir spectral function (not residual)
! 
logical                       :: lband, lband_same
integer                       :: nband
real*8,     save, allocatable :: band_coef(:,:,:), band_cent(:,:,:), band_width(:,:,:)   ! (nband, nspin, malf)
                                                         ! band_cent and band_width are of dimension of energy 
                                                         ! band_coef is dimensionless 
!
! 3d spin
!
logical                       :: lspin3d, lbfield3d, lbsite, lzfs
real*8                        :: d_xx, d_yy, d_zz        ! Dxx, Dyy, Dzz for zero-field splitting (ZFS)
real*8                        :: d_xx_1, d_yy_1, d_zz_1
real*8                        :: d_xx_2, d_yy_2, d_zz_2
real*8,     save, allocatable :: dbf3d(:,:)              ! dbf3d(ixyz,iorbs)
complex*16, save, allocatable :: pauli(:,:,:)            ! Pauli matrices
complex*16, save, allocatable :: sopr(:,:,:,:)           ! S(irho,jrho,iorbs,ixyz)
real*8,           allocatable :: sdot(:,:)               ! <Sx>, <Sy>, <Sz> for dots
!
! local temperature
!
logical                       :: linfcom, lrun_loc
integer                       :: nstep_loc
real*8                        :: energy_loc, entropy_loc, dinfcom_loc
!
! Truncation based on derivatives of ADOs
!
integer                       :: nnzeig, nnzeig_cf(2)
real*8,     save, allocatable :: dvalder(:), dvalder_cf(:,:)
complex*16, save, allocatable :: zmatder(:,:), zmatder_cf(:,:,:)
!
! delta voltage
!
logical                       :: ldeltav
real*8,     save, allocatable :: deltav_amps(:,:)
!
! Fano-function-fitting (FFF) scheme
!
integer,    save              :: nfff, npwfff, numfff
real*8,     save              :: rfff
integer,    save, allocatable :: mpfff(:), nfff1d(:), mfff1d(:)
real*8,     save, allocatable :: afff(:,:), bfff(:), temfff(:), gfff(:,:)
complex*16, save, allocatable :: cfff(:,:,:,:,:), dfff(:,:,:,:), thefff(:,:,:,:,:)
complex*16, save, allocatable :: ctheta(:,:,:,:) 
!
! Jacobi iteration
!
integer                       :: nnzjac   
real*8,     save, allocatable :: dvaljac(:)
complex*16, save, allocatable :: zmatjac(:,:)
real*8                        :: ejac
!
! Adiabatic HEOM
!
integer,    save              :: ntier_slow, ndrawer_slow
integer,    save, allocatable :: ns2f(:), icor_slow(:), icor_fast(:)
integer,    save, allocatable :: ms2f(:), mslow(:)
integer,    save, allocatable :: nnz_zgf(:), irow_zgf(:,:), icol_zgf(:,:)
integer,    save, allocatable :: nnz_sop(:), irow_sop(:,:), icol_sop(:,:)  ! sop stands for system operators
integer,    save, allocatable :: nnz_lvl(:), irow_lvl(:,:), icol_lvl(:,:)
integer,    save, allocatable :: nnz_lvr(:), irow_lvr(:,:), icol_lvr(:,:)
integer,    save, allocatable :: ind_lvl(:,:), ind_lvr(:,:), nnc_lvl(:,:), nnc_lvr(:,:)
integer,    save, allocatable :: nnz_lvlsop(:), irow_lvlsop(:,:), icol_lvlsop(:,:)
integer,    save, allocatable :: nnz_lvrsop(:), irow_lvrsop(:,:), icol_lvrsop(:,:)
integer,    save, allocatable :: npermu(:,:,:), npermu_fact(:,:)
integer,    save, allocatable :: nlvrow(:), nlvcol(:), nnrc2v(:,:), ncrc2v(:,:)
real*8,     save, allocatable :: dgama_drawer(:)
complex*16, save, allocatable :: zgefnl(:,:,:), zgefnr(:,:,:)
complex*16, save, allocatable :: zgelvl(:,:,:), zgelvr(:,:,:)              ! lv stands for Liouville
complex*16, save, allocatable :: zgefnl_sop(:,:,:), zgefnr_sop(:,:,:)
complex*16, save, allocatable :: zgelvl_sop(:,:,:), zgelvr_sop(:,:,:)
complex*16, save, allocatable :: cval_zgfl(:,:), cval_zgfr(:,:)
complex*16, save, allocatable :: cval_sopl(:,:), cval_sopr(:,:)
complex*16, save, allocatable :: cval_lvl(:,:), cval_lvr(:,:)
complex*16, save, allocatable :: cval_lvlsop(:,:), cval_lvrsop(:,:)
!
! system projector
!
integer,    save              :: nstate_proj
integer,    save, allocatable :: istate_proj(:)
complex*16, save, allocatable :: cproj(:,:)
complex*16, save, allocatable :: cproj_bar(:,:)   ! complement of cproj:  cproj+cproj_bar=I 
!
end module matmod
