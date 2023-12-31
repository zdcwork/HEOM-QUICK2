module hbmod
implicit none
!
logical                       :: lhb, lscba
integer                       :: npade_hb, ndrude_hb, nmode_hb, itype_hb
integer                       :: ncor_hb
integer                       :: nn_hb, nk_hb, np_hb, nbath_hb, nkp_hb
real*8                        :: dbeta_hb, dinvbeta_hb
real*8,     save, allocatable :: dwidth_hb(:), dcouple_hb(:,:), dcenter_hb(:)
complex*16, save, allocatable :: zcoef_hb(:), zcoefp_hb(:)
complex*16, save, allocatable :: zpeta_hb(:), zppole_hb(:), zpeigv_hb(:)
complex*16, save, allocatable :: cb_hb(:,:), cd_hb(:,:), cgamma_hb(:,:)
complex*16, save, allocatable :: qa_hb(:,:,:)
!
integer*8                     :: nunk_hb, lunk_hb_spa
complex*16, save, allocatable :: rho_hb(:,:,:,:), rho0_hb(:,:,:,:),     &
                                 rhotmp_hb(:,:,:,:), rhorhs_hb(:,:,:,:)
!
end module hbmod
