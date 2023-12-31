module sparsemod
implicit none
!
integer                       :: nnzero_spa   ! nonzero elements in each ams matrix
!
! amsall in CSR (compressed sparse row) format (3-array) variation
!
real*8,     save, allocatable :: ams_spa(:,:,:)
integer,    save, allocatable :: ja_spa(:,:,:)
integer,    save, allocatable :: ia_spa(:,:,:)
!
integer                       :: iread_spa    ! iread_spa = 0: cannot read sparsity information from file
                                              !           = 1: read-in sparsity information is different but compatible 
                                              !           = 2: read-in sparsity information is same as present job
!
! sparsity format of effective Hamiltonian 
!
integer                       :: nnz_hsys
integer,    save, allocatable :: irow_hsys(:), icol_hsys(:)
complex*16, save, allocatable :: cval_hsys(:)
!
integer                       :: nnztd_hsys
integer,    save, allocatable :: irowtd_hsys(:), icoltd_hsys(:)
complex*16, save, allocatable :: cvaltd_hsys(:)
!
! sparse representation of all ADOs
!
integer*8                     :: lunk_spa 
integer,    save, allocatable :: nnz_spa(:)
integer*8,  save, allocatable :: ind_spa(:)
integer,    save, allocatable :: irow_spa(:), icol_spa(:)
complex*16, save, allocatable :: rho_spa(:)
complex*16, save, allocatable :: rhs_spa(:)
complex*16, save, allocatable :: rho0_spa(:), rhotmp_spa(:), rhorhs_spa(:)
!
integer*8                     :: lunk_spa0 
integer,    save, allocatable :: nnz_spa0(:)
integer*8,  save, allocatable :: ind_spa0(:)
integer,    save, allocatable :: irow_spa0(:), icol_spa0(:)
complex*16, save, allocatable :: rho_spa0(:)
!
! sparse representation of all variables for system correlation functions
!
integer                       :: nnzcf_hsys
integer,    save, allocatable :: irowcf_hsys(:), icolcf_hsys(:)
complex*16, save, allocatable :: cvalcf_hsys(:)
!
integer*8                     :: lunkcf_spa
integer,    save, allocatable :: nnzcf_spa(:)
integer*8,  save, allocatable :: indcf_spa(:)
integer,    save, allocatable :: irowcf_spa(:), icolcf_spa(:)
complex*16, save, allocatable :: brhoh_spa(:), brhoa_spa(:)
!
complex*16, save, allocatable :: brho0_spa(:), brho1_spa(:)
complex*16, save, allocatable :: bdrhoh_spa(:), bdrhoa_spa(:)
complex*16, save, allocatable :: rhop_spa(:), rhoq_spa(:)
complex*16, save, allocatable :: brhorhs_spa(:), rhoprhs_spa(:)
!
end module sparsemod
