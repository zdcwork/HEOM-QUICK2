module corrfuncmod
implicit none
!
integer, save                 :: icont_cf, nresume_cf
logical, save                 :: lresume_cf, lplus_cf
!
integer, save                 :: iorbs_dos, ispin_dos, maxit_dos
real*8,  save                 :: dt_dos, tmax_dos, freq_dos, crit_dos
logical, save                 :: ljw_dos, lfreq_dos, lread_dos
!
integer                       :: a_ams_isgn, a_ams_iorbs, a_ams_ispin
integer                       :: b_ams_isgn, b_ams_iorbs, b_ams_ispin
complex*16, save, allocatable :: a_ams(:,:)
complex*16, save, allocatable :: b_ams(:,:), b_ams_dag(:,:)
!
! time propagation
!
complex*16, save, allocatable :: brhoh(:,:,:), brhoa(:,:,:)     ! for Brho
complex*16, save, allocatable :: brho0(:,:,:), brhorhs(:,:,:)
complex*16, save, allocatable :: rhoprhs(:,:,:)
complex*16, save, allocatable :: brho1(:,:,:), rhop(:,:,:)
complex*16, save, allocatable :: bdrhoh(:,:,:), bdrhoa(:,:,:)   ! for Brho^\dag
!
! linear problem by iterative methods
!
real*8,     save, allocatable :: qmrvecs_cf(:,:)
!
end module corrfuncmod
