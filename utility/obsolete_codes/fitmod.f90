module fitmod
implicit none
!
logical                       :: lfitexp, lfitjob, detail
integer                       :: nfitexp
integer*8                     :: noriexp
real*8                        :: fitdpre
real*8,           allocatable :: fitb(:), fitgama(:), fitx(:)
real*8,           allocatable :: amatsu(:), ab(:), agama(:)
!
end module fitmod
