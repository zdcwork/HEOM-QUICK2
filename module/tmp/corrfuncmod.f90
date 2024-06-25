module corrfuncmod
implicit none
!
integer*8                    :: cfnunk, cfntable, cfnoprunk, cfnoprtable
integer,   save, allocatable :: cfindextable(:), cfoprtable(:)
integer*8, save, allocatable :: cfifirst(:), cfilast(:)
integer*8, save, allocatable :: cfnfirst(:), cfnlast(:)
integer*8, save, allocatable :: cfioprfirst(:), cfioprlast(:)
integer*8, save, allocatable :: cfnoprfirst(:), cfnoprlast(:)
!
real*8                       :: memorycf


end module corrfuncmod
