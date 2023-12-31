module tmpmatmod
implicit none
!
integer,    save, allocatable :: iatmp1(:), iatmp2(:)
!
real*8,     save, allocatable :: denmbr(:,:), denmbi(:,:), denmcr(:,:), denmci(:,:)
!
real*8,     save, allocatable :: dmtmp1(:,:), dmtmp2(:,:)
real*8,     save, allocatable :: dmtmp3(:,:), dmtmp4(:,:)
real*8,     save, allocatable :: dmtmp5(:,:), dmtmp6(:,:)
!
complex*16, save, allocatable :: cmtmp1(:,:), cmtmp2(:,:)
complex*16, save, allocatable :: cmtmp3(:,:), cmtmp4(:,:)
complex*16, save, allocatable :: cmtmp5(:,:), cmtmp6(:,:)
complex*16, save, allocatable :: cmtmp7(:,:), cmtmp8(:,:), cmtmp9(:,:)
complex*16, save, allocatable :: cmtmpa(:,:), cmtmpb(:,:)
complex*16, save, allocatable :: cmtmp10(:,:), cmtmp11(:,:), cmtmp12(:,:), cmtmp13(:,:)
complex*16, save, allocatable :: cmtmp14(:,:), cmtmp15(:,:), cmtmp16(:,:), cmtmp17(:,:)
complex*16, save, allocatable :: cmtmp18(:,:)
complex*16, save, allocatable :: zmtmp1(:,:), zmtmp2(:,:)  
!
real*8,     save, allocatable :: dmtmp1_sp(:,:), dmtmp2_sp(:,:)
!
! Liouville space matrix
!
real*8,     save, allocatable :: dlmat1(:,:), dlmat2(:,:), dlmat3(:,:), dlmat4(:,:),  &
                                 dlmat5(:,:), dlmat6(:,:)
!
! Liouville space vector
!
real*8,     save, allocatable :: dlvec1(:), dlvec2(:), dlvec3(:), dlvec4(:)
real*8,     save, allocatable :: dlvec5(:), dlvec6(:), dlvec7(:), dlvec8(:)
!
complex*16, save, allocatable :: zlmtmp1(:,:), zlmtmp2(:,:)
complex*16, save, allocatable :: zlmtmp3(:,:), zlmtmp4(:,:)
!
end module tmpmatmod
