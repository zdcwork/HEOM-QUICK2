module tmpmatmod_omp
implicit none
!
real*8,     save, allocatable :: dmtmp1_omp(:,:,:), dmtmp2_omp(:,:,:),  &
                                 dmtmp3_omp(:,:,:), dmtmp4_omp(:,:,:),  &
                                 dmtmp5_omp(:,:,:), dmtmp6_omp(:,:,:),  &    
                                 dmtmp7_omp(:,:,:), dmtmp8_omp(:,:,:)
complex*16, save, allocatable :: cmtmp1_omp(:,:,:), cmtmp2_omp(:,:,:),  &
                                 cmtmp3_omp(:,:,:), cmtmp4_omp(:,:,:),  &
                                 cmtmp5_omp(:,:,:), cmtmp6_omp(:,:,:),  &
                                 cmtmp7_omp(:,:,:), cmtmp8_omp(:,:,:),  &
                                 cmtmpa_omp(:,:,:), cmtmpb_omp(:,:,:)
complex*16, save, allocatable :: zmtmp1_omp(:,:,:), zmtmp2_omp(:,:,:)
complex*16, save, allocatable :: zlmtmp1_omp(:,:,:), zlmtmp2_omp(:,:,:), zlmtmp3_omp(:,:,:)
!
end module tmpmatmod_omp
