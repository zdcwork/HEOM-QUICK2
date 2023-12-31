subroutine allocate_tmpmem_omp
use matmod_omp
use tmpmatmod_omp
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
integer :: nrho2
!
allocate(dmtmp1_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp2_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp3_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp4_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp5_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp6_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp7_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(dmtmp8_omp(nrho,nrho,nthreads_omp), STAT=istat)
!
dmtmp1_omp = 0.d0
dmtmp2_omp = 0.d0
dmtmp3_omp = 0.d0
dmtmp4_omp = 0.d0
dmtmp5_omp = 0.d0
dmtmp6_omp = 0.d0
dmtmp7_omp = 0.d0
dmtmp8_omp = 0.d0
!
allocate(cmtmp1_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp2_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp3_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp4_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp5_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp6_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp7_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmp8_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmpa_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(cmtmpb_omp(nrho,nrho,nthreads_omp), STAT=istat)
!
cmtmp1_omp = czero
cmtmp2_omp = czero
cmtmp3_omp = czero
cmtmp4_omp = czero
cmtmp5_omp = czero
cmtmp6_omp = czero
cmtmp7_omp = czero
cmtmp8_omp = czero
cmtmpa_omp = czero
cmtmpb_omp = czero
!
allocate(zmtmp1_omp(nrho,nrho,nthreads_omp), STAT=istat)
allocate(zmtmp2_omp(nrho,nrho,nthreads_omp), STAT=istat)
!
zmtmp1_omp = czero
zmtmp2_omp = czero
!
!if (ltrun_der) then
!   nrho2 = nrho * nrho
!   allocate(zlmtmp1_omp(nrho2,nrho2,nthreads_omp), STAT=istat)
!   allocate(zlmtmp2_omp(nrho2,nrho2,nthreads_omp), STAT=istat)
!   allocate(zlmtmp3_omp(nrho2,nrho2,nthreads_omp), STAT=istat)
!   zlmtmp1_omp = czero
!   zlmtmp2_omp = czero
!   zlmtmp3_omp = czero
!end if
!
end subroutine allocate_tmpmem_omp
!
subroutine free_tmpmem_omp
use tmpmatmod_omp
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(dmtmp1_omp, dmtmp2_omp, dmtmp3_omp, dmtmp4_omp, &
           dmtmp5_omp, dmtmp6_omp, dmtmp7_omp, dmtmp8_omp, &
           STAT=istat)
!
deallocate(cmtmp1_omp, cmtmp2_omp, cmtmp3_omp, cmtmp4_omp, &
           cmtmp5_omp, cmtmp6_omp, cmtmp7_omp, cmtmp8_omp, &
           STAT=istat)
deallocate(cmtmpa_omp, cmtmpb_omp, STAT=istat)
!
deallocate(zmtmp1_omp, zmtmp2_omp, STAT=istat)
!
!if (ltrun_der) then
!   deallocate(zlmtmp1_omp, zlmtmp2_omp, zlmtmp3_omp, STAT=istat)
!end if
!
end subroutine free_tmpmem_omp
