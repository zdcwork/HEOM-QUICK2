subroutine allocate_cf
use matmod
use matmod_omp
use corrfuncmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer    :: istat
real*8     :: dtmp1, dtmp2
integer*8  :: ltmp1, lnum
!
memocf = 0.d0
if (lsparse) then
   lnum = lunkcf_spa
else
   lnum = nunk * nrho**2 
end if
dtmp1  = dble(lnum * 16) / dble(1024**2)
!
allocate(a_ams(nrho,nrho), STAT=istat)
allocate(b_ams(nrho,nrho), b_ams_dag(nrho,nrho), STAT=istat)
!
if (lsparse) then
   allocate(brhoh_spa(lunkcf_spa), brhoa_spa(lunkcf_spa), STAT=istat)
   memocf = memocf + 2.d0 * dtmp1   
   if (.not. lfreq_dos) then
      allocate(brho0_spa(lunkcf_spa), brho1_spa(lunkcf_spa), STAT=istat)
      allocate(rhop_spa(lunkcf_spa), rhoprhs_spa(lunkcf_spa), STAT=istat)
      allocate(rhoq_spa(lunkcf_spa), brhorhs_spa(lunkcf_spa), STAT=istat)
      memocf = memocf + 6.d0 * dtmp1
      if (ljw_dos) then
         allocate(bdrhoh_spa(lunkcf_spa), bdrhoa_spa(lunkcf_spa), STAT=istat)
         memocf = memocf + 2.d0 * dtmp1
      end if
   else
      ltmp1 = lunkcf_spa * 2 * 2
      dtmp2 = dble(ltmp1 * 8) / dble(1024**2)
      allocate(qmrvecs_cf(ltmp1,9), STAT=istat)
      memocf = memocf + 9.d0 * dtmp2
   end if
else
   allocate(brhoh(nrho,nrho,nunk), STAT=istat)
   allocate(brhoa(nrho,nrho,nunk), STAT=istat)
   memocf = memocf + 2.d0 * dtmp1
   if (.not. lfreq_dos) then
      allocate(brho0(nrho,nrho,nunk), STAT=istat)
      allocate(brho1(nrho,nrho,nunk), STAT=istat)
      allocate(brhorhs(nrho,nrho,nunk), STAT=istat)
      allocate(rhoprhs(nrho,nrho,nunk), STAT=istat)
      allocate(rhop(nrho,nrho,nunk), STAT=istat)
      memocf = memocf + 5.d0 * dtmp1
      if (ljw_dos) then
         allocate(bdrhoh(nrho,nrho,nunk), STAT=istat)
         allocate(bdrhoa(nrho,nrho,nunk), STAT=istat)
         memocf = memocf + 2.d0 * dtmp1
      end if
   else 
      ltmp1 = (nunk * 2) * nrho**2 * 2   ! real and imaginary parts for both brhoh and brhoa
      dtmp2 = dble(ltmp1 * 8) / dble(1024**2)  
      allocate(qmrvecs_cf(ltmp1,9), STAT=istat)   
      memocf = memocf + 9.d0 * dtmp2
   end if
end if
!
if (istat .eq. 0) then
   write(6,*)'allocate_cf: memory allocated for corrfunc ', memocf, ' MB '
else 
   write(6,*)'allocate_cf: fail memory allocation! '
   stop
end if
!
call allocate_tmpmem_omp
!
return
end subroutine allocate_cf
!
subroutine free_cf
use matmod
use matmod_omp
use corrfuncmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: istat
!
deallocate(a_ams, b_ams, b_ams_dag, STAT=istat)
!
if (lsparse) then
   deallocate(brhoh_spa, brhoa_spa, STAT=istat)
   deallocate(nnzcf_spa, indcf_spa, irowcf_spa, icolcf_spa, STAT=istat)
   deallocate(irowcf_hsys, icolcf_hsys, cvalcf_hsys, STAT=istat)
else
   deallocate(brhoh, brhoa, STAT=istat)
end if
!
if (lfreq_dos) then
   deallocate(qmrvecs_cf, STAT=istat)
else 
   if (lsparse) then
      deallocate(brho0_spa, brho1_spa, brhorhs_spa, rhop_spa, rhoq_spa, rhoprhs_spa, STAT=istat)
      if (ljw_dos) then
         deallocate(bdrhoh_spa, bdrhoa_spa, STAT=istat)
      end if
   else
      deallocate(brho0, brho1, brhorhs, rhop, rhoprhs, STAT=istat)
      if (ljw_dos) then
         deallocate(bdrhoh, bdrhoa, STAT=istat)
      end if
   end if
end if
!
if (ltrun_der) then
   deallocate(zmatder_cf, dvalder_cf, STAT=istat)
end if
!
call free_tmpmem_omp
!
end subroutine free_cf
