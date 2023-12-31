subroutine initialize2(length0)
use random
use matmod
use tmpmatmod
use hbmod
implicit none
!
include '../include/sizes'
include '../include/common'
!
integer*8, intent(in) :: length0(*)
integer               :: istat, ibath
integer               :: ni, nj
integer               :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
integer*8             :: lni, lnj, lnk, lunk, ltmp1
logical               :: lsmall 
real*8                :: dtmp1
!
if (lsparse) then
   call initialize2_spa(length0)
   return
end if
!
allocate(hs(nrho,nrho), dhs(nrho,nrho), rhosys(nrho,nrho), STAT=istat)
allocate(rho(nrho,nrho,nunk), STAT=istat)
allocate(rsdm(norbs,norbs,nspin), STAT=istat)
!
write(6,*)
if (istat == 0) then
 write(6,*)' OK! Memory allocation for rho successful '
 write(6,*)' memory allocated : ', dble(nunk * 16 * nrho**2) / dble(1024**2), ' MB'
else
 write(6,*)' error! Memory allocation for rho failed  '
 write(6,*)' memory needed    : ', dble(nunk * 16 * nrho**2) / dble(1024**2), ' MB'
 stop
end if
call flush(6)
memory = memory + dble(nunk * 16 * nrho**2) / dble(1024**2)
!
! Initialize density matrix (and auxiliary density matrix) here, and
! note that normalization constraint for rho^0 should be satisfied at this stage.
!
open(unit=15, file='variables.sav', form='unformatted', status='unknown')
rewind(15)
read(15)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5  ! corresponding to norbs, nspin, ncor, jtier, and nalf, respectively
read(15)lunk                               ! number of unknowns up to jtier (jtier <= ntier0)
!
if (ntmp1 .ne. norbs .or. ntmp2 .ne. nspin .or. ntmp3 .ne. ncor .or. &
    ntmp5 .ne. nalf  .or. ntmp4 .le. 1) then 
  write(6,*)
  write(6,*)' error input file <variables.sav> '
  stop
end if
!
lsmall = .false.
if (ntmp4 < ntier0) lsmall = .true.
!
lnj = 0
do ni=1,min0(ntmp4,ntier0)
   lnj = lnj + length0(ni)
end do
if (lsmall) then
   if (lnj .ne. lunk) then
      write(6,*)
      write(6,*)'initialize2: error! incompatible <variables.sav> ', ntmp4, ntier0, lnj, lunk
      stop
   end if
   do lni=1,lunk
      read(15)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
   end do
   rho(1:nrho, 1:nrho, lunk+1:nunk) = czero
else 
   if (lnj .ne. nunk) then
      write(6,*)
      write(6,*)'initialize2: error! incompatible <variables.sav> ', ntmp4, ntier0, lnj, nunk
      stop
   end if
   do lni=1,nunk
      read(15)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
   end do
end if
close(15)
!
write(6,*)
write(6,*)'initialize2: <variables.sav> read successfully '
call flush(6)
!
! Initialize system Hamiltonian hs here, and note that hs should be a Hermitian matrix
!
hs(1:nrho,1:nrho) = czero
!
! Heat bath
!
if (lhb) then
   ltmp1 = nunk * nrho**2 * nbath_hb * 2  ! both real and imaginary parts
   dtmp1  = dble(ltmp1 * 8) / dble(1024**2)
   memohb = 0.d0
!
   allocate(rho_hb(nrho,nrho,nunk,nbath_hb), STAT=istat)
   memohb = memohb + dtmp1 
   if (istat .ne. 0) then
      write(6,*)'initialize2: error! memory allocation fail 2 '
      write(6,*)'memory desired ', memohb, ' MB'
      stop
   end if
!
   open(unit=65, file='rho_hb.data', form='binary', status='unknown')
   rewind(65)
   read(65)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5  ! norbs, nspin, ncor, jtier, and nalf
   if (ntmp1 .ne. norbs .or. ntmp2 .ne. nspin .or. ntmp3 .ne. ncor .or. &
       ntmp5 .ne. nalf  .or. ntmp4 .le. 1) then 
      write(6,*)
      write(6,*)'initialize2: error! wrong input file <rho_hb.data> 1 '
      write(6,*)norbs, nspin, ncor, ntier, nalf
      write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5
      stop
   end if
   read(65)ntmp1, ntmp2, ntmp3                ! nmode_hb, ncor_hb, nbath_hb
   if (ntmp1 .ne. nmode_hb .or. ntmp2 .ne. ncor_hb .or. ntmp3 .ne. nbath_hb) then
      write(6,*)
      write(6,*)'initialize2: error! wrong input file <rho_hb.data> 2 '
      write(6,*)nmode_hb, ncor_hb, nbath_hb
      write(6,*)ntmp1, ntmp2, ntmp3
      stop
   end if
   read(65)lunk                               ! number of unknowns up to jtier (jtier <= ntier)
!
   lsmall = .false.
   if (ntmp4 < ntier) lsmall = .true.
   lnj = 0
   do ni=1,min0(ntmp4,ntier)
      lnj = lnj + length0(ni)
   end do
   if (lsmall) then
      if (lnj .ne. lunk) then
         write(6,*)'initialize2: error! incompatible <rho_hb.data> '
         write(6,*)ntmp4, ntier, lnj, lunk
         stop
      end if
      do ibath=1,nbath_hb     
         do lni=1,lunk
            read(65)((rho_hb(ni,nj,lni,ibath), ni=1,nrho), nj=1,nrho)
         end do
         rho_hb(1:nrho,1:nrho,lunk+1:nunk,ibath) = czero
      end do
   else 
      if (lnj .ne. nunk) then
         write(6,*)'initialize2: error! incompatible <rho_hb.data> '
         write(6,*)ntmp4, ntier, lnj, nunk
         stop
      end if
      do ibath=1,nbath_hb
         do lni=1,nunk
            read(65)((rho_hb(ni,nj,lni,ibath), ni=1,nrho), nj=1,nrho)
         end do
      end do
   end if
   close(65)
!
   write(6,*)
   write(6,*)'initialize2: <rho_hb.data>   read successfully '
   call flush(6)
end if
!
return
end subroutine initialize2
