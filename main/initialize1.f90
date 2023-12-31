subroutine initialize1
use random
use matmod
use tmpmatmod
use hbmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer             :: istat
integer             :: ni, nj, itier
integer*8           :: lni, ltmp1
logical             :: lherm
real*8              :: dtmp1
!
namelist / readini / lreadini, lreadmat, lreadfile
!
if (lsparse) then
   call initialize1_spa
   return
end if
!
allocate(rho(nrho,nrho,nunk), STAT=istat)
allocate(rhodiag(nrho), STAT=istat)
allocate(hs(nrho,nrho), dhs(nrho,nrho), rhosys(nrho,nrho), STAT=istat)
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
!
! Initialize density matrix (and auxiliary density matrix) here, and
! note that normalization constraint for rho^0 should be satisfied at this stage.
!
rho(1:nrho,1:nrho,1:nunk) = czero
rho(nrho,nrho,1) = cunity
!
!do ni=1,nrho
!   rho(ni,ni,1) = cunity / dble(nrho)
!end do
!
if (norbs .eq. 2 .and. nspin .eq. 2) then
   if (iguess .eq. 0) goto 100
!
   rho(1:nrho,1:nrho,1:nunk) = czero
   write(6,*)
!
   if (iguess .eq. -1) then  ! empty state
      rho(1,1,1) = cunity
      write(6,*)'initialize1: initial rho is set to empty state '
!
   else if (iguess .eq. 1) then  ! S=0 (singlet state)
      rho( 7, 7,1) =  dcmplx(5.d-1, 0.d0)
      rho(10,10,1) =  dcmplx(5.d-1, 0.d0)
      rho( 7,10,1) = -dcmplx(5.d-1, 0.d0)
      rho(10, 7,1) = -dcmplx(5.d-1, 0.d0)
      write(6,*)'initialize1: initial rho is set to singlet state (S=0) '
!
   else if (iguess .eq. 2) then  ! S=1, Ms=0 (triplet state)
      rho( 7, 7,1) =  dcmplx(5.d-1, 0.d0)
      rho(10,10,1) =  dcmplx(5.d-1, 0.d0)
      rho( 7,10,1) =  dcmplx(5.d-1, 0.d0)
      rho(10, 7,1) =  dcmplx(5.d-1, 0.d0)
      write(6,*)'initialize1: initial rho is set to triplet state (S=1, Ms=0) '
!
   else if (iguess .eq. 3) then  ! S=1, Ms=1 (triplet state)
      rho(6,6,1) = cunity
      write(6,*)'initialize1: initial rho is set to triplet state (S=1, Ms=1) '
!
   else if (iguess .eq. 4) then  ! S=1, Ms=-1 (triplet state)
      rho(11,11,1) = cunity
      write(6,*)'initialize1: initial rho is set to triplet state (S=1, Ms=-1) '
!
   else if (iguess .eq. 5) then  ! S=1, Ms=1,-1 (triplet state) 1:1 mixture of Ms=1 and Ms=-1
      rho( 6, 6,1) =  dcmplx(5.d-1, 0.d0)
      rho(11,11,1) =  dcmplx(5.d-1, 0.d0)
      write(6,*)'initialize1: initial rho is set to triplet state (S=1, Ms=1,-1) '
      write(6,*)'initialize1: Ms=1, -1 are mixed at the level of density matrix  '
!
   else if (iguess .eq. 6) then  ! S=1, Ms=1,-1 (triplet state) mixing wavefunction
      rho( 6, 6,1) =  dcmplx(5.d-1, 0.d0)
      rho(11,11,1) =  dcmplx(5.d-1, 0.d0)
      rho( 6,11,1) =  dcmplx(5.d-1, 0.d0)
      rho(11, 6,1) =  dcmplx(5.d-1, 0.d0)
      write(6,*)'initialize1: initial rho is set to triplet state (S=1, Ms=1,-1) '
      write(6,*)'initialize1: Ms=1, -1 are mixed at the level of wavefunction    '
!
   else
      write(6,*)'initialize1: error! unknown iguess = ', iguess
      stop
   end if
   call flush(6)
100 continue
end if
!
if (grandom) then
!   dtmp1 = 1.d-2
!   do lni=1,nunk
!     do ni=1,nrho
!       do nj=1,nrho
!         rho(ni,nj,lni) = dcmplx(random_normal(), random_normal()) * dtmp1
!       end do
!     end do
!   end do
   call init_random_seed()
   rho = czero
   do nj=1,nrho
      do ni=1,nj-1
         rho(ni,nj,1) = dcmplx(random_normal(), random_normal())
         rho(nj,ni,1) = dconjg(rho(ni,nj,1))
      end do
      rho(nj,nj,1) = cunity / dble(nrho)
   end do
   write(6,*)
   write(6,*)'initialize1: rho is initialized randomly '
   call cmatout(nrho, nrho, rho(1,1,1), dmtmp1)
   call checkhermicity(rho(1,1,1), nrho, cmtmp1, nrho, lherm, dtmp1)
   if (.not. lherm) then
      write(6,*)'initialize1: error! rho is not hermitian ', dtmp1
      stop
   end if
   open(unit=68, file='rhoini.data', form='unformatted')
   write(68)nrho**2
   write(68)((rho(ni,nj,1), ni=1,nrho), nj=1,nrho)
   close(68)
   write(6,*)'initialize1: initial rho written to <rhoini.data> '
   call flush(6)
end if
!
rewind(5)
lreadini = .false.
lreadmat = .false.
lreadfile = .false.
read(5, readini, end=101)
101 continue
if (lreadini) then
   if (lreadfile) then
      open(unit=68, file='rhoini.data', form='unformatted')
      read(68)ni
      if (ni .ne. nrho**2) then
         write(6,*)'initialize1: error read from <rhoini.data> '
         write(6,*)ni, nrho**2
         stop
      end if
      read(68)((rho(ni,nj,1), ni=1,nrho), nj=1,nrho)
      close(68)
      write(6,*)'initialize1: initial rho read from <rhoini.data> '
      call cmatout(nrho, nrho, rho(1,1,1), dmtmp1)
      call checkhermicity(rho(1,1,1), nrho, cmtmp1, nrho, lherm, dtmp1)
      if (.not. lherm) then
         write(6,*)'initialize1: error! rho is not hermitian ', dtmp1
         stop
      end if
   else if (lreadmat) then
      read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,nrho), nj=1,nrho)  ! real part in order of (1,1), (2,1), (1,2), (2,2)  
      if (istat .ne. 0) then
         write(6,*)'initialize1: error reading real part of rho '
         stop
      end if
      read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,nrho), nj=1,nrho)  ! imag part in order of (1,1), (2,1), (1,2), (2,2) 
      if (istat .ne. 0) then
         write(6,*)'initialize1: error reading imag part of rho '
         stop
      end if
      do nj=1,nrho
         do ni=1,nrho
            rho(ni,nj,1) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
         end do
      end do
      write(6,*)'initialize1: rho read from input file in matrix form '
      call cmatout(nrho, nrho, rho(1,1,1), dmtmp1)
      call flush(6)
      call checkhermicity(rho(1,1,1), nrho, cmtmp1, nrho, lherm, dtmp1)
      if (.not. lherm) then
         write(6,*)'initialize1: error! rho is not hermitian ', dtmp1
         stop
      end if
   else
      read(5,*,iostat=istat) (rhodiag(ni), ni=1,nrho)
      write(6,*)'initialize1: read rhodiag from input file', istat
      write(6,*)'initialize1: ', (rhodiag(ni), ni=1,nrho)
      call flush(6)
      rho(1:nrho,1:nrho,1:nunk) = czero
      do ni=1,nrho
         rho(ni,ni,1) = dcmplx(rhodiag(ni), 0.d0)
      end do
   end if
endif
!
! hermitianize rho(*,*,1)
!
do ni=1,nrho
  do nj=1,ni
    rho(ni,nj,1) = 5.d-1 * (rho(ni,nj,1) + dconjg(rho(nj,ni,1)))
    rho(nj,ni,1) = dconjg( rho(ni,nj,1) )
  end do
end do
!
! Initialize system Hamiltonian hs as zero
!
hs(1:nrho,1:nrho) = czero
!
memory = memory + dble(nunk * 16 * nrho**2) / dble(1024**2)
!
! Heat bath
! 
if (lhb) then
   ltmp1 = nunk * nrho**2 * nbath_hb * 2  ! both real and imaginary parts
   dtmp1  = dble(ltmp1 * 8) / dble(1024**2)
   memohb = 0.d0
!
   allocate(rho_hb(nrho,nrho,nunk,nbath_hb), STAT=istat)
   rho_hb(1:nrho,1:nrho,1:nunk,1:nbath_hb) = czero
   memohb = memohb + dtmp1 
   if (istat .ne. 0) then
      write(6,*)'initialize1: error! memory allocation fail 2 '
      write(6,*)'memory desired ', memohb, ' MB'
      stop
   end if
end if
!
return
!
end subroutine initialize1
