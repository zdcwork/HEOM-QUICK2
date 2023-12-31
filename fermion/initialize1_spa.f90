subroutine initialize1_spa
use matmod
use tmpmatmod
use random
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer        :: istat, nnz, ni, nj
integer        :: ncount
integer*8      :: lni
logical        :: lherm
real*8         :: dtmp1
complex*16     :: ctmp1
!
complex*16, allocatable :: cmat1(:,:)
!
namelist / readini / lreadini, lreadmat, lreadfile
!
write(6,*)
write(6,*)'initialize_spa1: nnz(rhosys) = ', nnz_spa(1)
write(6,*)'initialize_spa1: (irow, icol) of nonzero elements in rhosys '
lni = ind_spa(1)
do ni=1,nnz_spa(1)
   write(6,*)ni, irow_spa(lni-1+ni), icol_spa(lni-1+ni)
end do
call flush(6)
!
allocate(rsdm(norbs,norbs,nspin), STAT=istat)
!
rho_spa(1)          = cunity
rho_spa(2:lunk_spa) = czero
!
ncount = 0 
nnz = nnz_spa(1)
lni = ind_spa(1)
do ni=1,nnz
   if (irow_spa(lni-1+ni) .eq. icol_spa(lni-1+ni)) ncount = ncount + 1
end do
ctmp1 = cunity / dble(ncount)
write(6,*)
write(6,*)'initialize1_spa: nonzero diagonal elements in rhosys ', ncount
call flush(6)
!
do ni=1,nnz
!   if (irow_spa(lni-1+ni) .eq. icol_spa(lni-1+ni)) then
!      rho_spa(lni-1+ni) = ctmp1
!   end if
   if (irow_spa(lni-1+ni) .eq. nrho .and. &
       icol_spa(lni-1+ni) .eq. nrho) then
       rho_spa(lni-1+ni) = cunity
   else 
       rho_spa(lni-1+ni) = czero
   end if
end do
!
! special case for double dot (norbs=2 and nspin=2)
!
if (norbs .eq. 2 .and. nspin .eq. 2) then  
   if (iguess .eq. 0) goto 100
!
   rho_spa(1:lunk_spa) = czero
   write(6,*)
!
   if (iguess .eq. -1) then   ! empty state
      loop_a: do ni=1,nnz
         if (irow_spa(ni) .eq. 1 .and. icol_spa(ni) .eq. 1) then
             rho_spa(ni) = cunity
             exit loop_a
         end if
      end do loop_a
      write(6,*)'initialize1_spa: initial rho is set to empty state '
!
   else if (iguess .eq. 1) then ! S=0 (singlet state)
      ncount = 0
      loop_b: do ni=1,nnz
         if (irow_spa(ni) .eq. 7 .and. icol_spa(ni) .eq. 7) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 10 .and. icol_spa(ni) .eq. 10) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 7 .and. icol_spa(ni) .eq. 10) then
             rho_spa(ni) = dcmplx(-5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 10 .and. icol_spa(ni) .eq. 7) then
             rho_spa(ni) = dcmplx(-5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (ncount .eq. 4) exit loop_b
      end do loop_b
      if (ncount .ne. 4) then
         write(6,*)'initialize1_spa: Warning! iguess, ncount ', iguess, ncount
      else 
         write(6,*)'initialize1_spa: initial rho is set to singlet state (S=0) '
      end if
!
   else if (iguess .eq. 2) then ! S=1, Ms=0 (triplet state)
      ncount = 0
      loop_c: do ni=1,nnz
         if (irow_spa(ni) .eq. 7 .and. icol_spa(ni) .eq. 7) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 10 .and. icol_spa(ni) .eq. 10) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 7 .and. icol_spa(ni) .eq. 10) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 10 .and. icol_spa(ni) .eq. 7) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (ncount .eq. 4) exit loop_c
      end do loop_c
      if (ncount .ne. 4) then
         write(6,*)'initialize1_spa: Warning! iguess, ncount ', iguess, ncount
      else
         write(6,*)'initialize1_spa: initial rho is set to triplet state (S=1, Ms=0) '
      end if
!
   else if (iguess .eq. 3) then ! S=1, Ms=1 (triplet state)
      loop_d: do ni=1,nnz
         if (irow_spa(ni) .eq. 6 .and. icol_spa(ni) .eq. 6) then
             rho_spa(ni) = cunity
             exit loop_d
         end if
      end do loop_d
      write(6,*)'initialize1_spa: initial rho is set to triplet state (S=1, Ms=1) '
!
   else if (iguess .eq. 4) then ! S=1, Ms=-1 (triplet state)
      loop_e: do ni=1,nnz
         if (irow_spa(ni) .eq. 11 .and. icol_spa(ni) .eq. 11) then
             rho_spa(ni) = cunity
             exit loop_e
         end if
      end do loop_e
      write(6,*)'initialize1_spa: initial rho is set to triplet state (S=1, Ms=-1) '
!
   else if (iguess .eq. 5) then ! S=1, Ms=1,-1 (triplet state), 1:1 mixture of Ms=1 and Ms=-1
      loop_f: do ni=1,nnz
         if (irow_spa(ni) .eq. 6 .and. icol_spa(ni) .eq. 6) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
         end if
         if (irow_spa(ni) .eq. 11 .and. icol_spa(ni) .eq. 11) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
         end if
      end do loop_f
      write(6,*)'initialize1_spa: initial rho is set to triplet state (S=1, Ms=1,-1) '
      write(6,*)'initialize1_spa: Ms=1,-1 are mixed at the level of density matrix   '
!
   else if (iguess .eq. 6) then ! S=1, Ms=1,-1 (triplet state) mixing wavefunction
      ncount = 0
      loop_g: do ni=1,nnz
         if (irow_spa(ni) .eq. 6 .and. icol_spa(ni) .eq. 6) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 11 .and. icol_spa(ni) .eq. 11) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 6 .and. icol_spa(ni) .eq. 11) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (irow_spa(ni) .eq. 11 .and. icol_spa(ni) .eq. 6) then
             rho_spa(ni) = dcmplx(5.d-1, 0.d0)
             ncount = ncount + 1
         end if
         if (ncount .eq. 4) exit loop_g
      end do loop_g
      if (ncount .ne. 4) then
         write(6,*)'initialize1_spa: Warning! iguess, ncount ', iguess, ncount
      else
         write(6,*)'initialize1_spa: initial rho is set to triplet state (S=1, Ms=1,-1) '
         write(6,*)'initialize1_spa: Ms=1,-1 are mixed at the level of wavefunction     '
      end if
!
   else 
      write(6,*)'initialize1_spa: error! unknown iguess = ', iguess
      stop
   end if
   call flush(6)
100 continue
end if
!
if (grandom) then
   call init_random_seed()
   rho_spa(1:lunk_spa) = czero
   do lni=ind_spa(1),ind_spa(1)+nnz_spa(1),1
      if (irow_spa(lni) .eq. icol_spa(lni)) then
         rho_spa(lni) = cunity / dble(nrho)
      else
         rho_spa(lni) = dcmplx(random_normal(), random_normal()) 
      end if
   end do
end if
!

allocate(cmat1(nrho,nrho), STAT=istat)
allocate(rhodiag(nrho), STAT=istat)
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
         write(6,*)'initialize1_spa: error read from <rhoini.data> '
         write(6,*)ni, nrho**2
         stop
      end if
      read(68)((cmat1(ni,nj), ni=1,nrho), nj=1,nrho)
      close(68)
      write(6,*)'initialize1_spa: initial rho read from <rhoini.data> '
      call cmatout(nrho, nrho, cmat1, dmtmp1)
      call checkhermicity(cmat1, nrho, cmtmp1, nrho, lherm, dtmp1)
      if (.not. lherm) then
         write(6,*)'initialize1_spa: error! rho is not hermitian ', dtmp1
         stop
      end if
      do lni=ind_spa(1),ind_spa(1)+nnz_spa(1),1
         rho_spa(lni) = cmat1(irow_spa(lni), icol_spa(lni))
      end do
   else if (lreadmat) then
      read(5,*,iostat=istat) ((dmtmp1(ni,nj), ni=1,nrho), nj=1,nrho)  ! real part in order of (1,1), (2,1), (1,2), (2,2)  
      if (istat .ne. 0) then
         write(6,*)'initialize1_spa: error reading real part of rho '
         stop
      end if
      read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,nrho), nj=1,nrho)  ! imag part in order of (1,1), (2,1), (1,2), (2,2) 
      if (istat .ne. 0) then
         write(6,*)'initialize1_spa: error reading imag part of rho '
         stop
      end if
      do nj=1,nrho
         do ni=1,nrho
            cmat1(ni,nj) = dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
         end do
      end do
      write(6,*)'initialize1_spa: rho read from input file in matrix form '
      call cmatout(nrho, nrho, cmat1, dmtmp1)
      call flush(6)
      call checkhermicity(cmat1, nrho, cmtmp1, nrho, lherm, dtmp1)
      if (.not. lherm) then
         write(6,*)'initialize1_spa: error! rho is not hermitian ', dtmp1
         stop
      end if
      do lni=ind_spa(1),ind_spa(1)+nnz_spa(1),1
         rho_spa(lni) = cmat1(irow_spa(lni), icol_spa(lni))
      end do
   else
      read(5,*,iostat=istat) (rhodiag(ni), ni=1,nrho)
      write(6,*)'initialize1_spa: read rhodiag from input file', istat
      write(6,*)'initialize1_spa: ', (rhodiag(ni), ni=1,nrho)
      call flush(6)
      rho_spa = czero
      do lni=ind_spa(1),ind_spa(1)+nnz_spa(1),1
         if (irow_spa(lni) .eq. icol_spa(lni)) then
            rho_spa(lni) = dcmplx(rhodiag(irow_spa(lni)), 0.d0)
         end if
      end do
   end if
endif
!
deallocate(cmat1, STAT=istat)

!
lni = ind_spa(1)
nnz = nnz_spa(1)
call zmat_herm_coo('u', nnz, rho_spa(lni), irow_spa(lni), icol_spa(lni), cmtmp1, nrho, nrho)
!
hs(1:nrho,1:nrho) = czero
!
return
end subroutine initialize1_spa
