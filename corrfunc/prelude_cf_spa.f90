subroutine prelude_cf_spa
use matmod
use tmpmatmod
use sparsemod
use corrfuncmod
use random
implicit none
include '../include/sizes'
include '../include/common'
!
! Note : The sparsity topology of Hamiltonian determines the sparsity of reduced density
!        matrix and all ADOs. 
!        With nonzero off-diagonal bath correlation function, H is replaced by some Heff,
!        which accounts for the off-diagonal hybridization function
!
integer                 :: isgnb, ispinb, iorbsb
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs, malf
integer                 :: na, nb, nc, mi, mj, nnz
integer*1               :: fact, itype
integer                 :: mindex(MAXTIER)
integer*8               :: lni, lnj, lnk, lnm, lnl, lunk_last, lstart
integer                 :: isgn, iorbs, ispin, ialf, iorbs2, nrho2, isame, icomp, icor
integer                 :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
integer                 :: ntmp10, ntmp11
logical                 :: lherm, lsame, lcomp, lexist, lexist1, lexist2, ltmp1, ltmp2, ltmp3, ltmp4
character*1             :: sidea, transa, transb
real*8                  :: dtmp1, dtmp2
real*8                  :: dsmall_spa
complex*16              :: ctmp1
real*8                  :: t_on
complex*16, allocatable :: work(:), cexphs(:,:), cexpheff(:,:), val_tmp(:)
real*8,     allocatable :: tmpval(:), rwork(:)
complex*16, allocatable :: hybrid(:,:,:)
complex*16              :: cmix1, cmix2
integer,    allocatable :: rowcoo(:), colcoo(:)
complex*16, allocatable :: valcoo(:)
character               :: matdescra(6)
integer,    allocatable :: nvec1(:), nvec2(:)
complex*16, allocatable :: cvec1(:)
complex*16, allocatable :: cmat1(:,:), cmat2(:,:)
!
namelist / coupling / readcp, readmat, lrcmplx
!
iorbsb = iorbs_dos
ispinb = ispin_dos
write(6,*)
write(6,*)'prelude_cf_spa: entering subroutine  '
write(6,*)'prelude_cf_spa: iorbsb, ispinb  ', iorbsb, ispinb
call flush(6)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'prelude_cf_spa: error entrance, lsparse = ', lsparse
   stop
end if
!
lexist = .false.
inquire(file='sparse_index_cf.data', exist=lexist1, err=999)   ! unit = 54
999 continue
inquire(file='sparse_info_cf.data', exist=lexist2, err=998)    ! unit = 55 
998 continue
if (lexist1 .and. lexist2) then
   lexist = .true.
end if
!
if (lexist) then
   write(6,*)
   write(6,*)'prelude_cf_spa: <sparse_index_cf.data> and <sparse_info_cf.data> found '
   write(6,*)'prelude_cf_spa: start reading '
   open(unit=54, file='sparse_index_cf.data', form='binary', status='old')
   rewind(54)
   read(54)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9, ntmp10, ntmp11
   if (ntmp1 .ne. ntier .or. ntmp2 .ne. ncor .or. ntmp3 .ne. norbs .or.     &
       ntmp4 .ne. nspin .or. ntmp5 .ne. nalf .or. ntmp6 .ne. numfff .or.    &
       ntmp7 .ne. ntier0 .or. ntmp8 .ne. ndrawer_slow .or. ntmp9 .ne. ncor_slow .or. &
       ntmp10 .ne. iorbs_dos .or. ntmp11 .ne. ispin_dos) then
       write(6,*)
       write(6,*)'prelude_cf_spa: <sparse_index_cf.data> incompatible with present job '
       write(6,*)'prelude_cf_spa: abort reading '
       write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9, ntmp10, ntmp11
       write(6,*)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow, iorbs_dos, ispin_dos
       call flush(6)
       close(54)
       lexist = .false.
       goto 15
   end if
   read(54)lni, ltmp1, ltmp2, ltmp3, ltmp4
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen) then
      write(6,*)'prelude_cf_spa: <sparse_index_cf.data> incompatible with present job '
      write(6,*)'prelude_cf_spa: abort reading '
      write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4
      write(6,*)nunk, offcor, lequileads, lsimple, lscreen
      call flush(6)
      close(54)
      lexist = .false.
      goto 15
   end if
!
   open(unit=55, file='sparse_info_cf.data', form='binary', status='old')
   rewind(55)
   read(55)lni, ltmp1, ltmp2, ltmp3, ltmp4
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen) then
       write(6,*)'prelude_cf_spa: <sparse_info_cf.data> incompatible with present job '
       write(6,*)'prelude_cf_spa: abort reading '
       write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4
       write(6,*)nunk, offcor, lequileads, lsimple, lscreen
       call flush(6)
       close(55)
       lexist = .false.
       goto 15
   end if
end if
!
15 continue
!
if (lproj) then
    dsmall_spa = dpico
else
    dsmall_spa = dsmall
end if
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), STAT=istat)
!
call refresh_rhosys_spa
cmat1(1:nrho,1:nrho) = rhosys(1:nrho,1:nrho)
call calchs(cmat1)
!
! cexphs = exp(-H)
!
allocate(cexphs(nrho,nrho), cexpheff(nrho,nrho), STAT=istat)
!
rewind(5)
do ni=1,5
   read(5,*)nj
end do
read(5,*)malf
!
allocate(hybrid(norbs,norbs,malf), STAT=istat)
!
cmat2(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
!
! --- consider possible change of system Hamiltonian during time evolution 
! --- assume Hamiltonian at time = t_on captures the sparsity topology
!if (igroundsteady .eq. 1 .and. .not. fixdot) then
if (igroundsteady .ne. 0 .and. .not. fixdot .and. lspdhs) then
    t_on = 1.d10   ! consider the asymptotic limit
    call calcdhs(t_on, cmat1)
    cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
    write(6,*)'prelude_cf_spa: consider the sparsity of dhs at time ', t_on
    !
    t_on = 1.d0    ! also consider a transient time (to be consistent with <prelude_td_spa.f90>)
    call calcdhs(t_on, cmat1)
    cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
    write(6,*)'prelude_cf_spa: consider the sparsity of dhs at time ', t_on
    !
    call flush(6)
end if
cmtmp2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho)
!
call checkhermicity(cmtmp2, nrho, cmtmp1, nrho, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)
   write(6,*)'prelude_cf_spa: error, H is not hermitian '
   stop
end if
!
if (lproj) then
    write(6,*)'prelude_cf_spa: projected system Hamiltonian '
else
    write(6,*)'prelude_cf_spa: system Hamiltonian '
end if
cmtmp1 = cmtmp2 * hbar
call cmatout(nrho, nrho, cmtmp1, dmtmp1)
call flush(6)
if (lproj) then
    !call calc_exp_hmat_diag(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexphs, 1, nprojf-nproji+1)
    call calc_exp_hmat_proj(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexphs, nstate_proj, istate_proj)
    cmtmp3(1:nrho,1:nrho) = cexphs(1:nrho,1:nrho)
else
    call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexphs)
    call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 11, cmtmp3)
    call checksparsity(cexphs, cmtmp3, nrho, nrho, dsmall, lsame, ntmp1)
    if (.not. lsame) then
       write(6,*)'prelude_cf_spa: error! sparsity analysis failed 1 ', ntmp1
       stop
    end if
end if
!call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexphs)
!
!lwork = nrho**2
!allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
!call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
!if (info .ne. 0) then
!   write(6,*)
!   write(6,*)'prelude_cf_spa: error! diagnalization failed 1 ', info
!   stop
!end if
!! comfine the absolute value of eigenvalues
!do ni=1,nrho
!   if ( tmpval(ni) .gt.  5.d0 ) tmpval(ni) =  5.d0
!   if ( tmpval(ni) .lt. -5.d0 ) tmpval(ni) = -5.d0
!end do
!cmtmp3(1:nrho,1:nrho) = czero
!do ni=1,nrho
!   cmtmp3(ni,ni) = dcmplx( dexp(-tmpval(ni)), 0.d0 )
!end do
!deallocate(work, rwork, tmpval, STAT=istat)
!call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
!           czero, cmtmp4, nrho)
!call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
!           czero, cmtmp3, nrho)
!cexphs(1:nrho,1:nrho) = cmtmp3(1:nrho,1:nrho)
!
hybrid(1:norbs,1:norbs,1:malf) = czero
do ialf=1,malf
   do ni=1,norbs
      hybrid(ni,ni,ialf) = dcmplx(linewidth(ialf,ni) / hbar, 0.d0)
   end do
end do
if (lspcor) then
    do ialf=1,malf
       do ispin=1,nspin
          do ni=1,norbs
             hybrid(ni,ni,ialf) = hybrid(ni,ni,ialf) + dcmplx(dlw_sp(ni,ispin,ialf) / hbar, 0.d0)
          end do
       end do
    end do
end if
if (lequileads) then
   do ialf=2,malf
      do ni=1,norbs
         hybrid(ni,ni,1) = hybrid(ni,ni,1) + hybrid(ni,ni,ialf)
      end do
   end do
end if
!
if (offcor) then
   cmtmp2(1:norbs,1:norbs) = czero
   if (megaflux) then
      hybrid(1:norbs,1:norbs,1:malf) = czero
      dtmp1 = 0.d0
      do ialf=1,malf
         hybrid(1:2,1:2,ialf) = dcmplx(linewidth(ialf,1) / hbar, 0.d0)
         dtmp1 = dtmp1 + linewidth(ialf,1)
      end do
      cmtmp2(1:norbs,1:norbs) = dcmplx(dtmp1, 0.d0)
   else
      hybrid(1:norbs,1:norbs,1:malf) = czero
!
      readcp  = .false. 
      readmat = .false. 
      lrcmplx = .false. 
      rewind(5)
      read(5, coupling, end=200)     
      200 continue
      if (readcp) then 
         do ialf=1,malf
            dmtmp3(1:norbs,1:norbs) = 0.d0
            if (readmat) then 
               do ni=1,norbs
                  if (lrcmplx) then
                     read(5,*,iostat=istat) (dmtmp2(ni,nj), nj=1,norbs), (dmtmp3(ni,nj), nj=1,norbs)
                  else
                     read(5,*,iostat=istat) (dmtmp2(ni,nj), nj=1,norbs)
                  end if
                  if (istat .ne. 0) then 
                     write(6,*)'prelude_cf_spa: error reading linewidth, ialf = ', ialf, ' iorbs = ', ni
                     stop
                  end if
               end do
            else 
               if (lrcmplx) then
                  read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs), ((dmtmp3(ni,nj), ni=1,norbs), nj=1,norbs)
               else
                  read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs)
               end if
               if (istat .ne. 0) then
                  write(6,*)'prelude_cf_spa: error reading linewidth matrix for ialf = ', ialf
                  stop
               end if
            end if
            hybrid(1:norbs,1:norbs,ialf) = dcmplx(dmtmp2(1:norbs,1:norbs), dmtmp3(1:norbs,1:norbs)) / hbar
            cmtmp2(1:norbs,1:norbs) = cmtmp2(1:norbs,1:norbs) + dcmplx(dmtmp2(1:norbs,1:norbs), dmtmp3(1:norbs,1:norbs)) / hbar
         end do  
      else
! hybrid function remains diagonal if readcp is turned off
         do ialf=1,malf
            do ni=1,norbs
               hybrid(ni,ni,ialf) = dcmplx(linewidth(ialf,ni) / hbar, 0.d0)
            end do
         end do
      end if
   end if
!
   if (lequileads) then
      do ialf=2,malf
         do ni=1,norbs
            hybrid(ni,ni,1) = hybrid(ni,ni,1) + hybrid(ni,ni,ialf)
         end do
      end do
   end if
!
! extract off-diagonal coupling
   cmtmp1(1:norbs,1:norbs) = czero
   do ni=1,norbs
      do nj=1,ni-1
         if (cdabs(cmtmp2(ni,nj)) .gt. dsmall_spa) then
            cmtmp1(ni,nj) = cmtmp2(ni,nj) * eye    ! (or -eye?) sign chosen arbitrarily
            cmtmp1(nj,ni) = dconjg(cmtmp1(ni,nj))  ! impose Hermicity
         end if
      end do
   end do  
! construct effective Hamiltonian corresponding to off-diagonal coupling
   !cmtmp2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho)
   cmtmp2(1:nrho,1:nrho) = czero
   do ni=1,norbs
      do nj=1,norbs
         if (cdabs(cmtmp1(ni,nj)) .gt. dsmall_spa) then 
            do ispin=1,nspin
               cmtmp3(1:nrho,1:nrho) = amsori(1:nrho,1:nrho,nj,ispin) * cunity
               call zamsmm1('L', 'C', 'N', ni, ispin, cmtmp1(ni,nj), cmtmp3, nrho, cunity, cmtmp2, nrho)
            end do
         end if
      end do
   end do
   if (lproj) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp2, nrho, czero, cmtmp1, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
   end if
   cmtmp2(1:nrho,1:nrho) = cmtmp2(1:nrho,1:nrho) + cmat2(1:nrho,1:nrho)
   call checkhermicity(cmtmp2, nrho, cmtmp3, nrho, lherm, dtmp1)
   if (.not. lherm) then
      write(6,*)
      write(6,*)'prelude_cf_spa: error, effective H is not hermitian '
      stop
   end if
!
   if (lproj) then
       !call calc_exp_hmat_diag(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexpheff, 1, nprojf-nproji+1)
       call calc_exp_hmat_proj(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexpheff, nstate_proj, istate_proj)
       cmtmp3(1:nrho,1:nrho) = cexpheff(1:nrho,1:nrho)
   else
       call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexpheff)
       call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 11, cmtmp3)
       call checksparsity(cexpheff, cmtmp3, nrho, nrho, dsmall, lsame, ntmp1)
       if (.not. lsame) then
          write(6,*)'prelude_cf_spa: error! sparsity analysis failed 2 ', ntmp1
          stop
       end if
   end if
   !call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexpheff)
!
!! diagonalize Heff
!   lwork = nrho**2
!   allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
!   call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
!   if (info .ne. 0) then
!      write(6,*)
!      write(6,*)'prelude_cf_spa: error! diagnalization failed 2 ', info
!      stop
!   end if
!! comfine the absolute value of eigenvalues
!   do ni=1,nrho
!      if ( tmpval(ni) .gt.  5.d0 ) tmpval(ni) =  5.d0
!      if ( tmpval(ni) .lt. -5.d0 ) tmpval(ni) = -5.d0
!   end do
!   cmtmp3(1:nrho,1:nrho) = czero
!   do ni=1,nrho
!      cmtmp3(ni,ni) = dcmplx( dexp(-tmpval(ni)), 0.d0 )
!   end do
!   deallocate(work, rwork, tmpval, STAT=istat)
!   call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
!              czero, cmtmp4, nrho)
!   call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
!              czero, cmtmp3, nrho)
!   cexpheff(1:nrho,1:nrho) = cmtmp3(1:nrho,1:nrho)
!
   cmtmp1(1:nrho,1:nrho) = cexphs(1:nrho,1:nrho) - cexpheff(1:nrho,1:nrho)
   call cmaxmat(nrho, nrho, cmtmp1, nrho, dtmp1)
   write(6,*)
   if (dabs(dtmp1) .lt. dsmall_spa) then
      write(6,*)'prelude_cf_spa: Heff is same as H '
   else
      write(6,*)'prelude_cf_spa: Heff is different from H '
   end if
end if
!
! zero-th tier : reduced density matrix
!
if (.not. offcor) then
   cmtmp3(1:nrho,1:nrho) = cexphs(1:nrho,1:nrho)
else 
   cmtmp3(1:nrho,1:nrho) = cexpheff(1:nrho,1:nrho)
end if
!
nrho2 = nrho**2
allocate(irowcf_hsys(nrho2), icolcf_hsys(nrho2), cvalcf_hsys(nrho2), STAT=istat)
cvalcf_hsys = czero
call zmat_dns2coo(nnzcf_hsys, irowcf_hsys, icolcf_hsys, cvalcf_hsys, cmtmp3, nrho, nrho, nrho, dsmall_spa)
if (nnzcf_hsys .le. 0 .or. nnzcf_hsys .gt. nrho2) then
   write(6,*)'prelude_cf_spa: error nnzcf_hsys ', nnzcf_hsys
   stop
end if
write(6,*)'prelude_cf_spa: nnzcf_hsys = ', nnzcf_hsys
call flush(6)
!
if (lexist) then
   read(54)ntmp1
   if (ntmp1 .le. 0 .or. ntmp1 .gt. nrho2) then
      write(6,*)'prelude_cf_spa: error nnzcf_hsys read from <sparse_index_cf.data> ', ntmp1
      write(6,*)'prelude_cf_spa: abort reading '
      close(54)
      close(55)
      lexist = .false.
   else
      allocate(nvec1(ntmp1), nvec2(ntmp1), cvec1(ntmp1), STAT=istat)
      do ni=1,ntmp1
         read(54)nvec1(ni), nvec2(ni), cvec1(ni)
      end do
      call zmat_coo2dns(ntmp1, nvec1, nvec2, cvec1, cmtmp1, nrho, nrho, nrho)
      call zmat_comp_sparsity(cmtmp1, nrho, nrho, nrho, cmtmp3, nrho, nrho, nrho, &
                              isame, icomp, na, nb, dsmall_spa) 
      deallocate(nvec1, nvec2, cvec1, STAT=istat)
      if (isame .ne. 1) then
         write(6,*)'prelude_cf_spa: incompatible <sparse_index_cf.data> '
         write(6,*)'prelude_cf_spa: abort reading '
         close(54)
         close(55)
         lexist = .false.
      end if 
   end if  
end if
!
if (lexist) then
   read(55)lunkcf_spa
   allocate(nnzcf_spa(nunk), indcf_spa(nunk), STAT=istat)
   allocate(irowcf_spa(lunkcf_spa), icolcf_spa(lunkcf_spa),          &
            brhoh_spa(lunkcf_spa), brhoa_spa(lunkcf_spa), STAT=istat)
   do lni=1,nunk
      read(54)nnzcf_spa(lni), indcf_spa(lni)
   end do
   close(unit=54)
   do lni=1,lunkcf_spa
      read(55)irowcf_spa(lni), icolcf_spa(lni)
   end do
   close(unit=55)
   write(6,*)'prelude_cf_spa: <sparse_index_cf.data> and <sparse_info_cf.data> read successfully'
   call flush(6)
   goto 101
end if
!
lunkcf_spa = 0
allocate(indcf_spa(nunk), nnzcf_spa(nunk), STAT=istat)
!
! exp(-Heff) * B * exp(-Heff)
sidea  = 'L'
transa = 'N'
transb = 'N'
if (lproj) then
    call zamsmm1_proj(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp1, nrho)
else
    call zamsmm1(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp1, nrho)
end if
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp1, nrho, czero, cmtmp7, nrho)
!
! exp(-Heff) * (B+) * exp(-Heff)
sidea  = 'R'
transa = 'C'
transb = 'N'
if (lproj) then
    call zamsmm1_proj(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp2, nrho)
else
    call zamsmm1(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp2, nrho)
end if
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp3, nrho, czero, cmtmp8, nrho)
do 
   cmix1 = dcmplx(random_normal(), random_normal())
   cmix2 = dcmplx(random_normal(), random_normal())
   if ( cdabs(cmix1) .gt. micron .and. cdabs(cmix2) .gt. micron ) exit
end do
!
! term_Heff =  exp(-Heff) * (B+ + B) * exp(-Heff)
cmtmp1(1:nrho,1:nrho) = cmix1 * cmtmp7(1:nrho,1:nrho) + cmix2 * cmtmp8(1:nrho,1:nrho)
!
nk = 0
do nj=1,nrho ! column first
   do ni=1,nrho
!      if (nj .eq. 1 .and. ni .eq. 1) then   ! enforce the first element to be recorded
!         nk = nk + 1
!         cycle
!      end if
      if ( cdabs(cmtmp1(ni,nj)) .gt. dsmall_spa ) then
         nk = nk + 1
      end if
   end do
end do
lunkcf_spa   = nk
indcf_spa(1) = 1
nnzcf_spa(1) = nk
allocate(irowcf_spa(lunkcf_spa), icolcf_spa(lunkcf_spa), STAT=istat)
na = 0
spa_rho: do nj=1,nrho
   do ni=1,nrho
      if (na .ge. nk) exit spa_rho
      if ( cdabs(cmtmp1(ni,nj)) .gt. dsmall_spa) then
      !if ( cdabs(cmtmp1(ni,nj)) .gt. dsmall .or. &
      !     (nj .eq. 1 .and. ni .eq. 1) ) then        ! (1,1) element must be included
         na = na + 1
         irowcf_spa(na) = ni
         icolcf_spa(na) = nj
      end if
   end do
end do spa_rho
!
open(unit=48, file='rhocoo.tmp', form='binary', status='unknown', access='sequential')
rewind(48)
do lni=1,lunkcf_spa
   write(48)irowcf_spa(lni), icolcf_spa(lni)
end do
open(unit=49, file='rhoval.tmp', form='binary', status='unknown', access='sequential')
rewind(49)
do lni=1,lunkcf_spa
   write(49)cmtmp1(irowcf_spa(lni),icolcf_spa(lni))
end do
!
lunk_last = lunkcf_spa
deallocate(irowcf_spa, icolcf_spa, STAT=istat)
!
! from second to terminal tier
!
nk = nrho**2
allocate(rowcoo(nk), colcoo(nk), valcoo(nk), STAT=istat)
!
do itier=2,ntier0
   allocate(irowcf_spa(lunk_last), icolcf_spa(lunk_last), val_tmp(lunk_last), STAT=istat)
   rewind(48)
   do lni=1,lunk_last
      read(48)irowcf_spa(lni), icolcf_spa(lni)
   end do
   rewind(49)
   do lni=1,lunk_last
      read(49)val_tmp(lni)
   end do
!
   do lni=nfirst(itier),nlast(itier)
      lnl = index_coef_ref(lni)
      zmtmp2(1:nrho,1:nrho) = czero
!
      lnm = ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1)
      mindex(1:itier-1) = indextable(lnm+1:lnm+itier-1)
!
      dtmp2 = dlarge
      do ni=1,itier-1
         nk = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
         if (psdfff .and. mmfff(nk) .gt. 1) then
             call findlowertier(itier, mindex, ni, lnk, itype, fact)
             if (itier .eq. 2) then
               lnj = 1
             else
               lnj = (lnk + 1 - ifirst(itier-1)) / (itier - 2) + nfirst(itier-1)
             end if
         else
             lnj    = indexcoef(lnl)
             itype  = itypecoef(lnl)
         end if
         isgn   = mpm(nk)
         iorbs  = morbs(nk)
         ispin  = mspin(nk)
         ialf   = ilead(nk)
         icor   = mcor(nk)
         lnl    = lnl + 1
         lstart = indcf_spa(lnj)
         nnz    = nnzcf_spa(lnj)
!
         sidea  = 'L'
         transa = 'N'
         transb = 'N'
         if (isgn .eq. 1) then
            transa = 'C'
         end if
         if (itype .ne. 0) then
            transb = 'C'
         end if   
!
         if (.not. offcor) then
            dtmp2 = min(dtmp2, cdabs(hybrid(iorbs,iorbs,ialf)))
            if ( cdabs(hybrid(iorbs,iorbs,ialf)) .lt. dsmall_spa ) cycle
! cmtmp4 : ams*rho
            sidea  = 'L'
            if (lprams) then
                call zamsmm_coo_proj(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                     icolcf_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
            else
                call zamsmm_coo(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                icolcf_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
            end if
! cmtmp5 : rho*ams
            sidea  = 'R'
            if (lprams) then
                call zamsmm_coo_proj(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                     icolcf_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
            else
                call zamsmm_coo(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                icolcf_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
            end if
!
! cmtmp1 : exp(-Heff) * (ams*rho)  
! cmtmp2 : exp(-Heff) * (rho*ams)
            call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, czero, cmtmp1, nrho)
            call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp5, nrho, czero, cmtmp2, nrho)
! avoid accidental cancellation between rho*ams and ams*rho
            do 
               cmix1 = dcmplx(random_normal(), random_normal())
               cmix2 = dcmplx(random_normal(), random_normal())
               if ( cdabs(cmix1) .gt. micron .and. cdabs(cmix2) .gt. micron ) exit
            end do
            zmtmp2(1:nrho,1:nrho) = zmtmp2(1:nrho,1:nrho) + cmix1 * cmtmp1(1:nrho,1:nrho) + cmix2 * &
                                    cmtmp2(1:nrho,1:nrho)
!
         else 
            dtmp1 = 0.d0 
            do iorbs2=1,norbs
               dtmp1 = dtmp1 + cdabs(hybrid(iorbs,iorbs2,ialf))
               if ( cdabs(hybrid(iorbs,iorbs2,ialf)) .lt. dsmall_spa ) cycle
               sidea  = 'L'
               if (lprams) then
                   call zamsmm_coo_proj(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                        icolcf_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
               else
                   call zamsmm_coo(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                   icolcf_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
               end if
               sidea  = 'R'
               if (lprams) then
                   call zamsmm_coo_proj(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                        icolcf_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
               else
                   call zamsmm_coo(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irowcf_spa(lstart), &
                                   icolcf_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
               end if
               call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, czero, cmtmp1, nrho)
               call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp5, nrho, czero, cmtmp2, nrho)
               do
                  cmix1 = dcmplx(random_normal(), random_normal())
                  cmix2 = dcmplx(random_normal(), random_normal())
                  if ( cdabs(cmix1) .gt. micron .and. cdabs(cmix2) .gt. micron ) exit
               end do
               zmtmp2(1:nrho,1:nrho) = zmtmp2(1:nrho,1:nrho) + cmix1 * cmtmp1(1:nrho,1:nrho) + cmix2 * &
                                       cmtmp2(1:nrho,1:nrho)
            end do ! end of loop over iorbs2
            dtmp2 = min(dtmp2, dtmp1)
         end if
      end do ! end of loop over ni
!
      if (dtmp2 .gt. dsmall_spa) then
          call zmat_dns2coo(na, rowcoo, colcoo, valcoo, zmtmp2, nrho, nrho, nrho, dsmall_spa)
      else
          na = 0
      end if
      lunkcf_spa = lunkcf_spa + na
      nnzcf_spa(lni) = na
      indcf_spa(lni) = indcf_spa(lni-1) + nnzcf_spa(lni-1)
      if (na .gt. 0) then
         do ni=1,na
            write(48)rowcoo(ni), colcoo(ni)
            write(49)valcoo(ni)
         end do
      end if
   end do 
!
   lunk_last = lunkcf_spa
   deallocate(irowcf_spa, icolcf_spa, val_tmp, STAT=istat)
end do
deallocate(rowcoo, colcoo, valcoo, STAT=istat)
!
allocate(irowcf_spa(lunk_last), icolcf_spa(lunk_last),           &
         brhoh_spa(lunk_last), brhoa_spa(lunk_last), STAT=istat)
rewind(48)
do lni=1,lunk_last
   read(48)irowcf_spa(lni), icolcf_spa(lni)
end do
close(unit=48, status='delete')
close(unit=49, status='delete')
!
do lni=1,lunk_last
   brhoh_spa(lni) = dcmplx( random_normal(), random_normal() )
   brhoa_spa(lni) = dcmplx( random_normal(), random_normal() )
end do
!
101 continue
!
deallocate(cexphs, cexpheff, STAT=istat)
deallocate(hybrid, STAT=istat)
deallocate(cmat1, cmat2, STAT=istat)
!
! account for memory of used arrays 
! including: indcf_spa, nnzcf_spa, irowcf_spa, icolcf_spa, brhoh_spa, brhoa_spa
!
dtmp1  = dble(nunk * (4 + 8) + lunkcf_spa * 2 * 4 + lunkcf_spa * 2 * 16) / dble(1024**2)
memory = memory + dtmp1
write(6,*)
write(6,*)'prelude_cf_spa: memory allocated for sparsity info ', dtmp1, ' MB '
call flush(6)
!
lni = nrho**2 * nunk 
write(6,*)
write(6,*)'prelude_cf_spa: total unknown variables = ', lni
write(6,*)'prelude_cf_spa: nonzero variables       = ', lunkcf_spa
write(6,*)'prelude_cf_spa: sparsity ratio for rho  = ', dble(lunkcf_spa) / dble(lni)
call flush(6)
!
! output files 
! 
if (.not. lexist) then
   open(unit=54, file='sparse_index_cf.data', form='binary')
   rewind(54)
   write(54)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow, iorbs_dos, ispin_dos
   write(54)nunk, offcor, lequileads, lsimple, lscreen
   write(54)nnzcf_hsys
   do ni=1,nnzcf_hsys
      write(54)irowcf_hsys(ni), icolcf_hsys(ni), cvalcf_hsys(ni)
   end do
   do lni=1,nunk
      write(54)nnzcf_spa(lni), indcf_spa(lni)
   end do
   close(54)
!
   open(unit=55, file='sparse_info_cf.data', form='binary')
   rewind(55)
   write(55)nunk, offcor, lequileads, lsimple, lscreen
   write(55)lunkcf_spa
   do lni=1,lunkcf_spa
      write(55)irowcf_spa(lni), icolcf_spa(lni)
   end do
   close(55)
end if
!
! Debug USE : compare brhoh and brhoh_spa, and brhoa and brhoa_spa, 
!             to verify the locations of nonzero elements
!
!lsame = .true.
!lcomp = .true.
!check1: do itier=2,ntier0
!   do lni=nfirst(itier),nlast(itier)
!      cmtmp1(1:nrho,1:nrho) = brhoh(1:nrho,1:nrho,lni) 
!      lstart = indcf_spa(lni)
!      call zmat_coo2dns(nnzcf_spa(lni), irowcf_spa(lstart), icolcf_spa(lstart), brhoh_spa(lstart), &
!                        cmtmp2, nrho, nrho, nrho)
!      call zmat_comp_sparsity(cmtmp1, nrho, nrho, nrho, cmtmp2, nrho, nrho, nrho, ni, nj, na, nb, dsmall_spa)
!      if (ni .eq. 0) lsame = .false.
!      if (nj .eq. 0) lcomp = .false.
!      if (.not. lcomp) exit check1
!   end do
!end do check1
!write(6,*)
!if (lsame) then
!   write(6,*)'prelude_cf_spa: brhoh and brhoh_spa have same sparsity topology '
!else
!   if (lcomp) then
!      write(6,*)'prelude_cf_spa: sparsity of brhoh is covered by brhoh_spa '
!   else
!      write(6,*)'prelude_cf_spa: attention! sparsity of brhoh differs from brhoh_spa '
!   end if
!end if
!call flush(6)
!!
!lsame = .true.
!lcomp = .true.
!check2: do itier=2,ntier0
!   do lni=nfirst(itier),nlast(itier)
!      cmtmp1(1:nrho,1:nrho) = brhoa(1:nrho,1:nrho,lni) 
!      lstart = indcf_spa(lni)
!      call zmat_coo2dns(nnzcf_spa(lni), irowcf_spa(lstart), icolcf_spa(lstart), brhoa_spa(lstart), &
!                        cmtmp2, nrho, nrho, nrho)
!      call zmat_comp_sparsity(cmtmp1, nrho, nrho, nrho, cmtmp2, nrho, nrho, nrho, ni, nj, na, nb, dsmall_spa)
!      if (ni .eq. 0) lsame = .false.
!      if (nj .eq. 0) lcomp = .false.
!      if (.not. lcomp) exit check2
!   end do
!end do check2
!if (lsame) then
!   write(6,*)'prelude_cf_spa: brhoa and brhoa_spa have same sparsity topology '
!else
!   if (lcomp) then
!      write(6,*)'prelude_cf_spa: sparsity of brhoa is covered by brhoa_spa '
!   else
!      write(6,*)'prelude_cf_spa: attention! sparsity of brhoa differs from brhoa_spa '
!   end if
!end if
!call flush(6)
!
111 continue
!
write(6,*)'prelude_cf_spa: leaving subroutine '
call flush(6)
!
return
end subroutine prelude_cf_spa
