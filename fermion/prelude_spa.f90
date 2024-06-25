subroutine prelude_spa
use matmod
use tmpmatmod
use sparsemod
use random
implicit none
include '../include/sizes'
include '../include/common'
!
! Note : The sparsity topology of Hamiltonian determines the sparsity of reduced density
!        matrix and all ADOs. 
!    (1) With nonzero off-diagonal bath correlation function, H is replaced by some Heff,
!        which accounts for the off-diagonal hybridization function
!    (2) For time-evolution case, if delta_H invovles new nonzero elements which are zero
!        in initial stationary-state H, these nonzero elements should be included in Heff
!
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs, malf
integer                 :: na, nb, nc, mi, mj, nnz, isame, icomp
integer*8               :: lni, lnj, lnk, lnm, lnl, lunk_last, lstart
integer*1               :: itype, fact
integer                 :: isgn, iorbs, ispin, ialf, iorbs2, nrho2, icor
integer                 :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
integer                 :: mindex(MAXTIER)
logical                 :: lherm, lsame, lcomp, lexist, lexist1, lexist2 
logical                 :: ltmp1, ltmp2, ltmp3, ltmp4, ltmp5
character*1             :: sidea, transa, transb
real*8                  :: dtmp1, dtmp2
real*8                  :: dsmall_spa
complex*16              :: ctmp1
real*8                  :: t_on
real*8                  :: cpu1, cpu2
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
!write(6,*)
!write(6,*)'prelude_spa: entering subroutine '
!call flush(6)
!
if (.not. lsparse) then
   write(6,*)
   write(6,*)'prelude_spa: error entrance, lsparse = ', lsparse
   stop
end if
!
! IMPORTANT!!! 
! If the sparsity topology of dhs is different from that of hs during the simulation 
! of real-time dynamics, the relevant nonzero elements should be explicitly accounted
! for in the following Heff (to determine the sparsity of reduced density matrix and 
! all ADOs)
!
!allocate(hs(nrho,nrho), STAT=istat)
!allocate(dhs(nrho,nrho), STAT=istat) 
!allocate(rhosys(nrho,nrho), STAT=istat)
!
! check whether <sparse_index.data> and <sparse_info.data> exist
! 
lexist = .false.
inquire(file='sparse_index.data', exist=lexist1, err=999)   ! unit = 52
999 continue
inquire(file='sparse_info.data', exist=lexist2, err=998)    ! unit = 53
998 continue
if (lexist1 .and. lexist2) then
   lexist = .true.
end if
!
if (lexist) then
   write(6,*)
   write(6,*)'prelude_spa: <sparse_index.data> and <sparse_info.data> found '
   write(6,*)'prelude_spa: start reading '
   open(unit=52, file='sparse_index.data', form='binary', status='old')
   rewind(52)
   read(52)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
   if (ntmp1 .ne. ntier .or. ntmp2 .ne. ncor .or. ntmp3 .ne. norbs .or.      &
       ntmp4 .ne. nspin .or. ntmp5 .ne. nalf .or. ntmp6 .ne. numfff .or.     &
       ntmp7 .ne. ntier0 .or. ntmp8 .ne. ndrawer_slow .or. ntmp9 .ne. ncor_slow) then
      write(6,*)
      write(6,*)'prelude_spa: <sparse_index.data> incompatible -- 1 '
      write(6,*)'prelude_spa: abort reading '
      write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
      write(6,*)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow
      call flush(6)
      close(52)
      lexist = .false.
      goto 15
   end if
   read(52)lni, ltmp1, ltmp2, ltmp3, ltmp4, ltmp5
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen .or. ltmp5 .ne. lproj) then
      write(6,*)'prelude_spa: <sparse_index.data> incompatible -- 2 '
      write(6,*)'prelude_spa: abort reading '
      write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4, ltmp5
      write(6,*)nunk, offcor, lequileads, lsimple, lscreen, lproj
      call flush(6)
      close(52)
      lexist = .false.
      goto 15
   end if
   read(52)ltmp1, ntmp1, ntmp2, ntmp3
   if (lproj .and. (ntmp1 .ne. nproji .or. ntmp2 .ne. nprojf) .and. (.not. lprocc)  .or. &
       lproj .and. (ltmp1 .ne. lprocc)                                              .or. &
       lproj .and. lprocc .and. (ntmp3 .ne. nprocc)) then
      write(6,*)'prelude_spa: <sparse_index.data> incompatible -- 3 '
      write(6,*)'prelude_spa: abort reading '
      write(6,*)ltmp1, ntmp1, ntmp2, ntmp3
      write(6,*)lprocc, nproji, nprojf, nprocc
      call flush(6)
      close(52)
      lexist = .false.
      goto 15
   end if
!
   open(unit=53, file='sparse_info.data', form='binary', status='old')
   rewind(53)
   read(53)lni, ltmp1, ltmp2, ltmp3, ltmp4, ltmp5
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen .or. ltmp5 .ne. lproj) then
      write(6,*)'prelude_spa: <sparse_info.data> incompatible -- 1 '
      write(6,*)'prelude_spa: abort reading '
      write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4, ltmp5
      write(6,*)nunk, offcor, lequileads, lsimple, lscreen, lproj
      call flush(6)
      close(53)
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
cmat1(1:nrho,1:nrho) = czero
!
!write(6,*)'prelude_spa: before calchs'
!call flush(6)
!lprths = .true.
lprths = .false.
call calchs(cmat1)
!write(6,*)'prelude_spa: after calchs'
!call flush(6)
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
    write(6,*)'prelude_spa: consider the sparsity of dhs at time ', t_on
    !
    t_on = 1.d0    ! also consider a transient time (to be consistent with <prelude_td_spa.f90>)
    call calcdhs(t_on, cmat1)
    cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
    write(6,*)'prelude_spa: consider the sparsity of dhs at time ', t_on
    !
    call flush(6)
end if
!if (igroundsteady .ne. 0 .and. .not. fixdot) then
!   if (igroundsteady .eq. 1) then
!      t_on = 1.d10
!   else if (igroundsteady .eq. 2) then
!      t_on = 1.d0   ! choose a time instant when the change to Hamiltonian is effective
!   else
!      write(6,*)'prelude_spa: error! unknown igroundsteady ', igroundsteady
!      stop
!   end if
!   call calcdhs(t_on, cmat1)
!   cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
!end if
cmtmp2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho)
! --- end of delta_Hamiltonian consideration
!
call checkhermicity(cmtmp2, nrho, cmtmp1, nrho, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)
    write(6,*)'prelude_spa: error, H is not hermitian '
    stop
end if
!
!call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 5, cmtmp3)
!write(6,*)
!write(6,*)'prelude_spa: approximate exp(-H)'
!call cmatout(nrho, nrho, cmtmp3, dmtmp1)
!call flush(6)
!
!call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexphs)
!call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 11, cmtmp3)
!if (lprams) then
if (lproj) then
    write(6,*)'prelude_spa: projected system Hamiltonian '
else
    write(6,*)'prelude_spa: system Hamiltonian '
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
       write(6,*)'prelude_spa: error! sparsity analysis failed 1 ', ntmp1
       stop
    end if
end if
!
!lwork = nrho**2
!allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
!call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
!if (info .ne. 0) then
!   write(6,*)
!   write(6,*)'prelude_spa: error! diagnalization failed 1 ', info
!   stop
!end if
!! comfine the absolute value of eigenvalues
!do ni=1,nrho
!   if ( tmpval(ni) .gt.  5.d0 ) tmpval(ni) =  5.d0
!   if ( tmpval(ni) .lt. -5.d0 ) tmpval(ni) = -5.d0
!end do
!!
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
!write(6,*)
!write(6,*)'prelude_spa: exp(-H)'
!call cmatout(nrho, nrho, cexphs, dmtmp1)
!call flush(6)
!!
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
      cmtmp2(1:norbs,1:norbs) = dcmplx(dtmp1,0.d0)
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
                     write(6,*)'prelude_spa: error reading linewidth, ialf = ', ialf, ' iorbs = ', ni
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
                  write(6,*)'prelude_spa: error reading linewidth matrix for ialf = ', ialf
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
      write(6,*)'prelude_spa: error, effective H is not hermitian '
      stop
   end if
!
   !if (lprams) then
   if (lproj) then
       !call calc_exp_hmat_diag(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexpheff, 1, nprojf-nproji+1)
       call calc_exp_hmat_proj(nrho, cmtmp2, nrho, dbeta(1)*hbar, cexpheff, nstate_proj, istate_proj)
       cmtmp3(1:nrho,1:nrho) = cexpheff(1:nrho,1:nrho)
   else
       call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 10, cexpheff)
       call calc_exp_hmat(nrho, cmtmp2, nrho, 0.31415926d0, 11, cmtmp3)
       call checksparsity(cexpheff, cmtmp3, nrho, nrho, dsmall, lsame, ntmp1)
       if (.not. lsame) then
          write(6,*)'prelude_spa: error! sparsity analysis failed 2 ', ntmp1
          stop
       end if
   end if
!
!! diagonalize Heff
!   lwork = nrho**2
!   allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
!   call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
!   if (info .ne. 0) then
!      write(6,*)
!      write(6,*)'prelude_spa: error! diagnalization failed 2 ', info
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
      write(6,*)'prelude_spa: Heff is same as H '
   else
      write(6,*)'prelude_spa: Heff is different from H '
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
allocate(irow_hsys(nrho2), icol_hsys(nrho2), cval_hsys(nrho2), STAT=istat)
cval_hsys = czero
call zmat_dns2coo(nnz_hsys, irow_hsys, icol_hsys, cval_hsys, cmtmp3, nrho, nrho, nrho, dsmall_spa)
if (nnz_hsys .le. 0 .or. nnz_hsys .gt. nrho2) then
   write(6,*)'prelude_spa: error nnz_hsys for the present job ', nnz_hsys
   stop
end if
write(6,*)'prelude_spa: nnz_hsys = ', nnz_hsys
call flush(6)
!
iread_spa = 0
if (lexist) then
   read(52)ntmp1
   if (ntmp1 .le. 0 .or. ntmp1 .gt. nrho2) then
      write(6,*)'prelude_spa: error nnz_hsys read from <sparse_index.data> ', ntmp1
      write(6,*)'prelude_spa: abort reading '
      close(52)
      close(53)
      lexist = .false.
   else
      allocate(nvec1(ntmp1), nvec2(ntmp1), cvec1(ntmp1), STAT=istat)
      do ni=1,ntmp1
         read(52)nvec1(ni), nvec2(ni), cvec1(ni)
      end do
      call zmat_coo2dns(ntmp1, nvec1, nvec2, cvec1, cmtmp1, nrho, nrho, nrho)
      call zmat_comp_sparsity(cmtmp1, nrho, nrho, nrho, cmtmp3, nrho, nrho, nrho, &
                              isame, icomp, na, nb, dsmall_spa) 
      deallocate(nvec1, nvec2, cvec1, STAT=istat)
      if (isame .eq. 1) then
         iread_spa = 2
      else if (icomp .eq. 1) then
         iread_spa = 1
      else
         write(6,*)'prelude_spa: incompatible <sparse_index.data> '
         write(6,*)'prelude_spa: abort reading '
         close(52)
         close(53)
         lexist = .false.
      end if 
   end if  
end if
!
if (lexist) then
   if (iread_spa .eq. 2) then  ! read sparsity info for present job
      read(53)lunk_spa
      allocate(nnz_spa(nunk), ind_spa(nunk), STAT=istat)
      allocate(irow_spa(lunk_spa), icol_spa(lunk_spa), rho_spa(lunk_spa), STAT=istat)
      do lni=1,nunk
         read(52)nnz_spa(lni), ind_spa(lni)
      end do
      close(unit=52)
      do lni=1,lunk_spa
         read(53)irow_spa(lni), icol_spa(lni)
      end do
      close(unit=53)
      write(6,*)'prelude_spa: <sparse_index.data> and <sparse_info.data> read successfully'
      call flush(6)
      goto 101
   else ! read sparsity info for a previous 'reference' job
        ! note added on Jan 24, 2020:
        ! the reference job should have same number of ADOs (nunk)
        ! and those ADOs must be more sparse than the present job (lunk_spa0 < lunk)
      read(53)lunk_spa0
      allocate(nnz_spa0(nunk), ind_spa0(nunk), STAT=istat)
      allocate(irow_spa0(lunk_spa0), icol_spa0(lunk_spa0), STAT=istat)
      do lni=1,nunk
         read(52)nnz_spa0(lni), ind_spa0(lni)
      end do
      close(unit=52)
      do lni=1,lunk_spa0
         read(53)irow_spa0(lni), icol_spa0(lni)
      end do
      close(unit=53)
      write(6,*)'prelude_spa: <sparse_index.data> and <sparse_info.data> read successfully '
      write(6,*)'prelude_spa: but the information cannot be used straightforwardly         '
      call flush(6)
   end if
end if
!
lunk_spa = 0
allocate(ind_spa(nunk), nnz_spa(nunk), STAT=istat)
nk = 0
do nj=1,nrho
   do ni=1,nrho
!      if (nj .eq. 1 .and. ni .eq. 1) then
!          nk = nk + 1
!          cycle
!      end if
      if ( cdabs(cmtmp3(ni,nj)) .gt. dsmall_spa ) then
         nk = nk + 1
      end if
   end do
end do
lunk_spa   = nk
ind_spa(1) = 1
nnz_spa(1) = nk
allocate(irow_spa(lunk_spa), icol_spa(lunk_spa), STAT=istat)
na = 0
spa_rho: do nj=1,nrho
   do ni=1,nrho
      if (na .ge. nk) exit spa_rho
      if ( cdabs(cmtmp3(ni,nj)) .gt. dsmall_spa) then
!      if ( cdabs(cmtmp3(ni,nj)) .gt. dsmall_spa .or.     &
!           (nj .eq. 1 .and. ni .eq. 1) ) then    ! (1,1) element is included (to allow for ltrend=F)
         na = na + 1
         irow_spa(na) = ni
         icol_spa(na) = nj
      end if
   end do
end do spa_rho
!
write(6,*)'prelude_spa: nnz(rhosys) = ', nnz_spa(1)
write(6,*)'prelude_spa: (irow, icol) of nonzero elements of rhosys '
lni = ind_spa(1)
do ni=1,na
   write(6,*)ni, irow_spa(lni-1+ni), icol_spa(lni-1+ni)
end do
call flush(6)
!
!open(unit=48, file='rhocoo.tmp', form='unformatted', status='unknown', access='sequential')
open(unit=48, file='rhocoo.tmp', form='binary', status='unknown', access='sequential')
rewind(48)
do lni=1,lunk_spa
   write(48)irow_spa(lni), icol_spa(lni)
end do
!open(unit=49, file='rhoval.tmp', form='unformatted', status='unknown', access='sequential')
open(unit=49, file='rhoval.tmp', form='binary', status='unknown', access='sequential')
rewind(49)
do lni=1,lunk_spa
   write(49)cmtmp3(irow_spa(lni),icol_spa(lni))
end do
!
lunk_last = lunk_spa
deallocate(irow_spa, icol_spa, STAT=istat)
!
! from second to terminal tier
!
nk = nrho**2
allocate(rowcoo(nk), colcoo(nk), valcoo(nk), STAT=istat)
!
write(6,*)
write(6,*)'prelude_spa: running through all the tiers '
call flush(6)
!
do itier=2,ntier0
   call cpu_time(cpu1)
   allocate(irow_spa(lunk_last), icol_spa(lunk_last), val_tmp(lunk_last), STAT=istat)
   rewind(48)
   do lni=1,lunk_last
      read(48)irow_spa(lni), icol_spa(lni)
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
!         lnj    = indexcoef(lnl)
!         itype  = itypecoef(lnl)
         lnl    = lnl + 1
         lstart = ind_spa(lnj)
         nnz    = nnz_spa(lnj)
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
                call zamsmm_coo_proj(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                     icol_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
            else
                call zamsmm_coo(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                icol_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
            end if
! cmtmp5 : rho*ams
            sidea  = 'R'
            if (lprams) then
                call zamsmm_coo_proj(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                     icol_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
            else
                call zamsmm_coo(sidea, transa, transb, iorbs, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                icol_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
            end if
!
! cmtmp1 : exp(-Heff) * (ams*rho)  
! cmtmp2 : exp(-Heff) * (rho*ams)
            matdescra(1) = 'g'
            matdescra(4) = 'f'
            call mkl_zcoomm('n', nrho, nrho, nrho, cunity, matdescra, val_tmp(1), irow_spa(1), icol_spa(1), &
                            nnz_spa(1), cmtmp4, nrho, czero, cmtmp1, nrho)
            call mkl_zcoomm('n', nrho, nrho, nrho, cunity, matdescra, val_tmp(1), irow_spa(1), icol_spa(1), &
                            nnz_spa(1), cmtmp5, nrho, czero, cmtmp2, nrho)
! avoid accidental cancellation between rho*ams and ams*rho
            do
               cmix1 = dcmplx(random_normal(), random_normal())
               cmix2 = dcmplx(random_normal(), random_normal())
               if (cdabs(cmix1) .gt. micron .and. cdabs(cmix2) .gt. micron) exit
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
                   call zamsmm_coo_proj(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                        icol_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
               else
                   call zamsmm_coo(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                   icol_spa(lstart), czero, cmtmp4, nrho, zmtmp1, nrho)
               end if
               sidea  = 'R'
               if (lprams) then
                   call zamsmm_coo_proj(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                        icol_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
               else
                   call zamsmm_coo(sidea, transa, transb, iorbs2, ispin, cunity, val_tmp(lstart), nnz, irow_spa(lstart), &
                                   icol_spa(lstart), czero, cmtmp5, nrho, zmtmp1, nrho)
               end if
               matdescra(1) = 'g'
               matdescra(4) = 'f'
               call mkl_zcoomm('n', nrho, nrho, nrho, cunity, matdescra, val_tmp(1), irow_spa(1), icol_spa(1), &
                               nnz_spa(1), cmtmp4, nrho, czero, cmtmp1, nrho)
               call mkl_zcoomm('n', nrho, nrho, nrho, cunity, matdescra, val_tmp(1), irow_spa(1), icol_spa(1), &
                               nnz_spa(1), cmtmp5, nrho, czero, cmtmp2, nrho)
               do
                  cmix1 = dcmplx(random_normal(), random_normal())
                  cmix2 = dcmplx(random_normal(), random_normal())
                  if (cdabs(cmix1) .gt. micron .and. cdabs(cmix2) .gt. micron) exit
               end do
               zmtmp2(1:nrho,1:nrho) = zmtmp2(1:nrho,1:nrho) + cmix1 * cmtmp1(1:nrho,1:nrho) + cmix2 * &
                                       cmtmp2(1:nrho,1:nrho)
            end do ! end of loop over iorbs2
            dtmp2 = min(dtmp2, dtmp1)
         end if
      end do ! end of loop over ni
!
!if (dtmp2 .lt. dsmall) then
!    write(99,*)itier, lni, dtmp2
!    call flush(99)
!end if
      if (dtmp2 .gt. dsmall_spa) then
          call zmat_dns2coo(na, rowcoo, colcoo, valcoo, zmtmp2, nrho, nrho, nrho, dsmall_spa)
      else
          na = 0
      end if
      lunk_spa = lunk_spa + na
      nnz_spa(lni) = na
      ind_spa(lni) = ind_spa(lni-1) + nnz_spa(lni-1)
      if (na .gt. 0) then
         do ni=1,na
            write(48)rowcoo(ni), colcoo(ni)
            write(49)valcoo(ni)
         end do
      end if
   end do 
!
   lunk_last = lunk_spa
   deallocate(irow_spa, icol_spa, val_tmp, STAT=istat)
!
   call cpu_time(cpu2)
   write(6,100)itier, lunk_spa, cpu2 - cpu1
   call flush(6)
end do
100 format('itier= ', I3, ' size= ', I20, 2x, f12.4)
deallocate(rowcoo, colcoo, valcoo, STAT=istat)
!
allocate(irow_spa(lunk_last), icol_spa(lunk_last), rho_spa(lunk_last), STAT=istat)
rewind(48)
do lni=1,lunk_last
   read(48)irow_spa(lni), icol_spa(lni)
end do
rewind(49)
do lni=1,lunk_last
   read(49)rho_spa(lni)
end do
close(unit=48, status='delete')
close(unit=49, status='delete')
!
101 continue
!
deallocate(cexphs, cexpheff, STAT=istat)
deallocate(hybrid, STAT=istat)
deallocate(cmat1, cmat2, STAT=istat)
!
! account for memory of used arrays 
! including: ind_spa, nnz_spa, irow_spa, icol_spa, rho_spa
!        and ind_spa0, nnz_spa0, irow_spa0, icol_spa0 (if iread_spa == 1)
!
dtmp1 = dble(nunk * (4 + 8) + lunk_spa * 2 * 4 + lunk_spa * 16) / dble(1024**2)
if (iread_spa .eq. 1) then
   dtmp1 = dtmp1 + dble(nunk * (4 + 8) + lunk_spa0 * 2 * 4) / dble(1024**2)
end if
memory = memory + dtmp1
write(6,*)'prelude_spa: memory allocated for sparsity info ', dtmp1, ' MB '
call flush(6)
!
lni = nrho**2 * nunk 
write(6,*)
write(6,*)'prelude_spa: total unknown variables = ', lni
write(6,*)'prelude_spa: nonzero variables       = ', lunk_spa
write(6,*)'prelude_spa: sparsity ratio for rho  = ', dble(lunk_spa) / dble(lni)
call flush(6)
!
! output files 
! 
if (iread_spa .ne. 2) then
   open(unit=52, file='sparse_index.data', form='binary')
   rewind(52)
   write(52)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow
   write(52)nunk, offcor, lequileads, lsimple, lscreen, lproj
   write(52)lprocc, nproji, nprojf, nprocc
   write(52)nnz_hsys
   do ni=1,nnz_hsys
      write(52)irow_hsys(ni), icol_hsys(ni), cval_hsys(ni)
   end do
   do lni=1,nunk
      write(52)nnz_spa(lni), ind_spa(lni)
   end do
   close(52)
!
   open(unit=53, file='sparse_info.data', form='binary')
   rewind(53)
   write(53)nunk, offcor, lequileads, lsimple, lscreen, lproj
   write(53)lunk_spa
   do lni=1,lunk_spa
      write(53)irow_spa(lni), icol_spa(lni)
   end do
   close(53)
end if
!
! Debug USE : compare rho and rho_spa, to verify the locations of nonzero elements
!
!lsame = .true.
!lcomp = .true.
!check1: do itier=2,ntier0
!   do lni=nfirst(itier),nlast(itier)
!      cmtmp1(1:nrho,1:nrho) = rho(1:nrho,1:nrho,lni) 
!      lstart = ind_spa(lni)
!      call zmat_coo2dns(nnz_spa(lni), irow_spa(lstart), icol_spa(lstart), rho_spa(lstart), &
!                        cmtmp2, nrho, nrho, nrho)
!      call zmat_comp_sparsity(cmtmp1, nrho, nrho, nrho, cmtmp2, nrho, nrho, nrho, ni, nj, na, nb, dsmall_spa)
!      if (ni .eq. 0) lsame = .false.
!      if (nj .eq. 0) lcomp = .false.
!      if (.not. lcomp) exit check1
!   end do
!end do check1
!write(6,*)
!if (lsame) then
!   write(6,*)'prelude_spa: rho and rho_spa have same sparsity topology '
!else
!   if (lcomp) then
!      write(6,*)'prelude_spa: sparsity of rho is covered by rho_spa '
!   else
!      write(6,*)'prelude_spa: attention! sparsity of rho differs from rho_spa '
!   end if
!end if
!call flush(6)
!
111 continue
!
write(6,*)'prelude_spa: leaving subroutine '
call flush(6)
!
return
end subroutine prelude_spa
