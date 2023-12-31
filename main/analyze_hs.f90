subroutine analyze_hs(zin)
use matmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in) :: zin(nrho,*)
integer :: info, istat, nwork, ncount
integer :: ni, nj, nk, nm, iorbs, ispin, nxyz, nnz
integer :: n1, n2
integer*8 :: lni, lnj
logical :: lherm
real*8  :: dtmp1, dtmp2, emin
complex*16 :: ctmp1
real*8,     allocatable   :: dnocc(:,:), ds2(:,:), dsz(:,:), dsztmp(:), ds12(:)
real*8,     allocatable   :: rwork(:), tmpval(:)
complex*16, allocatable   :: zwork(:,:), cvec1(:)
complex*16, allocatable   :: cmat1(:,:), cmat2(:,:), cmat3(:,:), cmat4(:,:)
!
!write(6,*)
!write(6,*)'analyze_hs: entering '
!call flush(6)
if (.not. lanahs) then
    write(6,*)'analyze_hs: error! wrong entry. lanahs=F found ', lanahs
    stop
end if
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), cmat3(nrho,nrho), cmat4(nrho,nrho), STAT=istat)
call checkhermicity(zin, nrho, cmat1, nrho, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)'analyze_hs: error! check Hermicity fails 1 ', dtmp1
    stop
end if
if (nspin .ne. 2) then
    write(6,*)'analyze_hs: error! code only works for nspin=2 ', nspin
    stop
end if
!
allocate(zwork(nrho,nrho), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
allocate(cvec1(nrho), STAT=istat)
allocate(dnocc(nspin,nrho), STAT=istat)
nxyz = 3
if (lspin3d) then
    allocate(ds2(nxyz,nrho), dsz(nxyz,nrho), dsztmp(nxyz), STAT=istat)
    if (norbs .eq. 2) allocate(ds12(nrho), STAT=istat)
end if
!
nwork = nrho**2
cmat1(1:nrho,1:nrho) = zin(1:nrho,1:nrho)
!
call zheev('V', 'u', nrho, cmat1, nrho, tmpval, zwork, nwork, rwork, info)
if (info .ne. 0) then
    write(6,*)'analyze_hs: error! diagonalization failed 1 ', info
    stop
end if
! sort eigenvalues from negative to positive
!write(6,*)
!write(6,*)'before sort'
!write(6,*)(tmpval(ni), ni=1,nrho)
!do ni=1,nrho-1
!   do nj=ni+1,nrho
!      if (tmpval(ni) .gt. tmpval(nj)) then
!          dtmp1      = tmpval(ni)
!          tmpval(ni) = tmpval(nj)
!          tmpval(nj) = dtmp1
!          cvec1(1:nrho)    = cmat1(1:nrho,ni)
!          cmat1(1:nrho,ni) = cmat1(1:nrho,nj)
!          cmat1(1:nrho,nj) = cvec1(1:nrho)
!      end if
!   end do
!end do
!write(6,*)'after  sort'
!write(6,*)(tmpval(ni), ni=1,nrho)
!call flush(6)
zwork(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho)
!
do nk=1,nrho
   ! rho_k = |psi_k> <psi_k|
   do nj=1,nrho
      do ni=1,nrho
         cmat2(ni,nj) = zwork(ni,nk) * dconjg(zwork(nj,nk))
      end do
   end do
   ! occupation number
   call getocc(dnocc(1,nk), cmat2, cmat3, cmat4)
   !
   if (lspin3d) then
       dsz(1:nxyz,nk) = 0.d0
       do iorbs=1,norbs
          call calcspinmoment(iorbs, dsztmp, cmat2)
          dsz(1:nxyz,nk) = dsz(1:nxyz,nk) + dsztmp(1:nxyz)
       end do
       call calctotalspin2(ds2(1,nk), cmat2)
       !
       if (norbs .eq. 2) then
           call calcspinproduct(1, 2, ds12(nk), cmat2)
       end if
   end if
end do
!
! output
! 
!write(6,*)'before output'
!write(6,*)(tmpval(ni), ni=1,nrho)
!call flush(6)
!
write(6,*)
write(6,*)'analyze_hs: analyze all eigenstates of system Hamiltonian '
if (lspin3d) then
    if (norbs .eq. 2) then
        write(6,*)' nk    eig      occ_u    occ_d     occ    <Sz>    <S^2>    <S1*S2> '
        do nk=1,nrho
           write(6,1000)nk, tmpval(nk)*hbar, dnocc(1,nk), dnocc(2,nk), dnocc(1,nk)+dnocc(2,nk),  &
                        dsz(nxyz,nk), ds2(1,nk)+ds2(2,nk)+ds2(3,nk), ds12(nk)
        end do
    else
        write(6,*)' nk    eig      occ_u    occ_d     occ    <Sz>    <S^2> '
        do nk=1,nrho
           write(6,1000)nk, tmpval(nk)*hbar, dnocc(1,nk), dnocc(2,nk), dnocc(1,nk)+dnocc(2,nk),  &
                        dsz(nxyz,nk), ds2(1,nk)+ds2(2,nk)+ds2(3,nk)
        end do
    end if
else
    write(6,*)' nk    eig      occ_u    occ_d     occ '
    do nk=1,nrho
       write(6,1000)nk, tmpval(nk)*hbar, dnocc(1,nk), dnocc(2,nk), dnocc(1,nk)+dnocc(2,nk)
    end do
end if
call flush(6)
1000 format(I4, 1x, f10.4, 1x, 7(f7.3, 1x))
!
! turn off lbzman and lhseig if lreadi=T 
!
if (lreadi) then
    lbzman = .false.
    lhseig = .false.
    write(6,*)
    write(6,*)'analyze_hs: lreadi=T found '
    write(6,*)'analyze_hs: override lbzman and lhseig '
    call flush(6)
end if
!
if (lbzman) then
    write(6,*)
    write(6,*)'analyze_hs: lanahs=T and lbzman=T found '
    write(6,*)'analyze_hs: use Boltzmann distribution to populate states '
    write(6,*)'analyze_hs: set lhseig = F '
    call flush(6)
    lhseig = .false.
    !
    if (lproj) then
        n1 = 1
        n2 = nstate_proj
        emin = tmpval(istate_proj(n1))
    else
        emin = tmpval(1)
        n1 = 1
        n2 = nrho
    end if
    cmat2(1:nrho,1:nrho) = czero
    do nm=n1,n2,1
       if (lproj) then
           nk = istate_proj(nm)
       else
           nk = nm
       end if
       dtmp1 = dexp(-dbeta(1) * (tmpval(nk) - emin) * hbar)
       do nj=1,nrho
          ctmp1 = dconjg(zwork(nj,nk)) * dtmp1
          do ni=1,nrho
             cmat2(ni,nj) = cmat2(ni,nj) + zwork(ni,nk) * ctmp1
          end do
       end do
    end do
    dtmp2 = 0.d0
    do ni=1,nrho
       dtmp2 = dtmp2 + dble(cmat2(ni,ni))
    end do
    cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) / dtmp2
!write(6,*)
!write(6,*)'cmat2 -- checkpoint 1 '
!call cmatout(nrho, nrho, cmat2, cmat1)
end if
!
! use eigenstate of hs to construct initial reduced density matrix
!
if (lhseig) then
    write(6,*)
    write(6,*)'analyze_hs: lanahs=T and lhseig=T found '
    if (nhseig .ge. 1 .and. nhseig .le. nrho) then
        write(6,*)'analyze_hs: use eigenstate No. ', nhseig, ' of system Hamiltonian '
        write(6,*)'            to build initial reduced density matrix '
        call flush(6)
    else
        write(6,*)'analyze_hs: error! wrong nhseig = ', nhseig
        stop
    end if
    if (lproj) then
        n1 = 0
        do nm=1,nstate_proj
           if (nhseig .eq. istate_proj(nm)) then
               n1 = nm
           end if
        end do
        if (n1 .eq. 0) then
            write(6,*)'analyze_hs: error! nhseig outside range of projector '
            write(6,*)'analyze_hs: nhseig, istate_proj '
            write(6,*)nhseig, (istate_proj(nm), nm=1,nstate_proj)
            stop
        end if
    end if
    if (lhsdeg) then
        write(6,*)'analyze_hs: consider degeneracy of system Hamiltonian '
        if (toldeg .lt. 0) then
            write(6,*)'analyze_hs: error! wrong toldeg = ', toldeg
            stop
        end if
        write(6,*)'analyze_hs: tolerance for degeneracy is toldeg = ', toldeg
        call flush(6)
    end if
    !
    ncount = 1
    do nj=1,nrho
       do ni=1,nrho
          cmat2(ni,nj) = zwork(ni,nhseig) * dconjg(zwork(nj,nhseig))
       end do
    end do
    write(6,*)
    write(6,*)'analyze_hs: orbital No. ', nhseig, ' adopted '
    call flush(6)
    if (lhsdeg) then
        do nk=1,nrho
           if (nk .eq. nhseig) cycle
           if ( dabs(tmpval(nk) - tmpval(nhseig))*hbar .gt. toldeg ) cycle
           if (lproj) then
               n1 = 0
               do nm=1,nstate_proj
                  if (nk .eq. istate_proj(nm)) then
                      n1 = nm
                      exit
                  end if
               end do
               if (n1 .eq. 0) cycle
           end if
           do nj=1,nrho
              do ni=1,nrho
                 cmat2(ni,nj) = cmat2(ni,nj) + zwork(ni,nk) * dconjg(zwork(nj,nk))
              end do
           end do
           ncount = ncount + 1
           write(6,*)'analyze_hs: orbital No. ', nk, ' adopted '
           call flush(6)
        end do
    end if
    cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) / dble(ncount)
end if
!
if (lbzman .or. lhseig) then
    call checkhermicity(cmat2, nrho, cmat1, nrho, lherm, dtmp1)
    if (.not. lherm) then
        write(6,*)'analyze_hs: error! cmat2 not Hermitian '
        stop
    end if
    dtmp1 = 0.d0
    do ni=1,nrho
       dtmp1 = dtmp1 + dble(cmat2(ni,ni))
    end do
    if (dabs(dtmp1-1.d0) .gt. dpico) then
        write(6,*)'analyze_hs: error! trace(cmat2)=1 not satisfied ', dtmp1
        stop
    end if
!write(6,*)
!write(6,*)'cmat2 '
!call cmatout(nrho, nrho, cmat2, cmat1)
    if (lsparse) then
        nnz = nnz_spa(1)
        lni = ind_spa(1)
        do ni=1,nnz
           lnj = lni - 1 + ni
           rho_spa(lnj) = cmat2(irow_spa(lnj), icol_spa(lnj))
        end do
        write(6,*)'analyze_hs: rho_spa updated '
        call flush(6)
    else
        rho(1:nrho,1:nrho,1) = cmat2(1:nrho,1:nrho)
        write(6,*)'analyze_hs: rho updated '
        call flush(6)
    end if
    call refresh_rhosys
!    write(6,*)'analyze_hs: rhosys: ' 
!    call cmatout(nrho, nrho, rhosys, cmat1)
!    call flush(6)
end if
!
deallocate(cmat1, cmat2, cmat3, cmat4, STAT=istat)
deallocate(zwork, rwork, tmpval, STAT=istat)
deallocate(cvec1, STAT=istat)
deallocate(dnocc, STAT=istat)
if (lspin3d) then
    deallocate(ds2, dsz, dsztmp, STAT=istat)
    if (norbs .eq. 2) deallocate(ds12, STAT=istat)
end if
!
return
end subroutine analyze_hs
!
subroutine getocc(docc, rhoinp, cmat1, cmat2)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in) :: rhoinp(nrho,*)
complex*16, intent(inout) :: cmat1(nrho,*), cmat2(nrho,*)
real*8, intent(out) :: docc(*)
integer :: ispin, iorbs, ni
!
docc(1:nspin) = czero
do ispin=1,nspin
   do iorbs=1,norbs
      call zamsmm1('l', 'n', 'n', iorbs, ispin, cunity, rhoinp, nrho, czero, cmat1, nrho)
      call zamsmm1('l', 'c', 'n', iorbs, ispin, cunity,  cmat1, nrho, czero, cmat2, nrho)
      do ni=1,nrho
         docc(ispin) = docc(ispin) + dble(cmat2(ni,ni))
      end do
   end do
end do
!
return
end subroutine getocc

