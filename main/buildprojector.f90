subroutine buildprojector
use matmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: info, istat, nwork, ncount
integer :: ni, nj, nk
logical :: lherm
real*8  :: docc0, tolocc
real*8  :: dtmp1, dtmp2
complex*16 :: ctmp1
real*8,     allocatable   :: dnocc(:,:)
real*8,     allocatable   :: rwork(:), tmpval(:)
complex*16, allocatable   :: zwork(:,:)
complex*16, allocatable   :: cmat1(:,:), cmat2(:,:), cmat3(:,:), cmat4(:,:)
!
!write(6,*)
!write(6,*)'buildprojector: entering '
!call flush(6)
if (.not. lproj) then
    write(6,*)'buildprojector: error! wrong entry. lproj=F found ', lproj
    stop
end if
tolocc = 1.d-6
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), cmat3(nrho,nrho), cmat4(nrho,nrho), STAT=istat)
!
cmat1(1:nrho,1:nrho) = czero
lprths = .false.
!lprths = .true. 
call calchs(cmat1)
cmat4(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
call checkhermicity(cmat4, nrho, cmat1, nrho, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)'buildprojector: error! check Hermicity fails 1 ', dtmp1
!    call cmatout(nrho, nrho, cmat4, cmat1(1,1))
    stop
end if
!
allocate(zwork(nrho,nrho), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
nwork = nrho**2
cmat1(1:nrho,1:nrho) = cmat4(1:nrho,1:nrho)
!
call zheev('V', 'u', nrho, cmat1, nrho, tmpval, zwork, nwork, rwork, info)
if (info .ne. 0) then
    write(6,*)'buildprojector: error! diagonalization failed 1 ', info
    stop
end if
zwork(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho)
!
cmat2(1:nrho,1:nrho) = czero
if (lprocc) then
    allocate(dnocc(nspin,nrho), STAT=istat)
    ncount = 0
    do nk=1,nrho
       ! rho_k = |psi_k> <psi_k|
       do nj=1,nrho
          do ni=1,nrho
             cmat1(ni,nj) = zwork(ni,nk) * dconjg(zwork(nj,nk))
          end do
       end do
       ! occupation number
       call getocc(dnocc(1,nk), cmat1, cmat3, cmat4)
       if (nspin .eq. 2) then
           dnocc(1,nk) = dnocc(1,nk) + dnocc(2,nk)
       end if
       if (nk .eq. 1) then
           docc0 = dnocc(1,nk)
           cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + cmat1(1:nrho,1:nrho)
           ncount = ncount + 1
           cycle
       end if
! Note (2022-06-04): assume all occupation numbers are integers, 
! but in many cases fractional occupation numbers may show up due to mixing of
! degenerate states. We cannot deal with such cases in the present version
       if ( dabs(docc0 - dnocc(1,nk)) .le. dble(nprocc) + tolocc ) then
           cmat2(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho) + cmat1(1:nrho,1:nrho)
           ncount = ncount + 1
       end if
    end do
else
    do nk=nproji,nprojf,1
       ! rho_k = |psi_k> <psi_k|
       do nj=1,nrho
          do ni=1,nrho
             cmat2(ni,nj) = cmat2(ni,nj) + zwork(ni,nk) * dconjg(zwork(nj,nk))
          end do
       end do
    end do
end if
!
allocate(cproj(nrho,nrho), cproj_bar(nrho,nrho), STAT=istat)
cproj(1:nrho,1:nrho) = cmat2(1:nrho,1:nrho)
cmat1 = czero
do ni=1,nrho
   cmat1(ni,ni) = cunity
end do
cproj_bar(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) - cproj(1:nrho,1:nrho)
!
lprgen = .true.
ldopro = .true.
!
! output
! 
write(6,*)
write(6,*)'buildprojector: all eigenstates of system Hamiltonian '
write(6,*)' nk    eig  '
do nk=1,nrho
   write(6,1000)nk, tmpval(nk)*hbar
end do
call flush(6)
1000 format(I4, 1x, f10.4, 1x, 7(f7.3, 1x))
!
if (lprocc) then
    write(6,*)'buildprojector: use occupation number as criterion '
    write(6,*)' from nocc= ', max(int(docc0)-nprocc, 0), ' to ', int(docc0)+nprocc
    nstate_proj = ncount
    allocate(istate_proj(nstate_proj), STAT=istat)
    ncount = 0
    do nk=1,nrho
       if ( dabs(docc0 - dnocc(1,nk)) .le. dble(nprocc) + tolocc ) then
           ncount = ncount + 1
           istate_proj(ncount) = nk
       end if
    end do
    if (ncount .ne. nstate_proj) then
        write(6,*)'buildprojector: error! ncount != nstate_proj found '
        write(6,*)ncount, nstate_proj
        stop
    end if
    write(6,*)'buildprojector: projector includes the following eigenstates: '
    do ncount=1,nstate_proj
       nk = istate_proj(ncount)
       write(6,1000)nk, dnocc(1,nk)
    end do
else
    write(6,*)'buildprojector: projector from ', nproji, ' to ', nprojf
    nstate_proj = nprojf - nproji + 1
    allocate(istate_proj(nstate_proj), STAT=istat)
    ncount = 0
    do nk=nproji,nprojf,1
       ncount = ncount + 1
       istate_proj(ncount) = nk
    end do
    write(6,*)'buildprojector: projector includes the following eigenstates: '
    do ncount=1,nstate_proj
       nk = istate_proj(ncount)
       write(6,1000)nk, tmpval(nk)*hbar
    end do
end if
call flush(6)
!
call checkamsprojector
!
deallocate(cmat1, cmat2, cmat3, cmat4, STAT=istat)
deallocate(zwork, rwork, tmpval, STAT=istat)
if (allocated(dnocc)) deallocate(dnocc, STAT=istat)
!
return
end subroutine buildprojector
