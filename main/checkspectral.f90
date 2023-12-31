subroutine checkspectral
use matmod
use tmpmatmod
implicit none
include "../include/sizes"
include "../include/common"
!
! Supposedly dxip (proportional to spectral density function)
!  should be Hermitian and non-negative-definite.
! Here, dxip is real, so it is a symmetric matrix.
!
! note: check only ialf=1&2
!
integer :: ni, nj, info, istat
logical :: lcheck
real*8  :: tmpval(maxorb), dwork(5*maxorb)
real*8, allocatable :: dmat1(:,:)
!
if (lspcor) then
    write(6,*)'checkspectral: lspcor=T found '
    write(6,*)'checkspectral: code not ready yet. exiting... '
    call flush(6)
    return
end if
!
allocate(dmat1(norbs,norbs), STAT=istat)
!
lcheck = .true.
dxip(1:nvar2,1:nspin,1:nalf) = dble(zxip(1:nvar2,1:nspin,1:nalf))
if (offcor) then
  do nj=1,norbs
    do ni=1,norbs
      dmtmp1(ni,nj) = dxip((nj-1)*norbs+ni, 1, 1)
      if (nalf .ne. 1) dmtmp2(ni,nj) = dxip((nj-1)*norbs+ni, 1, 2)
    end do
  end do
else
  dmtmp1 = 0.d0
  if (nalf .ne. 1) dmtmp2 = 0.d0
  do ni=1,norbs
    dmtmp1(ni,ni) = dxip(ni, 1, 1)
    if (nalf .ne. 1) dmtmp2(ni,ni) = dxip(ni, 1, 2)
  end do
end if
dmtmp1 = dmtmp1 * 2.d0 / bandwidth(1)**2
if (nalf .ne. 1) dmtmp2 = dmtmp2 * 2.d0 / bandwidth(2)**2
!
! write system-bath coupling
! 
dmat1(1:norbs,1:norbs) = dmtmp1(1:norbs,1:norbs)
write(6,*)
write(6,*)' coupling matrix for sytem-lead_L '
call amatout(norbs, norbs, dmat1)  ! show full matrix
!
do ni=1,norbs 
  do nj=1,ni-1
    if (dabs(dmtmp1(ni,nj) - dmtmp1(nj,ni)) .ge. dnano) then
      write(6,*)' symmetry error ', ni, nj, dmtmp1(ni,nj), dmtmp1(nj,ni) 
      stop
    end if
  end do
end do
if (nalf .ne. 1) then
  dmat1(1:norbs,1:norbs) = dmtmp2(1:norbs,1:norbs)
  write(6,*)
  write(6,*)' coupling matrix for sytem-lead_R '
  call amatout(norbs, norbs, dmat1)
  do ni=1,norbs 
    do nj=1,ni-1
      if (dabs(dmtmp2(ni,nj) - dmtmp2(nj,ni)) .ge. dnano) then
        write(6,*)' symmetry error ', ni, nj, dmtmp2(ni,nj), dmtmp2(nj,ni) 
        stop
      end if
    end do
  end do
end if
call flush(6)
!
write(6,*)
write(6,*)' check nonnegative-definite property for SL '
if (.not. offcor) then
  do ni=1,norbs
    if (dmtmp1(ni,ni) .lt. 0.d0) then
      write(6,*)' error! negative diag. ', ni, dmtmp1(ni,ni)
      lcheck = .false.
    end if
    linewidth(1,ni) = dmtmp1(ni,ni)
  end do
else
  call dsyev('N', 'U', norbs, dmtmp1, nrho, tmpval, dwork, 5*maxorb, info)
  if (info .ne. 0) then
    write(6,*)
    write(6,*)' error! diagonalization failed in checkspectral ', info
    stop
  end if
  do ni=1,norbs
    if (tmpval(ni) .lt. 0.d0) then
      write(6,*)' error! negative eigen. ', ni, tmpval(ni)
      lcheck = .false.
    end if
    linewidth(1,ni) = tmpval(ni)
  end do
end if
call flush(6)
!
if (nalf .ne. 1) then
  write(6,*)' check nonnegative-definite property for SR '
  if (.not. offcor) then
    do ni=1,norbs
      if (dmtmp2(ni,ni) .lt. 0.d0) then
        write(6,*)' error! negative diag. ', ni, dmtmp2(ni,ni)
        lcheck = .false.
      end if
      linewidth(2,ni) = dmtmp2(ni,ni)
    end do
  else
    call dsyev('N', 'U', norbs, dmtmp2, nrho, tmpval, dwork, 5*maxorb, info)
    if (info .ne. 0) then
      write(6,*)
      write(6,*)' error! diagonalization failed in checkspectral ', info
      stop
    end if
    do ni=1,norbs
      if (tmpval(ni) .lt. 0.d0) then
        write(6,*)' error! negative eigen. ', ni, tmpval(ni)
        lcheck = .false.
      end if
      linewidth(2,ni) = tmpval(ni)
    end do
  end if
end if
call flush(6)
!
write(6,*)' checkspectral complete! '
if (lcheck) then
  write(6,*)' dxip is verified to be OK! linewidth are as follows '
  write(6,1001)(linewidth(1,ni), ni=1,norbs)
  if (nalf .ne. 1) write(6,1002)(linewidth(2,ni), ni=1,norbs)
  call flush(6)
else 
  write(6,*)' error! negative eigenvalue was found! check your input! '
  stop
end if
!
1001 format('  left ', 1x, 6(e14.6e2, 2x))
1002 format(' right ', 1x, 6(e14.6e2, 2x))
!
deallocate(dmat1, STAT=istat)
return
!
end subroutine checkspectral
