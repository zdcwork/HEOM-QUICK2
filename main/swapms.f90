subroutine swapms(nrho, norbs, nspin, irowin, icolin, iorbs, ispin, &
                  irowout, icolout, isgnout)
implicit none
!
integer, intent(in)  :: nrho, norbs, nspin, iorbs, ispin, irowin, icolin
integer, intent(out) :: irowout, icolout, isgnout
integer :: ni, nj, nk, len0, istat
integer :: isgrow, isgcol
integer, allocatable :: iocc0(:,:)
!
irowout = 0
icolout = 0
isgnout = 0   ! -1, reserves sign 
              !  1, conserves sign
isgrow = 0
isgcol = 0
!
allocate(iocc0(norbs,nspin), STAT=istat)
!
if (nspin .eq. 1) then
!
! row
!
! figure out occupation pattern for row index
  do ni=1,norbs
    len0 = nrho / 2**(norbs - ni + 1)
    if (mod(irowin, len0) .eq. 0) then
      if (mod(irowin / len0, 2) .eq. 0) then
       iocc0(ni,1) = 1
      else 
       iocc0(ni,1) = 0
      end if
    else
      if (mod(irowin / len0, 2) .eq. 0) then
       iocc0(ni,1) = 0
      else
       iocc0(ni,1) = 1
      end if
    end if
  end do
  if (iocc0(iorbs-1,1) .eq. 1 .and. iocc0(iorbs,1) .eq. 1) then
    isgrow = -1
  else
    isgrow = 1
  end if
! swap occupation number for (iorb-1)th and (iorb)th orbitals
  nj = iocc0(iorbs-1,1)
  iocc0(iorbs-1,1) = iocc0(iorbs,1)
  iocc0(iorbs,1) = nj
! calculate row index corresponding to new occupation pattern
  nj = 1
  do ni=norbs,1,-1
    len0 = nrho / 2**(norbs - ni + 1)
    if (iocc0(ni,1) .eq. 1) nj = nj + len0
  end do
  irowout = nj
!
! column
!
! figure out occupation pattern for col index
  do ni=1,norbs
    len0 = nrho / 2**(norbs - ni + 1)
    if (mod(icolin, len0) .eq. 0) then
      if (mod(icolin / len0, 2) .eq. 0) then
       iocc0(ni,1) = 1
      else 
       iocc0(ni,1) = 0
      end if
    else
      if (mod(icolin / len0, 2) .eq. 0) then
       iocc0(ni,1) = 0
      else
       iocc0(ni,1) = 1
      end if
    end if
  end do
  if (iocc0(iorbs-1,1) .eq. 1 .and. iocc0(iorbs,1) .eq. 1) then
    isgcol = -1
  else
    isgcol = 1
  end if
! swap occupation number for (iorb-1)th and (iorb)th orbitals
  nj = iocc0(iorbs-1,1)
  iocc0(iorbs-1,1) = iocc0(iorbs,1)
  iocc0(iorbs,1) = nj
! calculate col index corresponding to new occupation pattern
  nj = 1
  do ni=norbs,1,-1
    len0 = nrho / 2**(norbs - ni + 1)
    if (iocc0(ni,1) .eq. 1) nj = nj + len0
  end do
  icolout = nj
!
! sign
!
  isgnout = isgrow * isgcol    
!
else 
!
! row
!
  do ni=1,norbs
    len0 = nrho / 4**(norbs - ni + 1)
    if (mod(irowin, len0) .eq. 0) then
      if (mod(irowin / len0, 4) .eq. 1) then       ! vacant
        iocc0(ni,1) = 0
        iocc0(ni,2) = 0
      else if (mod(irowin / len0, 4) .eq. 2) then  ! spin up only
        iocc0(ni,1) = 1
        iocc0(ni,2) = 0
      else if (mod(irowin / len0, 4) .eq. 3) then  ! spin down only
        iocc0(ni,1) = 0
        iocc0(ni,2) = 1                         
      else                                         ! doubly occupied 
        iocc0(ni,1) = 1
        iocc0(ni,2) = 1
      end if
    else
      if (mod(irowin / len0, 4) .eq. 1) then
        iocc0(ni,1) = 1
        iocc0(ni,2) = 0
      else if (mod(irowin / len0, 4) .eq. 2) then
        iocc0(ni,1) = 0
        iocc0(ni,2) = 1
      else if (mod(irowin / len0, 4) .eq. 3) then
        iocc0(ni,1) = 1
        iocc0(ni,2) = 1
      else
        iocc0(ni,1) = 0
        iocc0(ni,2) = 0
      end if
    end if
  end do
!
!  if (iocc0(1,ispin) .eq. 1 .and. iocc0(iorbs,ispin) .eq. 1) then
  if (iocc0(iorbs-1,ispin) .eq. 1 .and. iocc0(iorbs,ispin) .eq. 1) then
    isgrow = -1
  else
    isgrow = 1
  end if
!
!  nj = iocc0(1,ispin)
!  iocc0(1,ispin) = iocc0(iorbs,ispin)
  nj = iocc0(iorbs-1,ispin)
  iocc0(iorbs-1,ispin) = iocc0(iorbs,ispin)
  iocc0(iorbs,ispin) = nj
!
  nj = 1
  do ni=norbs,1,-1
    len0 = nrho / 4**(norbs - ni + 1)
    if (iocc0(ni,1) .eq. 1 .and. iocc0(ni,2) .eq. 0) then
      nj = nj + len0
    else if (iocc0(ni,1) .eq. 0 .and. iocc0(ni,2) .eq. 1) then
      nj = nj + len0 * 2
    else if (iocc0(ni,1) .eq. 1 .and. iocc0(ni,2) .eq. 1) then
      nj = nj + len0 * 3
    end if
  end do
  irowout = nj
!
! column
!
  do ni=1,norbs
    len0 = nrho / 4**(norbs - ni + 1)
    if (mod(icolin, len0) .eq. 0) then
      if (mod(icolin / len0, 4) .eq. 1) then       ! vacant
        iocc0(ni,1) = 0
        iocc0(ni,2) = 0
      else if (mod(icolin / len0, 4) .eq. 2) then  ! spin up only
        iocc0(ni,1) = 1
        iocc0(ni,2) = 0
      else if (mod(icolin / len0, 4) .eq. 3) then  ! spin down only
        iocc0(ni,1) = 0
        iocc0(ni,2) = 1                         
      else                                         ! doubly occupied 
        iocc0(ni,1) = 1
        iocc0(ni,2) = 1
      end if
    else
      if (mod(icolin / len0, 4) .eq. 1) then
        iocc0(ni,1) = 1
        iocc0(ni,2) = 0
      else if (mod(icolin / len0, 4) .eq. 2) then
        iocc0(ni,1) = 0
        iocc0(ni,2) = 1
      else if (mod(icolin / len0, 4) .eq. 3) then
        iocc0(ni,1) = 1
        iocc0(ni,2) = 1
      else
        iocc0(ni,1) = 0
        iocc0(ni,2) = 0
      end if
    end if
  end do
!
!  if (iocc0(1,ispin) .eq. 1 .and. iocc0(iorbs,ispin) .eq. 1) then
  if (iocc0(iorbs-1,ispin) .eq. 1 .and. iocc0(iorbs,ispin) .eq. 1) then
    isgcol = -1
  else
    isgcol = 1
  end if
!
!  nj = iocc0(1,ispin)
!  iocc0(1,ispin) = iocc0(iorbs,ispin)
  nj = iocc0(iorbs-1,ispin)
  iocc0(iorbs-1,ispin) = iocc0(iorbs,ispin)
  iocc0(iorbs,ispin) = nj
!
  nj = 1
  do ni=norbs,1,-1
    len0 = nrho / 4**(norbs - ni + 1)
    if (iocc0(ni,1) .eq. 1 .and. iocc0(ni,2) .eq. 0) then
      nj = nj + len0
    else if (iocc0(ni,1) .eq. 0 .and. iocc0(ni,2) .eq. 1) then
      nj = nj + len0 * 2
    else if (iocc0(ni,1) .eq. 1 .and. iocc0(ni,2) .eq. 1) then
      nj = nj + len0 * 3
    end if
  end do
  icolout = nj
!
! sign
!
  isgnout = isgrow * isgcol    
!
end if
!
deallocate(iocc0, STAT=istat)
end subroutine swapms
