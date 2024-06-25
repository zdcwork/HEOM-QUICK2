subroutine pre_atrx_bicg
use matmod
use auxmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, nl, iorbs, ispin, dim1, isgn
integer :: istat
!
dim1 = nrho**2
!
! new auxiliary matrices
!
allocate(lenams(norbs,nspin,nsgn,nlr,nnt), STAT=istat)
!
lenams = 0
do ispin=1,nspin
  do iorbs=1,norbs
    dmtmp1(1:nrho, 1:nrho) = amsall(1:nrho, 1:nrho, iorbs, ispin)
    call h2lmln(nrho, dmtmp1, nrho, dlmat1, dim1) 
    do nj=1,dim1
      do ni=1,dim1
        if (dabs(dlmat1(ni,nj)) .ge. dnano) then
          lenams(iorbs,ispin,1:nsgn,1,1) = lenams(iorbs,ispin,1:nsgn,1,1) + 1
        end if
      end do
    end do
    call h2lmrn(nrho, dmtmp1, nrho, dlmat1, dim1) 
    do nj=1,dim1
      do ni=1,dim1
        if (dabs(dlmat1(ni,nj)) .ge. dnano) then
          lenams(iorbs,ispin,1:nsgn,2,1) = lenams(iorbs,ispin,1:nsgn,2,1) + 1
        end if
      end do
    end do
    call h2lmlt(nrho, dmtmp1, nrho, dlmat1, dim1)
    do nj=1,dim1
      do ni=1,dim1
        if (dabs(dlmat1(ni,nj)) .ge. dnano) then
          lenams(iorbs,ispin,1:nsgn,1,2) = lenams(iorbs,ispin,1:nsgn,1,2) + 1
        end if
      end do
    end do
    call h2lmrt(nrho, dmtmp1, nrho, dlmat1, dim1)
    do nj=1,dim1
      do ni=1,dim1
        if (dabs(dlmat1(ni,nj)) .ge. dnano) then
          lenams(iorbs,ispin,1:nsgn,2,2) = lenams(iorbs,ispin,1:nsgn,2,2) + 1
        end if
      end do
    end do
  end do
end do
!
namsmax = 0
do ispin=1,nspin
  do iorbs=1,norbs
    do isgn=1,nsgn
      do ni=1,nlr
        do nj=1,nnt
          namsmax = max(namsmax, lenams(iorbs,ispin,isgn,ni,nj))
        end do
      end do
    end do
  end do
end do
write(6,*)
write(6,*)' maximal non-zero elements for lams ', namsmax
write(6,*)' nrho * nams                        ', nrho * nams
write(6,*)' # of all elements                  ', (nrho*nrho)**2
call flush(6)
!
allocate(lrowams(namsmax,norbs,nspin,nsgn,nlr,nnt), STAT=istat)
allocate(lcolams(namsmax,norbs,nspin,nsgn,nlr,nnt), STAT=istat)
allocate(lvalams(namsmax,norbs,nspin,nsgn,nlr,nnt), STAT=istat)
allocate(lrtmp1(namsmax), lctmp1(namsmax), lvtmp1(namsmax), STAT=istat)
!
write(6,*)
istat = 0
do ispin=1,nspin
  do iorbs=1,norbs
    do isgn=1,nsgn
      call h2laln(iorbs,ispin,isgn,nk)
      if (nk .ne. lenams(iorbs,ispin,isgn,1,1)) then
        istat = 1
        write(6,*)' error counting in h2laln ', iorbs, ispin, nk, lenams(iorbs,ispin,isgn,1,1)
      end if
      call h2larn(iorbs,ispin,isgn,nk)
      if (nk .ne. lenams(iorbs,ispin,isgn,2,1)) then
        istat = 1
        write(6,*)' error counting in h2larn ', iorbs, ispin, nk, lenams(iorbs,ispin,isgn,2,1)
      end if
      call h2lalt(iorbs,ispin,isgn,nk)
      if (nk .ne. lenams(iorbs,ispin,isgn,1,2)) then
        istat = 1
        write(6,*)' error counting in h2lalt ', iorbs, ispin, nk, lenams(iorbs,ispin,isgn,1,2)
      end if
      call h2lart(iorbs,ispin,isgn,nk)
      if (nk .ne. lenams(iorbs,ispin,isgn,2,2)) then
        istat = 1
        write(6,*)' error counting in h2lart ', iorbs, ispin, nk, lenams(iorbs,ispin,isgn,2,2)
      end if
    end do
  end do
end do
if (istat .ne. 0) then
  write(6,*)' error! check h2lams! '
  stop
else
  write(6,*)' h2lams completed succesfully '
  call flush(6)
end if
!
end subroutine pre_atrx_bicg
