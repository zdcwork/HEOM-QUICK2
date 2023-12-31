subroutine h2laln(korbs, kspin, ksgn, ncount)
use auxmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn
integer, intent(out) :: ncount
integer :: ni, nj, nk, nl, n1, n2
!
ncount = 0
do nk=1,nams
  do nl=1,nrho
    ncount = ncount + 1
    if (ksgn .eq. 2) then              ! ams
      lrowams(ncount,korbs,kspin,ksgn,1,1) = (nl - 1) * nrho + rowams(nk,korbs,kspin)
      lcolams(ncount,korbs,kspin,ksgn,1,1) = (nl - 1) * nrho + colams(nk,korbs,kspin)
    else if (ksgn .eq. 1) then         ! ams^t
      lrowams(ncount,korbs,kspin,ksgn,1,1) = (nl - 1) * nrho + colams(nk,korbs,kspin)
      lcolams(ncount,korbs,kspin,ksgn,1,1) = (nl - 1) * nrho + rowams(nk,korbs,kspin)
    else
      write(6,*)
      write(6,*)' unknown ksgn in h2laln', ksgn
      stop
    end if
    lvalams(ncount,korbs,kspin,ksgn,1,1) = valams(nk,korbs,kspin)
  end do
end do  
end subroutine h2laln
!
subroutine h2larn(korbs, kspin, ksgn, ncount)
use auxmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn
integer, intent(out) :: ncount
integer :: ni, nj, nk, nl, n1, n2
!
ncount = 0
do nk=1,nams
  do nl=1,nrho
    ncount = ncount + 1
    if (ksgn .eq. 2) then
      lrowams(ncount,korbs,kspin,ksgn,2,1) = (colams(nk,korbs,kspin) - 1) * nrho + nl
      lcolams(ncount,korbs,kspin,ksgn,2,1) = (rowams(nk,korbs,kspin) - 1) * nrho + nl
    else if (ksgn .eq. 1) then
      lrowams(ncount,korbs,kspin,ksgn,2,1) = (rowams(nk,korbs,kspin) - 1) * nrho + nl
      lcolams(ncount,korbs,kspin,ksgn,2,1) = (colams(nk,korbs,kspin) - 1) * nrho + nl
    else
      write(6,*)
      write(6,*)' unknown ksgn in h2larn', ksgn
      stop
    end if
    lvalams(ncount,korbs,kspin,ksgn,2,1) = valams(nk,korbs,kspin)
  end do
end do  
end subroutine h2larn
!
subroutine h2lalt(korbs, kspin, ksgn, ncount)
use auxmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn
integer, intent(out) :: ncount
integer :: ni, nj, nk, nl, n1, n2
!
ncount = 0
do nk=1,nams
  do nl=1,nrho
    ncount = ncount + 1
    if (ksgn .eq. 2) then
      lrowams(ncount,korbs,kspin,ksgn,1,2) = (nl - 1) * nrho + rowams(nk,korbs,kspin)
      lcolams(ncount,korbs,kspin,ksgn,1,2) = (colams(nk,korbs,kspin) - 1) * nrho + nl
    else if (ksgn .eq. 1) then
      lrowams(ncount,korbs,kspin,ksgn,1,2) = (nl - 1) * nrho + colams(nk,korbs,kspin)
      lcolams(ncount,korbs,kspin,ksgn,1,2) = (rowams(nk,korbs,kspin) - 1) * nrho + nl
    else
      write(6,*)
      write(6,*)' unknown ksgn in h2lalt', ksgn
      stop
    end if
    lvalams(ncount,korbs,kspin,ksgn,1,2) = valams(nk,korbs,kspin)
  end do
end do  
end subroutine h2lalt
!
subroutine h2lart(korbs, kspin, ksgn, ncount)
use auxmod
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: korbs, kspin, ksgn
integer, intent(out) :: ncount
integer :: ni, nj, nk, nl, n1, n2
!
ncount = 0
do nk=1,nams
  do nl=1,nrho
    ncount = ncount + 1
    if (ksgn .eq. 2) then
      lrowams(ncount,korbs,kspin,ksgn,2,2) = (colams(nk,korbs,kspin) - 1) * nrho + nl
      lcolams(ncount,korbs,kspin,ksgn,2,2) = (nl - 1) * nrho + rowams(nk,korbs,kspin)
    else if (ksgn .eq. 1) then
      lrowams(ncount,korbs,kspin,ksgn,2,2) = (rowams(nk,korbs,kspin) - 1) * nrho + nl
      lcolams(ncount,korbs,kspin,ksgn,2,2) = (nl - 1) * nrho + colams(nk,korbs,kspin)
    else
      write(6,*)
      write(6,*)' unknown ksgn in h2lart', ksgn
      stop
    end if
    lvalams(ncount,korbs,kspin,ksgn,2,2) = valams(nk,korbs,kspin)
  end do
end do  
end subroutine h2lart
