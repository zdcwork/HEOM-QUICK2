subroutine ball2drawerlite(itier, length0, length2)
use matmod
implicit none
!
! purpose : generate all possible cases for tier #itier,
!           and write to array indextable consecutively.
!
include '../include/sizes'
include '../include/common'
!
integer, intent(in)    :: itier
integer*8, intent(in)  :: length0(*), length2(*)
!
integer   :: i, j, k
integer   :: nj, nk, nl, istat, nballs, idraw
integer*8 :: ni, ncount, ipost
integer*8 :: lni, lnj, lnk
integer, allocatable :: ipoint(:), iwork(:)
!
! special treatment for itier=1
!
if (itier .eq. 1) then
  ifirst(1)     = length0(itier)
  ilast(1)      = length0(itier)
  nfirst(1)     = length0(itier)
  nlast(1)      = length0(itier)
  indextable(1) = 0              ! no drawer 
!----
  noprfirst(1)  = length2(itier)
  noprlast(1)   = length2(itier)
  ioprfirst(1)  = length2(itier)
  ioprlast(1)   = length2(itier)
  oprtable(1)   = 0              ! no operator
  return
end if
!
nballs = itier - 1
!
! calculate starting position in array indextable
!
ifirst(itier) = ilast(nballs) + 1                        ! first position of itier in indextable
ilast (itier) = ilast(nballs) + length0(itier) * nballs  ! last  ...
!
nfirst(itier) = nlast(nballs) + 1 
nlast (itier) = nlast(nballs) + length0(itier)
!
noprfirst(itier) = noprlast(nballs) + 1
noprlast (itier) = noprlast(nballs) + length2(itier)     
!
ioprfirst(itier) = ioprlast(nballs) + 1
ioprlast (itier) = ioprlast(nballs) + length2(itier) * nballs
!
allocate(ipoint(nballs), iwork(nballs), STAT=istat)
!
do i=1,nballs-1
  ipoint(i) = i
end do
ipoint(nballs) = nballs - 1
!
ncount = 0
!
write(6,*)' possible cases for '
write(6,*)' # of balls ', nballs, ' # of operators ', nopr
call flush(6)
!
ni = ioprfirst(itier)
ipost = -1
!
point: do 
!
! check pointer overflow
! 
   nj = nballs
   check: do
     ipoint(nj) = ipoint(nj) + 1
     if ( (ipoint(nj)+nballs-nj) .le. nopr ) then
       exit check
     else
       nj = nj - 1 
       if (nj >= 1) then
         cycle check
       else
         exit point
       endif
     endif
   end do check
!
   k = ipoint(nj)
   do i=nj+1,nballs
     k = k + 1
     ipoint(i) = k
   end do
!
! see whether its complement has already been recorded.
!
   if (ncount .gt. 0) then
     call findoprindex(ioprfirst(itier), ni-1, nballs, ipoint, nballs, iwork, ipost, 1)
   end if
!
! record
!
   if (ipost .lt. 0) then
     do nk=1,nballs
      oprtable(ni) = ipoint(nk)
      ni = ni + 1
     enddo
     ncount = ncount + 1
     if (length2(itier) .le. 1000) then
       write(6,100)ncount, (ipoint(nl), nl=1,nballs)
     endif
   end if
100 format(6x, i3, '| ', 4x, 20(2x, I3))
!
end do point
call flush(6)
!
! check
!
if (ncount .ne. length2(itier)) then
 write(6,*)
 write(6,*)' error building oprtable at tier ', itier
 write(6,*)ncount, length2(itier)
 stop
endif
!------------------
ncount = 0
ni = ifirst(itier)
!
do lni=noprfirst(itier),noprlast(itier)
  lnj = ioprfirst(itier) + (lni - noprfirst(itier)) * nballs - 1
  iwork(1:nballs)    = oprtable(lnj+1:lnj+nballs)
  ipoint(1:nballs-1) = 1
  ipoint(nballs)     = 0
  point2: do
    nj = nballs
    check2: do
      ipoint(nj) = ipoint(nj) + 1
      if ( ipoint(nj) .le. nvar0 ) then
        exit check2
      else
        nj = nj - 1
        if (nj >= 1) then
          cycle check2
        else
          exit point2
        end if
      end if
    end do check2
    ipoint(nj+1:nballs) = 1 
    do nk=1,nballs
      call look4drawer(jpm(iwork(nk)), jalf(iwork(nk)), jorbs(iwork(nk)), jspin(iwork(nk)), &
                       kmats(ipoint(nk)), idraw)
      indextable(ni) = idraw
      ni = ni + 1
    end do
    ncount = ncount + 1
  end do point2
end do
!
! check length of indextable 
! 
if (ncount .ne. length0(itier)) then
  write(6,*)
  write(6,*)' error building indextable at tier ', itier
  write(6,*)ncount, length0(itier)
  stop
end if
!
deallocate(ipoint, iwork, STAT=istat)
!
end subroutine ball2drawerlite
