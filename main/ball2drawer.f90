subroutine ball2drawer(itier, length0, length2)
use matmod
implicit none
!
! purpose : generate all possible cases for tier #itier,
!           and write to array indextable consecutively.
!
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)    :: itier
integer*8, intent(inout) :: length2(*)
integer*8, intent(inout) :: length0(*)
!
integer   :: i, j, k
integer   :: ncor0
integer   :: nj, nk, nl, istat, nballs, idraw, nlist, ilist
integer*8 :: ni, ncount, ipost
integer*8 :: lni, lnj, lnk
real*8    :: dtmp1, dtmp2
integer   :: nvec1(maxtier)
integer, allocatable :: ipoint(:), iwork(:)
integer, allocatable :: listncor(:,:)
logical   :: lrecord
!
if (itier .le. ntier) then
    ncor0 = ncor
else
    ncor0 = ncor_slow
    if (lad_fast) then
        ncor0 = ncor_fast
    end if
end if
if (lomit .and. itier .gt. 2) then
    dtmp1 = dcbmax**(itier-1) * dfact_omit
end if
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
nlist  = ncor0**nballs
!
! calculate starting position in array indextable
!
ifirst(itier) = ilast(nballs) + 1                        ! first position of itier in indextable
nfirst(itier) = nlast(nballs) + 1 
!
noprfirst(itier) = noprlast(nballs) + 1
!noprlast (itier) = noprlast(nballs) + length2(itier)     
!
ioprfirst(itier) = ioprlast(nballs) + 1
!ioprlast (itier) = ioprlast(nballs) + length2(itier) * nballs
!
allocate(ipoint(nballs), iwork(nballs), STAT=istat)
allocate(listncor(nballs,nlist), STAT=istat)
!
ipoint(1:nballs - 1) = 1
ipoint(nballs)       = 0
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
     if ( (ipoint(nj)) .le. nopr ) then
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
   ipoint(nj+1:nballs) = ipoint(nj)
!
   if (lscreen) then
       do nk=1,nballs
          if (jomit(ipoint(nk)) .eq. 1) then
              cycle point
          end if
       end do
   end if
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
!     if (length2(itier) .le. 1000) then
     if (length2(itier) .le. 100) then
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
if (.not. lscreen .and. ncount .ne. length2(itier)) then
    write(6,*)
    write(6,*)' error building oprtable at tier ', itier
    write(6,*)ncount, length2(itier)
    stop
endif
length2(itier)  = ncount
noprlast(itier) = noprlast(nballs) + length2(itier)     
ioprlast(itier) = ioprlast(nballs) + length2(itier) * nballs
!
ilist = 0
ipoint(1:nballs-1) = 1
ipoint(nballs)     = 0
point1: do
  nj = nballs
  check1: do
    ipoint(nj) = ipoint(nj) + 1
    if ( ipoint(nj) .le. ncor0 ) then
      exit check1
    else
      nj = nj - 1
      if (nj >= 1) then
        cycle check1
      else
        exit point1
      end if
    end if
  end do check1
  ipoint(nj+1:nballs) = 1 
  ilist = ilist + 1
  listncor(1:nballs,ilist) = ipoint(1:nballs)
end do point1
if (ilist .ne. nlist) then
  write(6,*)
  write(6,*)' error! error in <ball2drawer>! ilist, nlist : ', ilist, nlist
  stop
end if
!
!write(6,*)
!write(6,*)' arrangement of ipoint : ', nballs
!do nk=1,nlist
!  write(6,814)nk, (listncor(nl,nk), nl=1,nballs)
!end do
!814 format(I6, 2x, 4(I4,2x))
!call flush(6)
!
!------------------
ncount = 0
ni = ifirst(itier)
do lni=noprfirst(itier),noprlast(itier)
!do lni=noprfirst(itier),noprfirst(itier)+length2(itier)-1
  lnj = ioprfirst(itier) + (lni - noprfirst(itier)) * nballs - 1
  iwork(1:nballs) = oprtable(lnj+1:lnj+nballs)
  do ilist=1,nlist
     if (itier .le. ntier) then
         ipoint(1:nballs) = listncor(1:nballs,ilist)
     else
         if (lad_fast) then
             do nl=1,nballs
                ipoint(nl) = icor_fast(listncor(nl,ilist))
             end do
         else
             do nl=1,nballs
                ipoint(nl) = icor_slow(listncor(nl,ilist))
             end do
         end if
     end if
     lrecord = .true.
!
     check2: do nk=1,nballs-1
       do nl=nk+1,nballs
          if (iwork(nk) .eq. iwork(nl)) then
             if (lsdraw) then
                if (ipoint(nk) .gt. ipoint(nl)) then
                   lrecord = .false.
                   exit check2  
                end if
             else 
                if (ipoint(nk) .ge. ipoint(nl)) then   ! exclude same drawer (same operator and same memory component)
                   lrecord = .false.
                   exit check2  
                end if
             end if
          end if
       end do
     end do check2
!
     if (lrecord) then
         do nk=1,nballs
            call look4drawer(jpm(iwork(nk)), jalf(iwork(nk)), jorbs(iwork(nk)), jspin(iwork(nk)),  &
                             ipoint(nk), idraw)
            nvec1(nk) = idraw
         end do
         if (lomit .and. itier .gt. 2) then
             dtmp2 = 1.d0
             do nk=1,nballs
                idraw = nvec1(nk)
                dtmp2 = dtmp2 * cdabs(cb(morbs(idraw),mspin(idraw),mcor(idraw),ilead(idraw),mpm(idraw)))
             end do
             if (dtmp2 .gt. dtmp1) then
                 do nk=1,nballs
                    indextable(ni) = nvec1(nk)
                    ni = ni + 1
                 end do
                 ncount = ncount + 1
             end if
         else
             do nk=1,nballs
                indextable(ni) = nvec1(nk)
                ni = ni + 1
             end do
             ncount = ncount + 1
         end if
     end if
  end do
end do
length0(itier) = ncount
nlast(itier)   = nlast(nballs) + length0(itier)
ilast(itier)   = ilast(nballs) + length0(itier) * nballs
deallocate(ipoint, iwork, listncor, STAT=istat)
!
end subroutine ball2drawer
