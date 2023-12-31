subroutine ball2drawer_filter(itier, length0, length2)
use matmod
implicit none
!
! purpose : generate all possible cases for tier #itier 
!           with the preset filter applied, 
!           and write to array indextable consecutively.
!
include '../include/sizes'
include '../include/common'
!
integer,   intent(in)    :: itier
integer*8, intent(in)    :: length2(*)
integer*8, intent(inout) :: length0(*)
!
integer   :: i, j, k
integer   :: ni, nj, nk, nl, istat, nballs, idraw
integer   :: nmem, npick, npickmin
integer*8 :: ncount, ipost
integer*8 :: lni, lnj, lnk, lnl
integer*8 :: lena, lenb, lenc
integer*8 :: ilist, nlist
integer, allocatable :: ipoint(:), iwork(:), iveca(:), ivecb(:), ivecc(:)
integer, allocatable :: listncor(:,:)
logical   :: lrecord
integer*8            :: ifrac, ifrac2
external             :: ifrac, ifrac2
!
if (.not. lfilter) then
   write(6,*)'ball2drawer_filter: error entrance. lfilter = ', lfilter
   stop
end if
if (lscreen) then
    write(6,*)'ball2drawer_filter: error! lscreen=T found '
    write(6,*)'ball2drawer_filter: code unavailable yet   '
    stop
end if
!
! special treatment for itier=1 (no filter)
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
nmem   = min(nfilter_long, ncor)
if (nmem .ne. nfilter_long) then
   write(6,*)'ball2drawer_filter: warning! nfilter_long > ncor found '
   write(6,*)'ncor = ', ncor, ' nfilter_long = ', nfilter_long
   write(6,*)'nmem set to ncor = ', nmem
   call flush(6)
end if
!
npickmin = nfilter_count   ! npick goes from npickmin to nballs
if (npickmin > nballs) then
   write(6,*)'ball2drawer_filter: error! nballs < nfilter_count found '
   write(6,*)'nballs = ', nballs, ' nfilter_count = ', nfilter_count
   stop
end if
!
! calculate starting position in array indextable
!
ifirst(itier) = ilast(nballs) + 1                        ! first position of itier in indextable
nfirst(itier) = nlast(nballs) + 1 
!
noprfirst(itier) = noprlast(nballs) + 1
noprlast (itier) = noprlast(nballs) + length2(itier)     
!
ioprfirst(itier) = ioprlast(nballs) + 1
ioprlast (itier) = ioprlast(nballs) + length2(itier) * nballs
!
allocate(ipoint(nballs), iwork(nballs), STAT=istat)
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
lni = ioprfirst(itier)
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
! see whether its complement has already been recorded.
!
   if (ncount .gt. 0) then
      call findoprindex(ioprfirst(itier), lni-1, nballs, ipoint, nballs, iwork, ipost, 1)
   end if
!
! record
!
   if (ipost .lt. 0) then
      do nk=1,nballs
         oprtable(lni) = ipoint(nk)
         lni = lni + 1
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
!
ncount = 0
lnl    = ifirst(itier)
!
pick: do npick=npickmin,nballs
   lena = ifrac2(nballs, npick) / ifrac(nballs - npick)
   lenb = nmem**npick 
   lenc = 1
   allocate(iveca(lena*npick), ivecb(lenb*npick), STAT=istat)
   if (npick < nballs) then
      lenc = (ncor - nmem)**(nballs - npick)
      allocate(ivecc(lenc*(nballs-npick)), STAT=istat)
   end if
!
   call buildcombination(nballs, npick, iveca, lni)
   if (lni .ne. lena) then
      write(6,*)'ball2drawer_filter: error! at npick = ', npick
      write(6,*)'lena != lni ', lena, lni
      stop
   end if
!
   call assigndrawer(npick, nmem, ivecb, lnj)
   if (lnj .ne. lenb) then
      write(6,*)'ball2drawer_filter: error! at npick = ', npick
      write(6,*)'lenb != lnj ', lenb, lnj
      stop
   end if
!
   if (npick < nballs) then
      call assigndrawer(nballs-npick, ncor-nmem, ivecc, lnk)
      if (lnk .ne. lenc) then
         write(6,*)'ball2drawer_filter: error! at npick = ', npick
         write(6,*)'lenc != lnk ', lenc, lnk
         stop
      end if
   end if
!
   nlist = lena * lenb * lenc
   allocate(listncor(nballs,nlist), STAT=istat)
   listncor(1:nballs,1:nlist) = -1
   ilist = 0
   do lni=1,lena
      do lnj=1,lenb
         do lnk=1,lenc
            ilist = ilist + 1
            do ni=1,npick
               listncor(iveca((lni-1)*npick+ni),ilist) = kmats( ivecb((lnj-1)*npick+ni) )
            end do
            if (npick == nballs) cycle
            nj = 0
            do ni=1,nballs
               if (listncor(ni,ilist) > 0) cycle
               nj = nj + 1
               listncor(ni,ilist) = kmats( ivecc((lnk-1)*(nballs-npick)+nj) + nmem )
            end do
            if (nj .ne. nballs-npick) then
               write(6,*)'ball2drawer_filter: error! nballs-npick, nj ', nballs-npick, nj
               stop
            end if
         end do
      end do
   end do
!
   do lni=noprfirst(itier),noprlast(itier)
      lnj = ioprfirst(itier) + (lni - noprfirst(itier)) * nballs - 1
      iwork(1:nballs) = oprtable(lnj+1:lnj+nballs)
      do ilist=1,nlist
         ipoint(1:nballs) = listncor(1:nballs,ilist)
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
                    if (ipoint(nk) .ge. ipoint(nl)) then
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
               indextable(lnl) = idraw
               lnl = lnl + 1
            end do
            ncount = ncount + 1
         end if
      end do
   end do
!
   deallocate(iveca, ivecb, STAT=istat)
   deallocate(listncor, STAT=istat)
   if (npick < nballs) deallocate(ivecc, STAT=istat)
end do pick
!
length0(itier) = ncount
nlast(itier)   = nlast(nballs) + length0(itier)
ilast(itier)   = ilast(nballs) + length0(itier) * nballs
deallocate(ipoint, iwork, STAT=istat)
!
end subroutine ball2drawer_filter
