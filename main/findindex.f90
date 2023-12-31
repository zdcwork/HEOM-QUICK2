subroutine findindex(istart, iend, lrec, record, lwork, work, ipost, iop1, fact)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! purpose : find the exact position of a given record (or its complement)
!           in the indextable 
!           search from istart to iend in indextable
! input   
!                istart : starting position in indextable to search through
!                  iend : ending   position in indextable to search through
!                  lrec : length of record (number of balls)
! (-istart+iend+1)/lrec : number of records within the search range
!                record : the record 
!                  iop1 : 0, if search the input record
!                         1, if search the complement of the input record
!                 lwork : length of work space
!                  work : work space
!
! output
!                 ipost : exact position found if search succeeded
!                         points to the position PRIOR to the desired record, i.e., 
!                         if the record occupies 5th to 7th cells of indextable, then
!                         ipost=4.
!                         -1, if search failed
!                 fact  : 
!                          1, if even # of swaps made when sorting the record 
!                             (for iop1 = 0, fact is always equal to 1)
!                         -1, if odd ... 
! note 
!       the input record should follow the standard sequence upon entrance.
!
integer*8, intent(in)  :: istart, iend
integer, intent(in)    :: lrec, record(*), lwork
integer*1, intent(in)  :: iop1
integer, intent(inout) :: work(lwork)
integer*8, intent(out) :: ipost
integer*1, intent(out) :: fact
!
integer :: ni, nj, nk, igel, iswap
integer*8 :: low, high, lni, lnj
!
ipost = 0
if (lwork .lt. lrec) then
 write(6,*)
 write(6,*)' error! lwork < lrec in findindex ', lwork, lrec
 stop
end if
!
if (iop1 .eq. 0) then
  work(1:lrec) = record(1:lrec)
else if (iop1 .eq. 1) then
  do ni=1,lrec
    work(ni) = iredex(record(lrec-ni+1)) ! sequence is reversed for hermitian conjugate 
  end do
else
  write(6,*) 
  write(6,*)' error! unknown iop1 found in findindex ', iop1
end if
!
! sort before search if necessary
!
fact = 1
if (iop1 .eq. 1) then
  do ni=1,lrec-1
    do nj=ni+1,lrec
      if (work(ni) .gt. work(nj)) then
        nk       = work(ni)
        work(ni) = work(nj)
        work(nj) = nk
        fact     = fact * (-1)
      end if
    end do
  end do
end if
!
low  = 1
high = (iend - istart + 1) / lrec
lni  = (high + low) / 2
!
if (high .lt. low) then
 ipost = -1
 goto 100
end if
!
find: do while (low <= high)
        lnj = istart + (lni - 1) * lrec - 1
        call compareindex(lnj, lrec, work, igel)
!        write(6,*)' compare ', low, high, lni, igel
        if (low .eq. high .and. igel .ne. 0) then
         ipost = -1
         goto 100
        end if
        if (igel .gt. 0) then
          high = max(lni - 1, 1)
        else if (igel .lt. 0) then
          low  = min(lni + 1, (iend - istart + 1) / lrec)
        else
          ipost = lnj 
!          write(6,*)' found ', ipost
          goto 100
        end if
        lni = (high + low) / 2
      end do find
ipost = -1
100 continue
!-------debug
!write(6,*)'record, work'
!write(6,*)(record(ni),ni=1,lrec), '||', (work(nj),nj=1,lrec), '||', ipost
!write(6,*)'istart, iend, low, high, lni'
!write(6,*)istart, iend, low, high, lni
!call flush(6)
!-------end of debug
end subroutine findindex
!--------------------------------------------------------
subroutine compareindex(lni, nball, mindex, iout)
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer*8, intent(in) :: lni
integer,   intent(in) :: nball, mindex(*)
integer,   intent(out) :: iout
integer*8 :: lnj
integer   :: ni
!
iout = 0
do ni=1,nball
 if (indextable(lni + ni) .gt. mindex(ni)) then
  iout = 1
  return
 else if (indextable(lni + ni) .lt. mindex(ni)) then
  iout = -1
  return
 end if
end do
iout = 0
return
end subroutine compareindex
