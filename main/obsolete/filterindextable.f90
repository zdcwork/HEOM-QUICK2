subroutine filterindextable(length0)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8, intent(inout) :: length0(*)
integer                  :: itier, i, j, k, istat, nballs, ncount, idraw
logical                  :: lcount
integer                  :: ntmp1, ntmp2
integer*8                :: lni, lnj, lnk
real*8                   :: cpu1, cpu2
!
if (ntier .le. 2) return
open(unit=12, file='indextable.tmp', form='unformatted', status='unknown')
rewind(12)
do lni=ifirst(1), ilast(2)
  write(12)indextable(lni)
end do
!
do itier=3,ntier
  call cpu_time(cpu1)
  nballs = itier - 1
  ncount = 0
  do lni=nfirst(itier), nlast(itier)
    lcount = .true.
    do i=1,nballs
      idraw     = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs + i )
      iatmp1(i) = mcor1(idraw)
      iatmp2(i) = idraw
    end do
    do i=1,nballs
      if (iatmp1(i) .le. ndrude) goto 100
    end do
    ntmp1 = 0
    do i=1,nballs-1
      do j=i+1,nballs
        ntmp2 = abs(iatmp1(i) - iatmp1(j))
        ntmp1 = max(ntmp1, ntmp2)
      end do
    end do
    if (ntmp1 .gt. nfilter) lcount = .false.
    100 continue
    if (lcount .or. itier .eq. 3) then
      ncount = ncount + 1
      do k=1,nballs
        write(12)iatmp2(k)
      end do
    end if
  end do
  call cpu_time(cpu2)
  write(6,*)
  write(6,*)' length of indextable at itier ', itier
  write(6,*)' actual ', ncount, ' estimated ', length0(itier)
  write(6,*)' cpu time elapsed: ', cpu2 - cpu1
  call flush(6)
  length0(itier) = ncount
end do
close(12)
!
end subroutine filterindextable
