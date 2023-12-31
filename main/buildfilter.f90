subroutine buildfilter
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, istat, itier
integer                 :: ntmp1, ntmp2
integer                 :: index_nk(MAXTIER)
integer*8               :: lni, lnj
!
if (.not. lfilter) then
   write(6,*)'buildfilter: error entrance. lfilter = ', lfilter
   stop
end if
!
allocate(filtertable(nunk), STAT=istat)
write(6,*)'buildfilter: memory for filtertable (MB)  ', dble(nunk * 4) / dble(1024 * 1024)
call flush(6)
!
do lni=1,nunk
   filtertable(lni) = 0
end do
!
do itier=2,ntier
   do lni=nfirst(itier),nlast(itier)
      do ni=1,itier-1
         index_nk(ni) = mcor( indextable(ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni) )
      end do
      ntmp1 = 0 
      do ni=1,itier-1
         if (kpmats(index_nk(ni)) .le. nfilter_long) ntmp1 = ntmp1 + 1
      end do
      filtertable(lni) = ntmp1   
   end do
end do
!
end subroutine buildfilter
