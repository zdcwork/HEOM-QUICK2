subroutine clearindextable(length0)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8, intent(inout) :: length0(*)
integer*8 :: lni, lnj, lnk
integer   :: itier, istat
integer   :: nballs, nmem, npick
integer, allocatable :: tmpvec(:)
integer*8              :: ifrac, ifrac2
external               :: ifrac, ifrac2
!
if (ntier0 .le. 2) return
!
lnk = 1
do itier=2,ntier0
  lnk = lnk + length0(itier) * (itier - 1)
end do
!
deallocate(indextable, STAT=istat)
allocate(indextable(lnk), STAT=istat)
!
do itier=2,ntier0
  nfirst(itier) = nlast(itier - 1) + 1
  nlast (itier) = nlast(itier - 1) + length0(itier)
  ifirst(itier) = ilast(itier - 1) + 1
  ilast (itier) = ilast(itier - 1) + length0(itier) * (itier - 1)
end do
!
!open(unit=12, file='indextable.tmp', form='unformatted', status='unknown')
open(unit=12, file='indextable.tmp', form='binary', status='unknown')
rewind(12)
do lni=1,lnk 
  read(12)indextable(lni)
end do
close(unit=12, status="delete")
!

!open(99, file='index1.dat')
!rewind(99)
!itier=4
!do lni=nfirst(itier),nlast(itier)
!   lnj = ifirst(itier) - 1  + (lni - nfirst(itier)) * (itier - 1)
!   write(99,999)lni, (indextable(lnj+istat), istat=1,itier-1)
!end do
!close(99)
!999 format(I8, 2x, 10(I2, 1x))
!stop

!nballs = 3
!nmem   = 5
!lni    = nmem**nballs 
!lnk    = lni * nballs
!allocate(tmpvec(lnk), STAT=istat)
!write(6,*)
!write(6,*)'clearindextable: debugging '
!call assigndrawer(nballs, nmem, tmpvec, lnj)
!write(6,*)'clearindextable: lni ', lni, ' lnj ', lnj
!call flush(6)
!deallocate(tmpvec, STAT=istat)
!stop

!nballs = 5
!npick  = 3
!lni    = ifrac2(nballs,npick)/ifrac(nballs-npick)
!lnk    = lni * npick
!allocate(tmpvec(lnk), STAT=istat)
!write(6,*)
!write(6,*)'clearindextable: debugging '
!call buildcombination(nballs, npick, tmpvec, lnj)
!write(6,*)'clearindextable: lni ', lni, ' lnj ', lnj
!call flush(6)
!deallocate(tmpvec, STAT=istat)
!stop

!
end subroutine clearindextable
