subroutine examindextable(length0)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8, intent(inout) :: length0(*)
integer :: itier, i, j, k, istat, nballs, ncount, nscreen
integer*8 :: lni, lnj, lnk
real*8  :: cpu1, cpu2
integer, allocatable :: vec1(:)
!
if (ntier0 .le. 2) return
!
allocate(vec1(ntier0-1), STAT=istat)
!
!open(unit=12, file='indextable.tmp', form='unformatted', status='unknown')
open(unit=12, file='indextable.tmp', form='binary', status='unknown')
rewind(12)
do lni=ifirst(1), ilast(2)
  write(12)indextable(lni)
end do
!
nscreen = 0
do itier=3,ntier0
  call cpu_time(cpu1)
  nballs = itier - 1
  ncount = 0
  do lni=nfirst(itier), nlast(itier)
    do i=1,nballs
      iatmp1(i) = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs + i )
    end do
    istat = 0
! screen out the cases where `opr' are same while `ncor' are not in ascending order
    exam: do k=1,nballs-1
      if (iatmp1(k+1) .lt. iatmp1(k)) then    ! i <= j <= k <= ... is admitted  
        istat = -1
        exit exam
      end if
    end do exam
!
    if (lsimple) then  ! screen out all ADOs with duplicate 'opr'
       if (lwalf) then
          do i=1,nballs
             vec1(i) = mopr(iatmp1(i))
          end do 
          exam2: do j=1,nballs
             do i=1,j-1  ! (execute (j-1)-1 times)
                !if (i .lt. 0 .or. i .ge. j) cycle exam2
                if (vec1(i) .eq. vec1(j)) then
                   istat = -1
                   nscreen = nscreen + 1
                   exit exam2
                end if
             end do
          end do exam2
       else
          do i=1,nballs
             vec1(i) = msopr(iatmp1(i))
          end do 
          exam3: do j=1,nballs
             do i=1,j-1
                ! if (i .lt. 0 .or. i .ge. j) cycle exam3
                if (vec1(i) .eq. vec1(j)) then
                   istat = -1
                   nscreen = nscreen + 1
                   exit exam3
                end if
             end do
          end do exam3
       end if
    end if 
!
    if (istat .eq. 0) then
      ncount = ncount + 1
      do k=1,nballs
        write(12)iatmp1(k)
      end do
    end if
  end do
  call cpu_time(cpu2)
  write(6,*)
  write(6,*)' length of indextable at itier ', itier
  write(6,*)' actual ', ncount, ' estimated ', length0(itier)
!
  if (lsimple) then
     write(6,*)'examindextable: nscreen = ', nscreen
  end if
  write(6,*)' cpu time elapsed: ', cpu2 - cpu1
  call flush(6)
  length0(itier) = ncount
end do
close(12)
!
deallocate(vec1, STAT=istat)
return
!
end subroutine examindextable
