subroutine sortindextable3
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: itier
real*8  :: cpu1, cpu2
!
write(6,*)
write(6,*)' entering sortindextable3 ...        '
!write(6,*)' this could be very time-consuming   '
call flush(6)
!
do itier=2,ntier0
   call cpu_time(cpu1)
   call quicksort_index(itier)
   call cpu_time(cpu2)
   write(6,*)' indextable on tier ', itier, ' is sorted by ', cpu2 - cpu1, ' secs'
   call flush(6)
end do
end subroutine sortindextable3
!
subroutine quicksort_index(itier)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
! quicksort code from Numerical Recipes in Fortran 90
!
integer, intent(in) :: itier
integer   :: ni, nj, nk, nl, nr, nn
integer   :: istat, nballs
logical   :: lswap
integer*8 :: lni, lnj, lnk, lnl, lnr, lnn
integer*8 :: ini, inj, ink, inl, inr, inn
integer   :: ndirect = 20 
integer   :: nstack, jstack
integer*8, allocatable :: istack(:)
integer   :: na(MAXTIER), nb(MAXTIER)
!
nballs = itier - 1
lnn    = nlast(itier) - nfirst(itier) + 1
nstack = MAXTIER * 2 * int(dlog(dble(lnn)) / dlog(2.d0)) 
jstack = 0
lnl    = nfirst(itier)
lnr    = nlast(itier)
!
allocate(istack(nstack), STAT=istat)
!
loopout: do
   if (lnr - lnl < ndirect) then   ! insertion sort when subarray small enough
      do lnj=lnl+1,lnr
         inj = ifirst(itier) - 1 + (lnj - nfirst(itier)) * nballs
         na(1:nballs) = indextable(inj+1:inj+nballs)
         do lni=lnj-1,lnl,-1
            ini = ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs
            call swapornot(nballs, indextable(ini+1), na(1), lswap)
            if (.not. lswap) exit
            indextable(ini+1+nballs:ini+nballs*2) = indextable(ini+1:ini+nballs)
         end do
         if (lswap) then
            lni = lnl - 1
            ini = ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs
         end if
         indextable(ini+1+nballs:ini+nballs*2) = na(1:nballs)
      end do     
      if (jstack == 0) then
         deallocate(istack, STAT=istat)
         return
      end if
      lnr = istack(jstack)
      lnl = istack(jstack - 1)
      jstack = jstack - 2
   else
      lnk = (lnl + lnr) / 2
      ink = ifirst(itier) - 1 + (lnk - nfirst(itier)) * nballs
      inl = ifirst(itier) - 1 + (lnl - nfirst(itier)) * nballs
      inr = ifirst(itier) - 1 + (lnr - nfirst(itier)) * nballs
! swap kth and (l+1)th
      na(1:nballs) = indextable(ink+1:ink+nballs)
      indextable(ink+1:ink+nballs) = indextable(inl+1+nballs:inl+nballs*2)
      indextable(inl+1+nballs:inl+nballs*2) = na(1:nballs)
! if arr(lnl) > arr(lnr), swap
      call swapornot(nballs, indextable(inl+1), indextable(inr+1), lswap)
      if (lswap) then
         na(1:nballs) = indextable(inl+1:inl+nballs)
         indextable(inl+1:inl+nballs) = indextable(inr+1:inr+nballs)
         indextable(inr+1:inr+nballs) = na(1:nballs)
      end if
! if arr(lnl+1) > arr(lnr), swap
      call swapornot(nballs, indextable(inl+1+nballs), indextable(inr+1), lswap)
      if (lswap) then
         na(1:nballs) = indextable(inl+1+nballs:inl+nballs*2)
         indextable(inl+1+nballs:inl+nballs*2) = indextable(inr+1:inr+nballs)
         indextable(inr+1:inr+nballs) = na(1:nballs)
      end if
! if arr(lnl) > arr(lnl+1), swap
      call swapornot(nballs, indextable(inl+1), indextable(inl+1+nballs), lswap)
      if (lswap) then
         na(1:nballs) = indextable(inl+1+nballs:inl+nballs*2)
         indextable(inl+1+nballs:inl+nballs*2) = indextable(inl+1:inl+nballs)
         indextable(inl+1:inl+nballs) = na(1:nballs)
      end if
!
      lni = lnl + 1
      ini = ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs
      lnj = lnr
      inj = ifirst(itier) - 1 + (lnj - nfirst(itier)) * nballs
!
      na(1:nballs) = indextable(ini+1:ini+nballs)
!
      do 
         do 
            lni = lni + 1
            ini = ifirst(itier) - 1 + (lni - nfirst(itier)) * nballs
            call swapornot(nballs, na(1), indextable(ini+1), lswap)   !  if (a > arr(i)) lswap == .true.
            if (.not. lswap) exit
         end do
         do
            lnj = lnj - 1
            inj = ifirst(itier) - 1 + (lnj - nfirst(itier)) * nballs
            call swapornot(nballs, indextable(inj+1), na(1), lswap)
            if (.not. lswap) exit
         end do
         if (lnj < lni) exit   
         nb(1:nballs) = indextable(ini+1:ini+nballs)
         indextable(ini+1:ini+nballs) = indextable(inj+1:inj+nballs)
         indextable(inj+1:inj+nballs) = nb(1:nballs)         
      end do
      inj = ifirst(itier) - 1 + (lnj - nfirst(itier)) * nballs
      inl = ifirst(itier) - 1 + (lnl - nfirst(itier)) * nballs
      indextable(inl+1+nballs:inl+nballs*2) = indextable(inj+1:inj+nballs)
      indextable(inj+1:inj+nballs) = na(1:nballs)
!      
      jstack = jstack + 2
      if (jstack > nstack) then
         write(6,*)
         write(6,*)'quicksort_index: ERROR! nstack too small! '
         stop
      end if
      if (lnr - lni + 1 >= lnj - 1) then 
         istack(jstack)     = lnr
         istack(jstack - 1) = lni
         lnr = lnj - 1
      else 
         istack(jstack)     = lnj - 1
         istack(jstack - 1) = lnl 
         lnl = lni
      end if
!
   end if
end do loopout
!
deallocate(istack, STAT=istat)
return
!
end subroutine quicksort_index
