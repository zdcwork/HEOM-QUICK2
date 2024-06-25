subroutine add2list(list0, dim0, len0, list1, dim1, ninput, nindex)
implicit none
!
! add an integer to an ordered list
!
integer, intent(in)  :: dim0, dim1, len0, ninput
integer, intent(in)  :: list0(dim0)
integer, intent(out) :: nindex, list1(dim1)
integer              :: mi, mj
!
mj = len0 + 1
if (len0 .gt. dim0 .or. mj .gt. dim1) then
    write(6,*)
    write(6,*)'add2list: error size ', dim0, dim1, len0
    stop
end if
!
compare0: do mi=1,len0
   if (list0(mi) .gt. ninput) then
      mj = mi
      exit compare0
   end if
end do compare0
list1(mj) = ninput
if (mj .gt. 1) then
   list1(1:mj-1) = list0(1:mj-1)
end if
if (mj .lt. len0+1) then
   list1(mj+1:len0+1) = list0(mj:len0)
end if
nindex = mj
return
end subroutine add2list
!
subroutine remove4list(list0, dim0, len0, list1, dim1, nindex)
implicit none
!
! remove an element from an ordered list
! 
integer, intent(in)  :: dim0, dim1, len0, nindex
integer, intent(in)  :: list0(dim0)
integer, intent(out) :: list1(dim1)
integer              :: mi, mj
!
if (nindex .le. 0 .or. nindex .gt. len0) then
   write(6,*)
   write(6,*)'remove4list: error dimension ', len0, nindex
   stop
end if
if (len0 .gt. dim0 .or. len0-1 .gt. dim1) then
   write(6,*)
   write(6,*)'remove4list: error size ', dim0, dim1, len0
   stop
end if
if (len0 .eq. 1) then
   write(6,*)'remove4list: the resultant list is empty '
   return
end if
mj = 0
do mi=1,len0
   if (mi .eq. nindex) cycle
   mj = mj + 1
   list1(mj) = list0(mi)
end do
return
end subroutine remove4list
!
subroutine remove4list2(list, nlen, ninput, nindex)
implicit none
!
! remove an element from an ordered list, and store
! in the same list 
!
! upon input the list has a length of nlen, and
! upon output the length of list changes to (nlen-1)
! 
! note: the code does not check whether the order of list
!
integer, intent(in)    :: nlen, ninput
integer, intent(out)   :: nindex
integer, intent(inout) :: list(*)
!
integer                :: ni
!
if (nlen .le. 0) then
    write(6,*)
    write(6,*)'remove4list2: error! nlen ', nlen
    stop
end if
!
nindex = 0
if (ninput .gt. list(nlen)) then
    return
end if
!
do ni=1,nlen
   if (list(ni) .eq. ninput) then
       nindex = ni
       exit
   end if
end do
if (nindex .ne. 0) then
    do ni=nindex,nlen-1
       list(ni) = list(ni+1)
    end do
end if
!
return
end subroutine remove4list2
