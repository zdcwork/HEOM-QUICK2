subroutine findsubvec(vec0, len0, vecinp, nlen, vec1, len1, vindex) 
implicit none
!
!  vec0 is an ordered vector of len0 (ascending order), and
!  vecinp is a vector containing nlen elements
!
!  The subroutine finds an ordered sub-vector, vec1, of len1, 
!  which consists of all the elements of vec0 that are not in vecinp
!  (len1 = len0 - nlen)
! 
!  vindex(ni) stores the position of vecinp(ni) in a sub-vector of vec0,
!  in which vecinp(1:ni-1) were removed
!
!  Warning! this code does not check whether vec0 is ordered
!
integer, intent(in)  :: len0, nlen, vec0(*), vecinp(*)
integer, intent(out) :: len1, vec1(*), vindex(*)
!
integer              :: ni, nj
integer              :: lentmp
!
if (len0 .le. 0 .or. nlen .le. 0 .or. len0 .lt. nlen) then
    write(6,*)
    write(6,*)'findsubvec: error! len0, nlen ', len0, nlen
    stop
end if
!
vec1(1:len0) = vec0(1:len0)
vindex(1:len0) = 0
len1 = len0 - nlen
!
do ni=1,nlen
!
! locate the index number for vecinp(ni)    
   lentmp = len0 - ni + 1
   call remove4list2(vec1, lentmp, vecinp(ni), vindex(ni))
   if (vindex(ni) .eq. 0) then
       write(6,*)
       write(6,*)'findsubvec: error! cannot find vecinp(', ni, ')', vecinp(ni)
       write(6,*)'vec1: ', (vec1(nj), nj=1,lentmp)
       stop
   end if
end do
!
return
end subroutine findsubvec
