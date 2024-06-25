subroutine piksrt(n,arr)
   integer, parameter    :: dp = SELECTED_REAL_KIND(12, 60)
   integer               :: n
   real(dp)              :: arr(n)
   integer               :: i,j
   real(dp)              :: a

   do j=2, n
      a=arr(j)
      do i=j-1,1,-1
         if (arr(i).le.a) go to 10
         arr(i+1)=arr(i)
      end do
      i=0
10    arr(i+1)=a
   end do

   return

end subroutine
!++++++++++++++++++++++++++++
