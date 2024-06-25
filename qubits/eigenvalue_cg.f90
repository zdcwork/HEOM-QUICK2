subroutine eigenvalue_cg(Ns,rho,lambda)
   implicit none
   integer(kind=4), intent(in)  :: Ns
   complex(kind=8), intent(in)  :: rho(Ns,Ns)
   complex(kind=8), intent(out) :: lambda(Ns)
   integer(kind=4), parameter   :: matz=0

   integer(kind=4)              :: ierr                  
   real(kind=8)                 :: ar(Ns,Ns),ai(Ns,Ns)	                  
   real(kind=8)                 :: wr(Ns)
   real(kind=8)                 :: wi(Ns)
   real(kind=8)                 :: xr(Ns,Ns)
   real(kind=8)                 :: xi(Ns,Ns)
   integer(kind=4)              :: i,j
		
   do i=1, Ns
      do j=1, Ns
         ar(i,j)=dble(rho(i,j))
         ai(i,j)=dimag(rho(i,j))
      end do
   end do		
		
   call cg(Ns,Ns,ar,ai,wr,wi,matz,xr,xi,ierr)
   !call piksrt(Ns,wr)

   do i=1, Ns
	  lambda(i)=wr(i)+(0,1)*wi(i)
   end do

   return
	  
end subroutine
!+++++++++++++++++++++++++++++++