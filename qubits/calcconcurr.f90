!++++++entanglement_Concurrence of two-qubit+++++++
subroutine  calcconcurr(Ns,rho,concurrence)
   implicit none
   integer, parameter          :: matz=0     
   integer(kind=4), intent(in) :: Ns
   complex(kind=8), intent(in) :: rho(Ns,Ns)
   real(kind=8), intent(out)   :: concurrence
   complex(kind=8)             :: rho1(Ns,Ns),rho1cgj(Ns,Ns)
   complex(kind=8)             :: entangle_mat(Ns,Ns)
   integer(kind=4)             :: ierr                  
   real(kind=8)                :: ar(Ns,Ns),ai(Ns,Ns)                  
   real(kind=8)                :: wr(Ns)
   real(kind=8)                :: wi(Ns)
   real(kind=8)                :: xr(Ns,Ns)
   real(kind=8)                :: xi(Ns,Ns)
   real(kind=8)                :: sqrt_sum
   integer(kind=4)             :: i,j,l,m,l1,m1
   real(kind=8)                :: sigmayy1(4,4)

   sigmayy1=0.0d0
   sigmayy1(1,4)=-1.0d0
   sigmayy1(2,3)=1.0d0
   sigmayy1(3,2)=1.0d0
   sigmayy1(4,1)=-1.0d0

   do i=1, Ns
      do j=1, Ns
         rho1(i,j)=rho(i,j)
         rho1cgj(i,j)=dconjg(rho(i,j))
      end do
   end do

   entangle_mat=matmul(rho1,matmul(sigmayy1,matmul(rho1cgj,sigmayy1)))

   do l=1, Ns
      do m=1, Ns
         ar(l,m)=dble(entangle_mat(l,m))
         ai(l,m)=dimag(entangle_mat(l,m))
      end do
   end do

!! eigenvalues and eigenvectors of a complex general matrix.
!! matz=0: do not output eigenvectors		

   call cg(Ns,Ns,ar,ai,wr,wi,matz,xr,xi,ierr)     
!just use wr--the real parts of eigenvalues

!do l=1,Ns
!write(6,*)l, wr(l)
!end do
do l=1,Ns
   if (wr(l) .lt. 0.d0) then
      write(6,*)'calcconcurr: warning! negative eigenvalue found ', l, wr(l)
      wr(l) = 0.d0
   end if
end do

!! sort the eigenvalues
   call piksrt(Ns,wr)
   do l1=1, Ns
      wr(l1)=dsqrt(wr(l1))
   end do
!          
   sqrt_sum=0.0d0
   do m1=1, Ns-1
      sqrt_sum=sqrt_sum+wr(m1)
   end do
!  
   concurrence=max(0.0d0,wr(Ns)-sqrt_sum)

   return

end subroutine
!+++++++++++++++++++++++++++
