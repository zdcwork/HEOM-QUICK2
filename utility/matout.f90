subroutine cmatout(m,n,cmat,dwork)
implicit none
!
integer, intent(in)    :: n, m
complex*16, intent(in) :: cmat(m,*)
real*8,  intent(inout) :: dwork(m,*)
!
dwork(1:m,1:n) = dble(cmat(1:m,1:n))
write(6,*)
write(6,*)' real part '
call amatout(m,n,dwork)
call flush(6)
!
dwork(1:m,1:n) = dimag(cmat(1:m,1:n))
write(6,*)
write(6,*)' imag part '
call amatout(m,n,dwork)
call flush(6)
!
end subroutine cmatout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Print asymmetric matrices                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine amatout(m,n,mat)
implicit none
!
integer, intent(in) :: n,m
real*8, intent(in)  :: mat(m,*)
integer :: i,ii,j
integer :: ir,row,range
integer, parameter :: numcol=5
real*8, dimension(numcol) :: s
!
do i=1,n,numcol
 ii=min0(i+numcol-1,n)
 write(6,1000) (ir,ir=i,ii)
 range=min0(n-i+1,numcol)
 do row=1,m
  do j=1,range
   s(j)=mat(row,i+j-1)
  enddo
  write(6,1001) row,(s(j),j=1,range)
 enddo
enddo
!
1000 format(5(11X,I3))
1001 format(I4,9D14.6)
!
return
end subroutine amatout
!
subroutine symout(n,mat, title)
!
implicit double precision(a-h,o-z)
!
integer, intent(in) :: n
real*8, intent(in)  :: mat(*)
character*(*) title
!
integer :: i,ii,j
integer :: ir,row,range
integer, parameter :: numcol=5
real*8, dimension(numcol) :: s
!
ind (i,j) = ((max0(i,j)-1)*max0(i,j))/2+min0(i,j)
!
write (6,'(/10X,(A))') title
do i=1,n,numcol
 ii=min0(i+numcol-1,n)
 write(6,1000) (ir,ir=i,ii)
 range=min0(n-i+1,numcol)
 do row=1,n
  do j=1,range
   s(j)=mat(ind(row,i+j-1))
  enddo
  write(6,1001) row,(s(j),j=1,range)
 enddo
enddo
!
1000 format(5(11X,I3))
1001 format(I4,9D14.6)
!
return
end subroutine symout
