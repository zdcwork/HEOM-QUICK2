subroutine checkhermicity(cin, nmat, cwork, ndim, lherm, derr)
implicit none
!
integer,    intent(in)    :: nmat, ndim
complex*16, intent(in)    :: cin(nmat,*)
complex*16, intent(inout) :: cwork(ndim,*)
logical,    intent(out)   :: lherm
real*8,     intent(out)   :: derr
integer :: ni, nj
real*8  :: dtmp1
real*8, parameter :: dpico = 1.d-12
!
if (ndim < nmat) then
   write(6,*)'checkhermicity: error! ndim < nmat ', ndim, nmat
   stop
end if
!
lherm = .true.
do nj=1,nmat
   do ni=1,nmat
      cwork(ni,nj) = cin(ni,nj) - dconjg(cin(nj,ni))
   end do
end do
call cmaxmat(nmat, nmat, cwork, nmat, derr)
if (derr > dpico) then
    lherm = .false.
end if
!
return
end subroutine checkhermicity
