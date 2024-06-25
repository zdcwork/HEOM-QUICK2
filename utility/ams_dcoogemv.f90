subroutine ams_dcoogemv(transa, ldim, lval, lrow, lcol, length, dvin, dvout)
!
! Calculate A * Vin = Vout,
! A is a ldim * ldim sparse matrix with 'length' nonzero elements, 
! dvin and dvout are input and output vectors of dimension ldim.
! By using the sparsity of creation/annihilation operator matrix 
! (in Liouville space), the computational cost reduces from 
! ldim * ldim to length
!
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent (in) :: transa
integer, intent (in)     :: ldim, length
integer, intent (in)     :: lrow(*), lcol(*)
real*8, intent (in)      :: lval(*), dvin(*)
real*8, intent (out)     :: dvout(*)
!
integer :: ni, nj, nk
logical :: lsame
external :: lsame
logical :: nota, conja
!
nota  = lsame(transa, 'N')
conja = lsame(transa, 'T') 
!
if ( (.not. nota) .and. (.not. conja) ) then
   write(6,*)'ams_dcoogemv: error! nota = ', nota, ' conja = ', conja
   stop
end if
!
dvout(1:ldim) = 0.d0
!
if (nota) then
   do nk=1,length
      ni = lrow(nk)
      nj = lcol(nk)
      dvout(ni) = dvout(ni) + lval(nk) * dvin(nj)
   end do
else 
   do nk=1,length
      ni = lcol(nk)
      nj = lrow(nk)
      dvout(ni) = dvout(ni) + lval(nk) * dvin(nj)
   end do
end if
!
end subroutine ams_dcoogemv
