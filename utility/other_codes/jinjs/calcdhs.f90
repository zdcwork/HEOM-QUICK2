subroutine calcdhs(tt)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8, intent(in) :: tt
integer :: ni, nj
real*8  :: dshift, del, der
real*8  :: dtmp1, dtmp2
!
dhs=czero
return
call getengyshift(tt, 1, 1, del)
if (nalf .eq. 2) then
  call getengyshift(tt, 2, 1, der)
end if
!
if (nalf .eq. 1) then
!  dshift = 5.d-1 * del  ! U(t) = 5.d-1 * V(t)
  dshift = 0.d0         ! U(t) = 0
!  dshift = del           ! U(t) = V(t)
else if (nalf .eq. 2) then
  if (norbs .eq. 1 .and. nspin .eq. 2) then
    if (dabs(epara4) .ge. dnano) then
      dshift = 5.d-1 * (del - der) * epara4
    else
      dshift = 5.d-1 * (del + der)
    end if
  else
    dshift = 5.d-1 * (del + der)
  end if
end if
!
! dH_s = dshift * ( \sum_{ms} c^\dag_{ms} c_{ms} )
!
dhs = czero
!
!if (nalf .eq. 1) return  ! one lead case ---> fix dot level energy
!
do nj=1,nspin
  do ni=1,norbs
    cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, ni, nj), 0.d0)
    call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(dshift, 0.d0), cmtmp1, nrho, &
               cmtmp1, nrho, cunity, dhs, nrho)
  end do
end do
!
! second-order correction : ec * N(t)**2
!
!if (nspin .ne. 1) then
!  if (capa .ge. dnano) then
!    dtmp1 = capa * dshift / qe0   ! N(t)
!    dtmp2 = ec * dtmp1**2
!    do ni=1,nrho 
!      dhs(ni,ni) = dhs(ni,ni) + dtmp2
!    end do 
!  end if
!end if
!
end subroutine calcdhs
