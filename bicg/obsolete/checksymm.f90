subroutine checksymm
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer    :: itier, itype
integer    :: index_nk(MAXTIER)
integer*8  :: lni, lnj, lnk, lnl, lnm
integer    :: ni, nj, nk, ntmp1, ntmp2, ntmp3, ntmp4, dim1, istat
integer    :: i, j, k
real*8     :: cpu1, cpu2, dtmp1, dtmp2
complex*16 :: cgamma_nk
!
integer, allocatable :: icol(:)
!
if (igroundsteady .eq. 0) then
  cmtmp3(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho)
else if (igroundsteady .eq. 1) then
  call calcdhs(1.d20)
  cmtmp3(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho) + dhs(1:nrho, 1:nrho)
else
  write(6,*)
  write(6,*)' error! unknown igroundsteady in checksymm ', igroundsteady
  stop
end if
dim1 = nrho**2
!
! itier=1, nballs=0,  no balls, i.e., N = 0 in the derivation
!
itier = 1
!
do j=1,nrho
  do i=1,nrho
    dmtmp3(i,j) =  dble(-eye * cmtmp3(i,j))
    dmtmp4(i,j) = dimag(-eye * cmtmp3(i,j))
  end do
end do
call h2lmln(nrho, dmtmp3, nrho, dlmat1, dim1)
call h2lmln(nrho, dmtmp4, nrho, dlmat2, dim1)
call h2lmrn(nrho, dmtmp3, nrho, dlmat3, dim1)
call h2lmrn(nrho, dmtmp4, nrho, dlmat4, dim1)
dlmat1(1:dim1, 1:dim1) = dlmat1(1:dim1, 1:dim1) - dlmat3(1:dim1, 1:dim1)
dlmat2(1:dim1, 1:dim1) = dlmat2(1:dim1, 1:dim1) - dlmat4(1:dim1, 1:dim1)
!
write(6,*)
write(6,*)' Liouville matrix of -i[H, *] '
write(6,*)'        real part             '
call amatout(dim1, dim1, dlmat1)
write(6,*)'        imag part             '
call amatout(dim1, dim1, dlmat2)
call flush(6)
!
! check columns: both real and imaginary parts
!
allocate(icol(dim1), STAT=istat)
!
k = 0
do i=1,dim1
  dtmp1 = 0.d0
  dtmp2 = 0.d0
  do j=1,dim1
    dtmp1 = dtmp1 + dabs(dlmat1(j,i))
    dtmp2 = dtmp2 + dabs(dlmat2(j,i))
  end do
  if (dtmp1 .le. dnano .and. dtmp2 .le. dnano) then
    k = k + 1
    icol(k) = i
  end if
end do
ntmp1 = k
!
write(6,*)
write(6,*)' # of zero columns for Liouville matrix of -iH ', ntmp1
!write(6,*)(icol(k), k=1,ntmp1)
do i=1,ntmp1
  call jmod(icol(i), nrho, j, k)
  write(6,*)' row = ', j, ' col = ', k
end do
call flush(6)
!
deallocate(icol, STAT=istat)
!
end subroutine checksymm
