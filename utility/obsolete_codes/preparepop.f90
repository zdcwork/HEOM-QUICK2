subroutine preparepop
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, istat, info, i, j
integer :: lwork, imodms
complex*16, allocatable :: work(:), ctilde(:,:,:,:)
real*8,     allocatable :: rwork(:)
complex*16 :: ctmp1, ctmp2
external   :: imodms
!
if (igroundsteady .eq. 0) then
  dhs(1:nrho, 1:nrho) = czero
else if (igroundsteady .eq. 1) then
  call calcdhs(1.d20)
end if
cmtmp1(1:nrho, 1:nrho) = hs(1:nrho, 1:nrho) + dhs(1:nrho, 1:nrho)
vech(1:nrho, 1:nrho)   = cmtmp1(1:nrho, 1:nrho)
!
lwork = nrho*nrho
allocate(work(lwork), rwork(3*nrho-2), STAT=istat)
allocate(ctilde(nrho,nrho,norbs,nspin), STAT=istat)
!
call zheev('V', 'u', nrho, vech, nrho, valh, work, lwork, rwork, info)
if (info .ne. 0) then
  write(6,*)
  write(6,*)' diagonalization fail in preparepop ', info
  stop
end if
!
do nj=1,nrho
  do ni=1,nrho
    vinvh(ni,nj) = dconjg(vech(nj,ni))
  end do
end do
!
do ni=1,norbs
  do nj=1,nspin
    cmtmp2(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, ni, nj), 0.d0)
    call zgemm('c', 'n', nrho, nrho, nrho, cunity, vech, nrho, cmtmp2, nrho, &
               czero, cmtmp3, nrho)
    call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, vech, nrho, &
               czero, cmtmp2, nrho)
    ctilde(1:nrho, 1:nrho, ni, nj) = cmtmp2(1:nrho, 1:nrho)
  end do
end do
!
do ni=1,nvar
  ctmp1 = cgamma(imodms(ni)) 
  if (igroundsteady .eq. 1) then
    ctmp1 = ctmp1 + eye * dipm(ni) * engyshift(ilead(ni))
  end if
  do j=1,nrho
    do i=1,nrho
      ctmp2 = -cunity / (ctmp1 + eye * (valh(j) - valh(i)))
      if (dipm(ni) .le. 0.d0) then
        cmtmp2(i,j) = ctmp2 * ctilde(i,j, morbs(ni), mspin(ni))
      else
        cmtmp2(i,j) = ctmp2 * dconjg(ctilde(j,i, morbs(ni), mspin(ni)))
      end if
    end do
  end do
  call zgemm('n', 'n', nrho, nrho, nrho, cunity, vech, nrho, cmtmp2, nrho, &
             czero, cmtmp3, nrho)
  call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp3, nrho, vech, nrho, &
             czero, cmtmp2, nrho)
  vpop(1:nrho, 1:nrho, ni) = cmtmp2(1:nrho, 1:nrho)
end do
!
deallocate(work, rwork, STAT=istat)
deallocate(ctilde, STAT=istat)
!
end subroutine preparepop
