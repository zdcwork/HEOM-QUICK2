subroutine prepareresi
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, istat, info, i, j, iorbs, ispin
integer :: lwork, imodms
complex*16, allocatable :: work(:), ctilde(:,:,:,:)
real*8,     allocatable :: rwork(:)
complex*16 :: ctmp1, ctmp2, ctmp3, ctmp4
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
do ispin=1,nspin
  do iorbs=1,norbs
    do nk=1,nsgn
      do nj=1,nalf
        cmtmp1 = czero
        cmtmp2 = czero
        do j=1,nrho
          do i=1,nrho
            do ni=1,nmcorr
              ctmp1 = gmcorr(ni,nj,nk)
              if (igroundsteady .eq. 1) then
                ctmp1 = ctmp1 + eye * dpm(nk) * engyshift(nj)
              end if
              ctmp2 = -cunity / (ctmp1 + eye * (valh(j) - valh(i))) 
              cmtmp1(i,j) = cmtmp1(i,j) + ctmp2 * cbcorr(ni,nj,nk)  ! => vb
              cmtmp2(i,j) = cmtmp2(i,j) + ctmp2 * cdcorr(ni,nj,nk)  ! => vd
            end do
          end do
        end do
        do j=1,nrho
          do i=1,nrho
            if (dpm(nk) .le. 0.d0) then
              cmtmp1(i,j) = cmtmp1(i,j) * ctilde(i,j, iorbs, ispin)
              cmtmp2(i,j) = cmtmp2(i,j) * ctilde(i,j, iorbs, ispin)
            else
              cmtmp1(i,j) = cmtmp1(i,j) * dconjg(ctilde(j,i, iorbs, ispin))
              cmtmp2(i,j) = cmtmp2(i,j) * dconjg(ctilde(j,i, iorbs, ispin))
            end if
          end do
        end do
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, vech, nrho, cmtmp1, nrho, &
                   czero, cmtmp3, nrho)
        call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp3, nrho, vech, nrho, &
                   czero, cmtmp1, nrho)
        call zgemm('n', 'n', nrho, nrho, nrho, cunity, vech, nrho, cmtmp2, nrho, &
                   czero, cmtmp3, nrho)
        call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp3, nrho, vech, nrho, &
                   czero, cmtmp2, nrho)
        vb(1:nrho, 1:nrho, nj, nk, iorbs, ispin) = cmtmp1(1:nrho, 1:nrho)
        vd(1:nrho, 1:nrho, nj, nk, iorbs, ispin) = cmtmp2(1:nrho, 1:nrho)
      end do
    end do
  end do
end do
!
deallocate(work, rwork, STAT=istat)
deallocate(ctilde, STAT=istat)
!
end subroutine prepareresi
