subroutine preparesteady(itrun, iresi)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in) :: itrun, iresi
!
integer :: ni, nj, nk, istat, info, i, j, iorbs, ispin
integer :: ialf, isgn, imats
integer :: dim1
integer :: lwork
complex*16, allocatable :: work(:), ctilde(:,:,:,:)
real*8,     allocatable :: rwork(:)
complex*16 :: ctmp1, ctmp2, ctmp3, ctmp4
!
write(6,*)
write(6,*)' entering preparesteady '
call flush(6)
!
dim1 = nrho**2
!
dhs = czero
if (igroundsteady .eq. 1) call calcdhs(1.d20)
cmtmp1 = hs + dhs
vech   = cmtmp1
!
lwork = nrho*nrho
allocate(work(lwork), rwork(3*nrho-2), STAT=istat)
allocate(ctilde(nrho,nrho,norbs,nspin), STAT=istat)
!
call zheev('V', 'u', nrho, vech, nrho, valh, work, lwork, rwork, info)
if (info .ne. 0) then
  write(6,*)
  write(6,*)' diagonalization fail in preparesteady ', info
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
if (itrun .eq. 1) then
  do isgn=1,nsgn
    do ispin=1,nspin
      do iorbs=1,norbs
        do ialf=1,nalf
          do imats=0,nmats,1
            ctmp1 = cgama(iorbs,ispin,imats+1,ialf,isgn)
            if (igroundsteady .eq. 1) then
              if (.not. leadspin) then
                if (isgn .eq. 1) then
                  ctmp1 = ctmp1 + eye * engyshift(ialf,1)
                else
                  ctmp1 = ctmp1 - eye * engyshift(ialf,1)
                end if
              else
                if (isgn .eq. 1) then
                  ctmp1 = ctmp1 + eye * engyshift(ialf,ispin)
                else
                  ctmp1 = ctmp1 - eye * engyshift(ialf,ispin)
                end if
              end if
            end if
            do j=1,nrho
              do i=1,nrho
                ctmp2 = -cunity / (ctmp1 + eye * (valh(j) - valh(i)))
                if (isgn .eq. 2) then
                  cmtmp2(i,j) = ctmp2 * ctilde(i,j,iorbs,ispin)
                else
                  cmtmp2(i,j) = ctmp2 * dconjg(ctilde(j,i,iorbs,ispin))
                end if
              end do
            end do 
            cpop(1:nrho,1:nrho,imats+1,ialf,iorbs,ispin,isgn) = cmtmp2(1:nrho,1:nrho) 
          end do
        end do
      end do
    end do
  end do
!
! make lpopr and lpopi (sum over all drawers)
!
  lpopr = 0.d0
  lpopi = 0.d0
  do isgn=1,nsgn
    do ispin=1,nspin
      do iorbs=1,norbs
        do ialf=1,nalf
          do imats=0,nmats,1
!
            cmtmp1(1:nrho, 1:nrho) = cpop(1:nrho,1:nrho,imats+1,ialf,iorbs,ispin,isgn) * &
                                     cb(iorbs,ispin,imats+1,ialf,isgn)
            cmtmp2(1:nrho, 1:nrho) = cpop(1:nrho,1:nrho,imats+1,ialf,iorbs,ispin,isgn) * &
                                     cd(iorbs,ispin,imats+1,ialf,isgn)
! (1) + ams_left * B_left 
            if (isgn .eq. 1) then
              call zamsmm('l', 'n', 'n', iorbs, ispin, cunity, cmtmp1, nrho, czero, cmtmp3, nrho)
            else
              call zamsmm('l', 'c', 'n', iorbs, ispin, cunity, cmtmp1, nrho, czero, cmtmp3, nrho)
            end if
            dmtmp1 =  dble(cmtmp3)
            dmtmp2 = dimag(cmtmp3)
            call h2lmln(nrho, dmtmp1, nrho, dlmat1, dim1)
            call h2lmln(nrho, dmtmp2, nrho, dlmat2, dim1)
            lpopr = lpopr + dlmat1
            lpopi = lpopi + dlmat2
! (2) - ams_left * D_right (Ai=0)
            do j=1,nrho
              do i=1,nrho
                if (isgn .eq. 1) then
                  dmtmp1(i,j) = amsall(i,j,iorbs,ispin)
                else
                  dmtmp1(i,j) = amsall(j,i,iorbs,ispin)
                end if
              end do
            end do
            dmtmp2 = 0.d0
            dmtmp3 =  dble(cmtmp2)
            dmtmp4 = dimag(cmtmp2)
            call h2lmlrn(nrho, dmtmp1, dmtmp3, nrho, dlmat3, dim1)  ! ArBr
            call h2lmlrn(nrho, dmtmp2, dmtmp4, nrho, dlmat4, dim1)  ! AiBi
            dlmat1 = dlmat3 - dlmat4                                ! ArBr - AiBi (real part)
            call h2lmlrn(nrho, dmtmp1, dmtmp4, nrho, dlmat5, dim1)  ! ArBi
            call h2lmlrn(nrho, dmtmp2, dmtmp3, nrho, dlmat6, dim1)  ! AiBr
            dlmat2 = dlmat5 + dlmat6                                ! ArBi + AiBr (imag part)
            lpopr  = lpopr - dlmat1                                 ! note the minus sign
            lpopi  = lpopi - dlmat2
! (3) - B_left * ams_right (Bi=0)
            do j=1,nrho
              do i=1,nrho
                if (isgn .eq. 1) then
                  dmtmp3(i,j) = amsall(i,j,iorbs,ispin)
                else
                  dmtmp3(i,j) = amsall(j,i,iorbs,ispin)
                end if
              end do
            end do
            dmtmp4 = 0.d0
            dmtmp1 =  dble(cmtmp1)
            dmtmp2 = dimag(cmtmp1)
            call h2lmlrn(nrho, dmtmp1, dmtmp3, nrho, dlmat3, dim1)  ! ArBr
            call h2lmlrn(nrho, dmtmp2, dmtmp4, nrho, dlmat4, dim1)  ! AiBi
            dlmat1 = dlmat3 - dlmat4                                ! ArBr - AiBi (real part)
            call h2lmlrn(nrho, dmtmp1, dmtmp4, nrho, dlmat5, dim1)  ! ArBi
            call h2lmlrn(nrho, dmtmp2, dmtmp3, nrho, dlmat6, dim1)  ! AiBr
            dlmat2 = dlmat5 + dlmat6                                ! ArBi + AiBr (imag part)
            lpopr  = lpopr - dlmat1                                 ! note the minus sign
            lpopi  = lpopi - dlmat2
! (4) + D_right * ams_right
            if (isgn .eq. 1) then
              call zamsmm('r', 'n', 'n', iorbs, ispin, cunity, cmtmp2, nrho, czero, cmtmp3, nrho)
            else
              call zamsmm('r', 'c', 'n', iorbs, ispin, cunity, cmtmp2, nrho, czero, cmtmp3, nrho)
            end if 
            dmtmp1 =  dble(cmtmp3)
            dmtmp2 = dimag(cmtmp3)
            call h2lmrn(nrho, dmtmp1, nrho, dlmat1, dim1)
            call h2lmrn(nrho, dmtmp2, nrho, dlmat2, dim1)
            lpopr = lpopr + dlmat1
            lpopi = lpopi + dlmat2
!
          end do
        end do
      end do
    end do
  end do
!
end if
!
deallocate(work, rwork, STAT=istat)
deallocate(ctilde, STAT=istat)
!
write(6,*)
write(6,*)' leaving preparesteady '
call flush(6)
!
end subroutine preparesteady
