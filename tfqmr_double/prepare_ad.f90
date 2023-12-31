subroutine prepare_ad
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                   :: nrho2, info, istat, nwork, ni, nj, ncount
integer                   :: iorbs, iorbsp, ispin, ialf, isgn, icor, idraw, isopr, im
integer                   :: ntmp1, ntmp2, ntmp3, mi, mj, iorbs1, iorbs2, ij, ki
integer                   :: jsgn, jsgnb, jorbs1, jorbs2, jjspin, jjalf, jm, jcor, jmats, jdraw, jsopr, jjorbs
character*1               :: transa, transb
integer,    allocatable   :: ipiv(:)
real*8,     allocatable   :: rwork(:), tmpval(:), tmpval2(:)
complex*16, allocatable   :: zwork(:,:), zwork2(:,:)
complex*16, allocatable   :: cmat1(:,:), cmat2(:,:)
complex*16, allocatable   :: zlsys(:,:), zlself(:,:), zself1(:,:), zself2(:,:) 
complex*16, allocatable   :: zlu(:,:), zlu1(:,:), zlu2(:,:)
complex*16, allocatable   :: zlv1(:), zlv2(:), zlv3(:), zlv4(:), zlv5(:), zlv6(:)
logical                   :: lherm
real*8                    :: dcut
real*8                    :: dtmp0, dtmp1, dtmp2, dtmp3
complex*16                :: zcoefl, zcoefr, zgama, zcb, zcd, zw, zgamaj
complex*16                :: ctmp0, ctmp1, ctmp2, ctmp3
!
if (.not. lad) then
    write(6,*)
    write(6,*)'prepare_ad: error entrance ', lad
    stop
end if
nrho2 = nrho**2
nwork = nrho2**2
!
!dcut = dsmall
dcut = dfemto
write(6,*)'prepare_ad: cutoff adopted = ', dcut
call flush(6)
!
allocate(cmat1(nrho,nrho), cmat2(nrho,nrho), STAT=istat)
allocate(zlsys(nrho2,nrho2), STAT=istat)
if (lscba_ad) then
    !if (offcor) then
    !    write(6,*)
    !    write(6,*)'prepare_ad: offcor=T found for lscba_ad=T '
    !    write(6,*)'prepare_ad: code not ready yet '
    !    stop
    !end if
    allocate(zlself(nrho2,nrho2), zlu(nrho2,nrho2), STAT=istat)
    allocate(zself1(nrho2,nrho2), zlu1(nrho2,nrho2), STAT=istat)
    allocate(zself2(nrho2,nrho2), zlu2(nrho2,nrho2), STAT=istat)
    allocate(zlv1(nrho2), zlv2(nrho2), zlv3(nrho2), zlv4(nrho2), zlv5(nrho2), zlv6(nrho2), STAT=istat)
end if
!
if (igroundsteady .eq. 0) then
    cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
else if (igroundsteady .eq. 1) then
    cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
else 
!    write(6,*)
!    write(6,*)'prepare_ad: error! wrong igroundsteady ', igroundsteady
!    stop
    cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
    if (ltd2st) then
        write(6,*)
        write(6,*)'prepare_ad: ltd2st=T found, use chemical potential shift '
        write(6,*)'            at t=infty for time-dependent phase of C(t)  '
    end if
end if
!
if (.not. allocated(zmatder)) allocate(zmatder(nrho2,nrho2), STAT=istat)
if (.not. allocated(dvalder)) allocate(dvalder(nrho2), STAT=istat)
!
! diagonalize L
! 
call zhil2liou('l', 'n', nrho, cmat1, nrho, zlmtmp1, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zlmtmp2, nrho2)
zlsys(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
zmatder(1:nrho2,1:nrho2) = zlsys(1:nrho2,1:nrho2)
call checkhermicity(zmatder, nrho2, zlmtmp1, nrho2, lherm, dtmp1)
if (.not. lherm) then
    write(6,*)'prepare_ad: error! L is not Hermitian '
    stop
end if
!
allocate(ipiv(nrho2), STAT=istat)
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
allocate(zwork2(nrho2,nrho2), tmpval2(nrho2), STAT=istat)
call zheev('V', 'u', nrho2, zmatder, nrho2, dvalder, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
    write(6,*)'prepare_ad: error! diagonalization failed 1 ', info
    stop
end if
zwork2(1:nrho2,1:nrho2) = zmatder(1:nrho2,1:nrho2)
tmpval2(1:nrho2) = dvalder(1:nrho2)
!
!write(6,*)'prepare_ad: check zero eigenvalues and eigenvectors of L '
!do nj=1,nrho2
!   dtmp0 = 0.d0
!   do ni=1,nrho2
!      dtmp0 = dtmp0 + cdabs(zwork2(ni,nj))**2 
!   end do
!   write(6,*)nj, tmpval2(nj), dtmp0
!   call flush(6)
!end do
!stop
!
! sort and screen out eigenvectors with zero eigenvalues
!
zwork = czero
ncount = 0
do ni=1,nrho2
   if (dabs(dvalder(ni)) .le. dpico) cycle
   ncount = ncount + 1
   tmpval(ncount) = dvalder(ni)
   zwork(1:nrho2,ncount) = zmatder(1:nrho2,ni)
end do
write(6,*)'prepare_ad: ', ncount, ' out of ', nrho2, ' eigenvalues are nonzero '
call flush(6)
zmatder(1:nrho2,1:nrho2) = czero
zmatder(1:nrho2,1:ncount) = zwork(1:nrho2,1:ncount)
dvalder(1:nrho2) = 0.d0
dvalder(1:ncount) = tmpval(1:ncount)
nnzeig = ncount
!
! calculate the approximated generating functional for each drawer
! (both left and right coefficient matrices)
!
if (.not. allocated(zgefnl)) allocate(zgefnl(nrho,nrho,nvar), STAT=istat)
if (.not. allocated(zgefnr)) allocate(zgefnr(nrho,nrho,nvar), STAT=istat)
if (.not. allocated(zgelvl)) allocate(zgelvl(nrho2,nrho2,nvar), STAT=istat)
if (.not. allocated(zgelvr)) allocate(zgelvr(nrho2,nrho2,nvar), STAT=istat)
if (.not. allocated(zgefnl_sop)) then     
    allocate(zgefnl_sop(nrho,nrho,nsopr), STAT=istat)
end if
if (.not. allocated(zgefnr_sop)) then
    allocate(zgefnr_sop(nrho,nrho,nsopr), STAT=istat)
end if
if (.not. allocated(zgelvl_sop)) then     
    allocate(zgelvl_sop(nrho2,nrho2,nsopr), STAT=istat)
end if
if (.not. allocated(zgelvr_sop)) then
    allocate(zgelvr_sop(nrho2,nrho2,nsopr), STAT=istat)
end if
!
! Calculate the adiabatic integral for all the modes (drawers),
! in fact only those fast modes will be used explicitly 
!
zgefnl(1:nrho,1:nrho,1:nvar) = czero
zgefnr(1:nrho,1:nrho,1:nvar) = czero
!
if (lscba_ad) goto 90
!
do idraw=1,nvar
   isgn  = mpm(idraw)
   iorbs = morbs(idraw)
   ispin = mspin(idraw)
   icor  = mcor(idraw)
   ialf  = ilead(idraw)
   im    = mmfff(idraw) - 1
   call dfactorial(im, dtmp1)
!
   if (offcor) then  !  ctheta does not depend on iorbsp
       do iorbsp=1,norbs
          iorbs2 = (iorbsp - 1) * norbs + iorbs
          zcoefl = -eye * cb(iorbs2,ispin,icor,ialf,isgn) * dtmp1
          zcoefr = -eye * cd(iorbs2,ispin,icor,ialf,isgn) * dtmp1
          zgama  = cgama(iorbs2,ispin,icor,ialf,isgn)
          if (igroundsteady .eq. 1) then
              zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
          else if (igroundsteady .gt. 1 .and. ltd2st) then
              zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
          end if
! zgama should NOT be zero
          if (cdabs(zgama) .lt. dpico) then
              write(6,*)
              write(6,*)'prepare_ad: error! zero zgama found ', zgama
              write(6,*)isgn, iorbs, iorbsp, ispin, ialf, icor
              stop
          end if
          if (isgn .eq. 1) then
              do mj=1,nrho
                 do mi=1,nrho
                    cmat1(mi,mj) = dconjg(zamsall(mj,mi,iorbsp,ispin))
                 end do
              end do
          else
              cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbsp,ispin)
          end if
          cmat2(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * (-cunity / zgama)**(im + 1)
          zgefnl(1:nrho,1:nrho,idraw) = zgefnl(1:nrho,1:nrho,idraw) + zcoefl * cmat2(1:nrho,1:nrho)
          zgefnr(1:nrho,1:nrho,idraw) = zgefnr(1:nrho,1:nrho,idraw) + zcoefr * cmat2(1:nrho,1:nrho)
!
          call zgemv('c', nrho2, nnzeig, cunity, zmatder, nrho2, cmat1(1,1), 1, czero, cmat2(1,1), 1)
          ntmp3 = 0
          loop1: do mj=1,nrho
             do mi=1,nrho
                ntmp3 = ntmp3 + 1
                if (ntmp3 .gt. nnzeig) exit loop1
                ctmp1 = cunity / (eye * dvalder(ntmp3) - zgama)**(im + 1) - (-cunity / zgama)**(im + 1)
                cmat2(mi,mj) = cmat2(mi,mj) * ctmp1
             end do
          end do loop1
          call zgemv('n', nrho2, nnzeig, cunity, zmatder, nrho2, cmat2(1,1), 1, czero, cmat1(1,1), 1)
          zgefnl(1:nrho,1:nrho,idraw) = zgefnl(1:nrho,1:nrho,idraw) + zcoefl * cmat1(1:nrho,1:nrho)
          zgefnr(1:nrho,1:nrho,idraw) = zgefnr(1:nrho,1:nrho,idraw) + zcoefr * cmat1(1:nrho,1:nrho)
       end do
!
   else
       zcoefl = -eye * cb(iorbs,ispin,icor,ialf,isgn) * dtmp1
       zcoefr = -eye * cd(iorbs,ispin,icor,ialf,isgn) * dtmp1
       zgama  = cgama(iorbs,ispin,icor,ialf,isgn)
       if (igroundsteady .eq. 1) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       else if (igroundsteady .gt. 1 .and. ltd2st) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       end if
! zgama should NOT be zero
       if (cdabs(zgama) .lt. dpico) then
           write(6,*)
           write(6,*)'prepare_ad: error! zero zgama found ', zgama
           write(6,*)isgn, iorbs, ispin, ialf, icor
           stop
       end if
       if (isgn .eq. 1) then
           do mj=1,nrho
              do mi=1,nrho
                 cmat1(mi,mj) = dconjg(zamsall(mj,mi,iorbs,ispin)) 
              end do
           end do
       else
           cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbs,ispin)
       end if
       cmat2(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * (-cunity / zgama)**(im + 1)
       zgefnl(1:nrho,1:nrho,idraw) = zgefnl(1:nrho,1:nrho,idraw) + zcoefl * cmat2(1:nrho,1:nrho)
       zgefnr(1:nrho,1:nrho,idraw) = zgefnr(1:nrho,1:nrho,idraw) + zcoefr * cmat2(1:nrho,1:nrho)
!
       call zgemv('c', nrho2, nnzeig, cunity, zmatder, nrho2, cmat1(1,1), 1, czero, cmat2(1,1), 1)
       ntmp3 = 0
       loop0: do mj=1,nrho
          do mi=1,nrho
             ntmp3 = ntmp3 + 1
             if (ntmp3 .gt. nnzeig) exit loop0
             ctmp1 = cunity / (eye * dvalder(ntmp3) - zgama)**(im + 1) - (-cunity / zgama)**(im + 1)
             cmat2(mi,mj) = cmat2(mi,mj) * ctmp1
          end do
       end do loop0
       call zgemv('n', nrho2, nnzeig, cunity, zmatder, nrho2, cmat2(1,1), 1, czero, cmat1(1,1), 1)
       zgefnl(1:nrho,1:nrho,idraw) = zgefnl(1:nrho,1:nrho,idraw) + zcoefl * cmat1(1:nrho,1:nrho)
       zgefnr(1:nrho,1:nrho,idraw) = zgefnr(1:nrho,1:nrho,idraw) + zcoefr * cmat1(1:nrho,1:nrho)
   end if
end do
!
90 continue
!
! calculate the Liouville matrices zgelvl and zgelvr 
!
zgelvl(1:nrho2,1:nrho2,1:nvar) = czero
zgelvr(1:nrho2,1:nrho2,1:nvar) = czero
!
if (lscba_ad) goto 100
!
do idraw=1,nvar
   isgn  = mpm(idraw)
   iorbs = morbs(idraw)
   ispin = mspin(idraw)
   icor  = mcor(idraw)
   ialf  = ilead(idraw)
   im    = mmfff(idraw) - 1
   call dfactorial(im, dtmp1)
!
   if (offcor) then
       do iorbsp=1,norbs
          iorbs2 = (iorbsp - 1) * norbs + iorbs
          zcoefl = -eye * cb(iorbs2,ispin,icor,ialf,isgn) * dtmp1
          zcoefr = -eye * cd(iorbs2,ispin,icor,ialf,isgn) * dtmp1
          zgama  = cgama(iorbs2,ispin,icor,ialf,isgn)
          if (igroundsteady .eq. 1) then
              zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
          else if (igroundsteady .gt. 1 .and. ltd2st) then
              zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
          end if
! zgama should NOT be zero
          if (cdabs(zgama) .lt. dpico) then
              write(6,*)
              write(6,*)'prepare_ad: error! zero zgama found ', zgama
              write(6,*)isgn, iorbs, iorbsp, ispin, ialf, icor
              stop
          end if
          if (isgn .eq. 1) then
              do mj=1,nrho
                 do mi=1,nrho
                    cmat1(mi,mj) = dconjg(zamsall(mj,mi,iorbsp,ispin))
                 end do
              end do
          else
              cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbsp,ispin)
          end if
          !
          do mi=1,nrho2
             ctmp1 = cunity / (eye * tmpval2(mi) - zgama)**(im + 1)
             do mj=1,nrho2
                zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * ctmp1
             end do
          end do
          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2)
          !
          cmat2(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * zcoefl
          cmat1(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * zcoefr
          do ij=1,nrho2
             do mi=1,nrho
                do ni=1,nrho
                   ntmp1 = (ni - 1) * nrho + mi
                   do ki=1,nrho
                      ntmp2 = (ni - 1) * nrho + ki
                      ntmp3 = (ki - 1) * nrho + mi
                      zgelvl(ij,ntmp1,idraw) = zgelvl(ij,ntmp1,idraw) + zlmtmp2(ij,ntmp2) * cmat2(ki,mi)
                      zgelvr(ij,ntmp1,idraw) = zgelvr(ij,ntmp1,idraw) + zlmtmp2(ij,ntmp3) * cmat1(ni,ki)
                   end do
                end do
             end do
          end do
       end do
   else
       zcoefl = -eye * cb(iorbs,ispin,icor,ialf,isgn) * dtmp1
       zcoefr = -eye * cd(iorbs,ispin,icor,ialf,isgn) * dtmp1
       zgama  = cgama(iorbs,ispin,icor,ialf,isgn)
       if (igroundsteady .eq. 1) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       else if (igroundsteady .gt. 1 .and. ltd2st) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       end if
! zgama should NOT be zero
       if (cdabs(zgama) .lt. dpico) then
           write(6,*)
           write(6,*)'prepare_ad: error! zero zgama found ', zgama
           write(6,*)isgn, iorbs, ispin, ialf, icor
           stop
       end if
       transa = 'n'
       if (isgn .eq. 1) transa = 'c'
       cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbs,ispin)
!
       do mi=1,nrho2
          ctmp1 = cunity / (eye * tmpval2(mi) - zgama)**(im + 1)
          do mj=1,nrho2
             zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * ctmp1
          end do
       end do
       call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2)
!
       !cmat2(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * zcoefl
       !cmat1(1:nrho,1:nrho) = cmat1(1:nrho,1:nrho) * zcoefr
       !do ij=1,nrho2
       !   do mi=1,nrho
       !      do ni=1,nrho
       !         ntmp1 = (ni - 1) * nrho + mi
       !         do ki=1,nrho
       !            ntmp2 = (ni - 1) * nrho + ki
       !            ntmp3 = (ki - 1) * nrho + mi
       !            zgelvl(ij,ntmp1,idraw) = zgelvl(ij,ntmp1,idraw) + zlmtmp2(ij,ntmp2) * cmat2(ki,mi)
       !            zgelvr(ij,ntmp1,idraw) = zgelvr(ij,ntmp1,idraw) + zlmtmp2(ij,ntmp3) * cmat1(ni,ki)
       !         end do
       !      end do
       !   end do
       !end do
       call zhil2liou('l', transa, nrho, cmat1, nrho, zlmtmp1, nrho2)
       call zgemm('n', 'n', nrho2, nrho2, nrho2, zcoefl, zlmtmp2, nrho2, zlmtmp1, nrho2, cunity, zgelvl(1,1,idraw), nrho2) 
       call zhil2liou('r', transa, nrho, cmat1, nrho, zlmtmp1, nrho2)
       call zgemm('n', 'n', nrho2, nrho2, nrho2, zcoefr, zlmtmp2, nrho2, zlmtmp1, nrho2, cunity, zgelvr(1,1,idraw), nrho2) 
   end if
end do
!
100 continue
!
if (lscba_ad) then
    do idraw=1,nvar
       isgn  = mpm(idraw)
       iorbs = morbs(idraw)
       ispin = mspin(idraw)
       icor  = mcor(idraw)
       ialf  = ilead(idraw)
       im    = mmfff(idraw) - 1
       if (im .gt. 2) then
           write(6,*)
           write(6,*)'prepare_ad: error! im > 2 found for idraw = ', idraw
           write(6,*)'prepare_ad: code not available yet '
           stop
       end if
       call dfactorial(im, dtmp0)
       iorbs2 = (iorbs - 1) * norbs + iorbs
       if (offcor) then
           zgama = cgama(iorbs2,ispin,icor,ialf,isgn)
       else
           zgama = cgama(iorbs,ispin,icor,ialf,isgn)
       end if
       if (igroundsteady .eq. 1) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       else if (igroundsteady .gt. 1 .and. ltd2st) then
           zgama = zgama + eye * dipm(idraw) * eleadinfty(ialf,ispin)
       end if
       zw = -eye * zgama  ! zw should not be zero (probably no need to check)
       zlself(1:nrho2,1:nrho2) = czero
       zself1(1:nrho2,1:nrho2) = czero
       zself2(1:nrho2,1:nrho2) = czero
       !
       do jsgn=1,nsgn
          jsgnb = nsgn - jsgn + 1
          do jjalf=1,nalf
             do jjspin=1,nspin
                ctmp0 = czero
                if (igroundsteady .eq. 1) then
                    ctmp0 = eye * dpm(jsgn) * eleadinfty(jjalf,jjspin)
                else if (igroundsteady .gt. 1 .and. ltd2st) then
                    ctmp0 = eye * dpm(jsgn) * eleadinfty(jjalf,jjspin)
                end if
                do jorbs1=1,norbs 
                   do jorbs2=1,norbs 
                      if (.not. offcor .and. jorbs1 .ne. jorbs2) cycle
                      jjorbs = (jorbs2 - 1) * norbs + jorbs1
                      cmtmp1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,jorbs1,jjspin) 
                      cmtmp2(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,jorbs2,jjspin) 
                      transa = 'n'
                      transb = 'c'
                      if (jsgn .eq. 1) then
                          transa = 'c'
                          transb = 'n'
                      end if
                      call zhil2liou('l', transb, nrho, cmtmp1, nrho, zlmtmp1, nrho2)
                      call zhil2liou('r', transb, nrho, cmtmp1, nrho, zlmtmp2, nrho2)
                      zlmtmp4(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
                      !
                      zlv1(1:nrho2) = czero
                      zlv2(1:nrho2) = czero
                      if (im .ge. 1) then
                          zlv3(1:nrho2) = czero
                          zlv4(1:nrho2) = czero
                      end if
                      if (im .ge. 2) then
                          zlv5(1:nrho2) = czero
                          zlv6(1:nrho2) = czero
                      end if
                      !
                      do jcor=1,ncor
                         jm = 0
                         jmats = jcor - ndrude - nmats
                         if (numfff .gt. 0 .and. jmats .gt. 0) then
                             jm = mfff1d(jmats) - 1
                         end if
                         call dfactorial(jm, dtmp1)
                         if (offcor) then
                             zgamaj = cgama(jjorbs,jjspin,jcor,jjalf,jsgn) + ctmp0 
                             do ni=1,nrho2
                                ctmp1 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 1) * dtmp1
                                zlv1(ni) = zlv1(ni) + cb(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp1
                                zlv2(ni) = zlv2(ni) + cd(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp1
                             end do
                         else
                             zgamaj = cgama(jorbs1,jjspin,jcor,jjalf,jsgn) + ctmp0 
                             do ni=1,nrho2
                                ctmp1 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 1) * dtmp1
                                zlv1(ni) = zlv1(ni) + cb(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp1
                                zlv2(ni) = zlv2(ni) + cd(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp1
                             end do
                         end if
                         if (im .ge. 1) then
                             call dfactorial(jm+1, dtmp2)
                             if (offcor) then
                                 do ni=1,nrho2
                                    ctmp2 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 2) * dtmp2 * eye
                                    zlv3(ni) = zlv3(ni) + cb(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp2
                                    zlv4(ni) = zlv4(ni) + cd(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp2
                                 end do
                             else
                                 do ni=1,nrho2
                                    ctmp2 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 2) * dtmp2 * eye
                                    zlv3(ni) = zlv3(ni) + cb(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp2
                                    zlv4(ni) = zlv4(ni) + cd(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp2
                                 end do
                             end if
                         end if
                         if (im .ge. 2) then
                             call dfactorial(jm+2, dtmp3)
                             if (offcor) then
                                 do ni=1,nrho2
                                    ctmp3 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 3) * dtmp3 * (-1.d0)
                                    zlv5(ni) = zlv5(ni) + cb(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp3
                                    zlv6(ni) = zlv6(ni) + cd(jjorbs,jjspin,jcor,jjalf,jsgn) * ctmp3
                                 end do
                             else
                                 do ni=1,nrho2
                                    ctmp3 = (-1.d0 / (zgamaj + eye * (zw - tmpval2(ni))))**(jm + 3) * dtmp3 * (-1.d0)
                                    zlv5(ni) = zlv5(ni) + cb(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp3
                                    zlv6(ni) = zlv6(ni) + cd(jorbs1,jjspin,jcor,jjalf,jsgn) * ctmp3
                                 end do
                             end if
                         end if
                      end do
                      do mj=1,nrho2
                         do mi=1,nrho2
                            zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv1(mi)
                         end do
                      end do
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                      call zhil2liou('l', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, -eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zlself, nrho2)
                      do mj=1,nrho2
                         do mi=1,nrho2
                            zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv2(mi)
                         end do
                      end do
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                      call zhil2liou('r', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                      call zgemm('n', 'n', nrho2, nrho2, nrho2, eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zlself, nrho2)
                      if (im .ge. 1) then
                          do mj=1,nrho2
                             do mi=1,nrho2
                                zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv3(mi)
                             end do
                          end do
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                          call zhil2liou('l', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, -eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zself1, nrho2)
                          do mj=1,nrho2
                             do mi=1,nrho2
                                zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv4(mi)
                             end do
                          end do
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                          call zhil2liou('r', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zself1, nrho2)
                      end if
                      if (im .ge. 2) then
                          do mj=1,nrho2
                             do mi=1,nrho2
                                zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv5(mi)
                             end do
                          end do
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                          call zhil2liou('l', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, -eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zself2, nrho2)
                          do mj=1,nrho2
                             do mi=1,nrho2
                                zlmtmp1(mi,mj) = dconjg(zwork2(mj,mi)) * zlv6(mi)
                             end do
                          end do
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zwork2, nrho2, zlmtmp1, nrho2, czero, zlmtmp2, nrho2) 
                          call zhil2liou('r', transa, nrho, cmtmp2, nrho, zlmtmp1, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp2, nrho2, zlmtmp1, nrho2, czero, zlmtmp3, nrho2)
                          call zgemm('n', 'n', nrho2, nrho2, nrho2, eye, zlmtmp4, nrho2, zlmtmp3, nrho2, cunity, zself2, nrho2)
                      end if
                   end do
                end do
             end do
          end do
       end do
       !
       zlmtmp1(1:nrho2,1:nrho2) = czero
       do ni=1,nrho2
          zlmtmp1(ni,ni) = zw
       end do
       zlmtmp1(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlsys(1:nrho2,1:nrho2) - zlself(1:nrho2,1:nrho2)
       call zgetrf(nrho2, nrho2, zlmtmp1, nrho2, ipiv, info)
       if (info .ne. 0) then
           write(6,*)
           write(6,*)'prepare_ad: error with zgetrf for calculation of zlu ', info
           write(6,*)'prepare_ad: idraw = ', idraw
           stop
       end if
       call zgetri(nrho2, zlmtmp1, nrho2, ipiv, zwork, nwork, info)
       if (info .ne. 0) then
           write(6,*)
           write(6,*)'prepare_ad: error with zgetri for calculation of zlu ', info
           write(6,*)'prepare_ad: idraw = ', idraw
           stop
       end if
       zlu(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) * eye
       !
       if (im .ge. 1) then
           zlmtmp1(1:nrho2,1:nrho2) = -zself1(1:nrho2,1:nrho2)
           do ni=1,nrho2
              zlmtmp1(ni,ni) = zlmtmp1(ni,ni) + cunity
           end do
           call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp1, nrho2,     zlu, nrho2, czero, zlmtmp2, nrho2)
           call zgemm('n', 'n', nrho2, nrho2, nrho2,    eye,     zlu, nrho2, zlmtmp2, nrho2, czero,    zlu1, nrho2)
       end if
       if (im .ge. 2) then
           call zgemm('n', 'n', nrho2, nrho2, nrho2,    eye,    zlu1, nrho2, zlmtmp2, nrho2,  czero,    zlu2, nrho2)
           call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity,  zself2, nrho2,     zlu, nrho2,  czero, zlmtmp3, nrho2)
           call zgemm('n', 'n', nrho2, nrho2, nrho2,   -eye,     zlu, nrho2, zlmtmp3, nrho2, cunity,    zlu2, nrho2)
           call zgemm('n', 'n', nrho2, nrho2, nrho2, cunity, zlmtmp1, nrho2,    zlu1, nrho2,  czero, zlmtmp3, nrho2)
           call zgemm('n', 'n', nrho2, nrho2, nrho2,    eye,     zlu, nrho2, zlmtmp3, nrho2, cunity,    zlu2, nrho2)
       end if
       !
       if (im .eq. 0) then
           zlmtmp2(1:nrho2,1:nrho2) = zlu(1:nrho2,1:nrho2)
       else if (im .eq. 1) then
           zlmtmp2(1:nrho2,1:nrho2) = zlu1(1:nrho2,1:nrho2) * (-eye)
       else if (im .eq. 2) then
           zlmtmp2(1:nrho2,1:nrho2) = zlu2(1:nrho2,1:nrho2) * (-1.d0)
       else
           write(6,*)
           write(6,*)'prepare_ad: error! im > 2 found for calculation of zgelvl/zgelvr '
           write(6,*)'prepare_ad: code unavailable yet '
           stop
       end if
       !
       do iorbsp=1,norbs
          if (.not. offcor .and. iorbs .ne. iorbsp) cycle
          iorbs2 = (iorbsp - 1) * norbs + iorbs
          if (offcor) then
              zcoefl = -eye * cb(iorbs2,ispin,icor,ialf,isgn) 
              zcoefr = -eye * cd(iorbs2,ispin,icor,ialf,isgn) 
          else
              zcoefl = -eye * cb(iorbsp,ispin,icor,ialf,isgn) 
              zcoefr = -eye * cd(iorbsp,ispin,icor,ialf,isgn) 
          end if
          transa = 'n'
          if (isgn .eq. 1) transa = 'c'
          cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbsp,ispin) 
          call zhil2liou('l', transa, nrho, cmat1, nrho, zlmtmp1, nrho2)
          call zgemm('n', 'n', nrho2, nrho2, nrho2, zcoefl, zlmtmp2, nrho2, zlmtmp1, nrho2, cunity, zgelvl(1,1,idraw), nrho2) 
          call zhil2liou('r', transa, nrho, cmat1, nrho, zlmtmp1, nrho2)
          call zgemm('n', 'n', nrho2, nrho2, nrho2, zcoefr, zlmtmp2, nrho2, zlmtmp1, nrho2, cunity, zgelvr(1,1,idraw), nrho2) 
          !
          if (isgn .eq. 1) then
              do nj=1,nrho
                 do ni=1,nrho
                    cmat1(ni,nj) = dconjg(zamsall(nj,ni,iorbsp,ispin))
                 end do
              end do
          else
              cmat1(1:nrho,1:nrho) = zamsall(1:nrho,1:nrho,iorbsp,ispin) 
          end if
          call zgemv('n', nrho2, nrho2, cunity, zlmtmp2, nrho2, cmat1(1,1), 1, czero, cmat2(1,1), 1)
          zgefnl(1:nrho,1:nrho,idraw) = zgefnl(1:nrho,1:nrho,idraw) + zcoefl * cmat2(1:nrho,1:nrho)
          zgefnr(1:nrho,1:nrho,idraw) = zgefnr(1:nrho,1:nrho,idraw) + zcoefr * cmat2(1:nrho,1:nrho)
       end do
    end do 
end if
!
zgefnl_sop(1:nrho,1:nrho,1:nsopr) = czero
zgefnr_sop(1:nrho,1:nrho,1:nsopr) = czero
zgelvl_sop(1:nrho2,1:nrho2,1:nsopr) = czero
zgelvr_sop(1:nrho2,1:nrho2,1:nsopr) = czero
do isgn=1,nsgn
   do ispin=1,nspin
      do iorbs=1,norbs
         call look4sysopr(isgn, iorbs, ispin, isopr)
         do ialf=1,nalf
            loop_s: do icor=1,ncor
               call look4drawer(isgn, ialf, iorbs, ispin, icor, idraw)
               if (.not. lad_fast .and. mslow(idraw) .eq. 1) cycle loop_s  ! only sum up the fast modes
               if (      lad_fast .and. mslow(idraw) .eq. 0) cycle loop_s  ! only sum up the slow modes
               !if (mslow(idraw) .eq. 1) cycle loop_s                         
               ! only sum up the fast modes
               do nj=1,nrho
                  do ni=1,nrho
                     zgefnl_sop(ni,nj,isopr) = zgefnl_sop(ni,nj,isopr) + zgefnl(ni,nj,idraw)
                     zgefnr_sop(ni,nj,isopr) = zgefnr_sop(ni,nj,isopr) + zgefnr(ni,nj,idraw)
                  end do
               end do
               do nj=1,nrho2
                  do ni=1,nrho2
                     zgelvl_sop(ni,nj,isopr) = zgelvl_sop(ni,nj,isopr) + zgelvl(ni,nj,idraw)
                     zgelvr_sop(ni,nj,isopr) = zgelvr_sop(ni,nj,isopr) + zgelvr(ni,nj,idraw)
                  end do
               end do
            end do loop_s
         end do
      end do
   end do
end do
!
! analyze the sparsity of zgefnl/zgefnr 
!
mi = 0
mj = 0
do idraw=1,nvar
   ntmp1 = 0
   ntmp2 = 0
   do nj=1,nrho
      do ni=1,nrho
         if ( cdabs(zgefnl(ni,nj,idraw)) .gt. dcut ) ntmp1 = ntmp1 + 1
         if ( cdabs(zgefnr(ni,nj,idraw)) .gt. dcut ) ntmp2 = ntmp2 + 1
      end do
   end do
   mi = max(mi, ntmp1)
   mj = max(mj, ntmp2)
end do
write(6,*)'prepare_ad: max nonzero elements in zgefnl    ', mi
write(6,*)'prepare_ad: max nonzero elements in zgefnr    ', mj
call flush(6)
!
mi = max(mi, mj)
if (.not. allocated(nnz_zgf))                   allocate(nnz_zgf(nvar),      STAT=istat)
if (.not. allocated(irow_zgf)  .and. mi .gt. 0) allocate(irow_zgf(mi,nvar),  STAT=istat)
if (.not. allocated(icol_zgf)  .and. mi .gt. 0) allocate(icol_zgf(mi,nvar),  STAT=istat)
if (.not. allocated(cval_zgfl) .and. mi .gt. 0) allocate(cval_zgfl(mi,nvar), STAT=istat)
if (.not. allocated(cval_zgfr) .and. mi .gt. 0) allocate(cval_zgfr(mi,nvar), STAT=istat)
!
mi = 0
mj = 0
do idraw=1,nvar
   ntmp1 = 0
   ntmp2 = 0
   do nj=1,nrho2
      do ni=1,nrho2
         if ( cdabs(zgelvl(ni,nj,idraw)) .gt. dcut ) ntmp1 = ntmp1 + 1
         if ( cdabs(zgelvr(ni,nj,idraw)) .gt. dcut ) ntmp2 = ntmp2 + 1
      end do
   end do
   mi = max(mi, ntmp1)
   mj = max(mj, ntmp2)
end do
write(6,*)'prepare_ad: max nonzero elements in zgelvl    ', mi
write(6,*)'prepare_ad: max nonzero elements in zgelvr    ', mj
call flush(6)
!
if (.not. allocated(nnz_lvl))                  allocate(nnz_lvl(nvar),     STAT=istat)
if (.not. allocated(irow_lvl) .and. mi .gt. 0) allocate(irow_lvl(mi,nvar), STAT=istat)
if (.not. allocated(icol_lvl) .and. mi .gt. 0) allocate(icol_lvl(mi,nvar), STAT=istat)
if (.not. allocated(cval_lvl) .and. mi .gt. 0) allocate(cval_lvl(mi,nvar), STAT=istat)
if (.not. allocated(nnz_lvr))                  allocate(nnz_lvr(nvar),     STAT=istat)
if (.not. allocated(irow_lvr) .and. mj .gt. 0) allocate(irow_lvr(mj,nvar), STAT=istat)
if (.not. allocated(icol_lvr) .and. mj .gt. 0) allocate(icol_lvr(mj,nvar), STAT=istat)
if (.not. allocated(cval_lvr) .and. mj .gt. 0) allocate(cval_lvr(mj,nvar), STAT=istat)
!
if (.not. allocated(ind_lvl)) allocate(ind_lvl(nrho2,nvar), STAT=istat)
if (.not. allocated(nnc_lvl)) allocate(nnc_lvl(nrho2,nvar), STAT=istat)
if (.not. allocated(ind_lvr)) allocate(ind_lvr(nrho2,nvar), STAT=istat)
if (.not. allocated(nnc_lvr)) allocate(nnc_lvr(nrho2,nvar), STAT=istat)
!
mi = 0
mj = 0
do isopr=1,nsopr
   ntmp1 = 0
   ntmp2 = 0
   do nj=1,nrho
      do ni=1,nrho
         if ( cdabs(zgefnl_sop(ni,nj,isopr)) .gt. dcut ) ntmp1 = ntmp1 + 1
         if ( cdabs(zgefnr_sop(ni,nj,isopr)) .gt. dcut ) ntmp2 = ntmp2 + 1
      end do
   end do
   mi = max(mi, ntmp1)
   mj = max(mj, ntmp2)
end do
write(6,*)'prepare_ad: max nonzero elements in zgefnl_sop ', mi
write(6,*)'prepare_ad: max nonzero elements in zgefnr_sop ', mj
write(6,*)'prepare_ad: total number of elements           ', nrho2
call flush(6)
!
mi = max(mi, mj)
if (.not. allocated(nnz_sop))                   allocate(nnz_sop(nsopr),      STAT=istat)
if (.not. allocated(irow_sop)  .and. mi .gt. 0) allocate(irow_sop(mi,nsopr),  STAT=istat)
if (.not. allocated(icol_sop)  .and. mi .gt. 0) allocate(icol_sop(mi,nsopr),  STAT=istat)
if (.not. allocated(cval_sopl) .and. mi .gt. 0) allocate(cval_sopl(mi,nsopr), STAT=istat)
if (.not. allocated(cval_sopr) .and. mi .gt. 0) allocate(cval_sopr(mi,nsopr), STAT=istat)
!
mi = 0
mj = 0
do isopr=1,nsopr
   ntmp1 = 0
   ntmp2 = 0
   do nj=1,nrho2
      do ni=1,nrho2
         if ( cdabs(zgelvl_sop(ni,nj,isopr)) .gt. dcut ) ntmp1 = ntmp1 + 1
         if ( cdabs(zgelvr_sop(ni,nj,isopr)) .gt. dcut ) ntmp2 = ntmp2 + 1
      end do
   end do
   mi = max(mi, ntmp1)
   mj = max(mj, ntmp2)
end do
write(6,*)'prepare_ad: max nonzero elements in zgelvl_sop ', mi
write(6,*)'prepare_ad: max nonzero elements in zgelvr_sop ', mj
write(6,*)'prepare_ad: total number of elements           ', nrho2**2
call flush(6)
!
if (.not. allocated(nnz_lvlsop))                  allocate(nnz_lvlsop(nsopr),     STAT=istat)
if (.not. allocated(irow_lvlsop) .and. mi .gt. 0) allocate(irow_lvlsop(mi,nsopr), STAT=istat)
if (.not. allocated(icol_lvlsop) .and. mi .gt. 0) allocate(icol_lvlsop(mi,nsopr), STAT=istat)
if (.not. allocated(cval_lvlsop) .and. mi .gt. 0) allocate(cval_lvlsop(mi,nsopr), STAT=istat)
if (.not. allocated(nnz_lvrsop))                  allocate(nnz_lvrsop(nsopr),     STAT=istat)
if (.not. allocated(irow_lvrsop) .and. mj .gt. 0) allocate(irow_lvrsop(mj,nsopr), STAT=istat)
if (.not. allocated(icol_lvrsop) .and. mj .gt. 0) allocate(icol_lvrsop(mj,nsopr), STAT=istat)
if (.not. allocated(cval_lvrsop) .and. mj .gt. 0) allocate(cval_lvrsop(mj,nsopr), STAT=istat)
!
do idraw=1,nvar
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         if ( cdabs(zgefnl(ni,nj,idraw)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_zgf (ntmp1,idraw) = ni
             icol_zgf (ntmp1,idraw) = nj
             cval_zgfl(ntmp1,idraw) = zgefnl(ni,nj,idraw)
             cval_zgfr(ntmp1,idraw) = zgefnr(ni,nj,idraw)
         end if
      end do
   end do
   nnz_zgf(idraw) = ntmp1
!
   ntmp1 = 0
   do nj=1,nrho2
      ind_lvl(nj,idraw) = -1
      nnc_lvl(nj,idraw) = 0
      do ni=1,nrho2
         if ( cdabs(zgelvl(ni,nj,idraw)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_lvl(ntmp1,idraw) = ni
             icol_lvl(ntmp1,idraw) = nj
             cval_lvl(ntmp1,idraw) = zgelvl(ni,nj,idraw)
             if (ind_lvl(nj,idraw) .lt. 0) then
                 ind_lvl(nj,idraw) = ntmp1
             end if
             nnc_lvl(nj,idraw) = nnc_lvl(nj,idraw) + 1
         end if
      end do
   end do
   nnz_lvl(idraw) = ntmp1
   ntmp1 = 0
   do nj=1,nrho2
      ntmp1 = ntmp1 + nnc_lvl(nj,idraw)
   end do
   if (ntmp1 .ne. nnz_lvl(idraw)) then
       write(6,*)
       write(6,*)'prepare_ad: error! sum of nnc_lvl is not equal to nnz_lvl '
       write(6,*)idraw, nnz_lvl(idraw), ntmp1
       stop
   end if
!
   ntmp1 = 0
   do nj=1,nrho2
      ind_lvr(nj,idraw) = -1
      nnc_lvr(nj,idraw) = 0
      do ni=1,nrho2
         if ( cdabs(zgelvr(ni,nj,idraw)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_lvr(ntmp1,idraw) = ni
             icol_lvr(ntmp1,idraw) = nj
             cval_lvr(ntmp1,idraw) = zgelvr(ni,nj,idraw)
             if (ind_lvr(nj,idraw) .lt. 0) then
                 ind_lvr(nj,idraw) = ntmp1
             end if
             nnc_lvr(nj,idraw) = nnc_lvr(nj,idraw) + 1
         end if
      end do
   end do
   nnz_lvr(idraw) = ntmp1
   ntmp1 = 0
   do nj=1,nrho2
      ntmp1 = ntmp1 + nnc_lvr(nj,idraw)
   end do
   if (ntmp1 .ne. nnz_lvr(idraw)) then
       write(6,*)
       write(6,*)'prepare_ad: error! sum of nnc_lvr is not equal to nnz_lvr '
       write(6,*)idraw, nnz_lvr(idraw), ntmp1
       stop
   end if
end do
!
!write(6,*)
!write(6,*)'ind_lvl, nnc_lvl '
!do idraw=1,nvar
!   do nj=1,nrho2
!      write(6,1000)idraw, nj, ind_lvl(nj,idraw), nnc_lvl(nj,idraw)
!   end do
!end do
!write(6,*)
!write(6,*)'ind_lvr, nnc_lvr '
!do idraw=1,nvar
!   do nj=1,nrho2
!      write(6,1000)idraw, nj, ind_lvr(nj,idraw), nnc_lvr(nj,idraw)
!   end do
!end do
!call flush(6)
!1000 format(1x, 4(I6, 2x))
!
do isopr=1,nsopr
   ntmp1 = 0
   do nj=1,nrho
      do ni=1,nrho
         if ( cdabs(zgefnl_sop(ni,nj,isopr)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_sop (ntmp1,isopr) = ni
             icol_sop (ntmp1,isopr) = nj
             cval_sopl(ntmp1,isopr) = zgefnl_sop(ni,nj,isopr)
             cval_sopr(ntmp1,isopr) = zgefnr_sop(ni,nj,isopr)
         end if
      end do
   end do
   nnz_sop(isopr) = ntmp1
!
   ntmp1 = 0
   do nj=1,nrho2
      do ni=1,nrho2
         if ( cdabs(zgelvl_sop(ni,nj,isopr)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_lvlsop(ntmp1,isopr) = ni
             icol_lvlsop(ntmp1,isopr) = nj
             cval_lvlsop(ntmp1,isopr) = zgelvl_sop(ni,nj,isopr)
         end if
      end do
   end do
   nnz_lvlsop(isopr) = ntmp1
!
   ntmp1 = 0
   do nj=1,nrho2
      do ni=1,nrho2
         if ( cdabs(zgelvr_sop(ni,nj,isopr)) .gt. dcut ) then
             ntmp1 = ntmp1 + 1
             irow_lvrsop(ntmp1,isopr) = ni
             icol_lvrsop(ntmp1,isopr) = nj
             cval_lvrsop(ntmp1,isopr) = zgelvr_sop(ni,nj,isopr)
         end if
      end do
   end do
   nnz_lvrsop(isopr) = ntmp1
end do
!
if (.not. allocated(nlvrow)) allocate(nlvrow(nrho2), STAT=istat)
if (.not. allocated(nlvcol)) allocate(nlvcol(nrho2), STAT=istat)
if (.not. allocated(nnrc2v)) allocate(nnrc2v(nrho,nrho), STAT=istat)
if (.not. allocated(ncrc2v)) allocate(ncrc2v(nrho,nrho), STAT=istat)
ntmp1 = 0
do nj=1,nrho
   do ni=1,nrho
      ntmp1 = ntmp1 + 1
      nlvrow(ntmp1) = ni
      nlvcol(ntmp1) = nj
      nnrc2v(ni,nj) = ntmp1
      ncrc2v(nj,ni) = ntmp1
   end do
end do
!
deallocate(ipiv, STAT=istat)
deallocate(zwork, rwork, tmpval, cmat1, cmat2, STAT=istat)
deallocate(zwork2, tmpval2, STAT=istat)
deallocate(zlsys, STAT=istat)
if (lscba_ad) then
    deallocate(zlself, zself1, zself2, zlu, zlu1, zlu2, STAT=istat)
    deallocate(zlv1, zlv2, zlv3, zlv4, zlv5, zlv6, STAT=istat)
end if
!
if (lcop_ad) then
    write(6,*)
    write(6,*)'prepare_ad: lcop_ad=T '
    write(6,*)'prepare_ad: use zgelvl, zgelvr, zgelvl_sop, zgelvr_sop '
else
    write(6,*)
    write(6,*)'prepare_ad: lcop_ad=F '
    write(6,*)'prepare_ad: use zgefnl, zgefnr, zgefnl_sop, zgefnr_sop '
end if
call flush(6)
!
!open(99, file='zgelvl.dat', form='formatted', status='unknown')
!rewind(99)
!do idraw=1,nvar
!   do ni=1,nnz_lvl(idraw)
!      write(99,*)idraw, ni, irow_lvl(ni,idraw), icol_lvl(ni,idraw), cval_lvl(ni,idraw)
!   end do
!end do
!close(99)
!
write(6,*)
write(6,*)'prepare_ad: completed '
call flush(6)
!
!stop
!
return
end subroutine prepare_ad
!
subroutine dfactorial(nin, dout)
implicit none
!
integer, intent(in)  :: nin
real*8,  intent(out) :: dout
integer              :: ni
!
dout = 1.d0
if (nin .lt. 0) then
    write(6,*)
    write(6,*)'dfactorial: error! input is negative ', nin
    stop
end if
!
do ni=1,nin
   dout = dout * dble(ni)
end do
!
return
end subroutine dfactorial
