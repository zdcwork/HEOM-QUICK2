subroutine checkrho
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs
integer                 :: nmax_element, nzero, ncount
integer                 :: index_nk(MAXTIER)
integer,   allocatable  :: nnonzero(:)
real*8,    allocatable  :: dnonzero(:,:)
logical                 :: lcheck
integer                 :: isgn, iorbs, ispin, ialf, idraw
integer                 :: jopr(maxtier), ntmp1, ntmp2, ntmp3, ntmp4
integer*8               :: lnumado(maxtier)
real*8                  :: dvalado(maxtier)
integer*8               :: lni, lnj, lnk
real*8                  :: dtmp1, dtmp2, dtmp3, dtmp4
integer,    allocatable :: tmpmcor(:,:)
real*8,     allocatable :: tmpval(:), rwork(:), tmpmm0(:)
complex*16, allocatable :: work(:)
integer                 :: mcortmp0(maxtier-1)
integer*8               :: lnum5, lnum6, lnum7, lnum8, ltmp1
real*8                  :: d5, d6, d7, d8
real*8                  :: per5, per6, per7, per8
!
if (lsparse) then
   call checkrho_spa
   return
end if
!
d5 = 1.d-5
d6 = 1.d-6
d7 = 1.d-7
d8 = 1.d-8
lcheck = .true.
write(6,*)
write(6,*)'checkrho: check density matrix rho0 '
do ni=1,nrho
  if (dble(rho(ni,ni,1)) .lt. 0.d0) then
    write(6,*)'checkrho: Warning, negative diag. ', ni, dble(rho(ni,ni,1))
  end if
end do
call flush(6)
!
! rho0 should be a positive definite matrix
!
cmtmp1(1:nrho, 1:nrho) = rho(1:nrho, 1:nrho, 1)
!
lwork = nrho**2
allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
call zheev('N', 'u', nrho, cmtmp1, nrho, tmpval, work, lwork, rwork, info)
if (info .ne. 0) then
   write(6,*)
   write(6,*)'checkrho: error! diagnalization failed in checkrho ', info
   stop
end if
!
do ni=1,nrho
   if (tmpval(ni) .lt. 0.d0) then
      write(6,*)'checkrho: Warning, negative eigen.', ni, tmpval(ni)
      lcheck = .false.
   end if
end do
call flush(6)
deallocate(work ,rwork, tmpval, STAT=istat)
!
write(6,*)
write(6,*)'checkrho: positivite-definiteness check complete! '
if (lcheck) then
  write(6,*)'checkrho: rho is found OK! '
else
  write(6,*)'checkrho: Warning! rho is not positive definite! '
end if
call flush(6)
!
! check sparsity of rho(*,*,*)
!
dmtmp1(1:nrho,1:nrho) = 0.d0
do lni=1,nunk
   do nj=1,nrho
      do ni=1,nrho
         dmtmp1(ni,nj) = max(dmtmp1(ni,nj), cdabs(rho(ni,nj,lni)))
      end do
   end do
end do
write(6,*)
write(6,*)'checkrho: maximal values for elements of all density matrices'
call amatout(nrho,nrho,dmtmp1)
call flush(6)
!
lnum5 = 0
lnum6 = 0
lnum7 = 0
lnum8 = 0
do nj=1,nrho
   do ni=1,nrho
      if (dmtmp1(ni,nj) .ge. d5) then
         lnum5 = lnum5 + 1
      end if
      if (dmtmp1(ni,nj) .ge. d6) then
         lnum6 = lnum6 + 1
      end if
      if (dmtmp1(ni,nj) .ge. d7) then
         lnum7 = lnum7 + 1
      end if
      if (dmtmp1(ni,nj) .ge. d8) then
         lnum8 = lnum8 + 1
      end if
   end do
end do
write(6,*)'checkrho: overall elements ', nrho**2
write(6,1003)d5, lnum5, dble(lnum5) / dble(nrho*nrho) * 100.d0
write(6,1003)d6, lnum6, dble(lnum6) / dble(nrho*nrho) * 100.d0
write(6,1003)d7, lnum7, dble(lnum7) / dble(nrho*nrho) * 100.d0
write(6,1003)d8, lnum8, dble(lnum8) / dble(nrho*nrho) * 100.d0
1003 format('# of elements larger than ', e12.4e2, I8, ' percent ', f8.3)
call flush(6)
!
!!!!!!! debug
! check the first ADO
!
!call cmatout(nrho,nrho,rho(1,1,2),dmtmp1)
!call cmatout(nrho,nrho,rho(1,1,3),dmtmp1)
!!
!do ni=1,nrho
!   do nj=1,nrho
!!      cmtmp1(ni,nj) = rho(ni,nj,2) + dconjg(rho(nj,ni,2))
!!      cmtmp1(ni,nj) = dcmplx(dble(rho(ni,nj,2)), 0.d0)
!!      cmtmp1(ni,nj) = eye * (rho(ni,nj,2) - dconjg(rho(nj,ni,2)))
!!      cmtmp1(ni,nj) = rho(ni,nj,3) + dconjg(rho(nj,ni,3))
!      cmtmp1(ni,nj) = eye * (rho(ni,nj,3) - dconjg(rho(nj,ni,3)))
!   end do
!end do
!allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
!call zheev('N', 'u', nrho, cmtmp1, nrho, tmpval, work, lwork, rwork, info)
!if (info .ne. 0) then
!    write(6,*)
!    write(6,*)'checkrho: error! diagnalization failed 2 ', info
!    stop
!end if
!write(6,*)'checkrho: xxx eigenvalues of ADO + H.c. xxx '
!do ni=1,nrho
!   write(6,*)ni, tmpval(ni)
!end do
!deallocate(work ,rwork, tmpval, STAT=istat)
!call flush(6)
!
!!!!!!! debug
!
! check level of significance for rho at each iter
!
write(6,*)
write(6,*)'checkrho: check level of significance of rho at each tier'
do itier=2,ntier
   lnj   = nlast(itier) - nfirst(itier) + 1
   lnum5 = 0
   lnum6 = 0
   lnum7 = 0
   lnum8 = 0
   do lni=nfirst(itier),nlast(itier)
      call cmaxmat(nrho,nrho,rho(1,1,lni),nrho,dtmp1)     
      if (dtmp1 .le. d5) lnum5 = lnum5 + 1
      if (dtmp1 .le. d6) lnum6 = lnum6 + 1
      if (dtmp1 .le. d7) lnum7 = lnum7 + 1
      if (dtmp1 .le. d8) lnum8 = lnum8 + 1
   end do
   per5 = dble(lnum5) / dble(lnj) * 100.d0
   per6 = dble(lnum6) / dble(lnj) * 100.d0
   per7 = dble(lnum7) / dble(lnj) * 100.d0
   per8 = dble(lnum8) / dble(lnj) * 100.d0
   write(6,*)'itier = ', itier, ' # of ADOs ', lnj
   write(6,1002)d5, lnum5, per5
   write(6,1002)d6, lnum6, per6
   write(6,1002)d7, lnum7, per7
   write(6,1002)d8, lnum8, per8
end do
1002 format('# of ADOs smaller than ', e12.4e2, I8, ' percent ', f8.3)
call flush(6)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (lchkbig) then 
   itier = ntier
   dchkbig = max(dchkbig, d8)
   open(unit=44, file='ado_big.data', form='formatted', status='unknown')
   rewind(44)
   write(6,*)
   write(6,*)'checkrho: writing significant ADOs at tier', itier, 'to <ado_big.data> '
   write(6,*)'checkrho: dchkbig = ', dchkbig
   write(6,*)'checkrho: lni, iball, mpm, morbs, mspin, ilead, mcor '
!
   do lni=nfirst(itier),nlast(itier)
      call cmaxmat(nrho, nrho, rho(1,1,lni), nrho, dtmp1)
      if (dtmp1 .lt. dchkbig) cycle
      do ni=1,itier-1
         idraw = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
         write(44,1004)lni, ni, mpm(idraw), morbs(idraw), mspin(idraw), ilead(idraw), mcor(idraw)
      end do
      call flush(44)
   end do
   close(44)
   write(6,*)'checkrho: writing to <ado_big.data> finished'
   call flush(6)
!
   lnumado(1:ntier) = 0
   dvalado(1:ntier) = 0.d0
!   jtier  = 4     
!   iballs = 2   ! check ADOs which have less than $iballs memory components within the longest $jtier memory time
   jtier  = nchklong
   iballs = nchkcount 
   ltmp1 = 0 
   dtmp1 = 0.d0
   do lni=nfirst(itier),nlast(itier)
      ntmp2 = 0
      do ni=1,itier-1
         idraw = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
         if (mcor(idraw) .eq. 1 .or. mcor(idraw) .gt. (ncor - jtier + 1)) then  ! works only for case of one drude pole 
            ntmp2 = ntmp2 + 1
         end if
      end do
      call cmaxmat(nrho, nrho, rho(1,1,lni), nrho, dtmp2)
      if (ntmp2 .lt. iballs) then    
         ltmp1 = ltmp1 + 1
         dtmp1 = max(dtmp1, dtmp2)
      end if
      lnumado(ntmp2 + 1) = lnumado(ntmp2 + 1) + 1
      dvalado(ntmp2 + 1) = dvalado(ntmp2 + 1) + dtmp2
   end do
   do ni=1,ntier
      dvalado(ni) = dvalado(ni) / dble(lnumado(ni))
   end do
   write(6,*)
   write(6,*)'checkrho: test ADOs with less than ', iballs  
   write(6,*)'checkrho: components of the most   ', jtier, ' longest memory time'
   write(6,*)'checkrho: # of ADOs                ', ltmp1
   write(6,*)'checkrho: largest value            ', dtmp1
   write(6,*)
   write(6,*)'checkrho: relative significance of ADOs with ni long memory components'
   write(6,*)'checkrho: ni, lnum_ado, dval_ado'
   do ni=1,ntier
      write(6,1005)ni-1, lnumado(ni), dvalado(ni)
   end do
   call flush(6)
!
   ltmp1 = 0
   dtmp3 = dtmp1
   do lni=nfirst(itier),nlast(itier)
      call cmaxmat(nrho, nrho, rho(1,1,lni), nrho, dtmp2)
      if (dtmp2 .gt. dtmp1) then
         ltmp1 = ltmp1 + 1
         dtmp3 = max(dtmp2, dtmp3)
      end if
   end do
   write(6,*)
   write(6,*)'checkrho: # of ADOs even larger    ', ltmp1
   write(6,*)'checkrho: largest value            ', dtmp3
   call flush(6)
!
end if
1004 format(I10, 1x, I2, 1x, I1, 1x, I1, 1x, I1, 1x, I2, 1x, I4)
1005 format('checkrho: ', I2, 1x, I9, 1x, e12.3e2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! check contribution from "identical" operators
!
if (lchksame .and. ntier .ge. 3) then
  dtmp1 = -1.d9
  jtier = 0
  lnj   = 0
  do itier=3,ntier
    do lni=nfirst(itier),nlast(itier)
      do ni=1,itier-1
        idraw = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        call look4operator(mpm(idraw), morbs(idraw), mspin(idraw), ilead(idraw), jopr(ni))
      end do
      ntmp1 = 1
      do ni=2,itier-1
        ntmp1 = ntmp1 * (jopr(ni) - jopr(ni-1))
      end do
      if (ntmp1 .eq. 0) then  ! identical operators detected
        do ni=1,nrho
          do nj=1,nrho
            if (cdabs(rho(ni,nj,lni)) .gt. dtmp1) then
              dtmp1 = cdabs(rho(ni,nj,lni))
              jtier = itier
              lnj   = lni
            end if
          end do
        end do
      end if
    end do
  end do
  write(6,*)
  write(6,*)' maximal element from ''identical'' operators at tier ', jtier
  write(6,*)' absolute value ', dtmp1
  write(6,*)' isgn,  iorbs,  ispin,  ialf,  imats '
  do ni=1,jtier-1
    idraw = indextable( ifirst(jtier) - 1 + (lnj - nfirst(jtier)) * (jtier - 1) + ni )
    write(6,1000)mpm(idraw), morbs(idraw), mspin(idraw), ilead(idraw), mcor(idraw)
  end do
  call flush(6)
end if
1000 format(1x, 4(I2, 2x), I6)
!
if (lchkmm .and. ntier .ge. 3) then
  write(6,*)
  write(6,*)' check level of significance for the cross-memory terms '
! except for the drude terms, check whether the cross-memory term, with |icor1-icor2| = idiffchk
! note that icor_k >= ndrude
  write(6,*)' idiffchk = ', idiffchk
  if (idiffchk .le. 0) then
    write(6,*)' error! non-positive idiffchk : ', idiffchk
    stop
  end if
  allocate(tmpmm0(idiffchk), tmpmcor(idiffchk,2), STAT=istat)
  tmpmm0(1:idiffchk) = 0.d0
  do itier=3,ntier
    do lni=nfirst(itier),nlast(itier)
      iballs = itier - 1
      do ni=1,iballs
        idraw = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
        mcortmp0(ni) = mcor(idraw)
      end do
      ntmp1 = 0
      ntmp3 = 0
      ntmp4 = 0
      do ni=1,iballs-1
        do nj=ni+1,iballs
          if (mcortmp0(ni) .gt. ndrude .and. mcortmp0(nj) .gt. ndrude) then 
            ntmp2 = abs(mcortmp0(ni) - mcortmp0(nj))
            if (ntmp1 .lt. ntmp2) then
              ntmp1 = ntmp2
              ntmp3 = min(mcortmp0(ni), mcortmp0(nj))
              ntmp4 = max(mcortmp0(ni), mcortmp0(nj))
            end if
          end if 
        end do
      end do
      if (ntmp1 .ge. 1 .and. ntmp1 .le. idiffchk) then
        dtmp1 = 0.d0
        do ni=1,nrho
          do nj=1,nrho
            dtmp1 = max(dtmp1, cdabs(rho(ni,nj,lni)))
          end do
        end do
        if (tmpmm0(ntmp1) .lt. dtmp1) then
          tmpmm0(ntmp1)     = dtmp1
          tmpmcor(ntmp1,1) = ntmp3
          tmpmcor(ntmp1,2) = ntmp4
        end if
      end if
    end do
  end do
  write(6,*)
  write(6,*)' max. element of ADO with max|icor1-icor2| = mcor0 : '
  write(6,*)'   mcor0    max. element   icor1    icor2            ' 
  do ni=1,idiffchk
    write(6,1001)ni, tmpmm0(ni), tmpmcor(ni,1), tmpmcor(ni,2)
  end do
  call flush(6)
  deallocate(tmpmm0, tmpmcor, STAT=istat)
end if
1001 format(1x, I4, 2x, e13.4e3, 2x, 2(I4, 2x))
!
! check number of nonzero elements in ADOs
!
if (lchksparse) then
   write(6,*)
   write(6,*)'checkrho: check number of nonzero elements in ADOs'
   write(6,*)'checkrho: full size = ', nrho**2
   nk = 0
   do nj=1,nrho
      do ni=1,nrho
         if (cdabs(hs(ni,nj)) .gt. dpico) then
             nk = nk + 1
         end if
      end do
   end do
   write(6,*)'checkrho: # of nonzero elements in Hamiltonian = ', nk
   call flush(6)
   !
   do itier=2,ntier
      !write(6,*)'checkrho: itier = ', itier
      nmax_element = 0
      allocate(nnonzero(nrho**2), STAT=istat)
      allocate(dnonzero(nrho,nrho), STAT=istat)
      do ni=1,nrho**2
         nnonzero(ni) = 0
      end do
      do nj=1,nrho
         do ni=1,nrho
            dnonzero(ni,nj) = 0.d0
         end do
      end do
      !ncount = 0
      nzero = 0
      do lni=nfirst(itier),nlast(itier)
         nk = 0    
         do nj=1,nrho
            do ni=1,nrho
               !if (cdabs(rho(ni,nj,lni)) .gt. dnano) then
               if (cdabs(rho(ni,nj,lni)) .gt. dpico) then
                  nk = nk + 1
               end if
               dnonzero(ni,nj) = max(dnonzero(ni,nj), cdabs(rho(ni,nj,lni)))
            end do
         end do
         if (nk .gt. 0) then
            nmax_element = max(nmax_element, nk)
            nnonzero(nk) = nnonzero(nk) + 1
         else 
            nzero = nzero + 1
         end if
         !ncount = ncount + 1
      end do
      write(6,*)'checkrho: itier= ', itier, ' nmax_element= ', nmax_element
      write(6,*)'nnonzero( ', 0, ' )= ', nzero
      do ni=1,nmax_element
         if (nnonzero(ni) > 0) then
            write(6,*)'nnonzero( ', ni, ' )= ', nnonzero(ni)
         end if
      end do
      !call amatout(nrho,nrho,dnonzero)
      write(6,*)
      deallocate(nnonzero, STAT=istat)
      deallocate(dnonzero, STAT=istat)
   end do
end if
call flush(6)
!

if (lchkbigado) then
   do itier=2,ntier
      lnk   = nfirst(itier)
      dtmp4 = 0.d0
      do lni=nfirst(itier),nlast(itier)
         call cmaxmat(nrho, nrho, rho(1,1,lni), nrho, dtmp1)
         if (dtmp1 .gt. dtmp4) then
            lnk = lni
            dtmp4 = dtmp1
         end if
      end do
      write(6,*)
      write(6,*)'checkrho: largest ADO at tier ', itier
      write(6,*)'lnk = ', lnk
      write(6,*)'largest value = ', dtmp4
      call cmatout(nrho, nrho, rho(1,1,lnk), dmtmp1)
      write(6,*)
      do ni=1,itier-1
         index_nk(ni) = indextable(ifirst(itier)-1+(lnk-nfirst(itier))*(itier-1)+ni)
         nj = index_nk(ni)
         write(6,*)'draw No. ', ni, ' = ', nj
         write(6,*)'isgn, icor ', mpm(nj), mcor(nj)
         call flush(6)
      end do
   end do
end if
!

!write(6,*)' ams (iorb=1,ispin=1) '
!call printams(1,1)
!write(6,*)' ams (iorb=1,ispin=2) '
!call printams(1,2)
!call flush(6)

!nj = 0
!cmtmp1 = czero
!cmtmp3 = czero
if (lchkbigado) then
   do itier=2,ntier
      write(6,*)
      write(6,*)'checkrho: show all ADOs at tier ', itier
      do lni=nfirst(itier),nlast(itier)
         do ni=1,itier-1
            index_nk(ni) = indextable(ifirst(itier)-1+(lni-nfirst(itier))*(itier-1)+ni)
         end do
         write(6,*)
         write(6,*)'lni = ', lni
         write(6,*)'drawer ', (index_nk(ni), ni=1,itier-1)
         if (lsimple) then
            if (lwalf) then
               write(6,*)'iopr ', (mopr(index_nk(ni)), ni=1,itier-1)
            else
               write(6,*)'iops ', (msopr(index_nk(ni)), ni=1,itier-1)
            end if
         else
            write(6,*)'iopr ', (mopr(index_nk(ni)), ni=1,itier-1)
            write(6,*)'iops ', (msopr(index_nk(ni)), ni=1,itier-1)
         end if
         write(6,*)'isgn ', (mpm(index_nk(ni)), ni=1,itier-1)
         write(6,*)'icor ', (mcor(index_nk(ni)), ni=1,itier-1)
         write(6,*)'ispn ', (mspin(index_nk(ni)), ni=1,itier-1)
         write(6,*)'iorb ', (morbs(index_nk(ni)), ni=1,itier-1)
         write(6,*)'iled ', (ilead(index_nk(ni)), ni=1,itier-1)
         call cmatout(nrho, nrho, rho(1,1,lni), dmtmp1)
         call flush(6)
   !
   !      if (itier .eq. 2 .and. mpm(index_nk(1)) .eq. 1 .and. &
   !          mspin(index_nk(1)) .eq. 1) then
   !          cmtmp1(1:nrho,1:nrho) = cmtmp1(1:nrho,1:nrho) + rho(1:nrho,1:nrho,lni)
   !          nj = nj + 1
   !      end if
   !      if (itier .eq. 2 .and. mpm(index_nk(1)) .eq. 1 .and. &
   !          mspin(index_nk(1)) .eq. 2) then
   !          cmtmp3(1:nrho,1:nrho) = cmtmp3(1:nrho,1:nrho) + rho(1:nrho,1:nrho,lni)
   !          nj = nj + 1
   !      end if
   !
      end do
   end do
end if
!
!write(6,*)
!write(6,*)'sum: ispin=1', nj
!call cmatout(nrho, nrho, cmtmp1, dmtmp1)
!if (nspin .eq. 2) then
!   write(6,*)'sum: ispin=2', nj
!   call cmatout(nrho, nrho, cmtmp3, dmtmp1)
!end if
!call flush(6)
!!
!call zamsmm1('l', 'n', 'n', 1, 1,  cunity, cmtmp1, nrho, czero,  cmtmp2, nrho)
!call zamsmm1('r', 'n', 'n', 1, 1, -cunity, cmtmp1, nrho, cunity, cmtmp2, nrho)
!if (nspin .eq. 2) then
!   call zamsmm1('l', 'n', 'n', 1, 2,  cunity, cmtmp3, nrho, cunity, cmtmp2, nrho)
!   call zamsmm1('r', 'n', 'n', 1, 2, -cunity, cmtmp3, nrho, cunity, cmtmp2, nrho)
!end if
!!
!write(6,*)
!write(6,*)' c * rho+ - rho+ c '
!call cmatout(nrho, nrho, cmtmp2, dmtmp1)
!call flush(6)
!
return
end subroutine checkrho
