subroutine checkdos(isgnb, ispinb, iorbsb, rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)     :: isgnb, ispinb, iorbsb
complex*16, intent(in)  :: rhoinp(nrho,nrho,*)
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs, malf, itype
integer                 :: na, nb, nc, mi, mj
integer*8               :: lni, lnj, lnk, lnm, lnl
integer                 :: isgn, iorbs, ispin, ialf, iorbs2
logical                 :: lherm, lsame, lcomp
character*1             :: sidea, transa, transb
real*8                  :: dtmp1
complex*16, allocatable :: work(:), cexphs(:,:), cexpheff(:,:)
real*8,     allocatable :: tmpval(:), rwork(:)
complex*16, allocatable :: hybrid(:,:,:)
!
namelist / coupling / readcp, readmat, lrcmplx
!
write(6,*)
write(6,*)'checkdos: isgnb, iorbsb, ispinb ', isgnb, iorbsb, ispinb
call flush(6)
!
! cexphs = exp(-H)
!
allocate(cexphs(nrho,nrho), cexpheff(nrho,nrho), STAT=istat)
!
rewind(5)
do ni=1,5
   read(5,*)nj
end do
read(5,*)malf
!
allocate(hybrid(norbs,norbs,malf), STAT=istat)
!
cmtmp2(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
lwork = nrho**2
allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
if (info .ne. 0) then
   write(6,*)
   write(6,*)'checkdos: error! diagnalization failed 1 ', info
   stop
end if
cmtmp3(1:nrho,1:nrho) = czero
do ni=1,nrho
   cmtmp3(ni,ni) = dcmplx( dexp(-tmpval(ni)), 0.d0 )
end do
deallocate(work, rwork, tmpval, STAT=istat)
call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
           czero, cmtmp4, nrho)
call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
           czero, cmtmp3, nrho)
cexphs(1:nrho,1:nrho) = cmtmp3(1:nrho,1:nrho)
!
hybrid(1:norbs,1:norbs,1:malf) = czero
do ialf=1,malf
   do ni=1,norbs
      hybrid(ni,ni,ialf) = dcmplx(linewidth(ialf,ni) / hbar, 0.d0)
   end do
end do
if (lequileads) then
   do ialf=2,malf
      do ni=1,norbs
         hybrid(ni,ni,1) = hybrid(ni,ni,1) + hybrid(ni,ni,ialf)
      end do
   end do
end if
!
if (offcor) then
   if (megaflux) then
      write(6,*)
      write(6,*)'checkdos: warning, offcor with megaflux case not available yet '
      stop
   else
      cmtmp2(1:norbs,1:norbs) = czero
      hybrid(1:norbs,1:norbs,1:malf) = czero
!
      readcp  = .false. 
      readmat = .false. 
      lrcmplx = .false.
      rewind(5)
      read(5, coupling, end=200)     
      200 continue
      if (readcp) then 
         do ialf=1,malf
            dmtmp3(1:norbs,1:norbs) = 0.d0
            if (readmat) then 
               do ni=1,norbs
                  if (lrcmplx) then
                     read(5,*,iostat=istat) (dmtmp2(ni,nj), nj=1,norbs), (dmtmp3(ni,nj), nj=1,norbs)
                  else
                     read(5,*,iostat=istat) (dmtmp2(ni,nj), nj=1,norbs)
                  end if
                  if (istat .ne. 0) then 
                     write(6,*)'checkdos: error reading linewidth, ialf = ', ialf, ' iorbs = ', ni
                     stop
                  end if
               end do
            else 
               if (lrcmplx) then
                  read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs), ((dmtmp3(ni,nj), ni=1,norbs), nj=1,norbs)
               else
                  read(5,*,iostat=istat) ((dmtmp2(ni,nj), ni=1,norbs), nj=1,norbs)
               end if
               if (istat .ne. 0) then
                  write(6,*)'checkdos: error reading linewidth matrix for ialf = ', ialf
                  stop
               end if
            end if
            hybrid(1:norbs,1:norbs,ialf) = dcmplx(dmtmp2(1:norbs,1:norbs), dmtmp3(1:norbs,1:norbs)) / hbar
            cmtmp2(1:norbs,1:norbs) = cmtmp2(1:norbs,1:norbs) + dcmplx(dmtmp2(1:norbs,1:norbs), dmtmp3(1:norbs,1:norbs)) / hbar
         end do  
      else
         do ialf=1,malf
            do ni=1,norbs
               hybrid(ni,ni,ialf) = dcmplx(linewidth(ialf,ni) / hbar, 0.d0)
            end do
         end do
      end if
!
      if (lequileads) then
         do ialf=2,malf
            do ni=1,norbs
               hybrid(ni,ni,1) = hybrid(ni,ni,1) + hybrid(ni,ni,ialf)
            end do
         end do
      end if
!
! extract off-diagonal coupling
      cmtmp1(1:norbs,1:norbs) = czero
      do ni=1,norbs
         do nj=1,ni-1
            if (cdabs(cmtmp2(ni,nj)) .gt. dpico) then
               cmtmp1(ni,nj) = cmtmp2(ni,nj) * eye    ! (or -eye?) sign chosen arbitrarily
               cmtmp1(nj,ni) = dconjg(cmtmp1(ni,nj))  ! impose Hermicity
            end if
         end do
      end do  
! construct effective Hamiltonian corresponding to off-diagonal coupling
      cmtmp2(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
      do ni=1,norbs
         do nj=1,norbs
            if (cdabs(cmtmp1(ni,nj)) .gt. dpico) then 
               do ispin=1,nspin
                  cmtmp3(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,nj,ispin), 0.d0)
                  call zamsmm1('L', 'C', 'N', ni, ispin, cmtmp1(ni,nj), cmtmp3, nrho, cunity, cmtmp2, nrho)
               end do
            end if
         end do
      end do
      call checkhermicity(cmtmp2, nrho, cmtmp3, nrho, lherm, dtmp1)
      if (.not. lherm) then
         write(6,*)
         write(6,*)'checkdos: error, effective H is not hermitian '
         stop
      end if
! diagonalize Heff
      lwork = nrho**2
      allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
      call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
      if (info .ne. 0) then
         write(6,*)
         write(6,*)'checkdos: error! diagnalization failed 2 ', info
         stop
      end if
      cmtmp3(1:nrho,1:nrho) = czero
      do ni=1,nrho
         cmtmp3(ni,ni) = dcmplx( dexp(-tmpval(ni)), 0.d0 )
      end do
      deallocate(work, rwork, tmpval, STAT=istat)
      call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
                 czero, cmtmp4, nrho)
      call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
                 czero, cmtmp3, nrho)
      cexpheff(1:nrho,1:nrho) = cmtmp3(1:nrho,1:nrho)
   end if 
end if
!
!call checkhermicity(rhoinp(1,1,1), nrho, cmtmp1, nrho, lherm, dtmp1)
!if (lherm) then
!   write(6,*)'checkdos: rhoinp is Hermitian '
!else 
!   write(6,*)'checkdos: rhoinp is NOT Hermitian '
!end if
!
lsame = .true.
lcomp = .true.
lnk = 0
write(6,*)
if (.not. offcor) then
   cmtmp3(1:nrho,1:nrho) = cexphs(1:nrho,1:nrho)
   write(6,*)'checkdos: term_r = exp(-H) * (B + B+) * exp(-H) '
   write(6,*)'checkdos: term_h = exp(-H) * (ams*rho + rho*ams) '
else 
   cmtmp3(1:nrho,1:nrho) = cexpheff(1:nrho,1:nrho)
   write(6,*)'checkdos: term_r = exp(-Heff) * (B + B+) * exp(-Heff) '
   write(6,*)'checkdos: term_h = exp(-Heff) * (ams*rho + rho*ams)   '
end if
call flush(6)
!
sidea  = 'L'
transa = 'N'
transb = 'N'
if (isgnb .eq. 1) then
   transa = 'C'
end if
call zamsmm1(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp1, nrho)
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp1, nrho, czero, cmtmp5, nrho)
cmtmp1(1:nrho,1:nrho) = cmtmp5(1:nrho,1:nrho)  ! exp(-H)*B*exp(-H) or exp(-Heff)*B*exp(-Heff)
!
sidea  = 'R'
transa = 'C'
transb = 'C'
if (isgnb .eq. 1) then
   transa = 'N'
end if
call zamsmm1(sidea, transa, transb, iorbsb, ispinb, cunity, cmtmp3, nrho, czero, cmtmp2, nrho)
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp3, nrho, czero, cmtmp5, nrho)
cmtmp2(1:nrho,1:nrho) = cmtmp5(1:nrho,1:nrho)  ! exp(-H)*(B+)*exp(-H) or exp(-Heff)*(B+)*exp(-Heff)
!
dmtmp1 = 0.d0
dmtmp2 = 0.d0
do nj=1,nrho
   do ni=1,nrho 
      if ( cdabs(cmtmp1(ni,nj)) .gt. dpico .or. cdabs(cmtmp2(ni,nj)) .gt. dpico ) then 
          dmtmp1(ni,nj) = 1.d0
      end if
      if ( cdabs(rhoinp(ni,nj,1)) .gt. dpico ) then 
          dmtmp2(ni,nj) = 1.d0
      end if
   end do
end do
dmtmp1 = dmtmp1 - dmtmp2
call dmaxmat(nrho, nrho, dmtmp1, nrho, dtmp1)
write(6,*)
!write(6,*)' max_diff = ', dtmp1
if (dabs(dtmp1) .gt. dpico) then
   lcomp = .true.
   check9: do nj=1,nrho
      do ni=1,nrho
         if ( cdabs(rhoinp(ni,nj,1)) .gt. dpico ) then
            if ( cdabs(cmtmp1(ni,nj)) .lt. dpico .and. cdabs(cmtmp2(ni,nj)) .lt. dpico ) then 
               lcomp = .false.
               exit check9
            end if
         end if
      end do
   end do check9
   write(6,*)'checkdos: topology of rhoinp differs from term_r '
   if (lcomp) then
      write(6,*)'checkdos: rhoinp is compatiable with term_r '
   else
      write(6,*)'checkdos: rhoinp is not compatiable with term_r '
   end if
else 
   write(6,*)'checkdos: rhoinp has same topology as term_r '
end if
call flush(6)


!deallocate(cexphs, cexpheff, STAT=istat)
!deallocate(hybrid, STAT=istat)
!return 

!
write(6,*)
lsame = .true.
lcomp = .true.
lnk = 0
if (.not. offcor) then
   cmtmp3(1:nrho,1:nrho) = cexphs(1:nrho,1:nrho)
else 
   cmtmp3(1:nrho,1:nrho) = cexpheff(1:nrho,1:nrho)
end if
do itier=2,ntier
   do lni=nfirst(itier),nlast(itier)
!
      lnl = index_coef_ref(lni)
      na = 0
      dmtmp1(1:nrho,1:nrho) = 0.d0
      do nj=1,nrho
         do ni=1,nrho
            if ( cdabs(rhoinp(ni,nj,lni)) .gt. dpico ) then
               dmtmp1(ni,nj) = 1.d0
               na = na + 1
            end if
         end do
      end do
!
      nb = 0
      dmtmp2(1:nrho,1:nrho) = 0.d0
      do ni=1,itier-1
         nk    = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
         isgn  = mpm(nk)
         iorbs = morbs(nk)
         ispin = mspin(nk)
         ialf  = ilead(nk)
         lnj   = indexcoef(lnl)
         itype = itypecoef(lnl)
         lnl   = lnl + 1
!
         cmtmp2(1:nrho,1:nrho) = rhoinp(1:nrho,1:nrho,lnj)
         sidea  = 'L'
         transa = 'N'
         transb = 'N'
         if (isgn .eq. 1) then
            transa = 'C'
         end if
         if (itype .ne. 0) then
            transb = 'C'
         end if   

         if (.not. offcor) then
            if ( cdabs(hybrid(iorbs,iorbs,ialf)) .lt. dpico ) cycle
! cmtmp7 : ams*rho
            sidea  = 'L'
            call zamsmm1(sidea, transa, transb, iorbs, ispin, cunity, cmtmp2, nrho, czero, cmtmp7, nrho)
! cmtmp8 : rho*ams
            sidea  = 'R'
            call zamsmm1(sidea, transa, transb, iorbs, ispin, cunity, cmtmp2, nrho, czero, cmtmp8, nrho)
!! exp(-H)*B*exp(-H) * (ams*rho) and exp(-H)*B*exp(-H) * (rho*ams)
!            call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp7, nrho, czero, cmtmp5, nrho)
!            call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp8, nrho, czero, cmtmp6, nrho)
!! exp(-H)*B+*exp(-H) * (ams*rho) and exp(-H)*B+*exp(-H) * (rho*ams)
!            call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp7, nrho, czero, cmtmp3, nrho)
!            call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp8, nrho, czero, cmtmp4, nrho)
!
            call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp7, nrho, czero, cmtmp5, nrho)
            call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp8, nrho, czero, cmtmp6, nrho)

            do mj=1,nrho
               do mi=1,nrho
                  if ( dabs(dmtmp2(mi,mj)) .gt. dpico ) cycle
!                  if ( cdabs(cmtmp5(mi,mj)) .gt. dpico .or. cdabs(cmtmp6(mi,mj)) .gt. dpico .or.  &
!                       cdabs(cmtmp3(mi,mj)) .gt. dpico .or. cdabs(cmtmp4(mi,mj)) .gt. dpico ) then 
                  if ( cdabs(cmtmp5(mi,mj)) .gt. dpico .or. cdabs(cmtmp6(mi,mj)) .gt. dpico ) then
                     dmtmp2(mi,mj) = 1.d0
                  end if
               end do
            end do
!
         else 
            do iorbs2=1,norbs
               if ( cdabs(hybrid(iorbs,iorbs2,ialf)) .lt. dpico ) cycle
               sidea  = 'L'
               call zamsmm1(sidea, transa, transb, iorbs2, ispin, cunity, cmtmp2, nrho, czero, cmtmp7, nrho)
               sidea  = 'R'
               call zamsmm1(sidea, transa, transb, iorbs2, ispin, cunity, cmtmp2, nrho, czero, cmtmp8, nrho)
!
!               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp7, nrho, czero, cmtmp5, nrho)
!               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp8, nrho, czero, cmtmp6, nrho)
!               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp7, nrho, czero, cmtmp3, nrho)
!               call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp8, nrho, czero, cmtmp4, nrho)
               call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp7, nrho, czero, cmtmp5, nrho)
               call zhemm('l', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp8, nrho, czero, cmtmp6, nrho)
!
               do mj=1,nrho
                  do mi=1,nrho
                     if ( dabs(dmtmp2(mi,mj)) .gt. dpico ) cycle
!                     if ( cdabs(cmtmp5(mi,mj)) .gt. dpico .or. cdabs(cmtmp6(mi,mj)) .gt. dpico .or. &
!                          cdabs(cmtmp3(mi,mj)) .gt. dpico .or. cdabs(cmtmp4(mi,mj)) .gt. dpico ) then
                     if ( cdabs(cmtmp5(mi,mj)) .gt. dpico .or. cdabs(cmtmp6(mi,mj)) .gt. dpico ) then
                        dmtmp2(mi,mj) = 1.d0
                     end if
                  end do
               end do 
            end do ! end of loop over iorbs2
         end if
!
      end do ! end of loop over ni
!
      do nj=1,nrho
         do ni=1,nrho
            if ( dabs(dmtmp2(ni,nj)) .gt. dpico ) then
               nb = nb + 1
            end if
         end do
      end do
      lnk = lnk + nb
!
      dmtmp3(1:nrho,1:nrho) = dmtmp1(1:nrho,1:nrho) - dmtmp2(1:nrho,1:nrho)
      call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
      if (dabs(dtmp1) .gt. dpico) then
!
         lcomp = .true.
         check1: do nj=1,nrho
            do ni=1,nrho
               if (dabs(dmtmp1(ni,nj)) .gt. dpico) then
                  if (dabs(dmtmp2(ni,nj)) .lt. dpico) then 
                     lcomp = .false.
                     exit check1
                  end if
               end if
            end do
         end do check1
!
         if (.not. lcomp) then
            lsame = .false.
            write(6,*)'checkdos: ADO # ', lni, ' differs from term_h '
            write(6,*)'checkdos: nonzeros in ADO ', na, ' in term_h ', nb
            call flush(6)
         end if
      end if
!
   end do ! end of loop over lni
end do ! end of loop over itier
!
lni = nrho**2 * (nunk - 1)
write(6,*)
write(6,*)'checkdos: total ADO elements     = ', lni
write(6,*)'checkdos: nonzero ADO elements   = ', lnk
write(6,*)'checkdos: sparsity ratio         = ', dble(lnk) / dble(lni)
write(6,*)
if (lsame) then
   write(6,*)'checkdos: topology of all ADOs same as term_h '
else 
   write(6,*)'checkdos: topology of some ADOs differ from term_h '
end if
call flush(6)
!
deallocate(cexphs, cexpheff, STAT=istat)
deallocate(hybrid, STAT=istat)
end subroutine checkdos
