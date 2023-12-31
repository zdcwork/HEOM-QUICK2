subroutine calcdhs(tt, rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8,     intent(in) :: tt
complex*16, intent(in) :: rhoinp(nrho,*)
integer                :: ni, nj
integer                :: iorbs, ispin
logical                :: lams0
real*8                 :: dshift, del, der
real*8                 :: dtmp1, dtmp2
real*8                 :: occd1, occd2
complex*16             :: zshift
!
integer :: jtimes
data jtimes /0/
!
lams0  = lprams
!lprams = .false.
!
jtimes = jtimes + 1
dhs(1:nrho,1:nrho) = czero
!
! early return if fixdot=TRUE
! (system Hamiltonian does not change during time evolution)
!
if (fixdot) goto 1999
!
if (ldfield_hub) then
   if (dfieldtype .eq. 1) then
      dtmp1 = 0.d0
      dtmp2 = 0.d0
      if (tt .lt. tson) then
         dtmp1 = dedot1 * dexp(-gamadot * tt)
         dtmp2 = dedot2 * dexp(-gamadot * tt)
      else if (tt .lt. tsoff) then
         dtmp1 = -dedot1 * dexp(-gamadot * (tt - tson))
         dtmp2 = -dedot2 * dexp(-gamadot * (tt - tson))
      else 
         dtmp1 = dedot1 * dexp(-gamadot * (tt - tsoff))
         dtmp2 = dedot2 * dexp(-gamadot * (tt - tsoff))
      end if 
!
      do ispin=1,nspin
         do iorbs=1,norbs
            cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, ispin), 0.d0)
            if (iorbs .eq. 1) then
               call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(dtmp1, 0.d0), cmtmp1, nrho,   &
                          cmtmp1, nrho, cunity, dhs, nrho)
            else if (iorbs .eq. 2) then
               call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(dtmp2, 0.d0), cmtmp1, nrho,   &
                          cmtmp1, nrho, cunity, dhs, nrho)
            end if
         end do
      end do
   else 
      write(6,*)'calcdhs: error! unknown dfieldtype = ', dfieldtype
      write(6,*)'calcdhs:        for Hubbard system   '
      stop
   end if
   goto 412
end if
!
! Note: dot level shift is defined to be reversed to the lead level shift. 
!       (see the minus sign below and compared to "getenergyshift")
!
if (nspin .eq. 1 .and. norbs .eq. 1 .and. dabs(dedot) .ge. dnano) then
  if (dfieldtype .eq. 0) then
    dshift = dedot * (1.d0 - dexp(-tt / aD))
  else if (dfieldtype .eq. 1) then
    dshift = dedot * dsin(wD * tt)
  else
    write(6,*) 
    write(6,*)' error! unknown dfieldtype in calcdhs', dfieldtype
    stop
  end if
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0)
  call zgemm('c', 'n', nrho, nrho, nrho, dcmplx(dshift, 0.d0), cmtmp1, nrho, &
             cmtmp1, nrho, cunity, dhs, nrho)
end if

!if (nspin .eq. 1 .and. norbs .eq.1)then
!   if (dfieldtype .eq. 0) then
!    dshift = dedot * (1.d0 - dexp(-tt / aD))
!   else if (dfieldtype .eq. 1) then
!    dshift = dedot * dsin(wD * tt)
!   cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c1
!   call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
!             czero, cmtmp3, nrho)                                     ! n1
!   else if (dfieldtype .eq. 4) then
!      if (tt .ge. tson .and. tt .lt. tsoff) then
!        dhs = dhs + dedot1 * cmtmp3 
!   end if
!end if
!
if (nspin .eq. 2 .and. norbs .eq. 1) then
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c_up
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 2), 0.d0) ! c_down
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
              czero, cmtmp3, nrho)                                   ! n1_up
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
              czero, cmtmp4, nrho)                                   ! n1_down
!
  if ( dabs(egate) .gt. dpico .and. &
        tt .gt. ton_gate .and. tt .lt. toff_gate ) then
      dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + egate * (cmtmp3(1:nrho,1:nrho) + cmtmp4(1:nrho,1:nrho))
  end if
  if (jtimes .eq. 1) then
     write(6,*)'calcdhs: egate = ', egate
     call flush(6)
  end if
! 
  if (gatetype .ne. 0) then
     if (gatetype .eq. 1) then
        dtmp1 = gpara1 + gpara2 * dcos(gpara3 * tt + gpara4)
     else
        write(6,*)
        write(6,*)'calcdhs: error! unknown gatetype ', gatetype
        stop
     end if
     dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + dtmp1 * (cmtmp3(1:nrho,1:nrho) + cmtmp4(1:nrho,1:nrho))
  end if
!
  if (lspin3d) then
!apply external gate voltage on the dot
     if (dfieldtype .eq. 0) then
         dedot1 = 0.d0
     else if (dfieldtype .eq. 1) then
         if (tt .ge. wdot1) dedot1=0.d0 ! turn off efield if tt >= wdot1
         dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) &
                            + dedot1 * (cmtmp3 + cmtmp4)
     else if (dfieldtype .eq. 2) then
         dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) &
                            + dedot1 * dcos(wdot1 * tt) * (cmtmp3 + cmtmp4)
     end if
!apply external 3d magnetic field on the dot
     if (lbfield3d) then
        if (tt .gt. dpico) then
           cmtmp1(1:nrho, 1:nrho) = -dbf3d(1,1) * sopr(1:nrho, 1:nrho, 1, 1) &
                                    -dbf3d(2,1) * sopr(1:nrho, 1:nrho, 2, 1) &
                                    -dbf3d(3,1) * sopr(1:nrho, 1:nrho, 3, 1)
           dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + cmtmp1(1:nrho, 1:nrho)
         end if
      end if
  else
     cmtmp3(1:nrho,1:nrho)=czero  
     if (dfieldtype .eq. 1 .and. dabs(tflip) .ge. dnano) then
        zshift = dcmplx(tflip * dsin(wflip * tt), 0.d0)
     else if (dfieldtype .eq. 2 .and. dabs(tupdn) .ge. dnano) then
        zshift = dcmplx(tupdn * dcos(wflip * tt), 0.d0)
     else if (dfieldtype .eq. 3 .and. dabs(tflip) .ge. dnano) then
        zshift = tflip *  cdexp(eye * wflip * tt)
     else
        zshift = czero
     end if
     call zgemm('c', 'n', nrho, nrho, nrho, zshift, cmtmp1, nrho, cmtmp2, nrho, &
                 czero, cmtmp3, nrho)
     do nj=1,nrho
       do ni=1,nrho 
          dhs(ni,nj) = dhs(ni,nj) + cmtmp3(ni,nj) + dconjg(cmtmp3(nj,ni))
       end do
     end do
  end if
end if
!
if (nspin .eq. 2 .and. norbs .eq. 2) then
   cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c1up
   cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 2), 0.d0) ! c1down
   cmtmp3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c2up
   cmtmp4(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 2), 0.d0) ! c2down
   call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
              czero, cmtmp5, nrho)                                     ! n1up
   call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
              czero, cmtmp6, nrho)                                     ! n1down
   call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp3, nrho, &
              czero, cmtmp7, nrho)                                     ! n2up
   call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp4, nrho, &
              czero, cmtmp8, nrho)                                     ! n2down
   
   if (dfieldtype .eq. 0) then
! nothing to do
   else if (dfieldtype .eq. 1) then
! step type gate voltage applied on each dot, turn off when tt > wdot1 (wdot2)
      if (tt .ge. wdot1) dedot1 = 0.d0  !turn off efield1 when tt>=wdot1
      if (tt .ge. wdot2) dedot2 = 0.d0  !turn off efield2 when tt>=wdot2
      dhs(1:nrho, 1:nrho) = dhs(1:nrho, 1:nrho) + dedot1 * (cmtmp5 + cmtmp6) &
                                                + dedot2 * (cmtmp7 + cmtmp8)
! interdot coupling tdot12 is activated and kept constant
      dtmp1 = 0.d0
      if (tt .gt. dpico .and. dabs(tdot12) .gt. dpico) dtmp1 = tdot12
      cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,1), 0.d0) ! c2up
      call zamsmm1('l', 'c', 'n', 1, 1, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                    czero, cmtmp2, nrho)
      cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,2), 0.d0) ! c2down
      call zamsmm1('l', 'c', 'n', 1, 2, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                    cunity, cmtmp2, nrho)
      do nj=1,nrho
         do ni=1,nrho
            dhs(ni,nj) = dhs(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
         end do
      end do
! interdot exchange jdot12 is activated and kept constant
      if (lspin3d .and. dabs(jdot12) .gt. dpico .and. tt .gt. dpico) then
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,1,1), nrho,  &
                    sopr(1,1,1,2), nrho, czero,  cmtmp1, nrho)
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,2,1), nrho,  &
                    sopr(1,1,2,2), nrho, cunity, cmtmp1, nrho)
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,3,1), nrho,  &
                    sopr(1,1,3,2), nrho, cunity, cmtmp1, nrho)
         dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + jdot12 * cmtmp1(1:nrho,1:nrho)
      end if
   else if (dfieldtype .eq. 2) then
! cosine type external field applied on each dot
      dhs(1:nrho, 1:nrho) = dhs(1:nrho, 1:nrho) + dedot1 * dcos(wdot1 * tt) * (cmtmp5 + cmtmp6) &
                                                + dedot2 * dcos(wdot2 * tt) * (cmtmp7 + cmtmp8)

! interdot coupling tdot12 is activated and kept constant
      dtmp1 = 0.d0
      if (tt .gt. dpico .and. dabs(tdot12) .gt. dpico) dtmp1 = tdot12
      cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,1), 0.d0) ! c2up
      call zamsmm1('l', 'c', 'n', 1, 1, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                    czero, cmtmp2, nrho)
      cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,2,2), 0.d0) ! c2down
      call zamsmm1('l', 'c', 'n', 1, 2, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                    cunity, cmtmp2, nrho)
      do nj=1,nrho
         do ni=1,nrho
            dhs(ni,nj) = dhs(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
         end do
      end do
! interdot exchange jdot12 is activated and kept constant
      if (lspin3d .and. dabs(jdot12) .gt. dpico .and. tt .gt. dpico) then
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,1,1), nrho,  &
                    sopr(1,1,1,2), nrho, czero,  cmtmp1, nrho)
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,2,1), nrho,  &
                    sopr(1,1,2,2), nrho, cunity, cmtmp1, nrho)
         call zgemm('n', 'n', nrho, nrho, nrho, cunity, sopr(1,1,3,1), nrho,  &
                    sopr(1,1,3,2), nrho, cunity, cmtmp1, nrho)
         dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + jdot12 * cmtmp1(1:nrho,1:nrho)
      end if
   else if (dfieldtype .eq. 3) then
! other situations??
   end if
   if (lspin3d .and. lbfield3d) then
      if (tt .gt. dpico) then
!apply external 3d magnetic field on the dot
         cmtmp1(1:nrho,1:nrho) =    -dbf3d(1,1) * sopr(1:nrho,1:nrho,1,1)  &
                                    -dbf3d(2,1) * sopr(1:nrho,1:nrho,2,1)  &
                                    -dbf3d(3,1) * sopr(1:nrho,1:nrho,3,1)  &
                                    -dbf3d(1,2) * sopr(1:nrho,1:nrho,1,2)  &
                                    -dbf3d(2,2) * sopr(1:nrho,1:nrho,2,2)  &
                                    -dbf3d(3,2) * sopr(1:nrho,1:nrho,3,2)
            dhs(1:nrho,1:nrho) = dhs(1:nrho,1:nrho) + cmtmp1(1:nrho,1:nrho)
       end if
   end if
end if
!
if (nspin .eq. 1 .and. norbs .eq. 2) then
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c1
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp3, nrho)                                     ! n1
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)                                     ! n2
  if (dfieldtype .eq. 2) then
    dhs = dhs + dedot1 * dcos(wdot1 * tt) * cmtmp3 + &
                dedot2 * dcos(wdot2 * tt) * cmtmp4
  else if (dfieldtype .eq. 4) then
    if (tt .ge. tson .and. tt .lt. tsoff) then
        dhs = dhs + dedot1 * cmtmp3 + dedot2 * cmtmp4
    end if   
  end if
end if
!
if (nspin .eq. 1 .and. norbs .eq. 2 .and. dabs(pha) .ge. dnano)then
  cmtmp1(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 1, 1), 0.d0) ! c1
  cmtmp2(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, 2, 1), 0.d0) ! c2
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, &
             czero, cmtmp3, nrho)                                              ! n1
  call zgemm('c', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, cmtmp2, nrho, &
             czero, cmtmp4, nrho)                                              ! n2
  call zhemm('r','u', nrho, nrho, cunity, rho(1,1,1), nrho, cmtmp3, nrho, czero, &
              cmtmp1, nrho) 
  call zhemm('r','u', nrho, nrho, cunity, rho(1,1,1), nrho, cmtmp4, nrho, czero, & 
              cmtmp2, nrho) 
  occd1 = 0.d0
  occd2 = 0.d0
  do nj=1,nrho
    occd1 = occd1 + dble(cmtmp1(nj,nj))
    occd2 = occd2 + dble(cmtmp2(nj,nj))
  end do 
  dhs = dhs + pha * occd1 * cmtmp3 + pha * occd2 * cmtmp4
end if
!=========================================================
!
412 continue
!
if (lbfield .and. nspin .gt. 1 .and. .not. lspin3d) then
   cmtmp1 = czero
   cmtmp2 = czero
   do iorbs=1,norbs
      cmtmp3(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 1), 0.d0)
      call zamsmm1('l', 'c', 'n', iorbs, 1, cunity, cmtmp3, nrho, cunity, cmtmp1, nrho)
      cmtmp4(1:nrho, 1:nrho) = dcmplx(amsori(1:nrho, 1:nrho, iorbs, 2), 0.d0)
      call zamsmm1('l', 'c', 'n', iorbs, 2, cunity, cmtmp4, nrho, cunity, cmtmp2, nrho)
   end do
   cmtmp1 = cmtmp1 - cmtmp2          ! n_up - n_down
   do nj=1,nrho
      do ni=1,nrho
         dhs(ni,nj) = dhs(ni,nj) - dbfield * (1.d0 - dexp(-tt / tbfield)) * cmtmp1(ni,nj) * 0.5d0
      end do
   end do
end if
!
! Apply projector if needed: P*H*P
!
if (lprgen .and. ldopro) then
    call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, dhs, nrho, czero, cmtmp1, nrho)
    call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, czero, dhs, nrho)
end if
!
1999 continue
!
lprams = lams0
!
return
end subroutine calcdhs
