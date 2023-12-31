subroutine genhs
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer       :: istat, info, ni, nj, npt
integer       :: iorbs0, jorbs0, ispin0, jspin0, itype    ! itype = 0, U-term takes the Coulomb repulsion form
                                                          !    0.5 * \sum_{ij} U_{ij} n_i n_j, with n_i = \sum_s n_{is}
                                                          ! itype = 1, U-term takes the Hubbard interaction form
                                                          !    0.5 * \sum_{ij} U_{ij} \sum_{ss'} a+_{is}a+_{js'} a_{js'}a_{is} 
logical       :: lams0                                                         
real*8        :: dtmp1
character*120 :: cline
!
! general input format for system Hamiltonian
! 
! $hamil_sys  lgenhs=.true. $end
!
! hsys
! edot   ispin iorbs       para_edot             :  on-site (level) energy
! ecoup  ispin iorbs jorbs para_ecoup            :  inter-level coupling   (iorbs>jorbs is required)
! udot   itype iorbs       para_udot             :  on-site U term
! ucoup  itype iorbs jorbs para_ucoup            :  inter-level U coupling (iorbs>jorbs is required)
!
!
! for the single-site Anderson model, the e-e interaction term is U n_u n_d = 0.5 * U (n n - n_u - n_d)
! to recover this interaction term, if udot with itype=0 and the same U parameter is used, the on-site terms need
! to be modified, as follows
! edot 1 1  eup   - 0.5 * U
! edot 2 1  edown - 0.5 * U
!
if (.not. lgenhs) then
   write(6,*)
   write(6,*)'genhs: error calling subroutine'
   stop
end if
!
lams0  = lprams 
!lprams = .false.
!
hs(1:nrho,1:nrho) = czero
!
rewind(5)
find_sys: do
   read(5, '(A120)') cline
   istat = index(cline, 'hsys')
   if (istat > 0) then
      write(6,*)
      write(6,*)'genhs: start reading parameters for system Hamiltonian'
      exit find_sys
   end if
end do find_sys
!
if (istat > 0) then
   read_sys: do
      read(5, '(A120)') cline
      write(6,'(A120)') cline
!
      if ( index(cline, 'edot') .gt. 0 ) then  ! on-site energy
         npt = index(cline, 'edot') + 4
         read(cline(npt:120),*,iostat=istat) ispin0, iorbs0, dtmp1
         if (istat == 0) then
            if (ispin0 > nspin .or. iorbs0 > norbs) then
               write(6,*)'genhs: error! ispin, nspin, iorbs, norbs', ispin0, nspin, iorbs0, norbs
               stop
            end if
            dtmp1 = dtmp1 / hbar
            cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs0,ispin0), 0.d0)
            call zamsmm1('l', 'c', 'n', iorbs0, ispin0, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                        czero, cmtmp2, nrho)
            hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + cmtmp2(1:nrho,1:nrho)
         end if
!
      else if ( index(cline, 'ecoup') .gt. 0 ) then  ! hopping between sites (spin-conserved, iorbs0>jorbs0 reqired)
         npt = index(cline, 'ecoup') + 5
         read(cline(npt:120),*,iostat=istat) ispin0, iorbs0, jorbs0, dtmp1
         if (ispin0 > nspin .or. iorbs0 > norbs .or. jorbs0 > norbs) then
            write(6,*)'genhs: error! ispin, nspin, iorbs, jorbs, norbs', ispin0, nspin, iorbs0, jorbs0, norbs
            stop
         end if
         if (iorbs0 .le. jorbs0) then
            write(6,*)'genhs: error! iorbs > jorbs is required for ecoup', iorbs0, jorbs0
            stop
         end if
         if (istat == 0) then
            dtmp1 = dtmp1 / hbar
            cmtmp1(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,jorbs0,ispin0), 0.d0)
            call zamsmm1('l', 'c', 'n', iorbs0, ispin0, dcmplx(dtmp1,0.d0), cmtmp1, nrho, &
                        czero, cmtmp2, nrho)
            do nj=1,nrho
               do ni=1,nrho
                  hs(ni,nj) = hs(ni,nj) + cmtmp2(ni,nj) + dconjg(cmtmp2(nj,ni))
               end do
            end do
         end if
!
      else if ( index(cline, 'udot') .gt. 0 ) then  ! on-site U interaction (Coulomb or Hubbuard type)
         npt = index(cline, 'udot') + 4
         read(cline(npt:120),*,iostat=istat) itype, iorbs0, dtmp1
         if (iorbs0 > norbs) then
            write(6,*)'genhs: error! iorbs, norbs', iorbs0, norbs
            stop
         end if
         if (itype .ne. 0 .and. itype .ne. 1) then
            write(6,*)'genhs: error! unknown itype (0 for Coulomb and 1 for Hubbard)', itype
            stop
         end if
         if (istat == 0) then
            dtmp1 = dtmp1 / hbar
            if (itype == 0) then
               cmtmp3(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  cmtmp2(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs0,ispin0), 0.d0)
                  call zamsmm1('l', 'c', 'n', iorbs0, ispin0, cunity, cmtmp2, nrho, cunity, cmtmp3, nrho)
               end do
               call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.5d0 * dtmp1, 0.d0), cmtmp3, nrho, cmtmp3, nrho, &
                          czero, cmtmp1, nrho)         
            else if (itype == 1) then
               cmtmp1(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  cmtmp2(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs0,ispin0), 0.d0)
                  call zamsmm1('l', 'c', 'n', iorbs0, ispin0, cunity, cmtmp2, nrho, cunity, cmtmp1, nrho)
               end do
               cmtmp2(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  call zamsmm1('r', 'n', 'n', iorbs0, ispin0, cunity, cmtmp1, nrho,  czero, cmtmp3, nrho)
                  call zamsmm1('l', 'c', 'n', iorbs0, ispin0, cunity, cmtmp3, nrho, cunity, cmtmp2, nrho)
               end do
               cmtmp1(1:nrho,1:nrho) = 0.5d0 * dtmp1 * cmtmp2(1:nrho,1:nrho)
            else 
               write(6,*)'genhs: error! unknown itype for udot', itype
               stop
            end if
            hs(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + cmtmp1(1:nrho,1:nrho)
         end if
 ! 
      else if ( index(cline, 'ucoup') .gt. 0 ) then ! inter-level U-coupling (iorbs0>jorbs0 is required)
         npt = index(cline, 'ucoup') + 5
         read(cline(npt:120),*,iostat=istat) itype, iorbs0, jorbs0, dtmp1
         if (iorbs0 > norbs .or. jorbs0 > norbs) then
            write(6,*)'genhs: error! iorbs, jorbs, norbs', iorbs0, jorbs0, norbs
            stop
         end if
         if (iorbs0 .le. jorbs0) then
            write(6,*)'genhs: error! iorbs > jorbs is required for ucoup', iorbs0, jorbs0
            stop
         end if
         if (itype .ne. 0 .and. itype .ne. 1) then
            write(6,*)'genhs: error! unknown itype (0 for Coulomb and 1 for Hubbard)', itype
            stop
         end if
         if (istat == 0) then
            dtmp1 = dtmp1 / hbar
            if (itype == 0) then
               cmtmp3(1:nrho,1:nrho) = czero
               cmtmp4(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  cmtmp2(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,jorbs0,ispin0), 0.d0)
                  call zamsmm1('l', 'c', 'n', jorbs0, ispin0, cunity, cmtmp2, nrho, cunity, cmtmp3, nrho)
                  cmtmp2(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,iorbs0,ispin0), 0.d0)
                  call zamsmm1('l', 'c', 'n', iorbs0, ispin0, cunity, cmtmp2, nrho, cunity, cmtmp4, nrho)
               end do
               call zgemm('n', 'n', nrho, nrho, nrho, dcmplx(0.5d0 * dtmp1, 0.d0), cmtmp4, nrho, cmtmp3, nrho, &
                          czero, cmtmp1, nrho)         
!
            else if (itype == 1) then
               cmtmp1(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  cmtmp2(1:nrho,1:nrho) = dcmplx(amsori(1:nrho,1:nrho,jorbs0,ispin0), 0.d0)
                  call zamsmm1('l', 'c', 'n', jorbs0, ispin0, cunity, cmtmp2, nrho, cunity, cmtmp1, nrho)
               end do
               cmtmp2(1:nrho,1:nrho) = czero
               do ispin0=1,nspin
                  call zamsmm1('r', 'n', 'n', iorbs0, ispin0, cunity, cmtmp1, nrho,  czero, cmtmp3, nrho)
                  call zamsmm1('l', 'c', 'n', iorbs0, ispin0, cunity, cmtmp3, nrho, cunity, cmtmp2, nrho)
               end do
               cmtmp1(1:nrho,1:nrho) = 0.5d0 * dtmp1 * cmtmp2(1:nrho,1:nrho)
!
            else 
               write(6,*)'genhs: error! unknown itype for ucoup', itype
               stop
            end if
            do nj=1,nrho
               do ni=1,nrho
                  hs(ni,nj) = hs(ni,nj) + cmtmp1(ni,nj) + dconjg(cmtmp1(nj,ni))
               end do
            end do
         end if
!
      else 
         istat = -1
      end if 
!
      if (istat .ne. 0) then  
          exit read_sys
      end if
!
   end do read_sys
end if
!
1999 continue
!
lprams = lams0
!
return
end subroutine genhs
