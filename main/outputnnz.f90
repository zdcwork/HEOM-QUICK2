 subroutine outputnnz      
 use matmod   
 use tmpmatmod
 use sparsemod
 use random                
 implicit none              
 include '../include/sizes' 
 include '../include/common'
 
integer                 :: ni, nj, nk, info, istat, lwork, itier, jtier, iballs, malf
integer                 :: na, nb, nc, mi, mj, nnz, isame, icomp
integer*8               :: lni, lnj, lnk, lnm, lnl, lunk_last, lstart
integer*1               :: itype, fact
integer                 :: isgn, iorbs, ispin, ialf, iorbs2, nrho2, icor
integer                 :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
integer                 :: mindex(MAXTIER)
logical                 :: lsame, lexist, lexist1, lexist2, ltmp1, ltmp2, ltmp3, ltmp4
integer,    allocatable    :: nvec1(:), nvec2(:)
complex*16, allocatable     :: cvec1(:)
 
 nrho2 = nrho**2
 lexist = .false.
 inquire(file='sparse_index.data', exist=lexist1, err=999)   ! unit = 52
 999 continue
 inquire(file='sparse_info.data', exist=lexist2, err=998)    ! unit = 53
 998 continue
 if (lexist1 .and. lexist2) then
  lexist = .true.
 else
  return
 end if
 write(6,*)'outputnnz:',lexist
 
 if (lexist) then
   write(6,*)
   write(6,*)'outputnnz: <sparse_index.data> and <sparse_info.data> found '
   write(6,*)'outputnnz: start reading '
   open(unit=52, file='sparse_index.data', form='binary', status='old')
   rewind(52)
   read(52)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
   if (ntmp1 .ne. ntier .or. ntmp2 .ne. ncor .or. ntmp3 .ne. norbs .or.      &
       ntmp4 .ne. nspin .or. ntmp5 .ne. nalf .or. ntmp6 .ne. numfff .or.     &
       ntmp7 .ne. ntier0 .or. ntmp8 .ne. ndrawer_slow .or. ntmp9 .ne. ncor_slow) then
      write(6,*)
      write(6,*)'outputnnz: <sparse_index.data> incompatible with present job '
      write(6,*)'outputnnz: abort reading '
      write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
      write(6,*)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow
      call flush(6)
      close(52)
      lexist = .false.
      return
   end if
   read(52)lni, ltmp1, ltmp2, ltmp3, ltmp4
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen) then
      write(6,*)'outputnnz: <sparse_index.data> incompatible with present job '
      write(6,*)'outputnnz: abort reading '
      write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4
      write(6,*)nunk, offcor, lequileads, lsimple, lscreen
      call flush(6)
      close(52)
      lexist = .false.
      return
   end if
!
   open(unit=53, file='sparse_info.data', form='binary', status='old')
   rewind(53)
   read(53)lni, ltmp1, ltmp2, ltmp3, ltmp4
   if (lni .ne. nunk .or. ltmp1 .ne. offcor .or. ltmp2 .ne. lequileads .or. &
       ltmp3 .ne. lsimple .or. ltmp4 .ne. lscreen) then
      write(6,*)'outputnnz: <sparse_info.data> incompatible with present job '
      write(6,*)'outputnnz: abort reading '
      write(6,*)lni, ltmp1, ltmp2, ltmp3, ltmp4
      write(6,*)nunk, offcor, lequileads, lsimple, lscreen
      call flush(6)
      close(53)
      lexist = .false.
      return
   end if
end if

 
iread_spa = 2
!
if (lexist) then
   if (iread_spa .eq. 2) then  ! read sparsity info for present job
      read(53)lunk_spa
      close(unit=52)
      write(6,*)'outputnnz: <sparse_index.data> and <sparse_info.data> read successfully'
      call flush(6)
      goto 101
   else 
      write(6,*)'outputnnz: unknown error!'
      return
   end if
end if

101 continue

! output sparisty information

 open(unit = 80, file='RDO_and_ADO.data')
 rewind(80)
 lunk_last = lunk_spa
 
 do itier=1,2
!
   do lni=nfirst(itier),nlast(itier)
      call flush(80)
      lnl = index_coef_ref(lni)
!
      lnm = ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1)
!      write(80,*)'lni,lnl,lnm',lni,lnl,lnm
      ! mindex represents {j1,j2,...,jn} which is the subscript of rho^{n}_{j1,j2,...,jn}
      mindex(1:itier-1) = indextable(lnm+1:lnm+itier-1)
!
      write(80,*)'ntier-ADO balls  ', '   j1,    j2,    j3,    ....,    jn   '
      write(80,1001) (mindex(mi), mi = 1, itier-1)
      if (psdfff) then
          write(80,*) 'isgn  ', 'iorbs  ', 'ispin  ', 'ialf  ', 'icor  ', 'itype  ', 'mp'
         else
          write(80,*) 'isgn  ', 'iorbs  ', 'ispin  ', 'ialf  ', 'icor  ', 'itype'
         end if
      do ni=1,itier-1
         ! nk reads each ji in {j1,j2,...,jn}
         nk = indextable( ifirst(itier) - 1 + (lni - nfirst(itier)) * (itier - 1) + ni )
!         if (psdfff .and. mmfff(nk) .gt. 1) then
!             call findlowertier(itier, mindex, ni, lnk, itype, fact)
!             if (itier .eq. 2) then
!               lnj = 1
!             else
!               lnj = (lnk + 1 - ifirst(itier-1)) / (itier - 2) + nfirst(itier-1)
!             end if
!         else
!             lnj    = indexcoef(lnl)
             itype   = itypecoef(lnl)
!         end if
         !!!! if you want to know the details of each ball, print these following variables
         isgn    = mpm(nk)
         iorbs   = morbs(nk)
         ispin   = mspin(nk)
         ialf    = ilead(nk)
         icor    = mcor(nk)
!         lnj     = lnj
         itype    = itype
         if (psdfff) then
          write(80,1002) isgn, iorbs, ispin, ialf, icor, itype, mmfff(nk)
         else
          write(80,1003) isgn, iorbs, ispin, ialf, icor, itype
         end if
         !!! if you want to know the details of each ball, print these above variables
      end do
      lstart   = ind_spa(lni)
      nnz     = nnz_spa(lni)
      write(80,*)'num:', lni
      write(80,*)'      irow     ', '      icol     ', '       Re[rho]       ', '          Im[rho]  '
      do nj = 1, nnz
           write(80,1004) irow_spa(lstart - 1 + nj), icol_spa(lstart - 1 + nj), dble(rho_spa(lstart -1 + nj)), dimag(rho_spa(lstart -1 + nj))
      end do
      write(80,*)
  end do
 end do
 call flush(80)
 close(80)

 write(6,*)'outputnnz: done'
 write(6,*)
 
 1001 format(18x,32(I5,2x))
 1002 format(2(I4, 2x), 2(1x, I4, 2x), 3(I4, 2x))
 1003 format(2(I4, 2x), 2(1x, I4, 2x), 2(I4, 2x))
 1004 format(2(4x, I6, 4x),2(4x,F16.12,4x))
 end subroutine outputnnz
