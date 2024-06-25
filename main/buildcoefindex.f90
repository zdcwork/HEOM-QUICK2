subroutine buildcoefindex
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer   :: index_nk(MAXTIER), iopera(MAXTIER), index_nj(MAXTIER), index_np(MAXTIER)
logical   :: iexist1, iexist2, iexist3, iexist, ltmp1, ltmp2, ltmp3
logical   :: lslow1, lslow2, lfast1, lfast2
integer   :: itier, istat, mdim
integer*1 :: ifactfront, ifactrear, fact, iout
integer   :: ni, nj, nk, nl, nm, np, nq, nd
integer   :: mi, mj
integer   :: isgn, ialf, iorbs1, iorbs2, ispin, imats, iopr
integer   :: nfast, idraw, ndegen
integer   :: ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
integer   :: nvec_tmp1(MAXTIER), nvec_tmp2(MAXTIER)
integer*8 :: lni, lnj, lnk, lnl, lnm
real*8    :: cpu1, cpu2, cpua, cpub
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4
!
call cpu_time(cpua)
write(6,*)
write(6,*)' entering buildcoefindex     '
write(6,*)' this could take a while ... '
call flush(6)
!
allocate(index_coef_ref(nunk), STAT=istat)
index_coef_ref = 0
lnl = 1
!
if (lsimple) then
   allocate(ioprref(nunk), STAT=istat)
   ioprref = 0
   lnm = 1
end if
!
iexist = .false.
!
! Be careful with reading, make sure the existing coefindex.data is 
! compatible with the ongoing job.
!
inquire(file='coefindex.data', exist=iexist1, err=999)
999 continue
inquire(file='auxindex.data', exist=iexist3, err=997)
997 continue
if (iexist1 .and. iexist3) then
  iexist = .true.
end if
!
if (lsimple .and. iexist) then
   inquire(file='oprindex.data', exist=iexist2, err=996)
   996 continue
   if (.not. iexist2) iexist = .false.
end if
!
if (iexist) then
   write(6,*)
   write(6,*)' coefindex.data found, start to read '
   open(unit=18, file='coefindex.data', form='binary', status='old')
   rewind(18)
   read(18)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
   if (ntmp1 .ne. ntier .or. ntmp2 .ne. ncor .or. ntmp3 .ne. norbs .or.     &
       ntmp4 .ne. nspin .or. ntmp5 .ne. nalf .or. ntmp6 .ne. numfff .or.    &
       ntmp7 .ne. ntier0 .or. ntmp8 .ne. ndrawer_slow .or. ntmp9 .ne. ncor_slow) then
       write(6,*) 
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ntmp1, ntmp2, ntmp3, ntmp4, ntmp5, ntmp6, ntmp7, ntmp8, ntmp9
       write(6,*)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   read(18)ltmp1, ltmp2, ltmp3, ntmp1, ntmp2, dtmp1
   if (ltmp1 .ne. lad) then
       write(6,*)
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)' ltmp1, lad ', ltmp1, lad
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   if (lad) then
       if (ltmp2 .ne. lad_fast .or. ltmp3 .ne. lset_fast .or.       &
           ntmp1 .ne. idegen_fast .or. ntmp2 .ne. ndegen_fast .or.  &
           dabs(dtmp1 - dratio_fast) .gt. dpico) then
           write(6,*)
           write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
           write(6,*)ltmp1, ltmp2, ltmp3, ntmp1, ntmp2, dtmp1
           write(6,*)lad, lad_fast, lset_fast, idegen_fast, ndegen_fast, dratio_fast
           call flush(6)
           close(18)
           iexist = .false.
           goto 15
       end if
   end if
   read(18)ltmp1, ntmp2, ntmp3
   if (ltmp1 .ne. lfilter .or.                                              &
       lfilter .and. (ntmp2 .ne. nfilter_count .or. ntmp3 .ne. nfilter_long)) then 
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ltmp1, ntmp2, ntmp3
       write(6,*)lfilter, nfilter_count, nfilter_long
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   read(18)ltmp1
   if (ltmp1 .ne. ltrun_der) then
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ltmp1
       write(6,*)ltrun_der
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   read(18)ltmp1, ltmp2
   if (ltmp1 .ne. lsimple .or. ltmp2 .ne. lwalf) then
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ltmp1, lsimple
       write(6,*)ltmp2, lwalf
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   read(18)ltmp1
   if (ltmp1 .ne. lscreen) then
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ltmp1, lscreen
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   read(18)ltmp1, dtmp1
   if (ltmp1 .ne. lomit) then
       write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
       write(6,*)ltmp1, lomit
       call flush(6)
       close(18)
       iexist = .false.
       goto 15
   end if
   if (lomit) then
       if ( dabs(dtmp1 - dfact_omit) .gt. dpico) then
           write(6,*)'buildcoefindex: coefindex.data found incompatible! abort reading '
           write(6,*)dtmp1, dfact_omit
           call flush(6)
           close(18)
           iexist = .false.
           goto 15
       end if
   end if
   read(18)lall
!
   open(unit=37, file='auxindex.data', form='binary', status='old')
   rewind(37)
   read(37)lni
   if (lni .ne. nunk) then
      write(6,*)'buildcoefindex: error length from <auxindex.data> '
      write(6,*)'length(read), length(expected) ', lni, nunk 
      stop
   end if
!
   if (lsimple) then
      open(unit=67, file='oprindex.data', form='binary', status='old')
      rewind(67)
      read(67)ntmp1, ntmp2
      if (ntmp1 .ne. ntier .or. ntmp2 .ne. ntier0) then
         write(6,*)'buildcoefindex: <oprindex.data> found incompatible! '
         write(6,*)ntmp1, ntier
         stop
      end if
      read(67)lopr
   end if
!
   goto 101
!
else
  write(6,*)
  write(6,*)' coefindex.data not found '
  call flush(6)
end if
!
15 continue
!
mdim = MAXTIER
open(unit=16, file='coefindex.tmp', form='binary', status='unknown')
rewind(16)
lall = 0
!
if (lsimple) then
   open(unit=66, file='oprindex.tmp', form='binary', status='unknown') 
   rewind(66)
   lopr = 0
end if
!
do itier=2,ntier0
   call cpu_time(cpu1)
   dtmp1 = dcbmax**(itier-1) * dfact_omit
   do lni=nfirst(itier), nlast(itier)
      index_coef_ref(lni) = lnl
      if (lsimple) ioprref(lni) = lnm
!
!  Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
      dtmp2 = 1.d0
      do ni=1,itier-1 
         index_nk(ni) = indextable( ifirst(itier) - 1 +  (lni - nfirst(itier)) * (itier - 1) + ni )
         nd = index_nk(ni)
         call look4operator(mpm(nd), morbs(nd), mspin(nd), ilead(nd), iopera(ni))
         dtmp2 = dtmp2 * cdabs(cb(morbs(nd),mspin(nd),mcor(nd),ilead(nd),mpm(nd)))
      end do
!
      if (lad .and. itier .eq. ntier) then
          lslow1 = .true.
          lfast1 = .true.
          do ni=1,itier-1   ! all slow modes? 
             if (mslow(index_nk(ni)) .ne. 1) then
                 lslow1 = .false.
                 exit
             end if
          end do
          do ni=1,itier-1
             if (mslow(index_nk(ni)) .eq. 1) then
                 lfast1 = .false.
                 exit
             end if
          end do
      end if
!
!  nearest lower tier contribution (N - 1)
!
      do ni=1,itier-1
!  same tier if psdfff scheme is invoked
         if (psdfff .and. mmfff(index_nk(ni)) .gt. 1) then
             call findsametier(itier, index_nk(1), ni, index_nk(ni)-1, lnj, iout, fact)
             if (iout .ne. 0) then
                 write(6,*)'buildcoefindex: unexpected iout ', iout, itier, lni
                 stop
             end if
             lnk = (lnj + 1 - ifirst(itier)) / (itier - 1) + nfirst(itier)
             ifactfront = fact
             ifactrear  = fact
             lnl  = lnl + 1
             lall = lall + 1
             write(16)lnk, iout, ifactfront, ifactrear
             cycle
         end if
         if (lomit .and. itier .gt. 3) then
             nd = index_nk(ni)
             dtmp4 = dtmp2 / cdabs(cb(morbs(nd),mspin(nd),mcor(nd),ilead(nd),mpm(nd)))
             if (dtmp4 .lt. dtmp1/dcbmax) cycle 
         end if
         call findlowertier(itier, index_nk(1), ni, lnj, iout, fact)
         if (itier .eq. 2) then
             lnk = 1
         else
             lnk = (lnj + 1 - ifirst(itier-1)) / (itier - 2) + nfirst(itier-1)
         end if
!
         if (mod(ni, 2) .eq. 0) then
             ifactfront = -1
         else
             ifactfront = 1
         end if
         if (mod(itier-1-ni, 2) .eq. 0) then 
             ifactrear = 1 
         else
             ifactrear = -1
         end if
         if (iout .eq. 1) then
             ifactfront = ifactfront * fact
             ifactrear  = ifactrear  * fact
         end if
         lnl  = lnl + 1
         lall = lall + 1
         write(16)lnk, iout, ifactfront, ifactrear
!        write(6,*)'lower', lni, ' --> ', lnk, ' itype= ', iout
      end do
!
!  nearest higher tier contribution (N + 1)  
!
      if (itier .lt. ntier) then
          do isgn=1,nsgn
             do ispin=1,nspin
                orbs0: do iorbs1=1,norbs
!
                   if (lsimple .and. .not. lwalf) then
                      call look4sysopr(isgn, iorbs1, ispin, iopr)
                      ltmp1 = .false. 
                      iter1: do nj=1,itier-1
                         if (iopr .eq. msopr(index_nk(nj))) then
                            ltmp1 = .true.
                            exit iter1
                         end if
                      end do iter1
                      lopr = lopr + 1
                      lnm = lnm + 1
                      if (ltmp1) then
                         write(66)lopr, -1
                         cycle orbs0
                      else 
                         write(66)lopr, iopr
                      end if
                   end if
!
                   alf0: do ialf=1,nalf
                      call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                      if (lscreen .and. jomit(iopr) .eq. 1) then
                          cycle alf0
                      end if
                      if (lsimple .and. lwalf) then
                          ltmp1 = .false. 
                          iter0: do nj=1,itier-1
                             if (iopr .eq. mopr(index_nk(nj))) then
                                 ltmp1 = .true.
                                 exit iter0
                             end if
                          end do iter0
                          lopr = lopr + 1
                          lnm = lnm + 1
                          if (ltmp1) then
                              write(66)lopr, -1
                              cycle alf0
                          else 
                              write(66)lopr, iopr
                          end if
                      end if
!
                      do imats=1,ncor
                         if (lfilter .and. itier .eq. ntier-1) then
                             ntmp1 = filtertable(lni)
                             if (kpmats(imats) .le. nfilter_long) ntmp1 = ntmp1 + 1
                             if (ntmp1 .lt. nfilter_count) cycle
                         end if
                         call look4drawer(isgn, ialf, iorbs1, ispin, imats, nq)
                         if (lomit) then
                             dtmp3 = dtmp2 * cdabs(cb(morbs(nq),mspin(nq),mcor(nq),ilead(nq),mpm(nq)))
                             if (dtmp3 .lt. dtmp1*dcbmax) cycle
                         end if
                         call findhighertier(itier, index_nk(1), nq, lnj, iout, ifactfront, ifactrear)
                         if (iout .lt. 0) then
                             write(6,*)
                             write(6,*)'buildcoefindex: error! failed to find higher tier ADO. F1 '
                             stop
                         end if
                         lnk  = (lnj + 1 - ifirst(itier+1)) / itier + nfirst(itier+1) 
                         lnl  = lnl + 1
                         lall = lall + 1
                         write(16)lnk, iout, ifactfront, ifactrear
                      end do
                   end do alf0
               end do orbs0
             end do
          end do
      end if
!
      if (lad) then
          if (itier .eq. ntier) then
              do isgn=1,nsgn
                 do ispin=1,nspin
                    do iorbs1=1,norbs
                       do ialf=1,nalf
                          call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                          if (lscreen .and. jomit(iopr) .eq. 1) cycle
                          !
                          do imats=1,ncor
                             call look4drawer(isgn, ialf, iorbs1, ispin, imats, nq)
                             lslow2 = .false.
                             lfast2 = .false.
                             if (mslow(nq) .eq. 1) then
                                 lslow2 = .true.
                             else
                                 lfast2 = .true.
                             end if
                             !
                             if (lad_fast .and. lfast1 .and. lfast2 .or.            &
                                 (.not. lad_fast) .and. lslow1 .and. lslow2 ) then
                                 if (itier + 1 .gt. ntier_ad) cycle
                                 call findhighertier(itier, index_nk, nq, lnj, iout, ifactfront, ifactrear)
                                 if (iout .lt. 0) then
                                     write(6,*)
                                     write(6,*)'buildcoefindex: error! failed to find higher tier ADO. F2 '
                                     write(6,*)itier, index_nk(1:itier-1), nq
                                     stop
                                 end if
                                 lnk  = (lnj + 1 - ifirst(itier+1)) / itier + nfirst(itier+1) 
                                 lnl  = lnl + 1
                                 lall = lall + 1
                                 write(16)lnk, iout, ifactfront, ifactrear
                                 cycle
                             end if
                             !
                             call add2list(index_nk, mdim, itier-1, index_nj, mdim, nq, mj)
                             ! store fl = (-1)**(mj-1) in memory, 
                             !       fr = (-1)**(itier-mj), where itier is length of index_nj,
                             !             and (itier-1) is the length of index_nk
                             !            is to be determined by fr = fl * (-1)**(itier-1)
                             ! do NOT store fr to save the use of memory
                             ifactfront = (-1)**(mj - 1)  ! fl
                             !
                             ! find the fastest mode
! note added on Jan 30, 2020
! there may be more than one fastest mode (same imats, different iopr)
! note added on Feb 04, 2020
! there may be many modes that have rates close enough to the fastest mode, and
! these fast modes constitute a set of degenerate modes, when applying adiabatic
! approximation to these modes, the symmetry of approximated form may be important 
                             mi    = 0
                             ntmp1 = 0 
                             do mj=1,itier
                                if (ntmp1 .lt. ms2f(index_nj(mj))) then
                                    ntmp1 = ms2f(index_nj(mj))
                                    mi = mj
                                end if
                             end do
                             nfast = mi
!!! debug
!ntmp1 = index_nj(mi)
!if (mslow(ntmp1) .eq. 1) then
!    write(6,*)
!    write(6,*)'buildcoefindex: error! fastest mode is slow '
!    write(6,*)ntmp1,  mslow(ntmp1)
!    stop
!end if
!!! debug
                             call findlowertier(itier+1, index_nj, mi, lnj, iout, fact)
                             ! store gr = (-1)**(itier-mj) * fact in memory, where itier is length of index_nj,
                             !       gl = (-1)**(mj - 1) * fact is to be determined by 
                             !       gl = gr * (-1)**(itier-1)
                             ! do NOT store gl to save the use of memory
                             ifactrear  = (-1)**(itier - mi) * fact   ! gr
                             lnk = (lnj + 1 - ifirst(itier)) / (itier - 1) + nfirst(itier) 
                             !
                             ! if the upper-tier ADO contains identical drawers, its value is zero,
                             ! record it by setting ifactfront to zero (instead of +1/-1)
                             do mi=1,itier-1
                                if (index_nj(mi) .eq. index_nj(mi+1)) then
                                    ifactfront = 0
                                    exit
                                end if
                             end do
                             lnl  = lnl + 1
                             lall = lall + 1
                             write(16)lnk, iout, ifactfront, ifactrear
!
                             ! cycle if the upper-tier ADO contains identical modes
                             if (ifactfront .eq. 0) cycle
!
                             ! treat degenerate cases 
                             if (idegen_fast .eq. 0) cycle
                             !
                             nvec_tmp1(1) = index_nj(nfast)
                             nvec_tmp2(1) = nfast
                             ntmp1 = 1
                             do mi=1,itier
                                if (mi .eq. nfast) cycle
                                idraw = index_nj(mi)
                                if (dgama_drawer(idraw) .lt. dgama_drawer(index_nj(nfast))*dratio_fast) cycle
                                ntmp1 = ntmp1 + 1
                                nvec_tmp1(ntmp1) = idraw
                                nvec_tmp2(ntmp1) = mi
                             end do
                             ndegen = ntmp1
                             if (ndegen .eq. 1) cycle
!
                             if (idegen_fast .eq. 1) then
                                 ntmp1 = index_nj(nfast)
                                 do mj=2,ndegen,1
                                    idraw = nvec_tmp1(mj)
                                    mi    = nvec_tmp2(mj)
                                    call findlowertier(itier+1, index_nj, mi, lnj, iout, fact)
                                    ifactrear  = (-1)**(itier - mi) * fact   ! gr
                                    lnk = (lnj + 1 - ifirst(itier)) / (itier - 1) + nfirst(itier) 
                                    lnl  = lnl + 1
                                    lall = lall + 1
                                    write(16)lnk, iout, ifactfront, ifactrear
                                 end do
                             else if (idegen_fast .eq. 2) then
                                 if (ndegen .gt. ndegen_fast) cycle
                                 call findsubvec(index_nj, itier, nvec_tmp1, ndegen, index_np, ntmp2, nvec_tmp2)
                                 call findanytier(ntmp2+1, index_np, lnj, iout, fact)
                                 if (lnj .eq. 1) then
                                     lnk = 1
                                 else 
                                     lnk = (lnj + 1 - ifirst(ntmp2+1)) / ntmp2 + nfirst(ntmp2+1) 
                                 end if
                                 lnl  = lnl + 1
                                 lall = lall + 1
                                 write(16)lnk, iout, fact, fact
                             end if
                          end do
                       end do
                    end do
                 end do
              end do
!
          else if (itier .gt. ntier .and. itier .lt. ntier_ad) then
              if (lad_fast) then
                  ntmp1 = ncor_fast
              else
                  ntmp1 = ncor_slow
              end if
              !
              do isgn=1,nsgn
                 do ispin=1,nspin
                    do iorbs1=1,norbs
                       do ialf=1,nalf
                          call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                          if (lscreen .and. jomit(iopr) .eq. 1) cycle
                          !
                          do imats=1,ntmp1
                             if (lad_fast) then
                                 call look4drawer(isgn, ialf, iorbs1, ispin, icor_fast(imats), nq)
                             else
                                 call look4drawer(isgn, ialf, iorbs1, ispin, icor_slow(imats), nq)
                             end if
                             !
                             call findhighertier(itier, index_nk(1), nq, lnj, iout, ifactfront, ifactrear)
                             if (iout .lt. 0) then
                                 write(6,*)
                                 write(6,*)'buildcoefindex: error! failed to find higher tier ADO. F3 '
                                 write(6,*)itier, index_nk(1:itier-1), nq
                                 stop
                             end if
                             lnk  = (lnj + 1 - ifirst(itier+1)) / itier + nfirst(itier+1) 
                             lnl  = lnl + 1
                             lall = lall + 1
                             write(16)lnk, iout, ifactfront, ifactrear
                          end do
                       end do
                    end do
                 end do
              end do
!
          end if
      end if
!
!
!  Truncation for the derivatives of terminal-tier ADOs
!
      if (ltrun_der .and. itier .eq. ntier) then
          do isgn=1,nsgn
             do ispin=1,nspin
                do iorbs1=1,norbs
                   do ialf=1,nalf
                      call look4operator(isgn, iorbs1, ispin, ialf, iopr)
                      if (lscreen .and. jomit(iopr) .eq. 1) cycle
                      !
                      do imats=1,ncor
                         call look4drawer(isgn, ialf, iorbs1, ispin, imats, nq)
                         call add2list(index_nk, mdim, ntier-1, index_nj, mdim, nq, mj)
                         ifactfront = (-1)**(mj - 1)
                         ifactrear  = (-1)**(ntier - mj)
                         lnl  = lnl + 1
                         lall = lall + 1
                         lnk  = nq
                         iout = 0
                         write(16)lnk, iout, ifactfront, ifactrear 
                         do mi=1,ntier
                            call findlowertier(ntier+1, index_nj, mi, lnj, iout, fact)
                            lnk = (lnj + 1 - ifirst(ntier)) / (ntier - 1) + nfirst(ntier)
                            ifactfront = (-1)**(mi - 1)     
                            ifactrear  = (-1)**(ntier - mi)
                            if (iout .eq. 1) then
                               ifactfront = ifactfront * fact
                               ifactrear  = ifactrear  * fact
                            end if
                            lnl  = lnl + 1
                            lall = lall + 1
                            write(16)lnk, iout, ifactfront, ifactrear
                         end do
                      end do
                   end do
                end do
             end do
          end do
      end if
!
   end do
   call cpu_time(cpu2)
   write(6,100)itier, nlast(itier)-nfirst(itier)+1, cpu2 - cpu1
   call flush(6)
end do 
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
close(16)
close(26)
if (lsimple) close(66)
!
101 continue
!
! the integer(kind=8) uses 8 byte of memory
! the integer(kind=4) uses 4 byte of memory (confirmed by test)
!
!dtmp1 = dble(lall * (8 + 4 * 3) + nunk * 8) / dble(1024**2)    ! nunk is for index_coef_ref
dtmp1 = dble(lall * (8 + 1 * 3) + nunk * 8) / dble(1024**2)    ! nunk is for index_coef_ref
write(6,*)
write(6,*)' total searches made      : ', lall
write(6,*)' estimated general memory : ', dtmp1, ' MB '
call flush(6)
!
allocate(indexcoef(lall), itypecoef(lall), ifactfcoef(lall), ifactrcoef(lall), STAT=istat)
if (lsimple) then
   allocate(ioprindex(lopr), STAT=istat)
end if
!
write(6,*)
if (istat .eq. 0) then
 write(6,*)' memory allocation for indexcoef successful! '
 memory = memory + dtmp1
else 
 write(6,*)' memory allocation for indexcoef failed!     '
 stop
end if
call flush(6)
!
! write to file (and memory if possible)
!
if (iexist) then
  do lni=1,lall
    read(18)lnj, iout, ifactfront, ifactrear
    indexcoef(lni)  = lnj
    itypecoef(lni)  = iout
    ifactfcoef(lni) = ifactfront
    ifactrcoef(lni) = ifactrear
  end do
  close(18)
!
  do lni=1,nunk
     read(37)lnl
     index_coef_ref(lni) = lnl
  end do
  if (lsimple) then
     do lni=1,nunk
        read(37)lnm
        ioprref(lni) = lnm
     end do
  end if
  close(37)
!
  if (lsimple) then
     do lni=1,lopr
        read(67)ioprindex(lni)
     end do
     close(67)
  end if
!
  goto 200
end if
!
! read into memory
!
open(unit=16, file='coefindex.tmp', form='binary', status='unknown')
rewind(16)
do lni=1,lall
  read(16) lnj, iout, ifactfront, ifactrear
  indexcoef(lni)  = lnj
  itypecoef(lni)  = iout
  ifactfcoef(lni) = ifactfront
  ifactrcoef(lni) = ifactrear
end do
close(unit=16, status="delete")
!
if (lsimple) then
   open(unit=66, file='oprindex.tmp', form='binary', status='unknown')
   rewind(66)
   do lni=1,lopr
      read(66) lnj, nj
      ioprindex(lni) = nj
   end do
   close(unit=66, status="delete")
end if
!
! write to file
!
!open(unit=17, file='coefindex.data', form='binary', status='unknown', asynchronous='yes')
open(unit=17, file='coefindex.data', form='binary', status='unknown')
rewind(17)
write(17)ntier, ncor, norbs, nspin, nalf, numfff, ntier0, ndrawer_slow, ncor_slow
write(17)lad, lad_fast, lset_fast, idegen_fast, ndegen_fast, dratio_fast
write(17)lfilter, nfilter_count, nfilter_long
write(17)ltrun_der
write(17)lsimple, lwalf
write(17)lscreen
write(17)lomit, dfact_omit
write(17)lall
do lni=1,lall
   write(17)indexcoef(lni), itypecoef(lni), ifactfcoef(lni), ifactrcoef(lni)
end do
close(17)
!
open(unit=37, file='auxindex.data', form='binary', status='unknown')
rewind(37)
write(37)nunk
do lni=1,nunk
   write(37)index_coef_ref(lni)
end do
if (lsimple) then
   do lni=1,nunk
      write(37)ioprref(lni)
   end do
end if
close(37)
!
if (lsimple) then
   open(unit=67, file='oprindex.data', form='binary', status='unknown')
   rewind(67)
   write(67)ntier, ntier0
   write(67)lopr
   do lni=1,lopr
      write(67)ioprindex(lni)
   end do
   close(67)
end if
!
200 continue
!
call cpu_time(cpub)
write(6,*)
write(6,*)' leaving buildcoefindex ', cpub - cpua
call flush(6)
!
return
end subroutine buildcoefindex
