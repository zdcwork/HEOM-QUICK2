subroutine getengyshift(tt, ialf, ispin, eshift)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)  :: ialf, ispin     ! 1 for left, 2 for right
integer, parameter   :: maxpul = 10
integer              :: ntmp1, ntmp2, numsec1, numsec2, istat, npt, ni, nj
logical              :: lfound, lfound2
character*120        :: cline
real*8,  intent(in)  :: tt
real*8,  intent(out) :: eshift
real*8               :: tip(maxpul), tfp(maxpul), ampp(maxpul), omega_second(ialf,ispin), tchar_second(ialf,ispin)
!
if (tt >= t_off) then
   eshift = 0.d0
   return 
end if
!
if (fieldtype .eq. 0) then
   eshift = amps(ialf,ispin) * (1.d0 - dexp(-tt / tchar(ialf,ispin)))
!
else if (fieldtype .eq. 1) then
   eshift = amps(ialf,ispin) * dsin(2.0 * pi * tt / tchar(ialf,ispin))
   if (leshift_dc) then
       eshift = eshift + amps_dc(ialf,ispin)
   end if
!
else if (fieldtype .eq. 2) then
   eshift = amps(ialf,ispin) * (1.d0 - dexp(-tt / tchar(ialf,ispin)))
!
! pulsed_energy_shift  ialf  npulse    !note: npulse must not exceed maxpul
! (ti(ipulse), tf(ipulse), amp(ipulse), ipulse=1,npulse)
!
   lfound = .false.
   rewind(5)
   find_pulse: do
      read(5, '(A120)', iostat=istat) cline
      if (istat .ne. 0) exit find_pulse
      npt = index(cline, 'pulsed_energy_shift')
      if (npt .gt. 0) then
          npt = npt + 19
!write(6,*)'npt=',npt
!call flush(6)
          read(cline(npt:120),*,iostat=istat) ntmp1, ntmp2
!write(6,*)'ntmp1, ntmp2', ntmp1, ntmp2
!write(6,*)'istat, ialf', istat, ialf
!call flush(6)
          if (istat .eq. 0 .and. ntmp1 .eq. ialf) then
              lfound = .true. 
              exit find_pulse
          end if
      end if
   end do find_pulse
   if (lfound) then
       if (ntmp2 .le. 0 .or. ntmp2 .gt. maxpul) then
           write(6,*)'getenergyshift: error pulse number found ', ntmp2
           write(6,*)'getenergyshift: maximal pulses allowed   ', maxpul
           stop
       end if
!write(6,*)'ntmp2 = ', ntmp2
!call flush(6)
       read(5,*,iostat=istat) (tip(ntmp1), tfp(ntmp1), ampp(ntmp1), ntmp1=1,ntmp2)
       if (istat .ne. 0) then
           write(6,*)'getengyshift: error reading primary pulse for lead ', ialf
           stop
       end if
!write(6,*)'after read '
!call flush(6)
       ampp(1:ntmp2) = ampp(1:ntmp2) / hbar
       find_interval: do ntmp1=1,ntmp2
           if (tt .ge. tip(ntmp1) .and. tt .le. tfp(ntmp1)) then
               eshift = eshift + ampp(ntmp1)
               exit find_interval
           end if
       end do find_interval
   end if
!
!fieldtype = 3
else if (fieldtype .eq. 3) then
   eshift = 0.d0
!
! pulsed_energy_shift  ialf  npulse    !note: npulse must not exceed maxpul
! (ti(ipulse), tf(ipulse), amp(ipulse), ipulse=1,npulse)
!
   lfound = .false.
   rewind(5)
   find_pulse_3: do
      read(5, '(A120)', iostat=istat) cline
      if (istat .ne. 0) exit find_pulse_3
      npt = index(cline, 'pulsed_energy_shift')
      if (npt .gt. 0) then
          npt = npt + 19
!write(6,*)'npt=',npt
!call flush(6)
          read(cline(npt:120),*,iostat=istat) ntmp1, ntmp2
!write(6,*)'ntmp1, ntmp2', ntmp1, ntmp2
!write(6,*)'istat, ialf', istat, ialf
!call flush(6)
          if (istat .eq. 0 .and. ntmp1 .eq. ialf) then
              lfound = .true. 
              exit find_pulse_3
          end if
      end if
   end do find_pulse_3
   if (lfound) then
       if (ntmp2 .le. 0 .or. ntmp2 .gt. maxpul) then
           write(6,*)'getenergyshift: error pulse number found ', ntmp2
           write(6,*)'getenergyshift: maximal pulses allowed   ', maxpul
           stop
       end if
!write(6,*)'ntmp2 = ', ntmp2
!call flush(6)
       read(5,*,iostat=istat) (tip(ntmp1), tfp(ntmp1), ampp(ntmp1), ntmp1=1,ntmp2)
       if (istat .ne. 0) then
           write(6,*)'getengyshift: error reading fieldtype=3 pulse for lead ', ialf
           stop
       end if
!write(6,*)'after read '
!call flush(6)
       ampp(1:ntmp2) = ampp(1:ntmp2) / hbar
       find_interval_3: do ntmp1=1,ntmp2
           if (tt .ge. tip(ntmp1) .and. tt .le. tfp(ntmp1)) then
               eshift = eshift + ampp(ntmp1) * dsin(2.0 * pi * tt /tchar(ialf,ispin)+ phase(ialf,ispin))
               exit find_interval_3
           end if
       end do find_interval_3
   end if
!
!   if (tt .ge. (tcetr(ialf,ispin)-0.5 * twidth(ialf,ispin)) .and. (tt .le. (tcetr(ialf,ispin)+0.5 * twidth(ialf,ispin)))) then
!      eshift = amps(ialf,ispin) * dsin(2.0 * pi * tt / tchar(ialf,ispin))! * dexp(-(tt-tcetr(ialf,ispin)) * (tt-tcetr(ialf,ispin)) / twidth(ialf,ispin) / twidth(ialf,ispin))
!   else 
!      eshift = 0.d0
!   end if 

   lfound2 = .false.
   find_second_pulse: do
      read(5, '(A120)', iostat=istat) cline
      if (istat .ne. 0) exit find_second_pulse
      npt = index(cline, 'second_pulsed_shift')
      if (npt .gt. 0) then
          npt = npt + 19
!write(6,*)'npt=',npt
!call flush(6)
          read(cline(npt:120),*,iostat=istat) numsec1, numsec2
!write(6,*)'numsec1, numsec2', numsec1, numsec2
!write(6,*)'istat, ialf', istat, ialf
!call flush(6)
          if (istat .eq. 0 .and. numsec1 .eq. ialf) then
              lfound2 = .true. 
              exit find_second_pulse
          end if
      end if
   end do find_second_pulse
   if (lfound2) then
       if (numsec2 .le. 0 .or. numsec2 .gt. maxpul) then
           write(6,*)'getenergyshift: error second pulse number found ', numsec2
           write(6,*)'getenergyshift: maximal pulses allowed   ', maxpul
           stop
       end if
!write(6,*)'ntmp2 = ', ntmp2
!call flush(6)
       read(5,*,iostat=istat) ((omega_second(ni,nj), nj=1,ispin), ni=1,ialf)
       if (istat .ne. 0) then
          write(6,*)'getenergyshift: error occurred while reading the second pulse oemga, ialf, ispin= ', ialf, ispin
       end if
       omega_second(1:ialf,1:ispin) = omega_second(1:ialf,1:ispin) / hbar
       tchar_second(1:ialf,1:ispin) = 2.d0 * pi / omega_second(1:ialf,1:ispin)
       read(5,*,iostat=istat) (tip(numsec1), tfp(numsec1), ampp(numsec1), numsec1=1,numsec2)
       if (istat .ne. 0) then
           write(6,*)'getengyshift: error reading second pulse for lead ', ialf
           stop
       end if
!write(6,*)'after read '
!call flush(6)
       ampp(1:numsec2) = ampp(1:numsec2) / hbar
       find_interval_second: do numsec1=1,numsec2
           if (tt .ge. tip(numsec1) .and. tt .le. tfp(numsec1)) then
               eshift = eshift + ampp(numsec1) * dsin(2.0 * pi * tt /tchar_second(ialf,ispin)+ phase(ialf,ispin))
               exit find_interval_second
           end if
       end do find_interval_second
   end if


   if (leshift_dc) then
       eshift = eshift + amps_dc(ialf,ispin)
   end if
else if (fieldtype .eq. -1) then
   eshift = 0.d0
else
   write(6,*)'getengyshift: error! unknown fieldtype =', fieldtype
   stop
end if
!
end subroutine getengyshift
