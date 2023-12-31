subroutine getengyshift(tt, ialf, ispin, eshift)
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)  :: ialf, ispin     ! 1 for left, 2 for right
integer, parameter   :: maxpul = 10
integer              :: ntmp1, ntmp2, istat, npt
logical              :: lfound
character*120        :: cline
real*8,  intent(in)  :: tt
real*8,  intent(out) :: eshift
real*8               :: tip(maxpul), tfp(maxpul), ampp(maxpul)
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
else if (fieldtype .eq. -1) then
   eshift = 0.d0
else
   write(6,*)'getengyshift: error! unknown fieldtype =', fieldtype
   stop
end if
!
end subroutine getengyshift
