!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! debug utilities developed by sun jian
! for benchmarking the time consuming of each part of the code
!
! this module contains:
!    subroutines:
!      sj_dbg_print_time
!         print out time infomation together with the msg 
!      sj_dbg_timer_start
!         start the timer (ntimer)
!      sj_dbg_timer_stop
!         stop the timer (ntimer) and report the elapsed time and msg (optional)
!
!    aliases to the subroutines listed above
!      sj_print_time
!         alias for sj_dbg_print_time
!      sj_cpu_time 
!         alias for sj_dbg_print_time, not recommended
!      sj_start_timer
!         alias for sj_dbg_start_timer
!      sj_stop_timer
!         alias for sj_dbg_stop_timer
!
!    parameters:
!      msg:    the message you want to print out after the time 
!      ntimer: the identifier of your timer, ranged from 1 to 255
!              required parameter for sj_dbg_timer_start and sj_dbg_timer_stop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sj_timer
   implicit none

   private  clock_start, cpu_start, clock_rate
   public   sj_cpu_time, sj_dbg_print_time, sj_dbg_timer_start, sj_dbg_timer_stop

   integer*8 :: clock_start(255)=-1
   real*8    :: cpu_start(255)=-1.0
   integer*8 :: clock_rate =0

   interface sj_print_time  ! a short alias 
      module procedure sj_dbg_print_time
   end interface
   interface sj_cpu_time    ! an alias for compatible consideration, not recommended
      module procedure sj_dbg_print_time
   end interface
   interface sj_timer_start ! a short alias 
      module procedure sj_dbg_timer_start
   end interface
   interface sj_timer_print ! a short alias 
      module procedure sj_dbg_timer_print
   end interface
   interface sj_timer_stop  ! a short alias 
      module procedure sj_dbg_timer_stop
   end interface

contains
   subroutine sj_dbg_print_time(msg)
      character(30) :: msg
      real*8    :: t
      integer*8 :: c_start =0, c_end
      integer*8 :: c_rate
      real*8    :: elapsed_time
      save c_start
      save c_rate
      
      if(c_start .eq. 0) then
         call system_clock(count=c_start, count_rate=c_rate)
         print '("sj_debug_clockstart= ",i)', c_start
      end if
      if(c_rate .eq. 0) then
         call system_clock(count_rate=c_rate)
         print '("sj_debug_clockrate= ",i)', c_rate
      end if

      call cpu_time(t)
      call system_clock(count=c_end, count_rate=c_rate)
      elapsed_time=real((c_end-c_start),8)/real(c_rate,8)
      print '("sj_debug_time: cpu=[",f15.5, "] real= [", f15.5, "] ", a)', t, elapsed_time, msg
   end subroutine 


   subroutine sj_dbg_timer_start(ntimer)
      integer,intent(in) :: ntimer

      real*8     :: t
      integer*8  :: clock
      real*8     :: elapsed_time

      if ((clock_start(ntimer) .lt. 0) .and. (cpu_start(ntimer) .lt. -1.0d-8) ) then
         call cpu_time(t)
         cpu_start(ntimer)=t
         call system_clock(count=clock, count_rate=clock_rate)
         clock_start(ntimer)=clock
         print '("sj_dbg_timer_start: ntimer=",i4)', ntimer
      else 
         print '("sj_dbg_timer_start: error! ntimer=",i," is already running!")', ntimer
         print '("sj_dbg_timer_start: error! clock_start= ",i," cpu_start= ", f15.5)', clock_start(ntimer), cpu_start(ntimer)
      end if

      if(clock_rate .eq. 0) then
         call system_clock(count_rate=clock_rate)
      end if

   end subroutine


   subroutine sj_dbg_timer_print(ntimer, msg)
      integer,intent(in)                :: ntimer
      character(30),intent(in),optional :: msg

      real*8        :: t
      integer*8     :: clock_end
      real*8        :: elapsed_time

      if ((clock_start(ntimer) .lt. 0) .or. (cpu_start(ntimer) .lt. -1.0d-8) ) then
         print '("sj_dbg_timer_print: error! ntimer=",i4," is not running!")', ntimer
         print '("sj_dbg_timer_print: error! clock_start= ",i4," cpu_start= ", f15.5)', clock_start(ntimer), cpu_start(ntimer)
         return
      end if

      call cpu_time(t)
      t=t-cpu_start(ntimer)
      call system_clock(count=clock_end, count_rate=clock_rate)
      elapsed_time=real((clock_end-clock_start(ntimer)),8)/real(clock_rate,8)
      if (present(msg)) then
         print '("sj_dbg_timer_print: ntimer=",i4 " cpu=[",f15.5, "] real= [", f15.5, "] ", a)', ntimer, t, elapsed_time, msg
      else 
         print '("sj_dbg_timer_print: ntimer=",i4 " cpu=[",f15.5, "] real= [", f15.5, "] "   )', ntimer, t, elapsed_time
      end if

   end subroutine
   

   subroutine sj_dbg_timer_stop(ntimer, msg)
      integer,intent(in)                :: ntimer
      character(30),intent(in),optional :: msg

      real*8        :: t
      integer*8     :: clock_end
      real*8        :: elapsed_time

      if ((clock_start(ntimer) .lt. 0) .or. (cpu_start(ntimer) .lt. -1.0d-8) ) then
         print '("sj_dbg_timer_stop : error! ntimer=",i4," is not running!")', ntimer
         print '("sj_dbg_timer_stop : error! clock_start= ",i4," cpu_start= ", f15.5)', clock_start(ntimer), cpu_start(ntimer)
         return
      end if

      call cpu_time(t)
      t=t-cpu_start(ntimer)
      call system_clock(count=clock_end, count_rate=clock_rate)
      elapsed_time= real((clock_end-clock_start(ntimer)),8)/real(clock_rate,8)
      if (present(msg)) then
         print '("SJ_DBG_TIMER_STOP : ntimer=",I4 " CPU=[",F15.5, "] REAL= [", F15.5, "] ", A)', ntimer, t, elapsed_time, msg
      else 
         print '("SJ_DBG_TIMER_STOP : ntimer=",I4 " CPU=[",F15.5, "] REAL= [", F15.5, "] "   )', ntimer, t, elapsed_time
      end if

      clock_start(ntimer) = -1
      cpu_start(ntimer) = -1.0

   end subroutine

end module sj_timer
