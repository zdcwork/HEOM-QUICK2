!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug utilities developed by SUN Jian
! for benchmarking the time consuming of each part of the code
!
! this module contains:
!    Subroutines:
!      SJ_DBG_PRINT_TIME
!         print out time infomation together with the msg 
!      SJ_DBG_TIMER_START
!         start the timer (nTimer)
!      SJ_DBG_TIMER_STOP
!         stop the timer (nTimer) and report the elapsed time and msg (optional)
!
!    Aliases to the subroutines listed above
!      SJ_PRINT_TIME
!         alias for SJ_DBG_PRINT_TIME
!      SJ_CPU_TIME 
!         alias for SJ_DBG_PRINT_TIME, not recommended
!      SJ_START_TIMER
!         alias for SJ_DBG_START_TIMER
!      SJ_STOP_TIMER
!         alias for SJ_DBG_STOP_TIMER
!
!    parameters:
!      msg:    the message you want to print out after the time 
!      nTimer: the identifier of your timer, ranged from 1 to 255
!              required parameter for SJ_DBG_TIMER_START and SJ_DBG_TIMER_STOP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SJ_DBG
   implicit none

   private  clock_start, cpu_start, clock_rate
   public   SJ_CPU_TIME, SJ_DBG_PRINT_TIME, SJ_DBG_TIMER_START, SJ_DBG_TIMER_STOP

   INTEGER :: clock_start(255)=-1
   REAL*8  :: cpu_start(255)=-1.0
   INTEGER :: clock_rate =0

   interface SJ_PRINT_TIME  ! a short alias 
      module procedure SJ_DBG_PRINT_TIME
   end interface
   interface SJ_CPU_TIME    ! an alias for compatible consideration, not recommended
      module procedure SJ_DBG_PRINT_TIME
   end interface
   interface SJ_TIMER_START ! a short alias 
      module procedure SJ_DBG_TIMER_START
   end interface
   interface SJ_TIMER_PRINT ! a short alias 
      module procedure SJ_DBG_TIMER_PRINT
   end interface
   interface SJ_TIMER_STOP  ! a short alias 
      module procedure SJ_DBG_TIMER_STOP
   end interface

contains
   subroutine SJ_DBG_PRINT_TIME(msg)
      character(30) :: msg
      real*8 :: t
      INTEGER :: c_start =0,c_end, c_rate
      REAL(KIND=8) :: elapsed_time
      save c_start
      save c_rate
      
      if(c_start .eq. 0) then
         call system_clock(count=c_start)
         print '("SJ_DEBUG_clockstart= ",I)', c_start
      end if
      if(c_rate .eq. 0) then
         call system_clock(COUNT_RATE=c_rate)
         print '("SJ_DEBUG_clockRate= ",I)', c_rate
      end if

      call cpu_time(t)
      call system_clock(count=c_end)
      elapsed_time=REAL((c_end-c_start),8)/REAL(c_rate,8)
      print '("SJ_DEBUG_TIME: CPU=[",f15.5, "] REAL= [", f15.5, "] ", A)', t, elapsed_time, msg
   end subroutine 


   subroutine SJ_DBG_TIMER_START(nTimer)
      INTEGER,intent(in) :: nTimer

      real*8   :: t
      INTEGER  :: clock
      REAL*8   :: elapsed_time

      if ((clock_start(nTimer) .lt. 0) .and. (cpu_start(nTimer) .lt. -1.0d-8) ) then
         call cpu_time(t)
         cpu_start(nTimer)=t
         call system_clock(count=clock)
         clock_start(nTimer)=clock
         print '("SJ_DBG_TIMER_START: nTimer=",I4)', nTimer
      else 
         print '("SJ_DBG_TIMER_START: Error! nTimer=",I," is already running!")', nTimer
         print '("SJ_DBG_TIMER_START: Error! clock_start= ",I," cpu_start= ", F15.5)', clock_start(nTimer), cpu_start(nTimer)
      end if

      if(clock_rate .eq. 0) then
         call system_clock(COUNT_RATE=clock_rate)
      end if

   end subroutine


   subroutine SJ_DBG_TIMER_PRINT(nTimer, msg)
      INTEGER,intent(in)                :: nTimer
      character(30),intent(in),optional :: msg

      real*8        :: t
      INTEGER       :: clock_end
      REAL*8        :: elapsed_time

      if ((clock_start(nTimer) .lt. 0) .or. (cpu_start(nTimer) .lt. -1.0d-8) ) then
         print '("SJ_DBG_TIMER_PRINT: Error! nTimer=",I4," is not running!")', nTimer
         print '("SJ_DBG_TIMER_PRINT: Error! clock_start= ",I4," cpu_start= ", F15.5)', clock_start(nTimer), cpu_start(nTimer)
         return
      end if

      call cpu_time(t)
      t=t-cpu_start(nTimer)
      call system_clock(count=clock_end)
      elapsed_time=REAL((clock_end-clock_start(nTimer)),8)/REAL(clock_rate,8)
      if (present(msg)) then
         print '("SJ_DBG_TIMER_PRINT: nTimer=",I4 " CPU=[",f15.5, "] REAL= [", f15.5, "] ", A)', nTimer, t, elapsed_time, msg
      else 
         print '("SJ_DBG_TIMER_PRINT: nTimer=",I4 " CPU=[",f15.5, "] REAL= [", f15.5, "] "   )', nTimer, t, elapsed_time
      end if

   end subroutine
   

   subroutine SJ_DBG_TIMER_STOP(nTimer, msg)
      INTEGER,intent(in)                :: nTimer
      character(30),intent(in),optional :: msg

      real*8        :: t
      INTEGER       :: clock_end
      REAL*8        :: elapsed_time

      if ((clock_start(nTimer) .lt. 0) .or. (cpu_start(nTimer) .lt. -1.0d-8) ) then
         print '("SJ_DBG_TIMER_STOP : Error! nTimer=",I4," is not running!")', nTimer
         print '("SJ_DBG_TIMER_STOP : Error! clock_start= ",I4," cpu_start= ", F15.5)', clock_start(nTimer), cpu_start(nTimer)
         return
      end if

      call cpu_time(t)
      t=t-cpu_start(nTimer)
      call system_clock(count=clock_end)
      elapsed_time=REAL((clock_end-clock_start(nTimer)),8)/REAL(clock_rate,8)
      if (present(msg)) then
         print '("SJ_DBG_TIMER_STOP : nTimer=",I4 " CPU=[",f15.5, "] REAL= [", f15.5, "] ", A)', nTimer, t, elapsed_time, msg
      else 
         print '("SJ_DBG_TIMER_STOP : nTimer=",I4 " CPU=[",f15.5, "] REAL= [", f15.5, "] "   )', nTimer, t, elapsed_time
      end if

      clock_start(nTimer) = -1
      cpu_start(nTimer) = -1.0

   end subroutine

end module
