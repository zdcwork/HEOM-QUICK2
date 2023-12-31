subroutine get_omp_env
use omp_lib
use matmod_omp
implicit none
!
include '../include/sizes'
!
integer :: tid
!
!$omp parallel default(shared) private(tid) 
tid = omp_get_thread_num()
!
if (tid .eq. 0) then
   nthreads_omp = omp_get_num_threads()
   nprocs_omp   = omp_get_num_procs()
   maxt_omp     = omp_get_max_threads()
!
   write(6,*)'get_omp_env: number of threads = ', nthreads_omp
   write(6,*)'get_omp_env: number of procs   = ', nprocs_omp
   write(6,*)'get_omp_env: max thread number = ', maxt_omp
!
   if (maxprocs < nthreads_omp .or. maxprocs < nprocs_omp) then
      write(6,*)'get_omp_env: error thread number '
      write(6,*)'get_omp_env: maxprocs, nthreads, nprocs ',     &
                maxprocs, nthreads_omp, nprocs_omp
      stop
   end if
end if
!$omp end parallel
write(6,*)'get_omp_env: check nthreads_omp ', nthreads_omp
if (nthreads_omp > 1) then
   write(6,*)'get_omp_env: The code is running in parallel mode '
else
   write(6,*)'get_omp_env: The code is running in serial mode '
end if
call flush(6)
!
end subroutine get_omp_env
