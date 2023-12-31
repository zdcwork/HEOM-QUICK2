subroutine zpout_pfd(N, ns, z, zcoef, sqza)
!
! purpose : find poles and coefficients for PFD scheme
!
! f(z) = 0.5 + sum_i zcoef * 2 * z / (z^2 - zpole^2)
! with z = beta * w
! 
implicit none
include "../include/sizes"
include "../include/common"
!
integer, intent (in)     :: N, ns
complex*16, intent (out) :: z(*), zcoef(*), sqza(N,*)
integer                  :: mm, nn, nl, nm, np, nu, pp, istat
integer                  :: LDZZ                     ! dimension lengrth of ZZ array
integer                  :: LDvl, LDvr, info, lwork  ! parameter used for zgeev
character*1              :: jobvl, jobvr                  
real*8,     allocatable  :: delta(:,:)               ! delta function
real*8,     allocatable  :: rwork(:)                                                                                            
complex*16, allocatable  :: ZZ(:,:)                  ! ZZ function
complex*16, allocatable  :: sqz(:)                   ! alpha and gama are pole including fore-exp and up-exp cofficients 
complex*16, allocatable  :: vl(:,:), vr(:,:)         ! eigenvector of ZZ matrix
complex*16, allocatable  :: work(:)                  ! zgeev workspace
!                                                                                                                                                                
! delta function defined   
allocate(sqz(N), delta(N,N), ZZ(N,N), STAT=istat)
! 
do mm=1,N     
  do nn=1,N  
    delta(mm,nn) = 0.d0      
    ZZ(mm,nn)    = czero
  end do                 
  delta(mm,mm) = 1.d0
end do                                                                                                                                                          
! ZZ establish 
do mm=1,N  
  do nn=1,N 
    if (mm+1 .le. N) ZZ(mm,nn) = dble(2 * mm * (2 * mm - 1)) * delta(nn,mm+1) 
    ZZ(mm,nn) = ZZ(mm,nn) - dble(2 * N * (2 * N - 1)) * delta(mm, N)         
  end do                                                                    
end do    
!
jobvl = 'N' 
jobvr = 'N' 
LDZZ  = N   
LDvl  = N   
LDvr  = N   
lwork = 2 * N   
allocate(vl(N,N), vr(N,N), STAT=istat)         
allocate(work(lwork), rwork(lwork), STAT=istat) 
call zgeev(jobvl, jobvr, N, ZZ, LDZZ, z, vl, LDvl, vr, LDvr, work, lwork, rwork, info) 
if (info .ne. 0) then 
  write(6,*)'zpout_pfd: error! diagonalization of ZZ matrix failed' 
  stop
end if  
!
do pp=1,N
  sqz(pp) = cdsqrt(z(pp)) * 2.d0
end do     
! 
do pp=1,N
   if ( dimag(sqz(pp)) .gt. 0.d0 ) then
      sqza(pp,1) = sqz(pp) 
   else 
      sqza(pp,1) = -sqz(pp)
   end if   
   sqza(pp,2) = dconjg(sqza(pp,1))
end do
!
do pp=1,N
   zcoef(pp) = -cunity
end do
!
deallocate(delta, rwork, STAT=istat)
deallocate(ZZ, sqz, vl, vr, work, STAT=istat)
!
end subroutine
