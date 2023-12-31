subroutine calc_compress(dkinfo, denergy, dentropy)
!
! rho and dot{rho} (rhorhs) are in memory 
! output information compressibility (which has the dimension of inverse temperature)
!
! careful: Hamiltonian is time-independent (fixdot = .true.), 
!          otherwise its variation with respect to time needs to be considered
!
! Reference: EPL 85, 40004 (2009) 
!
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8, intent(out) :: dkinfo, denergy, dentropy
integer             :: ni, nj
real*8              :: dtmp1, dtmp2, dtmp3
integer             :: lwork, info, istat
complex*16, allocatable :: work(:)
real*8,     allocatable :: tmpval(:), rwork(:)
!
if (lsparse) then
   write(6,*)
   write(6,*)'calc_compress: calc_compress in SPARSE mode not available yet '
   return  
end if
!
if (.not. fixdot) then
   write(6,*)
   write(6,*)'calc_compress: error! fixdot = .false. not available yet'
   return
end if
!
! denergy = trace(rho * Hamiltonian)
!
cmtmp1(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
!
call zhemm('l', 'u', nrho, nrho, cunity, cmtmp1, nrho, rho(1,1,1), nrho, &
           czero, cmtmp2, nrho)
call zhemm('l', 'u', nrho, nrho, cunity, cmtmp1, nrho, rhorhs(1,1,1), nrho, &
           czero, cmtmp3, nrho)
denergy = 0.d0
dtmp1 = 0.d0
do ni=1,nrho
   denergy = denergy + dble( cmtmp2(ni,ni) )
   dtmp1 = dtmp1 + dble( cmtmp3(ni,ni) )      ! dE/dt
end do
!
! diagonalization of rho
!
cmtmp2(1:nrho,1:nrho) = rho(1:nrho,1:nrho,1)
!
lwork = nrho**2
allocate(work(lwork), rwork(3*nrho-2), tmpval(nrho), STAT=istat)
call zheev('V', 'u', nrho, cmtmp2, nrho, tmpval, work, lwork, rwork, info)
if (info .ne. 0) then
   write(6,*)
   write(6,*)'calc_compress: error! diagnalization failed ', info
   stop
end if
do ni=1,nrho
   if (tmpval(ni) .lt. 0.d0) then
      write(6,*)
      write(6,*)'calc_compress: negative eigenvalue found for rho ', ni, tmpval(ni)
      write(6,*)'calc_compress: code stop '
      stop
   end if
end do
!
! check
!
!cmtmp3(1:nrho,1:nrho) = czero
!do ni=1,nrho
!   cmtmp3(ni,ni) = dcmplx(tmpval(ni),0.d0)
!end do
!!
!call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
!           czero, cmtmp4, nrho)
!call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
!           czero, cmtmp5, nrho)
!cmtmp5(1:nrho,1:nrho) = cmtmp5(1:nrho,1:nrho) - rho(1:nrho,1:nrho,1)
!call cmaxmat(nrho, nrho, cmtmp5, nrho, dtmp2)
!write(6,*)
!write(6,*)'calc_compress: precision of diagnoalization ', dtmp2
!call flush(6)
!
! entropy = -trace( rho * ln(rho) )
!
! eigenvalue=0, xlnx = 0
dentropy = 0.d0
do ni=1,nrho
   if ( dabs(tmpval(ni)) .ge. dpico ) then
       dentropy = dentropy - tmpval(ni) * dlog(tmpval(ni))
   end if 
end do
!
cmtmp3(1:nrho,1:nrho) = 0.d0
do ni=1,nrho
   cmtmp3(ni,ni) = dcmplx(dlog(tmpval(ni)), 0.d0)
end do
call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, cmtmp2, nrho, &
           czero, cmtmp4, nrho)
call zgemm('n', 'c', nrho, nrho, nrho, cunity, cmtmp4, nrho, cmtmp2, nrho, &
           czero, cmtmp3, nrho)  ! cmtmp3 = dlog(rho)
do ni=1,nrho
   cmtmp3(ni,ni) = cmtmp3(ni,ni) + cunity
end do
call zhemm('r', 'u', nrho, nrho, cunity, cmtmp3, nrho, rhorhs(1,1,1), nrho, &
           czero, cmtmp4, nrho)
dtmp2 = 0.d0
do ni=1,nrho
   dtmp2 = dtmp2 - dble(cmtmp4(ni,ni))  ! dS/dt
end do
!
dkinfo = dtmp2 / dtmp1
!
deallocate(work, rwork, tmpval, STAT=istat)
!
return
end subroutine calc_compress
