subroutine liouamsmat_omp(ialr, irhont, iorbs, ispin, isgn, alpha, zrho, lrho, &
                          const, zout, lout, ithread)
use matmod
use tmpmatmod_omp
implicit none
include '../include/sizes'
include '../include/common'
!
! purpose : calculation in the Liouville space associated with
!           Y=A*rho or Y=rho*A or Y=A*rho^\dag or Y=rho^\dag*A. However, the 
!           Liouville coefficient matrix L_A should be transposed. 
!           A and rho are complex.
!
!           zout = const * zout + alpha * L^t_(-i*ams) * zrho
!
!           ialr   = 1 : A to the left of rho
!                  = 2 : A to the right of rho
!          
!           irhont = 0 : rho is in normal form
!                 != 0 : rho is in its Hermitian conjugated form
!
integer, intent(in)       :: ialr, irhont, iorbs, ispin, isgn, lrho, lout, ithread
real*8,     intent(in)    :: alpha
complex*16, intent(in)    :: const
complex*16, intent(in)    :: zrho(lrho,*)
complex*16, intent(inout) :: zout(lout,*)
!
integer    :: ni, nj, nn, nm, n1, n2
real*8     :: dtmp0
complex*16 :: ctmp0
!
zmtmp1_omp(1:nrho,1:nrho,ithread) = czero
ctmp0  = eye * alpha
!
if (ialr .eq. 1 .and. irhont .eq. 0) then  ! Y=A*rho
  do n1=1,nams
    if (isgn .eq. 1) then 
      ni = rowams(n1,iorbs,ispin)
      nm = colams(n1,iorbs,ispin)
    else
      ni = colams(n1,iorbs,ispin)
      nm = rowams(n1,iorbs,ispin)
    end if
    dtmp0 = valams(n1,iorbs,ispin)
    do nj=1,nrho
       zmtmp1_omp(ni,nj,ithread) = zmtmp1_omp(ni,nj,ithread) + dtmp0 * zrho(nm,nj)
    end do
  end do
!
else if (ialr .eq. 1 .and. irhont .ne. 0) then  ! Y=A*rho^\dag
  do n2=1,nams
    if (isgn .eq. 1) then
      nj = rowams(n2,iorbs,ispin)
      nm = colams(n2,iorbs,ispin)
    else
      nj = colams(n2,iorbs,ispin)
      nm = rowams(n2,iorbs,ispin)
    end if
    dtmp0 = valams(n2,iorbs,ispin)
    do ni=1,nrho
       zmtmp1_omp(ni,nj,ithread) = zmtmp1_omp(ni,nj,ithread) - dtmp0 * dconjg(zrho(nm,ni))
    end do
  end do
!
else if (ialr .eq. 2 .and. irhont .eq. 0) then  ! Y=rho*A
   do n2=1,nams
     if (isgn .eq. 1) then
       nj = colams(n2,iorbs,ispin)
       nn = rowams(n2,iorbs,ispin)
     else
       nj = rowams(n2,iorbs,ispin)
       nn = colams(n2,iorbs,ispin)
     end if
     dtmp0 = valams(n2,iorbs,ispin)
     do ni=1,nrho
        zmtmp1_omp(ni,nj,ithread) = zmtmp1_omp(ni,nj,ithread) + dtmp0 * zrho(ni,nn)
     end do
   end do
!
else if (ialr .eq. 2 .and. irhont .ne. 0) then  ! Y=rho^\dag *A
   do n1=1,nams
      if (isgn .eq. 1) then
         ni = colams(n1,iorbs,ispin)
         nn = rowams(n1,iorbs,ispin)
      else
         ni = rowams(n1,iorbs,ispin)
         nn = colams(n1,iorbs,ispin)
      end if
      dtmp0 = valams(n1,iorbs,ispin)
      do nj=1,nrho
         zmtmp1_omp(ni,nj,ithread) = zmtmp1_omp(ni,nj,ithread) - dtmp0 * dconjg(zrho(nj,nn))
      end do
   end do  
!
else
   write(6,*)'liouamsmat_omp: error! unknown ialr = ', ialr, ' irhont = ', irhont
   stop
end if
!
zout(1:nrho, 1:nrho) = zout(1:nrho, 1:nrho) * const +                 &
                       zmtmp1_omp(1:nrho, 1:nrho, ithread) * ctmp0
!
end subroutine liouamsmat_omp
