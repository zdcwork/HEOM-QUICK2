subroutine lioumatmat(ialr, irhont, ireala, iimaga, nrho, zamat, lamat, zrho, lrho, &
                      const, zout, lout) 
implicit none
!
! purpose : calculation in the Liouville space associated with
!           Y=A*rho or Y=rho*A or Y=A*rho^\dag or Y=rho^\dag*A. However, the 
!           Liouville coefficient matrix L_A should be transposed. 
!           A and rho are complex.
!
!           zout = const * zout + L^t_zamat * zrho
!
!           ialr   = 1 : A to the left of rho
!                  = 2 : A to the right of rho
!          
!           irhont = 0 : rho is in normal form
!                 != 0 : rho is in its Hermitian conjugated form
!
!  This version disregards the two inputs 'ireala' and 'iimaga':
!           ireala = 0 : zamat has zero real part
!                 != 0 : zamat has nonzero real part
!
!           iimaga = 0 : zamat has zero imag part
!                 != 0 : zamat has nonzero imag part
!
integer, intent(in)       :: ialr, irhont, ireala, iimaga, lamat, nrho, lrho, lout
complex*16, intent(in)    :: const
complex*16, intent(in)    :: zamat(lamat,*), zrho(lrho,*)
complex*16, intent(inout) :: zout(lout,*)
!
integer    :: ni, nj, nn, nm
complex*16 :: ctmp0
!
zout(1:nrho, 1:nrho) = zout(1:nrho, 1:nrho) * const
!
if (ialr .eq. 1 .and. irhont .eq. 0) then  ! Y=A*rho
  !do nj=1,nrho
  !  do ni=1,nrho
  !    if (ireala .ne. 0) then
  !      do nm=1,nrho
  !        zout(ni,nj) = zout(ni,nj) + dble(zamat(nm,ni)) * zrho(nm,nj)
  !      end do
  !    end if
  !    if (iimaga .ne. 0) then
  !      do nm=1,nrho
  !        zout(ni,nj) = zout(ni,nj) + dcmplx(0.d0, -dimag(zamat(nm,ni))) * zrho(nm,nj)
  !      end do
  !    end if
  !  end do
  !end do
  !return
  ! 
   do nj=1,nrho
      do ni=1,nrho
         do nm=1,nrho
            zout(ni,nj) = zout(ni,nj) + dconjg(zamat(nm,ni)) * zrho(nm,nj)
         end do
      end do
   end do
   return
!
else if (ialr .eq. 1 .and. irhont .ne. 0) then  ! Y=A*rho^\dag
  !do nj=1,nrho
  !  do ni=1,nrho
  !    if (ireala .ne. 0) then
  !      do nm=1,nrho
  !        zout(ni,nj) = zout(ni,nj) + dble(zamat(nm,nj)) * dconjg(zrho(nm,ni))
  !      end do
  !    end if
  !    if (iimaga .ne. 0) then
  !      do nm=1,nrho
  !        zout(ni,nj) = zout(ni,nj) + dimag(zamat(nm,nj)) * dcmplx(dimag(zrho(nm,ni)), &
  !                                                                  dble(zrho(nm,ni)))
  !      end do
  !    end if
  !  end do
  !end do
  !return
   do nj=1,nrho
      do ni=1,nrho
         do nm=1,nrho
            zout(ni,nj) = zout(ni,nj) + zamat(nm,nj) * dconjg(zrho(nm,ni))
         end do
      end do
   end do
   return
!
else if (ialr .eq. 2 .and. irhont .eq. 0) then  ! Y=rho*A
   !do nj=1,nrho
   !  do ni=1,nrho
   !    if (ireala .ne. 0) then
   !      do nn=1,nrho
   !        zout(ni,nj) = zout(ni,nj) + dble(zamat(nj,nn)) * zrho(ni,nn)
   !      end do
   !    end if
   !    if (iimaga .ne. 0) then
   !      do nn=1,nrho
   !        zout(ni,nj) = zout(ni,nj) + dcmplx(0.d0, -dimag(zamat(nj,nn))) * zrho(ni,nn)
   !      end do
   !    end if
   !  end do
   !end do
   !return
    do nn=1,nrho
       do nj=1,nrho
          ctmp0 = dconjg(zamat(nj,nn))
          do ni=1,nrho
             zout(ni,nj) = zout(ni,nj) + ctmp0 * zrho(ni,nn)
          end do
       end do
    end do
    return
!
else if (ialr .eq. 2 .and. irhont .ne. 0) then  ! Y=rho^\dag *A
   !do nj=1,nrho
   !  do ni=1,nrho
   !    if (ireala .ne. 0) then
   !      do nn=1,nrho
   !        zout(ni,nj) = zout(ni,nj) + dble(zamat(ni,nn)) * dconjg(zrho(nj,nn))
   !      end do
   !    end if
   !    if (iimaga .ne. 0) then
   !      do nn=1,nrho
   !        zout(ni,nj) = zout(ni,nj) + dimag(zamat(ni,nn)) * dcmplx(dimag(zrho(nj,nn)), &
   !                                                                  dble(zrho(nj,nn)))
   !      end do
   !    end if
   !  end do
   !end do
   !return
    do nn=1,nrho
       do nj=1,nrho
          ctmp0 = dconjg(zrho(nj,nn))
          do ni=1,nrho
             zout(ni,nj) = zout(ni,nj) + ctmp0 * zamat(ni,nn)
          end do
       end do
    end do
    return
!
else
   write(6,*)
   write(6,*)' error in lioumatmat, unknown ialr= ', ialr
   stop
end if
!
end subroutine lioumatmat
