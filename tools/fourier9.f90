program fourier
!
! purpose: do Fourier trnasformation for a signal
!
! note on 23Jun2008: Fourier and inverse Fourier transform (discrete FT)
!
   implicit none
   integer,   parameter :: nout = 3000
!   integer,   parameter :: nout = 20000
   integer*8, parameter :: ndata=20000000
   real*8, parameter  :: hbar = 0.65822d0
   real*8, parameter  :: pi = 3.1415926535897932384626d0
!
   integer :: ntimes, istat
   common /isignals/ ntimes
   real*8  :: t(ndata), signal(ndata), delta
   common /signals/ t, signal, delta
   real*8  :: da(ndata,2), freqm(ndata,2), dm(ndata,2)
!
!   real*8  :: freq(ndata), frsignal(ndata), fisignal(ndata)
   real*8  :: freq(nout), frsignal(nout), fisignal(nout)
   real*8, external :: fourier1
!
   integer :: i,j,k, n, nfreq, itimes, nread
   real*8  :: dt0, dt1, omega, dtmp1, dtmp2, tlas
   real*8  :: x1,x2,x3,f1,f2,f3,a1,b1,c1,xmax,fmax
   real*8  :: aa,bb,cc,dd,xx,yy,zz,y0
!
   logical :: ijobfreq, lfirst
!
   integer, parameter :: mgauleg1=8
double precision :: dpt(mgauleg1), dw(mgauleg1)
data dpt    /                &
             -0.96028986d0,  &
             -0.79666648d0,  &
             -0.52553241d0,  &
             -0.18343464d0,  &
              0.18343464d0,  &
              0.52553241d0,  &
              0.79666648d0,  &
              0.96028986d0 /
data dw     /                &
              0.10122854d0,  &
              0.22238103d0,  &
              0.31370665d0,  &
              0.36268378d0,  &
              0.36268378d0,  &
              0.31370665d0,  &
              0.22238103d0,  &
              0.10122854d0 /

!
!   write(*,*)' ijobfreq (true or false)'
!   read(*,*)ijobfreq
   ijobfreq = .false.
!
!   write(*,*)' ntimes, itimes, tlas '
!   read(*,*) ntimes, itimes, tlas
   itimes = 1
   write(*,*)' ntimes, tlas '
   read(*,*) ntimes, tlas
   if (ntimes>ndata) then
      write(6,*) 'Too many data!', ntimes, ndata
      stop
   end if
   write(6,*) 'number of points:',ntimes,' per ',itimes
   write(6,*) 'tlas:', tlas
!
   lfirst = .true.
   nread = ntimes
   open(unit=47, file='signal.dat', form='formatted', status='old')
   rewind(47)
   i = 1
   do j = 1, ntimes
      read(47, *, iostat=istat) t(i), signal(i)
      if (istat .eq. 0) then
        signal(i) = signal(i)*dexp(-tlas*t(i)/hbar)
      else
        if (lfirst) then
          write(6,*) 
          write(6,*)' input data ending at ', i
          lfirst = .false.
          nread = i
        end if
!
        if (i .le. 2) then
          write(6,*) ' error reading signal.dat '
          stop
        else
          t(i) = 2.d0 * t(i-1) - t(i-2)
          signal(i) = 0.d0
        end if
      end if
!
      if (itimes==1.or.mod(j,itimes)==1) then
         i = i + 1
         if (i==3) then
            delta = t(2) - t(1)
         elseif(i>3) then
            dt1 = t(i-1) - t(i-2)
            if (dabs(dt1-delta)>1.d-5) then
               write(6,*) 'The time step is not even!', i, dt1, delta
               write(6,*) 'The present code will stop!'
               stop
            end if
         end if
      end if
   end do
   close (47)
   write(6,*)
   write(6,*)' check dephasing constant, end of signal '
   write(6,*)' without ', signal(ntimes)*dexp(tlas*t(ntimes)/hbar), ' with ', signal(ntimes)
   ntimes = i-1
   if (mod(ntimes,2)/=0) then
      write(6,*) 'Warning: the input data number  : ',ntimes
      write(6,*) 'the number of data actually used: ',ntimes-1
      ntimes = ntimes -1
   end if
   write(6,*) 'the number of data actually used: ',ntimes
   dt1 = t(1)
   write(6,*) 'initial time :', t(1)
   write(6,*) 'final time   :', t(ntimes)
   do i = 1, ntimes
      t(i) = t(i) - dt1
!     write(6,*) 'ti,signa',i,t(i),signal(i)
!     signal(i) = signal(i)*dexp(-0.1d0*t(i))
   end do
   delta = t(2) - t(1)
   write(6,*) 'time step    :', delta
!   do n = 0, ntimes/2
   do n = 0, min(ntimes/2,nout-1)
      i = n + 1
      frsignal(i) = 0.d0
      fisignal(i) = 0.d0
      omega = 2.d0*pi*n/(ntimes*delta)
      freq(i) = omega
      !do j = 2, ntimes
      do j = 2, nread
         !frsignal(i) = frsignal(i) + signal(j)*dcos(omega*t(j))
         !fisignal(i) = fisignal(i) + signal(j)*dsin(omega*t(j))  ! exp^(iwt)
         aa = t(j-1)
         bb = t(j)
         y0 = signal(j-1)
         cc = (bb - aa) * 5.d-1
         dd = (bb + aa) * 5.d-1
         zz = (signal(j) - y0) / (bb - aa)
         do k=1,mgauleg1
            xx = cc * dpt(k) + dd
            yy = zz * (xx - aa) + y0
            frsignal(i) = frsignal(i) + yy * dcos(omega*xx) * dw(k) * 5.d-1
            fisignal(i) = fisignal(i) + yy * dsin(omega*xx) * dw(k) * 5.d-1
         end do
      end do
!     write(6,*) 'n,omega,fr,fi',n,omega, frsignal(i), fisignal(i)
   end do
   nfreq = min(ntimes/2 + 1, nout)
   write(6,*) 'maxiam frequency(eV):', freq(nfreq)*hbar
   write(6,*) 'frequency step(eV)  :', freq(2)*hbar
!
   if (ijobfreq) then
     open(99,file='sigfreq.dat')
     rewind(99)
     do i=1,nfreq
       read(99,*)dtmp1, frsignal(i), fisignal(i)
     end do
     close(99)
     open(99,file='sigtime.dat')
!
     do i = 1, ntimes
        dtmp1 = 0.d0
        do n = 1, nfreq
           omega = freq(n)
           dtmp2 = frsignal(n)*dcos(omega*t(i))    &
                   + fisignal(n)*dsin(omega*t(i))      ! exp^{iwt}
!                  - fisignal(n)*dsin(omega*t(i))      ! exp^(-iwt)
           dtmp1 = dtmp1 + dtmp2
           if (n/=1.and.n/=nfreq) dtmp1 = dtmp1 + dtmp2
        end do
        dtmp1 = dtmp1/dble(ntimes)/delta
!        write(99,'(4f18.10)') t(i), dtmp1
!        write(99,'(4f18.10)') t(i), dtmp1*dexp(tlas*t(i)/hbar)
        write(99,518) t(i), dtmp1*dexp(tlas*t(i)/hbar)

!        if (dabs(dtmp1-signal(i))>1.d-5) then
!          write(6,*) 'Fourier transformation may fail!'
!          write(6,*) 'i,dtmp1, s',i,dtmp1,signal(i)
!        end if
     end do
!
     stop
   endif
   
!   write(6,*)'nfreq, ntimes', nfreq, ntimes
!   stop 

   j = 0
!   do i = 1, nfreq
   do i=1, min(nfreq,nout)
!     frsignal(i) = frsignal(i)/dexp(-0.01d0*freq(i)*freq(i)/4.d0)
!     fisignal(i) = fisignal(i)/dexp(-0.01d0*freq(i)*freq(i)/4.d0)
      frsignal(i) = frsignal(i)*delta
      fisignal(i) = fisignal(i)*delta
      dtmp1 = dsqrt(frsignal(i)**2 + fisignal(i)**2)  
!      if (freq(i)*hbar>=50.d0) cycle
!      write(6,'(4f20.6)') freq(i)*hbar, frsignal(i), fisignal(i), dtmp1
      write(6,518) freq(i)*hbar, frsignal(i), fisignal(i), dtmp1
!      f1 = f2
!      f2 = f3
!      f3 = dtmp1
!      if (i>=3) then
!         if (f2>f1.and.f2>f3) then
!            x1 = freq(i-2)
!            x2 = freq(i-1)
!            x3 = freq(i)
!            a1 = ((f3-f1)/(x3-x1)-(f2-f1)/(x2-x1))/(x3-x2)
!            b1 = (f2-f1)/(x2-x1)-a1*x2
!            j = j + 1
!            freqm(j,1) =  -(b1-a1*x1)/(2.d0*a1)
!            dm(j,1)    = (freqm(j,1)-x1)*(a1*freqm(j,1)+b1)+f1
!            da(j,1)    = a1
!!            fmax=-golden(x1, x2, x3, fourier1, 1.d-5, xmax)
!            freqm(j,2) = xmax
!            dm(j,2) = fmax
!            x1 = fmax
!            x2 = -fourier1(xmax+0.001d0)
!            x3 = -fourier1(xmax-0.001d0)
!            da(j,2) = (x3+x2-2.d0*x1)*1.d6
!!           write(6,*) 'fmax, dm', fmax, dm(j,1)
!!           write(6,*) 'x1,x2,x3',x1,x2,x3
!!           write(6,*) 'f1,f2,f3',f1,f2,f3
!!           write(6,*) 'maxium:', -(b1-a1*x1)/(2.d0*a1)
!!           write(6,*) (x2-x1)*(a1*x2+b1)+f1-f2
!!           write(6,*) (x3-x1)*(a1*x3+b1)+f1-f3
!!           write(6,*) (freqm(j)-x1)*(a1*freqm(j)+b1)+f1
!         end if
!      end if
   end do
!   write(6,*) 'Peaks:' 
!   do i = 1, j
!      write(6,'(1x,i4,6f16.10)') i, freqm(i,2)*hbar, dm(i,2), da(i,1)
!   end do
!
518 format(f18.10, 2x, 3(e18.10e3, 2x))

end program fourier
!
real*8 function fourier1(omega)
! 
   implicit none
   integer, parameter :: ndata=50000
   real*8, parameter  :: hbar = 0.65822d0
   real*8, parameter  :: pi = 3.1415926535897932384626d0
!
   real*8, intent(in) :: omega
!
   integer :: ntimes
   common /isignals/ ntimes
   real*8  :: t(ndata), signal(ndata), delta
   common /signals/ t, signal, delta
!
   integer :: i,j,k
   real*8  :: dtmp1, dtmp2

   dtmp1 = 0.d0
   dtmp2 = 0.d0
   do j = 1, ntimes
      dtmp1 = dtmp1 + signal(j)*dcos(omega*t(j))
      dtmp2 = dtmp2 + signal(j)*dsin(omega*t(j))
   end do
!  dtmp1 = dtmp1/dexp(-0.01d0*omega*omega/4.d0)
!  dtmp2 = dtmp2/dexp(-0.01d0*omega*omega/4.d0)
   dtmp1 = dtmp1*delta
   dtmp2 = dtmp2*delta
   fourier1 =-dsqrt(dtmp1*dtmp1+dtmp2*dtmp2)
!
   return
end function fourier1

