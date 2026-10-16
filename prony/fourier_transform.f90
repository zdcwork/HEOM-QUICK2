! Fortran example.
! 1D complex to complex, and real to conjugate-even
subroutine fourier_transform(isgn_tmp,T_tmp,ialf,iorbs) 
!
! use MKL_DFT_TYPE
! use MKL_DFTI
 use fft_module
 use cubic_spline_interp
 implicit none
!
include "../include/sizes"
include "../include/common"
!
 integer, intent(in) :: ialf, iorbs
 real*8, intent(in) :: T_tmp, isgn_tmp
 real*8 :: T, dww, dtt
 integer :: isgn, isgn_t
 integer :: istat
 integer :: i, j, n_middle, p_nsample
 complex*16, allocatable :: XC_in(:), XC_out(:)
 complex*16, allocatable :: XC_remain(:), XC_interp(:)
 complex*16 :: ctmp
 real*8, allocatable :: tt_input(:), X(:), fermi(:)
 real*8, allocatable :: tt_remain(:)
 real*8, allocatable :: tt_interp(:)
! type(DFTI_DESCRIPTOR), pointer :: My_Desc_Handle
!
 integer, allocatable :: num_coloum(:), ntmp(:)
 real*8  :: dtmp, dtmp1, dtmp2, sam_freq
 complex*16, parameter :: ceye = dcmplx(0.d0,1.0d0)
 logical :: lfile_exist
 character(len=128) :: buffer
 character(len=128) :: ch_tmp
 character(len=2)   :: ch_ialf, ch_isgn, ch_iorbs
 character(len=128) :: title
 namelist / fourier  / lread_data, lcal_fermi,&
                       tol, tf_interp,&
                       linterp, num_interp
!
 i = 1
 p_nsample = -1
 lread_data = .TRUE.
 lcal_fermi = .TRUE.
 !
 isgn = int(isgn_tmp)
 T = T_tmp
 !scaling = 1.d0
 tol = 1.0d-6
 linterp = .FALSE.
 num_interp = 5001
 tf_interp = 1.d0/T
 !
 open(unit=5,status='unknown')
 rewind(5)
 read(UNIT=5,NML=fourier,END=100)
 100 continue
 close(5)
 write(6,'(" fourier_transform: invoke for isgn",I2," ialf",I2," iorbs",I2)') isgn,ialf,iorbs
 !
 !if (dabs(scaling-1.d0) .gt. tol) then
 !  write(6,*) 'scaling factor is', scaling
 !end if
 !if (scaling .ne. 1.d0) then
 !  write(6,*) 'warning: please set scaling = 1.d0',&
 !             ' for (inverse) fast Fourier transform (FFT)'
 !  stop
 !end if
 write(6,*)
 !
 !T = T * scaling
 !
 if (lread_data) then
   inquire(file='FT_INPUT.data',exist=lfile_exist)
   if (lfile_exist) then
     open(unit=10,file='FT_INPUT.data')
   else
     write(6,*) 'fourier_transform: does not find the input file FT_INPUT.data'
     write(6,*) 'fourier_transform: please check your input files'
     stop
   end if
   !
   rewind(10)
   do while(.not. is_iostat_end(istat))
      read(10,'(a)', iostat=istat) buffer
      p_nsample = p_nsample + 1
   end do
   !
   write(6,*) 'fourier_transform: the size of input data', p_nsample
   allocate(X(p_nsample), &
            num_coloum(p_nsample), ntmp(len(buffer)))
   rewind(10)
   ! check the number of coloum
   num_coloum = 0
   do i=1, p_nsample
      read(10,'(a)') buffer
      ntmp = 0
      if (buffer(1:1) .ne. char(32)) then
         num_coloum(i) = 1
      end if
      do j=1, len(buffer)
         if(buffer(j:j) .eq. char(32) .or.&
            buffer(j:j) .eq. char(9)  .or.&
            buffer(j:j) .eq. char(44)) then
           ntmp(j) = 1
         end if
      end do
      do j=1, len(buffer)-1
         if(ntmp(j) .eq. 1 .and. ntmp(j+1) .eq. 0) then
            num_coloum(i) = num_coloum(i) + 1
         end if
      end do
      !write(6,'(24(I1))') ntmp(1:24)
   end do
   if (maxval(num_coloum) .ne. minval(num_coloum)) then
      write(6,*) 'fourier_transform: the coloum of input data is not consist'
      write(6,*) 'fourier_transform: max number of coloum', maxval(num_coloum)
      write(6,*) 'fourier_transform: min number of coloum', minval(num_coloum)
      stop
   end if
   !
   write(6,*) 'fourier_transform: the couloum of input data', maxval(num_coloum)
   rewind(10)
   !
   !
   allocate(XC_in(p_nsample))
   if (maxval(num_coloum) .eq. 2) then ! real input
     write(6,*) 'fourier_transform: the type of input data is real'
     do i=1, p_nsample
        dtmp1 = 0.d0
        read(10,*) X(i), dtmp1
        XC_in(i) = dcmplx(dtmp1, 0.0d0)
     end do
   else if (maxval(num_coloum) .eq. 3) then ! complex input
     write(6,*) 'fourier_transform: the type of input data is complex'
     do i=1, p_nsample
        dtmp1 = 0.d0
        dtmp2 = 0.d0
        read(10,*) X(i), dtmp1, dtmp2
        XC_in(i) = dcmplx(dtmp1, dtmp2)
     end do
   else
     write(6,*) 'fourier_transform: wrong coloum'
     write(6,*) 'fourier_transform: max value of input coloum', maxval(num_coloum)
     stop
   end if
   rewind(10)
   close(10)
   !
   ! check the gap among two neighbour elements
   dtmp1 = (X(2) - X(1))
   do i=3, p_nsample
    if (dabs((X(i)-X(i-1))-dtmp1) .gt. tol) then
      write(6,*) 'the gap between two elements are not consistent! ', &
                 dabs((X(i)-X(i-1))-dtmp1), i
      stop
    end if 
   end do
 !else
 !  dww = (wf - wi) / p_nsample
 !  dtt = (tf - ti) / p_nsample
 end if
!
 !X = X * scaling
 !XC_in = XC_in * scaling 
! calculate fermi function
 if (lcal_fermi) then
   write(6,*) 'fourier_transform: calculating fermi function'
   write(6,*) 'fourier_transform: k_B * temperature is', T
   write(6,*) 'fourier_transform: sign of fermi function', isgn
   allocate(fermi(p_nsample))
   !dw = (X(p_nsample) - X(1)) / (p_nsample - 1)
   open(unit=7,file='FT_INPUT_FERMI.data')
   rewind(7)
   do i=1, p_nsample
      !ww(i) = dble(i-1) * dw + X(1)
      if ((X(i)*isgn) .ge. 0) then
         fermi(i) = 1.d0-&
                    1.d0/(1.d0+dexp(-isgn/T*X(i)))
      else
         fermi(i) = 1.d0/(1.d0+dexp(isgn/T*X(i)))
      end if
      !write(7,*) X(i), fermi(i)
      XC_in(i) = XC_in(i) * fermi(i) 
      write(7,'(4E17.8,2x)') X(i), fermi(i), dble(XC_in(i)), dimag(XC_in(i))
   end do
   write(6,*) 'fourier_transform: fermi function calculation done'
   close(7)
   deallocate(fermi)
   write(6,*)
   !
   if (isgn .eq. -1) then
     isgn_t = 2
   else if (isgn .eq. 1) then
     isgn_t = 1
   end if 
   write(ch_isgn,'(I2)') isgn_t
   write(ch_ialf,'(I2)') ialf
   write(ch_iorbs,'(I2)') iorbs
   buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
       '-IALF-'//trim(adjustl(ch_ialf))//&
     '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
   title = ('FT_INPUT_FERMI-'//trim(adjustl(buffer)))
   call system('mv ./FT_INPUT_FERMI.data ./'//trim(adjustl(title)))
 else
   goto 250
 end if 
 !
 250 continue
 deallocate(ntmp,num_coloum)
 !
 ! prelude to IFFT
 write(6,*) 'fourier_transform: performing inverse FFT'!,&
 allocate(tt_input(p_nsample), XC_out(p_nsample)) 
 tt_input = 0.d0 
 dww = dabs(X(2)-X(1)) 
 dtmp1 = 2.d0*pi/(dww*p_nsample) 
 call fftfreq(p_nsample,dww,tt_input) 
 call fftshift_1d(tt_input) 
 tt_input = tt_input * 2.d0 * pi 
 dtt = tt_input(2) - tt_input(1) 
 if (dabs(dtt - dtmp1) .gt. tol) then 
   write(6,*) 'fourier_transform: warning! check dt ', dtt, dtmp1 
 end if 
 write(6,*) 'information of frequency and time domain' 
 write(6,'(" from w[1]= ", F12.6," to w[-1]= ", F12.6)')& 
       minval(X), maxval(X) 
 write(6,'(" with dw = ", E14.6)') dww 
 write(6,'(" from t[1]= ", F12.6," to t[-1]= ", F12.6)')& 
       minval(tt_input), maxval(tt_input) 
 write(6,'(" with dt = ", E14.6)') dtt 
 write(6,*) 
 !
 ! Perform a complex to complex IFFT
! istat = DftiCreateDescriptor(My_Desc_Handle, DFTI_DOUBLE,& 
!         DFTI_COMPLEX, 1, p_nsample) 
! istat = DftiSetValue(My_Desc_Handle, DFTI_PLACEMENT,& 
!         DFTI_NOT_INPLACE) 
! istat = DftiCommitDescriptor(My_Desc_Handle) 
! istat = DftiComputeBackward(My_Desc_Handle, XC_in, XC_out) 
! istat = DftiFreeDescriptor(My_Desc_Handle) 
 XC_out = XC_out * dww / pi 
 call fftshift_1d(XC_out) 
 n_middle = p_nsample - p_nsample/2  
 allocate(tt_remain(n_middle),XC_remain(n_middle)) 
 tt_remain = tt_input(int(p_nsample/2)+1:p_nsample) 
 XC_remain = XC_OUT(int(p_nsample/2)+1:p_nsample) 
 do i=1, n_middle 
   if(mod(i,2) .eq. 0) then 
     XC_remain(i) = -1.0d0 * XC_remain(i) 
   end if 
 end do 
 !
 ! interpolate C(t)
 if (linterp) then 
   write(6,*) 'fourier_transform: enter interpolation moudle' 
   if (tf_interp .gt. maxval(tt_remain)) then 
     tf_interp = maxval(tt_remain) 
   end if 
   dtmp1 = tf_interp / dble(num_interp) 
   allocate(tt_interp(num_interp),XC_interp(num_interp)) 
   tt_interp = 0.d0 
   XC_interp = 0.d0 
   do i=1, num_interp 
     tt_interp(i) = dble(i-1) * dtmp1 
   end do 
   write(6,*) 'the number of interpolated points', num_interp 
   write(6,'(" from t[1]= ", F12.6," to t[-1]= ", F12.6)')& 
         minval(tt_interp), maxval(tt_interp) 
   write(6,'(" with dt = ", E14.6)') dtmp1 
   call cspline_interp_complex(tt_remain,XC_remain,tt_interp,XC_interp) 
   open(unit=14,file='IFT.data') 
   rewind(14) 
   write(14,*)  '     time      real      imag      abs' 
   do i=1, num_interp 
      write(14,200) tt_interp(i),& 
                    dreal(XC_interp(i)),& 
                    dimag(XC_interp(i)),& 
                    cdabs(XC_interp(i)) 
   end do 
   call flush(14) 
   close(14) 
 else 
   open(unit=14,file='IFT.data') 
   rewind(14) 
   write(14,*)  '     time      real      imag      abs' 
   do i=1, n_middle 
      write(14,200) tt_remain(i),& 
                    dreal(XC_remain(i)),& 
                    dimag(XC_remain(i)),& 
                    cdabs(XC_remain(i)) 
   end do 
   call flush(14) 
   close(14) 
 end if !
 call system('rm ./FT_INPUT.data ')
 !
 200 format(4(1x,E23.16))
!
end subroutine fourier_transform

