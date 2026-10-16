subroutine vasp2heom(ialf, iorbs, T_tmp)
  implicit none
  integer, intent(in) :: ialf, iorbs
  real*8, intent(in)  :: T_tmp
  integer  :: i, n, m, ntmp1, ntmp2, istat
  integer  :: sys_nband, env_nband, num_ww, chosen_spin
  real*8   :: dtmp1, dtmp2, dtmp3
  real*8   :: sys_fermi, sub_fermi, fermi_distance
  real*8   :: ww_i, ww_f, ww_step, width, heom_width
  integer,allocatable :: chosen_env_bands(:)
  real*8, allocatable :: ww(:), hyb(:)
  real*8, parameter   :: pi = 3.1415926d0
  !real*8, allocatable :: sys_eig(:,:), sys_occ(:,:)
  real*8, allocatable :: env_eig(:,:), env_occ(:,:)
  real*8, allocatable :: hopping(:,:), sum_hopping(:)
  character(len=128)  :: buffer, ch_tmp, title
  character(len=2)    :: ch_ialf, ch_isgn, ch_iorbs
  logical :: lfile_exist1, lfile_exist2
  !
  integer :: j, npade, num_tt
  real*8  :: T, sigma(2), crit
  real*8  :: tt_f, tt_step, tt
  real*8, parameter :: hbar = 0.658211928d0
  !real*8, parameter :: hbar = 1.d0
  complex*16 :: ctmp1, ctmp2
  complex*16, parameter :: cunit = (1.0d0, 0.0d0), &
                           eye = (0.0d0, 1.0d0), &
                           czero = (0.0d0, 0.0d0)
  complex*16, allocatable :: band_center(:), coupling(:)
  complex*16, allocatable :: zpeigv(:), zpcoef(:), zppole(:,:)
  complex*16, allocatable :: eta_tmp(:,:,:), gama_tmp(:,:,:)
  complex*16, allocatable :: ct_tot(:,:), ct(:,:,:)
  data sigma /1, -1/
  namelist / vasp / npade, tt_f, num_tt, crit,&
                    ww_i, ww_f, num_ww, width, chosen_spin,&
                    fermi_distance
                    
  if (ialf .lt. 1) then
    return
  end if
  ww_i = -10
  ww_f =  10  
  num_ww = 10001
  allocate(ww(num_ww), hyb(num_ww))
  ww = 0.d0
  hyb = 0.d0
  width = 0.03d0
  chosen_spin = 1
  heom_width = 0.3d0
  !
  npade = 650
  T = T_tmp
  tt_f = 200.d0
  num_tt = 5001
  crit = 1.d-8
  fermi_distance = 10.d0
  open(unit=5,status='unknown')
  rewind(5)
  read(5, vasp, end=93)
  93 continue
  close(5)
  if (mod(num_tt,2) .eq. 0) then
    write(6,*) 'warning: num_tt should be odd'
    num_tt = num_tt + 1
  end if
  tt_step = tt_f / (num_tt - 1)
  allocate(zpeigv(npade), zpcoef(npade), zppole(npade,2))
  zpeigv = czero
  zpcoef = czero
  zppole = czero
  call psd(npade, 2, zpeigv, zpcoef, zppole)
  !
  !write(6,*) '*********************************************************'
  !write(6,*) '          generate hybdridization funtion     '
  !write(6,*) '*********************************************************'
  write(6,'("Width of Lorentzian functions to approximate delta functions",E13.5)') width
  write(6,'("Width of Lorentzian function in HEOM calculations", E13.5)') heom_width
  write(6,'("Number of pade poles", I5)') npade
  write(6,'("Temperature of slab", E13.5)') T
  write(6,'("Final time of C(t)", E13.5)') tt_f
  write(6,'("Number of time stpe ", I7)') num_tt
  write(6,'("Time stpe of C(t)", E13.5)') tt_step
  write(6,'("Crition for minimal hopping strength", E13.5)') crit
  write(6,'("Single occupied spin is", I4)') chosen_spin
  write(6,'("Distance between fermi level", I4)') fermi_distance
  write(6,*)
  !
  sys_fermi = 0.d0
  sub_fermi = 0.d0
  inquire(file='fermi_sys', exist=lfile_exist2)
  if (lfile_exist2) then
     title = 'fermi_sys'
     open(unit=63, file=title)
     rewind(63)
     read(63,'(A)') buffer
     read(63,'(A)') buffer
     read(buffer(1:23),'(A10,F12.7)') ch_tmp, sys_fermi
     close(63)
     write(6,*) 'Fermi level of system: ', sys_fermi
  end if
  !
  inquire(file='fermi_sub', exist=lfile_exist1)
  if (lfile_exist1) then
      title = 'fermi_sub'
      open(unit=64, file=title)
      rewind(64)
      read(64,'(A)') buffer
      read(64,'(A)') buffer
      read(buffer(1:23),'(A10,F12.7)') ch_tmp, sub_fermi
      close(64)
      write(6,*) 'Fermi level of slab: ', sub_fermi
  end if
  !
  if (chosen_spin .eq. 1) then
    inquire(file='to_be_read_env_bands_sp1', exist=lfile_exist1)
    inquire(file='to_be_read_sys_bands_sp1', exist=lfile_exist2)
    if (lfile_exist1 .and. lfile_exist2) then
      open(unit=10,file='to_be_read_sys_bands_sp1')
      rewind(10)
      sys_nband = -1
      istat = 0
      write(6,*) 'read to_be_read_sys_bands_sp1 file'
      do while(.not. is_iostat_end(istat))
        read(10,'(a)', iostat=istat) buffer
        sys_nband = sys_nband + 1
      end do
      close(10)
      !
      write(6,*) 'number of system bands', sys_nband
      !
      open(unit=11,file='to_be_read_env_bands_sp1')
      rewind(11)
      env_nband = -1
      istat = 0
      write(6,*) 'read to_be_read_env_bands_sp1 file'
      do while(.not. is_iostat_end(istat))
        read(11,'(a)', iostat=istat) buffer
        env_nband = env_nband + 1
      end do
      !
      write(6,*) 'number of environment bands', env_nband
      !
      allocate(chosen_env_bands(env_nband))
      allocate(env_eig(env_nband,2), env_occ(env_nband,2))
      rewind(11)
      do i=1, env_nband
        read(11,100) chosen_env_bands(i),&
                     env_occ(i,1),env_eig(i,1),&
                     env_occ(i,2),env_eig(i,2)
      end do
      close(11)
      !
      ww_i = minval(env_eig(:,1)) - 15.d0
      ww_f = maxval(env_eig(:,1)) + 15.d0
      open(unit=5,status='unknown')
      rewind(5)
      read(5, vasp, end=113)
      113 continue
      close(5)
      ww_step = (ww_f - ww_i) / num_ww
      do m=1, num_ww
         ww(m) = (m-1) * ww_step + ww_i
      end do
      !
      allocate(hopping(sys_nband,env_nband))
      allocate(sum_hopping(env_nband))
      hopping = 0.d0
      sum_hopping = 0.d0
      open(unit=12, file='spin-1-coupling-details')
      rewind(12)
      do m=1, 6
         read(12,'(a)') buffer
      end do
      !
      do m=1, env_nband 
        do n=1, sys_nband
          read(12,101) ntmp1, ntmp2, &
                       dtmp1, dtmp2, dtmp3, &
                       hopping(n,m)
          sum_hopping(m) = sum_hopping(m) + &
                           hopping(n,m)**2
        end do  
      end do
      close(12)
      !
      open(unit=13,file='hyb-sp1')
      rewind(13)
      do m=1, num_ww
        do n=1, env_nband
          if (dabs(env_eig(n,1)-sub_fermi) .gt. fermi_distance) then
            cycle
          end if
          hyb(m) = hyb(m) + sum_hopping(n) * &
          width / ((ww(m)-env_eig(n,1))**2+(width*width))
        end do
        write(13,102) ww(m) - sub_fermi, hyb(m)
      end do
      close(13)
      !
      write(ch_isgn,'(I2)') 1
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('FT_INPUT-'//trim(adjustl(buffer)))
      call system('cp ./hyb-sp1 '//trim(adjustl(title)))
      !
      write(ch_isgn,'(I2)') 2
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('FT_INPUT-'//trim(adjustl(buffer)))
      call system('cp ./hyb-sp1 '//trim(adjustl(title)))
      !
      write(6,103) (1.0d0/heom_width)*sum(sum_hopping)!/&
                  !(atan(ww_f/heom_width)-atan(ww_i/heom_width))
      !
      ! Pade decomposition
      allocate(band_center(env_nband), coupling(env_nband))
      allocate(eta_tmp(npade+1,env_nband,2), gama_tmp(npade+1,env_nband,2))
      allocate(ct_tot(num_tt,2), ct(num_tt,env_nband,2))
      band_center = env_eig(:,1) - sub_fermi
      coupling = sum_hopping
      eta_tmp = czero
      gama_tmp = czero
      ct_tot = czero
      ct = czero
      !
      do i=1, 2
        do m=1, env_nband
          if(cdabs(coupling(m)) .lt. crit) then
            cycle
          end if
          if (dabs(env_eig(m,1)-sub_fermi) .gt. fermi_distance) then
            cycle
          end if
          ! generate spectral results
          gama_tmp(1,m,i) = (width - eye * sigma(i) * band_center(m)) / hbar
          ctmp1 = eye * gama_tmp(1,m,i) / T * hbar
          ctmp2 = dcmplx(5.d-1, 0.d0)
          do n=1, npade
            ctmp2 = ctmp2 + zpcoef(n) / (ctmp1 + zppole(n,i)) + &
                       zpcoef(n) / (ctmp1 - zppole(n,i))
          end do
          eta_tmp(1,m,i) = 0.5d0 * coupling(m) * width * ctmp2 / hbar**2
          do n=2, npade + 1
            gama_tmp(n,m,i) = ( -1.d0 * eye * sigma(i) * T * zppole(n-1,i) ) / hbar
            eta_tmp(n,m,i) = 2.d0 * eye * zpcoef(n-1) * coupling(m) * 0.5d0 * width**2 * T &
                          / ( (zppole(n-1,i) * T - band_center(m))**2 + width**2 ) / hbar**2
          end do
          ! calculate C(t)
          do j=1, num_tt
            tt = (j-1) * tt_step
            do n=1, npade+1
              ct(j,m,i) = ct(j,m,i) + &
               eta_tmp(n,m,i) * cdexp(-1.d0 * gama_tmp(n,m,i) * tt)
            end do
            ct_tot(j,i) = ct_tot(j,i) + ct(j,m,i)
          end do
          write(6,'("calculating for sigma = ",I2, " and No. ", I5," band with eigenvalue ", E12.5, " and hopping ",E12.5)')&
          i, chosen_env_bands(m), env_eig(m,1), dreal(coupling(m)) 
        end do
      end do
      !
      open(unit=21, file='ct-sp1.data')
      rewind(21)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(21,'(5E21.12)') tt, & 
                              dreal(ct_tot(j,1)), dimag(ct_tot(j,1)),&
                              dreal(ct_tot(j,2)), dimag(ct_tot(j,2))
      end do
      close(21)
      !
      write(ch_isgn,'(I2)') 1
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('IFT-'//trim(adjustl(buffer)))
      open(unit=51, file=title)
      rewind(51)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(51,'(3(E23.16,1x))') tt*hbar, & 
                              dreal(ct_tot(j,1))/(2.d0*hbar**2),&
                              dimag(ct_tot(j,1))/(2.d0*hbar**2)
      end do
      close(51)
      !
      write(ch_isgn,'(I2)') 2
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('IFT-'//trim(adjustl(buffer)))
      open(unit=52, file=title)
      rewind(52)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(52,'(3(E23.16,1x))') tt, & 
                              dreal(ct_tot(j,2)),&
                              dimag(ct_tot(j,2))
      end do
      close(52)
      !
      deallocate(sum_hopping,hopping,chosen_env_bands)
      deallocate(env_occ,env_eig)
      deallocate(band_center,coupling,eta_tmp,gama_tmp,ct_tot,ct)
    else
      write(6,*) 'Do not find to_be_read_sys_bands_sp1 '&
                 'or to_be_read_env_bands_sp1'
    end if
  end if
  !
  ww = 0.d0
  hyb = 0.d0
  !
  zpeigv = czero
  zpcoef = czero
  zppole = czero
  call psd(npade, 2, zpeigv, zpcoef, zppole)
  !
  if (chosen_spin .eq. 2) then
    inquire(file='to_be_read_env_bands_sp2', exist=lfile_exist1)
    inquire(file='to_be_read_sys_bands_sp2', exist=lfile_exist2)
    if (lfile_exist1 .and. lfile_exist2) then
      open(unit=14,file='to_be_read_sys_bands_sp2')
      rewind(14)
      sys_nband = -1
      istat = 0
      write(6,*) 'read to_be_read_sys_bands_sp2 file'
      do while(.not. is_iostat_end(istat))
        read(14,'(a)', iostat=istat) buffer
        sys_nband = sys_nband + 1
      end do
      close(14)
      !
      write(6,*) 'number of system bands', sys_nband
      !
      open(unit=15,file='to_be_read_env_bands_sp2')
      rewind(15)
      env_nband = -1
      istat = 0
      write(6,*) 'read to_be_read_env_bands_sp2 file'
      do while(.not. is_iostat_end(istat))
        read(15,'(a)', iostat=istat) buffer
        env_nband = env_nband + 1
      end do
      !
      write(6,*) 'number of environment bands', env_nband
      !
      allocate(chosen_env_bands(env_nband))
      allocate(env_eig(env_nband,2), env_occ(env_nband,2))
      rewind(15)
      do i=1, env_nband
        read(15,100) chosen_env_bands(i),&
                     env_occ(i,1),env_eig(i,1),&
                     env_occ(i,2),env_eig(i,2)
      end do
      close(15)
      !
      ww_i = minval(env_eig(:,2)) - 15.d0
      ww_f = maxval(env_eig(:,2)) + 15.d0
      open(unit=5,status='unknown')
      rewind(5)
      read(5, vasp, end=123)
      123 continue
      close(5)
      ww_step = (ww_f - ww_i) / num_ww
      do m=1, num_ww
         ww(m) = (m-1) * ww_step + ww_i
      end do
      !
      allocate(hopping(sys_nband,env_nband))
      allocate(sum_hopping(env_nband))
      hopping = 0.d0
      sum_hopping = 0.d0
      open(unit=16, file='spin-2-coupling-details')
      rewind(16)
      do m=1, 6
         read(16,'(a)') buffer
      end do
      do m=1, env_nband 
        do n=1, sys_nband
          read(16,101) ntmp1, ntmp2, &
                       dtmp1, dtmp2, dtmp3, &
                       hopping(n,m)
          sum_hopping(m) = sum_hopping(m) + &
                           hopping(n,m)**2
        end do  
      end do
      close(16)
      !
      open(unit=17,file='hyb-sp2')
      rewind(17)
      do m=1, num_ww
        do n=1, env_nband
          if (dabs(env_eig(n,2)-sub_fermi) .gt. fermi_distance) then
            cycle
          end if
          hyb(m) = hyb(m) + sum_hopping(n) * &
          width / ((ww(m)-env_eig(n,2))**2+(width*width))
        end do
        write(17,102) ww(m) - sub_fermi, hyb(m)
      end do
      close(17)
      !
      write(ch_isgn,'(I2)') 1
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('FT_INPUT-'//trim(adjustl(buffer)))
      call system('cp ./hyb-sp2 '//trim(adjustl(title)))
      !
      write(ch_isgn,'(I2)') 2
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('FT_INPUT-'//trim(adjustl(buffer)))
      call system('cp ./hyb-sp2 '//trim(adjustl(title)))
      !
      write(6,104) (1.0d0/heom_width)*sum(sum_hopping)!/&
                  !(atan(ww_f/heom_width)-atan(ww_i/heom_width))
      !
      ! Pade decomposition
      allocate(band_center(env_nband), coupling(env_nband))
      allocate(eta_tmp(npade+1,env_nband,2), gama_tmp(npade+1,env_nband,2))
      allocate(ct_tot(num_tt,2), ct(num_tt,env_nband,2))
      band_center = env_eig(:,1) - sub_fermi
      coupling = sum_hopping
      eta_tmp = czero
      gama_tmp = czero
      ct_tot = czero
      ct = czero
      !
      do i=1, 2
        do m=1, env_nband
          if(cdabs(coupling(m)) .lt. crit) then
            cycle
          end if
          if (dabs(env_eig(m,2)-sub_fermi) .gt. fermi_distance) then
            cycle
          end if
          gama_tmp(1,m,i) = (width - eye * sigma(i) * band_center(m)) / hbar
          ctmp1 = eye * gama_tmp(1,m,i) / T * hbar
          ctmp2 = dcmplx(5.d-1, 0.d0)
          do n=1, npade
            ctmp2 = ctmp2 + zpcoef(n) / (ctmp1 + zppole(n,i)) + &
                       zpcoef(n) / (ctmp1 - zppole(n,i))
          end do
          eta_tmp(1,m,i) = 0.5d0 * coupling(m) * width * ctmp2 / hbar**2
          do n=2, npade + 1
            gama_tmp(n,m,i) = ( -1.d0 * eye * sigma(i) * T * zppole(n-1,i) ) / hbar
            eta_tmp(n,m,i) = 2.d0 * eye * zpcoef(n-1) * coupling(m) * 0.5d0 * width**2 * T &
                          / ( (zppole(n-1,i) * T - band_center(m))**2 + width**2 ) / hbar**2
          end do
          ! calculate C(t)
          do j=1, num_tt
            tt = (j-1) * tt_step
            do n=1, npade+1
              ct(j,m,i) = ct(j,m,i) + &
               eta_tmp(n,m,i) * cdexp(-1.d0 * gama_tmp(n,m,i) * tt)
            end do
            ct_tot(j,i) = ct_tot(j,i) + ct(j,m,i)
          end do
          write(6,'("calculating for sigma = ",I2, " and No. ", I5," band with eigenvalue ", E12.5, " and hopping ",E12.5)')&
          i, chosen_env_bands(m), env_eig(m,2), dreal(coupling(m)) 
        end do
      end do
      open(unit=22, file='ct-sp2.data')
      rewind(22)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(22,'(5E21.12)') tt, & 
                              dreal(ct_tot(j,1)), dimag(ct_tot(j,1)),&
                              dreal(ct_tot(j,2)), dimag(ct_tot(j,2))
      end do
      close(22)
      !
      write(ch_isgn,'(I2)') 1
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('IFT-'//trim(adjustl(buffer)))
      open(unit=53, file=title)
      rewind(53)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(53,'(3(E23.16,1x))') tt, & 
                              dreal(ct_tot(j,1)),&
                              dimag(ct_tot(j,1))
      end do
      close(53)
      !
      write(ch_isgn,'(I2)') 2
      write(ch_ialf,'(I2)') ialf
      write(ch_iorbs,'(I2)') iorbs
      buffer = ('ISGN-'//trim(adjustl(ch_isgn))//&
                '-IALF-'//trim(adjustl(ch_ialf))//&
                '-IORBS-'//trim(adjustl(ch_iorbs))//'.data')
      title = ('IFT-'//trim(adjustl(buffer)))
      open(unit=54, file=title)
      rewind(54)
      do j=1, num_tt
        tt = (j-1) * tt_step
        write(54,'(3(E23.16,1x))') tt, & 
                              dreal(ct_tot(j,2)),&
                              dimag(ct_tot(j,2))
      end do
      close(54)
      !
      deallocate(sum_hopping,hopping,chosen_env_bands)
      deallocate(env_occ,env_eig)
    else
      write(6,*) 'Do not find to_be_read_sys_bands_sp2 '&
                 'or to_be_read_env_bands_sp2'
    end if
  end if
  !
  100 format('No. of band',1x,I5,2x,&
             ',occ. spin-1',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6,2x,&
             ',occ. spin-2',2x,F9.6,2x,&
             ',eigenergy',2x,F11.6)
  101 format(I10, 2X, I10, 4E18.10)
  102 format(2(E18.10,2x))
  103 format("For Lorentzian type spectral function,",&
             " effective spin-1 coupling is", 2x, E18.10)
  104 format("For Lorentzian type spectral function,",&
             " effective spin-2 coupling is", 2x, E18.10)
  !
end subroutine vasp2heom
