subroutine grid(iacc, nrg, nqg, nrtype, nwtype)
  use sj_dbg
  use omp_lib
  !
  ! purpose: generate grid for numerical integration
  !
  ! input  : 
  !     iacc: the accuracy index for spherical grid
  !     nrg : the number for radial grid
  !     nqg : the number for spherical grid
  !  nrtype : the type of radial grid
  !  nwtype : the type of partion function and partion weight
  !
  implicit none
  integer, parameter :: NSAVAIL = 14
  include '../include/parameter'
  include '../include/common_data'
  include '../include/grid_data'
  !
  integer, intent(in) :: iacc, nrg, nqg, nrtype, nwtype
  !
  integer :: i, j, k, itest, nq, np, ip, iw
  integer :: ia, ir, is, inq, ns, ip0, np0
  integer :: istat
  integer :: nspher(NSAVAIL), nqtype(5,4)
  real*8  :: rrp
  real*8  :: cput1, cput2, cput3
  !
  real*8, allocatable :: angx(:,:), angy(:,:), angz(:,:), angw(:,:)
  real*8, allocatable :: rgq(:), wrgq(:)
  real*8, allocatable :: rabk(:), ramin(:)
  real*8, allocatable :: meshx(:), meshy(:), meshz(:), we(:)

  real*8, allocatable :: drp(:), wbkp(:), bkp(:), uip(:)
  real*8, allocatable :: rmaxshell(:), r1shell(:), r2shell(:), vmaxshell(:),  ratom(:)
  !
  integer, external :: nq_rad
  real*8, external  :: radbeck
  !
  !openmp variables
  real*8, allocatable :: ompmeshx(:,:), ompmeshy(:,:), ompmeshz(:,:), ompwe(:,:)
  integer,allocatable :: iwomp(:)
  ! for test
  !real*8  :: dtmp(maxat), x1, y1, z1, w1, alpha
  !
  data nspher/6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302/
  data nqtype/                                 &
    !            6,  18,  50, 110,  50,          &
  ! in orignial code, it should be 18, instead of 26 in nqtype
    6,  26,  50, 110,  50,          &
    6,  38,  86, 194,  86,          &
    26,  86, 146, 302, 146,          &
    86, 194, 302, 590, 302/
  !
  ! allocate memeory for spherical part
  !
  call cpu_time(cput1)
  write(6,*)
  open(unit=19, file='grid.tmp', form='unformatted')
  rewind 19
  if(nqg>0) then
    itest = 0
    do i = 1, NSAVAIL
      if (nqg==nspher(i)) itest = 1
    end do
    if (itest==0) then
      write(6,*) 'The input spherical grid is not proper!'
      stop
    end if
    if (nqg<86) then
      write(6,*) 'Warning: The spherical grid may be too small!'
    end if
    allocate(angx(nqg,1), angy(nqg,1), angz(nqg,1), angw(nqg,1), STAT = istat)
  else
    if (iacc==1) then
      allocate(angx(110,5), angy(110,5), angz(110,5), angw(110,5), STAT = istat)
    else if(iacc==2) then
      allocate(angx(194,5), angy(194,5), angz(194,5), angw(194,5), STAT = istat)
    else if(iacc==3) then
      allocate(angx(302,5), angy(302,5), angz(302,5), angw(302,5), STAT = istat)
    else if(iacc==4) then
      allocate(angx(590,5), angy(590,5), angz(590,5), angw(590,5), STAT = istat)
    else
      write(6,*) 'iacc should only be 1, 2, 3,and 4!'
      stop
    end if
  end if
  !
  allocate(rabk(Nat), ramin(Nat), STAT = istat)
  !
  allocate(ompmeshx(MAXPT,Nat), ompmeshy(MAXPT,Nat), ompmeshz(MAXPT,Nat), ompwe(MAXPT,Nat),iwomp(Nat), STAT = istat)
  if(istat/=0)then
    print*,"allocate arrays for openmp failed."
    print*,"try to use a smaller MAXPT."
    stop
  endif


  if (Nat==1) then
    ramin(1) = 1.d5
  else
    do i = 1, Nat
      rabk(i) = radbeck(Atn(i))
      if (i/=Nat) then
        ramin(i) = Rr(i,i+1)
      else
        ramin(i) = Rr(i,i-1)
      end if
      do j = 1, Nat
        if(i/=j) then
          ramin(i) = min(ramin(i), Rr(i,j))
        end if
      end do
      !     write(6,*) 'i,ramin', i, ramin(i)
      ramin(i) = 0.18d0*ramin(i)
    end do
  end if
  !
  ! generate spherical grid
  !
  call SJ_TIMER_START(5)
  if (nqg>0) then
    call sphgrid(nqg, angx, angy, angz, angw)
  else
    do i = 1, 5
      nq = nqtype(i,iacc)
      call sphgrid(nq, angx(:,i), angy(:,i), angz(:,i), angw(:,i))
    end do
  end if
  !
  ! generate grid and weight for numerical integration
  !

  Nptot = 0
  if (nqg/=0) ns = nqg
  write(6,*) 'Generate radial grids (openmp version)'
  write(6,*) 'nqg=', nqg
  call SJ_TIMER_STOP(5,"-----grid: generate sphgrid                          ")
  call SJ_TIMER_START(5)
  !!$omp shared(ompmeshx,ompmeshy,ompmeshz,ompwe,iwomp)       &
  !$omp parallel default(shared)                                      & 
  !$omp private(ia,np,ir,rrp,np0,is,inq,ns,iw,meshx,meshy,meshz,we)   &
  !$omp private(rgq,wrgq)                                             &
  !$omp private(drp, bkp, wbkp, uip)
  allocate(rgq(nrg), wrgq(nrg), STAT = istat)
  allocate(meshx(MAXPT), meshy(MAXPT), meshz(MAXPT), we(MAXPT), STAT = istat)
  allocate(drp(Nsblk*Nat), wbkp(Nsblk), bkp(Nsblk*Nat), uip(Nsblk), STAT = istat)
  if(istat/=0)then
    print*,"grid: allocate arrays failed"
    stop
  endif
  !$ if(omp_get_thread_num()==0) write(80,*)"grid num threads:",omp_get_num_threads()
  !$omp do schedule(guided) reduction(+:Nptot)
  do ia = 1, Nat
    !
    ! generate radial grid 
    !
    np = 0
    !radgrid output: rgq,wrgq
    call radgrid(nrg, Atn(ia), rgq, wrgq, nrtype)
    do ir = 1, nrg
      rrp = rgq(ir)
      np0 = np + 1
      if (nqg/=0) then
        do is = 1, nqg
          np = np + 1
          !           if (ir==nrg.and.is==1) npt = np
          if (np>MAXPT) then
            write(6,*) 'The number per atom is too large!'
            stop
          end if
          meshx(np) = rrp*angx(is,1) + Xyz(1,ia)
          meshy(np) = rrp*angy(is,1) + Xyz(2,ia)
          meshz(np) = rrp*angz(is,1) + Xyz(3,ia)
          we(np) = wrgq(ir)*angw(is,1)
        end do
      else
        inq = nq_rad(Atn(ia), rrp)
        ns = nqtype(inq, iacc)
        !        write(6,*) 'ir, rrp, ns', ir, rrp, ns
        do is = 1, ns
          np = np + 1
          if (np>MAXPT) then
            write(6,*) 'The number per atom is too large!'
            stop
          end if
          meshx(np) = rrp*angx(is,inq) + Xyz(1,ia)
          meshy(np) = rrp*angy(is,inq) + Xyz(2,ia)
          meshz(np) = rrp*angz(is,inq) + Xyz(3,ia)
          we(np) = wrgq(ir)*angw(is,inq)
        end do
      end if
      if (rrp>ramin(ia)) then
        call weight(ns, Nat, nwtype, ia, rrp, meshx(np0), meshy(np0), meshz(np0),    &
          we(np0), rabk, Xyz, Rr, drp, bkp, wbkp, uip)
      end if
    end do
    !
    !  write(6,'(1x, A, i3, i8)') 'Total number of points for atom  ', ia, np
    !  write(6,*) 'meshx, y, z', meshx(1), meshy(1), meshz(1)
    !  write(6,*) 'we', we(1)
    !  write(6,*) 'meshx, y, z', meshx(npt), meshy(npt), meshz(npt)
    !  write(6,*) 'we', we(npt)
    iw = 0
    do ip = 1, np  
      if (we(ip)>1.d-20) then
        iw = iw + 1
        ompmeshx(iw,ia) = meshx(ip)
        ompmeshy(iw,ia) = meshy(ip)
        ompmeshz(iw,ia) = meshz(ip)
        ompwe(iw,ia)    =    we(ip)
        !meshx(iw) = meshx(ip)
        !meshy(iw) = meshy(ip)
        !meshz(iw) = meshz(ip)
        !   we(iw) =    we(ip)
      end if
    end do
    !write(6,'(1x, A, i11)') 'Total number of points after sort', iw
    Nptot = Nptot + iw
    iwomp(ia)=iw
    !Sequential version output:
    !do i = 1, iw
    !   write(19) meshx(i), meshy(i), meshz(i), we(i)
    !end do
  end do
  !$omp end do
  deallocate(rgq, wrgq, STAT = istat)
  deallocate(drp, wbkp, bkp, uip, STAT = istat)
  deallocate(meshx, meshy, meshz, we, STAT = istat)
  !$omp end parallel
  !Parallel version output:
  do ia=1,Nat
    do i = 1, iwomp(ia)
      write(19) ompmeshx(i,ia), ompmeshy(i,ia), ompmeshz(i,ia), ompwe(i,ia)
      !write(6,*) ompmeshx(i,ia), ompmeshy(i,ia), ompmeshz(i,ia), ompwe(i,ia)
    end do
    write(6,'(1x, A, i3, i11)') 'Total number of points after sort for atom',ia, iwomp(ia)
  enddo
  close(19)
  call SJ_TIMER_STOP(5,"-----grid: generate radial grid                 ")
  ! test here
  !open(unit=19, file='grid.tmp', form='unformatted')
  !rewind 19
  !dtmp = 0.d0
  !alpha = 1.9d0
  !do i = 1, Nptot
  !   read(19) x1, y1, z1, w1
  !   do j = 1, Nat
  !   rrp = (x1-Xyz(1,j))**2 + (y1-Xyz(2,j))**2 + (z1-Xyz(3,j))**2
  !   dtmp(j) = dtmp(j) + ((z1-Xyz(3,j))**2)*dexp(-rrp*alpha)*w1
  !   end do
  !end do
  !do i = 1, Nat
  !write(6,*) 'dtmp', dtmp(i)
  !end do
  !write(6,*) pi*dsqrt(pi)/(2.d0*alpha**2*dsqrt(alpha))
  !
  write(6,*) 'Total number of points', Nptot
  deallocate(rabk, STAT = istat)
  deallocate(angx, angy, angz, angw, STAT = istat)

  !
  ! sort the points to batches and calculate the neighbor
  ! atoms as well as basis tails
  !
  allocate(meshx(MAXPT), meshy(MAXPT), meshz(MAXPT), we(MAXPT), STAT = istat)
  allocate(rmaxshell(Nshell), r1shell(Nshell), r2shell(Nshell), vmaxshell(Nshell),   &
    ratom(Nat), STAT = istat)
  call shellpeak(rmaxshell, r1shell, r2shell, vmaxshell, ratom)
  !
  call cpu_time(cput2)
  write(6,*) 'cpu time for generate grid', cput2-cput1
  call sortgrid(meshx, meshy, meshz, we, vmaxshell, r1shell, r2shell)
  call cpu_time(cput3)
  write(6,*) 'cpu time for sorting grid', cput3-cput2
  !
  deallocate(rmaxshell, r1shell, r2shell, vmaxshell, ratom, STAT = istat)
  deallocate(meshx, meshy, meshz, we, STAT = istat)
  !
end subroutine grid
