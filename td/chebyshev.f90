subroutine chebyshev
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                 :: ni, nj, nk, nstep, istat
integer                 :: iorbs, ialf, ispin, imats
integer*8               :: lni, lnj
logical                 :: fexist
real*8, external        :: bessj
real*8                  :: tt, dtmp1, dtmp2, dtmp3, dcoefb, ddt, tout
real*8                  :: cpu1, cpu2, cpu3, cpu4
real*8                  :: occ, occu, occd, jleft, jright
real*8                  :: jleftu, jleftd, jrightu, jrightd
real*8                  :: eshift(maxalf,maxspin)
real*8,     allocatable :: bessfunc(:,:), tcheby(:)
complex*16, allocatable :: sigmau(:), sigmad(:)
complex*16, allocatable :: zdm0(:,:,:), zdm1(:,:,:,:,:,:)
complex*16              :: ctmp1
!
write(6,*)
write(6,*)' entering chebyshev '
call flush(6)
call cpu_time(cpu1)
!
tt = dpico
!
if (lresume .and. icont .eq. 1) then
  tt = tresume 
end if
!
do ialf=1,nalf
   do ispin=1,nspin
      call getengyshift(tt, ialf, ispin, eshift(ialf,ispin))
   end do
end do
eshift = eshift * hbar
!
allocate(bessfunc(ncheby,jcheby), tcheby(jcheby), STAT=istat)
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
allocate(zdm0(nrho,nrho,jcheby),                                             &
         zdm1(nrho,nrho,norbs,nalf,nspin,jcheby), STAT=istat)
!
write(6,*)
write(6,*)' temporary memory allocate for current and occupation '
write(6,*)' memory size = ', dble(nrho**2*jcheby*(norbs*nalf*nspin+1)*16) /  &
                             dble(1024**2), ' MB'
!
ddt   = dt / dble(jcheby)
dtmp1 = 1.d0 / cheby_delta
dtmp2 = 2.d0 * dtmp1
nstep = 0
call cpu_time(cpu3)
!
do ni=1,jcheby
  tcheby(ni) = dble(ni) * ddt
end do
!
do lni=1,ncheby
  do ni=1,jcheby
    bessfunc(lni, ni) = bessj(lni-1, tcheby(ni) * cheby_delta)
  end do
end do
!
100 continue
!
! Calculate occupation and current | output % by rho at tt
! 
call calcocc(occu, occd, sigmau, sigmad)
call calcj(jleftu, jrightu, jleftd, jrightd)
occ   = occu + occd
jleft = jleftu + jleftd
if (nalf .gt. 1) jright = jrightu + jrightd
!
call cpu_time(cpu4)
if (nspin .eq. 1) then
   if (nalf .eq. 1) then
      write(6,518)tt, jleftu, occu, cpu4-cpu3
      if (nstep .eq. 0) then
        write(13,518)tt, eshift(1,1), jleft, occu
        write(16,518)tt, eshift(1,1), occu
      end if
   else if (nalf .eq. 2) then
      write(6,518)tt, jleftu, jrightu, jleftu+jrightu, occu, cpu4-cpu3
      if (nstep .eq. 0) then
        write(13,518)tt, eshift(1,1), eshift(2,1), jleft, jright, jleft+jright, occu
        write(16,518)tt, eshift(1,1), eshift(2,1), occu
      end if
   else ! multi-leads
      write(6,518)tt, (jleads(ialf), ialf=1,nalf), occu, cpu4-cpu3
      if (nstep .eq. 0) then
         write(13,518)tt, (eshift(ialf,1), ialf=1,nalf), (jleads(ialf), ialf=1,nalf), occu
         write(16,518)tt, (eshift(ialf,1), ialf=1,nalf), occu
      end if
   end if
else 
   if (nalf .eq. 1) then
      write(6,518)tt, jleftu, jleftd, jleftu+jleftd, occu, occd, cpu4-cpu3
      if (nstep .eq. 0) then
        write(13,518)tt, eshift(1,1), jleftu, jleftd, jleftu+jleftd, jleftu-jleftd, occ
        write(16,518)tt, eshift(1,1), occu, occd, occ
      end if
   else if (nalf .eq. 2) then
      write(6,518)tt, jleft, jright, occu, occd, cpu4-cpu3
      if (nstep .eq. 0) then
        write(13,518)tt, eshift(1,1), eshift(2,1), jleft, jright,                      &
                     jleftu, jleftd, jrightu, jrightd, occu, occd
        write(16,518)tt, eshift(1,1), eshift(2,1), occu, occd, occ
      end if
   else ! multi-leads
      write(6,518)tt, (jleads(ialf), ialf=1,nalf), occu, occd, cpu4-cpu3
      if (nstep .eq. 0) then
         write(13,518)tt, (eshift(ialf,1), ialf=1,nalf), (jleads(ialf), ialf=1,nalf),  &
                     (jt(ialf,1), ialf=1,nalf), (jt(ialf,2), ialf=1,nalf), occu, occd
         write(16,518)tt, (eshift(ialf,1), ialf=1,nalf), occu, occd, occ
      end if
   end if
end if
call cpu_time(cpu3)
call flush(6)
!
if (tt .ge. tmax) goto 200
!
lni  = 1
rho0 = rho
rho  = bessfunc(1,jcheby) * rho0  
!
do ni=1,jcheby
  zdm0(1:nrho,1:nrho,ni) = rho0(1:nrho,1:nrho,1) * bessfunc(1,ni)
  do ialf=1,nalf
    do ispin=1,nspin
      do iorbs=1,norbs
        cmtmp1 = czero
        do imats=1,ncor
          if (lscale) then
             dtmp3 = dbsqrt(iorbs,ispin,imats,ialf,1)
          else
             dtmp3 = 1.d0
          end if
          call look4drawer(1, ialf, iorbs, ispin, imats, nk)
          cmtmp1(1:nrho, 1:nrho) = cmtmp1(1:nrho, 1:nrho) +               &
                                   rho0(1:nrho, 1:nrho, nfirst(2)-1+nk) * dtmp3
        end do
        zdm1(1:nrho,1:nrho,iorbs,ialf,ispin,ni) = cmtmp1(1:nrho,1:nrho) * &
                                                  bessfunc(1,ni)
      end do
    end do
  end do
end do
!
! P_n-1 : rho0 
! P_n-2 : rhotmp
!
if (nstep .eq. 0) then
  write(6,*)
  write(6,*)' convergency test for the first time-step ', dt
end if
do 
  dcoefb = bessfunc(lni+1, jcheby) * 2.d0
  if (nstep .eq. 0) then
    write(6,*)lni, dble(rho0(1,1,1)), dble(rho0(1,1,1)) * dcoefb
    call flush(6)
  end if
  if (lni .eq. 1) then
    rhotmp = rho0                          ! P0
    call calcderiv_omp(tt)
    rho0   = rhorhs * dtmp1                ! P1
  else if (lni .ge. 2) then
    rhorhs = rhotmp
    rhotmp = rho0                          ! P1 (P_n-1)
    rho0   = rhorhs                        ! P0 (P_n-2)
    call calcderiv_omp(tt)
    rho0   = rho0 + rhorhs * dtmp2         ! P2 (P_n)
  end if
  rho = rho + rho0 * dcoefb
  lni = lni + 1
!
  do ni=1,jcheby
    zdm0(1:nrho,1:nrho,ni) = zdm0(1:nrho,1:nrho,ni) + rho0(1:nrho,1:nrho,1) * &
                             bessfunc(lni,ni) * 2.d0
    do ialf=1,nalf
      do ispin=1,nspin
        do iorbs=1,norbs
          cmtmp1 = czero
          do imats=1,ncor
            call look4drawer(1, ialf, iorbs, ispin, imats, nk)
            cmtmp1(1:nrho, 1:nrho) = cmtmp1(1:nrho, 1:nrho) +                                      &
                                     rho0(1:nrho, 1:nrho, nfirst(2)-1+nk)
          end do
          zdm1(1:nrho,1:nrho,iorbs,ialf,ispin,ni) = zdm1(1:nrho,1:nrho,iorbs,ialf,ispin,ni) +      &
                                                    cmtmp1(1:nrho,1:nrho) * bessfunc(lni,ni) * 2.d0
        end do
      end do
    end do
  end do
!
  if (lni .ge. ncheby) exit
end do
!
! Evaluate occupation and current (fine meshgrids in time) | output to data
!
do ni=1,jcheby
  tout = tt + tcheby(ni)
  call chebyocc(zdm0(1,1,ni), occu, occd, sigmau, sigmad)
  call chebyj(zdm1(1,1,1,1,1,ni), jleftu, jrightu, jleftd, jrightd)
  occ = occu + occd
  jleft = jleftu + jleftd
  if (nalf .ne. 1) jright = jrightu + jrightd
!
  if (nspin .eq. 1) then
     if (nalf .eq. 1) then
       write(13,518)tout, eshift(1,1), jleft, occu
       write(16,518)tout, eshift(1,1), occu
     else if (nalf .eq. 2) then
       write(13,518)tout, eshift(1,1), eshift(2,1), jleft, jright, jleft+jright, occu
       write(16,518)tout, eshift(1,1), eshift(2,1), occu
     else ! multi-leads
       write(13,518)tout, (eshift(ialf,1), ialf=1,nalf), (jleads(ialf), ialf=1,nalf), occu
       write(16,518)tout, (eshift(ialf,1), ialf=1,nalf), occu
     end if
  else 
     if (nalf .eq. 1) then
       write(13,518)tout, eshift(1,1), jleftu, jleftd, jleftu+jleftd, jleftu-jleftd, occ
       write(16,518)tout, eshift(1,1), occu, occd, occ
     else if (nalf .eq. 2) then
       write(16,518)tout, eshift(1,1), eshift(2,1), occu, occd, occ
       write(13,518)tout, eshift(1,1), eshift(2,1), jleft, jright,                      &
                    jleftu, jleftd, jrightu, jrightd, occu, occd
       write(16,518)tout, eshift(1,1), eshift(2,1), occu, occd, occ
     else ! multi-leads
       write(13,518)tout, (eshift(ialf,1), ialf=1,nalf), (jleads(ialf), ialf=1,nalf), occu
       write(16,518)tout, (eshift(ialf,1), ialf=1,nalf), occu
     end if
  end if
  call flush(13)
  call flush(16)
end do
!
nstep = nstep + 1
tt    = tt + dt
!
!if (lresume) then
!  if (nstep .gt. 0 .and. mod(nstep, nresume) .eq. 0) then
!    call resumejob(tt, 2)
!  end if
!end if
!
inquire(file="stoptd", exist=fexist)
if (fexist) then
  call resumejob(tt, 2)
  write(6,*)
  write(6,*)' STOPTD detected, calculation terminates '
  call flush(6)
  goto 200
end if 
!
goto 100
200 continue
!
deallocate(bessfunc, tcheby, STAT=istat)
deallocate(sigmau, sigmad, zdm0, zdm1, STAT=istat)
!
call cpu_time(cpu2)
write(6,*)
write(6,*)' leaving chebyshev ', cpu2 - cpu1
call flush(6)
!
518  format(f12.6, 2x, 30(e15.6e3, 2x))
1000 format(1x, I10, 2x, e20.6e3, 2x, e20.6e3)
!
end subroutine chebyshev
