subroutine buildspin3d
use matmod
use tmpmatmod
implicit none
!
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, istat, ixyz, irho, jrho, iorbs
real*8  :: dtmp1, dtmp2, dtmp3
character*120 :: cline
!
if (.not. lspin3d .or. nspin .ne. 2) then
   write(6,*)'buildspin3d: error entrance, lspin3d, nspin ', lspin3d, nspin
   stop
end if
!
! build Pauli matrices (with 1/2 prefactor)
!
allocate(pauli(nspin,nspin,3), STAT=istat)
!  x
pauli(1,1,1) = czero                            ! 1/2 * ( 0  1 )
pauli(1,2,1) = dcmplx(0.5d0, 0.d0)              !       ( 1  0 )
pauli(2,1,1) = dcmplx(0.5d0, 0.d0)
pauli(2,2,1) = czero
!  y
pauli(1,1,2) = czero                            ! 1/2 * ( 0 -i )
pauli(1,2,2) = dcmplx(0.d0, -0.5d0)             !       ( i  0 )
pauli(2,1,2) = dcmplx(0.d0,  0.5d0)
pauli(2,2,2) = czero
!  z
pauli(1,1,3) = dcmplx(0.5d0, 0.d0)              ! 1/2 * ( 1  0 )
pauli(1,2,3) = czero                            !       ( 0 -1 ) 
pauli(2,1,3) = czero
pauli(2,2,3) = dcmplx(-0.5d0, 0.d0)
!
allocate(sopr(nrho,nrho,3,norbs), STAT=istat)
!
! S_i = \sum_ss' c^\dag_is pauli_ss' c_is' 
!
do iorbs=1,norbs
   do ixyz=1,3
      cmtmp1(1:nrho,1:nrho) = czero
      do ni=1,nspin
         do nj=1,nspin
            if ( cdabs(pauli(ni,nj,ixyz)) .lt. dpico ) cycle
            cmtmp2(1:nrho,1:nrho) = dcmplx(amsall(1:nrho,1:nrho,iorbs,ni), 0.d0)
            cmtmp3(1:nrho,1:nrho) = dcmplx(amsall(1:nrho,1:nrho,iorbs,nj), 0.d0)
            call zgemm('c', 'n', nrho, nrho, nrho, pauli(ni,nj,ixyz), cmtmp2, nrho,  &
                       cmtmp3, nrho, cunity, cmtmp1, nrho)
         end do 
      end do
      sopr(1:nrho,1:nrho,ixyz,iorbs) = cmtmp1(1:nrho,1:nrho)
   end do
end do
!
! check commutation relation [Sx, Sy] = iSz 
! also, [Sy, Sz] = iSx, [Sz, Sx] = iSy
!
do iorbs=1,norbs
   cmtmp1(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,1,iorbs)
   cmtmp2(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,2,iorbs)
   cmtmp3(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,3,iorbs)
! [x, y]
   call zgemm('n', 'n', nrho, nrho, nrho,  cunity, cmtmp1, nrho, cmtmp2, nrho,  &
               czero, cmtmp4, nrho)
   call zgemm('n', 'n', nrho, nrho, nrho, -cunity, cmtmp2, nrho, cmtmp1, nrho,  &
              cunity, cmtmp4, nrho)
   cmtmp5(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) - eye * cmtmp3(1:nrho,1:nrho)
   call cmaxmat(nrho, nrho, cmtmp5, nrho, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(6,*)'buildspin3d: error! [Sx, Sy] = iSz violated for iorbs = ', iorbs
      write(6,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
! [y, z]
   call zgemm('n', 'n', nrho, nrho, nrho,  cunity, cmtmp2, nrho, cmtmp3, nrho,  &
               czero, cmtmp4, nrho)
   call zgemm('n', 'n', nrho, nrho, nrho, -cunity, cmtmp3, nrho, cmtmp2, nrho,  &
              cunity, cmtmp4, nrho)
   cmtmp5(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) - eye * cmtmp1(1:nrho,1:nrho)
   call cmaxmat(nrho, nrho, cmtmp5, nrho, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(6,*)'buildspin3d: error! [Sy, Sz] = iSx violated for iorbs = ', iorbs
      write(6,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
! [z, x]
   call zgemm('n', 'n', nrho, nrho, nrho,  cunity, cmtmp3, nrho, cmtmp1, nrho,  &
               czero, cmtmp4, nrho)
   call zgemm('n', 'n', nrho, nrho, nrho, -cunity, cmtmp1, nrho, cmtmp3, nrho,  &
              cunity, cmtmp4, nrho)
   cmtmp5(1:nrho,1:nrho) = cmtmp4(1:nrho,1:nrho) - eye * cmtmp2(1:nrho,1:nrho)
   call cmaxmat(nrho, nrho, cmtmp5, nrho, dtmp1)
   if (dtmp1 .ge. dpico) then
      write(6,*)'buildspin3d: error! [Sz, Sx] = iSy violated for iorbs = ', iorbs
      write(6,*)'buildspin3d: largest deviation = ', dtmp1
      stop
   end if
end do
write(6,*)'buildspin3d: spin operator constructed and verified '
call flush(6)
!
! read 3d magnetic field
!
if (lbfield3d) then
   allocate(dbf3d(3,norbs), STAT=istat)
   rewind(5)
   find: do
      read(5, '(A120)') cline
      istat = index(cline, '3d magnetic field')
      if (istat > 0) exit find
   end do find
   if (lbsite) then  ! site-specific magnetic field
      do iorbs=1,norbs
         read(5,*,iostat=istat) (dbf3d(ni,iorbs), ni=1,3)
         if (istat .ne. 0) then
            write(6,*)'buildspin3d: error reading 3d magnetic field for iorbs = ', iorbs
            stop
         end if
      end do
   else ! homogeneous magnetic field for all sites
      read(5,*,iostat=istat) (dbf3d(ni,1), ni=1,3)
      if (istat .ne. 0) then
         write(6,*)'buildspin3d: error reading 3d magnetic field! '
         stop
      end if
      do iorbs=2,norbs
         dbf3d(1:3,iorbs) = dbf3d(1:3,1)
      end do
   end if
!
   dbf3d(1:3,1:norbs) = dbf3d(1:3,1:norbs) / hbar
   write(6,*)
   if (lbsite) then
      write(6,*)'buildspin3d: site-specific magnetic field in (x,y,z) direction' 
      write(6,*)'buildspin3d: H(i) = -[dbf3d(1,i)*Sx(i) + dbf3d(2,i)*Sy(i) +   '
      write(6,*)'buildspin3d:                             dbf3d(3,i)*Sz(i)]    '
      do iorbs=1,norbs
         write(6,'(I3, 1x, 3(e15.6e3, 2x))') iorbs, (dbf3d(ni,iorbs)*hbar, ni=1,3)
      end do
   else 
      write(6,*)'buildspin3d: homogeneous magnetic field in (x,y,z) direction' 
      write(6,*)'buildspin3d: H = -[dbf3d(1)*Sx + dbf3d(2)*Sy + dbf3d(3)*Sz] '
      write(6,'(3(e15.6e3, 2x))') (dbf3d(ni,1)*hbar, ni=1,3)
   end if
   call flush(6)
end if
!
! read zero-field (crystal field) splitting parameters
! 
if (lzfs) then
    d_xx = d_xx / hbar
    d_yy = d_yy / hbar
    d_zz = d_zz / hbar
    d_xx_1 = d_xx_1 / hbar
    d_yy_1 = d_yy_1 / hbar
    d_zz_1 = d_zz_1 / hbar
    d_xx_2 = d_xx_2 / hbar
    d_yy_2 = d_yy_2 / hbar
    d_zz_2 = d_zz_2 / hbar
    write(6,*)
    write(6,*)'buildspin3d: zero field splitting activated           '
    write(6,*)'buildspin3d: H = Dxx * Sx^2 + Dyy * Sy^2 + Dzz * Sz^2 '
    dtmp1 = d_xx**2 + d_yy**2 + d_zz**2
    if (dtmp1 .gt. dpico) then
        write(6,*)'buildspin3d: Dxx = ', d_xx * hbar
        write(6,*)'buildspin3d: Dyy = ', d_yy * hbar
        write(6,*)'buildspin3d: Dzz = ', d_zz * hbar
    end if
    dtmp1 = d_xx_1**2 + d_yy_1**2 + d_zz_1**2
    if (dtmp1 .gt. dpico) then
        write(6,*)'buildspin3d: Dxx_1 = ', d_xx_1 * hbar
        write(6,*)'buildspin3d: Dyy_1 = ', d_yy_1 * hbar
        write(6,*)'buildspin3d: Dzz_1 = ', d_zz_1 * hbar
    end if
    dtmp1 = d_xx_2**2 + d_yy_2**2 + d_zz_2**2
    if (dtmp1 .gt. dpico) then
        write(6,*)'buildspin3d: Dxx_2 = ', d_xx_2 * hbar
        write(6,*)'buildspin3d: Dyy_2 = ', d_yy_2 * hbar
        write(6,*)'buildspin3d: Dzz_2 = ', d_zz_2 * hbar
    end if
    call flush(6)
end if
!
allocate(sdot(3,norbs), STAT=istat)
!
return
end subroutine buildspin3d
!
!
subroutine calcspinmoment(iorbs, dmoment, rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!               <Sx> = tr(Sx * rho), same for y, z
integer,    intent(in)  :: iorbs
real*8,     intent(out) :: dmoment(*)
complex*16, intent(in)  :: rhoinp(nrho,*)
!
integer :: ni, ixyz
!
do ixyz=1,3
   dmoment(ixyz) = 0.d0
   if (lprgen) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, sopr(1,1,ixyz,iorbs), nrho, czero, cmtmp1, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
   else
       cmtmp2(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,ixyz,iorbs)
   end if
   call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, rhoinp, nrho, czero, cmtmp1, nrho)
   do ni=1,nrho
      dmoment(ixyz) = dmoment(ixyz) + dble(cmtmp1(ni,ni))
   end do
end do
!
return
end subroutine calcspinmoment
!
!
subroutine calcspinproduct(iorbs, iorbs2, dout, rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
! < S_i * S_j > 
integer,   intent(in)  :: iorbs, iorbs2
real*8,    intent(out) :: dout
complex*16, intent(in) :: rhoinp(nrho,*)
!
integer :: ni, ixyz
!
cmtmp1(1:nrho,1:nrho) = czero
do ixyz=1,3
   if (lprgen) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, sopr(1,1,ixyz,iorbs), nrho, czero, cmtmp2, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp2, nrho, czero, cmtmp3, nrho)
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, sopr(1,1,ixyz,iorbs2), nrho, czero, cmtmp2, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp2, nrho, czero, cmtmp4, nrho)
   else
       cmtmp3(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,ixyz,iorbs)
       cmtmp4(1:nrho,1:nrho) = sopr(1:nrho,1:nrho,ixyz,iorbs2)
   end if
   call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp3, nrho, cmtmp4, nrho, cunity, cmtmp1, nrho)
end do
call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, rhoinp, nrho, czero, cmtmp2, nrho)
!
dout = 0.d0
do ni=1,nrho
   dout = dout + dble(cmtmp2(ni,ni))
end do
!
return
end subroutine calcspinproduct
!
!
subroutine calctotalspin2(dout, rhoinp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
! < Sx^2 >, < Sy^2 >, < Sz^2 > (here, Sx, Sy, Sz correspond to total spin)
!
real*8,     intent(out) :: dout(*)
complex*16, intent(in)  :: rhoinp(nrho,*)
!
integer :: ni, ixyz, iorbs
!
do ixyz=1,3
   cmtmp1(1:nrho,1:nrho) = czero
   do iorbs=1,norbs
      cmtmp1(1:nrho,1:nrho) = cmtmp1(1:nrho,1:nrho) + sopr(1:nrho,1:nrho,ixyz,iorbs)
   end do
   if (lprgen) then
       call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
       call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp2, nrho, czero, cmtmp1, nrho)
   end if
   call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp1, nrho, czero, cmtmp2, nrho)
   call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp2, nrho, rhoinp, nrho, czero, cmtmp1, nrho) 
   dout(ixyz) = 0.d0
   do ni=1,nrho
      dout(ixyz) = dout(ixyz) + dble(cmtmp1(ni,ni))
   end do
end do
!
return
end subroutine calctotalspin2
