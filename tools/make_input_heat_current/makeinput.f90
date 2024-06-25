program makeinput
implicit none
!
integer :: ni, nj, nk, istat 
integer :: init, ntier, ncor, norbs, nspin, ialf
integer :: maxit0
integer :: fieldtype
integer :: icont, nresume
logical :: lresume, offcor, megaflux, readcp, readmat
real*8  :: aoffL, aoffR, phioffL, phioffR
integer :: nalf
real*8  :: crit
real*8, allocatable  :: bandwidth(:), dinvbeta(:), linewidth(:,:), engyshift(:,:), tchar(:,:)
real*8, allocatable  :: dmtmp1(:,:)
real*8  :: tmax, dt
real*8  :: e1up, e1down, e2up, e2down, uu1, uu2, u12, t12
real*8  :: temp1, temp2, factor1, factor2
real*8,  parameter :: tconst=1.1603d4 ! eV -> Kelvin
real*8  :: eup, edown, uu, ac, edot
integer :: iorbs_dos, ispin_dos 
logical :: lfreq_dos
real*8  :: freq_dos, crit_dos
!
! iv curve
!
real*8  :: vmin, vmax, dv, vv
integer :: nv, iv
!
init      = 3 
ntier     = 5 
ncor      = 8  
norbs     = 1
nspin     = 2
nalf      = 2
!
allocate(bandwidth(nalf), dinvbeta(nalf), linewidth(nalf,norbs), engyshift(nalf,nspin), &
         STAT=istat)
allocate(tchar(nalf,nspin), STAT=istat)
allocate(dmtmp1(norbs,norbs), STAT=istat)
!
tmax      = 200.d0
dt        = 0.02d0
!
! temperature -> beta^(-1)
! 
!temp1     = 0.5d0
!temp2     = 0.5d0
!dinvbeta(1) = temp1 / tconst * 1.d3
!dinvbeta(2) = temp2 / tconst * 1.d3
dinvbeta(1:nalf) = 0.05d0

!
! WBL : bandwidth >> linewidth
!
linewidth(1:nalf,1:norbs) = 0.1d0
!
bandwidth(1:nalf) = 2.d0
!
vmin     = -0.5d0
vmax     =  0.5d0 
!
open(5, file='vin.tmp')
rewind(5)
read(5,*)nv
read(5,*)iv
close(5)
!
dv = (vmax - vmin) / dble(nv - 1)
vv = vmin + dble(iv - 1) * dv
!
engyshift(1:nalf,1:nspin) = 0.d0
!engyshift(1,1:nspin) =-vv*5.d-1-dv*5.d-1
!engyshift(2,1:nspin) = vv*5.d-1+dv*5.d-1
!engyshift(1,1:nspin) = -0.1d0
!engyshift(2,1:nspin) =  0.1d0
engyshift(2,1:nspin) =  1.d-4
!
! write to file
!
open(7, file='input.in', status='unknown')
rewind(7)
write(7,222)init
write(7,222)ntier
write(7,222)ncor
write(7,222)norbs
write(7,222)nspin
write(7,222)nalf
write(7,333)(bandwidth(ni), ni=1,nalf)
write(7,333)((linewidth(ni,nj), nj=1,norbs), ni=1,nalf) 
write(7,333)(dinvbeta(ni), ni=1,nalf)
write(7,333)((engyshift(ni,nj), nj=1,nspin), ni=1,nalf) 
write(7,333)tmax
write(7,333)dt
!
222 format(I6)
333 format(8(e15.6e3, 2x))
!
edot  =  0.0
eup   = -0.2d0
edown = -0.2d0 
uu    =  0.5d0
e1up = 0.d0
e2up = 0.d0
e1down = 0.d0
e2down = 0.d0
uu1 = 0.d0
uu2 = 0.d0
u12 = 0.d0
!
if (norbs .eq. 1 .and. nspin .eq. 2) then
  write(7,*)
  write(7,452)eup, edown, uu
end if
if (norbs .eq. 1 .and. nspin .eq. 1) then
  write(7,*)
  write(7,445)edot
end if
if (norbs .eq. 2 .and. nspin .eq. 2) then
  write(7,*)
  write(7,459)e1up, e1down, e2up, e2down
  write(7,460)uu1, uu2, u12, t12
end if
!
!if (init .ne. 1) then
!  maxit0 = 10000
!  crit   = 1.d-8
!  write(7,*)
!  write(7,454)maxit0, crit
!end if
!
! field
!
fieldtype = 0
tchar(1:nalf,1:nspin) = 1.d-20
!
write(7,*)
write(7,449) fieldtype
write(7,333) ((tchar(ni,nj), nj=1,nspin), ni=1,nalf)
!
icont   = 0
lresume = .true.
nresume = 50
write(7,*)
write(7,451)lresume, nresume, icont
!
! off-diagonal with phi (norbs=2, nspin=1)
!
offcor = .false.
megaflux = .false.
readcp = .true.
readmat = .true. 
aoffL = 1.d0
aoffR = 1.d0
phioffL = 0.d0
phioffR = 0.d0
!
if (offcor) then
   write(7,*)
   write(7,450)offcor
   if (megaflux) then
      if (.not. (norbs .eq. 2 .and. nspin .eq. 1)) then
         write(6,*)'error! megaflux=.true. is only for norbs=2 and nspin=1'
         stop
      end if
      write(7,*)
      write(7,455)megaflux, aoffL, aoffR
      write(7,456)phioffL, phioffR
   else 
      if (readcp) then
         write(7,*)
         write(7,457)readcp, readmat
         do ialf=1,nalf
            dmtmp1(1:norbs,1:norbs) = 0.d0
            do ni=1,norbs
               dmtmp1(ni,ni) = linewidth(ialf,ni)
            end do
            !
            ! assign off-diagonal elements here
            !
            if (readmat) then
               do ni=1,norbs
                  write(7,458) (dmtmp1(ni,nj), nj=1,norbs)
               end do
            else 
               write(7,458) ((dmtmp1(ni,nj), ni=1,norbs), nj=1,norbs)
            end if
         end do
      end if
   end if
end if
!
!write(7,*)
!write(7,*)'$hamil_sys   lgenhs=.true.  $end'                                                                                       
!write(7,*)
!write(7,*)' hsys                           '                                                                                        
!write(7,*)' edot 1 1 -0.314                '                                                                                        
!write(7,*)' edot 2 1 -0.314                '                                                                                        
!write(7,*)' udot 0 1  0.314                '                                                                                        
!write(7,*)' edot 1 2 -0.942                '                                                                                        
!write(7,*)' edot 2 2 -0.942                '                                                                                        
!write(7,*)' udot 0 2  0.942                '                                                                                        
!write(7,*)' edot 1 3 -1.884                '                                                                                        
!write(7,*)' edot 2 3 -1.884                '                                                                                        
!write(7,*)' udot 0 3  1.884                '                                                                                        
!write(7,*)' ecoup  1  2  1  -0.157         '
!write(7,*)' ecoup  2  2  1  -0.157         '
!write(7,*)' ecoup  1  3  2  -0.157         '
!write(7,*)' ecoup  2  3  2  -0.157         '
!
! control convergency criteria
!
write(7,*)
write(7,*)'$converge  crit=', 1.d-9, ' $end'
call flush(7)
!
iorbs_dos = 1
ispin_dos = 1
lfreq_dos = .true.
freq_dos  = vv
crit_dos = 1.d-8
write(7,*)
write(7,*)'$dos  ldos=T iorbs_dos=', iorbs_dos, ' ispin_dos=', ispin_dos
write(7,*)'      lfreq_dos=T  freq_dos=', freq_dos
write(7,*)'      crit_dos=', crit_dos, ' $end'
call flush(7)
!
!write(7,*)
!write(7,*)'$combineleads  lequileads=T  $end '
!call flush(7)
!
write(7,*)
write(7,*)'$jobinfo  lsparse=T  $end'
call flush(7)
!
close(7)
!
444 format(' $para1  eup=', f12.6, 2x, 'edown=', f12.6, 2x, 'uu=', f12.6, 2x, '$end')
445 format(' $para2  edot=', f12.6, 2x, '$end')
449 format(' $field  fieldtype=', I2, 2x, '$end')
450 format(' $bathcorr  offcor=', L7, 2x, '$end')
451 format(' $resume  lresume=', L7, 2x, 'nresume= ', I4, 2x, 'icont= ', I2, 2x, '$end')
452 format(' $para1  eup=', f12.6, 2x, 'edown=', f12.6, 2x, 'uu=', f12.6, 2x, '$end')
454 format(' $converge  maxit0=', I6, 2x, ' crit= ', e16.8e3, 2x, '$end')
455 format(' $flux0  megaflux=' L7, 2x, 'aoffL=', f12.6, 2x, 'aoffR=', f12.6, 2x)
456 format('         phioffL=', f12.6, 2x, 'phioffR=', f12.6, 2x, '$end')
457 format(' $coupling  readcp=', L7, 2x, 'readmat=', L7, 2x, '$end')
458 format(1x, 32(e15.6e3,2x))
459 format(' $para4 e1up=', f12.6, 2x, 'e1down=', f12.6, 2x, 'e2up=', f12.6, 2x, 'e2down=', f12.6)
460 format('         uu1=', f12.6, 2x, 'uu2=', f12.6, 2x, 'u12=', f12.6, 2x, 't12=', f12.6, 2x, '$end')
!
deallocate(bandwidth, dinvbeta, linewidth, engyshift, tchar, STAT=istat)
deallocate(dmtmp1, STAT=istat)
!
end program makeinput
