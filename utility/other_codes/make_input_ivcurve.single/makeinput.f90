program makeinput
implicit none
!
integer :: ni, nj, nk, istat, funits, runits
integer :: init, ntier, ncor, norbs, nspin
integer :: methodss, itrun, iresi, nmatsu, maxit0
integer :: fieldtype
integer :: icont, nresume
logical :: lresume, onelead, doubledot
logical :: grandom
integer, parameter :: nalf=2
real*8  :: crit
real*8  :: aL, aR, wL, wR
real*8  :: bandwidth(nalf), linewidth(nalf), dinvbeta(nalf), engyshift(nalf)
real*8  :: tmax, dt
real*8  :: temp1, temp2, factor1, factor2
real*8,  parameter :: tconst=1.1603d4 ! eV -> Kelvin
real*8  :: eup, edown, uu, ac, edot
real*8  :: engy01, engy02, t12, u12
real*8  :: gamma0
!
! iv curve
!
real*8  :: vmin, vmax, dv, vv
integer :: nv, iv
!
onelead = .true.
!
init      = 3
ntier     = 3
ncor      = 100
!ncor      = 150
norbs     = 1
nspin     = 2
!
tmax      = 5.d1
dt        = 2.d-2
!
! temperature -> beta^(-1)
! 
!temp1     = 0.5d0
!temp2     = 0.5d0
!dinvbeta(1) = temp1 / tconst * 1.d3
!dinvbeta(2) = temp2 / tconst * 1.d3
!
! WBL : bandwidth >> linewidth
!
gamma0 = 5.d-1
!linewidth(1) = 0.5d0 * gamma0
!linewidth(2) = 0.5d0 * gamma0
!bandwidth(1:2) = 20.d0 * gamma0
!dinvbeta(1:2) = 0.1d0 * gamma0
!
bandwidth(1:2) = 15.d0
linewidth(1:2) = 0.05d0
dinvbeta (1:2) = 2.d-1
!
nv       =  21
vmax     =  1.d1
vmin     =  0.d0
dv       = (vmax - vmin) / dble(nv - 1)
!
open(5, file='vin.tmp')
rewind(5)
read(5,*)iv
close(5)
!
vv = vmin + dble(iv - 1) * dv
!
engyshift(1) =   vv 
!engyshift(2) =  -vv * 5.d-1
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
write(7,333)(bandwidth(ni), ni=1,nalf)
write(7,333)(linewidth(ni), ni=1,nalf)
write(7,333)( dinvbeta(ni), ni=1,nalf)
write(7,333)(engyshift(ni), ni=1,nalf)
write(7,333)tmax
write(7,333)dt
!
222 format(I6)
333 format(6(e14.6e2, 2x))
!
if (onelead) then
  write(7,*)
  write(7,455)onelead
end if
!
eup   =  0.5d0
edown =  2.5d0
uu    =  4.d0 
if (norbs .eq. 1 .and. nspin .eq. 2) then
  write(7,*)
  write(7,457)eup, edown, uu
end if
!
edot  = 0.d0
if (norbs .eq. 1 .and. nspin .eq. 1) then
  write(7,*)
  write(7,445)edot
end if
!
doubledot = .true.
engy01    = 0.5d0
engy02    = -0.5d0
u12       = 0.d0
t12       = 1.d0
if (norbs .eq. 2 .and. nspin .eq. 1) then
  if (doubledot) then
    write(7,*)
    write(7,500)doubledot
  end if
  write(7,*)
  write(7,501)engy01, engy02, u12, t12
end if

!
methodss = 0
itrun    = 0
iresi    = 0
nmatsu   = 200
funits   = 1
runits   = 0
grandom  = .true.
!
if (grandom) then
  write(7,*)
  write(7,456)grandom
end if
!
if (init .ne. 1) then
  maxit0 = 3000
  crit   = 1.d-6
  write(7,*)
  write(7,454)maxit0, crit
end if
!
! field
!
fieldtype = 0
aL = 1.d-20
aR = 1.d-20
wL = 0.1d0
wR = 0.1d0
!
write(7,*)
write(7,446)methodss, itrun
if (iresi .eq. 1) then
  write(7,*)
  write(7,447)iresi, nmatsu
end if
write(7,*)
write(7,448)funits, runits
!
write(7,*)
if (fieldtype .eq. 0) then
  write(7,449)fieldtype, aL, aR
else if (fieldtype .eq. 1) then
  write(7,450)fieldtype, wL, wR
end if
!
icont   = 0
lresume = .true.
nresume = 20
write(7,*)
write(7,451)lresume, nresume, icont
!
close(7)
!
444 format(' $para1  eup=', f12.6, 2x, 'edown=', f12.6, 2x, 'uu=', f12.6, 2x, '$end')
445 format(' $para2  edot=', f12.6, 2x, '$end')
446 format(' $method  methodss=', I2, 2x, ' itrun=', I2, 2x, ' $end')
447 format(' $residue  iresi=', I2, 2x, ' nmatsu=', I5, 2x, ' $end')
448 format(' $units  funits=', I2, 2x, ' runits= ', I2, 2x, ' $end')
449 format(' $field  fieldtype=', I2, 2x, 'aL= ', e14.6e3, 2x, 'aR= ', e14.6e3, 2x, '$end')
450 format(' $field  fieldtype=', I2, 2x, 'wL= ', e14.6e3, 2x, 'wR= ', e14.6e3, 2x, '$end')
451 format(' $resume  lresume=', L7, 2x, 'nresume= ', I4, 2x, 'icont= ', I2, 2x, '$end')
452 format(' $para1  eup=', f12.6, 2x, 'edown=', f12.6, 2x, 'uu=', f12.6, 2x)
453 format('          ac=', f12.6, 2x, '$end')
454 format(' $converge  maxit0=', I6, 2x, ' crit= ', e16.8e3, 2x, '$end')
455 format(' $lead  onelead=', L7, 2x, '$end')
456 format(' $guess0  grandom=', L7, 2x, '$end')
457 format(' $para1  eup=', f12.6, 2x, 'edown=', f12.6, 2x, 'uu=', f12.6, 2x, ' $end')
500 format(' $ddots  doubledot=', L7, 2x, ' $end')
501 format(' $para3 engy01=', f12.6, 2x, 'engy02=', f12.6, 2x, 'u12=', f12.6,2x,'t12=', f12.6,2x,' $end')
!
end program makeinput
