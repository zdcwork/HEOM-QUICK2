subroutine solvebicg
use bicgmod
use matmod
use matmod_omp
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer*8 :: lni, lnj, lnk, ltmp1, ltmp2
integer   :: ni, nj, nk, nl
integer   :: istat, iconvg, iter, ierror
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6
real*8    :: dsav1, errini, erriter
real*8    :: dnorm1, dnorm2, ratio12
real*8    :: cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cpu7, cpu8
real*8    :: dprev1, dprev2, dprev3
real*8    :: ddotlong, cdotlong
real*8    :: alphak, betak, res1, res2, occu, occd, occ
external  :: ddotlong, cdotlong
complex*16, allocatable :: sigmau(:), sigmad(:)
!
call cpu_time(cpu3)
write(6,*)
write(6,*)' entering solvebicg, nrho = ', nrho
call flush(6)
!
ierror = -1
iconvg = 0
iter   = 0
ltmp1  = nunk * nrho**2
ltmp2  = ltmp1 * 2
dprev3 = crit
!
if (lhf .and. .not. lfixhs) then  
   write(6,*)'solvebicg: HF with BICG does not work. Abort...'
   stop
end if
!
write(6,*)
write(6,*)' length of bicg vector ', ltmp2
!
allocate(sigmau(norbs), sigmad(norbs), STAT=istat)
!
call pre_atrx_bicg
!
call allocate_tmpauxmat_omp
!
!---------------------debug
!call checkatrx
!stop
!---------------------debug
!
! In the present BICG scheme, only the lower half of reduced density matrix
! (including all off-diagonal elements, and the real parts of its diagonal elements)
! namely, symmetry of rho is explicitly treated in the bicg algorithm 
! (or the minimal number of independent unknowns are explicitly solved)
!
sk(1:nrho, 1:nrho, 1:nunk) = rho(1:nrho, 1:nrho, 1:nunk)
do ni=1,nrho
  do nj=ni+1,nrho
    sk(ni,nj,1) = czero
  end do
end do
do ni=1,nrho
  sk(ni,ni,1) = dcmplx(dble(sk(ni,ni,1)), 0.d0)
end do
!
call calc_ax_bicg_omp(nrho, sk, rcgrhs)
!
rk(1:nrho, 1:nrho, 1:nunk) = -rcgrhs(1:nrho, 1:nrho, 1:nunk)
rk(1,1,1) = rk(1,1,1) + cunity
!
dnorm1 = cdotlong(ltmp1, rk(1,1,1), rk(1,1,1))
dtmp1  = dsqrt(dnorm1)
write(6,*)'solvebicg: initial norm ', dtmp1
call flush(6)
errini = dtmp1
!
! fast exit if initial rho is close enough to the real solution
!
if (errini .le. dprev3) then
   iconvg = 1
   write(6,*)'solvebicg: initial norm is smaller than criterion', dprev3
   goto 200
end if
!
sk(1:nrho, 1:nrho, 1:nunk) = rk(1:nrho, 1:nrho, 1:nunk)
dk(1:nrho, 1:nrho, 1:nunk) = rk(1:nrho, 1:nrho, 1:nunk)
fk(1:nrho, 1:nrho, 1:nunk) = sk(1:nrho, 1:nrho, 1:nunk)
!
write(6,*)
call cpu_time(cpu1)
!
100 continue
!
call cpu_time(cpu5)
call calc_ax_bicg_omp(nrho, dk, rcgrhs)
call cpu_time(cpu6)
!
res1   = cdotlong(ltmp1, sk(1,1,1), rk(1,1,1))
res2   = cdotlong(ltmp1, fk(1,1,1), rcgrhs(1,1,1))
alphak = res1 / res2
dsav1  = res1
!
rho(1:nrho, 1:nrho, 1:nunk) = rho(1:nrho, 1:nrho, 1:nunk) + dk(1:nrho, 1:nrho, 1:nunk) * alphak
! 
! important: recover the lower half of rho by Hermicity
!  
do ni=1,nrho
  do nj=ni+1,nrho
    rho(ni,nj,1) = dconjg(rho(nj,ni,1))
  end do
end do
!
rk(1:nrho, 1:nrho, 1:nunk) = rk(1:nrho, 1:nrho, 1:nunk) - rcgrhs(1:nrho, 1:nrho, 1:nunk) * alphak
!
call cpu_time(cpu7)
call calc_atrx_bicg_omp(nrho, fk, rcgrhs)
call cpu_time(cpu8)
!
sk(1:nrho, 1:nrho, 1:nunk) = sk(1:nrho, 1:nrho, 1:nunk) - rcgrhs(1:nrho, 1:nrho, 1:nunk) * alphak
!
dnorm2  = cdotlong(ltmp1, rk(1,1,1), rk(1,1,1))
ratio12 = dnorm2 / dnorm1
!
erriter = dsqrt(dnorm2)
if (erriter .le. dprev3) iconvg = 1
!
call calcocc(occu, occd, sigmau, sigmad)
!
call cpu_time(cpu2)
if (nspin .eq. 1) then
!  write(6,518)iter, res1, res2, ratio12, occu, cpu2 - cpu1
  write(6,518)iter, res1, res2, erriter, occu, cpu2 - cpu1
else 
!  write(6,518)iter, res1, res2, ratio12, occu, occd, cpu2 - cpu1, cpu6 - cpu5, cpu8 - cpu7
  write(6,518)iter, res1, res2, erriter, occu, occd, cpu2 - cpu1, cpu6 - cpu5, cpu8 - cpu7
end if
call flush(6)
call cpu_time(cpu1)
518 format(I6, 2x, 10(e12.4e2, 2x))
!
if (iconvg .eq. 1) then
 write(6,*)
 write(6,*)' converged, final norm ', erriter
 call flush(6)
 ierror = 0
 goto 200
end if
!
res1 = cdotlong(ltmp1, sk(1,1,1), rk(1,1,1))
res2 = dsav1
!
betak = res1 / res2
!
dk(1:nrho, 1:nrho, 1:nunk) = rk(1:nrho, 1:nrho, 1:nunk) + dk(1:nrho, 1:nrho, 1:nunk) * betak
fk(1:nrho, 1:nrho, 1:nunk) = sk(1:nrho, 1:nrho, 1:nunk) + fk(1:nrho, 1:nrho, 1:nunk) * betak
!
111 format(I4, 2x, ' check ', 2x, 6(e15.6e3,2x))
call flush(6)
!
iter = iter + 1
if (iter .gt. maxit0) then
 write(6,*)
 write(6,*)' error! maxit0 reached without convergency ', iter
! stop
 goto 200
end if
!
goto 100
!
200 continue
!
write(6,*)
write(6,*)'solvebicg: steady state solved by BICG method '
!
call outsteady(iter)
!
call cpu_time(cpu4)
deallocate(sigmau, sigmad, STAT=istat)
!
write(6,*)
write(6,*)' leaving solvebicg ', cpu4 - cpu3
call flush(6)
end subroutine solvebicg
