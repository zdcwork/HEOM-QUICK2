subroutine calc_adagx_bicg
!
! purpose : calculate the "generalized" matrix multiplication A * x for BICG process,
!           where A is N by N sparse matrix (indexed in indexcoef), x is N by 1 matrix 
!           (one-dimensional array).
!
! The reason that this evaluation is called "generalized" is because some unknowns are
! taken complex conjugate and then multiplied by A(i,j).
!
! Particularly for this program, psi(itier=1) = 1 is a constant (not an unknown), and 
! its relevant coefficient should be filtered out from indexcoef, which would play the
! role of vector B in this linear problem A * x = B
!    
! Here, x is spanned by 3N (N=nunk) arrays, and the output is also 3N-lengthed array.
!
use bicgmod
use matmod
implicit none
include './include/sizes'
include './include/common'
!
integer   :: itier
integer   :: index_nk(MAXTIER)
integer*8 :: lni, lnj, lnk, lnl, lnm
integer   :: ni, nj, nk, ntmp1, ntmp2
real*8    :: dfock
real*8    :: cpu1, cpu2
complex*16 :: cgamma_nk
!
if (igroundsteady .eq. 0) then
 dfock = 0.d0
else if (igroundsteady .eq. 1) then
 dfock = 5.d-1 * (engyshift(1) + engyshift(2)) / hbar
else
 write(6,*)
 write(6,*)' error call calc_adagx_bicg ', igroundsteady
 stop
end if
!
if (memmod .eq. 1) then
 open(unit=17, file='coefindex.data', form='unformatted', status='unknown')
 rewind(17)
 read(17)ntmp1, ntmp2
 read(17)lnm
 if (lnm .ne. lall) then
  write(6,*)
  write(6,*)' error! bad file <coefindex.data> ', lnm, lall
  stop
 end if
end if
lnl = 1
!
!######################################################################################
!                                                                                     !
!                   Follow the sequence settled in indextable                         !
!                                                                                     !
!######################################################################################
!
zetarhsr(1:nunk) = 0.d0 
phirhsr (1:nunk) = 0.d0 
psirhsr (1:nunk) = 0.d0 
zetarhsi(1:nunk) = 0.d0 
phirhsi (1:nunk) = 0.d0 
psirhsi (1:nunk) = 0.d0 
!
! First tier, no balls, i.e., N = 0 in the derivation
!
itier = 1
!
!phirhsr(1) = phirhsr(1) - dfock * phitmpi(1) 
phirhsi(1) = phirhsi(1) - dfock * phitmpr(1)
!
!phirhsi(1) = phirhsi(1) + dfock * phitmpr(1)
phirhsr(1) = phirhsr(1) + dfock * phitmpi(1)
!
do ni=1,nvar
 if (dipm(ni) .le. 0.d0) then
!  zetarhsr(1) = zetarhsr(1) - phitmpr(nfirst(itier + 1) - 1 + ni)
  phirhsr(nfirst(itier + 1) - 1 + ni) = phirhsr(nfirst(itier + 1) - 1 + ni) - zetatmpr(1)
!  zetarhsi(1) = zetarhsi(1) - phitmpi(nfirst(itier + 1) - 1 + ni)
  phirhsi(nfirst(itier + 1) - 1 + ni) = phirhsi(nfirst(itier + 1) - 1 + ni) - zetatmpi(1)
!
 else
!  zetarhsr(1) = zetarhsr(1) - phitmpr(nfirst(itier + 1) - 1 + iredex(ni))
  phirhsr(nfirst(itier + 1) - 1 + iredex(ni)) = phirhsr(nfirst(itier + 1) - 1 + iredex(ni)) - &
                                                zetatmpr(1)
!  zetarhsi(1) = zetarhsi(1) + phitmpi(nfirst(itier + 1) - 1 + iredex(ni))
  phirhsi(nfirst(itier + 1) - 1 + iredex(ni)) = phirhsi(nfirst(itier + 1) - 1 + iredex(ni)) + &
                                                zetatmpi(1)
!
!  phirhsr(1) = phirhsr(1) - 2.d0 * zetatmpr(nfirst(itier + 1) - 1 + ni) + &
!                                    psitmpr(nfirst(itier + 1) - 1 + ni) 
  psirhsr (nfirst(itier + 1) - 1 + ni) = psirhsr (nfirst(itier + 1) - 1 + ni) + phitmpr(1)
  zetarhsr(nfirst(itier + 1) - 1 + ni) = zetarhsr(nfirst(itier + 1) - 1 + ni) - 2.d0 *        &
                                         phitmpr(1)
!
!  phirhsi(1) = phirhsi(1) - 2.d0 * zetatmpi(nfirst(itier + 1) - 1 + ni) + &
!                                    psitmpi(nfirst(itier + 1) - 1 + ni) 
  zetarhsi(nfirst(itier + 1) - 1 + ni) = zetarhsi(nfirst(itier + 1) - 1 + ni) - 2.d0 *        &
                                         phitmpi(1)
  psirhsi (nfirst(itier + 1) - 1 + ni) = psirhsi (nfirst(itier + 1) - 1 + ni) + phitmpi(1)
 end if
end do
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,ntier
  call cpu_time(cpu1)
  do lni=nfirst(itier), nlast(itier)
!
    cgamma_nk = czero
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
    do ni=1,itier-1 
     index_nk(ni) = indextable( ifirst(itier) - 1 +                            &
                   (lni - nfirst(itier)) * (itier - 1) + ni )
     cgamma_nk    = cgamma_nk + cgamma(index_nk(ni)) 
     if (igroundsteady .eq. 1) then
      cgamma_nk = cgamma_nk + eye * dipm(index_nk(ni)) * engyshift(ilead(index_nk(ni))) / hbar
     end if
    end do
!
! same tier contribution (N)
!
!     zetarhsr(lni) = zetarhsr(lni) + dble ( cgamma_nk ) * zetatmpr(lni) -       &
!                                     dimag( cgamma_nk ) * zetatmpi(lni)
     zetarhsr(lni) = zetarhsr(lni) + dble ( cgamma_nk ) * zetatmpr(lni)
     zetarhsi(lni) = zetarhsi(lni) - dimag( cgamma_nk ) * zetatmpr(lni)
!     zetarhsi(lni) = zetarhsi(lni) + dble ( cgamma_nk ) * zetatmpi(lni) +       &
!                                     dimag( cgamma_nk ) * zetatmpr(lni)
     zetarhsi(lni) = zetarhsi(lni) + dble ( cgamma_nk ) * zetatmpi(lni)
     zetarhsr(lni) = zetarhsr(lni) + dimag( cgamma_nk ) * zetatmpi(lni)
!
!     phirhsr(lni)  = phirhsr(lni)  + dble ( cgamma_nk ) * phitmpr(lni)  -       &
!                                   ( dimag( cgamma_nk ) + dfock ) * phitmpi(lni)
     phirhsr(lni) = phirhsr(lni) + dble ( cgamma_nk ) * phitmpr(lni)
     phirhsi(lni) = phirhsi(lni) - ( dimag(cgamma_nk) + dfock ) * phitmpr(lni)
!     phirhsi(lni)  = phirhsi(lni)  + dble ( cgamma_nk ) * phitmpi(lni)  +       &
!                                   ( dimag( cgamma_nk ) + dfock ) * phitmpr(lni)
     phirhsi(lni) = phirhsi(lni) + dble ( cgamma_nk ) * phitmpi(lni)
     phirhsr(lni) = phirhsr(lni) + ( dimag(cgamma_nk) + dfock ) * phitmpi(lni)
!
!     psirhsr(lni)  = psirhsr(lni)  + dble ( cgamma_nk ) * psitmpr(lni)  -       &
!                                     dimag( cgamma_nk ) * psitmpi(lni)
     psirhsr(lni) = psirhsr(lni) + dble ( cgamma_nk ) * psitmpr(lni)
     psirhsi(lni) = psirhsi(lni) - dimag( cgamma_nk ) * psitmpr(lni)
!     psirhsi(lni)  = psirhsi(lni)  + dble ( cgamma_nk ) * psitmpi(lni)  +       &
!                                     dimag( cgamma_nk ) * psitmpr(lni)
     psirhsi(lni) = psirhsi(lni) + dble ( cgamma_nk ) * psitmpi(lni)
     psirhsr(lni) = psirhsr(lni) + dimag( cgamma_nk ) * psitmpi(lni)
!
! nearest lower tier contribution (N - 1)
!
    do ni=1,itier-1
      if (memmod .eq. 0) then
       lnj = indexcoef((lnl - 1) * 2 + 1)
       lnk = indexcoef((lnl - 1) * 2 + 2)
       lnl = lnl + 1
      else if (memmod .eq. 1) then
       read(17)lnj, lnk
      else 
       call findlowertier(itier, index_nk(1), ni, lnj, 0)
       call findlowertier(itier, index_nk(1), ni, lnk, 1)
      end if
!
      if (dipm(index_nk(ni)) .ge. 0.d0) then
!       zetarhsr(lni) = zetarhsr(lni) + dble ( cb(index_nk(ni) ) * phitmpr(lnj) -   &
!                                       dimag( cb(index_nk(ni) ) * phitmpi(lnj) 
!       zetarhsi(lni) = zetarhsi(lni) + dble ( cb(index_nk(ni) ) * phitmpi(lnj) +   &
!                                       dimag( cb(index_nk(ni) ) * phitmpr(lnj) 
        phirhsr(lnj) = phirhsr(lnj) + dble ( cb(index_nk(ni)) ) * zetatmpr(lni)
        phirhsi(lnj) = phirhsi(lnj) - dimag( cb(index_nk(ni)) ) * zetatmpr(lni)
        phirhsi(lnj) = phirhsi(lnj) + dble ( cb(index_nk(ni)) ) * zetatmpi(lni)
        phirhsr(lnj) = phirhsr(lnj) + dimag( cb(index_nk(ni)) ) * zetatmpi(lni)
!
!       psirhsr(lni)  = psirhsr(lni) + dble ( cm(index_nk(ni)) ) * phitmpr(lnj) -   &
!                                      dimag( cm(index_nk(ni)) ) * phitmpi(lnj)
!       psirhsi(lni)  = psirhsi(lni) + dble ( cm(index_nk(ni)) ) * phitmpi(lnj) +   &
!                                      dimag( cm(index_nk(ni)) ) * phitmpr(lnj)
        phirhsr(lnj) = phirhsr(lnj) + dble ( cm(index_nk(ni)) ) * psitmpr(lni)
        phirhsi(lnj) = phirhsi(lnj) - dimag( cm(index_nk(ni)) ) * psitmpr(lni)
        phirhsi(lnj) = phirhsi(lnj) + dble ( cm(index_nk(ni)) ) * psitmpi(lni)
        phirhsr(lnj) = phirhsr(lnj) + dimag( cm(index_nk(ni)) ) * psitmpi(lni)
!
      else 
!       zetarhsr(lni) = zetarhsr(lni) + (-1.d0)**(itier - 1) *                      &
!                                    ( dble ( cd(index_nk(ni)) ) * phitmpr(lnk) +   &
!                                      dimag( cd(index_nk(ni)) ) * phitmpi(lnk) )
!       zetarhsi(lni) = zetarhsi(lni) + (-1.d0)**(itier - 1) *                      &
!                                    ( dimag( cd(index_nk(ni)) ) * phitmpr(lnk) -   &
!                                      dble ( cd(index_nk(ni)) ) * phitmpi(lnk) )
        phirhsr(lnk) = phirhsr(lnk) + (-1.d0)**(itier - 1) * dble ( cd(index_nk(ni)) ) * zetatmpr(lni)
        phirhsi(lnk) = phirhsi(lnk) + (-1.d0)**(itier - 1) * dimag( cd(index_nk(ni)) ) * zetatmpr(lni)
        phirhsr(lnk) = phirhsr(lnk) + (-1.d0)**(itier - 1) * dimag( cd(index_nk(ni)) ) * zetatmpi(lni)
        phirhsi(lnk) = phirhsi(lnk) - (-1.d0)**(itier - 1) * dble ( cd(index_nk(ni)) ) * zetatmpi(lni)
!
!       phirhsr(lni) = phirhsr(lni) + dble ( ck(index_nk(ni)) ) * zetatmpr(lnj) -   &
!                                     dimag( ck(index_nk(ni)) ) * zetatmpi(lnj)
!       phirhsi(lni) = phirhsi(lni) + dble ( ck(index_nk(ni)) ) * zetatmpi(lnj) +   &
!                                     dimag( ck(index_nk(ni)) ) * zetatmpr(lnj)
        zetarhsr(lnj) = zetarhsr(lnj) + dble ( ck(index_nk(ni)) ) * phitmpr(lni)
        zetarhsi(lnj) = zetarhsi(lnj) - dimag( ck(index_nk(ni)) ) * phitmpr(lni)
        zetarhsi(lnj) = zetarhsi(lnj) + dble ( ck(index_nk(ni)) ) * phitmpi(lni)
        zetarhsr(lnj) = zetarhsr(lnj) + dimag( ck(index_nk(ni)) ) * phitmpi(lni)
!
       if (itier .ne. 2) then
!        phirhsr(lni) = phirhsr(lni) - ( dble ( cd(index_nk(ni)) ) * psitmpr(lnj) - &
!                                        dimag( cd(index_nk(ni)) ) * psitmpi(lnj) )
!        phirhsi(lni) = phirhsi(lni) - ( dble ( cd(index_nk(ni)) ) * psitmpi(lnj) + &
!                                        dimag( cd(index_nk(ni)) ) * psitmpr(lnj) )
        psirhsr(lnj) = psirhsr(lnj) - dble ( cd(index_nk(ni)) ) * phitmpr(lni)
        psirhsi(lnj) = psirhsi(lnj) + dimag( cd(index_nk(ni)) ) * phitmpr(lni)
        psirhsi(lnj) = psirhsi(lnj) - dble ( cd(index_nk(ni)) ) * phitmpi(lni)
        psirhsr(lnj) = psirhsr(lnj) - dimag( cd(index_nk(ni)) ) * phitmpi(lni)
       end if
!       psirhsr(lni) = psirhsr(lni) + (-1.d0)**(itier - 2) * (                      &
!                                     dble ( cm(index_nk(ni)) ) * phitmpr(lnk) +    &
!                                     dimag( cm(index_nk(ni)) ) * phitmpi(lnk) ) 
!       psirhsi(lni) = psirhsi(lni) + (-1.d0)**(itier - 2) * (                      &
!                                     dimag( cm(index_nk(ni)) ) * phitmpr(lnk) -    &
!                                     dble ( cm(index_nk(ni)) ) * phitmpi(lnk) )
        phirhsr(lnk) = phirhsr(lnk) + (-1.d0)**(itier - 2) * dble ( cm(index_nk(ni)) ) * psitmpr(lni)
        phirhsi(lnk) = phirhsi(lnk) + (-1.d0)**(itier - 2) * dimag( cm(index_nk(ni)) ) * psitmpr(lni)
        phirhsr(lnk) = phirhsr(lnk) + (-1.d0)**(itier - 2) * dimag( cm(index_nk(ni)) ) * psitmpi(lni)
        phirhsi(lnk) = phirhsi(lnk) - (-1.d0)**(itier - 2) * dble ( cm(index_nk(ni)) ) * psitmpi(lni)
      end if
    end do
!
! nearest higher tier contribution (N + 1)
!
    do ni=1,nvar
      if (itier .lt. ntier) then
!
        if (memmod .eq. 0) then
         lnj = indexcoef((lnl - 1) * 2 + 1)
         lnk = indexcoef((lnl - 1) * 2 + 2)
         lnl = lnl + 1
        else if (memmod .eq. 1) then
         read(17)lnj, lnk
        else 
         call findhighertier(itier, index_nk(1), ni, lnj, 0)
         call findhighertier(itier, index_nk(1), ni, lnk, 1)
        end if
!
        if (dipm(ni) .ge. 0.d0) then
!         zetarhsr(lni) = zetarhsr(lni) + (-1.d0)**itier * phitmpr(lnk)
!         zetarhsi(lni) = zetarhsi(lni) - (-1.d0)**itier * phitmpi(lnk)
          phirhsr(lnk) = phirhsr(lnk) + (-1.d0)**itier * zetatmpr(lni)
          phirhsi(lnk) = phirhsi(lnk) - (-1.d0)**itier * zetatmpi(lni)
!
!         phirhsr(lni) = phirhsr(lni) - 2.d0 * zetatmpr(lnj) + psitmpr(lnj)
!         phirhsi(lni) = phirhsi(lni) - 2.d0 * zetatmpi(lnj) + psitmpi(lnj)
          zetarhsr(lnj) = zetarhsr(lnj) - 2.d0 * phitmpr(lni)
          psirhsr(lnj)  = psirhsr(lnj)  + phitmpr(lni)
          zetarhsi(lnj) = zetarhsi(lnj) - 2.d0 * phitmpi(lni)
          psirhsi(lnj)  = psirhsi(lnj)  + phitmpi(lni)
        else
!         zetarhsr(lni) = zetarhsr(lni) - phitmpr(lnj)
!         zetarhsi(lni) = zetarhsi(lni) - phitmpi(lnj)
          phirhsr(lnj) = phirhsr(lnj) - zetatmpr(lni)
          phirhsi(lnj) = phirhsi(lnj) - zetatmpi(lni)
        end if
      end if
    end do
!
  end do
!
  call cpu_time(cpu2)
!  write(6,100)itier, nlast(itier), cpu2 - cpu1
!  call flush(6)
100 format('itier = ', I3, ' length ', I8, 2x, f12.6)
end do 
!
if (memmod .eq. 1) close(unit=17)
!
end subroutine calc_adagx_bicg
