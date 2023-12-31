subroutine inibicg
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
real*8    :: cpu1, cpu2
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
rhsrbicg(1:nunk) = 0.d0
rhsibicg(1:nunk) = 0.d0
!
! Higher tiers, # of balls = itier - 1. Here, itier >= 2, i.e., N >= 1.
!
do itier=2,min(2,ntier)
  call cpu_time(cpu1)
  do lni=nfirst(itier), nlast(itier)
!
! Get the multi-component index {n_k}, where k goes through 1, ..., nvar.   
!
    do ni=1,itier-1 
     index_nk(ni) = indextable( ifirst(itier) - 1 +                            &
                    (lni - nfirst(itier)) * (itier - 1) + ni )
    end do
!
! same tier contribution (N)
!
! null
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
      if (itier .eq. 2 .and. dipm(index_nk(ni)) .le. 0.d0) then
       rhsrbicg(lni) = dble ( cd(index_nk(ni)) ) * dble ( psi(1) ) -   &
                       dimag( cd(index_nk(ni)) ) * dimag( psi(1) )
       rhsibicg(lni) = dble ( cd(index_nk(ni)) ) * dimag( psi(1) ) +   &
                       dimag( cd(index_nk(ni)) ) * dble ( psi(1) )
      end if
    end do
!
! nearest higher tier contribution (N + 1)
!
    do ni=1,nvar
      if (itier .lt. ntier) then
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
!write(6,*)
!write(6,*)' rhsrbicg '
!write(6,*)(rhsrbicg(lni)*hbar**2, lni=1,nunk)
!write(6,*)' rhsibicg '
!write(6,*)(rhsibicg(lni)*hbar**2, lni=1,nunk)
!stop
!
end subroutine inibicg
