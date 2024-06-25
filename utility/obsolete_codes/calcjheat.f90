subroutine calcjheat(hleftu, hrightu, hleftd, hrightd)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
real*8, intent(out) :: hleftu, hrightu, hleftd, hrightd
integer*8           :: lni
integer             :: ni, nj, nk, iorbs, ispin, ialf, imats, isgn
integer             :: ntmp1, ntmp2
!
! important : the heat current expression used below is WRONG, a correct
!             expression is yet to be achieved.
!
if (.not. (mfdjob .and. thermopower)) then
  write(6,*)
  write(6,*)' error! <calcjheat> is exclusively for mfdjob and thermopower '
  stop
end if
!
jheat = 0.d0
do ispin=1,nspin
  do iorbs=1,norbs
    cmtmp1(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, iorbs, ispin), 0.d0)
    do ialf=1,nalf
      cmtmp2 = czero
      do imats=1,ncor
        call look4drawer(1, ialf, iorbs, ispin, imats, ni)
        cmtmp2(1:nrho, 1:nrho) = cmtmp2(1:nrho, 1:nrho) + rho(1:nrho, 1:nrho, nfirst(2)-1+ni) * &
                                                          wkln(imats)           ! this is wrong!
      end do
      call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, cmtmp2, nrho, czero, cmtmp3, nrho)
      do ni=1,nrho
        jheat(ialf,ispin) = jheat(ialf,ispin) - 2.d0 * dimag(cmtmp3(ni,ni))
      end do
    end do
  end do
end do
606 continue
!
jheat = jheat * jconst
!
hleftu  = jheat(1,1)
hrightu = jheat(2,1)
if (nspin .eq. 1) then
  hleftd  = 0.d0
  hrightd = 0.d0
else
  hleftd  = jheat(1,2)
  hrightd = jheat(2,2)
end if
!
end subroutine calcjheat
