subroutine chebyj(rhoinp, jleftu, jrightu, jleftd, jrightd)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in)  :: rhoinp(nrho,nrho,norbs,nalf,nspin)
real*8,     intent(out) :: jleftu, jrightu, jleftd, jrightd
integer*8               :: lni
integer                 :: ni, nj, nk, iorbs, ispin, ialf, imats, isgn
!
jt = 0.d0
do ispin=1,nspin
  do iorbs=1,norbs
    cmtmp1(1:nrho, 1:nrho) = dcmplx(amsall(1:nrho, 1:nrho, iorbs, ispin), 0.d0)
    do ialf=1,nalf
      call zgemm('n', 'n', nrho, nrho, nrho, cunity, cmtmp1, nrho, rhoinp(1,1,iorbs,ialf,ispin), &
                 nrho, czero, cmtmp3, nrho)
      do ni=1,nrho
        jt(ialf,ispin) = jt(ialf,ispin) - 2.d0 * dimag(cmtmp3(ni,ni))
      end do
    end do
  end do
end do
606 continue
jt = jt * jconst
!
jleads = 0.d0
do ialf=1,nalf
   do ispin=1,nspin
      jleads(ialf) = jleads(ialf) + jt(ialf,ispin)
   end do
end do
!
jleftu  = jt(1,1)
if (nalf .ne. 1) jrightu = jt(2,1)
if (nspin .eq. 1) then
  jleftd  = 0.d0
  if (nalf .ne. 1) jrightd = 0.d0
else
  jleftd  = jt(1,2)
  if (nalf .ne. 1) jrightd = jt(2,2)
end if
!
end subroutine chebyj
