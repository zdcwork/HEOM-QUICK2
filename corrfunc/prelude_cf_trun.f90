subroutine prelude_cf_trun
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer                   :: nrho2, info, istat, nwork, ni, nj
integer                   :: nrho2d, isgn, ncount
integer,    allocatable   :: ipiv(:)
real*8,     allocatable   :: rwork(:), tmpval(:)
complex*16, allocatable   :: zwork(:,:)
complex*16, allocatable   :: cmat1(:,:)
logical                   :: lherm
real*8                    :: dtmp1
!
if (.not. ltrun_der) then
   write(6,*)
   write(6,*)'prelude_cf_trun: error entrance ', ltrun_der
   stop
end if
nrho2 = nrho**2
nwork = nrho2**2
allocate(cmat1(nrho,nrho), STAT=istat)
!
if (.not. allocated(zmatder)) allocate(zmatder(nrho2,nrho2), STAT=istat)
if (.not. allocated(dvalder)) allocate(dvalder(nrho2), STAT=istat)
!
! diagonalize L
! 
cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
call zhil2liou('l', 'n', nrho, cmat1, nrho, zlmtmp1, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zlmtmp2, nrho2)
zmatder(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
call checkhermicity(zmatder, nrho2, zlmtmp1, nrho2, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)'prelude_cf_trun: error! L is not Hermitian '
   stop
end if
!
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), STAT=istat)
call zheev('V', 'u', nrho2, zmatder, nrho2, dvalder, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'prelude_cf_trun: error! diagonalization failed 1 ', info
   stop
end if
deallocate(zwork, rwork, STAT=istat)
!
nrho2d = nrho2 * 2
nwork  = nrho2d**2
if (.not. allocated(zmatder_cf)) allocate(zmatder_cf(nrho2d,nrho2d,nsgn), STAT=istat)
if (.not. allocated(dvalder_cf)) allocate(dvalder_cf(nrho2d,nsgn), STAT=istat)
zmatder_cf = czero
dvalder_cf = 0.d0
!
call zhil2liou('l', 'n', nrho, cmat1, nrho, zlmtmp1, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zlmtmp2, nrho2)
do isgn=1,nsgn
   zmatder_cf(1:nrho2,1:nrho2,isgn) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
   zmatder_cf(nrho2+1:nrho2d,nrho2+1:nrho2d,isgn) = zmatder_cf(1:nrho2,1:nrho2,isgn)
   if (lfreq) then
      do ni=1,nrho2
         zmatder_cf(ni,nrho2+ni,isgn) = -eye * dipm(isgn) * dfreq
         zmatder_cf(nrho2+ni,ni,isgn) =  eye * dipm(isgn) * dfreq
      end do
   end if
end do
allocate(zwork(nrho2d,nrho2d), rwork(3*nrho2d-2), tmpval(nrho2d), STAT=istat)
do isgn=1,nsgn
   call zheev('V', 'u', nrho2d, zmatder_cf(1,1,isgn), nrho2d, dvalder_cf(1,isgn), &
              zwork(1,1), nwork, rwork, info)
   if (info .ne. 0) then
      write(6,*)'prelude_cf_trun: error! diagonalization failed 2 ', info, isgn
      stop
   end if
   zwork = czero
   tmpval = 0.d0
   ncount = 0
   do ni=1,nrho2d
      if (dabs(dvalder_cf(ni,isgn)) .le. dpico) cycle
      ncount = ncount + 1
      zwork(1:nrho2d,ncount) = zmatder_cf(1:nrho2d,ni,isgn)
      tmpval(ncount) = dvalder_cf(ni,isgn)
   end do  
   nnzeig_cf(isgn) = ncount
   zmatder_cf(1:nrho2d,1:nrho2d,isgn) = zwork(1:nrho2d,1:nrho2d)
   dvalder_cf(1:nrho2d,isgn) = tmpval(1:nrho2d)
   write(6,*)'prelude_cf_trun: ', ncount, ' out of ', nrho2d, ' eigenvalues are nonzero ', isgn
end do
call flush(6)
deallocate(zwork, rwork, tmpval, cmat1, STAT=istat)
!
!write(6,*)'prelude_cf_trun: leaving subroutine'
!call flush(6)
return
end subroutine prelude_cf_trun
