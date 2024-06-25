subroutine prelude_jacobi(rcgtmp)
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
complex*16, intent(in)    :: rcgtmp(nrho,*)
integer                   :: nrho2, info, istat, nwork, ni, nj, nk, nm, nn, ncount
integer,    allocatable   :: ipiv(:)
real*8,     allocatable   :: rwork(:), tmpval(:)
real*8,     allocatable   :: dmat1(:,:)
complex*16, allocatable   :: zwork(:,:)
complex*16, allocatable   :: cmat1(:,:)
logical                   :: lherm, lread_jacobi
real*8                    :: dfact, e_jacobi
real*8                    :: dtmp1, dtmp2
!
namelist / jacobi / lread_jacobi, e_jacobi
!
nrho2 = nrho**2
nwork = nrho2**2
!
allocate(cmat1(nrho,nrho), STAT=istat)
allocate(dmat1(nrho2,nrho2), STAT=istat)
!
if (igroundsteady .eq. 0) then
   cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho)
else if (igroundsteady .eq. 1) then
   call calcdhs(1.d20, rcgtmp(1,1))
   cmat1(1:nrho,1:nrho) = hs(1:nrho,1:nrho) + dhs(1:nrho,1:nrho)
else
   write(6,*)
   write(6,*)'prelude_jacobi: error! unknown igroundsteady ', igroundsteady
   stop
end if
!
if (.not. allocated(zmatjac)) allocate(zmatjac(nrho2,nrho2), STAT=istat)
if (.not. allocated(dvaljac)) allocate(dvaljac(nrho2), STAT=istat)
!
! diagonalize L
! 
call zhil2liou('l', 'n', nrho, cmat1, nrho, zlmtmp1, nrho2)
call zhil2liou('r', 'n', nrho, cmat1, nrho, zlmtmp2, nrho2)
zmatjac(1:nrho2,1:nrho2) = zlmtmp1(1:nrho2,1:nrho2) - zlmtmp2(1:nrho2,1:nrho2)
call checkhermicity(zmatjac, nrho2, zlmtmp1, nrho2, lherm, dtmp1)
if (.not. lherm) then
   write(6,*)'prelude_jacobi: error! L is not Hermitian '
   stop
end if
!
allocate(zwork(nrho2,nrho2), rwork(3*nrho2-2), tmpval(nrho2), STAT=istat)
call zheev('V', 'u', nrho2, zmatjac, nrho2, dvaljac, zwork(1,1), nwork, rwork, info)
if (info .ne. 0) then
   write(6,*)'prelude_jacobi: error! diagonalization failed 1 ', info
   stop
end if
!
! sort and screen out eigenvectors with zero eigenvalues
!
zwork = czero
ncount = 0
do ni=1,nrho2
   if (dabs(dvaljac(ni)) .le. dpico) cycle
   ncount = ncount + 1
   tmpval(ncount) = dvaljac(ni)
   zwork(1:nrho2,ncount) = zmatjac(1:nrho2,ni)
end do
write(6,*)'prelude_jacobi: ', ncount, ' out of ', nrho2, ' eigenvalues are nonzero '
call flush(6)
nnzjac = ncount
zmatjac(1:nrho2,1:nrho2) = czero
zmatjac(1:nrho2,1:ncount) = zwork(1:nrho2,1:ncount)
dvaljac(1:nrho2) = 0.d0
dvaljac(1:ncount) = tmpval(1:ncount)
!
!write(6,*)
!write(6,*)'prelude_jacobi: zmatjac '
!call cmatout(nrho2, nrho2, zmatjac, dmat1)
!call flush(6)
!
! determine ejac
!
if (nnzjac .gt. 0) then
    dtmp1 = dvaljac(1)
    dtmp2 = dvaljac(1)
    do ni=2,nnzjac
       dtmp1 = min(dtmp1, dvaljac(ni))
       dtmp2 = max(dtmp2, dvaljac(ni))
    end do
    dtmp1 = dtmp2 - dtmp1
else
    dtmp1 = 0.d0
end if
dtmp2 = 0.d0
do ni=1,nsgn
   do nj=1,nalf
      do nk=1,ncor
         do nm=1,nspin
            do nn=1,nvar2
               dtmp2 = max(dtmp2, cdabs(cb(nn,nm,nk,nj,ni)))
            end do
         end do
      end do
   end do
end do
dtmp2 = dsqrt(dble(ntier-1) * dtmp2)
dfact = 1.d0
ejac  = max(dtmp1, dtmp2) * dfact
write(6,*)
write(6,*)'prelude_jacobi: gap of system Liouvillian ', dtmp1*hbar
write(6,*)'prelude_jacobi: sqrt(ntier*max(cb))       ', dtmp2*hbar
write(6,*)'prelude_jacobi: suggested ejac            ', ejac*hbar
call flush(6)
!
lread_jacobi = .false.
e_jacobi = 1.d0
rewind(5)
read(5, jacobi, end=100)
100 continue
if (lread_jacobi) then
    ejac = e_jacobi / hbar
    write(6,*)'prelude_jacobi: ejac read from input ', ejac*hbar
    call flush(6)
end if
if (ejac .le. 0.d0) then
    write(6,*)'prelude_jacobi: error, negative ejac '
    stop
end if
!
deallocate(zwork, rwork, tmpval, cmat1, STAT=istat)
deallocate(dmat1, STAT=istat)
return
end subroutine prelude_jacobi
