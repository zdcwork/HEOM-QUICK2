subroutine prelude_psdfff
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer                 :: ni, nj, nk, istat
integer                 :: isgn, ialf, ispin
real*8                  :: dtmp0, dtmp1, dtmp2
complex*16              :: ctmp0, ctmp1, ctmp2
real*8, allocatable     :: wtmp(:), otmp(:,:)
complex*16, allocatable :: zaptmp(:,:,:,:), zamtmp(:,:,:,:)
complex*16, allocatable :: zbptmp(:,:,:,:), zbmtmp(:,:,:,:)
!
if (.not. psdfff) then
   write(6,*)
   write(6,*)'prelude_psdfff: error! psdfff = ', psdfff
   stop
end if
!
select case (itype_fff) 
  case (1)
    nfff = nfff1
    rfff = rfff1
    allocate(mpfff(nfff), STAT=istat)
    allocate(afff(nfff,nalf), bfff(nfff), STAT=istat)
    mpfff(1:nfff) = mpfff1(1:nfff)
    bfff(1:nfff)  = bfff1(1:nfff) 
    do ialf=1,nalf
       afff(1:nfff,ialf) = afff1(1:nfff) / (dinvbeta(ialf) / hbar * rfff)
    end do
  case (2)
    nfff = nfff2
    rfff = rfff2
    allocate(mpfff(nfff), STAT=istat)
    allocate(afff(nfff,nalf), bfff(nfff), STAT=istat)
    mpfff(1:nfff) = mpfff2(1:nfff)
    bfff(1:nfff)  = bfff2(1:nfff) 
    do ialf=1,nalf
       afff(1:nfff,ialf) = afff2(1:nfff) / (dinvbeta(ialf) / hbar * rfff)
    end do
  case (3)
    nfff = nfff3
    rfff = rfff3
    allocate(mpfff(nfff), STAT=istat)
    allocate(afff(nfff,nalf), bfff(nfff), STAT=istat)
    mpfff(1:nfff) = mpfff3(1:nfff)
    bfff(1:nfff)  = bfff3(1:nfff) 
    do ialf=1,nalf
       afff(1:nfff,ialf) = afff3(1:nfff) / (dinvbeta(ialf) / hbar * rfff)
    end do
!    
  case default
    write(6,*)
    write(6,*)'prelude_psdfff: error! unknown itype_fff ', itype_fff
    stop
end select 
write(6,*)
write(6,*)'prelude_psdfff: parameters adopted in the FFF scheme '
write(6,*)'prelude_psdfff: itype_fff = ', itype_fff
write(6,*)'prelude_psdfff: T0 / T = ', rfff, ' nfff = ', nfff
write(6,*)'prelude_psdfff: powers of monomials: '
write(6,'(8(1x, I4))') (mpfff(ni), ni=1,nfff)
write(6,*)'prelude_psdfff: i, a(i,ialf), b(i) '
do ni=1,nfff
   write(6,1000)ni, (afff(ni,ialf)/hbar, ialf=1,nalf), bfff(ni)
end do
call flush(6)
1000 format(I4, 8(2x, f14.8))
!
numfff = 0
npwfff = 0
do ni=1,nfff
   npwfff = max(npwfff, mpfff(ni) + 1)
   numfff = numfff + mpfff(ni) + 1
end do
write(6,*)
write(6,*)'prelude_psdfff: number of monomial-exponentials: ', numfff
write(6,*)'prelude_psdfff: highest order of monomial:       ', npwfff-1 
call flush(6)
!
allocate(nfff1d(numfff), mfff1d(numfff), STAT=istat)
allocate(temfff(nalf), gfff(numfff,nalf), STAT=istat)
allocate(cfff(npwfff,nfff,nspin,nalf,nsgn), STAT=istat)
allocate(dfff(nfff,nspin,nalf,nsgn), STAT=istat)
allocate(thefff(npwfff,nfff,nspin,nalf,nsgn), STAT=istat)
allocate(wtmp(nalf), otmp(nalf,nspin), STAT=istat)
allocate(zaptmp(nfff,nspin,nalf,nsgn), zamtmp(nfff,nspin,nalf,nsgn), STAT=istat)
allocate(zbptmp(nfff,nspin,nalf,nsgn), zbmtmp(nfff,nspin,nalf,nsgn), STAT=istat)
!
temfff(1:nalf) = dinvbeta(1:nalf) * rfff
write(6,*)'prelude_psdfff: T0 = ', (temfff(ialf), ialf=1,nalf)
call flush(6)
!
nk = 0
do ni=1,nfff
   do nj=1,mpfff(ni)+1
      nk = nk + 1
      nfff1d(nk) = ni
      mfff1d(nk) = nj
   end do
end do
if (nk .ne. numfff) then
    write(6,*)
    write(6,*)'prelude_psdfff: error! wrong numfff ', numfff, nk
    stop
end if
do ni=1,numfff
   nj = nfff1d(ni)
   gfff(ni,1:nalf) = 1.d0 / afff(nj,1:nalf) 
end do
!
wtmp(1:nalf) = bandwidth(1:nalf) / hbar
otmp(1:nalf,1:nspin) = bcenter(1:nalf,1:nspin) / hbar
!write(6,*)
!write(6,*)'prelude_psdfff: bandwidth = ', (wtmp(ialf)*hbar, ialf=1,nalf)
!write(6,*)'prelude_psdfff: bcenter   = ', (otmp(ialf,1)*hbar, ialf=1,nalf)
!write(6,*)
!write(6,*)'prelude_psdfff: gfff      = '
!do ni=1,numfff
!   write(6,*)ni, (gfff(ni,ialf), ialf=1,nalf)
!end do
!call flush(6)

do ialf=1,nalf
   do ispin=1,nspin
      do ni=1,nfff
         zaptmp(ni,ispin,ialf,1) = dcmplx(1.d0 / afff(ni,ialf) + wtmp(ialf),  otmp(ialf,ispin))
         zaptmp(ni,ispin,ialf,2) = dcmplx(1.d0 / afff(ni,ialf) + wtmp(ialf), -otmp(ialf,ispin))
         zamtmp(ni,ispin,ialf,1) = dcmplx(1.d0 / afff(ni,ialf) - wtmp(ialf),  otmp(ialf,ispin))
         zamtmp(ni,ispin,ialf,2) = dcmplx(1.d0 / afff(ni,ialf) - wtmp(ialf), -otmp(ialf,ispin))
      end do
   end do
end do
zbptmp = dconjg(zaptmp)  ! see Cui Lei's draft on December 13, 2018 
zbmtmp = -zamtmp
!
do isgn=1,nsgn
   do ialf=1,nalf
      do ispin=1,nspin
         do ni=1,nfff
            ctmp0 = zaptmp(ni,ispin,ialf,isgn) * zamtmp(ni,ispin,ialf,isgn) 
            dtmp0 = dsgn(isgn) * otmp(ialf,ispin) * afff(ni,ialf)
            select case (mpfff(ni))  ! power of monomial 
              case (0)
                ctmp1 = afff(ni,ialf) * ctmp0
                cfff(1,ni,ispin,ialf,isgn) = -cunity / ctmp1
              case (1)
                ctmp1 = dcmplx(1.d0, dtmp0)
                ctmp2 = afff(ni,ialf)**3 * ctmp0**2
                cfff(1,ni,ispin,ialf,isgn) = -ctmp1 / ctmp2  
                ctmp1 = 2.d0 * afff(ni,ialf)**2 * ctmp0
                cfff(2,ni,ispin,ialf,isgn) = -cunity / ctmp1
              case (2)
                ctmp1 = 4.d0 * afff(ni,ialf)**3 * ctmp0
                ctmp2 = dcmplx(2.d0, dtmp0)
                cfff(1,ni,ispin,ialf,isgn) = -( cunity / zamtmp(ni,ispin,ialf,isgn)**2 + ctmp2 / ctmp0 + &
                                                cunity / zaptmp(ni,ispin,ialf,isgn)**2 ) / ctmp1
                ctmp1 = afff(ni,ialf)**2 * ctmp0
                ctmp2 = dcmplx(1.d0, dtmp0)
                cfff(2,ni,ispin,ialf,isgn) = -( cunity / 8.d0 / ctmp1 + ctmp2 / 2.d0 / ctmp1**2 )
                ctmp1 = afff(ni,ialf)**3 * ctmp0
                cfff(3,ni,ispin,ialf,isgn) = -cunity / 8.d0 / ctmp1
              case default
                write(6,*)
                write(6,*)'prelude_psdfff: error! unknown mpfff for cfff ', ni, mpfff(ni)
                stop
            end select
!             
            ctmp0 = dcmplx(wtmp(ialf), -dsgn(isgn) * otmp(ialf,ispin)) / wtmp(ialf) * afff(ni,ialf)
            ctmp1 = zbptmp(ni,ispin,ialf,isgn) * zbmtmp(ni,ispin,ialf,isgn) * afff(ni,ialf)**2
            select case (mpfff(ni)) 
              case (0)
                dfff(ni,ispin,ialf,isgn) = -ctmp0 / ctmp1 
              case (1)
                dfff(ni,ispin,ialf,isgn) =  ctmp0 / ctmp1**2   !!! no minus sign !?
              case (2)
                dfff(ni,ispin,ialf,isgn) = -ctmp0 / ctmp1**3
              case default
                write(6,*)
                write(6,*)'prelude_psdfff: error! unknown mpfff for dfff ', ni, mpfff(ni)
                stop
            end select
         end do ! ni
      end do ! ispin
   end do ! ialf
end do ! isgn
!
thefff = czero
do isgn=1,nsgn
   do ialf=1,nalf
      do ispin=1,nspin
         do ni=1,nfff
            do nj=2,mpfff(ni)+1   ! skip nj=1
               thefff(nj,ni,ispin,ialf,isgn) = dble(nj-1) * cfff(nj,ni,ispin,ialf,isgn) / cfff(nj-1,ni,ispin,ialf,isgn)
            end do ! nj
         end do ! ni
      end do ! ispin
   end do ! ialf
end do ! isgn
!
! check symmetry of thefff
! 
do ialf=1,nalf
   do ispin=1,nspin
      do ni=1,nfff
         do nj=2,mpfff(ni)+1   ! skip nj=1
            ctmp1 = thefff(nj,ni,ispin,ialf,1)
            ctmp2 = dconjg(thefff(nj,ni,ispin,ialf,2))
            if ( cdabs(ctmp1-ctmp2) .gt. dpico ) then
               write(6,*)
               write(6,*)'prelude_psdfff: error! wrong symmetry with thefff '
               write(6,*)' difference ', ctmp1, ctmp2
               write(6,*)' nj, ni, ispin, ialf ', nj, ni, ispin, ialf
               stop
            end if
         end do ! nj
      end do ! ni
   end do ! ispin
end do ! ialf
!
! check maximal of thefff and minimal of cfff, to make sure:
! [1] thefff is finite 
! [2] there is no trivial cfff 
! so that the hierarchy does not have zero ADO and the coefficients are finite
dtmp1 = 0.d0
dtmp2 = 1.d20
do isgn=1,nsgn
   do ialf=1,nalf
      do ispin=1,nspin
         do ni=1,nfff
            do nj=1,mpfff(ni)+1   
               dtmp1 = max(dtmp1, cdabs(thefff(nj,ni,ispin,ialf,isgn)))
               dtmp2 = min(dtmp2, cdabs(cfff(nj,ni,ispin,ialf,isgn)))
            end do ! nj
         end do ! ni
      end do ! ispin
   end do ! ialf
end do ! isgn
write(6,*)
write(6,*)'prelude_psdfff: max(thefff) = ', dtmp1
if (dtmp1 .gt. 1.d0/dnano) then
    write(6,*)'prelude_psdfff: Warning! max(thefff) is unusually large '
end if
write(6,*)'prelude_psdfff: min(cfff)   = ', dtmp2
if (dtmp2 .lt. dnano) then
    write(6,*)'prelude_psdfff: Warning! min(cfff) is unusually small '
end if
call flush(6)
!
deallocate(wtmp, otmp, STAT=istat)
deallocate(zaptmp, zamtmp, STAT=istat)
deallocate(zbptmp, zbmtmp, STAT=istat)
!
return
end subroutine prelude_psdfff
