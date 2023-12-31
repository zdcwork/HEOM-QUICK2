subroutine zamsmm1(sidea, transa, transb, korbs, kspin, zalpha, b, ldb, zbeta, c, ldc)
!
! Creation and annihilation operators are in sparse matrix form
! only nams nonzero matrix elements, being either +1 or -1.
! This subroutine calculates creation/annihilation matrix multiplied
! by genenral complex matrix b, resulting in complex matrix c.
! The matrix multiplication only involves nonzero elements of a and a^dag.
! In doing so the computation cost reduces from N^3 to nams * N
!
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent(in)   :: sidea, transa, transb
integer, intent(in)       :: korbs, kspin, ldb, ldc
complex*16, intent(in)    :: zalpha, zbeta, b(ldb,*)
complex*16, intent(inout) :: c(ldc,*)
!
integer  :: ni, nj, nk
integer  :: mi, mj, mk
logical  :: lsame
external :: lsame
logical  :: nota, notb, conja, conjb, lefta, righta
!
if (ldb .ne. nrho .or. ldc .ne. nrho) then
   write(6,*)'zamsmm1: error matrix dimension'
   write(6,*)'zamsmm1: nrho = ', nrho, ' ldb = ', ldb, ' ldc = ', ldc
   stop
end if
!
lefta  = lsame(sidea,  'L')
righta = lsame(sidea,  'R')
nota   = lsame(transa, 'N')
notb   = lsame(transb, 'N')
conja  = lsame(transa, 'C')
conjb  = lsame(transb, 'C')
!
if ( (.not. nota) .and. (.not. conja) ) then
  write(6,*)
  write(6,*)' error input for transa in zamsmm1', transa
  stop
end if
!
if ( (.not. notb) .and. (.not. conjb) ) then
  write(6,*)
  write(6,*)' error input for transb in zamsmm1', transb
  stop
end if
!
if ( (.not. lefta) .and. (.not. righta) ) then
  write(6,*)
  write(6,*)' error input for sidea in zamsmm1', sidea
  stop
end if
!
if (notb) then
   do nj=1,nrho
      do ni=1,nrho
         zmtmp1(ni,nj) = b(ni,nj)
      end do
   end do
else
   do nj=1,nrho
      do ni=1,nrho
         zmtmp1(ni,nj) = dconjg(b(nj,ni))
      end do
   end do
end if
!
zmtmp2 = czero
!
if (lefta) then
   if (nota) then
      do nk=1,nams
         mi = rowams(nk,korbs,kspin)
         mk = colams(nk,korbs,kspin)
         do mj=1,nrho
            if (sgnams(nk,korbs,kspin) .eq. 0) then  
               zmtmp2(mi,mj) = zmtmp2(mi,mj) + zmtmp1(mk,mj)
            else
               zmtmp2(mi,mj) = zmtmp2(mi,mj) - zmtmp1(mk,mj)
            end if
         end do
      end do
   else 
      do nk=1,nams
         mi = colams(nk,korbs,kspin)
         mk = rowams(nk,korbs,kspin)
         do mj=1,nrho
            if (sgnams(nk,korbs,kspin) .eq. 0) then  
               zmtmp2(mi,mj) = zmtmp2(mi,mj) + zmtmp1(mk,mj)
            else
               zmtmp2(mi,mj) = zmtmp2(mi,mj) - zmtmp1(mk,mj)
            end if
         end do
      end do
   end if
else 
   if (nota) then
      do nk=1,nams
         mk = rowams(nk,korbs,kspin)
         mj = colams(nk,korbs,kspin)
         do mi=1,nrho
            if (sgnams(nk,korbs,kspin) .eq. 0) then  
               zmtmp2(mi,mj) = zmtmp2(mi,mj) + zmtmp1(mi,mk)
            else
               zmtmp2(mi,mj) = zmtmp2(mi,mj) - zmtmp1(mi,mk)
            end if
         end do
      end do
   else 
      do nk=1,nams
         mk = colams(nk,korbs,kspin)
         mj = rowams(nk,korbs,kspin)
         do mi=1,nrho
            if (sgnams(nk,korbs,kspin) .eq. 0) then  
               zmtmp2(mi,mj) = zmtmp2(mi,mj) + zmtmp1(mi,mk)
            else
               zmtmp2(mi,mj) = zmtmp2(mi,mj) - zmtmp1(mi,mk)
            end if
         end do
      end do
   end if
end if
!
c(1:nrho, 1:nrho) = zbeta * c(1:nrho, 1:nrho) + zalpha * zmtmp2(1:nrho, 1:nrho)
!
zmtmp1 = czero
zmtmp2 = czero
!
end subroutine zamsmm1
