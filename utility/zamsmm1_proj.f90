subroutine zamsmm1_proj(sidea, transa, transb, korbs, kspin, zalpha, b, ldb, zbeta, c, ldc)
!
! Creation and annihilation operators are in sparse matrix form
! Use projected annihilation operators created in <checkamsprojector.f90>
! This subroutine calculates creation/annihilation matrix multiplied
! by genenral complex matrix b, resulting in complex matrix c.
! The matrix multiplication only involves nonzero elements of a and a^dag.
! In doing so the computation cost reduces from N^3 to nnz_ams * N
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
complex*16 :: ctmp1
!
if (ldb .ne. nrho .or. ldc .ne. nrho) then
   write(6,*)'zamsmm1_proj: error matrix dimension'
   write(6,*)'zamsmm1_proj: nrho = ', nrho, ' ldb = ', ldb, ' ldc = ', ldc
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
  write(6,*)' error input for transa in zamsmm1_proj', transa
  stop
end if
!
if ( (.not. notb) .and. (.not. conjb) ) then
  write(6,*)
  write(6,*)' error input for transb in zamsmm1_proj', transb
  stop
end if
!
if ( (.not. lefta) .and. (.not. righta) ) then
  write(6,*)
  write(6,*)' error input for sidea in zamsmm1_proj', sidea
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
      do nk=1,nnz_ams(korbs,kspin)
         mi = row_ams(nk,korbs,kspin)
         mk = col_ams(nk,korbs,kspin)
         ctmp1 = cval_ams(nk,korbs,kspin)
         do mj=1,nrho
            zmtmp2(mi,mj) = zmtmp2(mi,mj) + ctmp1 * zmtmp1(mk,mj)
         end do
      end do
   else 
      do nk=1,nnz_ams(korbs,kspin)
         mi = col_ams(nk,korbs,kspin)
         mk = row_ams(nk,korbs,kspin)
         ctmp1 = dconjg(cval_ams(nk,korbs,kspin))
         do mj=1,nrho
            zmtmp2(mi,mj) = zmtmp2(mi,mj) + ctmp1 * zmtmp1(mk,mj)
         end do
      end do
   end if
else 
   if (nota) then
      do nk=1,nnz_ams(korbs,kspin)
         mk = row_ams(nk,korbs,kspin)
         mj = col_ams(nk,korbs,kspin)
         ctmp1 = cval_ams(nk,korbs,kspin)
         do mi=1,nrho
            zmtmp2(mi,mj) = zmtmp2(mi,mj) + ctmp1 * zmtmp1(mi,mk)
         end do
      end do
   else 
      do nk=1,nnz_ams(korbs,kspin)
         mk = col_ams(nk,korbs,kspin)
         mj = row_ams(nk,korbs,kspin)
         ctmp1 = dconjg(cval_ams(nk,korbs,kspin))
         do mi=1,nrho
            zmtmp2(mi,mj) = zmtmp2(mi,mj) + ctmp1 * zmtmp1(mi,mk)
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
end subroutine zamsmm1_proj
