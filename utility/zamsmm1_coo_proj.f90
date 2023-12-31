subroutine zamsmm1_coo_proj(sidea, transa, transb, korbs, kspin, &
                            zalpha, b, nnz_b, irow_b, icol_b,    & 
                            zbeta, c, nnz_c, irow_c, icol_c, zmtmp, lwork)
!
! Creation and annihilation operators are in sparse matrix form
! Use projected annihilation operators created in <checkamsprojector.f90>
! This subroutine calculates creation/annihilation matrix multiplied
! by genenral complex matrix b, resulting in complex matrix c.
! The matrix multiplication only involves nonzero elements of a and a^dag.
! In doing so the computation cost reduces from N^3 to nnz_ams * N
!
! input matrix b and output matrix c are both in COO sparse format 
!
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent(in)   :: sidea, transa, transb
integer, intent(in)       :: korbs, kspin, nnz_b, nnz_c, lwork
integer,    intent(in)    :: irow_b(*), icol_b(*), irow_c(*), icol_c(*)
complex*16, intent(in)    :: zalpha, zbeta, b(*)
complex*16, intent(inout) :: c(*), zmtmp(lwork,*)
!
integer                   :: ni, nj, nk
integer                   :: mi, mj, mk
integer                   :: irowvec(nnz_b), icolvec(nnz_b)     
logical                   :: lsame
external                  :: lsame
logical                   :: nota, notb, conja, conjb, lefta, righta
complex*16                :: bvec(nnz_b), cvec(nnz_c)
complex*16                :: ctmp1
!
if (nnz_c .eq. 0) return
if (nnz_c .lt. 0) then
   write(6,*)
   write(6,*)'zamsmm1_coo_proj: error, negative nnz_c found ', nnz_c
   stop
end if
if (nnz_b .eq. 0) then
   c(1:nnz_c) = c(1:nnz_c) * zbeta
   return
end if
if (nnz_b .lt. 0) then
   write(6,*)
   write(6,*)'zamsmm1_coo_proj: error, negative nnz_b found ', nnz_b
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
   write(6,*)'zamsmm1_coo_proj: error input for transa ', transa
   stop
end if
if ( (.not. notb) .and. (.not. conjb) ) then
   write(6,*)
   write(6,*)'zamsmm1_coo_proj: error input for transb ', transb
   stop
end if
if ( (.not. lefta) .and. (.not. righta) ) then
   write(6,*)
   write(6,*)'zamsmm1_coo_proj: error input for sidea ', sidea
   stop
end if
!
call zmat_a2b_coo(transb, nnz_b, irow_b, icol_b, b, irowvec, icolvec, bvec)
!
zmtmp(1:nrho,1:nrho) = czero
!
if (lefta) then
   if (nota) then
      do nk=1,nnz_ams(korbs,kspin)
         mi = row_ams(nk,korbs,kspin)
         mk = col_ams(nk,korbs,kspin)
         ctmp1 = cval_ams(nk,korbs,kspin)
         do ni=1,nnz_b
            if ( irowvec(ni) .ne. mk ) cycle
            mj = icolvec(ni)
            zmtmp(mi,mj) = zmtmp(mi,mj) + ctmp1 * bvec(ni)
         end do
      end do
   else 
      do nk=1,nnz_ams(korbs,kspin)
         mi = col_ams(nk,korbs,kspin)
         mk = row_ams(nk,korbs,kspin)
         ctmp1 = dconjg(cval_ams(nk,korbs,kspin))
         do ni=1,nnz_b
            if ( irowvec(ni) .ne. mk ) cycle
            mj = icolvec(ni)
            zmtmp(mi,mj) = zmtmp(mi,mj) + ctmp1 * bvec(ni)
         end do
      end do
   end if
else 
   if (nota) then
      do nk=1,nnz_ams(korbs,kspin)
         mk = row_ams(nk,korbs,kspin)
         mj = col_ams(nk,korbs,kspin)
         ctmp1 = cval_ams(nk,korbs,kspin)
         do ni=1,nnz_b
            if ( icolvec(ni) .ne. mk ) cycle
            mi = irowvec(ni)
            zmtmp(mi,mj) = zmtmp(mi,mj) + ctmp1 * bvec(ni)
         end do  
      end do
   else 
      do nk=1,nnz_ams(korbs,kspin)
         mk = col_ams(nk,korbs,kspin)
         mj = row_ams(nk,korbs,kspin)
         ctmp1 = dconjg(cval_ams(nk,korbs,kspin))
         do ni=1,nnz_b
            if ( icolvec(ni) .ne. mk ) cycle
            mi = irowvec(ni)
            zmtmp(mi,mj) = zmtmp(mi,mj) + ctmp1 * bvec(ni)
         end do  
      end do
   end if
end if
call zmat_extract_coo(nnz_c, irow_c, icol_c, cvec, zmtmp, nrho, nrho, nrho)
!
do ni=1,nnz_c
   c(ni) = zbeta * c(ni) + zalpha * cvec(ni)
end do
!
end subroutine zamsmm1_coo_proj
