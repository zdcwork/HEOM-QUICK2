subroutine zamsmm2_coo(sidea, transa, transb, korbs, kspin, zalpha, b, ldb, zbeta,  &
                       c, nnz_c, irow_c, icol_c, zmtmp, ldtmp)
!
! zalpha * op(a_ms) * op(B) + zbeta * C ==> C
! a_ms and the output C are in sparse COO format, input B is a dense matrix
!
use matmod
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent(in)   :: sidea, transa, transb
integer,    intent(in)    :: korbs, kspin, ldb, nnz_c, ldtmp
integer,    intent(in)    :: irow_c(*), icol_c(*)
complex*16, intent(in)    :: zalpha, zbeta, b(ldb,*)
complex*16, intent(inout) :: c(*), zmtmp(ldtmp,*)
!
integer                   :: ni, nj, nk
integer                   :: mi, mj, mk
logical                   :: lsame
external                  :: lsame
logical                   :: nota, conja, lefta, righta
complex*16                :: cvec(nnz_c)
!
if (nnz_c .eq. 0) return
if (nnz_c .lt. 0) then
   write(6,*)
   write(6,*)'zamsmm2_coo: error, negative nnz_c found ', nnz_c
   stop
end if
!
lefta  = lsame(sidea,  'L')
righta = lsame(sidea,  'R')
nota   = lsame(transa, 'N')
conja  = lsame(transa, 'C')
!
if ( (.not. nota) .and. (.not. conja) ) then
   write(6,*)
   write(6,*)'zamsmm2_coo: error input for transa ', transa
   stop
end if
if ( (.not. lefta) .and. (.not. righta) ) then
   write(6,*)
   write(6,*)'zamsmm2_coo: error input for sidea ', sidea
   stop
end if
!
call zmat_a2b_dns(transb, nrho, b, ldb, zmtmp, ldtmp)
cvec(1:nnz_c) = czero
!
if (lefta) then
   if (nota) then
      do ni=1,nnz_c
         mi = irow_c(ni)
         mj = icol_c(ni)
         do nk=1,nams
            if ( rowams(nk,korbs,kspin) .ne. mi ) cycle
            if ( sgnams(nk,korbs,kspin) .eq. 0 ) then
               cvec(ni) = cvec(ni) + zmtmp(colams(nk,korbs,kspin),mj)
            else
               cvec(ni) = cvec(ni) - zmtmp(colams(nk,korbs,kspin),mj)
            end if
         end do
      end do
   else 
      do ni=1,nnz_c
         mi = irow_c(ni)
         mj = icol_c(ni)
         do nk=1,nams
            if ( colams(nk,korbs,kspin) .ne. mi ) cycle
            if ( sgnams(nk,korbs,kspin) .eq. 0 ) then
               cvec(ni) = cvec(ni) + zmtmp(rowams(nk,korbs,kspin),mj)
            else
               cvec(ni) = cvec(ni) - zmtmp(rowams(nk,korbs,kspin),mj)
            end if
         end do
      end do
   end if
else 
   if (nota) then
      do ni=1,nnz_c
         mi = irow_c(ni)
         mj = icol_c(ni)
         do nk=1,nams
            if ( colams(nk,korbs,kspin) .ne. mj ) cycle
            if ( sgnams(nk,korbs,kspin) .eq. 0 ) then
               cvec(ni) = cvec(ni) + zmtmp(mi,rowams(nk,korbs,kspin))
            else
               cvec(ni) = cvec(ni) - zmtmp(mi,rowams(nk,korbs,kspin))
            end if
         end do
      end do
   else 
      do ni=1,nnz_c
         mi = irow_c(ni)
         mj = icol_c(ni)
         do nk=1,nams
            if ( rowams(nk,korbs,kspin) .ne. mj ) cycle
            if ( sgnams(nk,korbs,kspin) .eq. 0 ) then
               cvec(ni) = cvec(ni) + zmtmp(mi,colams(nk,korbs,kspin))
            else
               cvec(ni) = cvec(ni) - zmtmp(mi,colams(nk,korbs,kspin))
            end if
         end do
      end do
   end if
end if
!
do ni=1,nnz_c
   c(ni) = zbeta * c(ni) + zalpha * cvec(ni)
end do
!
end subroutine zamsmm2_coo
