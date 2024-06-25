subroutine zmat_zcoomm1(lefta, transa, transb, n, alpha, a, nnza, irowa, icola,  &
                        b, ldb, beta, c, nnzc, irowc, icolc,                     &
                        cvaltmp, irowtmp, icoltmp)
implicit none
! 
! alpha * op(A) * op(B) + beta * op(C) ==> C
!
! A and C are in sparse COO format, B is dense matrix
!
character*1, intent(in)    :: lefta, transa, transb
integer,     intent(in)    :: n, nnza, ldb, nnzc, irowa(*), icola(*), irowc(*), icolc(*)
complex*16,  intent(in)    :: a(*), b(ldb,*), alpha, beta
integer,     intent(inout) :: irowtmp(*), icoltmp(*)
complex*16,  intent(inout) :: c(*), cvaltmp(*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj
integer                    :: mi, mj, mk
!
if (nnzc .lt. 0 .or. nnza .lt. 0) then
   write(6,*)
   write(6,*)'zmat_zcoomm1: error! wrong nnz found ', nnza, nnzc
   stop
end if
if (nnzc .eq. 0) return
!
c(1:nnzc) = c(1:nnzc) * beta
if (nnza .eq. 0) return
!
call zmat_a2b_coo(transa, nnza, irowa, icola, a, irowtmp, icoltmp, cvaltmp)
cvaltmp(1:nnza) = alpha * cvaltmp(1:nnza)
!
if ( lsame(lefta, 'L') ) then
   if ( lsame(transb, 'N') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (irowtmp(nj) .ne. mi) cycle
            mk = icoltmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * b(mk, mj)
         end do
      end do 
   else if ( lsame(transb, 'T') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (irowtmp(nj) .ne. mi) cycle
            mk = icoltmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * b(mj, mk)
         end do
      end do 
   else if ( lsame(transb, 'R') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (irowtmp(nj) .ne. mi) cycle
            mk = icoltmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * dconjg(b(mk, mj))
         end do
      end do
   else if ( lsame(transb, 'C') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (irowtmp(nj) .ne. mi) cycle
            mk = icoltmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * dconjg(b(mj, mk))
         end do
      end do
   else
      write(6,*)
      write(6,*)'zmat_zcoomm1: error, unknown transb ', transb
      stop
   end if
!
else if ( lsame(lefta, 'R') ) then
   if ( lsame(transb, 'N') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (icoltmp(nj) .ne. mj) cycle
            mk = irowtmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * b(mi, mk)
         end do
      end do 
   else if ( lsame(transb, 'R') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (icoltmp(nj) .ne. mj) cycle
            mk = irowtmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * dconjg(b(mi, mk))
         end do
      end do 
   else if ( lsame(transb, 'T') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (icoltmp(nj) .ne. mj) cycle
            mk = irowtmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * b(mk, mi)
         end do
      end do 
   else if ( lsame(transb, 'C') ) then
      do ni=1,nnzc
         mi = irowc(ni)
         mj = icolc(ni)
         do nj=1,nnza
            if (icoltmp(nj) .ne. mj) cycle
            mk = irowtmp(nj)
            c(ni) = c(ni) + cvaltmp(nj) * dconjg(b(mk, mi))
         end do
      end do 
   else
      write(6,*)
      write(6,*)'zmat_zcoomm1: error, unknown transb ', transb
      stop
   end if
!
else
   write(6,*)
   write(6,*)'zmat_zcoomm1: error, unknown lefta ', lefta
   stop
end if
!
return
end subroutine zmat_zcoomm1
