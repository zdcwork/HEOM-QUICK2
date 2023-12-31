subroutine zmat_zcoomm2(lefta, transa, transb, n, alpha, a, nnza, irowa, icola,  &
                        b, nnzb, irowb, icolb, beta, c, nnzc, irowc, icolc,      &
                        cval1, irow1, icol1, cval2, irow2, icol2,                &
                        zmtmp, ldtmp)
implicit none
!
! alpha * op(A) * op(B) + beta * op(C) ==> C
! A, B, C are all in sparse COO format
!
! IMPORTANT: 
!    For COO storage, nonzeros with smaller column index appear earlier in the array
!
! Note added on Feb 20, 2020: 
! Cautious! It seems that this code could be wrong if transa='c' or transb='c', because 
! in such cases, after call of zmat_a2b_coo, nonzero elements do not follow 
! the rule that "nonzeros with smaller column index appear earlier in the array" 
!
!
character*1, intent(in)    :: lefta, transa, transb
integer,     intent(in)    :: n, nnza, nnzb, nnzc, irowa(*), icola(*), irowb(*), icolb(*), irowc(*), icolc(*)
integer,     intent(in)    :: ldtmp
complex*16,  intent(in)    :: a(*), b(*), alpha, beta
integer,     intent(inout) :: irow1(*), icol1(*), irow2(*), icol2(*)
complex*16,  intent(inout) :: c(*), cval1(*), cval2(*), zmtmp(ldtmp,*)
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj
integer                    :: mi, mj, mk
!
if (nnzc .lt. 0 .or. nnzb .lt. 0 .or. nnza .lt. 0) then
   write(6,*)
   write(6,*)'zmat_zcoomm2: error! wrong nnz found ', nnza, nnzb, nnzc
   stop
end if
if (nnzc .eq. 0) return
!
! early stop if transa='T' or 'C' or transb='T' or 'C' (added on Feb 20, 2020)
! the code may be wrong in such cases, see the note above
! 
if ( lsame(transa, 'T') .or. lsame(transa, 'C') .or.   &
     lsame(transb, 'T') .or. lsame(transb, 'C') ) then
    write(6,*)
    write(6,*)'zmat_zcoomm2: need to transpose matrix ', transa, transb
    write(6,*)'              code may not be reliable '
    stop
end if
!
c(1:nnzc) = c(1:nnzc) * beta
if (nnza .eq. 0 .or. nnzb .eq. 0) return
!
call zmat_a2b_coo(transa, nnza, irowa, icola, a, irow1, icol1, cval1)
call zmat_a2b_coo(transb, nnzb, irowb, icolb, b, irow2, icol2, cval2)
!
zmtmp(1:n,1:n) = dcmplx(0.d0, 0.d0)
!
if ( lsame(lefta, 'L') ) then
   !do ni=1,nnza
   !   mi = irow1(ni)
   !   mk = icol1(ni)
   !   do nj=1,nnzb
   !      if (irow2(nj) .ne. mk) cycle
   !      mj = icol2(nj)
   !      zmtmp(mi,mj) = zmtmp(mi,mj) + cval1(ni) * cval2(nj)
   !   end do
   !end do
   do ni=1,nnzb
      mj = icol2(ni)
      mk = irow2(ni)
      loop1: do nj=1,nnza
         if (icol1(nj) .gt. mk) exit loop1
         if (icol1(nj) .ne. mk) cycle
         mi = irow1(nj)
         zmtmp(mi,mj) = zmtmp(mi,mj) + cval1(nj) * cval2(ni)
      end do loop1
   end do
!
else if ( lsame(lefta, 'R') ) then
   !do ni=1,nnzb
   !   mi = irow2(ni)
   !   mk = icol2(ni)
   !   do nj=1,nnza
   !      if (irow1(nj) .ne. mk) cycle
   !      mj = icol1(nj)
   !      zmtmp(mi,mj) = zmtmp(mi,mj) + cval2(ni) * cval1(nj)
   !   end do
   !end do
   do ni=1,nnza
      mj = icol1(ni)
      mk = irow1(ni)
      loop2: do nj=1,nnzb
         if (icol2(nj) .gt. mk) exit loop2
         if (icol2(nj) .ne. mk) cycle
         mi = irow2(nj)
         zmtmp(mi,mj) = zmtmp(mi,mj) + cval2(nj) * cval1(ni)
      end do loop2
   end do
!
else
   write(6,*)
   write(6,*)'zmat_zcoomm2: error, unknown lefta ', lefta
   stop
end if
!
call zmat_extract_coo(nnzc, irowc, icolc, cval1, zmtmp, n, n, ldtmp)
call zmat_add_coo('n', nnzc, alpha, cval1, dcmplx(1.d0, 0.d0), c)
!
return
end subroutine zmat_zcoomm2
