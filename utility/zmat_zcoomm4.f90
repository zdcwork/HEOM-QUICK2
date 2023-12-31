subroutine zmat_zcoomm4(lefta, transa, transb, alpha, a, nnza, irowa, icola,     &
                        b, nnzb, irowb, icolb, beta, c, ldc, m, n, cval1, cval2)
implicit none
!
! alpha * op(A) * op(B) + beta * op(C) ==> C
! A, B are in sparse COO format, and C is a full matrix
!
! IMPORTANT: 
!    For COO storage, nonzeros with smaller column index appear earlier in the array
!
! Note: 
! (lefta, transa, transb) = (L, N, X) or (R, X, N) is faster,
! i.e., when the left matrix does not undergo transpose 
!
character*1, intent(in)    :: lefta, transa, transb
integer,     intent(in)    :: nnza, nnzb, irowa(*), icola(*), irowb(*), icolb(*)
integer,     intent(in)    :: m, n, ldc
complex*16,  intent(in)    :: a(*), b(*), alpha, beta
complex*16,  intent(inout) :: c(ldc,*), cval1(*), cval2(*)
logical                    :: lca, lcb
logical                    :: lsame
external                   :: lsame
integer                    :: ni, nj
integer                    :: mi, mj, mk
complex*16                 :: ctmp1
!
if (nnza .lt. 0 .or. nnzb .lt. 0 .or. ldc .le. 0) then
    write(6,*)
    write(6,*)'zmat_zcoomm4: error! wrong nnz found ', nnza, nnzb, ldc
    stop
end if
!
c(1:m,1:n) = c(1:m,1:n) * beta
if (nnza .eq. 0 .or. nnzb .eq. 0) return
!
if (lsame(transa, 'N') .or. lsame(transa, 'R')) then       ! no need to transpose
    lca = .false.  
else if (lsame(transa, 'C') .or. lsame(transa, 'T')) then  ! need to transpose
    lca = .true.
else
    write(6,*)
    write(6,*)'zmat_zcoomm4: error! unknown transa ', transa
    stop
end if
if (lsame(transb, 'N') .or. lsame(transb, 'R')) then       ! no need to transpose
    lcb = .false.  
else if (lsame(transb, 'C') .or. lsame(transb, 'T')) then  ! need to transpose
    lcb = .true.
else
    write(6,*)
    write(6,*)'zmat_zcoomm4: error! unknown transb ', transb
    stop
end if
!
if (lsame(transa, 'N') .or. lsame(transa, 'T')) then       ! no conjugation 
    cval1(1:nnza) = a(1:nnza)
else                                                       ! conjugation
    cval1(1:nnza) = dconjg(a(1:nnza))
end if
if (lsame(transb, 'N') .or. lsame(transb, 'T')) then       ! no conjugation
    cval2(1:nnzb) = b(1:nnzb)
else                                                       ! conjugation
    cval2(1:nnzb) = dconjg(b(1:nnzb))
end if
!
if ( lsame(lefta, 'L') ) then
    if (.not. lca .and. .not. lcb) then
        do ni=1,nnzb
           mj = icolb(ni)
           mk = irowb(ni)
           ctmp1 = alpha * cval2(ni)
           loop1: do nj=1,nnza
              if (icola(nj) .gt. mk) exit loop1
              if (icola(nj) .ne. mk) cycle
              mi = irowa(nj)
              c(mi,mj) = c(mi,mj) + cval1(nj) * ctmp1
           end do loop1
        end do
        !
    else if (.not. lca .and. lcb) then
        do ni=1,nnzb
           mj = irowb(ni)
           mk = icolb(ni)
           ctmp1 = alpha * cval2(ni)
           loop2: do nj=1,nnza
              if (icola(nj) .gt. mk) exit loop2
              if (icola(nj) .ne. mk) cycle
              mi = irowa(nj)
              c(mi,mj) = c(mi,mj) + cval1(nj) * ctmp1
           end do loop2
        end do
        !
    else if (lca .and. .not. lcb) then
        do ni=1,nnzb
           mj = icolb(ni)
           mk = irowb(ni)
           ctmp1 = alpha * cval2(ni)
           loop3: do nj=1,nnza
              if (irowa(nj) .ne. mk) cycle
              mi = icola(nj)
              c(mi,mj) = c(mi,mj) + cval1(nj) * ctmp1
           end do loop3
        end do
        !
    else  ! lca .and. lcb 
        do ni=1,nnzb
           mj = irowb(ni)
           mk = icolb(ni)
           ctmp1 = alpha * cval2(ni)
           loop4: do nj=1,nnza
              if (irowa(nj) .ne. mk) cycle
              mi = icola(nj)
              c(mi,mj) = c(mi,mj) + cval1(nj) * ctmp1
           end do loop4
        end do
        !
    endif
!
else if ( lsame(lefta, 'R') ) then
    if (.not. lca .and. .not. lcb) then
        do ni=1,nnza
           mj = icola(ni)
           mk = irowa(ni)
           ctmp1 = alpha * cval1(ni)
           loop5: do nj=1,nnzb
              if (icolb(nj) .gt. mk) exit loop5
              if (icolb(nj) .ne. mk) cycle
              mi = irowb(nj)
              c(mi,mj) = c(mi,mj) + cval2(nj) * ctmp1
           end do loop5
        end do
        !
    else if (lca .and. .not. lcb) then
        do ni=1,nnza
           mj = irowa(ni)
           mk = icola(ni)
           ctmp1 = alpha * cval1(ni)
           loop6: do nj=1,nnzb
              if (icolb(nj) .gt. mk) exit loop6
              if (icolb(nj) .ne. mk) cycle
              mi = irowb(nj)
              c(mi,mj) = c(mi,mj) + cval2(nj) * ctmp1
           end do loop6
        end do
        !
    else if (.not. lca .and. lcb) then
        do ni=1,nnza
           mj = icola(ni)
           mk = irowa(ni)
           ctmp1 = alpha * cval1(ni)
           loop7: do nj=1,nnzb
              if (irowb(nj) .ne. mk) cycle
              mi = icolb(nj)
              c(mi,mj) = c(mi,mj) + cval2(nj) * ctmp1
           end do loop7
        end do
        !
    else ! lca .and. lcb
        do ni=1,nnza
           mj = irowa(ni)
           mk = icola(ni)
           ctmp1 = alpha * cval1(ni)
           loop8: do nj=1,nnzb
              if (irowb(nj) .ne. mk) cycle
              mi = icolb(nj)
              c(mi,mj) = c(mi,mj) + cval2(nj) * ctmp1
           end do loop8
        end do
        !
    end if
!
else
   write(6,*)
   write(6,*)'zmat_zcoomm4: error, unknown lefta ', lefta
   stop
end if
!
return
end subroutine zmat_zcoomm4
