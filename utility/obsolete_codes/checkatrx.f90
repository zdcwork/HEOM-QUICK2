subroutine checkatrx
use matmod
use bicgmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer   :: ni, nj, nk, nl, istat, info, ntmp1, ntmp2, n1, n2
integer, allocatable :: ipiv(:), irow0(:)
integer*8 :: lni, lnj, lnk, lnl, ltmp1, ltmp2
real*8    :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, vl(1)
real*8, allocatable :: tmp1cg(:), tmp2cg(:), tmp3cg(:)
real*8, allocatable :: a(:,:), at(:,:), tmpvec(:,:), valr(:), vali(:), b(:,:), bt(:,:)
real*8     :: ddotlong, cdotlong
external   :: ddotlong, cdotlong
!
ltmp1 = nunk * nrho**2
ltmp2 = ltmp1 * 2
!
!call checksymm
!stop
!
allocate(tmp1cg(ltmp2), tmp2cg(ltmp2), tmp3cg(ltmp2), STAT=istat)
allocate(a(ltmp2,ltmp2), at(ltmp2,ltmp2), STAT=istat)
if (istat .ne. 0) then
  write(6,*)' allocation fail in checkatrx 1', istat
  stop
end if
!
do lni=1,ltmp2
  tmp1cg(1:ltmp2) = 0.d0
  tmp1cg(lni)     = 1.d0
  call matvec(rk, nrho, nunk, tmp1cg, ltmp2, 1)
  call calc_ax_bicg(nrho, rk, rcgrhs)
  call matvec(rcgrhs, nrho, nunk, tmp2cg, ltmp2, 0)
  a(1:ltmp2, lni) = tmp2cg(1:ltmp2)
  call calc_atrx_bicg(nrho, rk, rcgrhs)  
  call matvec(rcgrhs, nrho, nunk, tmp2cg, ltmp2, 0)
  at(1:ltmp2, lni) = tmp2cg(1:ltmp2)
end do 
!
dtmp3 = 0.d0
lnk   = 0
do lni=1,ltmp2
  do lnj=1,ltmp2
    if (dtmp3 .le. dabs(a(lni,lnj) - at(lnj,lni))) then
      lnk = lni
      lnl = lnj
      dtmp3 = dabs(a(lni,lnj) - at(lnj,lni))
    end if
  end do
end do
write(6,*)
write(6,*)' maximal diff (a-at) ', dtmp3
if (.not. (lnk .eq. ltmp2 .and. lnl .eq. ltmp2)) then
  write(6,*)' irow = ', lnk, ' icol = ', lnl
  write(6,*)a(lnk,lnl), at(lnl,lnk)
end if
!
sk(1:nrho, 1:nrho, 1:nunk) = rho(1:nrho, 1:nrho, 1:nunk)
nk = 0
do ni=1,nrho
  do nj=ni+1,nrho
    nk = nk + 1
    sk(ni,nj,1) = czero
  end do
end do
call matvec(sk, nrho, nunk, tmp1cg, ltmp2, 0)
!
do lni=1,ltmp2
  tmp2cg(lni) = 0.d0
  tmp3cg(lni) = 0.d0
  do lnj=1,ltmp2
    tmp2cg(lni) = tmp2cg(lni) +  a(lni,lnj) * tmp1cg(lnj)
    tmp3cg(lni) = tmp3cg(lni) + at(lni,lnj) * tmp1cg(lnj)
  end do
end do
!
call matvec(rk, nrho, nunk, tmp2cg, ltmp2, 1)
call matvec(dk, nrho, nunk, tmp3cg, ltmp2, 1)
!
dtmp1 = cdotlong(ltmp1, sk(1,1,1), rk(1,1,1))
dtmp2 = cdotlong(ltmp1, sk(1,1,1), dk(1,1,1))
!
write(6,*)
write(6,*)' rho^T * A   * rho = ', dtmp1
write(6,*)' rho^T * A^T * rho = ', dtmp2
!
dtmp4 = 0.d0
dtmp5 = 0.d0
do lni=1,ltmp2
  do lnj=1,ltmp2
    dtmp4 = dtmp4 +  a(lni,lnj)**2
    dtmp5 = dtmp5 + at(lni,lnj)**2
  end do
end do
write(6,*)
write(6,*)' sum(a, at) = ', dtmp4, dtmp5
call flush(6)
!
stop
!
! check if any variable is irrelevant (rows of at)
!
write(6,*)
write(6,*)' check rows of at'
ni = 0
do lni=1,ltmp2
  dtmp1 = 0.d0
  do lnj=1,ltmp2
    dtmp1 = dtmp1 + dabs(at(lni,lnj))**2
  end do
  write(6,*)lni, dtmp1
  if (dtmp1 .le. dnano) ni = ni + 1
end do
write(6,*)' # of zeros ', ni
call flush(6)
!
! check if any variable is irrelevant (columns of at)
!
write(6,*)
write(6,*)' check columns of at'
ni = 0
do lni=1,ltmp2
  dtmp1 = 0.d0
  do lnj=1,ltmp2
    dtmp1 = dtmp1 + dabs(at(lnj,lni))**2
  end do
  write(6,*)lni, dtmp1
  if (dtmp1 .le. dnano) ni = ni + 1
end do
write(6,*)' # of zeros ', ni
call flush(6)
!
! try to invert matrix a (the correct a should be obtained by transpose(at))
!
ntmp2 = 2 * nunk * nrho**2
if (dabs( dble(ntmp2) - dble(ltmp2) ) .ge.  1.d-9) stop
ntmp1 = ntmp2
!
ntmp2 = ntmp1 - nrho * (nrho - 1)   ! consider hermicity of rho(*,*,1) (off-diagonal elements)
ntmp2 = ntmp2 - nrho                ! consider hermicity of rho(*,*,1) (    diagonal elements, imaginary part)
!
allocate(ipiv(ntmp2), irow0(nrho*(nrho-1)), STAT=istat)
allocate(tmpvec(ntmp2,ntmp2), valr(ntmp2), vali(ntmp2), STAT=istat)
allocate(b(ntmp2,ntmp2), bt(ntmp2,ntmp2), STAT=istat)
!
!goto 60
!
! remove rows and columns from "at" in the aspect of hermicity of rho(*,*,1)
!
a(1:ntmp1, 1:ntmp1) = at(1:ntmp1, 1:ntmp1)
!     hermicity  (off-diagonal elements)
nk = 0
do ni=1,nrho
  do nj=1,nrho
    if (nj .gt. ni) then
      nk = nk + 1
      irow0(nk) = (nj - 1) * nrho + ni
      nk = nk + 1
      irow0(nk) = irow0(nk-1) + ntmp1 / 2
    end if
  end do
end do
!     hermicity  (diagonal elements)
do ni=1,nrho
  nk = nk + 1
  irow0(nk) = (ni - 1) * nrho + ni + ntmp1 / 2
end do
!
write(6,*)
write(6,*)' # of rows to remove : ', nk 
!
do ni=1,nk-1
  do nj=ni+1,nk
    if (irow0(ni) .lt. irow0(nj)) then
      n1 = irow0(ni)
      irow0(ni) = irow0(nj)
      irow0(nj) = n1
    end if
  end do
end do
write(6,*)(irow0(ni), ni=1,nk)
!
do ni=1,nk
  do nj=irow0(ni)+1,ntmp1
    at(nj-1,1:ntmp1) = at(nj,1:ntmp1)
  end do
end do
!
do ni=1,nk
  do nj=irow0(ni)+1,ntmp1 
    at(1:ntmp1,nj-1) = at(1:ntmp1,nj)
  end do
end do
!
60 continue
!
bt(1:ntmp2,1:ntmp2) = at(1:ntmp2,1:ntmp2)
!
do ni=1,ntmp2
  do nj=1,ntmp2
    b(ni,nj) = bt(nj,ni)
  end do
end do
!
call dgeev('n', 'v', ntmp2, b, ntmp2, valr, vali, vl, 1, tmpvec, ntmp2, bt, &
           ntmp2**2, info)
write(6,*)
write(6,*)' eigenvalues of a '
nl = 0
do ni=1,ntmp2
  write(6,*)ni, valr(ni), vali(ni)
  if ( valr(ni)**2 + vali(ni)**2 .le. 1.d-9) nl = nl + 1
end do
write(6,*)' # of zero eigenvalues ', nl
call flush(6)
!
stop
!
bt(1:ntmp2, 1:ntmp2) = b(1:ntmp2, 1:ntmp2)
call dgetrf(ntmp2, ntmp2, bt, ntmp2, ipiv, info)
call dgetri(ntmp2, bt, ntmp2, ipiv, b, ntmp2**2, info)
!
write(6,*)
write(6,*)' output flag from inversion of a ', info
call flush(6)
!
deallocate(ipiv, irow0, b, bt, STAT=istat)
deallocate(tmpvec, valr, vali, STAT=istat)
!
deallocate(a, at, STAT=istat)
deallocate(tmp1cg, tmp2cg, tmp3cg, STAT=istat)
!
end subroutine checkatrx
!
subroutine matvec(cmat, nrho, nunk, dvec, len0, iop)
implicit none
!
integer*8, intent(in) :: nunk, len0
integer,   intent(in) :: nrho, iop
complex*16, intent(inout) :: cmat(nrho,nrho,*)
real*8,     intent(inout) :: dvec(*)
!
integer :: ni, nj, nk
integer*8 :: ltmp1, lni, lnj
!
ltmp1 = len0 / 2
do lni=1,nunk
  lnj = (lni - 1) * nrho**2
  nk  = 0
  do nj=1,nrho
    do ni=1,nrho
      nk = nk + 1
      if (iop .eq. 0) then   ! mat2vec
        dvec(nk+lnj)       =  dble(cmat(ni,nj,lni))
        dvec(nk+lnj+ltmp1) = dimag(cmat(ni,nj,lni))
      else                   ! vec2mat
        cmat(ni,nj,lni) = dcmplx(dvec(nk+lnj), dvec(nk+lnj+ltmp1))
      end if
    end do 
  end do
end do
end subroutine matvec
