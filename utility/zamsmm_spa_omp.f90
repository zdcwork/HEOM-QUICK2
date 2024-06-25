subroutine zamsmm_spa_omp(sidea, transa, transb, korbs, kspin, zalpha, b, ldb, zbeta, c, ldc, &
                          dmtmp1, ldtmp1, dmtmp2, ldtmp2, dmtmp3, ldtmp3, dmtmp4, ldtmp4)
!
! use sparse BLAS to calculate ams*ADO 
!
use matmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
character*1, intent(in)   :: sidea, transa, transb
integer, intent(in)       :: korbs, kspin, ldb, ldc, ldtmp1, ldtmp2, ldtmp3, ldtmp4
complex*16, intent(in)    :: zalpha, zbeta, b(ldb,*)
complex*16, intent(inout) :: c(ldc,*)
real*8, intent(inout)     :: dmtmp1(ldtmp1,*), dmtmp2(ldtmp2,*), dmtmp3(ldtmp3,*), dmtmp4(ldtmp4,*)
!
integer     :: ni, nj, nk
integer     :: mi, mj, mk
logical     :: lsame
external    :: lsame
logical     :: nota, notb, conja, conjb, lefta, righta
integer     :: icase
character*1 :: trans_ams
character   :: matdescra(6)
real*8      :: dtmp1, dtmp2
!
if (ldb .ne. nrho .or. ldc .ne. nrho) then
   write(6,*)'zamsmm_spa_omp: error matrix dimension'
   write(6,*)'zamsmm_spa_omp: nrho = ', nrho, ' ldb = ', ldb, ' ldc = ', ldc
   stop
end if
!
lefta  = lsame(sidea,  'L')
righta = lsame(sidea,  'R')
nota   = lsame(transa, 'N')
conja  = lsame(transa, 'C')
notb   = lsame(transb, 'N')
conjb  = lsame(transb, 'C')
!
icase = -1
if (lefta  .and. notb ) icase = 1
if (lefta  .and. conjb) icase = 2
if (righta .and. notb ) icase = 3
if (righta .and. conjb) icase = 4
!
select case (icase)
case (1)
! dmtmp1 == alpha_r * b_r - alpha_i * b_i
! dmtmp2 == alpha_r * b_i + alpha_i * b_r
     dtmp1 = dble (zalpha)
     dtmp2 = dimag(zalpha)
     do nj=1,nrho
        do ni=1,nrho
           dmtmp1(ni,nj) = dtmp1 * dble (b(ni,nj)) - dtmp2 * dimag(b(ni,nj))
           dmtmp2(ni,nj) = dtmp1 * dimag(b(ni,nj)) + dtmp2 * dble (b(ni,nj))
        end do
     end do
     if (nota) then
        trans_ams = 'n'
     else
        trans_ams = 't'
     end if
case (2)
! dmtmp1 ==  alpha_r * b_r^t + alpha_i * b_i^t
! dmtmp2 == -alpha_r * b_i^t + alpha_i * b_r^t
     dtmp1 = dble (zalpha)
     dtmp2 = dimag(zalpha)
     do nj=1,nrho
        do ni=1,nrho
           dmtmp1(ni,nj) =  dtmp1 * dble (b(nj,ni)) + dtmp2 * dimag(b(nj,ni))
           dmtmp2(ni,nj) = -dtmp1 * dimag(b(nj,ni)) + dtmp2 * dble (b(nj,ni))
        end do
     end do
     if (nota) then
        trans_ams = 'n'
     else
        trans_ams = 't'
     end if
case (3)
! dmtmp1 ==  alpha_r * b_r^t - alpha_i * b_i^t
! dmtmp2 == -alpha_r * b_i^t - alpha_i * b_r^t
     dtmp1 = dble (zalpha)
     dtmp2 = dimag(zalpha)
     do nj=1,nrho
        do ni=1,nrho
           dmtmp1(ni,nj) =  dtmp1 * dble (b(nj,ni)) - dtmp2 * dimag(b(nj,ni))
           dmtmp2(ni,nj) = -dtmp1 * dimag(b(nj,ni)) - dtmp2 * dble (b(nj,ni))
        end do
     end do
     if (nota) then
        trans_ams = 't'
     else
        trans_ams = 'n'
     end if
case (4)
! dmtmp1 == alpha_r * b_r + alpha_i * b_i
! dmtmp2 == alpha_r * b_i - alpha_i * b_r
     dtmp1 = dble (zalpha)
     dtmp2 = dimag(zalpha)
     do nj=1,nrho
        do ni=1,nrho
           dmtmp1(ni,nj) = dtmp1 * dble (b(ni,nj)) + dtmp2 * dimag(b(ni,nj))
           dmtmp2(ni,nj) = dtmp1 * dimag(b(ni,nj)) - dtmp2 * dble (b(ni,nj))
        end do
     end do
     if (nota) then
        trans_ams = 't'
     else
        trans_ams = 'n'
     end if
case default
     write(6,*)
     write(6,*)'zamsmm_spa_omp: wrong case found ', icase
end select
!
! matdescra
matdescra(1) = 'g'
matdescra(4) = 'f'
!
! real part
call mkl_dcsrmm(trans_ams, nrho, nrho, nrho, 1.d0, matdescra, ams_spa(1,korbs,kspin),  &
                ja_spa(1,korbs,kspin), ia_spa(1,korbs,kspin), ia_spa(2,korbs,kspin),   &
                dmtmp1, nrho, 0.d0, dmtmp3, nrho)
! imaginary part
call mkl_dcsrmm(trans_ams, nrho, nrho, nrho, 1.d0, matdescra, ams_spa(1,korbs,kspin),  &
                ja_spa(1,korbs,kspin), ia_spa(1,korbs,kspin), ia_spa(2,korbs,kspin),   &
                dmtmp2, nrho, 0.d0, dmtmp4, nrho)
!
! beta * C 
!
! dmtmp1 == beta_r  * c_r - beta_i  * c_i
! dmtmp2 == beta_r  * c_i + beta_i  * c_r
!if (cdabs(zbeta) .lt. dpico) then 
if (cdabs(zbeta) .lt. dsmall) then 
   dmtmp1(1:nrho,1:nrho) = 0.d0
   dmtmp2(1:nrho,1:nrho) = 0.d0
else
   dtmp1 = dble (zbeta)
   dtmp2 = dimag(zbeta)
   do nj=1,nrho
      do ni=1,nrho
         dmtmp1(ni,nj) = dtmp1 * dble (c(ni,nj)) - dtmp2 * dimag(c(ni,nj))
         dmtmp2(ni,nj) = dtmp1 * dimag(c(ni,nj)) + dtmp2 * dble (c(ni,nj))
      end do
   end do
end if
!
if (lefta) then
   do nj=1,nrho
      do ni=1,nrho
         c(ni,nj) = dcmplx(dmtmp3(ni,nj), dmtmp4(ni,nj)) + dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
      end do
   end do
else
   do nj=1,nrho
      do ni=1,nrho
         c(ni,nj) = dcmplx(dmtmp3(nj,ni), -dmtmp4(nj,ni)) + dcmplx(dmtmp1(ni,nj), dmtmp2(ni,nj))
      end do
   end do
end if
!
end subroutine zamsmm_spa_omp
