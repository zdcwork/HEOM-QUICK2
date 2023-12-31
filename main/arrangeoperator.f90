subroutine arrangeoperator
use matmod
use tmpmatmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, istat
integer :: ispin, iorbs 
real*8  :: dtmp0
!
integer :: ijob(6)
integer :: info
!
if (nspin .eq. 1) then
  nams = 2**(norbs - 1)
else if (nspin .eq. 2) then
  nams = 2 * 4**(norbs - 1)
end if
write(6,*)
write(6,*)'arrangeoperator: # of nonzero element in c_ms : ', nams
call flush(6)
!
allocate(rowams(nams,norbs,nspin), STAT=istat)
allocate(colams(nams,norbs,nspin), STAT=istat)
allocate(sgnams(nams,norbs,nspin), STAT=istat)
allocate(valams(nams,norbs,nspin), STAT=istat)
!
do ispin=1,nspin
   do iorbs=1,norbs
      nk = 0
      do nj=1,nrho
         do ni=1,nrho
            dtmp0 = amsall(ni,nj,iorbs,ispin)
            if (dabs(dtmp0) .ge. dpico) then
               nk = nk + 1
               rowams(nk,iorbs,ispin) = ni
               colams(nk,iorbs,ispin) = nj
               valams(nk,iorbs,ispin) = dtmp0
               if (dtmp0 .gt. 0.d0) then
                  sgnams(nk,iorbs,ispin) = 0
               else 
                  sgnams(nk,iorbs,ispin) = 1
               end if
            end if
         end do
      end do
      if (nk .ne. nams) then
        write(6,*)
        write(6,*)'arrangeoperator: error when analyzing amsall ', iorbs, ispin
        stop
      end if
   end do
end do
!
if (.not. lsparse) return   
!
write(6,*)
write(6,*)'arrangeoperator: SPARSE mode is turned on '
nnzero_spa = nams
!
allocate(ams_spa(nnzero_spa,norbs,nspin), STAT=istat)
allocate(ja_spa(nnzero_spa,norbs,nspin), STAT=istat)
allocate(ia_spa(nrho+1,norbs,nspin), STAT=istat)
!
ijob(1) = 0   ! convert dense matrix to CSR format
ijob(2) = 1   ! one-based index for dense matrix
ijob(3) = 1   ! one-based index for sparse matrix
ijob(4) = 2   ! the whole dense matrix is used
ijob(5) = nnzero_spa
ijob(6) = 1   ! generate ams_spa, ia_spa, and ja_spa
!
do ispin=1,nspin
   do iorbs=1,norbs
      call mkl_ddnscsr(ijob, nrho, nrho, amsall(1,1,iorbs,ispin), nrho, ams_spa(1,iorbs,ispin), &
                       ja_spa(1,iorbs,ispin), ia_spa(1,iorbs,ispin), info)
      if (info .ne. 0) then
         write(6,*)'arrangeoperator: error compressing ams matrix ', iorbs, ispin
         stop
      end if
   end do
end do
!
write(6,*)
write(6,*)'arrangeoperator: amsall successfully converted to CSR format '
call flush(6)
!
end subroutine arrangeoperator
