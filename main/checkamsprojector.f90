subroutine checkamsprojector
use matmod
use sparsemod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, ntmp1
integer :: iorb, ispin, istat
!
if (.not. (lproj .and. lprgen)) then
    write(6,*)
    write(6,*)'checkamsprojector: error entry! '
    write(6,*)lproj, lprgen
    stop
end if
!
zamsall = czero
allocate(nnz_ams(norbs,nspin), STAT=istat)
!
ntmp1 = nams
do ispin=1,nspin
   do iorb=1,norbs
      cmtmp1(1:nrho,1:nrho) = amsall(1:nrho,1:nrho,iorb,ispin) * cunity
      call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, &
                 czero, cmtmp3, nrho)
      call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp3, nrho, &
                 czero, cmtmp2, nrho)
      zamsall(1:nrho,1:nrho,iorb,ispin) = cmtmp2(1:nrho,1:nrho)
      nk = 0
      do nj=1,nrho
         do ni=1,nrho
            !if (cdabs(cmtmp2(ni,nj)) .gt. dsmall) then
            if (cdabs(cmtmp2(ni,nj)) .gt. dpico) then
                nk = nk + 1
            end if
         end do
      end do
      nnz_ams(iorb,ispin) = nk
      ntmp1 = max(ntmp1, nk)
   end do
end do
!
allocate( row_ams(ntmp1,norbs,nspin), col_ams(ntmp1,norbs,nspin), &
         cval_ams(ntmp1,norbs,nspin), STAT=istat)
!
write(6,*)
write(6,*)'checkamsprojector: check sparsity of a_ms * cproj '
write(6,*)' nams = ', nams
write(6,*)' iorb, ispin, nonzeros '
call flush(6)
do ispin=1,nspin
   do iorb=1,norbs
      cmtmp1(1:nrho,1:nrho) = amsall(1:nrho,1:nrho,iorb,ispin) * cunity
      call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, &
                 czero, cmtmp3, nrho)
      call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp3, nrho, &
                 czero, cmtmp2, nrho)
!
! update:    rowams, colams, nnz_ams, cval_ams
! unchanged: sgnams, valams 
! note (2023-04-24): intuitively, nnz_ams should not exceed nams
!
      nk = 0
      do nj=1,nrho
         do ni=1,nrho
            !if (cdabs(cmtmp2(ni,nj)) .gt. dsmall) then
            if (cdabs(cmtmp2(ni,nj)) .gt. dpico) then
                nk = nk + 1
!                if (nk .gt. nams) cycle
                cval_ams(nk,iorb,ispin) = cmtmp2(ni,nj)
                row_ams (nk,iorb,ispin) = ni
                col_ams (nk,iorb,ispin) = nj
            end if
         end do
      end do
      write(6,1000)iorb, ispin, nk
      call flush(6)
   end do
end do
1000 format(2(1x, I4), 2x, I6)
!
!lprams = .true.
write(6,*)'checkamsprojector: annihilation operators projected ' 
if (lprams) then
    write(6,*)'checkamsprojector: utilize projected a_ms '
else
    write(6,*)'checkamsprojector: do not utilize projected a_ms '
end if
call flush(6)
!
! apply projector to occupancy operators
!
if (lprams) then
    do ispin=1,nspin
       do iorb=1,norbs
          cmtmp1(1:nrho,1:nrho) = znms(1:nrho,1:nrho,iorb,ispin)
          call zhemm('r', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp1, nrho, &
                     czero, cmtmp2, nrho)
          call zhemm('l', 'u', nrho, nrho, cunity, cproj, nrho, cmtmp2, nrho, &
                     czero, cmtmp1, nrho)
          znms(1:nrho,1:nrho,iorb,ispin) = cmtmp1(1:nrho,1:nrho)
       end do
    end do
    write(6,*)'checkamsprojector: utilize projected n_ms '
    call flush(6)
end if
!
if (ntmp1 .gt. nams) then
    write(6,*)'checkamsprojector: nnz_ams > nams found ', ntmp1, nams
    write(6,*)'checkamsprojector: Warning! projection on a_ms is ineffective '
    call flush(6)
end if
!
return
end subroutine checkamsprojector
