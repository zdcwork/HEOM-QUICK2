subroutine testdiis
implicit none
!
integer :: ni, nj, istat, info, i, j, ierror
integer :: norbs, lwork
integer, allocatable :: ipiv(:)
real*8              :: diff, summ, cpu1, cpu2, cpu3, cpu4
real*8, allocatable :: dtmp1(:,:), dtmp2(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:), xdiis(:), ctmp3(:), ctmp4(:,:)
complex*16, allocatable :: work(:), xlapck(:)

!
!open(unit=99, file='tmpmat.data', form='unformatted')
!rewind(99)
!read(99)norbs
!
norbs = 2000
allocate(dtmp1(norbs,norbs), dtmp2(norbs,norbs), STAT=istat)
allocate(ctmp1(norbs,norbs), ctmp4(norbs,norbs), STAT=istat)
allocate(ctmp2(norbs), xdiis(norbs), ctmp3(norbs), STAT=istat)
!!
!read(99)((dtmp1(j,i), i=1,norbs), j=1,norbs)
!read(99)((dtmp2(j,i), i=1,norbs), j=1,norbs)
!close(99)
!
do ni=1,norbs
 do nj=1,ni
  if (ni .eq. nj) then
    ctmp1(ni,nj) = dcmplx(1.d0, 1.d0)
  else if (nj .eq. ni - 1) then
    ctmp1(ni,nj) = dcmplx(2.d0, 5.d-1)
  else
    ctmp1(ni,nj) = dcmplx(1.d0 / dble(ni), 1.d0 / dble(nj))
  end if
! 
  if (ni .ne. nj) ctmp1(nj,ni) = -ctmp1(ni,nj)
!
 enddo
enddo
!
!ctmp1(1:norbs,1:norbs) = dcmplx( dtmp1(1:norbs, 1:norbs), &
!                                 dtmp2(1:norbs, 1:norbs) )
do nj=1,norbs
 ctmp2(nj) = dcmplx(dble(nj), 5.d-1)
end do
!
! solve linear problem ::  ctmp1 * x = ctmp2
!
lwork = norbs**2
allocate(ipiv(norbs), work(lwork), xlapck(norbs), STAT=istat)
!
xlapck(1:norbs) = ctmp2(1:norbs)
!
ctmp4(1:norbs, 1:norbs) = ctmp1(1:norbs, 1:norbs)
call cpu_time(cpu1)
call zgesv(norbs, 1, ctmp4, norbs, ipiv, xlapck, norbs, info)
call cpu_time(cpu2)
!
write(6,*)
write(6,*)'                  norbs = ', norbs
write(6,*)' out from zgesv,   info = ', info
!write(6,*)' xlapck(1), xlapck(norbs) ', xlapck(1), xlapck(norbs)
!write(6,*)ctmp1(1,1), ctmp1(norbs, norbs)
call flush(6)
!
ctmp3(1:norbs) = dcmplx(0.d0, 0.d0)
call zgemv('normal', norbs, norbs, 1.d0, ctmp1, norbs, xlapck, 1, 0.d0, ctmp3, 1)
!
!write(6,*)
!write(6,*)' check zgesv result'
!do ni=1,norbs
! write(6,*)ctmp3(ni), ctmp2(ni)
!enddo
!call flush(6)
!
ctmp4(1:norbs, 1:norbs) = ctmp1(1:norbs, 1:norbs)
!
call cpu_time(cpu3)
!call diis0(norbs, norbs, ctmp4, ctmp2, xdiis, ierror)
call bicg0(norbs, ctmp4, xdiis, ctmp2, ierror)
call cpu_time(cpu4)
!
write(6,*)
write(6,*)' out from diis, ierror = ', ierror
call flush(6)
!
!stop
!
ctmp3(1:norbs) = dcmplx(0.d0, 0.d0)
call zgemv('normal', norbs, norbs, 1.d0, ctmp1, norbs, xdiis, 1, 0.d0, ctmp3, 1)
write(6,*)
write(6,*)' check diis result'
do ni=1,norbs
 write(6,*)ctmp3(ni), ctmp2(ni)
enddo
call flush(6)
!
!stop
!
diff = 0.d0
summ = 0.d0
write(6,*)' xdiis,    xlapck '
do ni=1,norbs
 diff = diff + cdabs(xdiis(ni) - xlapck(ni))
 summ = summ + cdabs(xdiis(ni) + xlapck(ni))
! write(6,*)xdiis(ni), xlapck(ni)
end do
write(6,*)
write(6,*)' diff(xdiis, xlapck) = ', diff 
write(6,*)' summ(xdiis, xlapck) = ', summ
write(6,*)' ratio               = ', diff / summ
call flush(6)
!
write(6,*)
write(6,*)' cpu time for norbs =     ', norbs
write(6,*)' time from calling lapack ', cpu2 - cpu1
write(6,*)' time from calling diis   ', cpu4 - cpu3
write(6,*)
write(6,*)' test done '
call flush(6)
!
deallocate(ipiv, work, xlapck,  STAT=istat)
deallocate(ctmp2, ctmp3, ctmp4, xdiis, STAT=istat)
deallocate(dtmp1, dtmp2, ctmp1, STAT=istat)
!
end subroutine testdiis
