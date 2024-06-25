subroutine checkanticomm
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ni, nj, nk, nl
integer :: iorbs, ispin
real*8 :: dtmp1
!
ams(1:nrho, 1:nrho) = 0.d0
do ni=1,nrho
  ams(ni,ni) = 1.d0
end do
!
do ni=1,nspin
  do nj=1,nspin
    do nk=1,norbs
      do nl=1,nk
        dmtmp1(1:nrho, 1:nrho) = amsori(1:nrho, 1:nrho, nk, ni)
        dmtmp2(1:nrho, 1:nrho) = amsori(1:nrho, 1:nrho, nl, nj)
!
! {c1, c2}
!
        call dgemm('n', 'n', nrho, nrho, nrho, 1.d0, dmtmp1, nrho, dmtmp2, nrho, &
                   0.d0, dmtmp3, nrho)
        call dgemm('n', 'n', nrho, nrho, nrho, 1.d0, dmtmp2, nrho, dmtmp1, nrho, &
                   1.d0, dmtmp3, nrho)
        call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
        if (dtmp1 .ge. dnano) then
          write(6,*)
          write(6,*)' error! {c1, c2} ', dtmp1
          write(6,*)' 1 => iorbs = ', nk, ' ispin = ', ni
          write(6,*)' 2 => iorbs = ', nl, ' ispin = ', nj
          call printams(nk, ni)
          call printams(nl, nj)
          stop
        end if
!
! {c1, c2+} 
!
        call dgemm('n', 'c', nrho, nrho, nrho, 1.d0, dmtmp1, nrho, dmtmp2, nrho, &
                   0.d0, dmtmp3, nrho)
        call dgemm('c', 'n', nrho, nrho, nrho, 1.d0, dmtmp2, nrho, dmtmp1, nrho, &
                   1.d0, dmtmp3, nrho)
        if (nl .eq. nk .and. ni .eq. nj) then
          dmtmp3(1:nrho, 1:nrho) = dmtmp3(1:nrho, 1:nrho) - ams(1:nrho, 1:nrho)
        end if
        call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
        if (dtmp1 .ge. dnano) then
          write(6,*)
          write(6,*)' error! {c1, c2+} ', dtmp1
          write(6,*)' 1 => iorbs = ', nk, ' ispin = ', ni
          call printams(nk, ni)
          write(6,*)' 2 => iorbs = ', nl, ' ispin = ', nj
          call printams(nl, nj)
          stop
        end if
      end do
    end do
  end do
end do
!
write(6,*)
write(6,*)' anti-commutation relation check pass! '
call flush(6)
!
amsall = amsori
!
! build occupancy operator for each orbital and each spin direction
!
do ispin=1,nspin
   do iorbs=1,norbs
      call dgemm('c', 'n', nrho, nrho, nrho, 1.d0, amsori(1,1,iorbs,ispin), nrho, &
                 amsori(1,1,iorbs,ispin), nrho, 0.d0, dmtmp1, nrho)
      znms(1:nrho,1:nrho,iorbs,ispin) = dmtmp1(1:nrho,1:nrho) * cunity
   end do
end do
!
goto 200
!
amsall = dabs(amsori)
do ni=1,nspin
  do nj=1,nspin
    do nk=1,norbs
      do nl=1,nk
        dmtmp1(1:nrho, 1:nrho) = amsall(1:nrho, 1:nrho, nk, ni)
        dmtmp2(1:nrho, 1:nrho) = amsall(1:nrho, 1:nrho, nl, nj)
!
        if (ni .eq. nj .and. nk .eq. nl) then
! {c1, c1+} 
          call dgemm('n', 'c', nrho, nrho, nrho, 1.d0, dmtmp1, nrho, dmtmp2, nrho, &
                     0.d0, dmtmp3, nrho)
          call dgemm('c', 'n', nrho, nrho, nrho, 1.d0, dmtmp2, nrho, dmtmp1, nrho, &
                     1.d0, dmtmp3, nrho)
          dmtmp3 = dmtmp3 - ams
          call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
          if (dtmp1 .ge. dnano) then
            write(6,*)
            write(6,*)' error! {c1, c1+} ', dtmp1
            write(6,*)' iorbs = ', nk, ' ispin = ', ni
            call printams(nk, ni)
            stop
          end if
        else
! [c1, c2]
          call dgemm('n', 'n', nrho, nrho, nrho,  1.d0, dmtmp1, nrho, dmtmp2, nrho, &
                     0.d0, dmtmp3, nrho)
          call dgemm('n', 'n', nrho, nrho, nrho, -1.d0, dmtmp2, nrho, dmtmp1, nrho, &
                     1.d0, dmtmp3, nrho)
          call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
          if (dtmp1 .ge. dnano) then
            write(6,*)
            write(6,*)' error! [c1, c2] ', dtmp1
            write(6,*)' 1 => iorbs = ', nk, ' ispin = ', ni
            write(6,*)' 2 => iorbs = ', nl, ' ispin = ', nj
            call printams(nk, ni)
            call printams(nl, nj)
            stop
          end if
! [c1, c2+] 
          call dgemm('n', 'c', nrho, nrho, nrho,  1.d0, dmtmp1, nrho, dmtmp2, nrho, &
                     0.d0, dmtmp3, nrho)
          call dgemm('c', 'n', nrho, nrho, nrho, -1.d0, dmtmp2, nrho, dmtmp1, nrho, &
                     1.d0, dmtmp3, nrho)
          call dmaxmat(nrho, nrho, dmtmp3, nrho, dtmp1)
          if (dtmp1 .ge. dnano) then
            write(6,*)
            write(6,*)' error! [c1, c2+] ', dtmp1
            write(6,*)' 1 => iorbs = ', nk, ' ispin = ', ni
            call printams(nk, ni)
            write(6,*)' 2 => iorbs = ', nl, ' ispin = ', nj
            call printams(nl, nj)
            stop
          end if
        end if
      end do
    end do
  end do
end do
write(6,*)'      commutation relation check pass! '
call flush(6)
!
200 continue
!
zamsall(1:nrho,1:nrho,1:norbs,1:nspin) = amsall(1:nrho,1:nrho,1:norbs,1:nspin) * cunity
!
end subroutine checkanticomm
