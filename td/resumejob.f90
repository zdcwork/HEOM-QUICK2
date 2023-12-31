subroutine resumejob(tt, iop)
use matmod
use sparsemod
use hbmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer, intent(in)    :: iop
real*8,  intent(inout) :: tt
character*10           :: chtmp0
integer :: ni, nj, nk, nl, nm, nn, istat, ncount
logical :: ltmp1 
integer*8 :: lni, lnj
real*8  :: dtmp1, dtmp2
!
if (iop .eq. 1) then   ! read
  open(unit=23, file='TAPE.resume', form='binary', status='old', iostat=istat)
  if (istat .ne. 0) then
     write(6,*)
     write(6,*)'resumejob: error! TAPE.resume does not exist! '
     stop
  end if
  rewind(23)
  read(23)tt
  read(23)ni, nj, nk, nl, nm, nn
  if (ni .ne. ntier .or. nj .ne. norbs .or. nk .ne. nspin .or. & 
      nl .ne. ncor .or. nm .ne. nalf .or. nn .ne. numfff) then
     write(6,*)
     write(6,*)'resumejob: error! bad TAPE.resume to read 1! '
     write(6,*)ni, nj, nk, nl, nm, nn
     write(6,*)ntier, norbs, nspin, ncor, nalf, numfff
     stop
  end if
  read(23)lni
  if (lni .ne. nunk) then
     write(6,*)
     write(6,*)'resumejob: error! bad TAPE.resume to read 2! '
     write(6,*)lni, nunk
     stop
  end if
  read(23)ltmp1
  if (ltmp1 .ne. lhb) then
     write(6,*)
     write(6,*)'resumejob: error! bad TAPE.resume to read 3! '
     write(6,*)ltmp1, lhb
     stop
  end if
  if (lhb) then
     read(23)ni, nj, nk, nl
     if (ni .ne. npade_hb .or. nj .ne. ndrude_hb .or. nk .ne. nmode_hb .or. nl .ne. itype_hb) then
        write(6,*)
        write(6,*)'resumejob: error! bad TAPE.resume to read 4! '
        write(6,*)ni, nj, nk, nl
        write(6,*)npade_hb, ndrude_hb, nmode_hb, itype_hb
        stop
     end if
  end if
!
  if (lsparse) then
     read(23)lni
     if (lni .ne. lunk_spa) then
        write(6,*)
        write(6,*)'resumejob: error! bad TAPE.resume to read 3! '
        write(6,*)lni, lunk_spa
        stop
     end if
     do lni=1,lunk_spa
        read(23)rho_spa(lni)
     end do
  else
     do lni=1,nunk
        read(23)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
  end if
  close(23)
!
else if (iop .eq. 2) then ! write
  ncount = nfbase + int((tt + dnano) / dt) / nresume
!  write(6,*)'resumejob: writing at ', ncount
!  call flush(6)
  write(chtmp0,'(I10)')ncount
  open(unit=23, file='TAPE_'//trim(adjustl(chtmp0))//'.resume', form='binary', &
       status='unknown', iostat=istat)
  rewind(23)
  write(23)tt
  write(23)ntier, norbs, nspin, ncor, nalf, numfff
  write(23)nunk
  write(23)lhb
  if (lhb) then
     write(23)npade_hb, ndrude_hb, nmode_hb, itype_hb
  end if
  if (lsparse) then
     write(23)lunk_spa
     do lni=1,lunk_spa
        write(23)rho_spa(lni)
     end do
  else
     do lni=1,nunk
       write(23)((rho(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
  end if
  close(23)
!
  write(6,*)
  write(6,*)'resumejob: Resume file updated at tt = ', tt
  call flush(6)
!
else
  write(6,*)
  write(6,*)'resumejob: error! unknown iop = ', iop
  stop
end if
!
if (.not. lsparse) then
   do ni=1,ntier
     dtmp1 = 0.d0
     lnj   = -1
     do lni=nfirst(ni), nlast(ni)
       do nk=1,nrho
         do nj=1,nrho
           if (dtmp1 .le. cdabs(rho(nj,nk,lni))) then
              dtmp1 = cdabs(rho(nj,nk,lni))
              lnj   = lni
           end if
         end do
       end do
     end do
     write(6,103)lnj, lnj-nfirst(ni), dtmp1
   end do
end if
call flush(6)
103 format(2x, I12, 2x, I12, 2x, e18.7e3)
!
end subroutine resumejob
