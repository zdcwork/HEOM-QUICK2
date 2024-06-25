subroutine resume_cf(tt, iop)
use matmod
use corrfuncmod
use sparsemod
implicit none
include '../include/sizes'
include '../include/common'
!
! Format of binary file <TAPE.resume_cf> 
!   time, dt_dos
!   ntier, norbs, nspin, ncor, nalf
!   nunk  
!   iorbs_dos, ispin_dos
!   ljw_dos
!   lunkcf_spa (lsparse = .true.)    
!   brhoh
!   brhoa
!   bdrhoh (ljw_dos = .true.)
!   bdrhoa (ljw_dos = .true.)
!
integer, intent(in)    :: iop
real*8,  intent(inout) :: tt
integer                :: ni, nj, nk, nl, nm, nn, istat
logical                :: la
integer*8              :: lni, lnj
real*8                 :: dtmp1, dtmp2
!
if (iop .eq. 1) then   ! read
  open(unit=39, file='TAPE.resume_cf', form='unformatted', status='old', iostat=istat)
  if (istat .ne. 0) then
     write(6,*)
     write(6,*)'resume_cf: error! TAPE.resume_cf does not exist! '
     stop
  end if
  rewind(39)
  read(39)tt, dtmp1
  if (dabs(dtmp1 - dt_dos) .ge. dpico) then
! allow change of dt_dos when resuming the corrfunc job 
! be careful about the subsequent fourier transform (uneven time step is used)
     write(6,*)'resume_cf: warning! change in dt_dos identified from TAPE.resume_cf'
     write(6,*)'resume_cf: dt_read, dt_dos ', dtmp1, dt_dos
     write(6,*)'resume_cf: be careful about fourier transform to frequency domain  '
  end if
!
  read(39)ni, nj, nk, nl, nm
  if (ni .ne. ntier .or. nj .ne. norbs .or. nk .ne. nspin .or. & 
      nl .ne. ncor .or. nm .ne. nalf) then
     write(6,*)'resume_cf: error reading line 2 of TAPE.resume_cf'
     write(6,*)ni, nj, nk, nl, nm
     write(6,*)ntier, norbs, nspin, ncor, nalf
     stop
  end if
!
  read(39)lni
  if (lni .ne. nunk) then
     write(6,*)
     write(6,*)'resume_cf: error reading line 3 of TAPE.resume_cf'
     write(6,*)lni, nunk
     stop
  end if
!
  read(39)ni, nj
  if (ni .ne. iorbs_dos .or. nj .ne. ispin_dos) then
     write(6,*)
     write(6,*)'resume_cf: error reading line 4 of TAPE.resume_cf'
     write(6,*)ni, nj
     write(6,*)iorbs_dos, ispin_dos
     stop
  end if
!
  read(39)la
  if (la .ne. ljw_dos) then
     write(6,*)
     write(6,*)'resume_cf: error reading line 5 of TAPE.resume_cf'
     write(6,*)la, ljw_dos
     stop
  end if
!
  if (lsparse) then
     read(39)lni
     if (lni .ne. lunkcf_spa) then
        write(6,*)
        write(6,*)'resume_cf: error, wrong lunkcf_spa from <TAPE.resume_cf> '
        write(6,*)'resume_cf: ', lni, lunkcf_spa
        stop
     end if
     do lni=1,lunkcf_spa
        read(39)brhoh_spa(lni)
     end do
     do lni=1,lunkcf_spa
        read(39)brhoa_spa(lni)
     end do
     if (ljw_dos) then
        do lni=1,lunkcf_spa
           read(39)bdrhoh_spa(lni)
        end do
        do lni=1,lunkcf_spa
           read(39)bdrhoa_spa(lni)
        end do
     end if
  else
     do lni=1,nunk
        read(39)((brhoh(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
     do lni=1,nunk
        read(39)((brhoa(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
     if (ljw_dos) then
        do lni=1,nunk
           read(39)((bdrhoh(ni,nj,lni), ni=1,nrho), nj=1,nrho)
        end do
        do lni=1,nunk
           read(39)((bdrhoa(ni,nj,lni), ni=1,nrho), nj=1,nrho)
        end do
     end if
  end if
!
  close(39)
!
else if (iop .eq. 2) then ! write
  open(unit=39, file='TAPE.resume_cf', form='unformatted', status='unknown', iostat=istat)
  rewind(39)
  write(39)tt, dt_dos
  write(39)ntier, norbs, nspin, ncor, nalf
  write(39)nunk
  write(39)iorbs_dos, ispin_dos
  write(39)ljw_dos
  if (lsparse) then
     write(39)lunkcf_spa
     do lni=1,lunkcf_spa
        write(39)brhoh_spa(lni)
     end do
     do lni=1,lunkcf_spa
        write(39)brhoa_spa(lni)
     end do
     if (ljw_dos) then
        do lni=1,lunkcf_spa
           write(39)bdrhoh_spa(lni)
        end do
        do lni=1,lunkcf_spa
           write(39)bdrhoa_spa(lni)
        end do
     end if
  else
     do lni=1,nunk
       write(39)((brhoh(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
     do lni=1,nunk
       write(39)((brhoa(ni,nj,lni), ni=1,nrho), nj=1,nrho)
     end do
     if (ljw_dos) then
        do lni=1,nunk
          write(39)((bdrhoh(ni,nj,lni), ni=1,nrho), nj=1,nrho)
        end do
        do lni=1,nunk
          write(39)((bdrhoa(ni,nj,lni), ni=1,nrho), nj=1,nrho)
        end do
     end if
  end if
  close(39)
  write(6,*)'resume_cf: TAPE.resume_cf updated at tt = ', tt
!
else
  write(6,*)'resume_cf: error! unknown iop for resume_cf ', iop
  stop
end if
!
call flush(6)
!
end subroutine resume_cf
