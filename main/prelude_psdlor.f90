subroutine prelude_psdlor
use matmod
implicit none
include "../include/sizes"
include "../include/common"
!
integer :: ni, nj, nk, ntmp0
integer :: istat, ialf
logical :: lexist
real*8  :: dtmp0
!
if (.not. psdlor) then
   write(6,*)
   write(6,*)'prelude_psdlor: Warning! error entrance, psdlor = ', psdlor
   call flush(6)
   return 
end if
!
inquire(file='psd_lor.data', exist=lexist)
if (.not. lexist) then
   write(6,*)
   write(6,*)'prelude_psdlor: Error! file <psd_lor.data> does not exist'
   stop
end if
!
ntmp0 = nmats
open(unit=42, file='psd_lor.data', form='formatted', status='old')
rewind(42)
read(42,*)nmats, nlor
!
if (nmats .le. 0 .or. nlor .le. 0) then
   write(6,*)'prelude_psdlor: Error! nmats, nlor ', nmats, nlor
   stop
end if
!
allocate(dlor_coef (nlor,nalf), STAT=istat)
allocate(dlor_cent (nlor,nalf), STAT=istat)
allocate(dlor_width(nlor,nalf), STAT=istat)
!
rewind(42)
do ialf=1,nalf
   read(42,*)ni, nj
   do nk=1,nlor
      read(42,*)dlor_coef(nk,ialf), dlor_width(nk,ialf), dlor_cent(nk,ialf)
   end do
end do
close(42)
!
write(6,*)'prelude_psdlor: <psd_lor.data> read successfully! '
write(6,*)'prelude_psdlor: Pade point number changed from ', ntmp0
write(6,*)'prelude_psdlor:                             to ', nmats 
call flush(6)
!
! sort and symmetrize the residual reservoir spectral function
! 
if (lsymlor) then
   do ialf=1,nalf
      do ni=1,nlor-1      ! sort center: from negative to positive
         do nj=ni+1,nlor  
            if (dlor_cent(ni,ialf) .gt. dlor_cent(nj,ialf)) then
               dtmp0               = dlor_cent(ni,ialf)            
               dlor_cent(ni,ialf)  = dlor_cent(nj,ialf)
               dlor_cent(nj,ialf)  = dtmp0
               dtmp0               = dlor_coef(ni,ialf)            
               dlor_coef(ni,ialf)  = dlor_coef(nj,ialf)
               dlor_coef(nj,ialf)  = dtmp0
               dtmp0               = dlor_width(ni,ialf)            
               dlor_width(ni,ialf) = dlor_width(nj,ialf)
               dlor_width(nj,ialf) = dtmp0
            end if
         end do
      end do
!
      do ni=1,nlor/2
         nj = nlor - ni + 1
         dtmp0 = 5.d-1 * (dabs(dlor_cent(ni,ialf)) + dabs(dlor_cent(nj,ialf)))
         dlor_cent(ni,ialf) = -dtmp0
         dlor_cent(nj,ialf) =  dtmp0
         dtmp0 = 5.d-1 * (dabs(dlor_width(ni,ialf)) + dabs(dlor_width(nj,ialf)))
         dlor_width(ni,ialf) =  dtmp0
         dlor_width(nj,ialf) =  dtmp0
         dtmp0 = 5.d-1 * (dabs(dlor_coef(ni,ialf)) + dabs(dlor_coef(nj,ialf)))
         if (dlor_coef(ni,ialf) .lt. dlor_coef(nj,ialf)) then
            dlor_coef(ni,ialf) = -dtmp0
            dlor_coef(nj,ialf) =  dtmp0
         else    
            dlor_coef(ni,ialf) =  dtmp0
            dlor_coef(nj,ialf) = -dtmp0
         end if
      end do
   end do
end if
!
100 continue
!
end subroutine prelude_psdlor
