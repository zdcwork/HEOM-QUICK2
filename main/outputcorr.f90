subroutine outputcorr
use matmod
use tmpmatmod
implicit none
include '../include/sizes'
include '../include/common'
!
integer :: ialf, isgn, iorbs, ispin, imats, iorbs2, ifff,i
real*8  :: hbarscale
!
open(unit=79, file='res_corr.data', status='unknown')
!
i=0
hbarscale=1.d0
if (funits .eq. 2) then
   hbarscale = hbar
end if
!
write(79,*) nsgn
write(79,*) nalf
write(79,*) ncor
write(79,*) nspin
write(79,*) nvar2
!
iorbs2 = 1
ispin = 1
ialf = 1
isgn = 1
!
!do iorbs2=1, nvar2 
!  do ispin=1, nspin
    if(psdfff) then
      do imats=1, ncor
        if(imats .le. ndrude + nmats) then
!          do ialf=1, nalf
!            do isgn=1, nsgn
               i = i + 1
               write(79,100) i, iorbs2, ispin, imats, ialf, isgn, dble(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, &
                             dimag(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, & 
                             dimag(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, &
                             dimag(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, 0
!            end do
!          end do
        else 
          ifff = mfff1d(imats - ndrude - nmats) - 1
!          do ialf=1, nalf
!            do isgn=1, nsgn
               i = i + 1
               write(79,100) i, iorbs2, ispin, imats, ialf, isgn, dble(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, &
                             dimag(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, & 
                             dimag(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, &
                             dimag(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, ifff
!            end do
!          end do
        end if
      end do
    else
     do imats=1, ncor
!       do ialf=1, nalf
!         do isgn=1, nsgn
            i = i + 1
            write(79,100) i, iorbs2, ispin, imats, ialf, isgn, dble(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, &
                             dimag(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, & 
                             dimag(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, &
                             dimag(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, 0
!         end do
!        end do
      end do
    end if
!  end do
!end do
!
close(79)

100 format(I3,1x,5(I3,1x),6(1x, f16.10),2I)
101 format(I3,1x,5(I3,1x),6(1x, f16.10),2I)

end subroutine outputcorr
