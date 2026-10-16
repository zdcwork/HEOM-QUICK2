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
!hbarscale=1.d0
!if (funits .eq. 2) then
   hbarscale = hbar
!end if
!
write(79,*) nsgn
write(79,*) nalf
write(79,*) ncor
write(79,*) nspin
write(79,*) nvar2
!
do iorbs2=1, nvar2 
  do ispin=1, nspin
    if(psdfff) then
      do imats=1, ncor
        do ialf=1, nalf
          ! In diff mode, PSD-FFF uses compact per-lead Fano tail layout:
          ! Fano terms begin after ndrude+nmats_alf(ialf).
          ifff = 0
          if (lpsdfff_alf(ialf) .and. imats .gt. ndrude+nmats_alf(ialf) .and. &
              imats .le. ncor_alf(ialf)) then
             ifff = mfff1d(imats - ndrude - nmats_alf(ialf)) - 1
          else if (.not. (ldiff_spe_dec .or. ldiff_spe_dec_sp) .and. imats .gt. ndrude+nmats) then
             ifff = mfff1d(imats - ndrude - nmats) - 1
          end if
            do isgn=1, nsgn
               i = i + 1
               write(79,100) i, iorbs2, ispin, imats, ialf, isgn, dble(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, &
                             dimag(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, & 
                             dimag(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, &
                             dimag(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, ifff
            end do
        end do
      end do
    else
     do imats=1, ncor
       do ialf=1, nalf
         do isgn=1, nsgn
            i = i + 1
            write(79,100) i, iorbs2, ispin, imats, ialf, isgn, dble(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, &
                             dimag(cb(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, & 
                             dimag(cd(iorbs2,ispin,imats,ialf,isgn))*hbarscale**2, dble(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, &
                             dimag(cgama(iorbs2,ispin,imats,ialf,isgn))*hbarscale, 0
         end do
        end do
      end do
    end if
  end do
end do
!write(79,*) ' orb, spin, mats, alf, sgn'
!
!write(79,*) 'hbarscale = ', hbarscale
close(79)

100 format(I4,1x,5(I4,1x),6(1x, E16.8),I4)
101 format(I4,1x,5(I4,1x),6(1x, E16.8),I4)

end subroutine outputcorr
