subroutine calcams(norbs, nspin, iorbs, ispin, lmat, amsout)
implicit none
!
! purpose : construct annilation matrix a_{ms} in occupation number representation
!
! input 
!  norbs  : number of energy levels
!  nspin  : number of spin degree of freedom involved
!  iorbs  : the m-th energy level on which an electron is annihilated
!  ispin  : for spin system (nspin = 2),
!           1, spin up   electron
!           2, spin down electron
!   lmat  : dimension of matrix a_{ms}, lmat = (nspin*2)**norbs
!
! output 
! amsout  : the output matrix a_{ms}
!--------------------------------------------------------------------------
!                                 ALGORITHM
!
! (1) basis set
!
!     for spinless systems, each orbital can be vacant or occupied, denoted by
!     0 and 1, respectively; whereas for spin systems, each orbital can be vacant, 
!     populated by spin-up electron, spin-down electron, or doubly occupied, 
!     denoted by 0, u, d and 2, respectively.
!
!     when forming a complete set for a  multi-orbital system, the first orbital
!     changes its population most frequently, while the last orbital changes its 
!     population status only once.
!
!     e.g. for norbs=2, nspin=1, the matrix representation for 
!          c_1 = |00> <10| + |01> <11| (iorbs=1, ispin=1) is as follows, 
!
!                <00| <10| <01| <11|         
!              ---------------------
!         |00> |   0    1    0    0   
!         |10> |   0    0    0    0
!         |01> |   0    0    0    1
!         |11> |   0    0    0    0 
!
!     e.g. for norbs=1, nspin=2, the matrix representation for 
!          c_d = |0> <d| + |u> <2| (iorbs=1, ispin=2) is as follows, 
!
!                <0|  <u|  <d|  <2|         
!              ---------------------
!         |0> |   0    0    1    0   
!         |u> |   0    0    0    1
!         |d> |   0    0    0    0
!         |2> |   0    0    0    0 
!
integer, intent(in)  :: norbs, nspin, iorbs, ispin, lmat
real*8,  intent(out) :: amsout(lmat,*)
integer :: ni, nj, len0, irow0, irow1, icol0, icol1
!
! check input
!
if (lmat .ne. (nspin * 2)**norbs) then
 write(6,*)
 write(6,*)' error input for calcams ', norbs, nspin, lmat
 stop
end if
if (nspin .eq. 2 .and. ispin .ne. 1 .and. ispin .ne. 2) then
 write(6,*)
 write(6,*)' error, invalid spin index in calcams ', nspin, ispin
 stop
end if
amsout(1:lmat, 1:lmat) = 0.d0
!
if (nspin .eq. 1) then      ! spinless system
 len0 = lmat / 2**(norbs - iorbs + 1)
 do ni=1,2**(norbs-iorbs)
   irow0 = (ni - 1) * len0 * 2 + 1
   icol0 = (ni - 1) * len0 * 2 + 1 + len0
   do nj=1,len0
     amsout(irow0+nj-1, icol0+nj-1) = 1.d0
   end do
 end do
!
else if (nspin .eq. 2) then ! spin system
 len0 = lmat / 4**(norbs - iorbs + 1)
 if (ispin .eq. 1) then     ! spin up annihilation 
   do ni=1,4**(norbs-iorbs)
     irow0 = (ni - 1) * len0 * 4 + 1
     irow1 = (ni - 1) * len0 * 4 + 1 + len0 * 2
     icol0 = (ni - 1) * len0 * 4 + 1 + len0
     icol1 = (ni - 1) * len0 * 4 + 1 + len0 * 3
     do nj=1,len0
       amsout(irow0+nj-1, icol0+nj-1) = 1.d0
       amsout(irow1+nj-1, icol1+nj-1) = 1.d0
     end do
   end do
 else                       ! spin down annihilation
   do ni=1,4**(norbs-iorbs) 
     irow0 = (ni - 1) * len0 * 4 + 1
     irow1 = (ni - 1) * len0 * 4 + 1 + len0 
     icol0 = (ni - 1) * len0 * 4 + 1 + len0 * 2
     icol1 = (ni - 1) * len0 * 4 + 1 + len0 * 3
     do nj=1,len0
       amsout(irow0+nj-1, icol0+nj-1) = 1.d0
       amsout(irow1+nj-1, icol1+nj-1) = 1.d0
     end do
   enddo
 end if
!
else
 write(6,*)
 write(6,*)' error, unknown nspin input for calcams ', nspin
 stop
end if
end subroutine calcams
