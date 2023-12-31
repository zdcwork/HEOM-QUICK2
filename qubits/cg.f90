!
subroutine cg(nm,n,ar,ai,wr,wi,matz,zr,zi,ierr)          !,fv1,fv2,fv3
   implicit none
   integer            :: n,nm
   real(kind=8)       :: ar(nm,n),ai(nm,n)
   real(kind=8)       :: wr(n),wi(n)
   real(kind=8)       :: zr(nm,n),zi(nm,n)
   real(kind=8)       :: fv1(n),fv2(n),fv3(n)
   integer            :: matz
   integer            :: ierr
   integer            :: is1,is2

   if (n .le. nm) go to 10
   ierr=10 * n
   go to 50

10 call cbal(nm,n,ar,ai,is1,is2,fv1)

   call corth(nm,n,is1,is2,ar,ai,fv2,fv3)

   if (matz .ne. 0) go to 20

   call comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)

   go to 50

20 call comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)

   if (ierr .ne. 0) go to 50

   call cbabk2(nm,n,is1,is2,fv1,n,zr,zi)

50 return

end subroutine
!++++++++++++++++++++++++++++++++
