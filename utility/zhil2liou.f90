subroutine zhil2liou(lefta, transa, dimh, zhmat, lhmat, zlmat, llmat)
implicit none
!
integer,     intent(in)  :: dimh, lhmat, llmat
character*1, intent(in)  :: lefta, transa
complex*16,  intent(in)  :: zhmat(lhmat,*)
complex*16,  intent(out) :: zlmat(llmat,*)
!
integer    :: dim1, ni, nj, nm, nn, n1, n2
logical    :: lsame
external   :: lsame
complex*16 :: czero
!
if (dimh .gt. lhmat .or. dimh**2 .gt. llmat) then
   write(6,*)
   write(6,*)'zhil2liou: error! dimh, lhmat, llmat ', dimh, lhmat, llmat
   stop
end if
czero = dcmplx(0.d0, 0.d0)
dim1  = dimh**2
zlmat(1:dim1,1:dim1) = czero
!
if ( lsame(lefta, 'L') ) then
   if ( lsame(transa, 'N') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nm=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nj - 1) * dimh + nm
               zlmat(n1,n2) = zhmat(ni,nm)
            end do
         end do
      end do
   else if ( lsame(transa, 'R') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nm=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nj - 1) * dimh + nm
               zlmat(n1,n2) = dconjg(zhmat(ni,nm))
            end do
         end do
      end do
   else if ( lsame(transa, 'T') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nm=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nj - 1) * dimh + nm
               zlmat(n1,n2) = zhmat(nm,ni)
            end do
         end do
      end do
   else if ( lsame(transa, 'C') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nm=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nj - 1) * dimh + nm
               zlmat(n1,n2) = dconjg(zhmat(nm,ni))
            end do
         end do
      end do
   else
      write(6,*)
      write(6,*)'zhil2liou: error! unknown transa ', transa
      stop
   end if
!
else if ( lsame(lefta, 'R') ) then
   if ( lsame(transa, 'N') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nn=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nn - 1) * dimh + ni
               zlmat(n1,n2) = zhmat(nn,nj)
            end do
         end do
      end do
   else if ( lsame(transa, 'R') ) then
      do nj=1,dimh
         do ni=1,dimh
            do nn=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nn - 1) * dimh + ni
               zlmat(n1,n2) = dconjg(zhmat(nn,nj))
            end do
         end do
      end do
   else if ( lsame(transa, 'T') ) then
      do nj=1,dimh 
         do ni=1,dimh
            do nn=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nn - 1) * dimh + ni
               zlmat(n1,n2) = zhmat(nj,nn)
            end do
         end do
      end do
   else if ( lsame(transa, 'C') ) then
      do nj=1,dimh 
         do ni=1,dimh
            do nn=1,dimh
               n1 = (nj - 1) * dimh + ni
               n2 = (nn - 1) * dimh + ni
               zlmat(n1,n2) = dconjg(zhmat(nj,nn))
            end do
         end do
      end do
   else
      write(6,*)
      write(6,*)'zhil2liou: error! unknown transa ', transa
      stop
   end if
!
else
   write(6,*)
   write(6,*)'zhil2liou: error! unknown lefta ', lefta
   stop
end if
!
return
end subroutine zhil2liou


!
subroutine zh2lmlr(dimh, zahmat, lda, zbhmat, ldb, zlmat, llmat)
implicit none
!
! Y = A * X * B  => LY = L_AB * LX  ---> to check!!!
!
integer, intent(in)     :: dimh, lda, ldb, llmat
complex*16, intent(in)  :: zahmat(lda,*), zbhmat(ldb,*)
complex*16, intent(out) :: zlmat(llmat,*)
!
integer                 :: dim1, ni, nj, nn, nm, n1, n2
!
dim1 = dimh**2
zlmat(1:dim1, 1:dim1) = dcmplx(0.d0,0.d0)
do nj=1,dimh
   do ni=1,dimh
      n1 = (nj - 1) * dimh + ni
      do nn=1,dimh
         do nm=1,dimh
            n2 = (nn - 1) * dimh + nm
            zlmat(n1,n2) = zahmat(ni,nm) * zbhmat(nn,nj)
         end do
      end do
   end do
end do
!
return
end subroutine zh2lmlr

