subroutine cbabk2(nm,n,low,igh,scale,m,zr,zi)
   implicit none
   integer       i,j,k,m,n,ii,nm,igh,low
   real(kind=8)  scale(n),zr(nm,m),zi(nm,m)
   real(kind=8)  s

      if (m .eq. 0) go to 200
      if (igh .eq. low) go to 120

      do 110 i = low, igh
         s = scale(i)
         !left hand eigenvectors are back transformed if the foregoing statement is replaced by s=1.0d0/scale(i)
         do 100 j = 1, m
            zr(i,j) = zr(i,j) * s
            zi(i,j) = zi(i,j) * s
  100    continue

  110 continue
      !for i=low-1 step -1 until 1, igh+1 step 1 until n do
  120 do 140 ii = 1, n
         i = ii
         if (i .ge. low .and. i .le. igh) go to 140
         if (i .lt. low) i = low - ii
         k = scale(i)
         if (k .eq. i) go to 140

         do 130 j = 1, m
            s = zr(i,j)
            zr(i,j) = zr(k,j)
            zr(k,j) = s
            s = zi(i,j)
            zi(i,j) = zi(k,j)
            zi(k,j) = s
  130    continue

140   continue

200   return

end subroutine
!+++++++++++++++++++++++++++
