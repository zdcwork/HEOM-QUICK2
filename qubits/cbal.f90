subroutine cbal(nm,n,ar,ai,low,igh,scale)
   implicit none
   integer              i,j,k,l,m,n,jj,nm,igh,low,iexc
   real(kind=8)         ar(nm,n),ai(nm,n),scale(n)
   real(kind=8)         c,f,g,r,s,b2,radix
   logical              noconv

      radix = 16.0d0

      b2 = radix * radix
      k = 1
      l = n
      go to 100
      ! in-line procedure for row and column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50

      do 30 i = 1, l
         f = ar(i,j)
         ar(i,j) = ar(i,m)
         ar(i,m) = f
         f = ai(i,j)
         ai(i,j) = ai(i,m)
         ai(i,m) = f
   30 continue

      do 40 i = k, n
         f = ar(j,i)
         ar(j,i) = ar(m,i)
         ar(m,i) = f
         f = ai(j,i)
         ai(j,i) = ai(m,i)
         ai(m,i) = f
   40 continue

   50 go to (80,130), iexc

      ! search for rows isolating an eigenvalue and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
      ! for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj

         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (ar(j,i) .ne. 0.0d0 .or. ai(j,i) .ne. 0.0d0) go to 120
  110    continue

         m = l
         iexc = 1
         go to 20
  120 continue

      go to 140
      ! search for columns isolating an eigenvalue and push them left ..........
  130 k = k + 1

  140 do 170 j = k, l

         do 150 i = k, l
            if (i .eq. j) go to 150
            if (ar(i,j) .ne. 0.0d0 .or. ai(i,j) .ne. 0.0d0) go to 170
  150    continue

         m = k
         iexc = 2
         go to 20
  170 continue
      ! now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0d0
      ! iterative loop for norm reduction ..........
  190 noconv = .false.

      do 270 i = k, l
         c = 0.0d0
         r = 0.0d0

         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + dabs(ar(j,i)) + dabs(ai(j,i))
            r = r + dabs(ar(i,j)) + dabs(ai(i,j))
  200    continue
         ! guard against zero c or r due to underflow ..........
         if (c .eq. 0.0d0 .or. r .eq. 0.0d0) go to 270
         g = r / radix
         f = 1.0d0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
         ! now balance ..........
  240    if ((c + r) / f .ge. 0.95d0 * s) go to 270
         g = 1.0d0 / f
         scale(i) = scale(i) * f
         noconv = .true.

         do 250 j = k, n
            ar(i,j) = ar(i,j) * g
            ai(i,j) = ai(i,j) * g
  250    continue

         do 260 j = 1, l
            ar(j,i) = ar(j,i) * f
            ai(j,i) = ai(j,i) * f
  260    continue

  270 continue

      if (noconv) go to 190

  280 low = k
      igh = l

      return

end subroutine
!++++++++++++++++++++++++++++++
