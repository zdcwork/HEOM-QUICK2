subroutine corth(nm,n,low,igh,ar,ai,ortr,orti)
   implicit none
   integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
   real(kind=8) ar(nm,n),ai(nm,n),ortr(igh),orti(igh)
   real(kind=8) f,g,h,fi,fr,scale,pythag

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200

      do 180 m = kp1, la
         h = 0.0d0
         ortr(m) = 0.0d0
         orti(m) = 0.0d0
         scale = 0.0d0
         !scale column (algol tol then not needed)
         do 90 i = m, igh
   90    scale = scale + dabs(ar(i,m-1)) + dabs(ai(i,m-1))

         if (scale .eq. 0.0d0) go to 180
         mp = m + igh
         !for i=igh step -1 until m do
         do 100 ii = m, igh
            i = mp - ii
            ortr(i) = ar(i,m-1) / scale
            orti(i) = ai(i,m-1) / scale
            h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
  100    continue

         g = dsqrt(h)
         f = pythag(ortr(m),orti(m))
         if (f .eq. 0.0d0) go to 103
         h = h + f * g
         g = g / f
         ortr(m) = (1.0d0 + g) * ortr(m)
         orti(m) = (1.0d0 + g) * orti(m)
         go to 105

  103    ortr(m) = g
         ar(m,m-1) = scale
         !form (i-(u*ut)/h) * a
  105    do 130 j = m, n
            fr = 0.0d0
            fi = 0.0d0
            !for i=igh step -1 until m do
            do 110 ii = m, igh
               i = mp - ii
               fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
               fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
  110       continue

            fr = fr / h
            fi = fi / h

            do 120 i = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(i) + fi * orti(i)
               ai(i,j) = ai(i,j) - fr * orti(i) - fi * ortr(i)
  120       continue

  130    continue
         !form (i-(u*ut)/h)*a*(i-(u*ut)/h)
         do 160 i = 1, igh
            fr = 0.0d0
            fi = 0.0d0
            !for j=igh step -1 until m do
            do 140 jj = m, igh
               j = mp - jj
               fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
               fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
  140       continue

            fr = fr / h
            fi = fi / h

            do 150 j = m, igh
               ar(i,j) = ar(i,j) - fr * ortr(j) - fi * orti(j)
               ai(i,j) = ai(i,j) + fr * orti(j) - fi * ortr(j)
  150       continue

  160    continue

         ortr(m) = scale * ortr(m)
         orti(m) = scale * orti(m)
         ar(m,m-1) = -g * ar(m,m-1)
         ai(m,m-1) = -g * ai(m,m-1)
  180 continue

  200 return

end subroutine
!++++++++++++++++++++++++++++++
