!
!     complex division, (cr,ci) = (ar,ai)/(br,bi)
!
subroutine cdiv(ar,ai,br,bi,cr,ci)
   implicit none
   real(kind=8) ar,ai,br,bi,cr,ci
   real(kind=8) s,ars,ais,brs,bis

   s = dabs(br) + dabs(bi)
   ars = ar/s
   ais = ai/s
   brs = br/s
   bis = bi/s
   s = brs**2 + bis**2
   cr = (ars*brs + ais*bis)/s
   ci = (ais*brs - ars*bis)/s

   return

end subroutine
!+++++++++++++++++++++++++
