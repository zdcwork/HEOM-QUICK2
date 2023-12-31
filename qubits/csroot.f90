subroutine csroot(xr,xi,yr,yi)
   implicit none
   real(kind=8)  xr,xi,yr,yi
   real(kind=8)  s,tr,ti,pythag

      tr = xr
      ti = xi
      s = dsqrt(0.5d0*(pythag(tr,ti) + dabs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)

      return

end subroutine
!+++++++++++++++++++++++++
