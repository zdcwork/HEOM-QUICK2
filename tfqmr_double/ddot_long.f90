      double precision function ddot_long(n,dx,incx,dy,incy) 
!                                                                       
!     forms the dot product of two vectors.                             
!     uses unrolled loops for increments equal to one.                  
!     jack dongarra, linpack, 3/11/78.                                  
!                                                                       
      double precision dx(1),dy(1),dtemp 
!      integer i,incx,incy,ix,iy,m,mp1,n 
      integer*8 i,ix,iy,m,mp1,n 
      integer incx, incy
!                                                                       
      ddot_long = 0.0d0 
      dtemp = 0.0d0 
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!        code for unequal increments or equal increments                
!          not equal to 1                                               
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        dtemp = dtemp + dx(ix)*dy(iy) 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      ddot_long = dtemp 
      return 
!                                                                       
!        code for both increments equal to 1                            
!                                                                       
!                                                                       
!        clean-up loop                                                  
!                                                                       
   20 m = mod(n,5) 
      if( m .eq. 0 ) go to 40 
      do 30 i = 1,m 
        dtemp = dtemp + dx(i)*dy(i) 
   30 continue 
      if( n .lt. 5 ) go to 60 
   40 mp1 = m + 1 
      do 50 i = mp1,n,5 
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +             &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue 
   60 ddot_long = dtemp 
      return 
      END                                           
