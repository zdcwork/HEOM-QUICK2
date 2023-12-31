      double precision function dnrm2_long ( n, dx, incx) 
!      integer i, incx, ix, j, n, next 
      integer*8 i, ix, j, n, next 
      integer incx
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one 
      data   zero, one /0.0d0, 1.0d0/ 
!                                                                       
!     euclidean norm of the n-vector stored in dx() with storage        
!     increment incx .                                                  
!     if    n .le. 0 return with result = 0.                            
!     if n .ge. 1 then incx must be .ge. 1                              
!                                                                       
!           c.l.lawson, 1978 jan 08                                     
!     modified to correct problem with negative increment, 8/21/90.     
!     modified to correct failure to update ix, 1/25/92.                
!                                                                       
!     four phase method     using two built-in constants that are       
!     hopefully applicable to all machines.                             
!         cutlo = maximum of  dsqrt(u/eps)  over all known machines.    
!         cuthi = minimum of  dsqrt(v)      over all known machines.    
!     where                                                             
!         eps = smallest no. such that eps + 1. .gt. 1.                 
!         u   = smallest positive no.   (underflow limit)               
!         v   = largest  no.            (overflow  limit)               
!                                                                       
!     brief outline of algorithm..                                      
!                                                                       
!     phase 1    scans zero components.                                 
!     move to phase 2 when a component is nonzero and .le. cutlo        
!     move to phase 3 when a component is .gt. cutlo                    
!     move to phase 4 when a component is .ge. cuthi/m                  
!     where m = n for x() real and m = 2*n for complex.                 
!                                                                       
!     values for cutlo and cuthi..                                      
!     from the environmental parameters listed in the imsl converter    
!     document the limiting values are as follows..                     
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)                         
!                   thus cutlo = 2**(-51) = 4.44089e-16                 
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.          
!                   thus cuthi = 2**(63.5) = 1.30438e19                 
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.             
!                   thus cutlo = 2**(-33.5) = 8.23181d-11               
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19                    
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /                        
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /                        
      data cutlo, cuthi / 8.232d-11,  1.304d19 / 
!                                                                       
      if(n .gt. 0) go to 10 
         dnrm2_long  = zero 
         go to 300 
!                                                                       
   10 assign 30 to next 
      sum = zero 
      i = 1 
      if( incx .lt. 0 )i = (-n+1)*incx + 1 
      ix = 1 
!                                                 begin main loop       
   20    go to next,(30, 50, 70, 110) 
   30 if( dabs(dx(i)) .gt. cutlo) go to 85 
      assign 50 to next 
      xmax = zero 
!                                                                       
!                        phase 1.  sum is zero                          
!                                                                       
   50 if( dx(i) .eq. zero) go to 200 
      if( dabs(dx(i)) .gt. cutlo) go to 85 
!                                                                       
!                                prepare for phase 2.                   
      assign 70 to next 
      go to 105 
!                                                                       
!                                prepare for phase 4.                   
!                                                                       
  100 continue 
      ix = j 
      assign 110 to next 
      sum = (sum / dx(i)) / dx(i) 
  105 xmax = dabs(dx(i)) 
      go to 115 
!                                                                       
!                   phase 2.  sum is small.                             
!                             scale to avoid destructive underflow.     
!                                                                       
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75 
!                                                                       
!                     common code for phases 2 and 4.                   
!                     in phase 4 sum is large.  scale to avoid overflow.
!                                                                       
  110 if( dabs(dx(i)) .le. xmax ) go to 115 
         sum = one + sum * (xmax / dx(i))**2 
         xmax = dabs(dx(i)) 
         go to 200 
!                                                                       
  115 sum = sum + (dx(i)/xmax)**2 
      go to 200 
!                                                                       
!                                                                       
!                  prepare for phase 3.                                 
!                                                                       
   75 sum = (sum * xmax) * xmax 
!                                                                       
!                                                                       
!     for real or d.p. set hitest = cuthi/n                             
!     for complex      set hitest = cuthi/(2*n)                         
!                                                                       
   85 hitest = cuthi/float( n ) 
!                                                                       
!                   phase 3.  sum is mid-range.  no scaling.            
!                                                                       
      do 95 j = ix,n 
      if(dabs(dx(i)) .ge. hitest) go to 100 
         sum = sum + dx(i)**2 
         i = i + incx 
   95 continue 
      dnrm2_long = dsqrt( sum ) 
      go to 300 
!                                                                       
  200 continue 
      ix = ix + 1 
      i = i + incx 
      if( ix .le. n ) go to 20 
!                                                                       
!              end of main loop.                                        
!                                                                       
!              compute square root and adjust for scaling.              
!                                                                       
      dnrm2_long = xmax * dsqrt(sum) 
  300 continue 
      return 
      END                                           
