      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info) 
      integer ldx,n,k,job,info 
      double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1), &
     &                 xb(1)                                            
!                                                                       
!     dqrsl applies the output of dqrdc to compute coordinate           
!     transformations, projections, and least squares solutions.        
!     for k .le. min(n,p), let xk be the matrix                         
!                                                                       
!            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))              
!                                                                       
!     formed from columnns jpvt(1), ... ,jpvt(k) of the original        
!     n x p matrix x that was input to dqrdc (if no pivoting was        
!     done, xk consists of the first k columns of x in their            
!     original order).  dqrdc produces a factored orthogonal matrix q   
!     and an upper triangular matrix r such that                        
!                                                                       
!              xk = q * (r)                                             
!                       (0)                                             
!                                                                       
!     this information is contained in coded form in the arrays         
!     x and qraux.                                                      
!                                                                       
!     on entry                                                          
!                                                                       
!        x      double precision(ldx,p).                                
!               x contains the output of dqrdc.                         
!                                                                       
!        ldx    integer.                                                
!               ldx is the leading dimension of the array x.            
!                                                                       
!        n      integer.                                                
!               n is the number of rows of the matrix xk.  it must      
!               have the same value as n in dqrdc.                      
!                                                                       
!        k      integer.                                                
!               k is the number of columns of the matrix xk.  k         
!               must nnot be greater than min(n,p), where p is the      
!               same as in the calling sequence to dqrdc.               
!                                                                       
!        qraux  double precision(p).                                    
!               qraux contains the auxiliary output from dqrdc.         
!                                                                       
!        y      double precision(n)                                     
!               y contains an n-vector that is to be manipulated        
!               by dqrsl.                                               
!                                                                       
!        job    integer.                                                
!               job specifies what is to be computed.  job has          
!               the decimal expansion abcde, with the following         
!               meaning.                                                
!                                                                       
!                    if a.ne.0, compute qy.                             
!                    if b,c,d, or e .ne. 0, compute qty.                
!                    if c.ne.0, compute b.                              
!                    if d.ne.0, compute rsd.                            
!                    if e.ne.0, compute xb.                             
!                                                                       
!               note that a request to compute b, rsd, or xb            
!               automatically triggers the computation of qty, for      
!               which an array must be provided in the calling          
!               sequence.                                               
!                                                                       
!     on return                                                         
!                                                                       
!        qy     double precision(n).                                    
!               qy conntains q*y, if its computation has been           
!               requested.                                              
!                                                                       
!        qty    double precision(n).                                    
!               qty contains trans(q)*y, if its computation has         
!               been requested.  here trans(q) is the                   
!               transpose of the matrix q.                              
!                                                                       
!        b      double precision(k)                                     
!               b contains the solution of the least squares problem    
!                                                                       
!                    minimize norm2(y - xk*b),                          
!                                                                       
!               if its computation has been requested.  (note that      
!               if pivoting was requested in dqrdc, the j-th            
!               component of b will be associated with column jpvt(j)   
!               of the original matrix x that was input into dqrdc.)    
!                                                                       
!        rsd    double precision(n).                                    
!               rsd contains the least squares residual y - xk*b,       
!               if its computation has been requested.  rsd is          
!               also the orthogonal projection of y onto the            
!               orthogonal complement of the column space of xk.        
!                                                                       
!        xb     double precision(n).                                    
!               xb contains the least squares approximation xk*b,       
!               if its computation has been requested.  xb is also      
!               the orthogonal projection of y onto the column space    
!               of x.                                                   
!                                                                       
!        info   integer.                                                
!               info is zero unless the computation of b has            
!               been requested and r is exactly singular.  in           
!               this case, info is the index of the first zero          
!               diagonal element of r and b is left unaltered.          
!                                                                       
!     the parameters qy, qty, b, rsd, and xb are not referenced         
!     if their computation is not requested and in this case            
!     can be replaced by dummy variables in the calling program.        
!     to save storage, the user may in some cases use the same          
!     array for different parameters in the calling sequence.  a        
!     frequently occuring example is when one wishes to compute         
!     any of b, rsd, or xb and does not need y or qty.  in this         
!     case one may identify y, qty, and one of b, rsd, or xb, while     
!     providing separate arrays for anything else that is to be         
!     computed.  thus the calling sequence                              
!                                                                       
!          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)         
!                                                                       
!     will result in the computation of b and rsd, with rsd             
!     overwriting y.  more generally, each item in the following        
!     list contains groups of permissible identifications for           
!     a single callinng sequence.                                       
!                                                                       
!          1. (y,qty,b) (rsd) (xb) (qy)                                 
!                                                                       
!          2. (y,qty,rsd) (b) (xb) (qy)                                 
!                                                                       
!          3. (y,qty,xb) (b) (rsd) (qy)                                 
!                                                                       
!          4. (y,qy) (qty,b) (rsd) (xb)                                 
!                                                                       
!          5. (y,qy) (qty,rsd) (b) (xb)                                 
!                                                                       
!          6. (y,qy) (qty,xb) (b) (rsd)                                 
!                                                                       
!     in any group the value returned in the array allocated to         
!     the group corresponds to the last member of the group.            
!                                                                       
!     linpack. this version dated 08/14/78 .                            
!     g.w. stewart, university of maryland, argonne national lab.       
!                                                                       
!     dqrsl uses the following functions and subprograms.               
!                                                                       
!     blas daxpy,dcopy,ddot_long                                             
!     fortran dabs,min0,mod                                             
!                                                                       
!     internal variables                                                
!                                                                       
      integer i,j,jj,ju,kp1 
      double precision ddot_long,t,temp 
      logical cb,cqy,cqty,cr,cxb 
!                                                                       
!                                                                       
!     set info flag.                                                    
!                                                                       
      info = 0 
!                                                                       
!     determine what is to be computed.                                 
!                                                                       
      cqy = job/10000 .ne. 0 
      cqty = mod(job,10000) .ne. 0 
      cb = mod(job,1000)/100 .ne. 0 
      cr = mod(job,100)/10 .ne. 0 
      cxb = mod(job,10) .ne. 0 
      ju = min0(k,n-1) 
!                                                                       
!     special action when n=1.                                          
!                                                                       
      if (ju .ne. 0) go to 40 
         if (cqy) qy(1) = y(1) 
         if (cqty) qty(1) = y(1) 
         if (cxb) xb(1) = y(1) 
         if (.not.cb) go to 30 
            if (x(1,1) .ne. 0.0d0) go to 10 
               info = 1 
            go to 20 
   10       continue 
               b(1) = y(1)/x(1,1) 
   20       continue 
   30    continue 
         if (cr) rsd(1) = 0.0d0 
      go to 250 
   40 continue 
!                                                                       
!        set up to compute qy or qty.                                   
!                                                                       
         if (cqy) call dcopy(n,y,1,qy,1) 
         if (cqty) call dcopy(n,y,1,qty,1) 
         if (.not.cqy) go to 70 
!                                                                       
!           compute qy.                                                 
!                                                                       
            do 60 jj = 1, ju 
               j = ju - jj + 1 
               if (qraux(j) .eq. 0.0d0) go to 50 
                  temp = x(j,j) 
                  x(j,j) = qraux(j) 
                  t = -ddot_long(int8(n-j+1),x(j,j),1,qy(j),1)/x(j,j) 
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1) 
                  x(j,j) = temp 
   50          continue 
   60       continue 
   70    continue 
         if (.not.cqty) go to 100 
!                                                                       
!           compute trans(q)*y.                                         
!                                                                       
            do 90 j = 1, ju 
               if (qraux(j) .eq. 0.0d0) go to 80 
                  temp = x(j,j) 
                  x(j,j) = qraux(j) 
                  t = -ddot_long(int8(n-j+1),x(j,j),1,qty(j),1)/x(j,j) 
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1) 
                  x(j,j) = temp 
   80          continue 
   90       continue 
  100    continue 
!                                                                       
!        set up to compute b, rsd, or xb.                               
!                                                                       
         if (cb) call dcopy(k,qty,1,b,1) 
         kp1 = k + 1 
         if (cxb) call dcopy(k,qty,1,xb,1) 
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1) 
         if (.not.cxb .or. kp1 .gt. n) go to 120 
            do 110 i = kp1, n 
               xb(i) = 0.0d0 
  110       continue 
  120    continue 
         if (.not.cr) go to 140 
            do 130 i = 1, k 
               rsd(i) = 0.0d0 
  130       continue 
  140    continue 
         if (.not.cb) go to 190 
!                                                                       
!           compute b.                                                  
!                                                                       
            do 170 jj = 1, k 
               j = k - jj + 1 
               if (x(j,j) .ne. 0.0d0) go to 150 
                  info = j 
!           ......exit                                                  
                  go to 180 
  150          continue 
               b(j) = b(j)/x(j,j) 
               if (j .eq. 1) go to 160 
                  t = -b(j) 
                  call daxpy(j-1,t,x(1,j),1,b,1) 
  160          continue 
  170       continue 
  180       continue 
  190    continue 
         if (.not.cr .and. .not.cxb) go to 240 
!                                                                       
!           compute rsd or xb as required.                              
!                                                                       
            do 230 jj = 1, ju 
               j = ju - jj + 1 
               if (qraux(j) .eq. 0.0d0) go to 220 
                  temp = x(j,j) 
                  x(j,j) = qraux(j) 
                  if (.not.cr) go to 200 
                     t = -ddot_long(int8(n-j+1),x(j,j),1,rsd(j),1)/x(j,j) 
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1) 
  200             continue 
                  if (.not.cxb) go to 210 
                     t = -ddot_long(int8(n-j+1),x(j,j),1,xb(j),1)/x(j,j) 
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1) 
  210             continue 
                  x(j,j) = temp 
  220          continue 
  230       continue 
  240    continue 
  250 continue 
      return 
      END                                           
