      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job) 
      integer ldx,n,p,job 
      integer jpvt(1) 
      double precision x(ldx,1),qraux(1),work(1) 
!                                                                       
!     dqrdc uses householder transformations to compute the qr          
!     factorization of an n by p matrix x.  column pivoting             
!     based on the 2-norms of the reduced columns may be                
!     performed at the users option.                                    
!                                                                       
!     on entry                                                          
!                                                                       
!        x       double precision(ldx,p), where ldx .ge. n.             
!                x contains the matrix whose decomposition is to be     
!                computed.                                              
!                                                                       
!        ldx     integer.                                               
!                ldx is the leading dimension of the array x.           
!                                                                       
!        n       integer.                                               
!                n is the number of rows of the matrix x.               
!                                                                       
!        p       integer.                                               
!                p is the number of columns of the matrix x.            
!                                                                       
!        jpvt    integer(p).                                            
!                jpvt contains integers that control the selection      
!                of the pivot columns.  the k-th column x(k) of x       
!                is placed in one of three classes according to the     
!                value of jpvt(k).                                      
!                                                                       
!                   if jpvt(k) .gt. 0, then x(k) is an initial          
!                                      column.                          
!                                                                       
!                   if jpvt(k) .eq. 0, then x(k) is a free column.      
!                                                                       
!                   if jpvt(k) .lt. 0, then x(k) is a final column.     
!                                                                       
!                before the decomposition is computed, initial columns  
!                are moved to the beginning of the array x and final    
!                columns to the end.  both initial and final columns    
!                are frozen in place during the computation and only    
!                free columns are moved.  at the k-th stage of the      
!                reduction, if x(k) is occupied by a free column        
!                it is interchanged with the free column of largest     
!                reduced norm.  jpvt is not referenced if               
!                job .eq. 0.                                            
!                                                                       
!        work    double precision(p).                                   
!                work is a work array.  work is not referenced if       
!                job .eq. 0.                                            
!                                                                       
!        job     integer.                                               
!                job is an integer that initiates column pivoting.      
!                if job .eq. 0, no pivoting is done.                    
!                if job .ne. 0, pivoting is done.                       
!                                                                       
!     on return                                                         
!                                                                       
!        x       x contains in its upper triangle the upper             
!                triangular matrix r of the qr factorization.           
!                below its diagonal x contains information from         
!                which the orthogonal part of the decomposition         
!                can be recovered.  note that if pivoting has           
!                been requested, the decomposition is not that          
!                of the original matrix x but that of x                 
!                with its columns permuted as described by jpvt.        
!                                                                       
!        qraux   double precision(p).                                   
!                qraux contains further information required to recover 
!                the orthogonal part of the decomposition.              
!                                                                       
!        jpvt    jpvt(k) contains the index of the column of the        
!                original matrix that has been interchanged into        
!                the k-th column, if pivoting was requested.            
!                                                                       
!     linpack. this version dated 08/14/78 .                            
!     g.w. stewart, university of maryland, argonne national lab.       
!                                                                       
!     dqrdc uses the following functions and subprograms.               
!                                                                       
!     blas daxpy,ddot_long,dscal,dswap,dnrm2_long                                 
!     fortran dabs,dmax1,min0,dsqrt                                     
!                                                                       
!     internal variables                                                
!                                                                       
      integer j,jj,jp,l,lp1,lup,maxj,pl,pu 
      double precision maxnrm,dnrm2_long,tt 
      double precision ddot_long,nrmxl,t 
      logical negj,swapj 
!                                                                       
!                                                                       
      pl = 1 
      pu = 0 
      if (job .eq. 0) go to 60 
!                                                                       
!        pivoting has been requested.  rearrange the columns            
!        according to jpvt.                                             
!                                                                       
         do 20 j = 1, p 
            swapj = jpvt(j) .gt. 0 
            negj = jpvt(j) .lt. 0 
            jpvt(j) = j 
            if (negj) jpvt(j) = -j 
            if (.not.swapj) go to 10 
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1) 
               jpvt(j) = jpvt(pl) 
               jpvt(pl) = j 
               pl = pl + 1 
   10       continue 
   20    continue 
         pu = p 
         do 50 jj = 1, p 
            j = p - jj + 1 
            if (jpvt(j) .ge. 0) go to 40 
               jpvt(j) = -jpvt(j) 
               if (j .eq. pu) go to 30 
                  call dswap(n,x(1,pu),1,x(1,j),1) 
                  jp = jpvt(pu) 
                  jpvt(pu) = jpvt(j) 
                  jpvt(j) = jp 
   30          continue 
               pu = pu - 1 
   40       continue 
   50    continue 
   60 continue 
!                                                                       
!     compute the norms of the free columns.                            
!                                                                       
      if (pu .lt. pl) go to 80 
      do 70 j = pl, pu 
         qraux(j) = dnrm2_long(int8(n),x(1,j),1) 
         work(j) = qraux(j) 
   70 continue 
   80 continue 
!                                                                       
!     perform the householder reduction of x.                           
!                                                                       
      lup = min0(n,p) 
      do 200 l = 1, lup 
         if (l .lt. pl .or. l .ge. pu) go to 120 
!                                                                       
!           locate the column of largest norm and bring it              
!           into the pivot position.                                    
!                                                                       
            maxnrm = 0.0d0 
            maxj = l 
            do 100 j = l, pu 
               if (qraux(j) .le. maxnrm) go to 90 
                  maxnrm = qraux(j) 
                  maxj = j 
   90          continue 
  100       continue 
            if (maxj .eq. l) go to 110 
               call dswap(n,x(1,l),1,x(1,maxj),1) 
               qraux(maxj) = qraux(l) 
               work(maxj) = work(l) 
               jp = jpvt(maxj) 
               jpvt(maxj) = jpvt(l) 
               jpvt(l) = jp 
  110       continue 
  120    continue 
         qraux(l) = 0.0d0 
         if (l .eq. n) go to 190 
!                                                                       
!           compute the householder transformation for column l.        
!                                                                       
            nrmxl = dnrm2_long(int8(n-l+1),x(l,l),1) 
            if (nrmxl .eq. 0.0d0) go to 180 
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l)) 
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1) 
               x(l,l) = 1.0d0 + x(l,l) 
!                                                                       
!              apply the transformation to the remaining columns,       
!              updating the norms.                                      
!                                                                       
               lp1 = l + 1 
               if (p .lt. lp1) go to 170 
               do 160 j = lp1, p 
                  t = -ddot_long(int8(n-l+1),x(l,l),1,x(l,j),1)/x(l,l) 
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1) 
                  if (j .lt. pl .or. j .gt. pu) go to 150 
                  if (qraux(j) .eq. 0.0d0) go to 150 
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2 
                     tt = dmax1(tt,0.0d0) 
                     t = tt 
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2 
                     if (tt .eq. 1.0d0) go to 130 
                        qraux(j) = qraux(j)*dsqrt(t) 
                     go to 140 
  130                continue 
                        qraux(j) = dnrm2_long(int8(n-l),x(l+1,j),1) 
                        work(j) = qraux(j) 
  140                continue 
  150             continue 
  160          continue 
  170          continue 
!                                                                       
!              save the transformation.                                 
!                                                                       
               qraux(l) = x(l,l) 
               x(l,l) = -nrmxl 
  180       continue 
  190    continue 
  200 continue 
      return 
      END                                           
