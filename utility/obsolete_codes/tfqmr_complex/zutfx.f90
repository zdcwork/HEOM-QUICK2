!                                                                       
! modified on 18 Aug, 2011:                                             
! change data type from 'integer' to 'integer * 8' for nlen and ndim    
!                                                                       
!********************************************************************** 
!                                                                       
!     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal        
!     All rights reserved.                                              
!                                                                       
!     This code is part of a copyrighted package.  For details, see the 
!     file "cpyrit.doc" in the top-level directory.                     
!                                                                       
!     ***************************************************************** 
!     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE 
!                             COPYRIGHT NOTICE                          
!     ***************************************************************** 
!                                                                       
!********************************************************************** 
!                                                                       
!     This file contains the routine for the TFQMR algorithm.           
!                                                                       
!********************************************************************** 
!                                                                       
      SUBROUTINE ZUTFX (NDIM,NLEN,NLIM,VECS,TOL,INFO) 
!                                                                       
!     Purpose:                                                          
!     This subroutine uses the TFQMR algorithm to solve linear systems. 
!     It runs the  algorithm to convergence  or until  a user-specified 
!     limit on the number of iterations is reached.                     
!                                                                       
!     The  code is  set up  to solve  the system  A x = b  with initial 
!     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system, 
!     and  it  is  connected with the  original system as follows.  Let 
!     B y = c be the original unpreconditioned system to be solved, and 
!     let y_0 be an arbitrary initial guess for its solution.  Then:    
!          A x = b, where  A = M_1^{-1} B M_2^{-1},                     
!                          x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0). 
!     Here M = M_1 M_2 is the preconditioner.                           
!                                                                       
!     To recover the final iterate y_n for  the original system B y = c 
!     from the final iterate x_n for the preconditioned system A x = b, 
!     set                                                               
!               y_n = y_0 + M_2^{-1} x_n.                               
!                                                                       
!     The algorithm was first described  in the RIACS Technical Report  
!     91.18, `A Transpose-Free  Quasi-Minimal  Residual Algorithm  for  
!     Non-Hermitian Linear Systems`, by Roland Freund, September 1991,  
!     which subsequently appeared in  SIAM J. Sci. Comput., 14 (1993),  
!     pp. 470--482.                                                     
!                                                                       
!     Parameters:                                                       
!     For a description of  the parameters, see the file "zutfx.doc" in 
!     the current directory.                                            
!                                                                       
!     External routines used:                                           
!     double precision dlamch(ch)                                       
!        LAPACK routine, computes machine-related constants.            
!     double precision dznrm2(n,x,incx)                                 
!        BLAS-1 routine, computes the 2-norm of x.                      
!     subroutine zaxpby(n,z,a,x,b,y)                                    
!        Library routine, computes z = a * x + b * y.                   
!     double precision zdotu(n,x,incx,y,incy)                           
!        BLAS-1 routine, computes y^H * x.                              
!     subroutine zrandn(n,x,seed)                                       
!        Library routine, fills x with random numbers.                  
!                                                                       
!     Noel M. Nachtigal                                                 
!     April 13, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC CDABS, DBLE, DSQRT, MAX0 
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN 
      DOUBLE COMPLEX ZDOTU 
      DOUBLE PRECISION DLAMCH, DZNRM2 
!                                                                       
      INTEGER INFO(4), NLIM 
      INTEGER*8 NDIM, NLEN 
      DOUBLE COMPLEX VECS(NDIM,9) 
      DOUBLE PRECISION TOL 
!                                                                       
!     Miscellaneous parameters.                                         
!                                                                       
      DOUBLE COMPLEX ZONE, ZZERO 
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0)) 
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO 
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0) 
!                                                                       
!     Local variables, permanent.                                       
!                                                                       
      INTEGER IERR, N, RETLBL, TF, TRES, VF 
      SAVE    IERR, N, RETLBL, TF, TRES, VF 
      DOUBLE COMPLEX ALPHA, BETA, ETA, RHO 
      SAVE           ALPHA, BETA, ETA, RHO 
      DOUBLE PRECISION COS, VAR, R0, RESN, TAU, UCHK, UNRM 
      SAVE             COS, VAR, R0, RESN, TAU, UCHK, UNRM 
!                                                                       
!     Local variables, transient.                                       
!                                                                       
      INTEGER INIT, REVCOM 
      DOUBLE COMPLEX ZTMP 
      DOUBLE PRECISION DTMP 
!                                                                       
!     Initialize some of the permanent variables.                       
!                                                                       
      DATA RETLBL /0/ 
!                                                                       
!     Check the reverse communication flag to see where to branch.      
!        REVCOM   RETLBL      Comment                                   
!           0        0    first call, go to label 10                    
!           1       30    returning from AXB, go to label 30            
!           1       40    returning from AXB, go to label 40            
!           1       60    returning from AXB, go to label 60            
!           1       70    returning from AXB, go to label 70            
!                                                                       
      REVCOM  = INFO(2) 
      INFO(2) = 0 
      IF (REVCOM.EQ.0) THEN 
         N = 0 
         IF (RETLBL.EQ.0) GO TO 10 
      ELSE IF (REVCOM.EQ.1) THEN 
         IF (RETLBL.EQ.30) THEN 
            GO TO 30 
         ELSE IF (RETLBL.EQ.40) THEN 
            GO TO 40 
         ELSE IF (RETLBL.EQ.60) THEN 
            GO TO 60 
         ELSE IF (RETLBL.EQ.70) THEN 
            GO TO 70 
         END IF 
      END IF 
      IERR = 1 
      GO TO 90 
!                                                                       
!     Check whether the inputs are valid.                               
!                                                                       
   10 IERR = 0 
      IF (NDIM.LT.1)    IERR = 2 
      IF (NLEN.LT.1)    IERR = 2 
      IF (NLIM.LT.1)    IERR = 2 
      IF (NLEN.GT.NDIM) IERR = 2 
      IF (IERR.NE.0) GO TO 90 
!                                                                       
!     Extract from INFO the output units TF and VF, the true residual   
!     flag TRES, and the left starting vector flag INIT.                
!                                                                       
      VF   = MAX0(INFO(1),0) 
      INIT = VF / 100000 
      VF   = VF - INIT * 100000 
      TRES = VF / 10000 
      VF   = VF - TRES * 10000 
      TF   = VF / 100 
      VF   = VF - TF * 100 
!                                                                       
!     Check the convergence tolerance.                                  
!                                                                       
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E')) 
!                                                                       
!     Start the trace messages and convergence history.                 
!                                                                       
      IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') 0, 0, DONE, DONE 
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') 0, 0, DONE, DONE 
!                                                                       
!     Set x_0 = 0 and compute the norm of the initial residual.         
!                                                                       
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,2),ZZERO,VECS(1,5)) 
      CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1)) 
      R0 = DZNRM2(NLEN,VECS(1,5),1) 
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 90 
!                                                                       
!     Check whether the auxiliary vector must be supplied.              
!                                                                       
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,3),1) 
!                                                                       
!     Initialize the variables.                                         
!                                                                       
      N    = 1 
      RESN = DONE 
      RHO  = ZONE 
      VAR  = DZERO 
      ETA  = ZZERO 
      TAU  = R0 * R0 
      IERR = 8 
      CALL ZAXPBY (NLEN,VECS(1,8),ZZERO,VECS(1,8),ZZERO,VECS(1,8)) 
      CALL ZAXPBY (NLEN,VECS(1,4),ZZERO,VECS(1,4),ZZERO,VECS(1,4)) 
      CALL ZAXPBY (NLEN,VECS(1,6),ZZERO,VECS(1,6),ZZERO,VECS(1,6)) 
!                                                                       
!     This is one step of the TFQMR algorithm.                          
!     Compute \beta_{n-1} and \rho_{n-1}.                               
!                                                                       
   20 ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,5),1) 
      BETA = ZTMP / RHO 
      RHO  = ZTMP 
!                                                                       
!     Compute y_{2n-1}, v_{n-1}, and A y_{2n-1}.                        
!                                                                       
      CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,8)) 
      CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,5),BETA,VECS(1,6)) 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,6),VECS(1,9))                                 
!                                                                       
      INFO(2) = 1 
      INFO(3) = 6 
      INFO(4) = 9 
      RETLBL  = 30 
      RETURN 
   30 CALL ZAXPBY (NLEN,VECS(1,4),BETA,VECS(1,4),ZONE,VECS(1,9)) 
!                                                                       
!     Compute \sigma{n-1} and check for breakdowns.                     
!                                                                       
      ZTMP = ZDOTU(NLEN,VECS(1,3),1,VECS(1,4),1) 
      IF ((CDABS(ZTMP).EQ.DZERO).OR.(CDABS(RHO).EQ.DZERO)) THEN 
         IERR = 8 
         GO TO 90 
      END IF 
!                                                                       
!     Compute \alpha_{n-1}, d_{2n-1} and w_{2n}.                        
!                                                                       
      ALPHA = RHO / ZTMP 
      ZTMP  = VAR * ETA / ALPHA 
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7)) 
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,9)) 
!                                                                       
!     Compute \varepsilon_{2n-1}^2, \eta_{2n-1}^2, c_{2n-1}^2, and      
!     \tau_{2n-1}^2.                                                    
!                                                                       
      DTMP = DZNRM2(NLEN,VECS(1,5),1) 
      DTMP = DTMP * DTMP 
      VAR  = DTMP / TAU 
      COS  = DONE / ( DONE + VAR ) 
      TAU  = DTMP * COS 
      ETA  = ALPHA * COS 
!                                                                       
!     Compute x_{2n-1} and the upper bound for its residual norm.       
!                                                                       
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7)) 
!                                                                       
!     Compute the residual norm upper bound.                            
!     If the scaled upper bound is within one order of magnitude of the 
!     target convergence norm, compute the true residual norm.          
!                                                                       
      UNRM = DSQRT(DBLE(2*N) * TAU) / R0 
      UCHK = UNRM 
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN)) GO TO 50 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,1),VECS(1,9))                                 
!                                                                       
      INFO(2) = 1 
      INFO(3) = 1 
      INFO(4) = 9 
      RETLBL  = 40 
      RETURN 
   40 CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9)) 
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0 
      UCHK = RESN 
!                                                                       
!     Output the trace messages and convergence history.                
!                                                                       
   50 IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN 
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N-1, UNRM, RESN 
!                                                                       
!     Check for convergence or termination.  Stop if:                   
!         1. algorithm converged;                                       
!         2. the residual norm upper bound is smaller than the computed 
!     residual norm by a factor of at least 100.                        
!                                                                       
      IF (RESN.LE.TOL) THEN 
         IERR = 0 
         GO TO 90 
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN 
         IERR = 4 
         GO TO 90 
      END IF 
!                                                                       
!     Compute y_{2n}, A y_{2n}, d_{2n}, and w_{2n+1}.                   
!                                                                       
      CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,6),-ALPHA,VECS(1,4)) 
      ZTMP = VAR * COS 
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,6),ZTMP,VECS(1,7)) 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,6),VECS(1,8))                                 
!                                                                       
      INFO(2) = 1 
      INFO(3) = 6 
      INFO(4) = 8 
      RETLBL  = 60 
      RETURN 
   60 CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,5),-ALPHA,VECS(1,8)) 
!                                                                       
!     Compute \varepsilon_{2n}^2, \eta_{2n}^2, c_{2n}^2, and            
!     \tau_{2n}^2.                                                      
!                                                                       
      DTMP = DZNRM2(NLEN,VECS(1,5),1) 
      DTMP = DTMP * DTMP 
      VAR  = DTMP / TAU 
      COS  = DONE / ( DONE + VAR ) 
      TAU  = DTMP * COS 
      ETA  = ALPHA * COS 
!                                                                       
!     Compute x_{2n}.                                                   
!                                                                       
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ETA,VECS(1,7)) 
!                                                                       
!     Compute the residual norm upper bound.                            
!     If the scaled upper bound is within one order of magnitude of the 
!     target convergence norm, compute the true residual norm.          
!                                                                       
      UNRM = DSQRT(DBLE(2*N+1) * TAU) / R0 
      UCHK = UNRM 
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 80 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,1),VECS(1,9))                                 
!                                                                       
      INFO(2) = 1 
      INFO(3) = 1 
      INFO(4) = 9 
      RETLBL  = 70 
      RETURN 
   70 CALL ZAXPBY (NLEN,VECS(1,9),ZONE,VECS(1,2),-ZONE,VECS(1,9)) 
      RESN = DZNRM2(NLEN,VECS(1,9),1) / R0 
      UCHK = UNRM 
!                                                                       
!     Output the trace messages and convergence history.                
!                                                                       
   80 IF (VF.NE.0) WRITE (VF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN 
      IF (TF.NE.0) WRITE (TF,'(2I8,2E11.4)') N, 2*N, UNRM, RESN 
!                                                                       
!     Check for convergence or termination.  Stop if:                   
!         1. algorithm converged;                                       
!         2. the residual norm upper bound is smaller than the computed 
!     residual norm by a factor of at least 100;                        
!         3. algorithm exceeded the iterations limit.                   
!                                                                       
      IF (RESN.LE.TOL) THEN 
         IERR = 0 
         GO TO 90 
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN 
         IERR = 4 
         GO TO 90 
      ELSE IF (N.GE.NLIM) THEN 
         IERR = 4 
         GO TO 90 
      END IF 
!                                                                       
!     Update the running counter.                                       
!                                                                       
      N = N + 1 
      GO TO 20 
!                                                                       
!     Done.                                                             
!                                                                       
   90 NLIM    = N 
      RETLBL  = 0 
      INFO(1) = IERR 
!                                                                       
      RETURN 
      END                                           
