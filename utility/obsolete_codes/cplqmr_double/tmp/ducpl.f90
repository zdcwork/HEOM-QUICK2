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
!     This  file  contains  the  routines  for  the  QMR  algorithm for 
!     unsymmetric  matrices,  using  the  coupled  two-term  recurrence 
!     variant of the look-ahead Lanczos algorithm.                      
!                                                                       
!********************************************************************** 
!                                                                       
      SUBROUTINE DUCPL (NDIM,NLEN,NLIM,MAXPQ,MAXVW,M,MVEC,NORMS,        &
     &                  DWK,IDX,IWK,VECS,TOL,INFO)                      
!                                                                       
!     Purpose:                                                          
!     This subroutine uses the QMR algorithm based  on the coupled two- 
!     term  variant of  the look-ahead Lanczos  process to solve linear 
!     systems.  It runs the algorithm  to convergence or  until a user- 
!     specified limit on the number of iterations is reached.           
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
!     The algorithm  was first described  in the RIACS Technical Report 
!     92.15, `An Implementation of the QMR Method Based on Coupled Two- 
!     Term Recurrences`,  June 1992.  A good reference  for the details 
!     of  this implementation  is the  RIACS  Technical  Report  92.19, 
!     `Implementation  details of the coupled  QMR algorithm`,  by R.W. 
!     Freund and N.M. Nachtigal, October 1992.                          
!                                                                       
!     Parameters:                                                       
!     For a description of  the parameters, see the file "ducpl.doc" in 
!     the current directory.                                            
!                                                                       
!     External routines used:                                           
!     subroutine daxpby(n,z,a,x,b,y)                                    
!        Library routine, computes z = a * x + b * y.                   
!     double precision ddot(n,x,incx,y,incy)                            
!        BLAS-1 routine, computes y^H * x.                              
!     double precision dlamch(ch)                                       
!        LAPACK routine, computes machine-related constants.            
!     double precision dnrm2(n,x,incx)                                  
!        BLAS-1 routine, computes the 2-norm of x.                      
!     subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)                   
!        LINPACK routine, computes the QR factorization of x.           
!     subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)      
!        LINPACK routine, applies the QR factorization of x.            
!     subroutine drandn(n,x,seed)                                       
!        Library routine, fill x with random numbers.                   
!     subroutine drotg(a,b,cos,sin)                                     
!        BLAS-1 routine, computes the Givens rotation which rotates the 
!        vector [a; b] into [ sqrt(a**2 + b**2); 0 ].                   
!     subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)         
!        LINPACK routine, computes the SVD of x.                        
!     subroutine ducpl1(ndim,nlen,m,n,k,kstar,l,lstar,mk,mkstar,nl,     
!        nlstar,vf,ierr,adjust,norm1,norm2,dwk,iwk,idx,vecs)            
!        Low-level routine, rebuilds vectors PQ in CPL.                 
!     subroutine ducpl2(ndim,nlen,m,n,k,kstar,l,lstar,mk,mkstar,nl,     
!        nlstar,vf,ierr,adjust,norm1,norm2,dwk,iwk,idx,vecs)            
!        Low-level routine, rebuilds vectors VW in CPL.                 
!     double precision function ducpll(i,n)                             
!        User-supplied routine, computes inner recurrence coefficients. 
!     double precision ducplom(n)                                       
!        User-supplied routine, specifies the QMR scaling factors.      
!     double precision function ducplu(i,n)                             
!        User-supplied routine, computes inner recurrence coefficients. 
!                                                                       
!     Noel M. Nachtigal                                                 
!     March 1, 1992                                                     
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC DABS, DBLE, DMAX1, DSQRT, MAX0, MIN0, MOD 
      EXTERNAL DAXPBY, DDOT, DLAMCH, DNRM2, DQRDC, DQRSL, DRANDN, DROTG 
      EXTERNAL DSVDC, DUCPL1, DUCPL2, DUCPLL, DUCPLO, DUCPLU 
      DOUBLE PRECISION DDOT, DLAMCH, DNRM2, DUCPLL, DUCPLO, DUCPLU 
!                                                                       
!      INTEGER INFO(4), M, MAXPQ, MAXVW, MVEC, NDIM, NLEN, NLIM 
      INTEGER INFO(4), M, MAXPQ, MAXVW, MVEC, NLIM 
      integer*8 ndim, nlen
      INTEGER IDX(6,NLIM+2), IWK(M,13) 
      DOUBLE PRECISION DWK(M,8*M+18), NORMS(2), TOL 
      DOUBLE PRECISION VECS(NDIM,5*MVEC+4-1) 
!                                                                       
!     Common block variables.                                           
!                                                                       
!                                                                       
!     Common block DUCPLX.                                              
!                                                                       
      DOUBLE PRECISION NORMA 
      COMMON /DUCPLX/NORMA 
!                                                                       
!     Miscellaneous parameters.                                         
!                                                                       
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO 
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0) 
!                                                                       
!     Local variables, permanent.                                       
!                                                                       
      LOGICAL IBUILT, INNER, RERUN 
      SAVE    IBUILT, INNER, RERUN 
      INTEGER IEND, IERR, K, KBLKSZ, KSTAR, L, LBLKSZ, LSTAR, MK 
      SAVE    IEND, IERR, K, KBLKSZ, KSTAR, L, LBLKSZ, LSTAR, MK 
      INTEGER MKMAX, MKSTAR, MPQBLT, MVWBLT, N, NL, NLMAX, NLSTAR, NMAX 
      SAVE    MKMAX, MKSTAR, MPQBLT, MVWBLT, N, NL, NLMAX, NLSTAR, NMAX 
      INTEGER NP1, NQMR, NUMCHK, PQBEG, RETLBL, TF, TRES, VF, VWBEG 
      SAVE    NP1, NQMR, NUMCHK, PQBEG, RETLBL, TF, TRES, VF, VWBEG 
      DOUBLE PRECISION ADJUST, MAXOMG, NORM1, NORM2, TMAX, TMIN, TNRM 
      SAVE             ADJUST, MAXOMG, NORM1, NORM2, TMAX, TMIN, TNRM 
      DOUBLE PRECISION R0, RESN, UCHK, UNRM 
      SAVE             R0, RESN, UCHK, UNRM 
!                                                                       
!     Local variables, transient.                                       
!                                                                       
      INTEGER I, IBASE, IBLKSZ, IJ, INIT, J, KEND, LEND, MI, MIP1 
      INTEGER NI, NIP1, REVCOM 
      DOUBLE PRECISION DTMP, DTMP1, DTMP2, DTMP3, DTMP4, DTMP5 
      DOUBLE PRECISION IDTMP1, IDTMP2, IDTMP3, IDTMP4 
!                                                                       
!     Initialize some of the permanent variables.                       
!                                                                       
      DATA RETLBL /0/ 
!                                                                       
!     Check the reverse communication flag to see where to branch.      
!        REVCOM   RETLBL      Comment                                   
!           0        0    first call, go to label 10                    
!           1      200    returning from AXB, go to label 200           
!           1      350    returning from AXB, go to label 350           
!           1      720    returning from AXB, go to label 720           
!           2      370    returning from ATXB, go to label 370          
!                                                                       
write(6,*)'ducpl: here1'
call flush(6)
      REVCOM  = INFO(2) 
      INFO(2) = 0 
write(6,*)'ducpl: here2, revcom, retlbl', REVCOM, RETLBL
call flush(6)
      IF (REVCOM.EQ.0) THEN 
         N      = 0 
         NQMR   = 0 
         MPQBLT = 0 
         MVWBLT = 0 
         IF (RETLBL.EQ.0) GO TO 10 
      ELSE IF (REVCOM.EQ.1) THEN 
         IF (RETLBL.EQ.200) THEN 
            GO TO 200 
         ELSE IF (RETLBL.EQ.350) THEN 
            GO TO 350 
         ELSE IF (RETLBL.EQ.720) THEN 
            GO TO 720 
         END IF 
      ELSE IF (REVCOM.EQ.2) THEN 
         IF (RETLBL.EQ.370) GO TO 370 
      END IF 
      IERR = 1 
      GO TO 750 
!                                                                       
!     Check whether the inputs are valid.                               
!                                                                       
   10 IERR = 0 
      IF (NDIM.LT.1)          IERR = 2 
      IF (NLEN.LT.1)          IERR = 2 
      IF (NLIM.LT.1)          IERR = 2 
      IF (MAXPQ.LT.1)         IERR = 2 
      IF (MAXVW.LT.1)         IERR = 2 
      IF (NLEN.GT.NDIM)       IERR = 2 
      IF (MVEC.LT.MAXPQ+1)    IERR = 2 
      IF (MVEC.LT.MAXVW+1)    IERR = 2 
      IF (M.LT.MAXVW+MAXPQ+2) IERR = 2 
      IF (IERR.NE.0) GO TO 750 
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
!     Extract the norms.                                                
!                                                                       
      NORMA = DZERO 
      NORM1 = DABS(NORMS(1)) 
      NORM2 = DABS(NORMS(2)) 
!                                                                       
!     Set the adjustment parameters.                                    
!                                                                       
      NUMCHK = 25 
      ADJUST = DTEN 
!                                                                       
!     Extract and check the various tolerances and norm estimates.      
!                                                                       
      TNRM = DLAMCH('E') * DTEN 
      TMIN = DSQRT(DSQRT(DLAMCH('S'))) 
      TMAX = DONE / TMIN 
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E')) 
!                                                                       
!     Start the trace messages and convergence history.                 
!                                                                       
      IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') 0, DONE, DONE 
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') 0, DONE, DONE 
!                                                                       
!     Set up wrapped indices.  The following indices are used:          
!        IDX(1,I) = indices used for the work arrays row-dimensioned M; 
!        IDX(2,I) = indices used for v_i;                               
!        IDX(3,I) = indices used for w_i;                               
!        IDX(4,I) = indices used for p_i;                               
!        IDX(5,I) = indices used for q_i;                               
!        IDX(6,I) = indices used for s_i.                               
!                                                                       
      DO 20 I = 1, NLIM+2 
         IDX(1,I) = MOD(I-1,M) + 1 
         IDX(2,I) = 4 + 1 + 0*MVEC + MOD(I-1,MVEC) 
         IDX(3,I) = 4 + 1 + 1*MVEC + MOD(I-1,MVEC) 
         IDX(4,I) = 4 + 1 + 2*MVEC + MOD(I-1,MVEC) 
         IDX(5,I) = 4 + 1 + 3*MVEC + MOD(I-1,MVEC) 
         IDX(6,I) = 4 + 1 + 4*MVEC + MOD(I-1,MVEC-1) 
   20 END DO 
!                                                                       
!     Set x_0 = 0 and compute the norm of the initial residual.         
!                                                                       
      CALL DAXPBY (NLEN,VECS(1,IDX(2,1)),DONE,VECS(1,2),DZERO,VECS(1,   &
     &IDX(2,1)))                                                        
      CALL DAXPBY (NLEN,VECS(1,1),DZERO,VECS(1,1),DZERO,VECS(1,1)) 
      R0 = DNRM2(NLEN,VECS(1,IDX(2,1)),1) 
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 750 
!                                                                       
!     Check whether the auxiliary vector must be supplied.              
!                                                                       
      IF (INIT.EQ.0) CALL DRANDN (NLEN,VECS(1,3),1) 
      CALL DAXPBY (NLEN,VECS(1,IDX(3,1)),DONE,VECS(1,3),DZERO,VECS(1,   &
     &IDX(3,1)))                                                        
!                                                                       
!     Scale the first pair of Lanczos vectors and check for invariant   
!     subspaces.                                                        
!                                                                       
      DTMP1 = R0 
      DTMP2 = DNRM2(NLEN,VECS(1,IDX(3,1)),1) 
      IF (DTMP1.LT.TNRM) IERR = IERR + 16 
      IF (DTMP2.LT.TNRM) IERR = IERR + 32 
      IF (IERR.NE.0) GO TO 750 
      DWK(IDX(1,1),8*M+10)   = DONE 
      DWK(IDX(1,1),IDX(1,1)) = DDOT(NLEN,VECS(1,IDX(2,1)),1,VECS(1,IDX(3&
     &,1)),1) / ( DTMP1 * DTMP2 )                                       
      IF ((DTMP1.GE.TMAX).OR.(DTMP1.LE.TMIN)) THEN 
         DTMP = DONE / DTMP1 
         CALL DAXPBY (NLEN,VECS(1,IDX(2,1)),DTMP,VECS(1,IDX(2,1)),DZERO,&
     &VECS(1,IDX(2,1)))                                                 
         DTMP1 = DONE 
      END IF 
      IF ((DTMP2.GE.TMAX).OR.(DTMP2.LE.TMIN)) THEN 
         DTMP = DONE / DTMP2 
         CALL DAXPBY (NLEN,VECS(1,IDX(3,1)),DTMP,VECS(1,IDX(3,1)),DZERO,&
     &VECS(1,IDX(3,1)))                                                 
         DTMP2 = DONE 
      END IF 
      DWK(IDX(1,1),8*M+11) = DONE / DTMP1 
      DWK(IDX(1,1),8*M+12) = DONE / DTMP2 
!                                                                       
!     Initialize the counters.                                          
!                                                                       
      K      = 0 
      L      = 1 
      N      = 1 
      NMAX   = 0 
      KSTAR  = 1 
      LSTAR  = 1 
      MKMAX  = 0 
      NLMAX  = 1 
      PQBEG  = 1 
      VWBEG  = 1 
      IWK(IDX(1,0+1),2) = 1 
      IWK(IDX(1,1+1),2) = 1 
      IWK(IDX(1,0+1),1) = 1 
      IWK(IDX(1,1+1),1) = 1 
      MPQBLT = 1 
      MVWBLT = 1 
      MK     = IWK(IDX(1,K+1),1) 
      NL     = IWK(IDX(1,L+1),2) 
      IWK(IDX(1,N),3) = NL 
      NLSTAR = IWK(IDX(1,LSTAR+1),2) 
      IWK(IDX(1,N),5) = NLSTAR 
!                                                                       
!     Set up the QMR iteration.                                         
!                                                                       
      RESN   = DONE 
      DWK(IDX(1,1),8*M+8) = DUCPLO(1) 
      DWK(IDX(1,1),8*M+4) = DWK(IDX(1,1),8*M+8) * R0 
      MAXOMG = DONE / DWK(IDX(1,1),8*M+8) 
!                                                                       
!     This is one step of the coupled two-term Lanczos algorithm.       
!     V_n,  W_n,  P_{n-1}, Q_{n-1}, D_{n-1}, E_{n-1}, F_{n-1}, L_{n-1}, 
!     U_{n-1}, and D_{nn} are given.                                    
!                                                                       
!     Note that the previous step was  not necessarily step N-1,  as it 
!     could have been a restart.                                        
!     Except at the first step, the following hold:                     
!       NL.LE.N, MK.LE.N-1, N-NL.LE.MAXVW-1, N-MK.LE.MAXPQ.             
!                                                                       
!     The references in the comments  to step numbers correspond to the 
!     description in Section 4 of the report `Implementation details of 
!     the coupled QMR algorithm`, by Freund and Nachtigal, RIACS report 
!     92.19, October 1992.                                              
!                                                                       
   30 IWK(IDX(1,N),9)   = L 
      IWK(IDX(1,N),7)   = N 
      IWK(IDX(1,N),8)  = K 
      IWK(IDX(1,N),12)  = K 
      IWK(IDX(1,N),10) = KSTAR 
      IWK(IDX(1,N),13) = KSTAR 
      IWK(IDX(1,N),11)  = LSTAR 
write(6,*)'ducpl: right after 30'
call flush(6)
!                                                                       
!     Build p_n and q_n.                                                
!     There are three possible cases:                                   
!       N > NMAX   : the PQ sequence is being built new.                
!       N > MKMAX  : the PQ sequence is being rebuilt, build as if new. 
!       N <= MKMAX : the VW sequence is being rebuilt, just rerun PQ.   
!     The first two are identical in terms of code.                     
!     The code is inlined because of the matrix multiplications.        
!                                                                       
!     Initialize Lanczos step variables.                                
!                                                                       
      NP1        = N + 1 
      DWK(IDX(1,N),8*M+18) = DZERO 
      KBLKSZ     = N - MK 
      IBUILT     = .FALSE. 
      IERR       = 0 
!                                                                       
!     Clear current column of U.                                        
!                                                                       
      DO 40 I = 1, M 
         DWK(IDX(1,I),6*M+IDX(1,N)) = DZERO 
   40 END DO 
write(6,*)'ducpl: right after 40'
call flush(6)
      DWK(IDX(1,N),6*M+IDX(1,N)) = DONE 
!                                                                       
!     The first step is different, deal with it here.                   
!                                                                       
      IF (N.EQ.1) THEN 
         NP1     = 2 
         MKSTAR  = 1 
         DWK(IDX(1,1),8*M+16) = DONE 
         DWK(IDX(1,1),8*M+17) = DONE 
         DWK(IDX(1,1),8*M+14)  = DWK(IDX(1,1),8*M+11) 
         DWK(IDX(1,1),8*M+15)  = DWK(IDX(1,1),8*M+12) 
         INNER   = .FALSE. 
         CALL DAXPBY (NLEN,VECS(1,IDX(4,1)),DONE,VECS(1,IDX(2,1)),DZERO,&
     &VECS(1,IDX(4,1)))                                                 
         CALL DAXPBY (NLEN,VECS(1,IDX(5,1)),DONE,VECS(1,IDX(3,1)),DZERO,&
     &VECS(1,IDX(5,1)))                                                 
         GO TO 340 
      END IF 
!                                                                       
!     Step 1:                                                           
!     Update D^{(n-1)} to D^{(n)}.                                      
!                                                                       
      DO 60 I = NL, N-1 
         DTMP = DZERO 
         DO 50 J = MAX0(NLSTAR,IWK(IDX(1,I),3)), N-1 
            DTMP = DTMP + DWK(IDX(1,I),IDX(1,J)) * DWK(IDX(1,J),2*M+    &
     &IDX(1,N-1))                                                       
   50    CONTINUE 
         DTMP     = ( DWK(IDX(1,I),M+IDX(1,N-1)) - DTMP ) / DWK(IDX(1,N)&
     &,2*M+IDX(1,N-1))                                                  
         DWK(IDX(1,I),IDX(1,N)) = DTMP 
         DWK(IDX(1,N),IDX(1,I)) = DTMP * DWK(IDX(1,N),8*M+10) /         &
     &DWK(IDX(1,I),8*M+10)                                              
   60 END DO 
write(6,*)'ducpl: right after 60'
call flush(6)
!                                                                       
!     Step 2:                                                           
!     Compute k^\star.                                                  
!                                                                       
      I = KSTAR 
      DO 70 J = I+1, K 
         IF (IWK(IDX(1,J+1),1).LE.NL-1) KSTAR = J 
   70 END DO 
write(6,*)'ducpl: right after 70'
call flush(6)
      MKSTAR = IWK(IDX(1,KSTAR+1),1) 
!                                                                       
!     Check memory allocation  (MKSTAR.GE.PQBEG).  Update  the range of 
!     valid indices:  the new vectors  p_n and q_n  will overwrite  the 
!     vectors at location N-MVEC.                                       
!                                                                       
      IF (MKSTAR.LT.PQBEG) THEN 
         IERR = 64 
         GO TO 750 
      END IF 
      PQBEG = MAX0(PQBEG,N-MVEC+1) 
!                                                                       
!     Step 3:                                                           
!     Compute F_{n,1:n-1} = F_{n,m_{k^\star}:n-1}.                      
!                                                                       
      DO 90 I = MKSTAR, N-1 
         DTMP = DZERO 
         DO 80 J = MAX0(NL,IWK(IDX(1,I),5)), I+1 
            DTMP = DTMP + DWK(IDX(1,N),IDX(1,J)) * DWK(IDX(1,J),2*M+    &
     &IDX(1,I))                                                         
   80    CONTINUE 
         DWK(IDX(1,N),M+IDX(1,I)) = DTMP 
         DTMP     = DABS(DTMP) / DWK(IDX(1,I),8*M+16) 
         NORMA    = DMAX1(NORMA,DTMP) 
         NORM1    = DMAX1(ADJUST*NORMA,NORM1) 
         NORM2    = DMAX1(ADJUST*NORMA,NORM2) 
   90 END DO 
write(6,*)'ducpl: right after 90'
call flush(6)
!                                                                       
!     Step 4:                                                           
!     Check whether E_k is nonsingular.                                 
!                                                                       
      IF (KBLKSZ.EQ.1) THEN 
         INNER = DABS(DWK(IDX(1,MK),5*M+IDX(1,MK))).EQ.DZERO 
      ELSE 
         DO 110 I = MK, N-1 
            DO 100 J = MK, N-1 
               DWK(I-MK+1,4*M+J-MK+1) = DWK(IDX(1,I),5*M+IDX(1,J)) 
  100       CONTINUE 
  110    CONTINUE 
         CALL DSVDC (DWK(1,4*M+1),M,KBLKSZ,KBLKSZ,DWK(1,8*M+1),DWK(1,8*M&
     &+3),DZERO,0,DZERO,0,DWK(1,8*M+2),0,I)                             
         IF (I.NE.0) THEN 
            IERR = -I 
            GO TO 750 
         END IF 
         DTMP = DZERO 
         IF (DWK(1,8*M+1).NE.DZERO) DTMP = DWK(KBLKSZ,8*M+1) / DWK(1,8*M&
     &+1)                                                               
         INNER = DTMP.LT.(DBLE(NLEN) * DLAMCH('E')) 
      END IF 
      KEND = K 
      IF (INNER) KEND = K - 1 
write(6,*)'ducpl: here5'
call flush(6)
!                                                                       
!     Step 5:                                                           
!     Compute U_{m_i:m_{i+1}-1,n}, for i = k^\star, ..., kend, kend = k 
!     or kend = k-1.                                                    
!                                                                       
      IF (KEND.EQ.K) IWK(IDX(1,K+1+1),1) = N 
write(6,*)'ducpl: here5.1'
call flush(6)
      DO 150 I = KSTAR, KEND 
         MI     = IWK(IDX(1,I+1),1) 
         MIP1   = IWK(IDX(1,I+1+1),1) 
         IBLKSZ = MIP1 - MI 
         IF (IBLKSZ.EQ.1) THEN 
            DWK(IDX(1,MI),6*M+IDX(1,N)) = DWK(IDX(1,N),M+IDX(1,MI)) /   &
     &DWK(IDX(1,MI),5*M+IDX(1,MI)) * DWK(IDX(1,MI),8*M+10) / DWK(IDX(1,N&
     &),8*M+10)                                                         
         ELSE 
            DO 130 J = MI, MIP1-1 
               DWK(J-MI+1,8*M+1) = DWK(IDX(1,N),M+IDX(1,J)) * DWK(IDX(1,&
     &J),8*M+10) / DWK(IDX(1,N),8*M+10)                                 
               DO 120 IJ = MI, MIP1-1 
                  DWK(J-MI+1,4*M+IJ-MI+1) = DWK(IDX(1,J),5*M+IDX(1,IJ)) 
  120          CONTINUE 
  130       CONTINUE 
write(6,*)'ducpl: before dqrdc'
call flush(6)
            CALL DQRDC (DWK(1,4*M+1),M,IBLKSZ,IBLKSZ,DWK(1,8*M+3),0,    &
     &DZERO,0)                                                          
write(6,*)'ducpl: after  dqrdc'
call flush(6)
            CALL DQRSL (DWK(1,4*M+1),M,IBLKSZ,IBLKSZ,DWK(1,8*M+3),DWK(1,&
     &8*M+1),                                                           &
     &                  DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),   &
     &DZERO,100,J)                                                      
write(6,*)'ducpl: after  dqrsl'
call flush(6)
            DO 140 J = MI, MIP1-1 
               DWK(IDX(1,J),6*M+IDX(1,N)) = DWK(J-MI+1,8*M+2) 
  140       CONTINUE 
         END IF 
  150 END DO 
write(6,*)'ducpl: right after 150'
call flush(6)
!                                                                       
!     Step 6:                                                           
!     Build the part common to both inner and regular vectors.          
!                                                                       
      DWK(IDX(1,N),8*M+14) = DWK(IDX(1,N),8*M+11) 
      DWK(IDX(1,N),8*M+15) = DWK(IDX(1,N),8*M+12) 
      CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,IDX(2,N)),DZERO,VECS(1,3)) 
      CALL DAXPBY (NLEN,VECS(1,4),DONE,VECS(1,IDX(3,N)),DZERO,VECS(1,4)) 
      DO 160 I = MKSTAR, MK-1 
         DTMP = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+14) /     &
     &DWK(IDX(1,N),8*M+14)                                              
         CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,3),-DTMP,VECS(1,IDX(4,I&
     &)))                                                               
         DTMP = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+15) /     &
     &DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+10)
         CALL DAXPBY (NLEN,VECS(1,4),DONE,VECS(1,4),-DTMP,VECS(1,IDX(5,I&
     &)))                                                               
  160 END DO 
write(6,*)'ducpl: right after 160'
call flush(6)
!                                                                       
!     Check whether PQ is being rerun.                                  
!                                                                       
      RERUN = (N.LE.NMAX).AND.(N.LE.MKMAX) 
write(6,*)'ducpl: here6, rerun', rerun
call flush(6)
      IF (RERUN) THEN 
!                                                                       
!     Check whether E_k has become singular (this should not happen).   
!                                                                       
         IF ((IWK(IDX(1,K+1+1),1).EQ.N).AND.INNER) THEN 
            IERR = 128 
            GO TO 750 
         END IF 
!                                                                       
!     Determine whether this was an inner vector or not.                
!                                                                       
         INNER = IWK(IDX(1,K+1+1),1).NE.N 
         IF (VF.NE.0) THEN 
            IF (INNER) THEN 
               WRITE (VF,'(A17,I5,A14)') 'Rerunning vector ',N,         &
     &' (PQ) as inner'                                                  
            ELSE 
               WRITE (VF,'(A17,I5,A16)') 'Rerunning vector ',N,         &
     &' (PQ) as regular'                                                
            END IF 
         END IF 
!                                                                       
!     If building an inner vector, set the coefficients U_{m_k:n-1,n}.  
!                                                                       
         IF (INNER) THEN 
            DO 170 I = MK, N-1 
               DWK(IDX(1,I),6*M+IDX(1,N)) = DUCPLU(I,N) 
  170       CONTINUE 
         END IF 
      ELSE 
!                                                                       
!     PQ is not being rerun, determine whether it is being rebuilt.     
!                                                                       
         IF ((VF.NE.0).AND.(N.LE.NMAX)) THEN 
            WRITE (VF,'(A15,I8)') 'Rebuilding P&Q:', N 
         END IF 
!                                                                       
!     If E_k is singular, build an inner vector.                        
!                                                                       
         IF (INNER) THEN 
            IF (VF.NE.0) WRITE (VF,'(A31)')                             &
     &'... moment matrix E is singular'                                 
            GO TO 300 
         END IF 
      END IF 
!                                                                       
!     Either E_k is nonsingular, or PQ is being rerun.  For the latter, 
!     just finish building  the vectors.  For the former,  the check in 
!     Step 7 (G_{m_{k-1}:n-1,n-1})  could be done.  However,  the look- 
!     ahead strategy requires  the smaller of the checks  in Step 7 and 
!     Step 9, in case the norm  estimates need updating.  Hence, Step 7 
!     and Step 9  are switched,  and their  checks  are later performed 
!     together.                                                         
!                                                                       
!     Leave the common part of the vectors in the temporary vectors.    
!                                                                       
      CALL DAXPBY (NLEN,VECS(1,IDX(4,N)),DONE,VECS(1,3),DZERO,VECS(1,   &
     &IDX(4,N)))                                                        
      CALL DAXPBY (NLEN,VECS(1,IDX(5,N)),DONE,VECS(1,4),DZERO,VECS(1,   &
     &IDX(5,N)))                                                        
!                                                                       
!     Step 8:                                                           
!     Build regular vectors and compute A p_n and q_n^T A p_n.          
!                                                                       
      DO 180 I = MK, N-1 
         DTMP = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+14) /     &
     &DWK(IDX(1,N),8*M+14)                                              
         CALL DAXPBY (NLEN,VECS(1,IDX(4,N)),DONE,VECS(1,IDX(4,N)),-DTMP,&
     &VECS(1,IDX(4,I)))                                                 
         DTMP = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+15) /     &
     &DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+10)
         CALL DAXPBY (NLEN,VECS(1,IDX(5,N)),DONE,VECS(1,IDX(5,N)),-DTMP,&
     &VECS(1,IDX(5,I)))                                                 
  180 END DO 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,IDX(4,N)),VECS(1,IDX(2,NP1)))                 
!                                                                       
  190 INFO(2) = 1 
      INFO(3) = IDX(4,N) 
      INFO(4) = IDX(2,NP1) 
      RETLBL  = 200 
      RETURN 
write(6,*)'ducpl: before 200'
call flush(6)
  200 DWK(IDX(1,N),5*M+IDX(1,N)) = DDOT(NLEN,VECS(1,IDX(5,N)),1,VECS(1, &
     &IDX(2,NP1)),1) * DWK(IDX(1,N),8*M+14) * DWK(IDX(1,N),8*M+15)      
!                                                                       
!     Step 9:                                                           
!     Compute the last column of E_k^{-1}.                              
!                                                                       
      IF (KBLKSZ.EQ.1) THEN 
         DWK(IDX(1,MK),7*M+IDX(1,K)) = DONE / DWK(IDX(1,MK),5*M+IDX(1,MK&
     &))                                                                
      ELSE 
         DO 210 I = 1, KBLKSZ-1 
            DWK(I,8*M+1) = DZERO 
  210    CONTINUE 
         DWK(KBLKSZ,8*M+1) = DONE 
         CALL DQRSL (DWK(1,4*M+1),M,KBLKSZ,KBLKSZ,DWK(1,8*M+3),DWK(1,8*M&
     &+1),                                                              &
     &               DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),DZERO,&
     &100,J)                                                            
         DO 220 I = MK, N-1 
            DWK(IDX(1,I),7*M+IDX(1,K)) = DWK(I-MK+1,8*M+2) 
  220    CONTINUE 
      END IF 
!                                                                       
!     Check for zero norms.                                             
!                                                                       
      DWK(IDX(1,N),8*M+16) = DNRM2(NLEN,VECS(1,IDX(4,N)),1) * DWK(IDX(1,&
     &N),8*M+14)                                                        
      DWK(IDX(1,N),8*M+17) = DNRM2(NLEN,VECS(1,IDX(5,N)),1) * DWK(IDX(1,&
     &N),8*M+15)                                                        
      IF (DWK(IDX(1,N),8*M+16).EQ.DZERO) GO TO 360 
      IF (DWK(IDX(1,N),8*M+17).EQ.DZERO) GO TO 360 
!                                                                       
!     If PQ is being rerun, then skip to the end.                       
!                                                                       
      IF (RERUN) GO TO 360 
!                                                                       
!     Step 7:                                                           
!     Build the full column G_{m_{k-1}:n-1,n-1} and compute its norm.   
!                                                                       
      DTMP1 = DZERO 
      DTMP2 = DZERO 
      DO 240 I = IWK(IDX(1,K-1+1),1), N-1 
         DTMP = DZERO 
         DO 230 J = MAX0(I,NLSTAR), N 
            DTMP = DTMP + DWK(IDX(1,I),6*M+IDX(1,J)) * DWK(IDX(1,J),2*M+&
     &IDX(1,N-1))                                                       
  230    CONTINUE 
         DTMP  = DABS(DTMP) 
         DTMP1 = DTMP1 + DTMP * DWK(IDX(1,I),8*M+16) 
         DTMP2 = DTMP2 + DTMP * DWK(IDX(1,I),8*M+17) / DWK(IDX(1,I),8*M+&
     &10)                                                               
  240 END DO 
      DTMP1 = DTMP1 / DWK(IDX(1,N-1),8*M+16) 
      DTMP2 = DTMP2 * DWK(IDX(1,N-1),8*M+10) / DWK(IDX(1,N-1),8*M+17) 
!                                                                       
!     Step 9:                                                           
!     Build the 2nd term for the next step, G_{m_k:n-1,n}, regular.     
!     Compute the norm of G_{m_k:n-1,n}.                                
!                                                                       
      DTMP3 = DZERO 
      DTMP4 = DZERO 
      DTMP5 = DWK(IDX(1,N),5*M+IDX(1,N)) * DWK(IDX(1,N),2*M+IDX(1,N-1)) &
     &* DWK(IDX(1,N-1),8*M+10) / DWK(IDX(1,N),8*M+10)                   
      DO 250 I = MK, N-1 
         DTMP = DABS(DTMP5 * DWK(IDX(1,I),7*M+IDX(1,K))) 
         DTMP3 = DTMP3 + DTMP * DWK(IDX(1,I),8*M+16) 
         DTMP4 = DTMP4 + DTMP * DWK(IDX(1,I),8*M+17) / DWK(IDX(1,I),8*M+&
     &10)                                                               
  250 END DO 
      DTMP3 = DTMP3 / DWK(IDX(1,N),8*M+16) 
      DTMP4 = DTMP4 * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,N),8*M+17) 
!                                                                       
!     Steps 7 and 9:                                                    
!     Check G_{m_{k-1}:n-1,n-1} and G_{m_k:n-1,n}.                      
!                                                                       
      DTMP  = DMAX1(DTMP1,DTMP2,DTMP3,DTMP4) 
      INNER = DTMP.GT.NORM1 
      IF (.NOT.INNER) GO TO 360 
      DWK(IDX(1,N),8*M+18) = DTMP 
!                                                                       
!     If G_{m_{k-1}:n-1,n-1} is bad, build inner vectors.               
!                                                                       
      IF (DMAX1(DTMP1,DTMP2).GT.NORM1) GO TO 300 
!                                                                       
!     If G_{m_k:n-1,n} is bad, check the inner vectors.                 
!     This only applies if m_k > 1.                                     
!                                                                       
      IF (MK.LE.1) GO TO 300 
!                                                                       
!     Build the inner vectors to compute the 2nd term at the next step. 
!     Get the coefficients U_{m_k:n-1,n} in a temporary location.       
!     Build the inner vectors, their norms are needed.                  
!                                                                       
      IBUILT = .TRUE. 
      DO 260 I = MK, N-1 
         DWK(IDX(1,I),6*M+IDX(1,NP1)) = DUCPLU(I,N) 
         DTMP       = DWK(IDX(1,I),6*M+IDX(1,NP1)) * DWK(IDX(1,I),8*M+14&
     &) / DWK(IDX(1,N),8*M+14)                                          
         CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,3),-DTMP,VECS(1,IDX(4,I&
     &)))                                                               
         DTMP       = DWK(IDX(1,I),6*M+IDX(1,NP1)) * DWK(IDX(1,I),8*M+15&
     &) / DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M&
     &+10)                                                              
         CALL DAXPBY (NLEN,VECS(1,4),DONE,VECS(1,4),-DTMP,VECS(1,IDX(5,I&
     &)))                                                               
  260 END DO 
      DWK(IDX(1,NP1),8*M+16) = DNRM2(NLEN,VECS(1,3),1) * DWK(IDX(1,N),  &
     &8*M+14)                                                           
      DWK(IDX(1,NP1),8*M+17) = DNRM2(NLEN,VECS(1,4),1) * DWK(IDX(1,N),  &
     &8*M+15)                                                           
      IF (DWK(IDX(1,NP1),8*M+16).EQ.DZERO) GO TO 310 
      IF (DWK(IDX(1,NP1),8*M+17).EQ.DZERO) GO TO 310 
!                                                                       
!     Build the 2nd term for the next step, G_{m_{k-1}:m_k-1,n}, inner. 
!     Compute the norm of G_{m_{k-1}:m_k-1,n}.                          
!                                                                       
      DTMP = DZERO 
      DO 270 J = MKSTAR, MK-1 
         DTMP = DTMP + DWK(IDX(1,J),6*M+IDX(1,N)) * DWK(IDX(1,J),5*M+   &
     &IDX(1,MK)) / DWK(IDX(1,J),8*M+10)                                 
  270 END DO 
      DO 280 J = MK, N-1 
         DTMP = DTMP + DWK(IDX(1,J),6*M+IDX(1,NP1)) * DWK(IDX(1,J),5*M+ &
     &IDX(1,MK)) / DWK(IDX(1,J),8*M+10)                                 
  280 END DO 
      DTMP1 = DZERO 
      DTMP2 = DZERO 
      DTMP3 = DMAX1(DTMP3,DTMP4) 
      DTMP  = DWK(IDX(1,N),M+IDX(1,MK)) - DTMP * DWK(IDX(1,N),8*M+10) 
      DTMP5 = DTMP * DWK(IDX(1,MK),2*M+IDX(1,MK-1)) * DWK(IDX(1,MK-1),  &
     &8*M+10) / DWK(IDX(1,MK),8*M+10)                                   
      DO 290 I = IWK(IDX(1,K-1+1),1), MK-1 
         DTMP = DABS(DTMP5 * DWK(IDX(1,I),7*M+IDX(1,K-1))) 
         DTMP1 = DTMP1 + DTMP * DWK(IDX(1,I),8*M+16) 
         DTMP2 = DTMP2 + DTMP * DWK(IDX(1,I),8*M+17) / DWK(IDX(1,I),8*M+&
     &10)                                                               
  290 END DO 
      DTMP1 = DTMP1 / DWK(IDX(1,NP1),8*M+16) 
      DTMP2 = DTMP2 * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,NP1),8*M+17) 
!                                                                       
!     Compare the inner and regular versions of the 2nd term at the next
!     step.  Build the vector corresponding to the smaller term.        
!                                                                       
      INNER = DTMP3.GT.DMAX1(DTMP1,DTMP2) 
      IF (.NOT.INNER) GO TO 360 
!                                                                       
!     Build inner vectors.                                              
!     Check whether the P&Q block has to be forced to close.            
!                                                                       
  300 IF (VF.NE.0) WRITE (VF,'(A7,I5,A14)') 'Vector ',N,' (PQ) is inner' 
      IF (N-MK.EQ.MAXPQ) THEN 
         CALL DUCPL1 (NDIM,NLEN,M,N,K,KSTAR,L,LSTAR,MK,MKSTAR,NL,NLSTAR,&
     &                VF,IERR,ADJUST,NORM1,NORM2,DWK,IWK,IDX,VECS)      
         IF (IERR.NE.0) GO TO 750 
         INNER  = .FALSE. 
         KBLKSZ = N - MK 
         RERUN  = .TRUE. 
         NP1    = N + 1 
         GO TO 190 
      END IF 
!                                                                       
!     The temporary vectors  contain either just partial inner vectors, 
!     or the  completed ones,  depending on whether  IBUILT is  TRUE or 
!     FALSE.  In either case, replace the regular vectors.              
!                                                                       
  310 CALL DAXPBY (NLEN,VECS(1,IDX(4,N)),DONE,VECS(1,3),DZERO,VECS(1,   &
     &IDX(4,N)))                                                        
      CALL DAXPBY (NLEN,VECS(1,IDX(5,N)),DONE,VECS(1,4),DZERO,VECS(1,   &
     &IDX(5,N)))                                                        
!                                                                       
!     Step 11:                                                          
!     Get the coefficients U_{m_k:n-1,n} and build inner vectors.       
!                                                                       
      IF (IBUILT) THEN 
         DO 320 I = MK, N-1 
            DWK(IDX(1,I),6*M+IDX(1,N)) = DWK(IDX(1,I),6*M+IDX(1,NP1)) 
  320    CONTINUE 
         DWK(IDX(1,N),8*M+16) = DWK(IDX(1,NP1),8*M+16) 
         DWK(IDX(1,N),8*M+17) = DWK(IDX(1,NP1),8*M+17) 
      ELSE 
         DO 330 I = MK, N-1 
            DWK(IDX(1,I),6*M+IDX(1,N)) = DUCPLU(I,N) 
            DTMP     = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+14)&
     & / DWK(IDX(1,N),8*M+14)                                           
            CALL DAXPBY (NLEN,VECS(1,IDX(4,N)),DONE,VECS(1,IDX(4,N)),-  &
     &DTMP,VECS(1,IDX(4,I)))                                            
            DTMP     = DWK(IDX(1,I),6*M+IDX(1,N)) * DWK(IDX(1,I),8*M+15)&
     & / DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+&
     &10)                                                               
            CALL DAXPBY (NLEN,VECS(1,IDX(5,N)),DONE,VECS(1,IDX(5,N)),-  &
     &DTMP,VECS(1,IDX(5,I)))                                            
  330    CONTINUE 
         DWK(IDX(1,N),8*M+16) = DNRM2(NLEN,VECS(1,IDX(4,N)),1) *        &
     &DWK(IDX(1,N),8*M+14)                                              
         DWK(IDX(1,N),8*M+17) = DNRM2(NLEN,VECS(1,IDX(5,N)),1) *        &
     &DWK(IDX(1,N),8*M+15)                                              
      END IF 
!                                                                       
!     Step 11:                                                          
!     Compute A p_n and q_n^T A p_n.                                    
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,IDX(4,N)),VECS(1,IDX(2,NP1)))                 
!                                                                       
  340 INFO(2) = 1 
      INFO(3) = IDX(4,N) 
      INFO(4) = IDX(2,NP1) 
      RETLBL  = 350 
      RETURN 
write(6,*)'ducpl: before 350'
call flush(6)
  350 DWK(IDX(1,N),5*M+IDX(1,N)) = DDOT(NLEN,VECS(1,IDX(5,N)),1,VECS(1, &
     &IDX(2,NP1)),1) * DWK(IDX(1,N),8*M+14) * DWK(IDX(1,N),8*M+15)      
!                                                                       
!     Step 10:                                                          
!     If regular vectors were built, update the counters.               
!                                                                       
  360 IF (.NOT.INNER) THEN 
         K  = K + 1 
         MK = IWK(IDX(1,K+1),1) 
      END IF 
!                                                                       
!     Update counters.                                                  
!                                                                       
      MPQBLT = MAX0(MPQBLT,MK-IWK(IDX(1,K-1+1),1)) 
      MKMAX  = MAX0(MK,MKMAX) 
      IWK(IDX(1,N),6) = MKSTAR 
      IWK(IDX(1,N),4) = MK 
!                                                                       
!     Check for invariant subspaces.                                    
!                                                                       
      IF (DWK(IDX(1,N),8*M+16).LT.TNRM) IERR = IERR + 16 
      IF (DWK(IDX(1,N),8*M+17).LT.TNRM) IERR = IERR + 32 
      IF (IERR.NE.0) GO TO 750 
!                                                                       
!     Step 13:                                                          
!     Compute A^T q_n.                                                  
!     Have the caller carry out ATXB, then return here.                 
!        CALL ATXB (VECS(1,IDX(5,N)),VECS(1,IDX(3,NP1)))                
!                                                                       
      INFO(2) = 2 
      INFO(3) = IDX(5,N) 
      INFO(4) = IDX(3,NP1) 
      RETLBL  = 370 
      RETURN 
!                                                                       
!     Update the norm estimates.                                        
!                                                                       
write(6,*)'ducpl: before 370'
call flush(6)
  370 IF (N.LE.NUMCHK) THEN 
         DTMP1 = DNRM2(NLEN,VECS(1,IDX(2,NP1)),1) * DWK(IDX(1,N),8*M+14)&
     & / DWK(IDX(1,N),8*M+16)                                           
         DTMP2 = DNRM2(NLEN,VECS(1,IDX(3,NP1)),1) * DWK(IDX(1,N),8*M+15)&
     & / DWK(IDX(1,N),8*M+17)                                           
         NORMA = DMAX1(DTMP1,DTMP2,NORMA) 
      END IF 
write(6,*)'ducpl: right after 370'
call flush(6)
      DTMP  = DABS(DWK(IDX(1,N),5*M+IDX(1,N))) / ( DWK(IDX(1,N),8*M+16) &
     &* DWK(IDX(1,N),8*M+17) )                                          
      NORMA = DMAX1(DTMP,NORMA) 
      NORM1 = DMAX1(ADJUST*NORMA,NORM1) 
      NORM2 = DMAX1(ADJUST*NORMA,NORM2) 
!                                                                       
!     Build v_{n+1} and w_{n+1}.                                        
!     There are three cases:                                            
!       N > NMAX   : the VW sequence is being built new.                
!       N >= NLMAX : the VW sequence is being rebuilt, build as if new. 
!       N < NLMAX  : the PQ sequence is being rebuilt, just rerun VW.   
!     The first two are identical in terms of code.                     
!                                                                       
      IWK(IDX(1,N),12)   = K 
      IWK(IDX(1,N),13)  = KSTAR 
      DWK(IDX(1,N),8*M+13) = DZERO 
      IBUILT     = .FALSE. 
      LBLKSZ     = N - NL + 1 
!                                                                       
!     Clear current column of L.                                        
!                                                                       
      DO 380 I = 1, M 
         DWK(IDX(1,I),2*M+IDX(1,N)) = DZERO 
  380 END DO 
write(6,*)'ducpl: right after 380'
call flush(6)
!                                                                       
!     Step 14:                                                          
!     Update E^{(n-1)} to E^{(n)}.                                      
!                                                                       
      DO 400 I = MK, N-1 
         DTMP = DZERO 
         DO 390 J = MAX0(MKSTAR,IWK(IDX(1,I),4)), N-1 
            DTMP = DTMP + DWK(IDX(1,J),6*M+IDX(1,N)) * DWK(IDX(1,J),5*M+&
     &IDX(1,I)) / DWK(IDX(1,J),8*M+10)                                  
  390    CONTINUE 
         DTMP     = DWK(IDX(1,N),M+IDX(1,I)) - DTMP * DWK(IDX(1,N),8*M+ &
     &10)                                                               
         DWK(IDX(1,N),5*M+IDX(1,I)) = DTMP 
         DWK(IDX(1,I),5*M+IDX(1,N)) = DTMP * DWK(IDX(1,I),8*M+10) /     &
     &DWK(IDX(1,N),8*M+10)                                              
         DTMP1    = DABS(DWK(IDX(1,N),5*M+IDX(1,I))) / (DWK(IDX(1,N),8*M&
     &+17) * DWK(IDX(1,I),8*M+16))                                      
         DTMP2    = DABS(DWK(IDX(1,I),5*M+IDX(1,N))) / (DWK(IDX(1,I),8*M&
     &+17) * DWK(IDX(1,N),8*M+16))                                      
         NORMA    = DMAX1(NORMA,DTMP1,DTMP2) 
         NORM1    = DMAX1(ADJUST*NORMA,NORM1) 
         NORM2    = DMAX1(ADJUST*NORMA,NORM2) 
  400 END DO 
write(6,*)'ducpl: right after 400'
call flush(6)
!                                                                       
!     Step 15:                                                          
!     Compute l^\star.                                                  
!                                                                       
      I = LSTAR 
      DO 410 J = I+1, L 
         IF (IWK(IDX(1,J+1),2).LE.MK) LSTAR = J 
  410 END DO 
write(6,*)'ducpl: right after 410'
call flush(6)
      NLSTAR = IWK(IDX(1,LSTAR+1),2) 
!                                                                       
!     Check memory allocation  (NLSTAR.GE.VWBEG).  Update  the range of 
!     valid indices: the new vectors v_{n+1} and w_{n+1} will overwrite 
!     the vectors at location N+1-MVEC.                                 
!                                                                       
      IF (NLSTAR.LT.VWBEG) THEN 
         IERR = 64 
         GO TO 750 
      END IF 
      VWBEG = MAX0(VWBEG,NP1-MVEC+1) 
!                                                                       
!     Step 16:                                                          
!     Compute F_{1:n,n} = F_{n_{l^\star}:n,n}.                          
!                                                                       
      DO 430 I = NLSTAR, N 
         DTMP = DZERO 
         DO 420 J = MAX0(MK,IWK(IDX(1,I),6)), I 
            DTMP = DTMP + DWK(IDX(1,J),6*M+IDX(1,I)) * DWK(IDX(1,J),5*M+&
     &IDX(1,N)) / DWK(IDX(1,J),8*M+10)                                  
  420    CONTINUE 
         DWK(IDX(1,I),M+IDX(1,N)) = DTMP * DWK(IDX(1,I),8*M+10) 
         DTMP     = DABS(DTMP) / DWK(IDX(1,N),8*M+16) 
         NORMA    = DMAX1(NORMA,DTMP) 
         NORM1    = DMAX1(ADJUST*NORMA,NORM1) 
         NORM2    = DMAX1(ADJUST*NORMA,NORM2) 
  430 END DO 
write(6,*)'ducpl: right after 430'
call flush(6)
!                                                                       
!     Step 17:                                                          
!     Check whether D_l is nonsingular.                                 
!                                                                       
      IF (LBLKSZ.EQ.1) THEN 
         INNER = DABS(DWK(IDX(1,NL),IDX(1,NL))).EQ.DZERO 
      ELSE 
         DO 450 I = NL, N 
            DO 440 J = NL, N 
               DWK(I-NL+1,4*M+J-NL+1) = DWK(IDX(1,I),IDX(1,J)) 
  440       CONTINUE 
  450    CONTINUE 
         CALL DSVDC (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,8*M+1),DWK(1,8*M&
     &+3),DZERO,0,DZERO,0,DWK(1,8*M+2),0,I)                             
         IF (I.NE.0) THEN 
            IERR = -I 
            GO TO 750 
         END IF 
         DTMP = DZERO 
         IF (DWK(1,8*M+1).NE.DZERO) DTMP = DWK(LBLKSZ,8*M+1) / DWK(1,8*M&
     &+1)                                                               
         INNER = DTMP.LT.(DBLE(NLEN) * DLAMCH('E')) 
      END IF 
      LEND = L 
      IF (INNER) LEND = L - 1 
!                                                                       
!     Step 18:                                                          
!     Compute L_{n_i:n_{i+1}-1,n}, for i = l^\star, ..., lend, lend = l 
!     or lend = l-1.                                                    
!                                                                       
      IF (LEND.EQ.L) IWK(IDX(1,L+1+1),2) = NP1 
      DO 490 I = LSTAR, LEND 
         NI     = IWK(IDX(1,I+1),2) 
         NIP1   = IWK(IDX(1,I+1+1),2) 
         IBLKSZ = NIP1 - NI 
         IF (IBLKSZ.EQ.1) THEN 
            DWK(IDX(1,NI),2*M+IDX(1,N)) = DWK(IDX(1,NI),M+IDX(1,N)) /   &
     &DWK(IDX(1,NI),IDX(1,NI))                                          
         ELSE 
            DO 470 J = NI, NIP1-1 
               DWK(J-NI+1,8*M+1) = DWK(IDX(1,J),M+IDX(1,N)) 
               DO 460 IJ = NI, NIP1-1 
                  DWK(J-NI+1,4*M+IJ-NI+1) = DWK(IDX(1,J),IDX(1,IJ)) 
  460          CONTINUE 
  470       CONTINUE 
            CALL DQRDC (DWK(1,4*M+1),M,IBLKSZ,IBLKSZ,DWK(1,8*M+3),0,    &
     &DZERO,0)                                                          
            CALL DQRSL (DWK(1,4*M+1),M,IBLKSZ,IBLKSZ,DWK(1,8*M+3),DWK(1,&
     &8*M+1),                                                           &
     &                  DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),   &
     &DZERO,100,J)                                                      
            DO 480 J = NI, NIP1-1 
               DWK(IDX(1,J),2*M+IDX(1,N)) = DWK(J-NI+1,8*M+2) 
  480       CONTINUE 
         END IF 
  490 END DO 
write(6,*)'ducpl: right after 490'
call flush(6)
!                                                                       
!     Step 19:                                                          
!     Build the part common to both inner and regular vectors.          
!     Assume that V(NP1) = A p_n and W(NP1) = A^T q_n.                  
!                                                                       
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DZERO 
      DO 500 I = NLSTAR, NL-1 
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+11) /     &
     &DWK(IDX(1,N),8*M+14)                                              
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1)),- &
     &DTMP,VECS(1,IDX(2,I)))                                            
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+12) /     &
     &DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+10)
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1)),- &
     &DTMP,VECS(1,IDX(3,I)))                                            
  500 END DO 
write(6,*)'ducpl: right after 500'
call flush(6)
!                                                                       
!     Check whether VW is being rerun.                                  
!                                                                       
      RERUN = (N.LE.NMAX).AND.(N.LT.NLMAX) 
      IF (RERUN) THEN 
!                                                                       
!     Check whether D_l has become singular (this should not happen).   
!                                                                       
         IF ((IWK(IDX(1,L+1+1),2).EQ.NP1).AND.INNER) THEN 
            IERR = 128 
            GO TO 750 
         END IF 
!                                                                       
!     Determine whether this was an inner vector or not.                
!                                                                       
         INNER = IWK(IDX(1,L+1+1),2).NE.NP1 
         IF (VF.NE.0) THEN 
            IF (INNER) THEN 
               WRITE (VF,'(A17,I5,A14)') 'Rerunning vector ',NP1,       &
     &' (VW) as inner'                                                  
            ELSE 
               WRITE (VF,'(A17,I5,A16)') 'Rerunning vector ',NP1,       &
     &' (VW) as regular'                                                
            END IF 
         END IF 
!                                                                       
!     If building an inner vector, set the coefficients L_{n_l:n,n}.    
!                                                                       
         IF (INNER) THEN 
            DO 510 I = NL, N 
               DWK(IDX(1,I),2*M+IDX(1,N))  = DUCPLL(I,N) 
  510       CONTINUE 
         END IF 
      ELSE 
!                                                                       
!     VW is not being rerun, determine whether it is being rebuilt.     
!                                                                       
         IF ((VF.NE.0).AND.(N.LE.NMAX)) THEN 
            WRITE (VF,'(A15,I8)') 'Rebuilding V&W:', NP1 
         END IF 
!                                                                       
!     If E_k is singular, build an inner vector.                        
!                                                                       
         IF (INNER) THEN 
            IF (VF.NE.0) WRITE (VF,'(A31)')                             &
     &'... moment matrix D is singular'                                 
            GO TO 630 
         END IF 
      END IF 
!                                                                       
!     Either D_l is nonsingular, or VW is being rerun.  For the latter, 
!     just finish building  the vectors.  For the former,  the check in 
!     Step 20  (H_{n_{l-1}:n,n}) could be done. However, the look-ahead 
!     strategy  requires  the smaller of the checks in Step 20 and Step 
!     22, in case the norm estimates need updating.  Hence, Step 20 and 
!     Step 22  are switched,  and  their  checks  are  later  performed 
!     together.                                                         
!                                                                       
!     Save the common part of the vectors.                              
!                                                                       
      CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,IDX(2,NP1)),DZERO,VECS(1,3&
     &))                                                                
      CALL DAXPBY (NLEN,VECS(1,4),DONE,VECS(1,IDX(3,NP1)),DZERO,VECS(1,4&
     &))                                                                
!                                                                       
!     Step 21.                                                          
!     Build regular vectors.                                            
!                                                                       
      DO 520 I = NL, N 
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+11) /     &
     &DWK(IDX(1,N),8*M+14)                                              
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1)),- &
     &DTMP,VECS(1,IDX(2,I)))                                            
         DTMP = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+12) /     &
     &DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+10)
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1)),- &
     &DTMP,VECS(1,IDX(3,I)))                                            
  520 END DO 
write(6,*)'ducpl: right after 520'
call flush(6)
!                                                                       
!     Compute scale factors for the new vectors.                        
!                                                                       
  530 DTMP3      = DNRM2(NLEN,VECS(1,IDX(2,NP1)),1) 
      DTMP4      = DNRM2(NLEN,VECS(1,IDX(3,NP1)),1) 
      DTMP1      = DWK(IDX(1,N),8*M+14) * DTMP3 
      DTMP2      = DWK(IDX(1,N),8*M+15) * DTMP4 
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DTMP1 
      IF (DTMP1.LT.TNRM) IERR = IERR + 16 
      IF (DTMP2.LT.TNRM) IERR = IERR + 32 
      IF (IERR.NE.0) GO TO 670 
      DWK(IDX(1,NP1),8*M+10)     = DWK(IDX(1,N),8*M+10) * DTMP1 / DTMP2 
      DWK(IDX(1,NP1),IDX(1,NP1)) = DDOT(NLEN,VECS(1,IDX(3,NP1)),1,VECS(1&
     &,IDX(2,NP1)),1) / ( DTMP3 * DTMP4 )                               
      IF ((DTMP3.GE.TMAX).OR.(DTMP3.LE.TMIN)) THEN 
         DTMP = DONE / DTMP3 
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DTMP,VECS(1,IDX(2,NP1)),  &
     &DZERO,VECS(1,IDX(2,NP1)))                                         
         DTMP3 = DONE 
      END IF 
      IF ((DTMP4.GE.TMAX).OR.(DTMP4.LE.TMIN)) THEN 
         DTMP = DONE / DTMP4 
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DTMP,VECS(1,IDX(3,NP1)),  &
     &DZERO,VECS(1,IDX(3,NP1)))                                         
         DTMP4 = DONE 
      END IF 
      DWK(IDX(1,NP1),8*M+11) = DONE / DTMP3 
      DWK(IDX(1,NP1),8*M+12) = DONE / DTMP4 
!                                                                       
!     Step 22:                                                          
!     Compute the last column of D_l^{-1}.                              
!                                                                       
      IF (LBLKSZ.EQ.1) THEN 
         DWK(IDX(1,NL),3*M+IDX(1,L)) = DONE / DWK(IDX(1,NL),IDX(1,NL)) 
      ELSE 
         DO 540 I = 1, LBLKSZ-1 
            DWK(I,8*M+1) = DZERO 
  540    CONTINUE 
         DWK(LBLKSZ,8*M+1) = DONE 
         CALL DQRSL (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,8*M+3),DWK(1,8*M&
     &+1),                                                              &
     &               DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),DZERO,&
     &100,J)                                                            
         DO 550 I = NL, N 
            DWK(IDX(1,I),3*M+IDX(1,L)) = DWK(I-NL+1,8*M+2) 
  550    CONTINUE 
      END IF 
write(6,*)'ducpl: right after 550'
call flush(6)
!                                                                       
!     If VW is being rerun, then skip to the end.                       
!                                                                       
      IF (RERUN) GO TO 670 
!                                                                       
!     Step 20:                                                          
!     Build the full column H_{n_{l-1}:n,n} and compute its norm.       
!                                                                       
      DTMP1 = DZERO 
      DTMP2 = DZERO 
      DO 570 I = IWK(IDX(1,L-1+1),2), N 
         DTMP = DZERO 
         DO 560 J = MAX0(1,I-1), N 
            DTMP = DTMP + DWK(IDX(1,I),2*M+IDX(1,J)) * DWK(IDX(1,J),6*M+&
     &IDX(1,N))                                                         
  560    CONTINUE 
         DTMP  = DABS(DTMP) 
         DTMP1 = DTMP1 + DTMP 
         DTMP2 = DTMP2 + DTMP / DWK(IDX(1,I),8*M+10) 
  570 END DO 
      DTMP2 = DTMP2 * DWK(IDX(1,N),8*M+10) 
!                                                                       
!     Step 22:                                                          
!     Build the 2nd term for the next step, H_{n_l:n,n+1}, regular.     
!     Compute the norm of H_{n_l:n,n+1}.                                
!                                                                       
      DTMP3 = DZERO 
      DTMP4 = DZERO 
      DTMP5 = DWK(IDX(1,NP1),IDX(1,NP1)) * DWK(IDX(1,NP1),2*M+IDX(1,N)) &
     &* DWK(IDX(1,N),8*M+10) / DWK(IDX(1,NP1),8*M+10)                   
      DO 580 I = NL, N 
         DTMP  = DABS(DTMP5 * DWK(IDX(1,I),3*M+IDX(1,L))) 
         DTMP3 = DTMP3 + DTMP 
         DTMP4 = DTMP4 + DTMP / DWK(IDX(1,I),8*M+10) 
  580 END DO 
write(6,*)'ducpl: right after 580'
call flush(6)
      DTMP4 = DTMP4 * DWK(IDX(1,NP1),8*M+10) 
write(6,*)'ducpl: here 2.9'
call flush(6)
!                                                                       
!     Steps 20 and 22:                                                  
!     Check H_{n_{l-1}:n,n} and H_{n_l:n-1,n+1}.                        
!                                                                       
      DTMP  = DMAX1(DTMP1,DTMP2,DTMP3,DTMP4) 
      INNER = DTMP.GT.NORM2 
write(6,*)'ducpl: here 2.95, inner', inner
call flush(6)
      IF (.NOT.INNER) GO TO 670 
write(6,*)'ducpl: here 2.96'
call flush(6)
      DWK(IDX(1,N),8*M+13) = DTMP 
!                                                                       
!     If H_{n_{l-1}:n,n} is bad, build inner vectors.                   
!                                                                       
write(6,*)'ducpl: here3'
call flush(6)
      IF (DMAX1(DTMP1,DTMP2).GT.NORM2) GO TO 630 
write(6,*)'ducpl: here4'
call flush(6)
!                                                                       
!     If H_{n_l:n-1,n+1} is bad, check the inner vectors.               
!     This only applies if n_l > 1.                                     
!                                                                       
      IF (NL.LE.1) GO TO 630 
!                                                                       
!     Build the inner vectors to compute the 2nd term at the next step. 
!     Get the coefficients L_{n_l:n+1,n} in a temporary location.       
!     Build inner vectors in VECS(1,3) and VECS(1,4).                   
!                                                                       
      IBUILT = .TRUE. 
      DO 590 I = NL, N 
         DWK(IDX(1,I),2*M+IDX(1,NP1)) = DUCPLL(I,N) 
         DTMP       = DWK(IDX(1,I),2*M+IDX(1,NP1)) * DWK(IDX(1,I),8*M+11&
     &) / DWK(IDX(1,N),8*M+14)                                          
         CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,3),-DTMP,VECS(1,IDX(2,I&
     &)))                                                               
         DTMP       = DWK(IDX(1,I),2*M+IDX(1,NP1)) * DWK(IDX(1,I),8*M+12&
     &) / DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M&
     &+10)                                                              
         CALL DAXPBY (NLEN,VECS(1,4),DONE,VECS(1,4),-DTMP,VECS(1,IDX(3,I&
     &)))                                                               
  590 END DO 
write(6,*)'ducpl: right after 590'
call flush(6)
      IDTMP3       = DNRM2(NLEN,VECS(1,3),1) 
      IDTMP4       = DNRM2(NLEN,VECS(1,4),1) 
      IDTMP1       = DWK(IDX(1,N),8*M+14) * IDTMP3 
      IDTMP2       = DWK(IDX(1,N),8*M+15) * IDTMP4 
      DWK(IDX(1,NP1),2*M+IDX(1,NP1)) = IDTMP1 
      IF (IDTMP1.LT.TNRM) IERR = IERR + 16 
      IF (IDTMP2.LT.TNRM) IERR = IERR + 32 
      IF (IERR.NE.0) GO TO 640 
!                                                                       
!     Build the 2nd term for the next step, H_{n_{l-1}:n_l-1,n+1},      
!     inner.  Compute the norm of H_{n_{l-1}:n_l-1,n+1}.                
!                                                                       
      DTMP = DZERO 
      DO 600 J = NLSTAR, NL-1 
         DTMP = DTMP + DWK(IDX(1,NL),IDX(1,J)) * DWK(IDX(1,I),2*M+IDX(1,&
     &N))                                                               
  600 END DO 
      DO 610 J = NL, N 
         DTMP = DTMP + DWK(IDX(1,NL),IDX(1,J)) * DWK(IDX(1,I),2*M+IDX(1,&
     &NP1))                                                             
  610 END DO 
      DTMP1 = DZERO 
      DTMP2 = DZERO 
      DTMP3 = DMAX1(DTMP3,DTMP4) 
      DTMP  = ( DWK(IDX(1,NL),M+IDX(1,N)) - DTMP ) / DWK(IDX(1,NP1),2*M+&
     &IDX(1,NP1))                                                       
      DTMP5 = DTMP * DWK(IDX(1,NL),2*M+IDX(1,NL-1)) * DWK(IDX(1,NL-1),  &
     &8*M+10) / DWK(IDX(1,NL),8*M+10)                                   
      DO 620 I = IWK(IDX(1,L-1+1),2), NL-1 
         DTMP  = DABS(DWK(IDX(1,I),3*M+IDX(1,L-1)) * DTMP5) 
         DTMP1 = DTMP1 + DTMP 
         DTMP2 = DTMP2 + DTMP / DWK(IDX(1,I),8*M+10) 
  620 END DO 
write(6,*)'ducpl: right after 620'
call flush(6)
      DTMP2 = DTMP2 * DWK(IDX(1,N),8*M+10) 
!                                                                       
!     Compare the inner and regular versions of the 2nd term at the next
!     step.  Build the vector corresponding to the smaller term.        
!                                                                       
      INNER = DTMP3.GT.DMAX1(DTMP1,DTMP2) 
      IF (.NOT.INNER) GO TO 670 
!                                                                       
!     Build inner vectors.                                              
!     Check whether the V&W block has to be forced to close.            
!                                                                       
  630 IF (VF.NE.0) WRITE (VF,'(A7,I5,A14)') 'Vector ',NP1,              &
     &' (VW) is inner'                                                  
      IF (NP1-NL.EQ.MAXVW) THEN 
         CALL DUCPL2 (NDIM,NLEN,M,N,K,KSTAR,L,LSTAR,MK,MKSTAR,NL,NLSTAR,&
     &                VF,IERR,ADJUST,NORM1,NORM2,DWK,IWK,IDX,VECS)      
         IF (IERR.NE.0) GO TO 750 
         LBLKSZ = N - NL + 1 
         INNER  = .FALSE. 
         RERUN  = .TRUE. 
         NP1    = N + 1 
         GO TO 530 
      END IF 
!                                                                       
!     The temporary vectors  contain either just partial inner vectors, 
!     or the  completed ones,  depending on whether  IBUILT is  TRUE or 
!     FALSE.  In either case, replace the regular vectors.              
!                                                                       
  640 CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,3),DZERO,VECS(1, &
     &IDX(2,NP1)))                                                      
      CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,4),DZERO,VECS(1, &
     &IDX(3,NP1)))                                                      
!                                                                       
!     Step 24:                                                          
!     Get the coefficients L_{n_l:n+1,n} and build inner vectors.       
!                                                                       
      IF (IBUILT) THEN 
         DO 650 I = NL, NP1 
            DWK(IDX(1,I),2*M+IDX(1,N)) = DWK(IDX(1,I),2*M+IDX(1,NP1)) 
  650    CONTINUE 
         DTMP1 = IDTMP1 
         DTMP2 = IDTMP2 
         DTMP3 = IDTMP3 
         DTMP4 = IDTMP4 
      ELSE 
         DO 660 I = NL, N 
            DWK(IDX(1,I),2*M+IDX(1,N))  = DUCPLL(I,N) 
            DTMP      = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+11&
     &) / DWK(IDX(1,N),8*M+14)                                          
            CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1))&
     &,-DTMP,VECS(1,IDX(2,I)))                                          
            DTMP      = DWK(IDX(1,I),2*M+IDX(1,N)) * DWK(IDX(1,I),8*M+12&
     &) / DWK(IDX(1,N),8*M+15) * DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M&
     &+10)                                                              
            CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1))&
     &,-DTMP,VECS(1,IDX(3,I)))                                          
  660    CONTINUE 
         DTMP3      = DNRM2(NLEN,VECS(1,IDX(2,NP1)),1) 
         DTMP4      = DNRM2(NLEN,VECS(1,IDX(3,NP1)),1) 
         DTMP1      = DWK(IDX(1,N),8*M+14) * DTMP3 
         DTMP2      = DWK(IDX(1,N),8*M+15) * DTMP4 
         DWK(IDX(1,NP1),2*M+IDX(1,N)) = DTMP1 
         IF (DTMP1.LT.TNRM) IERR = IERR + 16 
         IF (DTMP2.LT.TNRM) IERR = IERR + 32 
      END IF 
      IF (IERR.NE.0) GO TO 670 
      DWK(IDX(1,NP1),8*M+10)     = DWK(IDX(1,N),8*M+10) * DTMP1 / DTMP2 
      DWK(IDX(1,NP1),IDX(1,NP1)) = DDOT(NLEN,VECS(1,IDX(3,NP1)),1,VECS(1&
     &,IDX(2,NP1)),1) / ( DTMP3 * DTMP4 )                               
      IF ((DTMP3.GE.TMAX).OR.(DTMP3.LE.TMIN)) THEN 
         DTMP = DONE / DTMP3 
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DTMP,VECS(1,IDX(2,NP1)),  &
     &DZERO,VECS(1,IDX(2,NP1)))                                         
         DTMP3 = DONE 
      END IF 
      IF ((DTMP4.GE.TMAX).OR.(DTMP4.LE.TMIN)) THEN 
         DTMP = DONE / DTMP4 
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DTMP,VECS(1,IDX(3,NP1)),  &
     &DZERO,VECS(1,IDX(3,NP1)))                                         
         DTMP4 = DONE 
      END IF 
      DWK(IDX(1,NP1),8*M+11) = DONE / DTMP3 
      DWK(IDX(1,NP1),8*M+12) = DONE / DTMP4 
!                                                                       
!     Step 23:                                                          
!     If regular vectors were built, update the counters.               
!                                                                       
  670 IF (.NOT.INNER) THEN 
         L  = L + 1 
         NL = IWK(IDX(1,L+1),2) 
      END IF 
write(6,*)'ducpl: right after 670'
call flush(6)
!                                                                       
!     Update counters.                                                  
!                                                                       
      MVWBLT = MAX0(MVWBLT,NL-IWK(IDX(1,L-1+1),2)) 
      NLMAX  = MAX0(NL,NLMAX) 
      IWK(IDX(1,N),5) = NLSTAR 
      IWK(IDX(1,N),3) = NL 
!                                                                       
!     Update the counter for steps taken.                               
!                                                                       
      NMAX = MAX0(NMAX,N) 
write(6,*)'ducpl: here 3.1'
call flush(6)
      IF ((N.LT.MKMAX).OR.(N.LT.NLMAX-1)) GO TO 740 
write(6,*)'ducpl: here 3.2'
call flush(6)
!                                                                       
!     The QMR code starts here.                                         
!     At this point, (N.GE.MKMAX).AND.(N.GE.NLMAX-1), so that  steps up 
!     to MIN(MKMAX,NLMAX-1) are guaranteed not to be rebuilt.  Also, no 
!     errors are allowed in IERR, other than  possibly having found one 
!     or both invariant subspaces, in which case all remaining iterates 
!     are computed.                                                     
!                                                                       
      IEND = MIN0(MKMAX,NLMAX-1) 
      IF (IERR.NE.0) IEND = N 
write(6,*)'ducpl: here 3.3, nqmr, iend, ierr', nqmr, iend, ierr
call flush(6)
  680 IF (NQMR.GT.IEND-1) GO TO 740 
write(6,*)'ducpl: here 3.4'
call flush(6)
!                                                                       
!     Update the QMR iteration counter.                                 
!                                                                       
      NQMR = NQMR + 1 
!                                                                       
!     Get the next scaling factor omega(i) and update MAXOMG.           
!                                                                       
      DWK(IDX(1,NQMR+1),8*M+8) = DUCPLO(NQMR+1) 
      MAXOMG      = DMAX1(MAXOMG,DONE/DWK(IDX(1,NQMR+1),8*M+8)) 
!                                                                       
!     Compute the starting index IBASE for the column of \hat{R}.       
!                                                                       
      IBASE       = MAX0(1,IWK(IDX(1,NQMR),5)-1) 
      DWK(IDX(1,IBASE),8*M+5) = DZERO 
!                                                                       
!     Multiply the new column by the previous omegas.                   
!                                                                       
      DO 690 J = IWK(IDX(1,NQMR),5), NQMR+1 
         DWK(IDX(1,J),8*M+5) = DWK(IDX(1,J),8*M+8) * DWK(IDX(1,J),2*M+  &
     &IDX(1,NQMR))                                                      
  690 END DO 
!                                                                       
!     Apply the previous rotations.                                     
!     The loop below explicitly implements a call to DROT:              
!     CALL DROT (1,DWK(IDX(1,J-1),8*M+5),1,DWK(IDX(1,J),8*M+5),1,DWK(IDX
!                                                                       
      DO 700 J = IBASE+1, NQMR 
         DTMP1     = DWK(IDX(1,J),8*M+5) 
         DTMP2     = DWK(IDX(1,J-1),8*M+5) 
         DWK(IDX(1,J-1),8*M+5) = DWK(IDX(1,J),8*M+9) * DTMP2 + DWK(IDX(1&
     &,J),8*M+6) * DTMP1                                                
         DWK(IDX(1,J),8*M+5)   = DWK(IDX(1,J),8*M+9) * DTMP1 - DWK(IDX(1&
     &,J),8*M+6) * DTMP2                                                
  700 END DO 
!                                                                       
!     Compute the rotation for the last element (this also applies it). 
!                                                                       
      CALL DROTG (DWK(IDX(1,NQMR),8*M+5),DWK(IDX(1,NQMR+1),8*M+5),      &
     &DWK(IDX(1,NQMR+1),8*M+9),DWK(IDX(1,NQMR+1),8*M+6))                
!                                                                       
!     Apply the new rotation to the right-hand side vector.             
!     Could be replaced with:                                           
!        DWK(IDX(1,NQMR+1),8*M+4) = DZERO                               
!        CALL DROT (1,DWK(IDX(1,NQMR),8*M+4),1,DWK(IDX(1,NQMR+1),8*M+4),
!                                                                       
      DWK(IDX(1,NQMR+1),8*M+4) = -DWK(IDX(1,NQMR+1),8*M+6) * DWK(IDX(1, &
     &NQMR),8*M+4)                                                      
      DWK(IDX(1,NQMR),8*M+4)   =  DWK(IDX(1,NQMR+1),8*M+9) * DWK(IDX(1, &
     &NQMR),8*M+4)                                                      
!                                                                       
!     Compute the next search direction s_i.                            
!     This is more complicated than it might have to be because storage 
!     for the vectors VECS(1,IDX(6,NQMR)) is minimized.                 
!                                                                       
      DTMP2 = DZERO 
      DTMP  = DWK(IDX(1,NQMR),8*M+14) 
      DO 710 J = IBASE, NQMR-1 
         DTMP1 = DWK(IDX(1,J),8*M+5) * DWK(IDX(1,J),8*M+7) / DTMP 
         IF (DTMP1.EQ.DZERO) GO TO 710 
         CALL DAXPBY (NLEN,VECS(1,IDX(6,NQMR)),DTMP2,VECS(1,IDX(6,NQMR))&
     &,-DTMP1,VECS(1,IDX(6,J)))                                         
         DTMP2 = DONE 
  710 END DO 
      CALL DAXPBY (NLEN,VECS(1,IDX(6,NQMR)),DTMP2,VECS(1,IDX(6,NQMR)),  &
     &DONE,VECS(1,IDX(4,NQMR)))                                         
      DTMP = DTMP / DWK(IDX(1,NQMR),8*M+5) 
!                                                                       
!     Compute the new QMR iterate, then scale the search direction.     
!                                                                       
      DTMP1 = DTMP * DWK(IDX(1,NQMR),8*M+4) 
      CALL DAXPBY (NLEN,VECS(1,1),DONE,VECS(1,1),DTMP1,VECS(1,IDX(6,NQMR&
     &)))                                                               
      DWK(IDX(1,NQMR),8*M+7) = DTMP 
      DTMP      = DABS(DTMP) 
      IF ((DTMP.GE.TMAX).OR.(DTMP.LE.TMIN)) THEN 
         DWK(IDX(1,NQMR),8*M+7) = DONE 
         CALL DAXPBY (NLEN,VECS(1,IDX(6,NQMR)),DTMP,VECS(1,IDX(6,NQMR)),&
     &DZERO,VECS(1,IDX(6,NQMR)))                                        
      END IF 
!                                                                       
!     Compute the residual norm upper bound.                            
!     If the scaled upper bound is within one order of magnitude of the 
!     target convergence norm, compute the true residual norm.          
!                                                                       
      UNRM = DSQRT(DBLE(NQMR+1)) * MAXOMG * DABS(DWK(IDX(1,NQMR+1),8*M+4&
     &)) / R0                                                           
      UCHK = UNRM 
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 730 
!                                                                       
!     Have the caller carry out AXB, then return here.                  
!        CALL AXB (VECS(1,1),VECS(1,3))                                 
!                                                                       
      INFO(2) = 1 
      INFO(3) = 1 
      INFO(4) = 3 
      RETLBL  = 720 
      RETURN 
write(6,*)'ducpl: before 720'
call flush(6)
  720 CALL DAXPBY (NLEN,VECS(1,3),DONE,VECS(1,2),-DONE,VECS(1,3)) 
      RESN = DNRM2(NLEN,VECS(1,3),1) / R0 
      UCHK = RESN 
!                                                                       
!     Output the convergence history.                                   
!                                                                       
  730 IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') NQMR, UNRM, RESN 
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') NQMR, UNRM, RESN 
!                                                                       
!     Check for convergence or termination.  Stop if:                   
!         1. algorithm converged;                                       
!         2. there is an error condition;                               
!         3. the residual norm upper bound is smaller than the computed 
!     residual norm by a factor of at least 100;                        
!         4. algorithm exceeded the iterations limit.                   
!                                                                       
      IF (RESN.LE.TOL) THEN 
         IERR = 0 
         GO TO 750 
      ELSE IF (IERR.NE.0) THEN 
         GO TO 750 
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN 
         IERR = 4 
         GO TO 750 
      ELSE IF (NQMR.GE.NLIM) THEN 
         IERR = 4 
         GO TO 750 
      END IF 
      GO TO 680 
!                                                                       
!     Update the running counter.                                       
!                                                                       
  740 IF (IERR.NE.0) GO TO 750 
write(6,*)'ducpl: after 740, ierr, n, nlim', ierr, n, nlim
call flush(6)
      IERR = 4 
      IF (N.GE.NLIM) GO TO 750 
      N = N + 1 
write(6,*)'ducpl: here4'
call flush(6)
      GO TO 30 
!                                                                       
!     Done.                                                             
!                                                                       
  750 RETLBL   = 0 
      NLIM     = NQMR 
      INFO(1)  = IERR 
      MAXPQ    = MPQBLT 
      MAXVW    = MVWBLT 
      NORMS(1) = NORM1 
      NORMS(2) = NORM2 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!********************************************************************** 
!                                                                       
      SUBROUTINE DUCPL1 (NDIM,NLEN,M,N,K,KSTAR,L,LSTAR,MK,MKSTAR,NL,    &
     &                   NLSTAR,VF,IERR,ADJUST,NORM1,NORM2,DWK,IWK,     &
     &                   IDX,VECS)                                      
!                                                                       
!     Purpose:                                                          
!     This subroutine rebuilds the data  for vectors p_n and q_n at the 
!     point where a block was  forced to close.  The routine  is called 
!     internally by the coupled Lanczos code,  and is  not meant  to be 
!     called directly by the user code.                                 
!                                                                       
!     Parameters:                                                       
!        See descriptions in the main routine DUCPL.                    
!                                                                       
!     External routines used:                                           
!     subroutine daxpby(n,z,a,x,b,y)                                    
!        Library routine, computes z = a * x + b * y.                   
!     double precision ddot(n,x,incx,y,incy)                            
!        BLAS-1 routine, computes y^H * x.                              
!     double precision dnrm2(n,x,incx)                                  
!        BLAS-1 routine, computes the 2-norm of x.                      
!     subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)                   
!        LINPACK routine, computes the QR factorization of x.           
!     subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)      
!        LINPACK routine, applies the QR factorization of x.            
!                                                                       
!     Noel M. Nachtigal                                                 
!     May 31, 1993                                                      
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC DMAX1, MAX0 
      EXTERNAL DAXPBY, DQRDC, DQRSL 
!                                                                       
      INTEGER IERR, K, KSTAR, L, LSTAR, M, MK, MKSTAR, NL, NLSTAR, N 
!      INTEGER IDX(6,*), IWK(M,13), NDIM, NLEN, VF 
      INTEGER IDX(6,*), IWK(M,13), VF 
      integer*8 ndim, nlen
      DOUBLE PRECISION ADJUST, DWK(M,8*M+18), NORM1, NORM2 
      DOUBLE PRECISION VECS(NDIM,*) 
!                                                                       
!     Miscellaneous parameters.                                         
!                                                                       
      DOUBLE PRECISION DONE, DZERO 
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0) 
!                                                                       
!     Local variables.                                                  
!                                                                       
      INTEGER I, IJ, J, KBLKSZ, NP1 
      DOUBLE PRECISION DTMP1, DTMP2 
!                                                                       
!     Find the index of the vector pair with the smallest pass value.   
!                                                                       
      IERR = 0 
      IF (VF.NE.0) WRITE (VF,'(A23)') 'PQ block did not close:' 
      J     = MK + 1 
      DTMP1 = DWK(IDX(1,J),8*M+18) 
      DO 100 I = MK+2, N 
         DTMP2 = DWK(IDX(1,I),8*M+18) 
         IF (DTMP2.GT.DZERO) THEN 
            IF ((DTMP1.EQ.DZERO).OR.(DTMP2.LT.DTMP1)) THEN 
               J     = I 
               DTMP1 = DTMP2 
            END IF 
         END IF 
  100 END DO 
      IF (DTMP1.EQ.DZERO) THEN 
         IF (VF.NE.0) WRITE (VF,'(A47)')                                &
     &'... no new norm estimates available (aborting).'                 
         IERR = 8 
         RETURN 
      END IF 
      NORM1 = ADJUST * DTMP1 
      NORM2 = DMAX1(NORM1,NORM2) 
      IF (VF.NE.0) WRITE (VF,'(A40,I5,2E11.4)')                         &
     &   '... updated norms, restarting from step:', IWK(IDX(1,J),7),   &
     &NORM1, NORM2                                                      
      IF (IWK(IDX(1,J),7).EQ.N) RETURN 
      N      = IWK(IDX(1,J),7) 
      L      = IWK(IDX(1,N),9) 
      NL     = IWK(IDX(1,L+1),2) 
      K      = IWK(IDX(1,N),8) 
      MK     = IWK(IDX(1,K+1),1) 
      LSTAR  = IWK(IDX(1,N),11) 
      NLSTAR = IWK(IDX(1,LSTAR+1),2) 
      KSTAR  = IWK(IDX(1,N),10) 
      MKSTAR = IWK(IDX(1,KSTAR+1),1) 
!                                                                       
!     Initialize local variables.                                       
!                                                                       
      DWK(IDX(1,N),6*M+IDX(1,N))   = DONE 
      NP1        = N + 1 
      DWK(IDX(1,N),8*M+18) = DZERO 
      KBLKSZ     = N - MK 
!                                                                       
!     Step 2:                                                           
!     Compute k^\star.                                                  
!                                                                       
      I = KSTAR 
      DO 110 J = I+1, K 
         IF (IWK(IDX(1,J+1),1).LE.NL-1) KSTAR = J 
  110 END DO 
      MKSTAR = IWK(IDX(1,KSTAR+1),1) 
!                                                                       
!     Compute U_{m_k:n-1,n}.  Save the old coefficients.                
!                                                                       
      IWK(IDX(1,K+1+1),1) = N 
      IF (KBLKSZ.EQ.1) THEN 
         DWK(IDX(1,MK),6*M+IDX(1,NP1)) = DWK(IDX(1,MK),6*M+IDX(1,N)) 
         DWK(IDX(1,MK),6*M+IDX(1,N))   = DWK(IDX(1,N),M+IDX(1,MK)) /    &
     &DWK(IDX(1,MK),5*M+IDX(1,MK)) * DWK(IDX(1,MK),8*M+10) / DWK(IDX(1,N&
     &),8*M+10)                                                         
      ELSE 
         DO 130 J = MK, N-1 
            DWK(J-MK+1,8*M+1) = DWK(IDX(1,N),M+IDX(1,J)) * DWK(IDX(1,J),&
     &8*M+10) / DWK(IDX(1,N),8*M+10)                                    
            DO 120 IJ = MK, N-1 
               DWK(J-MK+1,4*M+IJ-MK+1) = DWK(IDX(1,J),5*M+IDX(1,IJ)) 
  120       CONTINUE 
  130    CONTINUE 
         CALL DQRDC (DWK(1,4*M+1),M,KBLKSZ,KBLKSZ,DWK(1,8*M+3),0,DZERO,0&
     &)                                                                 
         CALL DQRSL (DWK(1,4*M+1),M,KBLKSZ,KBLKSZ,DWK(1,8*M+3),DWK(1,8*M&
     &+1),                                                              &
     &               DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),DZERO,&
     &100,J)                                                            
         DO 140 J = MK, N-1 
            DWK(IDX(1,J),6*M+IDX(1,NP1)) = DWK(IDX(1,J),6*M+IDX(1,N)) 
            DWK(IDX(1,J),6*M+IDX(1,N))   = DWK(J-MK+1,8*M+2) 
  140    CONTINUE 
      END IF 
!                                                                       
!     Convert inner vectors to regular vectors.                         
!                                                                       
      DO 150 I = MK, N-1 
         DTMP1 = DWK(IDX(1,I),6*M+IDX(1,N)) - DWK(IDX(1,I),6*M+IDX(1,NP1&
     &))                                                                
         DTMP2 = DTMP1 * DWK(IDX(1,I),8*M+14) / DWK(IDX(1,N),8*M+14) 
         CALL DAXPBY (NLEN,VECS(1,IDX(4,N)),DONE,VECS(1,IDX(4,N)),-DTMP2&
     &,VECS(1,IDX(4,I)))                                                
         DTMP2 = DTMP1 * DWK(IDX(1,I),8*M+15) / DWK(IDX(1,N),8*M+15) *  &
     &DWK(IDX(1,N),8*M+10) / DWK(IDX(1,I),8*M+10)                       
         CALL DAXPBY (NLEN,VECS(1,IDX(5,N)),DONE,VECS(1,IDX(5,N)),-DTMP2&
     &,VECS(1,IDX(5,I)))                                                
  150 END DO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!********************************************************************** 
!                                                                       
      SUBROUTINE DUCPL2 (NDIM,NLEN,M,N,K,KSTAR,L,LSTAR,MK,MKSTAR,NL,    &
     &                   NLSTAR,VF,IERR,ADJUST,NORM1,NORM2,DWK,IWK,     &
     &                   IDX,VECS)                                      
!                                                                       
!     Purpose:                                                          
!     This subroutine rebuilds the data for vectors v_{n+1} and w_{n+1} 
!     at the point where a block  was forced  to close.  The routine is 
!     called internally  by the coupled Lanczos code,  and is not meant 
!     to be called directly by the user code.                           
!                                                                       
!     Parameters:                                                       
!        See descriptions in the main routine DUCPL.                    
!                                                                       
!     External routines used:                                           
!     subroutine daxpby(n,z,a,x,b,y)                                    
!        Library routine, computes z = a * x + b * y.                   
!     subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)                   
!        LINPACK routine, computes the QR factorization of x.           
!     subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)      
!        LINPACK routine, applies the QR factorization of x.            
!                                                                       
!     Noel M. Nachtigal                                                 
!     May 31, 1993                                                      
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC DABS, DMAX1, MAX0 
      EXTERNAL DAXPBY, DQRDC, DQRSL 
!                                                                       
      INTEGER IERR, K, KSTAR, L, LSTAR, M, MK, MKSTAR, NL, NLSTAR, N 
!      INTEGER IDX(6,*), IWK(M,13), NDIM, NLEN, VF 
      INTEGER IDX(6,*), IWK(M,13), VF 
      integer*8 ndim, nlen
      DOUBLE PRECISION ADJUST, DWK(M,8*M+18), NORM1, NORM2 
      DOUBLE PRECISION VECS(NDIM,*) 
!                                                                       
!     Miscellaneous parameters.                                         
!                                                                       
      DOUBLE PRECISION DONE, DZERO 
      PARAMETER (DONE = 1.0D0,DZERO = 0.0D0) 
!                                                                       
!     Local variables.                                                  
!                                                                       
      INTEGER I, IJ, J, LBLKSZ, NP1 
      DOUBLE PRECISION DTMP1, DTMP2, SCALV, SCALW 
!                                                                       
!     Find the index of the vector pair with the smallest pass value.   
!                                                                       
      IERR = 0 
      IF (VF.NE.0) WRITE (VF,'(A23)') 'VW block did not close:' 
      J     = NL 
      DTMP1 = DWK(IDX(1,J),8*M+13) 
      DO 100 I = NL+1, N 
         DTMP2 = DWK(IDX(1,I),8*M+13) 
         IF (DTMP2.GT.DZERO) THEN 
            IF ((DTMP1.EQ.DZERO).OR.(DTMP2.LT.DTMP1)) THEN 
               J     = I 
               DTMP1 = DTMP2 
            END IF 
         END IF 
  100 END DO 
      IF (DTMP1.EQ.DZERO) THEN 
         IF (VF.NE.0) WRITE (VF,'(A47)')                                &
     &'... no new norm estimates available (aborting).'                 
         IERR = 8 
         RETURN 
      END IF 
      NORM2 = ADJUST * DTMP1 
      NORM1 = DMAX1(NORM1,NORM2) 
      IF (VF.NE.0) WRITE (VF,'(A40,I5,2E11.4)')                         &
     &      '... updated norms, restarting from step:', IWK(IDX(1,J),7),&
     & NORM1, NORM2                                                     
      IF (IWK(IDX(1,J),7).EQ.N) RETURN 
      N      = IWK(IDX(1,J),7) 
      L      = IWK(IDX(1,N),9) 
      NL     = IWK(IDX(1,L+1),2) 
      K      = IWK(IDX(1,N),12) 
      MK     = IWK(IDX(1,K+1),1) 
      LSTAR  = IWK(IDX(1,N),11) 
      NLSTAR = IWK(IDX(1,LSTAR+1),2) 
      KSTAR  = IWK(IDX(1,N),13) 
      MKSTAR = IWK(IDX(1,KSTAR+1),1) 
!                                                                       
!     Initialize local variables.                                       
!                                                                       
      NP1        = N + 1 
      DWK(IDX(1,N),8*M+13) = DZERO 
      LBLKSZ     = N - NL + 1 
!                                                                       
!     Step 15:                                                          
!     Compute l^\star.                                                  
!                                                                       
      I = LSTAR 
      DO 110 J = I+1, L 
         IF (IWK(IDX(1,J+1),2).LE.MK) LSTAR = J 
  110 END DO 
      NLSTAR = IWK(IDX(1,LSTAR+1),2) 
!                                                                       
!     Compute L_{n_l:n,n}.  Save the old coefficients.                  
!                                                                       
      IWK(IDX(1,L+1+1),2) = NP1 
      IF (LBLKSZ.EQ.1) THEN 
         DWK(IDX(1,NL),2*M+IDX(1,NP1)) = DWK(IDX(1,NL),2*M+IDX(1,N)) 
         DWK(IDX(1,NL),2*M+IDX(1,N))   = DWK(IDX(1,NL),M+IDX(1,N)) /    &
     &DWK(IDX(1,NL),IDX(1,NL))                                          
      ELSE 
         DO 130 J = NL, N 
            DWK(J-NL+1,8*M+1) = DWK(IDX(1,J),M+IDX(1,N)) 
            DO 120 IJ = NL, N 
               DWK(J-NL+1,4*M+IJ-NL+1) = DWK(IDX(1,J),IDX(1,IJ)) 
  120       CONTINUE 
  130    CONTINUE 
         CALL DQRDC (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,8*M+3),0,DZERO,0&
     &)                                                                 
         CALL DQRSL (DWK(1,4*M+1),M,LBLKSZ,LBLKSZ,DWK(1,8*M+3),DWK(1,8*M&
     &+1),                                                              &
     &               DZERO,DWK(1,8*M+1),DWK(1,8*M+2),DWK(1,8*M+1),DZERO,&
     &100,J)                                                            
         DO 140 J = NL, N 
            DWK(IDX(1,J),2*M+IDX(1,NP1)) = DWK(IDX(1,J),2*M+IDX(1,N)) 
            DWK(IDX(1,J),2*M+IDX(1,N))   = DWK(J-NL+1,8*M+2) 
  140    CONTINUE 
      END IF 
!                                                                       
!     Convert inner vectors to regular vectors.                         
!                                                                       
      SCALV = DWK(IDX(1,NP1),2*M+IDX(1,N)) * DWK(IDX(1,NP1),8*M+11) 
      SCALW = DWK(IDX(1,NP1),2*M+IDX(1,N)) * DWK(IDX(1,NP1),8*M+12) *   &
     &DWK(IDX(1,N),8*M+10) / DWK(IDX(1,NP1),8*M+10)                     
      DWK(IDX(1,NP1),2*M+IDX(1,N)) = DZERO 
      DO 150 I = NL, N 
         DTMP1 = DWK(IDX(1,I),2*M+IDX(1,N)) - DWK(IDX(1,I),2*M+IDX(1,NP1&
     &))                                                                
         DTMP2 = DTMP1 * DWK(IDX(1,I),8*M+11) / SCALV 
         CALL DAXPBY (NLEN,VECS(1,IDX(2,NP1)),DONE,VECS(1,IDX(2,NP1)),- &
     &DTMP2,VECS(1,IDX(2,I)))                                           
         DTMP2 = DTMP1 * DWK(IDX(1,I),8*M+12) / SCALW * DWK(IDX(1,N),8*M&
     &+10) / DWK(IDX(1,I),8*M+10)                                       
         CALL DAXPBY (NLEN,VECS(1,IDX(3,NP1)),DONE,VECS(1,IDX(3,NP1)),- &
     &DTMP2,VECS(1,IDX(3,I)))                                           
  150 END DO 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!********************************************************************** 
!                                                                       
      DOUBLE PRECISION FUNCTION DUCPLL(I,J) 
!                                                                       
!     Purpose:                                                          
!     Return the recurrence coefficients for inner vectors VW.          
!                                                                       
!     Parameters:                                                       
!     I = row index of the coefficient (input).                         
!     J = column index of the coefficient (input).                      
!                                                                       
!     Noel M. Nachtigal                                                 
!     April 17, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
!                                                                       
!     Common block DUCPLX.                                              
!                                                                       
      DOUBLE PRECISION NORMA 
      COMMON /DUCPLX/NORMA 
!                                                                       
!                                                                       
      INTEGER I, J 
!                                                                       
      IF ((I.LT.1).OR.(J.LT.1)) THEN 
         DUCPLL = 0.0D0 
      ELSE IF (I.EQ.J) THEN 
         DUCPLL = NORMA 
      ELSE 
         DUCPLL = 0.0D0 
      END IF 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!********************************************************************** 
!                                                                       
      DOUBLE PRECISION FUNCTION DUCPLO (I) 
!                                                                       
!     Purpose:                                                          
!     Return the scaling parameter OMEGA(I).                            
!                                                                       
!     Parameters:                                                       
!     I = the index of the parameter OMEGA (input).                     
!                                                                       
!     Noel M. Nachtigal                                                 
!     June 1, 1992                                                      
!                                                                       
!********************************************************************** 
!                                                                       
      INTEGER I 
!                                                                       
      DUCPLO = 1.0D0 
!                                                                       
      RETURN 
      END                                           
!                                                                       
!********************************************************************** 
!                                                                       
      DOUBLE PRECISION FUNCTION DUCPLU(I,J) 
!                                                                       
!     Purpose:                                                          
!     Return the recurrence coefficients for inner vectors PQ.          
!                                                                       
!     Parameters:                                                       
!     I = row index of the coefficient (input).                         
!     J = column index of the coefficient (input).                      
!                                                                       
!     Noel M. Nachtigal                                                 
!     April 17, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
      INTEGER I, J 
!                                                                       
      IF ((I.LT.1).OR.(J.LT.1)) THEN 
         DUCPLU = 0.0D0 
      ELSE 
         DUCPLU = 0.0D0 
      END IF 
!                                                                       
      RETURN 
      END                                           
