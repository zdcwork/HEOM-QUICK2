!                                                                       
! modified on 18 Aug, 2011:                                             
! change data type from 'integer' to 'integer * 8' for some variables   
!                                                                       
!*********************************************************************  
!                                                                       
!     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal        
!     All rights reserved.                                              
!                                                                       
!     This code is part of a copyrighted package.  For details, see the 
!     file `cpyrit.doc' in the top-level directory.                     
!                                                                       
!     ***************************************************************** 
!     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE 
!                             COPYRIGHT NOTICE                          
!     ***************************************************************** 
!                                                                       
!********************************************************************** 
!                                                                       
      SUBROUTINE ZAXPBY (N,ZZ,ZA,ZX,ZB,ZY) 
!                                                                       
!     Purpose:                                                          
!     This subroutine computes ZZ = ZA * ZX + ZB * ZY.  Several special 
!     cases are handled separately:                                     
!        ZA =  0.0, ZB =  0.0 => ZZ = 0.0                               
!        ZA =  0.0, ZB =  1.0 => ZZ = ZY  (this is COPY)                
!        ZA =  0.0, ZB = -1.0 => ZZ = -ZY                               
!        ZA =  0.0, ZB =   ZB => ZZ = ZB * ZY  (this is SCAL)           
!        ZA =  1.0, ZB =  0.0 => ZZ = ZX  (this is COPY)                
!        ZA =  1.0, ZB =  1.0 => ZZ = ZX + ZY                           
!        ZA =  1.0, ZB = -1.0 => ZZ = ZX - ZY                           
!        ZA =  1.0, ZB =   ZB => ZZ = ZX + ZB * ZY (this is AXPY)       
!        ZA = -1.0, ZB =  0.0 => ZZ = -ZX                               
!        ZA = -1.0, ZB =  1.0 => ZZ = -ZX + ZY                          
!        ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY                          
!        ZA = -1.0, ZB =   ZB => ZZ = -ZX + ZB * ZY                     
!        ZA =   ZA, ZB =  0.0 => ZZ = ZA * ZX  (this is SCAL)           
!        ZA =   ZA, ZB =  1.0 => ZZ = ZA * ZX + ZY  (this is AXPY)      
!        ZA =   ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY                      
!        ZA =   ZA, ZB =   ZB => ZZ = ZA * ZX + ZB * ZY                 
!     ZZ may be the same as ZX or ZY.                                   
!                                                                       
!     Parameters:                                                       
!     N  = the dimension of the vectors (input).                        
!     ZZ = the vector result (output).                                  
!     ZA = scalar multiplier for ZX (input).                            
!     ZX = one of the vectors (input).                                  
!     ZB = scalar multiplier for ZY (input).                            
!     ZY = the other vector (input).                                    
!                                                                       
!     Noel M. Nachtigal                                                 
!     March 23, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC DIMAG, DREAL 
!                                                                       
      INTEGER*8 N 
      DOUBLE COMPLEX ZA, ZB, ZX(N), ZY(N), ZZ(N) 
!                                                                       
!     Local variables.                                                  
!                                                                       
      INTEGER*8 I 
      DOUBLE PRECISION DAI, DAR, DBI, DBR 
!                                                                       
      IF (N.LE.0) RETURN 
!                                                                       
      DAI = DIMAG(ZA) 
      DAR = DREAL(ZA) 
      DBI = DIMAG(ZB) 
      DBR = DREAL(ZB) 
      IF ((DAR.EQ.0.0D0).AND.(DAI.EQ.0.0D0)) THEN 
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 0.0, ZB = 0.0 => ZZ = 0.0.                             
            DO 10 I = 1, N 
               ZZ(I) = (0.0D0,0.0D0) 
   10       CONTINUE 
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 0.0, ZB = 1.0 => ZZ = ZY (this is COPY).               
            DO 20 I = 1, N 
               ZZ(I) = ZY(I) 
   20       CONTINUE 
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 0.0, ZB = -1.0 => ZZ = -ZY.                            
            DO 30 I = 1, N 
               ZZ(I) = -ZY(I) 
   30       CONTINUE 
         ELSE 
!           ZA = 0.0, ZB = ZB => ZZ = ZB * ZY (this is SCAL).           
            DO 40 I = 1, N 
               ZZ(I) = ZB * ZY(I) 
   40       CONTINUE 
         END IF 
      ELSE IF ((DAR.EQ.1.0D0).AND.(DAI.EQ.0.0D0)) THEN 
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 1.0, ZB = 0.0 => ZZ = ZX (this is COPY).               
            DO 50 I = 1, N 
               ZZ(I) = ZX(I) 
   50       CONTINUE 
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 1.0, ZB = 1.0 => ZZ = ZX + ZY.                         
            DO 60 I = 1, N 
               ZZ(I) = ZX(I) + ZY(I) 
   60       CONTINUE 
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = 1.0, ZB = -1.0 => ZZ = ZX - ZY.                        
            DO 70 I = 1, N 
               ZZ(I) = ZX(I) - ZY(I) 
   70       CONTINUE 
         ELSE 
!           ZA = 1.0, ZB = ZB => ZZ = ZX + ZB * ZY (this is AXPY).      
            DO 80 I = 1, N 
               ZZ(I) = ZX(I) + ZB * ZY(I) 
   80       CONTINUE 
         END IF 
      ELSE IF ((DAR.EQ.-1.0D0).AND.(DAI.EQ.0.0D0)) THEN 
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = -1.0, ZB = 0.0 => ZZ = -ZX                             
            DO 90 I = 1, N 
               ZZ(I) = -ZX(I) 
   90       CONTINUE 
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = -1.0, ZB = 1.0 => ZZ = -ZX + ZY                        
            DO 100 I = 1, N 
               ZZ(I) = -ZX(I) + ZY(I) 
  100       CONTINUE 
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY.                      
            DO 110 I = 1, N 
               ZZ(I) = -ZX(I) - ZY(I) 
  110       CONTINUE 
         ELSE 
!           ZA = -1.0, ZB = ZB => ZZ = -ZX + ZB * ZY                    
            DO 120 I = 1, N 
               ZZ(I) = -ZX(I) + ZB * ZY(I) 
  120       CONTINUE 
         END IF 
      ELSE 
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = ZA, ZB = 0.0 => ZZ = ZA * ZX (this is SCAL).           
            DO 130 I = 1, N 
               ZZ(I) = ZA * ZX(I) 
  130       CONTINUE 
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = ZA, ZB = 1.0 => ZZ = ZA * ZX + ZY (this is AXPY)       
            DO 140 I = 1, N 
               ZZ(I) = ZA * ZX(I) + ZY(I) 
  140       CONTINUE 
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN 
!           ZA = ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY.                    
            DO 150 I = 1, N 
               ZZ(I) = ZA * ZX(I) - ZY(I) 
  150       CONTINUE 
         ELSE 
!           ZA = ZA, ZB = ZB => ZZ = ZA * ZX + ZB * ZY.                 
            DO 160 I = 1, N 
               ZZ(I) = ZA * ZX(I) + ZB * ZY(I) 
  160       CONTINUE 
         END IF 
      END IF 
!                                                                       
      RETURN 
      END                                           
