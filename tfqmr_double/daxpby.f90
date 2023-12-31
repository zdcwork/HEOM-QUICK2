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
      SUBROUTINE DAXPBY (N,DZ,DA,DX,DB,DY) 
!                                                                       
!     Purpose:                                                          
!     This subroutine computes DZ = DA * DX + DB * DY.  Several special 
!     cases are handled separately:                                     
!        DA =  0.0, DB =  0.0 => DZ = 0.0                               
!        DA =  0.0, DB =  1.0 => DZ = DY  (this is COPY)                
!        DA =  0.0, DB = -1.0 => DZ = -DY                               
!        DA =  0.0, DB =   DB => DZ = DB * DY  (this is SCAL)           
!        DA =  1.0, DB =  0.0 => DZ = DX  (this is COPY)                
!        DA =  1.0, DB =  1.0 => DZ = DX + DY                           
!        DA =  1.0, DB = -1.0 => DZ = DX - DY                           
!        DA =  1.0, DB =   DB => DZ = DX + DB * DY (this is AXPY)       
!        DA = -1.0, DB =  0.0 => DZ = -DX                               
!        DA = -1.0, DB =  1.0 => DZ = -DX + DY                          
!        DA = -1.0, DB = -1.0 => DZ = -DX - DY                          
!        DA = -1.0, DB =   DB => DZ = -DX + DB * DY                     
!        DA =   DA, DB =  0.0 => DZ = DA * DX  (this is SCAL)           
!        DA =   DA, DB =  1.0 => DZ = DA * DX + DY  (this is AXPY)      
!        DA =   DA, DB = -1.0 => DZ = DA * DX - DY                      
!        DA =   DA, DB =   DB => DZ = DA * DX + DB * DY                 
!     DZ may be the same as DX or DY.                                   
!                                                                       
!     Parameters:                                                       
!     N  = the dimension of the vectors (input).                        
!     DZ = the vector result (output).                                  
!     DA = scalar multiplier for DX (input).                            
!     DX = one of the vectors (input).                                  
!     DB = scalar multiplier for DY (input).                            
!     DY = the other vector (input).                                    
!                                                                       
!     Noel M. Nachtigal                                                 
!     March 23, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
!      INTEGER N 
      integer*8 n
      DOUBLE PRECISION DA, DB, DX(N), DY(N), DZ(N) 
!                                                                       
!     Local variables.                                                  
!                                                                       
!      INTEGER I 
      integer*8 i
!                                                                       
      IF (N.LE.0) RETURN 
!                                                                       
      IF (DA.EQ.0.0D0) THEN 
         IF (DB.EQ.0.0D0) THEN 
!           DA = 0.0, DB = 0.0 => DZ = 0.0.                             
            DO 10 I = 1, N 
               DZ(I) = 0.0D0 
   10       CONTINUE 
         ELSE IF (DB.EQ.1.0D0) THEN 
!           DA = 0.0, DB = 1.0 => DZ = DY (this is COPY).               
            DO 20 I = 1, N 
               DZ(I) = DY(I) 
   20       CONTINUE 
         ELSE IF (DB.EQ.-1.0D0) THEN 
!           DA = 0.0, DB = -1.0 => DZ = -DY.                            
            DO 30 I = 1, N 
               DZ(I) = -DY(I) 
   30       CONTINUE 
         ELSE 
!           DA = 0.0, DB = DB => DZ = DB * DY (this is SCAL).           
            DO 40 I = 1, N 
               DZ(I) = DB * DY(I) 
   40       CONTINUE 
         END IF 
      ELSE IF (DA.EQ.1.0D0) THEN 
         IF (DB.EQ.0.0D0) THEN 
!           DA = 1.0, DB = 0.0 => DZ = DX (this is COPY).               
            DO 50 I = 1, N 
               DZ(I) = DX(I) 
   50       CONTINUE 
         ELSE IF (DB.EQ.1.0D0) THEN 
!           DA = 1.0, DB = 1.0 => DZ = DX + DY.                         
            DO 60 I = 1, N 
               DZ(I) = DX(I) + DY(I) 
   60       CONTINUE 
         ELSE IF (DB.EQ.-1.0D0) THEN 
!           DA = 1.0, DB = -1.0 => DZ = DX - DY.                        
            DO 70 I = 1, N 
               DZ(I) = DX(I) - DY(I) 
   70       CONTINUE 
         ELSE 
!           DA = 1.0, DB = DB => DZ = DX + DB * DY (this is AXPY).      
            DO 80 I = 1, N 
               DZ(I) = DX(I) + DB * DY(I) 
   80       CONTINUE 
         END IF 
      ELSE IF (DA.EQ.-1.0D0) THEN 
         IF (DB.EQ.0.0D0) THEN 
!           DA = -1.0, DB = 0.0 => DZ = -DX                             
            DO 90 I = 1, N 
               DZ(I) = -DX(I) 
   90       CONTINUE 
         ELSE IF (DB.EQ.1.0D0) THEN 
!           DA = -1.0, DB = 1.0 => DZ = -DX + DY                        
            DO 100 I = 1, N 
               DZ(I) = -DX(I) + DY(I) 
  100       CONTINUE 
         ELSE IF (DB.EQ.-1.0D0) THEN 
!           DA = -1.0, DB = -1.0 => DZ = -DX - DY.                      
            DO 110 I = 1, N 
               DZ(I) = -DX(I) - DY(I) 
  110       CONTINUE 
         ELSE 
!           DA = -1.0, DB = DB => DZ = -DX + DB * DY                    
            DO 120 I = 1, N 
               DZ(I) = -DX(I) + DB * DY(I) 
  120       CONTINUE 
         END IF 
      ELSE 
         IF (DB.EQ.0.0D0) THEN 
!           DA = DA, DB = 0.0 => DZ = DA * DX (this is SCAL).           
            DO 130 I = 1, N 
               DZ(I) = DA * DX(I) 
  130       CONTINUE 
         ELSE IF (DB.EQ.1.0D0) THEN 
!           DA = DA, DB = 1.0 => DZ = DA * DX + DY (this is AXPY)       
            DO 140 I = 1, N 
               DZ(I) = DA * DX(I) + DY(I) 
  140       CONTINUE 
         ELSE IF (DB.EQ.-1.0D0) THEN 
!           DA = DA, DB = -1.0 => DZ = DA * DX - DY.                    
            DO 150 I = 1, N 
               DZ(I) = DA * DX(I) - DY(I) 
  150       CONTINUE 
         ELSE 
!           DA = DA, DB = DB => DZ = DA * DX + DB * DY.                 
            DO 160 I = 1, N 
               DZ(I) = DA * DX(I) + DB * DY(I) 
  160       CONTINUE 
         END IF 
      END IF 
!                                                                       
      RETURN 
      END                                           
