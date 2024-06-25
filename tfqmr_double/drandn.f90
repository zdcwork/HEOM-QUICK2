!********************************************************************** 
!                                                                       
      SUBROUTINE DRANDN (N,DX,SEED) 
!                                                                       
!     Purpose:                                                          
!     Fills the vector DX with random numbers  between 0 and 1.  If the 
!     SEED is given, it should be odd and positive.  The generator is a 
!     fairly unsophisticated one, from Pearson's  "Numerical methods in 
!     engineering and science" book.                                    
!                                                                       
!     Parameters:                                                       
!     N    = the dimension of the vector (input).                       
!     DX   = the vector to fill with random numbers (output).           
!     SEED = the seed for the generator (input).                        
!                                                                       
!     Noel M. Nachtigal                                                 
!     April 23, 1993                                                    
!                                                                       
!********************************************************************** 
!                                                                       
      INTRINSIC DBLE, IABS, MOD 
!                                                                       
!      INTEGER N, SEED 
      INTEGER SEED 
      integer*8 n
      DOUBLE PRECISION DX(N) 
!                                                                       
!     Local variables.                                                  
!                                                                       
!      INTEGER I, J 
      integer*8 i, j
!                                                                       
!     Local variables that are saved from one call to the next.         
!                                                                       
      DOUBLE PRECISION DMAX 
      INTEGER IM, IMAX, IS 
      SAVE DMAX, IM, IMAX, IS 
      DATA IM/0/ 
!                                                                       
!     Initialize the generator data.                                    
!                                                                       
      IF (IM.EQ.0) THEN 
         J  = 0 
         IM = 1 
         DO 10 I = 1, 31 
            J = J + 1 
            IF (IM*2.LE.IM) GO TO 20 
            IM = IM * 2 
   10    CONTINUE 
   20    IMAX = (IM-1) * 2 + 1 
         DMAX = DBLE(IMAX) 
         DO 30 I = 1, MOD(J,3) 
            J = J - 1 
            IM = IM / 2 
   30    CONTINUE 
         IM = IM + 5 
         IS = IABS(MOD(IM*30107,IMAX)) 
      END IF 
!                                                                       
!     Check whether we have a new seed.                                 
!                                                                       
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1 
!                                                                       
!     Here goes the rest.                                               
!                                                                       
      DO 40 I = 1, N 
         DX(I) = DBLE(IS) / DMAX 
         IS    = IABS(MOD(IM*IS,IMAX)) 
   40 END DO 
!                                                                       
      RETURN 
      END                                           
