      SUBROUTINE GAS_POTENTIAL(XI,XIDOT,FIN,FDN)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIN(3),FDN(3)
      REAL*8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), 
     &      DENS, DENS_DOT, FIN_ABS, FIN_ABS_DOT
C
      R12 = 0.0
      R12_DOT = 0.0
      DO 1 K = 1, 3 
      R12 = R12 + XI(K)**2
      R12_DOT = R12_DOT + XI(K)*XIDOT(K)
 1    CONTINUE
      R12 = SQRT(R12)
      R12_DOT = R12_DOT/R12      
C
      DO 2 K =1, 3
         THE(K) = XI(K)/R12
         THE_DOT(K) = XIDOT(K)/R12 - XI(K)*R12_DOT/R12**2
 2    CONTINUE
C
      YEAR = TIME/TWOPI
      DENS = DENS0*EXP(-YEAR/T_DEP)*R12**(-1.5)
      DENS_DOT = DENS0*(-1/T_DEP)*EXP(-YEAR/T_DEP)*R12**(-1.5) 
     &           + DENS0*(-1.5)*R12**(-2.5)*R12_DOT*EXP(-YEAR/T_DEP)
C
      FIN_ABS = TWOPI*DENS*(0.2*(R12/R_EDGE)**(2.5)
     &          + 32*(R12/R_EDGE)**(4.5)) 
C
      FIN_ABS_DOT = TWOPI*DENS_DOT*(0.2*(R12/R_EDGE)**(2.5)
     &                              + 32*(R12/R_EDGE)**(4.5))
     &              + TWOPI*DENS*(0.5*(R12/R_EDGE)**(1.5)
     &                              + 144*(R12/R_EDGE)**(3.5))
     &              *(R12_DOT/R_EDGE)
C      
      DO 3 K = 1, 3
         FIN(K) = FIN(K) + FIN_ABS*THE(K)
         FDN(K) = FDN(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 3    CONTINUE
C
      END
