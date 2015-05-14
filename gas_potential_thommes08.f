      SUBROUTINE GAS_POTENTIAL_THOMMES08(XI,XIDOT,FIRR,FD)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIRR(3),FD(3)
      REAL*8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), 
     &      DENS, DENS_DOT, VC2, VC2_DOT, 
     &      FIN_ABS, FIN_ABS_DOT
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
      DENS = DENS0*EXP(-YEAR/T_DEP)
      DENS_DOT = DENS*(-1/T_DEP)
C
      VC2 = TWOPI*DENS
      VC2_DOT = VC2*DENS_DOT/DENS  
C
      FIN_ABS = VC2/R12
      FIN_ABS_DOT = -VC2*R12_DOT/R12**2 + VC2_DOT/R12
C
      DO 3 K = 1, 3
         FIRR(K) = FIRR(K) + FIN_ABS*THE(K)
         FD(K) = FD(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 3    CONTINUE
C
      END
