      SUBROUTINE GAS_POTENTIAL_RED(XI,XIDOT,FIRR,FD)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIRR(3),FD(3)
      REAL*8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), 
     &     DENS, DENS_DOT, FIN_ABS, FIN_ABS_DOT,
     &     DENS_IN, DENS_IND_DOT
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
      DENS = DENS0*EXP(-YEAR/T_DEP)*R12**(-KG)
      DENS_DOT = DENS0*(-1/T_DEP)*EXP(-YEAR/T_DEP)*R12**(-KG)/TWOPI 
     &     + DENS0*(-KG)*R12**(-KG-1)*R12_DOT*EXP(-YEAR/T_DEP)
      DENS_IN = DENS * 0.2
      DENS_IN_DOT = DENS_DOT * 0.2
C
      FIN_ABS = TWOPI*( DENS*((0.5/(1+KG))*(R12/R_EDGE)**(1+KG)
     &     + (0.5625/(3+KG))*(R12/R_EDGE)**(3+KG))
     &   - DENS_IN*((0.75/(4-KG))*(R_IN/R12)**(4-KG)
     &           + (0.703125/(6-KG))*(R_IN/R12)**(6-KG)
     &           + (1.0/(2-KG))*(R_IN/R12)**(2-KG)) )
C
      FIN_ABS_DOT = TWOPI*( DENS_DOT *(0.5/(1+KG)*(R12/R_EDGE)**(1+KG)
     &     + 0.5625/(3+KG)*(R12/R_EDGE)**(3+KG))
     &   - DENS_IN_DOT *(0.75/(4-KG)*(R_IN/R12)**(4-KG)
     &            + 0.703125/(6-KG)*(R_IN/R12)**(6-KG)
     &            + 1.0/(2-KG)*(R_IN/R12)**(2-KG)) )
     & + TWOPI*DENS*(R12_DOT/R_EDGE)*(0.5*(R12/R_EDGE)**(KG)
     &     + 0.5625*(R12/R_EDGE)**(2+KG))
     & + TWOPI*DENS_IN*(R12_DOT/R_IN)*(0.75*(R_IN/R12)**(5-KG)
     &     + 0.703125*(R_IN/R12)**(7-KG)
     &     + 1.0/(R_IN/R12)**(3-KG))             
C
      DO 3 K = 1, 3
         FIRR(K) = FIRR(K) + FIN_ABS*THE(K)
         FD(K) = FD(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 3    CONTINUE
C
      END
