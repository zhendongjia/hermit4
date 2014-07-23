      SUBROUTINE GAS_POTENTIAL_ADD_INNER_POT(XI,XIDOT,FIRR,FD)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIRR(3),FD(3)
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
      DENS = DENS0*EXP(-YEAR/T_DEP)*R12**(-KG)
      DENS_DOT = DENS*(-1/T_DEP)/TWOPI + DENS*(-KG)*R12_DOT/R12
C
      FIN = DENS*TWOPI/2*R12**(KG-2)/(2-KG)
      FIN_DOT = FIN*DENS_DOT/DENS + FIN*(KG-2)*R12_DOT/R12
C
      DELT_R = 1.0
      FOUT = 0.0
      FOUT_DOT = 0.0
C
      IF(R12.GE.R_IN.AND.R12.LT.R_EDGE) THEN
         ROUT = R_EDGE + DELT_R
         ROUT_DOT = 0.0
         CALL GET_OUTER_GRAVITY(R12,R12_DOT,ROUT,ROUT_DOT,
     &                           DENS,DENS_DOT,FOUT,FOUT_DOT)
         FIN_ABS = FIN*(1-R_IN**(2-KG)) + FOUT
         FIN_ABS_DOT = FIN_DOT*(1-RIN**(2-KG)) + FOUT_DOT
      END IF
C
      IF(R12.LT.R_IN) THEN
         ROUT = R_EDGE + DELT_R
         ROUT_DOT = 0.0
         CALL GET_OUTER_GRAVITY(R12,R12_DOT,ROUT,ROUT_DOT,
     &                           DENS,DENS_DOT,FOUT,FOUT_DOT)
         FIN_ABS = FIN*(1-R12**(2-KG)) + FOUT
         FIN_ABS_DOT = FIN_DOT*(1-R12**(2-KG)) 
     &                 + FIN*R12**(1-KG)*(KG-2)*R12_DOT + FOUT_DOT
      END IF
C
      IF(R12.GE.R_EDGE) THEN
         ROUT = R12 + DELT_R
         ROUT_DOT = R12_DOT
         CALL GET_OUTER_GRAVITY(R12,R12_DOT,ROUT,ROUT_DOT,
     &                           DENS,DENS_DOT,FOUT,FOUT_DOT)
         FIN_ABS = FIN*(1-R_IN**(2-KG)+R_EDGE**(2-KG)-R12**(2-KG))
     &                  + FOUT
         FIN_ABS_DOT = FIN_DOT*(1-R_IN**(2-KG)+R_EDGE**(2-KG)
     &                                           -R12**(2-KG))
     &             + FIN*(KG-2)*R12**(1-KG)*R12_DOT 
     &             + FOUT_DOT
      END IF
C      
      DO 3 K = 1, 3
         FIRR(K) = FIRR(K) + FIN_ABS*THE(K)
         FD(K) = FD(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 3    CONTINUE
C
      END
