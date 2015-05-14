      SUBROUTINE GAS_POTENTIAL_S(XI,XIDOT,FIRR,FD)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIRR(3),FD(3)
      REAL*8 YEAR, R12, R12_DOT, THE(3), THE_DOT(3), 
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

      IF (R12 .LE. AS_LOWBER_BOUND) THEN
         FR = AS_F(1)
         FR_DOT = AS_FDOT(1)
      ELSE
         IDX = (R12 - AS_LOWER_BOUND) / AS_STEP + 1
         IF (IDX .GE. AS_LINE) THEN
            FR = AS_F(AS_LINE)
            FR_DOT = AS_FDOT(AS_LINE)
         ELSE
            F1 = R12 - AS_LOWER_BOUND - AS_STEP * (IDX - 1)
            F0 = AS_STEP - F1
            FR = (AS_F(IDX) * F0 + AS_F(IDX+1) * F1) / (F0 + F1)
            FR_DOT = (AS_FDOT(IDX) * F0 + AS_FDOT(IDX+1) * F1)
     &           / (F0 + F1)
         END IF
      END IF
C
      YEAR = TIME/TWOPI
      FIN_ABS = FR * EXP(-YEAR/T_DEP)
      FIN_ABS_DOT = FR_DOT * (R12_DOT/R12) * EXP(-YEAR/T_DEP)
     &     + FR * EXP(-YEAR/T_DEP) * (-1.0/T_DEP)/TWOPI
      FIN_ABS = FIN_ABS * DENS0
      FIN_ABS_DOT = FIN_ABS_DOT * DENS0
C     
      DO 3 K = 1, 3
         FIRR(K) = FIRR(K) + FIN_ABS*THE(K)
         FD(K) = FD(K) + FIN_ABS_DOT*THE(K) + FIN_ABS*THE_DOT(K)
 3    CONTINUE
C
      END
