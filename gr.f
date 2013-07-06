      SUBROUTINE GR(XI,XIDOT,FIN,FDN,I)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIN(3),FDN(3)
      REAL*8 SEMII, ECCI, R12, R12_DOT, V12, V12_DOT, V_C,
     &          SEMII_DOT, ECCI_DOT
C
      R12 = 0.0
      V12 = 0.0
      R12_DOT = 0.0
      V12_DOT = 0.0
      DO 1 K =1, 3
         R12 = R12 + XI(K)**2
         V12 = V12 + XIDOT(K)**2
         R12_DOT = R12_DOT + XI(K)*XIDOT(K)
         V12_DOT = V12_DOT + XIDOT(K)*FIN(K)
 1       CONTINUE
       R12 = SQRT(R12)
       R12_DOT = R12_DOT/R12
       V12 = SQRT(V12)
       V12_DOT = V12_DOT/V12
C
       SEMII = 2/R12 - V12**2/(1+BODY(I))
       SEMII = 1/SEMII
       SEMII_DOT = (-SEMII**2)*(-2*R12_DOT/R12**2 
     &                - 2*V12*V12_DOT/(1+ BODY(I)))
       ECCI = (1.0 - R12/SEMII)**2 +(R12_DOT*R12)**2/(SEMII*(1+BODY(I)))
       ECCI = SQRT(ECCI)
       ECCI_DOT = 0
C       ECCI_DOT = (0.5/ECCI)*(2*(1-R12/SEMII)*(-R12_DOT/SEMII+R12  
C     &                + R12*2*SEMII_DOT/SEMII**3 
C     &                + 2*R12_DOT*R12_DOT_DOT*R12**2/)
C     
       V_C = 173*365/TWOPI
       DO  2 K = 1, 3
          FIN(K) = FIN(K) - 3*SEMII*(1-ECCI**2)*XI(K)/(V_C**2*R12**5)
          FDN(K) = FDN(K) -3*SEMII_DOT*(1-ECCI**2)*XI(K)/(V_C**2*R12**5)
     &         - 3*SEMII*(1-2*ECCI*ECCI_DOT)*XI(K)/(V_C**2*R12**5)
     &         - 3*SEMII*(1-ECCI**2)*XIDOT(K)/(V_C**2*R12**5)
     &         +15*SEMII*(1-ECCI**2)*XI(K)*R12_DOT/(V_C**2*R12**6)
 2        CONTINUE
C
       END

