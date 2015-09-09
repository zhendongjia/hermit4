
      SUBROUTINE DAMPING(XI,XIDOT,FIRR,FD,I) 
      INCLUDE 'commonp.h'
      INTEGER I
      REAL*8 XI(3), XIDOT(3), FIRR(3), FD(3)
      REAL*8 YEAR, R12, R12_DOT, V12, V12_DOT, THE(3), THE_DOT(3),
     &       VCI(3), VCI_DOT(3), V_KEP, V_KEP_DOT,
     &       DELT_V(3), DELT_V_DOT(3), ABS_DELT_V, ABS_DELT_V_DOT,
     &       HI, DENS_S, DENS_G, HI_DOT, DENS_S_DOT, DENS_G_DOT,
     &       T_TIDAL1_DOT, T_TIDAL2_DOT
C
C
      R12 = 0.0
      R12_DOT = 0.0
      V12 = 0.0
      V12_DOT = 0.0
      DO 1 K = 1, 3 
      R12 = R12 + XI(K)**2
      R12_DOT = R12_DOT + XI(K)*XIDOT(K)
      V12 = V12 + XIDOT(K)**2
      V12_DOT = V12_DOT + XIDOT(K)*FIRR(K)
 1    CONTINUE
      R12 = SQRT(R12)
      R12_DOT = R12_DOT/R12
      V12 = SQRT(V12)
      V12_DOT = V12_DOT/V12
C
      DO 2 K =1, 3
         THE(K) = XI(K)/R12
         THE_DOT(K) = XIDOT(K)/R12 - XI(K)*R12_DOT/R12**2
 2     CONTINUE
C
      HI = 0.05*R12**1.25 * 0.5
      HI_DOT = HI*1.25*R12_DOT/R12
C
      V_KEP = R12**(-0.5) - 1.625*HI**2*R12**(-2.5)
      V_KEP_DOT = -0.5*R12**(-1.5)*R12_DOT 
     &     - 1.625*2*HI*HI_DOT*R12**(-2.5)
     &     - 1.625*HI**2*(-2.5)*R12**(-3.5)*R12_DOT
C      
      VCI(1) = - V_KEP*THE(2) 
      VCI(2) = V_KEP*THE(1)
      VCI(3) = V_KEP*THE(3)
C
      VCI_DOT(1) = - V_KEP_DOT*THE(2) - V_KEP*THE_DOT(2)
      VCI_DOT(2) = V_KEP_DOT*THE(1) + V_KEP*THE_DOT(1)
      VCI_DOT(3) = V_KEP_DOT*THE(3) + V_KEP*THE_DOT(3)
C
      DO 3 K = 1, 3
         DELT_V(K) = XIDOT(K) - VCI(K)
         DELT_V_DOT(K) = FIRR(K) - VCI_DOT(K)
 3    CONTINUE
C
      ABS_DELT_V = 0.0
      ABS_DELT_V_DOT = 0.0
      DO 4 K = 1, 3
         ABS_DELT_V = ABS_DELT_V + DELT_V(K)**2
         ABS_DELT_V_DOT = ABS_DELT_V_DOT + DELT_V(K)*DELT_V_DOT(K)
 4    CONTINUE
      ABS_DELT_V = SQRT(ABS_DELT_V)
      ABS_DELT_V_DOT = ABS_DELT_V_DOT/ABS_DELT_V
C
      YEAR = TIME/TWOPI
      DENS_S = DENS0*EXP(-YEAR/T_DEP)*R12**(-KG)
      DENS_S_DOT = DENS_S*(-1/T_DEP)/TWOPI + DENS_S*(-KG)*R12_DOT/R12 
      DENS_G = DENS_S/HI
      DENS_G_DOT = DENS_G*DENS_S_DOT/DENS_S - DENS_G*HI_DOT/HI
C
C
      T_TIDAL1(I) = (DENS_P/DENS_G)*RADIUS(I)/(ABS_DELT_V*CD)
      T_TIDAL1_DOT = T_TIDAL1(I)*(-DENS_G_DOT)/DENS_G
     &         + T_TIDAL1(I)*(-ABS_DELT_V_DOT)/ABS_DELT_V
C
      T_TIDAL2(I) = (1/BODY(I))*(1/(DENS_S*R12**2))*(HI/R12)**4
     &                   *R12**1.5*0.2
      T_TIDAL2_DOT = T_TIDAL2(I)*(-DENS_S_DOT)/DENS_S
     &           + T_TIDAL2(I)*(-2*R12_DOT)/R12
     &           + T_TIDAL2(I)*(4*HI_DOT)/HI
     &           + T_TIDAL2(I)*(-4*R12_DOT)/R12
     &           + T_TIDAL2(I)*R12_DOT/R12
     &           + T_TIDAL2(I)*(-V_KEP_DOT)/V_KEP    
C
C      T_UNIT = SQRT(1./2.)
C      M_UNIT = 2.0
C      V_UNIT = SQRT(2.)
C      T_RATIO = T_TIDAL1(I)/(TWOPI*T_DEP)
C      write(0,*) DENS_P/DENS_G, RADIUS(I), 1./(ABS_DELT_V*CD), T_TIDAL1(I)/(TWOPI*T_DEP)
C     write(0,*) 1/BODY(I),1/(DENS_S*R12**2),(HI/R12)**4, R12**1.5/V_UNIT, T_TIDAL2(I)/(TWOPI*T_DEP)
C      IF (T_RATIO.LT.1) THEN
C         WRITE (0,*) ABS_DELT_V*V_UNIT, R12, T_RATIO
C      END IF
C
      DO 5 K = 1, 3
         FIRR(K) = FIRR(K) + (-DELT_V(K))*(1/T_TIDAL1(I)+1/T_TIDAL2(I))
         FD(K) = FD(K) 
     &        + (-DELT_V_DOT(K))*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
     &        + (-DELT_V(K))*(-T_TIDAL1_DOT/T_TIDAL1(I)**2
     &           -T_TIDAL2_DOT/T_TIDAL2(I)**2)
 5    CONTINUE
C
      END


