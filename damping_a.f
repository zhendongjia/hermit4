
      SUBROUTINE DAMPING_A(XI,XIDOT,FIRR,FD,I) 
      INCLUDE 'commonp.h'
      INTEGER I
      REAL*8 XI(3), XIDOT(3), FIRR(3), FD(3)
      REAL*8 YEAR, R12, R12_DOT, V12, V12_DOT, THE(3), THE_DOT(3),
     &       VCI(3), VCI_DOT(3), V_KEP, V_KEP_DOT,
     &       DELT_V(3), DELT_V_DOT(3), ABS_DELT_V, ABS_DELT_V_DOT,
     &       HI, DENS_S, DENS_G, HI_DOT, DENS_S_DOT, DENS_G_DOT,
     &       T_TIDAL1_DOT, T_TIDAL2_DOT,
     &       SEMI, SEMI_V, SEMI_V_DOT, SEMI_DOT 
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
       SEMI_V = 2.0/R12 - V12**2
       SEMI_V_DOT = 2.0*(-1)/R12**2*R12_DOT
     &         - 2*V12*V12_DOT  
       SEMI = 1./SEMI_V
       SEMI_DOT = (-1)*SEMI*SEMI_V_DOT/SEMI_V
C     
      HI = 0.05*SEMI**1.25*0.5
      HI_DOT = HI*1.25*SEMI_DOT/SEMI
C
      V_KEP = SEMI**(-0.5) - 1.625*HI**2*SEMI**(-2.5)
      V_KEP_DOT = -0.5*SEMI**(-1.5)*SEMI_DOT 
     &     - 1.625*2*HI*HI_DOT*SEMI**(-2.5)
     &     - 1.625*HI**2*(-2.5)*SEMI**(-3.5)*SEMI_DOT
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
      DENS_S = DENS0*EXP(-YEAR/T_DEP)*SEMI**(-KG)
      DENS_S_DOT = DENS_S*(-1/T_DEP)/TWOPI + DENS_S*(-KG)*SEMI_DOT/SEMI 
      DENS_G = DENS_S/HI
      DENS_G_DOT = DENS_G*DENS_S_DOT/DENS_S - DENS_G*HI_DOT/HI
C
C
      T_TIDAL1(I) = (DENS_P/DENS_G)*RADIUS(I)/(ABS_DELT_V*CD)
      T_TIDAL1_DOT = T_TIDAL1(I)*(-DENS_G_DOT)/DENS_G
     &         + T_TIDAL1(I)*(-ABS_DELT_V_DOT)/ABS_DELT_V
C
      T_TIDAL2(I) = (1/BODY(I))*(1/(DENS_S*SEMI**2))*(HI/SEMI)**4
     &                   *SEMI**1.5*0.2
      T_TIDAL2_DOT = T_TIDAL2(I)*(-DENS_S_DOT)/DENS_S
     &           + T_TIDAL2(I)*(-2*SEMI_DOT)/SEMI
     &           + T_TIDAL2(I)*(4*HI_DOT)/HI
     &           + T_TIDAL2(I)*(-4*SEMI_DOT)/SEMI
     &           + T_TIDAL2(I)*1.5*SEMI_DOT/SEMI    
C
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


