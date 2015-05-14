      SUBROUTINE NBINT(I,IKS)
*
*
*       N-body integration with iteration.
*       ----------------------------------
*
      INCLUDE 'commonp.h'
      PARAMETER  (ONE24=1.0/24.0D0)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3),DV(3),FIN(3),FDN(3),
     &        X0T(3),V0T(3),AT3(3),BT2(3), GCM2_MAU2, GCM3_MAU3
      SAVE IPLOT
      DATA IPLOT /0/
      DATA ECRIT /1.00D0/
*
      IF (TWOGG.GT.0.AND.TIME.LE.TWOPI*T_S) THEN
         BODY(3) = (M_S - M_S_INT)*TIME/(TWOPI*T_S) + M_S_INT
         R_EDGE = (R_S - R_EDGE_INT)*TIME/(TWOPI*T_S) + R_EDGE_INT
      END IF
*          
      IESC = 0
*       Check regularization criterion for single particles.
      IF (STEP(I).LT.DTMIN.AND.I.LE.N.AND.IFIRST.EQ.1) THEN
*       See whether dominant body can be regularized.
          IF (IKS.EQ.0) THEN
              CALL SEARCH(I,IKS)
          END IF
      END IF
*
*       Form time-step factors and update T0.
      DT = TIME - T0(I)
      DTSQ = DT**2
      DT6 = 6.0/(DT*DTSQ)
      DT2 = 2.0/DTSQ
      DT12 = ONE12*DT
      DT60 = 0.2*DT12
*     DTSQ12 = ONE12*DTSQ
*     DT13 = ONE3*DT
      DT3 = DT*DTSQ
      DT4 = ONE24*DT3*DT
      DT3 = ONE6*DT3
      DT025 = 0.25*DT
      DT02 = 0.2*DT
      T0(I) = TIME
*
*       Predict body #I to order F3DOT and initialize scalars.
      DO 5 K = 1,3
*         X0T(K) = X(K,I)
*         V0T(K) = XDOT(K,I)
          XI(K) = X(K,I) + (D3(K,I)*DT02 + D2(K,I))*DT4
          XIDOT(K) = XDOT(K,I) + (D3(K,I)*DT025 + D2(K,I))*DT3
          FIN(K) = 0.0D0
          FDN(K) = 0.0D0
    5 CONTINUE
**
*       Treat c.m. body more carefully
      IF (I.GT.N) THEN
          CALL CMF(XI,XIDOT,FIN,FDN)
          J = I
*       Add indirect force due to body #I.
          GO TO 32
      END IF
*
*
*       Obtain force & derivative from non-zero mass particles.
      DO 10 J = IFIRST,NMASS
          IF (J.EQ.I) GO TO 8
          IF (BODY(J).LT.M_CRIT) GO TO 10
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIN(1) = FIN(1) + A1*DR3I
          FIN(2) = FIN(2) + A2*DR3I
          FIN(3) = FIN(3) + A3*DR3I
          FDN(1) = FDN(1) + (DV(1) - A1*DRDV)*DR3I
          FDN(2) = FDN(2) + (DV(2) - A2*DRDV)*DR3I
          FDN(3) = FDN(3) + (DV(3) - A3*DRDV)*DR3I
    8     IF (KZ(3).EQ.0) GO TO 10
*       Include indirect terms with option #3.
          XJ = X(1,J)
          YJ = X(2,J)
          ZJ = X(3,J)
          RJ2 = XJ**2 + YJ**2 + ZJ**2
          FJ = BODY(J)/(RJ2*SQRT (RJ2))
          FIN(1) = FIN(1) - XJ*FJ
          FIN(2) = FIN(2) - YJ*FJ
          FIN(3) = FIN(3) - ZJ*FJ
          RRD = 3.0*(XJ*XDOT(1,J) + YJ*XDOT(2,J) + ZJ*XDOT(3,J))/RJ2
          FDN(1) = FDN(1) - (XDOT(1,J) - RRD*XJ)*FJ
          FDN(2) = FDN(2) - (XDOT(2,J) - RRD*YJ)*FJ
          FDN(3) = FDN(3) - (XDOT(3,J) - RRD*ZJ)*FJ
   10 CONTINUE

      IF (IFIRST.EQ.1) GO TO 40
      J = NTOT
          A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*       Check the c.m. approximation.
          R2 = R(1)**2 + R(2)**2 + R(3)**2
          IF (RIJ2.GT.CMSEP2*R2.OR.LIST(1).EQ.0) GO TO 30
          CALL RESOLV
          J = 1
   20     A1 = X(1,J) - XI(1)
          A2 = X(2,J) - XI(2)
          A3 = X(3,J) - XI(3)
          DV(1) = XDOT(1,J) - XIDOT(1)
          DV(2) = XDOT(2,J) - XIDOT(2)
          DV(3) = XDOT(3,J) - XIDOT(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
   30     DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FIN(1) = FIN(1) + A1*DR3I
          FIN(2) = FIN(2) + A2*DR3I
          FIN(3) = FIN(3) + A3*DR3I
          FDN(1) = FDN(1) + (DV(1) - A1*DRDV)*DR3I
          FDN(2) = FDN(2) + (DV(2) - A2*DRDV)*DR3I
          FDN(3) = FDN(3) + (DV(3) - A3*DRDV)*DR3I
*       Include indirect terms with option #3.
   32     IF (KZ(3).EQ.0) GO TO 35
          IF (BODY(J).LT.M_CRIT) GO TO 35
          XJ = X(1,J)
          YJ = X(2,J)
          ZJ = X(3,J)
          RJ2 = XJ**2 + YJ**2 + ZJ**2
          FJ = BODY(J)/(RJ2*SQRT (RJ2))
          FIN(1) = FIN(1) - XJ*FJ
          FIN(2) = FIN(2) - YJ*FJ
          FIN(3) = FIN(3) - ZJ*FJ
          RRD = 3.0*(XJ*XDOT(1,J) + YJ*XDOT(2,J) + ZJ*XDOT(3,J))/RJ2
          FDN(1) = FDN(1) - (XDOT(1,J) - RRD*XJ)*FJ
          FDN(2) = FDN(2) - (XDOT(2,J) - RRD*YJ)*FJ
          FDN(3) = FDN(3) - (XDOT(3,J) - RRD*ZJ)*FJ
   35     IF (J.EQ.1) THEN
              J = J + 1
              GO TO 20
          END IF
*
*       Perform one iteration with dominant body of mass = 1.
 40   DO 50 ITER = 1,2
          RIN2 = 1.0/(XI(1)**2 + XI(2)**2 + XI(3)**2)
          RIN3 = RIN2*SQRT(RIN2)
          RD = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &                               XI(3)*XIDOT(3))*RIN2
*
*       Add solar contribution to the planetary perturbations.
          DO 42 K = 1,3
              FIRR(K) = FIN(K) - RIN3*XI(K)
              FD(K) = FDN(K) - (XIDOT(K) - RD*XI(K))*RIN3
   42     CONTINUE
*
      IF (TWOGG.GT.0.AND.I.EQ.3.AND.TIME.LE.T_S*TWOPI) GO TO 100
*
*       Add gas disk gravity
      IF (INT(G_P).EQ.1) THEN
         IF (BODY(I).GE.1E-4) THEN
                CALL GAS_POTENTIAL(XI,XIDOT,FIRR,FD)
         ELSE
                CALL GAS_POTENTIAL_FULL(XI, XIDOT, FIRR, FD)
         END IF
      END IF
      IF (INT(G_P).EQ.2) THEN
             IF (BODY(I).LE.1E-4) THEN
                CALL GAS_POTENTIAL_P(XI,XIDOT,FIRR,FD)
             END IF
             IF (BODY(I).GT.5E-4) THEN
                CALL GAS_POTENTIAL_J(XI, XIDOT, FIRR, FD)
             END IF
             IF (BODY(I).GT.1E-4.AND.BODY(I).LE.5E-4) THEN
                CALL GAS_POTENTIAL_S(XI, XIDOT, FIRR, FD)
             END IF
      END IF                  
*      
*       Add general relativity
      IF (G_R.GT.0.0) CALL GR(XI, XIDOT, FIRR, FD, I)
*
*       Add the tidal force by gas disk
      RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
      VI2 = XIDOT(1)**2 + XIDOT(2)**2 + XIDOT(3)**2
      RD = XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3)
      RI = SQRT(RI2)
      VI = SQRT(VI2)
      ZMB = 1.0 + BODY(I)
      SEMI = 2.0/RI - VI2/ZMB
      SEMI = 1/SEMI
      ECC2 = (1.0 - RI/SEMI)**2 + RD**2/(SEMI*ZMB)
      ECC = SQRT(ECC2)
      IF (G_D.GT.0) THEN
	IF (BODY(I).LT.1E-4) THEN
	   CALL DAMPING(XI, XIDOT, FIRR, FD, I)
 	END IF
      END IF
*
*      IF (G_D.GT.0.0.AND.(RI.LE.R_IN.OR.RI.GE.R_EDGE)) THEN 
*         IF (BODY(I).LT.1E-4) THEN
*                CALL DAMPING(XI, XIDOT, FIRR, FD, I)
*         END IF
*      END IF
*
*       Include the corrector after first or second iteration.
 100  DO 45 K = 1,3
       	      DF = 2.0*F(K,I) - FIRR(K)
	      FID = 6.0*FDOT(K,I)
	      SUM = FID + FD(K)
	      AT3(K) = 2.0*DF + DT*SUM
    	      BT2(K) = -3.0*DF - DT*(SUM + FID)
*
*       Replace standard corrector by Kokubo et al. (M.N. 297, 1067, 1998).
              FSUM = F(K,I) + 0.5*FIRR(K)
              DFD = FID - FD(K)
              XIDOT(K) = (DFD*DT12 + FSUM)*DT + X0DOT(K,I)
              V01 = 0.5*(X0DOT(K,I) + XIDOT(K))
*       Adopt \alpha = 7/6 as recommended by Kokubo & Makino (PASJ 56, 861).
              XI(K) = ((SUM*DT + 7.0*DF)*DT60 + V01)*DT + X0(K,I)
*
*             XI(K) = X0T(K) + (0.6*AT3(K) + BT2(K))*DTSQ12
*             XIDOT(K) = V0T(K) + (0.75*AT3(K) + BT2(K))*DT13
   45     CONTINUE
   50 CONTINUE
*
*       Copy corrected values and set new F, FDOT, D2 & D3.
      DO 60 K = 1,3
          X0(K,I) = XI(K)
          X0DOT(K,I) = XIDOT(K)
	  F(K,I) = 0.5*FIRR(K)
	  FDOT(K,I) = ONE6*FD(K)
          D3(K,I) = AT3(K)*DT6
          D2(K,I) = (3.0*AT3(K) + BT2(K))*DT2
*       NOTE: These are real derivatives!
   60 CONTINUE
*
*       Specify new time-step (standard criterion or fast expression).
      IF (KZ(5).EQ.0) THEN
          TTMP = TSTEP(FIRR,FD,D2(1,I),D3(1,I),ETA)
      ELSE
          TTMP = STEPI(FIRR,FD,D2(1,I),D3(1,I),ETA)
      END IF
      DT0 = TTMP
*
*       Select discrete value (increased by 2, decreased by 2 or unchanged).
      IF (TTMP.GT.2.0*STEP(I)) THEN
          IF (DMOD(TIME,2.0*STEP(I)).EQ.0.0D0) THEN 
              TTMP = MIN(2.0*STEP(I),DTMAX)
          ELSE
              TTMP = STEP(I) 
          END IF
      ELSE IF (TTMP.LT.STEP(I)) THEN
          TTMP = 0.5*STEP(I)
*       Check for further reduction by factor 2.
          IF (TTMP.GT.DT0) THEN
              TTMP = 0.5*TTMP
          END IF
      ELSE
          TTMP = STEP(I)
      END IF
*
*       Set new block step and update next time.
      STEP(I) = TTMP
      TNEXT(I) = STEP(I) + T0(I)
*
*       See whether any KS candidates are in the same block as body #I.
      IF (IKS.GT.0.AND.I.EQ.ICOMP) THEN
*       Accept same time, otherwise reduce STEP(ICOMP) and/or delay.
          IF (T0(JCOMP).EQ.T0(ICOMP)) THEN
              ICOMP = MIN(ICOMP,JCOMP)
              JCOMP = MAX(I,JCOMP)
          ELSE IF (T0(JCOMP) + STEP(JCOMP).LT.T0(ICOMP)) THEN
              STEP(ICOMP) = 0.5D0*STEP(ICOMP)
              TNEXT(ICOMP) = STEP(ICOMP) + T0(ICOMP)
              IKS = 0
          ELSE
              IKS = 0
          END IF
      END IF
*
*       Include optional diagnostic output for plotting.
      IF (NAME(I).EQ.KZ(8)) THEN
      IPLOT = IPLOT + 1
      IF (MOD(IPLOT,10).EQ.0.AND.STEP(I).GT.20.0*DTMIN) THEN
      J = KZ(8)
      RJ2 = 0.0
      VJ2 = 0.0
      DO 80 K = 1,3
          RJ2 = RJ2 + X(K,J)**2
          VJ2 = VJ2 + XDOT(K,J)**2
   80 CONTINUE
      SUNPL = 1.0 + BODY(J)
      SEMI = 2.0/SQRT(RJ2) - VJ2/SUNPL
      SEMI = 1.0/SEMI
      WRITE (10,85) TIME, SEMI, SQRT(RJ2)
   85 FORMAT (' ',F14.6,1P,2E14.6)
      END IF
      END IF
*
*     check if (escape/migrate too close to the host star) happen already
      RI2 = XI(1)**2 + XI(2)**2 + XI(3)**2
      VI2 = XIDOT(1)**2 + XIDOT(2)**2 + XIDOT(3)**2
      RD = XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3)
      RI = SQRT(RI2)
      VI = SQRT(VI2)
      ZMB = 1.0 + BODY(I)
      SEMI = 2.0/RI - VI2/ZMB
      SEMI = 1/SEMI
      ECC2 = (1.0 - RI/SEMI)**2 + RD**2/(SEMI*ZMB)
      ECC = SQRT(ECC2)
      V_ESC = SQRT(2/RI)
      IF ((ECC.GT.ECRIT.OR.RI.GT.R_ESC.OR.SEMI.LE.R_MHS).AND.KZ(9).GT.0) THEN
         WRITE (6,47) I, NAME(I), N-1, ECC, RI, SEMI
 47     FORMAT (' ESCAPE    I NAM N ECC R A ',3I4,2X,F9.4,2X,2F15.3)
          WRITE (0,*) TIME/TWOPI, NAME(I), ECC, RI, SEMI
          IESC = 1
      END IF
      RETURN
*
      END
