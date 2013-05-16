      SUBROUTINE BHINT(IQ)
*
*
*       Integration of Burdet-Heggie equations.
*       ---------------------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XI(6),VI(6),R0D2(3),R0D3(3),B0D(3),B0D2(3),FP(6),FD(6)
      REAL*8  RP(3),RPDOT(3)
      LOGICAL ITERM,ICOLL
*       Define collision distance and merged binary period (days).
      DATA  RCOLL,TMERGE /3.0D-04,10.0/
*
*
*       Define logical variables for termination and collision at TBLOCK.
      ITERM = .FALSE.
      ICOLL = .FALSE.
*
*       Advance the time until TBLOCK exceeded.
    1 TIME = T0(1) + STEP(1)
      IF (TIME.GT.TBLOCK) GO TO 50
*
*       Predict all BH elements to high order (R and RDOT to RDOT5).
      DT3 = ONE6*DTAU**3
      DT4 = 0.25*DT3*DTAU
      DO 5 K = 1,3
          RP(K) = ((RDOT3(K)*DTAU + RDOT2(K))*DTAU +
     &                             RDOT(K))*DTAU + R(K)
          RPDOT(K) = (3.0*RDOT3(K)*DTAU + 2.0*RDOT2(K))*DTAU + RDOT(K)
          B(K) = (0.5*BDOT2(K)*DTAU + BDOT(K))*DTAU + B(K)
          R(K) = RP(K) + (0.2*RDOT5(K)*DTAU + RDOT4(K))*DT4
          RDOT(K) = RPDOT(K) + (0.25*RDOT5(K)*DTAU + RDOT4(K))*DT3
    5 CONTINUE
      PPRED = (0.5*PDOT2*DTAU + PDOT)*DTAU + P
      P = PPRED + (0.25*PDOT4*DTAU + PDOT3)*DT3
*
*       Predict c.m. and perturbers and copy c.m. values.
      CALL RESOLV
      DO 10 K = 1,3
          XI(K) = X(K,ICOMP)
          XI(K+3) = X(K,JCOMP)
          VI(K) = XDOT(K,ICOMP)
          VI(K+3) = XDOT(K,JCOMP)
   10 CONTINUE
*
*       Save lower derivatives for corrector.
      P0D = PDOT
      P0D2 = PDOT2
      DO 15 K = 1,3
          R0D2(K) = 2.0*RDOT2(K)
          R0D3(K) = 6.0*RDOT3(K)
          B0D(K) = BDOT(K)
          B0D2(K) = BDOT2(K)
   15 CONTINUE
*
*       Obtain perturbation and regularized derivatives.
      CALL BHPERT(XI,VI,FP,FD)
*
*       Form time-step scalars for the corrector.
      DTU = DTAU
      DTSQ = DTU**2
      DT6 = 6.0/(DTU*DTSQ)
      DT2 = 2.0/DTSQ
      DT13 = ONE3*DTU
      DTSQ12 = ONE12*DTSQ
*
*       Apply the corrector for R(K).
      DO 20 K = 1,3
          DF = R0D2(K) - RDOT2(K)
          SUM = R0D3(K) + RDOT3(K)
          AT3 = 2.0D0*DF + SUM*DTU
          BT2 = -3.0D0*DF - (SUM + R0D3(K))*DTU
          R(K) = RP(K) + (0.6D0*AT3 + BT2)*DTSQ12
          RDOT(K) = RPDOT(K) + (0.75D0*AT3 + BT2)*DT13
          RDOT2(K) = 0.5D0*RDOT2(K)
          RDOT3(K) = ONE6*RDOT3(K)
          RDOT4(K) = (3.0D0*AT3 + BT2)*DT2
          RDOT5(K) = AT3*DT6
   20 CONTINUE
*
*       Include the corrector for B(K) (as in energy corrector).
      DO 25 K = 1,3
          DF = B0D(K) - BDOT(K)
          SUM = B0D2(K) + BDOT2(K)
          AT3 = 2.0D0*DF + SUM*DTU
          BT2 = -3.0D0*DF - (SUM + B0D2(K))*DTU
          B(K) = B(K) + (0.75D0*AT3 + BT2)*DT13
          BDOT3(K) = (3.0D0*AT3 + BT2)*DT2
          BDOT4(K) = AT3*DT6
   25 CONTINUE
*
*       Correct binding energy and save higher derivatives.
      DPD = P0D - PDOT
      SUM = P0D2 + PDOT2
      AT3 = 2.0D0*DPD + SUM*DTU
      BT2 = -3.0D0*DPD - (SUM + P0D2)*DTU
      P = PPRED + (0.75D0*AT3 + BT2)*DT13
      PDOT3 = (3.0D0*AT3 + BT2)*DT2
      PDOT4 = AT3*DT6
*
*       Prepare evaluation of regularized time derivatives.
   28 RR = 0.0
      RPR = 0.0
      DO 30 K = 1,3
          RR = RR + R(K)**2
          RPR = RPR + R(K)*RDOT(K)
   30 CONTINUE
      RR = SQRT(RR)
      RPR = RPR/RR
      TDOT2 = RPR
      TDOT3 = P*RR + BODY(NTOT)
      TDOT4 = PDOT*RR + P*RPR
*       Include all but the R*F'' term in TDOT5.
      TDOT5 = PDOT2*RR + 2.0*PDOT*RPR + P*TDOT3
*
*       Construct higher time derivatives (F'=R*dF/dt).
      DO 35 K = 1,3
          TDOT3 = TDOT3 + R(K)*FP(K)*RR
          TDOT4 = TDOT4 + (RDOT(K)*FP(K) + R(K)*RR*FD(K))*RR +
     &                                     R(K)*FP(K)*RPR
          TDOT5 = TDOT5 + (RDOT2(K)*FP(K) + 2.0*RR*RDOT(K)*FD(K))*RR
     &                  + 2.0*(RDOT(K)*FP(K) + RR*R(K)*FD(K))*RPR
     &                  + R(K)*FP(K)*TDOT3
   35 CONTINUE
*
*       Prescribe new time-step (only for binary) and update STEP & T0.
      IF (P.LT.0.0) THEN
          A2 = MIN(RMIN/BODY(NTOT),1.0/ABS(P))
          DTAU = ETAU*SQRT(A2)/(1.0 + 1000.0*GAMMA)**0.333
      END IF
      STEP(1) = ((((TDOT5*0.2*DTAU + TDOT4)*DTAU/24.0 +
     &              ONE6*TDOT3)*DTAU + 0.5*TDOT2)*DTAU + RR)*DTAU
      T0(1) = TIME
      NSTEPU = NSTEPU + 1
*
*       Activate logical variable on termination condition.
      IF (RR.GT.RR0) THEN
          ITERM = .TRUE.
      END IF
*
*       Include inelastic collision for periods < TMERGE days.
      IF (P.LT.0.0.AND.KZ(7).GT.0) THEN
          SEMI = -BODY(NTOT)/P
          TK = 365.0*TWOPI*SQRT(SEMI**3/BODY(NTOT))
      ELSE
          TK = 1000.0
      END IF
*
*       Delay collision orbit until end of block-step.
      IF ((RR.LT.RCOLL.OR.TK.LT.TMERGE).AND..NOT.ICOLL.AND.
     &    KZ(7).GT.0) THEN
          DTR = TBLOCK - TIME
          IF (DTR.LT.STEP(1)) THEN
              A2 = 1.0/RR
              A3 = A2*DTR
              A4 = 3.0*TDOT2**2*A2 - TDOT3
              DTU = ((ONE6*A4*A3 - 0.5*TDOT2)*A2*A3 + 1.0)*A3
              DTAU = MIN(DTU,DTAU)
              STEP(1) = TBLOCK - T0(1)
              ICOLL = .TRUE.
              GO TO 1
          END IF
      END IF
*
*       Check termination near end of block-step.
      IF (ITERM) THEN
          DTR = TBLOCK - TIME
          IF (DTR.LE.STEP(1)) THEN
              IF (DTR.LE.0.0D0) THEN
                  IQ = 1
                  GO TO 50
              ELSE
*       Invert remaining interval for integration to TIME = TBLOCK.
                  A2 = 1.0/RR
                  A3 = A2*DTR
                  A4 = 3.0*TDOT2**2*A2 - TDOT3
                  DTU = ((ONE6*A4*A3 - 0.5*TDOT2)*A2*A3 + 1.0)*A3
                  DTAU = MIN(DTU,DTAU)
                  STEP(1) = TBLOCK - T0(1)
                  GO TO 1
              END IF
          END IF
      END IF
*
*       Check for optional collision or coalescence of close binary.
      IF ((ICOLL.AND.(RR.LT.RCOLL.OR.TK.LT.TMERGE)).AND.
     &    KZ(7).GT.0) THEN
          ZMU = BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
          EB = 0.5*ZMU*P
*       Subtract current binding energy for conservation.
          BE(3) = BE(3) - EB
          WRITE (6,40)  TIME/TWOPI, NAME(ICOMP), NAME(JCOMP), N-1, RR,
     &                  EB, TK
   40     FORMAT (' COLLISION    YRS NAM N R EB TK ',
     &                           F13.5,3I4,1P,E10.2,E12.4,0P,F9.2)
          IF (N.LE.2) STOP
*       Ensure name of heaviest body for progenitor.
          I1 = ICOMP
          IF (BODY(JCOMP).GT.BODY(ICOMP)) I1 = JCOMP
          NAME(NTOT) = NZERO + NAME(I1)
*       Replace regularized components with current c.m.
          CALL REMOVE(JCOMP)
          CALL REMOVE(ICOMP)
*       Compensate for particle reduction by 2 and reset IFIRST.
          N = N + 1
          NMASS = NMASS + 1
          NPAIRS = NPAIRS - 1
          IFIRST = 2*NPAIRS + 1
*
*       Improve the solutions and total energy if perturber is present.
          IF (LIST(1).EQ.0) GO TO 48
*       Predict all bodies and initialize forces & time-steps of perturbers.
          CALL XVPRED(IFIRST,N)
          NP = LIST(1) + 1
*       Include improvement of new body (add and subtract 2).
          LIST(NP + 1) = N + 2
          DO 45 L = 2,NP+1
*       Note that perturbers have been moved up by 2.
              J = LIST(L) - 2
              DO 42 K = 1,3
                  X0DOT(K,J) = XDOT(K,J)
   42         CONTINUE
              CALL FPOLY1(J,J)
              CALL STEPS(J,J)
   45     CONTINUE
*
*       Define new total energy (suppressed; may not be necessary).
*         CALL ENERGY
*         BE(3) = ZKIN - POT
*
*       Return to start of INTGRT with current c.m. as new body.
   48     IQ = -1
          GO TO 50
      END IF
*
*       Continue until end of block-step even for termination.
      GO TO 1
*
*       Restore current block-time.
   50 TIME = TBLOCK
*
      RETURN
*
      END
