      SUBROUTINE OUTPUT
*
*
*       Energy check and output.
*       ------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XS(3),XSD(3),SEMI(NMAX),ECC(NMAX),A(4)
      REAL*8  TA_F,W_F,W
      DATA ECRIT,RESC /0.99D0,5.0D0/
*
*
*       Predict X & XDOT for all particles.
      CALL XVPRED(IFIRST,NTOT)
*
*       Obtain the total energy at current time.
      CALL ENERGY
*
*       Initialize c.m. terms.
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
      DO 12 K = 1,4
          A(K) = 0.0
   12 CONTINUE
*
*       Obtain c.m. & angular momentum integrals and Z-moment of inertia.
      AZ = 0.0D0
      ZM = 0.0D0
      ZMASS = 0.0D0
      DO 20 I = 1,N
          ZMASS = ZMASS + BODY(I)
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
          DO 18 K = 1,2
              A(K) = A(K) + BODY(I)*X(K,I)
              A(K+2) = A(K+2) + BODY(I)*XDOT(K,I)
   18     CONTINUE
          AZ = AZ + BODY(I)*(X(1,I)*XDOT(2,I) - X(2,I)*XDOT(1,I))
   20 CONTINUE
*
*       Correct for indirect terms.
      IF (KZ(3).GT.0) THEN
          AZ = AZ - (A(1)*A(4) - A(2)*A(3))/(ZMASS + 1.0)
      END IF
*
*       Form c.m. coordinates & velocities (vectors & scalars).
      DO 25 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   25 CONTINUE
*
      CMR(4) = SQRT(CMR(1)**2 + CMR(2)**2 + CMR(3)**2)
      CMRDOT(4) = SQRT(CMRDOT(1)**2 + CMRDOT(2)**2 + CMRDOT(3)**2)
*
*       Define crossing time and save single particle energy.
      ETOT = ZKIN - POT
      TCR = ZMASS**2.5/(2.0*ABS(ETOT))**1.5
      ETOT = ETOT + EBIN
*
*       Update energies and form the relative error (divide by ZKIN or ETOT).
      IF (TIME.LE.0.0D0) THEN
          DE = 0.0D0
          BE(1) = ETOT
          BE(3) = ETOT
      ELSE
          BE(2) = BE(3)
          BE(3) = ETOT
          DE = BE(3) - BE(2)
          DETOT = DETOT + DE
          DE = DE/MAX(ZKIN,ABS(ETOT))
          ERRTOT = ERRTOT + DE
      END IF
*
*       Print diagnostic information.
      WRITE (6,30)  TIME/TWOPI, NSTEPS, NSTEPU, BE(3), DE, AZ
   30 FORMAT (/,' YRS =',1P,E9.1,'  # =',0P,I10,I8,'  E =',F10.6,
     &          '  DE =',1P,E10.2,'  AZ =',0P,F12.8)
   40 IESC = 0
*     IF (KZ(6).GT.0) THEN
          DO 50 I = IFIRST,N
              RI2 = 0.0
              VI2 = 0.0
              RD = 0.0
              DO 45 K = 1,3
                  XS(K) = X(K,I)
                  XSD(K) = XDOT(K,I)
                  RI2 = RI2 + XS(K)**2
                  VI2 = VI2 + XSD(K)**2
                  RD = RD + XS(K)*XSD(K)
   45         CONTINUE
              RI = SQRT(RI2)
              ZMB = 1.0 + BODY(I)
              SEMI(I) = 2.0/RI - VI2/ZMB
              SEMI(I) = 1.0/SEMI(I)
              ECC2 = (1.0 - RI/SEMI(I))**2 + RD**2/(SEMI(I)*ZMB)
              ECC(I) = SQRT(ECC2)
*       Consider escape removal (or high ECC) but exclude encounter.
          IF (ECC(I).GT.ECRIT.AND.STEP(I).GT.20.0*DTMIN.AND.
     &        KZ(9).GT.0) THEN
              IF (IESC.EQ.0.AND.(RI.GT.RESC.OR.ECC(I).LT.1.0)) THEN
                  WRITE (6,47) I, NAME(I), N-1, ECC(I), RI, SEMI(I)
   47             FORMAT (' ESCAPE    I NAM N ECC R A ',
     &                                3I4,F9.4,2F8.3)
                  IF (IESC.EQ.0) THEN
                      IESC = I
                      ECCI = ECC(I)
                  END IF
              END IF
          END IF
*       Include optional output of elements.
          IF (KZ(6).GT.0.AND.IESC.EQ.0.AND.STEP(I).GT.20.0*DTMIN) THEN
             W_F = 0
             TA_F = 0
             W = GET_PRECESSION(X(:,I), XDOT(:,I), SEMI(I), ECC(I),
     &            W_F, TA_F)
             WRITE (6,48) I, NAME(I), ECC(I), RI, SEMI(I), STEP(I),
     &            W, W_F, TA_F, X(1,I), X(2,I), XDOT(1,I), XDOT(2,I),
     &            T_TIDAL1(I), T_TIDAL2(I)
   48         FORMAT (' ORBIT    I NAM ECC R A S W', 2I4, 2X, F15.4, 2X,
     &            F15.4, 2X, F15.4, 1P, E10.2, 0P, 3F15.3, 6E15.3)
          END IF
   50     CONTINUE
*     END IF
*
*       Check optional escaper removal (also ECC > ECRIT).
      IF (KZ(9).GT.0.AND.IESC.GT.0) THEN
          CALL REMOVE(IESC)
          IF (N.EQ.1) STOP
*       Update the total energy.
          CALL ENERGY
          BE(3) = ZKIN - POT + EBIN
*       Reduce the steps to compensate for force discontinuity.
          DO 55 J = IFIRST,NTOT
              STEP(J) = 0.25*STEP(J)
              TNEXT(J) = T0(J) + STEP(J)
   55     CONTINUE
*       Try again just in case.
          GO TO 40
      END IF
*
*       Include optional output of individual bodies.
      IF (KZ(4).GT.0) THEN
          CALL BODIES
      END IF
*
      TPRINT = TPRINT + DELTAT
      CALL CPUTIM(TCOMP)
      CPUTOT = CPUTOT + TCOMP - CPU0
      CPU0 = TCOMP
*
*       Save COMMON after energy check.
      IF (KZ(2).GE.1) CALL MYDUMP(1,2)
*
*       Check termination criterion.
      IF (TIME.GE.TCRIT) THEN
*       Terminate after optional COMMON save.
          WRITE (6,60)  TIME, CPUTOT/60.0, ERRTOT, DETOT
   60     FORMAT (//,9X,'END RUN',3X,'TIME =',F7.1,'  CPUTOT =',F7.1,
     &                  '  ERRTOT =',1P,E10.2,'  DETOT =',E10.2)
          IF (KZ(1).GT.0) CALL MYDUMP(1,1)
          STOP
      END IF
*
      RETURN
*
      END
