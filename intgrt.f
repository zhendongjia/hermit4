      SUBROUTINE INTGRT
*
*
*       N-body integrator flow control.
*       -------------------------------
*
      INCLUDE 'commonp.h'
      INTEGER  NXTLST(NMAX),LISTQ(NMAX),NL(40)
      LOGICAL LOOP
      SAVE IQ,LQ,LOOP
      DATA IQ,LQ,LOOP /0,11,.TRUE./
*
*
*       Enforce level search on return, except new and terminated KS.
      IF (IPHASE.NE.1.AND.IPHASE.NE.2) LOOP = .TRUE.
*
*       (Re-)Initialize end-point of integration times and set DTM.
    1 DTM = 1.0
      DO 2 I = IFIRST,NTOT
          TNEXT(I) = T0(I) + STEP(I)
          DTM = MIN(DTM,STEP(I))
    2 CONTINUE
      IPHASE = 0
      IKS = 0
      IQ = 0
*
*       Determine level for the smallest step (ignore extreme values).
      LQS = 40
      DO 3 L = 1,40
          IF (DTM.EQ.DTK(L)) THEN
              LQS = L
          END IF
    3 CONTINUE
*
*       Specify upper level for optimized membership.
      LQB = MAX(LQS - 4,1)
*       Enforce new block step search initially and on significant change.
      TLISTQ = TIME
*
*       Check updating new list of block steps with T0 + STEP =< TLISTQ.
    5 IF (TIME.GE.TLISTQ) THEN
*       Update interval by optimization at major times using sqrt(N).
          IF (DMOD(TLISTQ,2.0D0).EQ.0.0D0.OR.LOOP) THEN
              LOOP = .FALSE.
              DO 10 L = 1,40
                  NL(L) = 0
   10         CONTINUE
              DO 14 I = IFIRST,NTOT
*       Count steps at five different levels for the smallest values.
                  DO 12 L = LQB,LQS
                      IF (STEP(I).LT.DTK(L)) NL(L) = NL(L) + 1
   12             CONTINUE
   14         CONTINUE
              NLSUM = 0
*       Determine interval by summing smallest steps until near sqrt(N).
              NSQ = SQRT(FLOAT(N))
              LQ = LQS
              DO 15 L = LQS,LQB,-1
                  NLSUM = NLSUM + NL(L)
                  IF (NLSUM.LE.NSQ) LQ = L
   15         CONTINUE
*             WRITE (6,16)  TIME+TOFF,NQ,NLSUM,LQ,(NL(K),K=LQB,LQS)
*  16         FORMAT (' LEVEL CHECK:    T NQ NLSUM LQ NL  ',
*    &                                  F9.3,3I5,2X,7I4)
          END IF
*
*       Increase interval by optimized value.
          NQ = 0
          TMIN = 1.0D+10
   18     TLISTQ = TLISTQ + DTK(LQ)
          DO 20 I = IFIRST,NTOT
              IF (TNEXT(I).LE.TLISTQ) THEN
                  NQ = NQ + 1
                  LISTQ(NQ) = I
                  TMIN = MIN(TNEXT(I),TMIN)
              END IF
   20     CONTINUE
*       Increase interval in rare case of zero membership.
          IF (NQ.EQ.0) GO TO 18
*       Make a slight adjustment for high levels and small membership.
          IF (LQ.LE.15.AND.NQ.LE.2) LQ = MAX(LQ - 1,4)
      END IF
*
*       Find all particles in next block (TNEXT = TMIN).
      CALL INEXT(NQ,LISTQ,TMIN,NXTLEN,NXTLST)
*
*       Set new time and save block time (for regularization terminations).
      I = NXTLST(1)
      TIME = T0(I) + STEP(I)
      TBLOCK = TIME
      LI = 0
*
*     Re-determine list if current time exceeds boundary.
      IF (TIME.GT.TLISTQ) GO TO 5
*
*       Include commensurability test (may be suppressed if no problems).
*     IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
*         WRITE (6,30) I, NAME(I), NSTEPS, TIME, STEP(I), TIME/STEP(I)
*  30     FORMAT (' DANGER!   I NM # TIME STEP T/DT ',
*    &                        2I5,I11,F12.5,1P,E9.1,0P,F16.4)
*         STOP
*     END IF
*
*       Check for new regularization at end of block.
      IF (IKS.GT.0) THEN
          TIME = TPREV
          IPHASE = 1
          GO TO 100
      END IF
*
*       Check next output time.
      IF (TIME.GT.TPRINT) THEN
          TIME = TPREV
          IPHASE = 3
          GO TO 100
      END IF
*
*       See whether to advance regularized solution at first new time.
      IF (TIME.GT.TPREV.AND.NPAIRS.GT.0) THEN
          CALL BHINT(IQ)
          IF (IQ.LT.0) GO TO 1
      END IF
*
*       Predict all coordinates & velocities to order FDOT.
      DO 40 J = IFIRST,NTOT
          S = TIME - T0(J)
          S1 = 1.5*S
          S2 = 2.0*S
          X(1,J) = ((FDOT(1,J)*S + F(1,J))*S + X0DOT(1,J))*S + X0(1,J)
          X(2,J) = ((FDOT(2,J)*S + F(2,J))*S + X0DOT(2,J))*S + X0(2,J)
          X(3,J) = ((FDOT(3,J)*S + F(3,J))*S + X0DOT(3,J))*S + X0(3,J)
          XDOT(1,J) = (FDOT(1,J)*S1 + F(1,J))*S2 + X0DOT(1,J)
          XDOT(2,J) = (FDOT(2,J)*S1 + F(2,J))*S2 + X0DOT(2,J)
          XDOT(3,J) = (FDOT(3,J)*S1 + F(3,J))*S2 + X0DOT(3,J)
   40 CONTINUE
*
*       Save new time and increase # blocks.
      TPREV = TIME
      NBLOCK = NBLOCK + 1
      TMIN = 1.0D+10
*
*       Advance the pointer (<= NXTLEN) and select next particle index.
   50 LI = LI + 1
      IF (NXTLEN.EQ.0) GO TO 1
      IF (LI.GT.NXTLEN) GO TO 5
*
*       Select the next sequential particle.
      I = NXTLST(LI)
      TIME = T0(I) + STEP(I)
*
*       Perform the integration step for body #I.
      CALL NBINT(I,IKS)
      IF (IESC.GT.0) THEN 
         CALL REMOVE(I)
         CALL RM_FROM_LIST(I, LISTQ, NQ)
         CALL RM_FROM_LIST(I, NXTLST, NXTLEN)
         LI = LI - 1
         IF (N.LE.NLEAST) STOP
         CALL ENERGY
         BE(3) = ZKIN - POT + EBIN
          DO 55 J = IFIRST,NTOT
             STEP(J) = 0.25*STEP(J)
             TNEXT(J) = T0(J) + STEP(J)
 55          CONTINUE
            GO TO 50
        END IF         
*
*       Determine next block time.
      TMIN = MIN(TNEXT(I),TMIN)
*
      IF (LI.EQ.NXTLEN) THEN
          DO 35 L = 1,NXTLEN
              I = NXTLST(L)
              DO 32 K = 1,3
                  X(K,I) = X0(K,I)
                  XDOT(K,I) = X0DOT(K,I)
   32         CONTINUE
   35     CONTINUE
      END IF
*
*       Check termination of regularization.
      IF (IQ.NE.0.AND.LI.EQ.NXTLEN) THEN
          IPHASE = 2
          GO TO 100
      END IF
*
*       Increase counters and check timer & optional COMMON save.
      NSTEPS = NSTEPS + 1
      NTIMER = NTIMER + 1
      IF (NTIMER.LT.100000) GO TO 50
*
      NTIMER = 0
*     IF (NSTEPS.GE.100000) THEN
*         NSTEPS = 0
*         IF (KZ(1).GT.1) CALL MYDUMP(1,1)
*     END IF
*
*       Include facility for termination of run (create dummy file STOP).
      OPEN (99,FILE='STOP',STATUS='OLD',FORM='FORMATTED',IOSTAT=IO)
      IF (IO.EQ.0) THEN
          CLOSE (99)
          WRITE (6,70)
   70     FORMAT  (/,9X,'TERMINATION BY MANUAL INTERVENTION')
          CPU = 0.0
      END IF
*
*       Repeat cycle until elapsed computing time exceeds the limit.
      CALL CPUTIM(TCOMP)
      IF (TCOMP.LT.CPU) GO TO 50
*
*       Terminate run with optional COMMON save.
      IF (KZ(1).GT.0) THEN
          CPUTOT = CPUTOT + TCOMP - CPU0
          CALL MYDUMP(1,1)
          WRITE (6,80)  TIME/TWOPI, TCOMP, CPUTOT/60.0, ERRTOT, DETOT
   80     FORMAT (/,9X,'COMMON SAVED AT YRS =',1P,E9.1,
     &                 '  TCOMP =',0P,F7.1,
     &                 '  CPUTOT =',F6.1,'  ERRTOT =',F10.6,
     &                 '  DETOT =',F10.6)
      END IF
*
      STOP
*
  100 RETURN
*
      END
