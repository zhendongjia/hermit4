      SUBROUTINE STEPS(I1,I2)
*
*
*       Initialization of time-steps & prediction variables.
*       ----------------------------------------------------
*
      INCLUDE 'commonp.h'
*
*
*       Set new steps and initialize prediction variables.
      DO 40 I = I1,I2
*
*       Obtain time-steps using DT = ETA*F/FDOT (re-instated 9/98).
          FI2 = F(1,I)**2 + F(2,I)**2 + F(3,I)**2
          FD2 = FDOT(1,I)**2 + FDOT(2,I)**2 + FDOT(3,I)**2
*       Include precaution for small velocities (i.e. DT = ETA*TCR).
          IF (FD2.LT.FI2) FD2 = FI2/TCR**2
          DT = ETA*SQRT(FI2/FD2)
          IF (TIME.EQ.0.0D0) THEN
              DT = 0.1*DT
          ELSE IF (I.GT.N) THEN
              DT = 0.5*DT
          END IF
*
*       Initialize the time and obtain quantized step.
          T0(I) = TIME
*
*       Convert predicted step to nearest block time-step (truncated down).
          CALL STEPK(DT,DTN)
          IF (TIME.LE.0.0D0) THEN
              STEP(I) = DTN
          ELSE 
*       Reduce steps by factor 2 until commensurate with current time.
              STEP(I) = DTN
              ITER = 0
   10         IF (DMOD(TIME,STEP(I)).NE.0.0D0) THEN
                  STEP(I) = 0.5D0*STEP(I)
                  ITER = ITER + 1
                  IF (ITER.LT.16.OR.STEP(I).GT.DTK(40)) GO TO 10
                  STEP(I) = DTK(40)
                  WRITE (6,15)  I, ITER, TIME/STEP(I), DT, STEP(I)
   15             FORMAT (' WARNING!   I ITER T/STEP DT STEP ',
     &                                 I5,I4,F16.4,1P,2E9.1)
              END IF
          END IF
*
*       Initialize or update array for new block times.
          TNEXT(I) = T0(I) + STEP(I)
*
*       Set prediction variables (X0DOT set by START).
          DO 30 K = 1,3
              X0(K,I) = X(K,I)
              F(K,I) = 0.5D0*F(K,I)
              FDOT(K,I) = ONE6*FDOT(K,I)
              D2(K,I) = 0.0D0
              D3(K,I) = 0.0D0
   30     CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
