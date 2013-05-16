      SUBROUTINE SEARCH(I,IKS)
*
*
*       Close encounter search.
*       -----------------------
*
      INCLUDE 'commonp.h'
*
*
*       Increase counter for regularization attempts.
      NKSTRY = NKSTRY + 1
      FMAX = 0.0
      NCLOSE = 0
      JCOMP = 0
*
*       Find dominant neighbour by selecting all STEP(J) < 2*DTMIN.
      RJMIN2 = 1.0
      DO 6 J = IFIRST,NMASS
          IF (J.EQ.I.OR.STEP(J).GT.DTMIN) GO TO 6
          IF (BODY(J).LT.1E-6) GO TO 6
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          IF (RIJ2.LT.RMIN2) THEN
              NCLOSE = NCLOSE + 1
              JLIST(NCLOSE) = J
*       Remember index of every single body with small step inside 2*RMIN.
              FIJ = (BODY(I) + BODY(J))/RIJ2
              IF (FIJ.GT.FMAX) THEN
                  FMAX = FIJ
*       Save square distance and global index of dominant body.
                  RJMIN2 = RIJ2
                  JCOMP = J
              END IF
          END IF
    6 CONTINUE
*
      IF (RJMIN2.GT.RMIN2) GO TO 10
*
      RD = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Only select approaching particles (include nearly circular case).
      RIJMIN = SQRT(RJMIN2)
      IF (RD.GT.0.02*SQRT((BODY(I) + BODY(JCOMP))*RIJMIN)) GO TO 10
*
*       Evaluate vectorial perturbation due to the close bodies.
      CALL FPERT(I,JCOMP,NCLOSE,PERT)
*
*       Accept #I & JCOMP if the relative motion is dominant (GI < 0.2).
      GI = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      IF (GI.GT.0.2) GO TO 10
*
*       Save index and increase indicator to denote new regularization.
      ICOMP = I
      IKS = IKS + 1
*
   10 RETURN
*
      END
