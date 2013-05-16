      SUBROUTINE BODIES
*
*
*       Output of single bodies or binaries.
*       ------------------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  A(3)
*
*
*       Consider all interacting particles (singles or binaries).
      SMAX = 0.01
      DO 40 I = 1,NMASS
          IF (KZ(4).EQ.1) THEN
              WRITE (6,20)  NAME(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
   20         FORMAT (' BODIES    NAM X XD ',
     &             I5,2X,F15.3,2X,F15.3,2X,F15.3,2X,3F7.2)
              GO TO 40
          END IF
          JMIN = 0
          RJMIN2 = 2.0
          DO 30 J = I+1,NMASS
              A1 = X(1,I) - X(1,J)
              A2 = X(2,I) - X(2,J)
              A3 = X(3,I) - X(3,J)
              RIJ2 = A1**2 + A2**2 + A3**2
              IF (RIJ2.LT.RJMIN2) THEN
                  RJMIN2 = RIJ2
                  JMIN = J
              END IF
   30     CONTINUE
          IF (JMIN.LE.I) GO TO 40
          RIJMIN = SQRT(RJMIN2)
          VR2 = (XDOT(1,I) - XDOT(1,JMIN))**2 +
     &          (XDOT(2,I) - XDOT(2,JMIN))**2 +
     &          (XDOT(3,I) - XDOT(3,JMIN))**2
          EREL = 0.5*VR2 - (BODY(I) + BODY(JMIN))/RIJMIN
*       Only print significant binaries.
          IF (EREL.GT.-0.1) GO TO 40
          SEMI = -0.5*(BODY(I) + BODY(JMIN))/EREL
          ZN = SQRT((BODY(I) + BODY(JMIN))/SEMI**3)
          RD = (X(1,I) - X(1,JMIN))*(XDOT(1,I) - XDOT(1,JMIN)) +
     &         (X(2,I) - X(2,JMIN))*(XDOT(2,I) - XDOT(2,JMIN)) +
     &         (X(3,I) - X(3,JMIN))*(XDOT(3,I) - XDOT(3,JMIN))
          ECC2 = (1.0 - RIJMIN/SEMI)**2 +
     &                             RD**2/(SEMI*(BODY(I) + BODY(JMIN)))
          ECC = SQRT(ECC2)
          RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
          WRITE (6,35)  NAME(I), NAME(JMIN),
     &                  SEMI, ZN, RIJMIN, RI, ECC
   35     FORMAT ('   BINARY ',2I4,2F8.4,F9.4,F9.4,F9.4)
   40 CONTINUE
*
   70 RETURN
*
      END
