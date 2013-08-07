      SUBROUTINE REMOVE(I)
*
*
*       Particle removal.
*       -----------------
*
      INCLUDE 'commonp.h'
      REAL*8  A(6)
*
*
*       Move up all COMMON variables uless body #I is last.
      IF (I.GT.NTOT) GO TO 20
*
      DO 10 J = I,NTOT
          J1 = J + 1
          DO 5 K = 1,3
              X(K,J) = X(K,J1)
              X0(K,J) = X0(K,J1)
              X0DOT(K,J) = X0DOT(K,J1)
              XDOT(K,J) = XDOT(K,J1)
              F(K,J) = F(K,J1)
              FDOT(K,J) = FDOT(K,J1)
              D1(K,J) = D1(K,J1)
              D2(K,J) = D2(K,J1)
              D3(K,J) = D3(K,J1)
    5     CONTINUE
*
          BODY(J) = BODY(J1)
          NAME(J) = NAME(J1)
          RADIUS(J) = RADIUS(J1)
          STEP(J) = STEP(J1)
          T0(J) = T0(J1)
          TNEXT(J) = TNEXT(J1)
   10 CONTINUE
*
*       Reduce particle number.
      NTOT = NTOT - 1
      N = N - 1
      NMASS = NMASS - 1
*
   20 RETURN
*
      END
