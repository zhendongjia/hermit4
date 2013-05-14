      SUBROUTINE ENERGY
*
*
*       Total energy.
*       -------------
*
      INCLUDE 'commonp.h'
      REAL*8  SUM(4),VSUN(3)
*
*
*       Calculate the potential energy.
      ZKIN = 0.0D0
      POT = 0.0D0
      I = 1
   20 JMIN = I + 1
      IF (I.LE.2*NPAIRS) THEN
          JMIN = IFIRST
      END IF
      POTJ = 0.0D0
*
      DO 30 J = JMIN,N
          A1 = X(1,I) - X(1,J)
          A2 = X(2,I) - X(2,J)
          A3 = X(3,I) - X(3,J)
          POTJ = POTJ + BODY(J)/SQRT(A1*A1 + A2*A2 + A3*A3)
   30 CONTINUE
*
      POT = POT + BODY(I)*POTJ
      I = I + 1
      IF (I.LT.N) GO TO 20
*
      DO 35 K = 1,4
          SUM(K) = 0.0
   35 CONTINUE
*
*       Sum the kinetic energy and solar interaction.
      DO 40 I = IFIRST,NTOT
          ZKIN = ZKIN + BODY(I)*(XDOT(1,I)**2 + XDOT(2,I)**2 +
     &                                          XDOT(3,I)**2)
          RI2 = X(1,I)**2 + X(2,I)**2 + X(3,I)**2
          POT = POT + BODY(I)/SQRT(RI2)
          SUM(4) = SUM(4) + BODY(I)
          DO 36 K = 1,3
              SUM(K) = SUM(K) + BODY(I)*XDOT(K,I)
   36     CONTINUE
   40 CONTINUE
      ZKIN = 0.5D0*ZKIN
*
*       Define velocity of Sun (M_sun = 1).
      VS2 = 0.0
      DO 45 K = 1,3
          VSUN(K) = -SUM(K)
          VS2 = VS2 + VSUN(K)**2
   45 CONTINUE
*
*       Include indirect terms if option is active.
      IF (KZ(3).GT.0) THEN
          ZKIN = ZKIN - 0.5*VS2*(1.0 - SUM(4))
      END IF
*
*       Obtain the binding energy of any regularized pair.
      IF (IFIRST.GT.1) THEN
          CALL RESOLV
      ELSE
          EBIN = 0.0
      END IF
*       Total energy = ZKIN - POT + EBIN.
*
      RETURN
*
      END
