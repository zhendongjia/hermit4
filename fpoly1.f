      SUBROUTINE FPOLY1(I1,I2)
*
*
*       Force & first derivative.
*       -------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  A(9),F1(3),F1DOT(3),XI(3),XIDOT(3)
*
*
*       Loop over all bodies or one single body.
      DO 40 I = I1,I2
*
*       Initialize force & first derivative for body #I.
      DO 10 K = 1,3
          F1(K) = 0.0
          F1DOT(K) = 0.0
          XI(K) = X(K,I)
          XIDOT(K) = XDOT(K,I)
   10 CONTINUE
*
*       Obtain force & first derivative by summing over all bodies.
      DO 30 J = IFIRST,NMASS
          IF (J.EQ.I) GO TO 30
          DO 15 K = 1,3
              A(K) = X(K,J) - X(K,I)
              A(K+3) = XDOT(K,J) - XDOT(K,I)
   15     CONTINUE
*
          A(7) = 1.0/(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
          A(8) = BODY(J)*A(7)*SQRT(A(7))
          A(9) = 3.0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))*A(7)
*
          IF (KZ(3).EQ.0) GO TO 18
*       Include indirect terms.
          XJ = X(1,J)
          YJ = X(2,J)
          ZJ = X(3,J)
          RJ2 = XJ**2 + YJ**2 + ZJ**2
          FJ = BODY(J)/(RJ2*SQRT (RJ2))
          F1(1) = F1(1) - XJ*FJ
          F1(2) = F1(2) - YJ*FJ
          F1(3) = F1(3) - ZJ*FJ
          RRD = 3.0*(XJ*XDOT(1,J) + YJ*XDOT(2,J) + ZJ*XDOT(3,J))/RJ2
          F1DOT(1) = F1DOT(1) - (XDOT(1,J) - RRD*XJ)*FJ
          F1DOT(2) = F1DOT(2) - (XDOT(2,J) - RRD*YJ)*FJ
          F1DOT(3) = F1DOT(3) - (XDOT(3,J) - RRD*ZJ)*FJ
   18     DO 20 K = 1,3
              F1(K) = F1(K) + A(K)*A(8)
              F1DOT(K) = F1DOT(K) + (A(K+3) - A(K)*A(9))*A(8)
   20     CONTINUE
   30 CONTINUE
*
*       Add solar terms including body #I.
      SUNPL = 1.0 + BODY(I)
      RIN2 = 1.0/(XI(1)**2 + XI(2)**2 + XI(3)**2)
      RIN3 = SUNPL*RIN2*SQRT(RIN2)
      RD = 3.0*(XI(1)*XIDOT(1) + XI(2)*XIDOT(2) + XI(3)*XIDOT(3))*RIN2
      DO 32 K = 1,3
          F1(K) = F1(K) - RIN3*XI(K)
          F1DOT(K) = F1DOT(K) - (XIDOT(K) - RD*XI(K))*RIN3
   32 CONTINUE
*
*       Save F & FDOT for body #I (note D1 not used).
      DO 35 K = 1,3
          F(K,I) = F1(K)
          D1(K,I) = F1DOT(K)
          FDOT(K,I) = F1DOT(K)
   35 CONTINUE
   40 CONTINUE
*
      RETURN
*
      END
