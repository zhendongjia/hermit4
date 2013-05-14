      SUBROUTINE RESOLV
*
*
*       Physical coordinates and velocities.
*       ------------------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  Y(3),V(3)
*
*
*       Expand regularized interval to third order.
      RR = SQRT(R(1)**2 + R(2)**2 + R(3)**2)
      A2 = 1.0/RR
      A3 = A2*(TIME - T0(1))
      A4 = 3.0*TDOT2**2*A2 - TDOT3
      DTU = ((ONE6*A4*A3 - 0.5*TDOT2)*A2*A3 + 1.0)*A3
      DTU = MIN(DTU,DTAU)
*
*       Predict relative values and convert velocity to physical units.
      RR = 0.0
      DO 10 K = 1,3
          Y(K) = ((RDOT3(K)*DTU + RDOT2(K))*DTU + RDOT(K))*DTU + R(K)
          V(K) = (3.0*RDOT3(K)*DTU + 2.0*RDOT2(K))*DTU + RDOT(K)
          RR = RR + Y(K)**2
   10 CONTINUE
      RR = SQRT(RR)
      DO 12 K = 1,3
          V(K) = V(K)/RR
   12 CONTINUE
*
*       Update perturbers and c.m. to order FDOT2.
      NP1 = LIST(1) + 2
      LIST(NP1) = NTOT
      NP1 = LIST(1) + 1
      DO 25 L = 2,NP1
          J = LIST(L)
          DT = TIME - T0(J)
          DT3 = ONE6*DT
          DT4 = 0.5*DT3
          DO 20 K = 1,3
              X(K,J) = (((D2(K,J)*DT4 + FDOT(K,J))*DT + F(K,J))*DT +
     &                                  X0DOT(K,J))*DT + X0(K,J)
              XDOT(K,J) = ((D2(K,J)*DT3 + 3.0*FDOT(K,J))*DT +
     &                                    2.0*F(K,J))*DT + X0DOT(K,J)
*         X(K,J) = X(K,J) + D2(K,J)*DT4
*         XDOT(K,J) = XDOT(K,J) + D2(K,J)*DT3
   20     CONTINUE
   25 CONTINUE
*
*       Obtain the global coordinates and velocities.
      DO 30 K = 1,3
          X(K,ICOMP) = X(K,NTOT) + BODY(JCOMP)*Y(K)/BODY(NTOT)
          X(K,JCOMP) = X(K,ICOMP) - Y(K)
          XDOT(K,ICOMP) = XDOT(K,NTOT) + BODY(JCOMP)*V(K)/BODY(NTOT)
          XDOT(K,JCOMP) = XDOT(K,ICOMP) - V(K)
   30 CONTINUE
*
*       Predict current binding energy.
      ZMU = BODY(ICOMP)*BODY(JCOMP)/BODY(NTOT)
      PP = (((PDOT4*DTU/24.0 + ONE6*PDOT3)*DTU + 0.5*PDOT2)*DTU +
     &                                           PDOT)*DTU + P
      EBIN = 0.5*ZMU*PP
*
      RIJ2 = 0.0
      VIJ2 = 0.0
      DO 34 K = 1,3
      RIJ2 = RIJ2 + (X(K,1) - X(K,2))**2
      VIJ2 = VIJ2 + (XDOT(K,1) - XDOT(K,2))**2
   34 CONTINUE
      SEMI = 2.0/SQRT(RIJ2) - VIJ2/(BODY(ICOMP)+BODY(JCOMP))
      SEMI = 1.0/SEMI
      EB = -BODY(1)*BODY(2)/(2.0*SEMI)
      ERR = (EB - EBIN)/BE(3)
*     WRITE (6,36)  SQRT(RIJ2), SEMI, (EB - EBIN)/EBIN, ERR
*  36 FORMAT (' COMPARE   RIJ A DB DE ',1P,2E13.4,2E10.2)
      RETURN
*
      END
