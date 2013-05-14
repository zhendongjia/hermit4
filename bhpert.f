      SUBROUTINE BHPERT(XI,VI,FP,FD)
*
*
*       Physical perturbation in Burdet-Heggie regularization.
*       ------------------------------------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XI(6),VI(6),FP(6),FD(6)
*
*
*       Initialize individual force components.
      DO 5 K = 1,6
          FP(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE
*
*       Obtain planetary perturbation on each component.
      NP1 = LIST(1) + 1
      DO 10 L = 2,NP1
          K = LIST(L)
          A1 = X(1,K) - XI(1)
          A2 = X(2,K) - XI(2)
          A3 = X(3,K) - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          V1 = XDOT(1,K) - VI(1)
          V2 = XDOT(2,K) - VI(2)
          V3 = XDOT(3,K) - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*
          A1 = X(1,K) - XI(4)
          A2 = X(2,K) - XI(5)
          A3 = X(3,K) - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = XDOT(1,K) - VI(4)
          V2 = XDOT(2,K) - VI(5)
          V3 = XDOT(3,K) - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
   10 CONTINUE
*
*       Include solar perturbation (note: indirect terms cancel).
          SUN1 = 1.0 + BODY(1)
          SUN2 = 1.0 + BODY(2)
          A1 = - XI(1)
          A2 = - XI(2)
          A3 = - XI(3)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = SUN1/(RIJ2*SQRT(RIJ2))
          FP(1) = FP(1) + A1*A6
          FP(2) = FP(2) + A2*A6
          FP(3) = FP(3) + A3*A6
          V1 = - VI(1)
          V2 = - VI(2)
          V3 = - VI(3)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(1) = FD(1) + (V1 - A1*A9)*A6
          FD(2) = FD(2) + (V2 - A2*A9)*A6
          FD(3) = FD(3) + (V3 - A3*A9)*A6
*
          A1 = - XI(4)
          A2 = - XI(5)
          A3 = - XI(6)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
          A6 = SUN2/(RIJ2*SQRT(RIJ2))
          FP(4) = FP(4) + A1*A6
          FP(5) = FP(5) + A2*A6
          FP(6) = FP(6) + A3*A6
          V1 = - VI(4)
          V2 = - VI(5)
          V3 = - VI(6)
          A9 = 3.0D0*(A1*V1 + A2*V2 + A3*V3)/RIJ2
          FD(4) = FD(4) + (V1 - A1*A9)*A6
          FD(5) = FD(5) + (V2 - A2*A9)*A6
          FD(6) = FD(6) + (V3 - A3*A9)*A6
*
*       Form relative perturbation and physical derivative.
      FP2 = 0.0
      DO 20 K = 1,3
          FP(K) = FP(K) - FP(K+3)
          FD(K) = FD(K) - FD(K+3)
          FP2 = FP2 + FP(K)**2
   20 CONTINUE
*
*       Obtain current two-body separation.
      R2 = R(1)**2 + R(2)**2 + R(3)**2
      RR = SQRT(R2)
*
*       Evaluate useful scalars.
      PDOT = 0.0
      PDOT2 = 0.0
      RPR = 0.0
      RF = 0.0
      RDF = 0.0
      RFPR = 0.0
      RPR2 = 0.0
      RR2D = 0.0
*       Note that RDOT(K) and PDOT are regularized derivatives.
      DO 25 K = 1,3
          PDOT = PDOT + 2.0*RDOT(K)*FP(K)
          RPR = RPR + R(K)*RDOT(K)
          RF = RF + R(K)*FP(K)
          RDF = RDF + RDOT(K)*FP(K)
          RFPR = RFPR + R(K)*FD(K)
          RPR2 = RPR2 + RDOT(K)**2
   25 CONTINUE
      RPR = RPR/RR
      RFPR = RR*RFPR
*       Use RR*RPR in BDOT for scalar product R*RDOT.
      TDOT2 = RPR
      TDOT3 = P*RR + BODY(NTOT) + RF*RR
      TDOT4 = PDOT*RR + P*RPR + RF*RPR
*
*       Construct the derivatives of R and BH elements (F'=R*dF/dt).
      DO 30 K = 1,3
          BDOT(K) = -2.0*RDF*R(K) + RF*RDOT(K) + RR*RPR*FP(K)
          RDOT2(K) = P*R(K) + B(K) + R2*FP(K)
          RDOT3(K) = PDOT*R(K) + P*RDOT(K) + BDOT(K) +
     &                           (2.0*RPR*FP(K) + R2*FD(K))*RR
          PDOT2 = PDOT2 + 2.0*(RDOT2(K)*FP(K) + RDOT(K)*RR*FD(K))
          TDOT4 = TDOT4 + (RDOT(K)*FP(K) + R(K)*RR*FD(K))*RR
          RR2D = RR2D + R(K)*RDOT2(K)
   30 CONTINUE
*
*       Form second derivative of B by boot-strapping.
      DO 35 K = 1,3
          BDOT2(K) = -2.0*(RDOT2(K)*FP(K) + RDOT(K)*RR*FD(K))*R(K) -
     &                (RDF - RFPR)*RDOT(K) + RF*RDOT2(K) +
     &                (RPR2 + RR2D)*FP(K) + RPR*RR**2*FD(K)
   35 CONTINUE
*
*       Specify relative perturbation.
      GAMMA = SQRT(FP2)*R2/BODY(NTOT)
*
      RETURN
*
      END
