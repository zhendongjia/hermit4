      SUBROUTINE BHINIT
*
*
*       Initialization of Burdet-Heggie regularization.
*       -----------------------------------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XI(6),VI(6),FP(6),FD(6)
*
*
*       Introduce c.m. name and mass.
      ICOMP = 1
      JCOMP = 2
      NAME(NTOT) = NZERO + NAME(ICOMP)
      BODY(NTOT) = BODY(ICOMP) + BODY(JCOMP)
*
*       Define c.m. coordinates & velocities and set XDOT for components.
      DO 10 K = 1,3
          X(K,NTOT) = (BODY(ICOMP)*X(K,ICOMP) + BODY(JCOMP)*X(K,JCOMP))/
     &                                                        BODY(NTOT)
          X0DOT(K,NTOT) = (BODY(ICOMP)*X0DOT(K,ICOMP) + BODY(JCOMP)*
     &                                        X0DOT(K,JCOMP))/BODY(NTOT)
          XDOT(K,NTOT) = X0DOT(K,NTOT)
          XDOT(K,ICOMP) = X0DOT(K,ICOMP)
          XDOT(K,JCOMP) = X0DOT(K,JCOMP)
   10 CONTINUE
*
*       Obtain new force polynomials and time-steps for c.m. body.
*     CALL XVPRED(IFIRST,N)
      CALL FPOLY1(NTOT,NTOT)
      CALL STEPS(NTOT,NTOT)
*
*       Define relative coordinates and velocities in physical units.
      R2 = 0.0
      V2 = 0.0
      RRD = 0.0
      DO 20 K = 1,3
          R(K) = X(K,ICOMP) - X(K,JCOMP)
          RDOT(K) = XDOT(K,ICOMP) - XDOT(K,JCOMP)
          R2 = R2 + R(K)**2
          V2 = V2 + RDOT(K)**2
          RRD = RRD + R(K)*RDOT(K)
          XI(K) = X(K,ICOMP)
          XI(K+3) = X(K,JCOMP)
          VI(K) = XDOT(K,ICOMP)
          VI(K+3) = XDOT(K,JCOMP)
   20 CONTINUE
*
*       Introduce the Burdet-Heggie variables and convert dR/dt to R'.
      RR = SQRT(R2)
      RR0 = RR
      P = V2 - 2.0*BODY(NTOT)/RR
      DO 30 K = 1,3
          B(K) = BODY(NTOT)*R(K)/RR - V2*R(K) + RRD*RDOT(K)
          RDOT(K) = RR*RDOT(K)
   30 CONTINUE
*
*       Form perturber list (LIST(1) holds membership).
      CALL BHLIST
*
*       Initialize the BH elements and derivative.
      CALL BHPERT(XI,VI,FP,FD)
*
*       Absorb factorials in first and second derivatives.
      DO 40 K = 1,3
          RDOT2(K) = 0.5*RDOT2(K)
          RDOT3(K) = ONE6*RDOT3(K)
   40 CONTINUE
*
*       Assign a conservative step and convert to physical units.
      A2 = MIN(RMIN/BODY(NTOT),1.0/ABS(P))
      DTAU = ETAU*SQRT(A2)/(1.0 + 1000.0*GAMMA)**0.333
      STEP(1) = (((TDOT4/24.0D0*DTAU + ONE6*TDOT3)*DTAU +
     &                       0.5*TDOT2)*DTAU + RR)*DTAU
      T0(1) = TIME
      PDOT3 = 0.0
      PDOT4 = 0.0
*
*       Evaluate eccentricity and pericentre distance.
      RD2 = 0.0
      DO 45 K = 1,3
          RD2 = RD2 + RDOT(K)**2
   45 CONTINUE
      RD2 = RD2/RR**2
      A0 = 2.0/RR - RD2/(BODY(ICOMP) + BODY(JCOMP))
      A0 = 1.0/A0
      ECC2 = (1.0 - RR/A0)**2 + RRD**2/(A0*BODY(NTOT))
      ECC = SQRT(ECC2)
      PERI = A0*(1.0 - ECC)
*
      IF (KZ(10).GE.1) THEN
          SEMI = -BODY(NTOT)/P
          WRITE (6,50)  TIME/TWOPI, NAME(ICOMP), NAME(JCOMP), LIST(1),
     &                  ECC, RR, PERI, GAMMA
   50     FORMAT (' NEW BHREG    YRS NAM NP E R PER G ',
     &                           F10.5,3I4,F8.4,1P,3E10.2)
          CALL FLUSH(6)
      END IF
*
      RETURN
*
      END
