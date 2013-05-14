      SUBROUTINE BHTERM
*
*
*       Termination of Burdet-Heggie regularization.
*       --------------------------------------------
*
      INCLUDE 'commonp.h'
*
*
*       Obtain current global coordinates and velocities.
      CALL RESOLV
*
*       Check optional diagnostics.
      IF (KZ(10).GE.2) THEN
          SEMI = -BODY(NTOT)/P
          RR = SQRT(R(1)**2 + R(2)**2 + R(3)**2)
          WRITE (6,10)  TIME/TWOPI, LIST(1), RR, SEMI, GAMMA
   10     FORMAT (' END BHREG    YRS NP R A G ',F10.5,I4,1P,3E10.2)
      END IF
*
*       Update decision-making variables.
      NPAIRS = NPAIRS - 1
      IFIRST = 2*NPAIRS + 1
      NTOT = N
      LIST(1) = 0
*
*       Initialize primary velocities (X0 done in routine STEPS).
      DO 20 K = 1,3
          X0DOT(K,ICOMP) = XDOT(K,ICOMP)
          X0DOT(K,JCOMP) = XDOT(K,JCOMP)
   20 CONTINUE
*
*       Initialize force polynomials and prescribe new time-steps.
      CALL FPOLY1(ICOMP,JCOMP)
      CALL STEPS(ICOMP,JCOMP)
*
*       Specify new total energy using current variables.
      CALL XVPRED(IFIRST,NTOT)
      CALL ENERGY
      BE(3) = ZKIN - POT
*
      RETURN
*
      END
