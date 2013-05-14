      SUBROUTINE CMF(XI,XIDOT,FIN,FDN)
*
*
*       Force on c.m. particle.
*       -----------------------
*
      INCLUDE 'commonp.h'
      REAL*8  XI(3),XIDOT(3),DV(3),FIN(3),FDN(3),FF(6),FD(6)
*
*
*       Obtain the coordinates and velocities.
      RR = SQRT(R(1)**2 + R(2)**2 + R(3)**2)
      DO 1 K = 1,3
*       Note small phase error since R & RDOT known just before TBLOCK.
          X(K,ICOMP) = XI(K) + BODY(JCOMP)*R(K)/BODY(NTOT)
          X(K,JCOMP) = XI(K) - BODY(ICOMP)*R(K)/BODY(NTOT)
          XDOT(K,ICOMP) = XIDOT(K) + BODY(JCOMP)*RDOT(K)/(RR*BODY(NTOT))
          XDOT(K,JCOMP) = XIDOT(K) - BODY(ICOMP)*RDOT(K)/(RR*BODY(NTOT))
    1 CONTINUE
*
*       Initialize individual force components.
      DO 8 K = 1,6
          FF(K) = 0.0
          FD(K) = 0.0
    8 CONTINUE
*
      DO 10 J = IFIRST,NMASS
*       Sum over both components (makes hardly any difference).
      DO 5 I = 1,2
          K = 3*(I - 1)
          A1 = X(1,J) - X(1,I)
          A2 = X(2,J) - X(2,I)
          A3 = X(3,J) - X(3,I)
          DV(1) = XDOT(1,J) - XDOT(1,I)
          DV(2) = XDOT(2,J) - XDOT(2,I)
          DV(3) = XDOT(3,J) - XDOT(3,I)
          RIJ2 = A1*A1 + A2*A2 + A3*A3
*
          DR2I = 1.0/RIJ2
          DR3I = BODY(J)*DR2I*SQRT(DR2I)
          DRDV = 3.0*(A1*DV(1) + A2*DV(2) + A3*DV(3))*DR2I
*
          FF(1+K) = FF(1+K) + A1*DR3I
          FF(2+K) = FF(2+K) + A2*DR3I
          FF(3+K) = FF(3+K) + A3*DR3I
          FD(1+K) = FD(1+K) + (DV(1) - A1*DRDV)*DR3I
          FD(2+K) = FD(2+K) + (DV(2) - A2*DRDV)*DR3I
          FD(3+K) = FD(3+K) + (DV(3) - A3*DRDV)*DR3I
    5 CONTINUE
   10 CONTINUE
*
*       Form the mass-weighted c.m. force.
      DO 20 K = 1,3
          FIN(K) = (BODY(ICOMP)*FF(K) + BODY(JCOMP)*FF(K+3))/BODY(NTOT)
          FDN(K) = (BODY(ICOMP)*FD(K) + BODY(JCOMP)*FD(K+3))/BODY(NTOT)
   20 CONTINUE
*
*       Add optional indirect terms (only single summation required).
      IF (KZ(3).GT.0) THEN
      DO 30 J = IFIRST,NMASS
          XJ = X(1,J)
          YJ = X(2,J)
          ZJ = X(3,J)
          RJ2 = XJ**2 + YJ**2 + ZJ**2
          FJ = BODY(J)/(RJ2*SQRT (RJ2))
          FIN(1) = FIN(1) - XJ*FJ
          FIN(2) = FIN(2) - YJ*FJ
          FIN(3) = FIN(3) - ZJ*FJ
          RRD = 3.0*(XJ*XDOT(1,J) + YJ*XDOT(2,J) + ZJ*XDOT(3,J))/RJ2
          FDN(1) = FDN(1) - (XDOT(1,J) - RRD*XJ)*FJ
          FDN(2) = FDN(2) - (XDOT(2,J) - RRD*YJ)*FJ
          FDN(3) = FDN(3) - (XDOT(3,J) - RRD*ZJ)*FJ
   30 CONTINUE
      END IF
*
      RETURN
*
      END
