      SUBROUTINE GAS_POTENTIAL_3D(XI,XIDOT,FIRR,FD)
      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIRR(3),FD(3), T1, T2
C
      YEAR = TIME / TWOPI
      DENS_EDGE = DENS0 * EXP(-YEAR/T_DEP) * R_EDGE**(-KG)
      DENS_EDGE_DOT = DENS_EDGE * (-1/T_DEP) * (1/TWOPI)
C
      RP = SQRT(XI(1)**2 + XI(2)**2 + XI(3)**2)
      RP_DOT = (XI(1) * XIDOT(1) + XI(2) * XIDOT(2) 
     &                           + XI(3) * XIDOT(3)) / RP
C
      T1 = 0.5/(1+KG)
      T2 = 144.0/(3+KG)
      CONST = TWOPI * DENS_EDGE * (T1 * (RP/R_EDGE) 
     &          + T2 * (RP/R_EDGE)**3) 
C
C
      CONST_DOT = CONST * DENS_EDGE_DOT / DENS_EDGE
     &         + TWOPI * DENS_EDGE * (T1 * (1/R_EDGE)
     &          + T2 * (1/R_EDGE)**3 * 3 * RP**2 * RP_DOT) 
C
      FIRR(1) = FIRR(1) + CONST * XI(1) / RP
      FIRR(2) = FIRR(2) + CONST * XI(2) / RP
      FIRR(3) = FIRR(3) - 2 * CONST * XI(3) / RP
C
      FD(1) = FD(1) + CONST_DOT * XI(1) / RP
     &              + CONST * XIDOT(1) / RP
     &              - CONST * XI(1) * RP_DOT / RP**2
      FD(2) = FD(2) + CONST_DOT * XI(2) / RP
     &              + CONST * XIDOT(2) / RP
     &              - CONST * XI(2) * RP_DOT / RP**2      
      FD(3) = FD(3) - 2 * CONST_DOT * XI(3) / RP
     &              - 2 * CONST * XIDOT(3) / RP
     &     + 2 * CONST * XI(3) * RP_DOT / RP**2
C   
      END
