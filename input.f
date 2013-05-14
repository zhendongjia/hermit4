      SUBROUTINE INPUT
*
*
*       Parameter input.
*       ----------------
*
      INCLUDE 'commonp.h'
*
*
*       Options.
*       --------------------------------------------------
*       1   Common dump on unit #1 on 'touch STOP' or TCRIT.
*       2   Common dump on unit #2 at TCRIT.
*       3   Indirect force terms.
*       4   Output of individual bodies (=1) or binaries (>1).
*       5   Alternative time-step expression (not recommended).
*       6   Output of eccentricity and semi-major axis.
*       7   Inelastic collision (comps replaced by c.m.).
*       8   Output on fort.8 of TIME, SEMI, R for I = KZ(8).
*       9   Removal of escapers or ECC > ECRIT (see OUTPUT).
*      10   Diagnostic output of NEW BHREG and END BHREG.
*       --------------------------------------------------
*
*       Read & print the main input parameters.
      READ (5,*)  N, NRAND, ETA, DELTAT, TCRIT
      READ (5,*)  (KZ(J),J=1,10)
      NMASS = N
      NTOT = N
      NZERO = N
*
      WRITE (6,5)  N, ETA
    5 FORMAT (//,5X,'N =',I3,'  ETA =',F8.5)
      WRITE (6,10)  (KZ(J),J=1,10)
   10 FORMAT (/,5X,'OPTIONS:   ',10I4,/)
*
      READ (5,*)  DTMIN, RMIN, ETAU, GMIN
      WRITE (6,15)
   15 FORMAT (//,5X,'DTMIN     RMIN      ETAU      GMIN')
      WRITE (6,20)  DTMIN, RMIN, ETAU, GMIN
   20 FORMAT (/,2X,1P,4E10.1)
*
*       Convert time unit to 2*pi for planets (YRS = TIME/TWOPI).
      DELTAT = TWOPI*DELTAT
      TCRIT = TWOPI*TCRIT
*
*       Set random number skip for routine DATA.
      IDUM1 = NRAND
      RMIN2 = RMIN**2
      CMSEP2 = GMIN**(-0.666667)
*
      RETURN
*
      END
