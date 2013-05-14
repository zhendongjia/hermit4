      SUBROUTINE ZERO
*
*
*       Initialization of global scalars. 
*       ---------------------------------
*
      INCLUDE 'commonp.h'
*
*
*       Initialize parameters & counters and set useful constants.
      TIME = 0.0D0
      TPRINT = 0.0D0
      CPUTOT = 0.0
      ERRTOT = 0.0D0
      DETOT = 0.0D0
      DTMAX = 1.0
      MODEL = 0
      NSTEPS = 0
      NTIMER = 0
      NKSTRY = 0
      NSTEPU = 0
      NPAIRS = 0
      IPHASE = 0
      IFIRST = 1
      LIST(1) = 0
*
*       Set fractional constants & two PI.
      ONE3 = 1.0D0/3.0D0
      ONE6 = 1.0D0/6.0D0
      ONE9 = 1.0D0/9.0D0
      ONE12 = 1.0D0/12.0D0
      TWOPI = 8.0D0*ATAN(1.0D0)
*
      RETURN
*
      END
