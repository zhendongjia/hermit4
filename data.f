      SUBROUTINE DATA
*
*
*       Initial conditions.
*       -------------------
*
      INCLUDE 'commonp.h'
      REAL*4  RAN2
*
*
*       Read initial conditions.
      DO 1 I = 1,N
          READ (5,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
    1 CONTINUE

*       Read parameters about tidal force by gas disk
      READ (5,*) G_D, DENS_ORI, DENS_P, T_TIDAL1, T_TIDAL2
      WRITE (6,5) G_D, DENS_ORI, DENS_P, T_TIDAL1, T_TIDAL2
 5    FORMAT (/, 2X, 1P, 6E10.1)

*
*       Initialize the portable random number generator (range: 0 to 1).
      KDUM = -1
      RN1 = RAN2(KDUM)
*       Skip the first random numbers (IDUM1 specified at input).
      DO 10 K = 1,IDUM1
          RN1 = RAN2(KDUM)
   10 CONTINUE
*
*       Initialize centre of mass terms.
      DO 25 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   25 CONTINUE
*
      ZMASS = 0.0
      RMAX2 = 0.0
      DO 40 I = 1,NMASS
          ZMASS = ZMASS + BODY(I)
          RI2 = 0.0
          DO 35 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
              RI2 = RI2 + X(K,I)**2
   35     CONTINUE
          RMAX2 = MAX(RI2,RMAX2)
   40 CONTINUE
*
*       Define nominal crossing time (cf. initial step and zero velocity).
      TCR = TWOPI
*
*       Save random number sequence in COMMON for future use.
      IDUM1 = KDUM
*
*       Define maximum quantized time-step from RMAX2.
      RK2 = 1.0
      DO 50 K = 1,8
          IF (RK2.LT.RMAX2) THEN
              DTMAX = 2.0D0*DTMAX
              RK2 = 4.0*RK2
          END IF
   50 CONTINUE
*
      RETURN
*
      END
