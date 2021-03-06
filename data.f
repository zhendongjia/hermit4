
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
      READ (5,*) G_P, G_D, G_R, T_DEP, R_EDGE, R_IN, DENS0,
     &             DENS_P, M_CRIT, R_ESC, KG, EJ, CD, R_MHS, THREE_D,
     &             NLEAST, TWOGG
      WRITE (6,5) G_P, G_D, G_R, T_DEP, R_EDGE, R_IN, DENS0,
     &             DENS_P, M_CRIT, R_ESC, KG, EJ, CD, R_MHS, THREE_D,
     &             NLEAST, TWOGG
 5    FORMAT (/, 2X, 1P, 10E10.1, F10.2, 3F10.4, 3I10)
      GCM2_MAU2 = 1.125D-7
      GCM3_MAU3 = 1.7D6
      DENS0 = DENS0*GCM2_MAU2
      DENS_P = DENS_P*GCM3_MAU3
      DO 2 I = 1, N
      RADIUS(I) = (BODY(I)*3/(2*TWOPI*DENS_P))**0.333
 2    CONTINUE
*
      IF (TWOGG.GT.0) THEN 
         READ(5,*) M_S, M_S_INT, T_S, R_S, R_EDGE_INT
      END IF
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
