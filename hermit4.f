*
*             H E R M I T 4
*             *************
*
*       N-body code with Hermite block-steps, iteration & regularization.
*       -----------------------------------------------------------------
*
*       Method of Kokubo, Yoshinaga & Makino M.N. 297, 1067, 1998
*       Burdet - Heggie two-body regularization included 2/04
*       ---------------------------------------------------------
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      PROGRAM HERMIT4
*
      INCLUDE 'commonp.h'
*
*
*       Initialize the timer.
      CALL CPUTIM(CPU0)
*
      CALL LOAD_DATA
*       Read start/restart indicator & CPU time.
      READ (5,*)  KSTART, TCOMP
*
      IF (KSTART.EQ.1) THEN
*
*       Read input parameters, perform initial setup and obtain output.
          CPU = TCOMP
          CALL START
          CALL OUTPUT
      ELSE
*
*       Read previously saved COMMON variables from tape/disc on unit 1.
          CALL MYDUMP(0,1)
*       Safety indicator preventing repeated restarts set in routine CHECK.
          CPU = TCOMP
          CPU0 = 0.0
          IF (KSTART.EQ.3) READ (5,*) J,K
          IF (J.GT.0) KZ(J) = K
      END IF
*
*       Advance solutions until next output or change of procedure.
    1 CALL INTGRT
*
      IF (IPHASE.EQ.1) THEN
          CALL BHREG
      ELSE IF (IPHASE.EQ.2) THEN
          CALL BHTERM
      ELSE IF (IPHASE.EQ.3) THEN
          CALL OUTPUT
      END IF
*
*       Continue integration.
      GO TO 1
*
      END
