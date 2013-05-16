      SUBROUTINE MYDUMP(I,J)
*
*
*       COMMON save or read.
*       --------------------
*
      PARAMETER  (NMAX=1000,NA=65*NMAX+26,NP=2*NMAX+84)
      REAL*4  A,E,F,P
      INTEGER  IR,NAME
*
*
      COMMON/NBODY/  A(NA)
      COMMON/PARAMS/ E(64)
      COMMON/BHREG1/ F(104)
      COMMON/BHNAME/ NAME(209)
      COMMON/RAND2/  IR(99)
      COMMON/BLOCKS/ P(NP)
*
*
*       Open unit #J by reading dummy and rewinding.
      REWIND J
      READ (J,ERR=10,END=10)  DUMMY
   10 REWIND J
*
*       Read or save all COMMON variables (valid for tape or disc).
      IF (I.EQ.0) THEN
          READ (J)   A, E, F, NAME, IR, P
      ELSE
          WRITE (J)  A, E, F, NAME, IR, P
          END FILE J
          CLOSE (UNIT=J)
      END IF
*
      RETURN
*
      END
