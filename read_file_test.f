      IMPLICIT NONE

      INTEGER MAXLINE
      REAL*8 AP_VAR, AP_F, AP_FDOT
      REAL*8 AP_LOWER_BOUND, AP_UPPER_BOUND, AP_STEP
      INTEGER AP_LINE

      PARAMETER (MAXLINE=32768)
      COMMON/AP_F_FDOT/ AP_VAR(MAXLINE), AP_F(MAXLINE),
     &     AP_FDOT(MAXLINE), AP_LOWER_BOUND, AP_UPPER_BOUND,
     &     AP_STEP, AP_LINE

      REAL*8 X, F, FDOT

      CALL READ_AP_F_FDOT()
      WRITE(*,*) AP_LOWER_BOUND, AP_UPPER_BOUND
      DO WHILE(.TRUE.)
         READ (*, *) X
         if (X.LT.AP_LOWER_BOUND .OR. X.GT.AP_UPPER_BOUND) EXIT
         CALL GET_AP_F_FDOT(X, F, FDOT)
         WRITE(*, *) X, F, FDOT
      END DO
      END

      SUBROUTINE GET_AP_F_FDOT(X, F, FDOT)
      INTEGER MAXLINE
      REAL*8 AP_VAR, AP_F, AP_FDOT
      REAL*8 AP_LOWER_BOUND, AP_UPPER_BOUND, AP_STEP
      INTEGER AP_LINE

      PARAMETER (MAXLINE=65536)
      COMMON/AP_F_FDOT/ AP_VAR(MAXLINE), AP_F(MAXLINE),
     &     AP_FDOT(MAXLINE), AP_LOWER_BOUND, AP_UPPER_BOUND,
     &     AP_STEP, AP_LINE
      REAL*8 X, F, FDOT, P0, P1
      INTEGER IDX
      IF (X.LE.AP_LOWER_BOUND) THEN
         F = AP_F(1)
         FDOT = AP_FDOT(1)
         RETURN
      END IF
      IF (X.GE.AP_UPPER_BOUND) THEN
         F = AP_F(AP_LINE)
         FDOT = AP_FDOT(AP_LINE)
         RETURN
      END IF
      IDX = (X - AP_LOWER_BOUND) / AP_STEP + 1
      P0 = AP_VAR(IDX+1) - X
      P1 = X - AP_VAR(IDX)
      F = (AP_F(IDX) * P0 + AP_F(IDX+1) * P1) / (P0 + P1)
      FDOT = (AP_FDOT(IDX) * P0 + AP_FDOT(IDX+1) * P1) / (P0 + P1)
      RETURN
      END


      SUBROUTINE READ_AP_F_FDOT()
      INTEGER MAXLINE
      REAL*8 AP_VAR, AP_F, AP_FDOT
      REAL*8 AP_LOWER_BOUND, AP_UPPER_BOUND, AP_STEP
      INTEGER AP_LINE

      PARAMETER (MAXLINE=65536)
      COMMON/AP_F_FDOT/ AP_VAR(MAXLINE), AP_F(MAXLINE),
     &     AP_FDOT(MAXLINE), AP_LOWER_BOUND, AP_UPPER_BOUND,
     &     AP_STEP, AP_LINE

      OPEN(10, FILE='ap_f_fdot', status='old')

      AP_LINE=1
      DO WHILE(.TRUE.)
         READ(10, *, end=100) AP_VAR(AP_LINE), AP_F(AP_LINE),
     &        AP_FDOT(AP_LINE)
         AP_LINE = AP_LINE + 1
      END DO
 100  AP_LINE = AP_LINE - 1
      AP_LOWER_BOUND = AP_VAR(1)
      AP_UPPER_BOUND = AP_VAR(AP_LINE)
      AP_STEP = (AP_UPPER_BOUND - AP_LOWER_BOUND) / (AP_LINE - 1)
      RETURN
      END
-
