      SUBROUTINE LOAD_DATA
      INCLUDE 'commonp.h'

      OPEN(10, FILE='../ap_f_fdot', status='old')
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


      OPEN(20, FILE='../aj_f_fdot', status='old')
      AJ_LINE=1
      DO WHILE(.TRUE.)
         READ(20, *, end=200) AJ_VAR(AJ_LINE), AJ_F(AJ_LINE),
     &        AJ_FDOT(AJ_LINE)
         AJ_LINE = AJ_LINE + 1
      END DO
 200  AJ_LINE = AJ_LINE - 1
      AJ_LOWER_BOUND = AJ_VAR(1)
      AJ_UPPER_BOUND = AJ_VAR(AJ_LINE)
      AJ_STEP = (AJ_UPPER_BOUND - AJ_LOWER_BOUND) / (AJ_LINE - 1)

      IF (N.GT.2) THEN
      OPEN(30, FILE='../as_f_fdot', status='old')
      AS_LINE=1
      DO WHILE(.TRUE.)
         READ(30, *, end=300) AS_VAR(AS_LINE), AS_F(AS_LINE),
     &        AS_FDOT(AS_LINE)
         AS_LINE = AS_LINE + 1
      END DO
 300  AS_LINE = AS_LINE - 1
      AS_LOWER_BOUND = AS_VAR(1)
      AS_UPPER_BOUND = AS_VAR(AS_LINE)
      AS_STEP = (AS_UPPER_BOUND - AS_LOWER_BOUND) / (AS_LINE - 1)
      END IF


      END
