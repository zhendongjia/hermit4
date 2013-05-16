      REAL*8 FUNCTION GET_PRECESSION(X,XDOT,SEMI,ECC,W_F,TA_F)
      IMPLICIT NONE

      REAL*8 X(3),XDOT(3),SEMI,ECC
      REAL*8 R12,H,RDOT12,SIN_F,COS_F,SIN_W_F,COS_W_F,COS_W,SIN_W,W
      REAL*8 W_F,TA_F
      R12 = SQRT(X(1)**2 + X(2)**2)
      H = X(1)*XDOT(2) - X(2)*XDOT(1)
      RDOT12 = (X(1)*XDOT(1) + X(2)*XDOT(2))/R12

      SIN_F = SEMI*(1-ECC**2)*RDOT12/(H*ECC)
      COS_F = SEMI*(1-ECC**2)/(R12*ECC) - 1/ECC

      SIN_W_F = X(2)/R12
      COS_W_F = X(1)/R12 
      
*     W_F = ASIN(SIN_W_F)
*     TA_F = ASIN(SIN_F)

      COS_W = COS_W_F*COS_F + SIN_W_F*SIN_F
      SIN_W = SIN_W_F*COS_F - COS_W_F*SIN_F
      
      W_F = ATAN2(SIN_W_F,COS_W_F)
      W_F = W_F*180/ACOS(-1.)
      TA_F = ATAN2(SIN_F,COS_F)
      TA_F = TA_F*180/ACOS(-1.)
      W = ATAN2(SIN_W,COS_W)
      W = W*180/ACOS(-1.)
*     TA_F = TA_F*180/ACOS(-1.)
*     W_F = W_F*180/ACOS(-1.)
*     W = W_F - TA_F

      GET_PRECESSION = W

      END
