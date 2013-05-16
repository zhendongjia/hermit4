      SUBROUTINE GAS_POTENTIAL(XI,XIDOT,FIN,FDN,DENS0)
      

      INCLUDE 'commonp.h'
      REAL*8 XI(3),XIDOT(3),FIN(3),FDN(3),DENS0
      REAL*8 R12,R_DOT12,COS_THE,SIN_THE,
     &     COS_THE_DOT,SIN_THE_DOT,DENS,DENS_DOT,FIN_ABS,
     &     FIN_ABS_DOT
      

      R12 = SQRT(XI(1)**2 + XI(2)**2)
      R_DOT12 = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2))/R12

      COS_THE = XI(1)/R12
      COS_THE_DOT = XIDOT(1)/R12 - XI(1)*R_DOT12/R12**2

      SIN_THE = XI(2)/R12
      SIN_THE_DOT = XIDOT(2)/R12 - XI(2)*R_DOT12/R12**2
         
      DENS = DENS0*EXP(-TIME/T_DEP)*R12**(-1.5)
      DENS_DOT = DENS0*(-1/T_DEP)*EXP(-TIME/T_DEP)*R12**(-1.5) 
     &           + DENS0*(-1.5)*R12**(-2.5)*R_DOT12*EXP(-TIME/T_DEP)
      
      FIN_ABS = 2*3.14*DENS*(0.2*(R12/R_EDGE)**(2.5)
     &          + 32*(R12/R_EDGE)**(4.5)) 


      FIN(1) = FIN(1) + FIN_ABS*COS_THE
      FIN(2) = FIN(2) + FIN_ABS*SIN_THE



      FIN_ABS_DOT = 2*3.14*DENS_DOT*(0.2*(R12/R_EDGE)**(2.5)
     &                              + 32*(R12/R_EDGE)**(4.5))
     &              +2*3.14*DENS*(0.5*(R12/R_EDGE)**(1.5)
     &                              + 144*(R12/R_EDGE)**(3.5))
     &              *(R_DOT12/R_EDGE)

      FDN(1) = FDN(1) + FIN_ABS_DOT*COS_THE + FIN_ABS*COS_THE_DOT
      FDN(2) = FDN(2) + FIN_ABS_DOT*SIN_THE + FIN_ABS*SIN_THE_DOT


      RETURN


      END
