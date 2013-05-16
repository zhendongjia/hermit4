
      SUBROUTINE DAMPING(XI,XIDOT,FIN,FDN,I)
      
      INCLUDE 'commonp.h'
      INTEGER I
      REAL*8 XI(3),XIDOT(3),FIN(3),FDN(3)
      REAL*8 R12,R12_DOT,V,V_DOT,COSTHE,SINTHE,VC(3),VC_DOT(3),
     &       DELT_V(3),DELT_V_DOT(3),F_TIDAL(3),F_TIDAL_DOT(3),
     &       RP,M_SUN,YR_S,AU_CM,H,DENS_S,DENS_G,H_DOT,DENS_S_DOT,
     &       DENS_G_DOT,V_E,V_AU_YR,T_TIDAL1_DOT,
     &       T_TIDAL2_DOT,V_AU_YR_DOT,V_T_TIDAL_DOT,
     &       ABS_DELT_V, ABS_DELT_V_DOT

      
      R12 = SQRT(XI(1)**2+XI(2)**2)
      V = SQRT(XIDOT(1)**2+XIDOT(2)**2)
      R12_DOT = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2))/R12
    

      COSTHE = XI(1)/R12
      SINTHE = XI(2)/R12

      VC(1) = -R12**(-0.5)*SINTHE
      VC(2) = R12**(-0.5)*COSTHE
      VC(3) = 0

      VC_DOT(1) = 1.5*R12**(-1.5)*R12_DOT*SINTHE - R12**(-1.5)*XIDOT(2)
      VC_DOT(2) = R12**(-1.5)*XIDOT(1) - 1.5*R12**(-1.5)*R12_DOT*COSTHE
      VC_DOT(3) = 0

      DELT_V(1) = XIDOT(1) - VC(1)
      DELT_V(2) = XIDOT(2) - VC(2)
      DELT_V(3) = XIDOT(3) - VC(3)

      ABS_DELT_V = SQRT(DELT_V(1)**2+DELT_V(2)**2+DELT_V(3)**2)
            

*in g
      M_SUN = 1.989E33     
*get planet radius in cm unit
      RP = 4*3.14*DENS_P/3
      RP = BODY(I)*M_SUN/RP
      RP = RP**(0.333)
*unit change from s to yr
      YR_S = 3E7
*speed of sound in cm/yr
*      V_S = 30000*YR_S
*unit change from au to cm
      AU_CM = 1.5E13
*disk scale height in au
      H = 0.05*R12**1.25
*gas disk surface density in g/cm^2
      DENS_S = DENS_ORI*EXP(-TIME/T_DEP)
*gas disk density in g/cm^3
      DENS_G = DENS_S/(H*AU_CM)

      H_DOT = H*1.25*R12_DOT/R12
      DENS_S_DOT = DENS_S*(-1/T_DEP)
      DENS_G_DOT = DENS_G*DENS_S_DOT/DENS_S - DENS_G*H_DOT/H 


*earth velocity in km/s
      V_E = 30
*change planet velocity from nbody unit to au/year
      V_AU_YR = V*V_E*1E5*YR_S/AU_CM

      T_TIDAL1(I)=(DENS_P/DENS_G)*(RP/AU_CM)*(8/3)/(ABS_DELT_V*V_AU_YR)
      
      T_TIDAL2(I) = 1/BODY(I)
      T_TIDAL2(I) = T_TIDAL2(I)*M_SUN/(DENS_S*(R12*AU_CM)**2)
      T_TIDAL2(I) = T_TIDAL2(I)*(H/R12)**4
      T_TIDAL2(I) = T_TIDAL2(I)*R12/V_AU_YR

*      WRITE (0,*) T_TIDAL1, T_TIDAL2

      F_TIDAL(1) = -DELT_V(1)*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
      F_TIDAL(2) = -DELT_V(2)*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))
      F_TIDAL(3) = -DELT_V(3)*(1/T_TIDAL1(I) + 1/T_TIDAL2(I))

      FIN(1) = FIN(1) + F_TIDAL(1)
      FIN(2) = FIN(2) + F_TIDAL(2)
      FIN(3) = FIN(3) + F_TIDAL(3)

      V_DOT = (XIDOT(1)*FIN(1) + XIDOT(2)*FIN(2))/V
      V_AU_YR_DOT = V_AU_YR*V_DOT/V

      DELT_V_DOT(1) = FIN(1) - VC_DOT(1)
      DELT_V_DOT(2) = FIN(2) - VC_DOT(2)
      
      ABS_DELT_V_DOT = DELT_V(1)*DELT_V_DOT(1)+DELT_V(2)*DELT_V_DOT(2)  

      T_TIDAL1_DOT = T_TIDAL1(I)*(-DENS_G_DOT)/DENS_G
     &               + T_TIDAL1(I)*(-ABS_DELT_V_DOT)/ABS_DELT_V   
      T_TIDAL2_DOT = T_TIDAL2(I)*(-DENS_S_DOT/DENS_S - 2*R12_DOT/R12)  
     &               + T_TIDAL2(I)*(4*H_DOT/H -4*R12_DOT/R12)
     &               + T_TIDAL2(I)*(R12_DOT/R12 -V_AU_YR_DOT/V_AU_YR)

      V_T_TIDAL_DOT = T_TIDAL1_DOT/T_TIDAL1(I)**2  
     &      + T_TIDAL2_DOT/T_TIDAL2(I)**2
      
      F_TIDAL_DOT(1) = F_TIDAL(1)*DELT_V_DOT(1)/DELT_V(1)
     &                 + DELT_V(1)*V_T_TIDAL_DOT
      F_TIDAL_DOT(2) = F_TIDAL(2)*DELT_V_DOT(2)/DELT_V(2)
     &                 + DELT_V(2)*V_T_TIDAL_DOT

      FDN(1) = FDN(1) + F_TIDAL_DOT(1)
      FDN(2) = FDN(2) + F_TIDAL_DOT(2)
      

      RETURN      
      END


