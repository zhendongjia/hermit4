      SUBROUTINE GET_OUTER_GRAVITY(R12,R12_DOT,ROUT,ROUT_DOT,
     &                              DENS,DENS_DOT,FOUT,FOUT_DOT)
      INCLUDE 'commonp.h'
      REAL*8 R12,R12_DOT,ROUT,ROUT_DOT,DENS,DENS_DOT
      REAL*8 FOUT,FOUT_DOT
C      
      FOUT = TWOPI*DENS*((0.5/(1+KG))*(R12/ROUT)**(1+KG)
     &          + (144.0/(3+KG))*(R12/ROUT)**(3+KG)) 
      FOUT_DOT = FOUT*DENS_DOT/DENS
     &      + TWOPI*DENS*(0.5*(R12/ROUT)**(KG)
     &            + 144.0*(R12/ROUT)**(2+KG))
     &       *(R12_DOT/ROUT - R12*ROUT_DOT/ROUT**2)
C
      END
