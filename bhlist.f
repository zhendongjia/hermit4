      SUBROUTINE BHLIST
*
*
*       Perturber selection.
*       --------------------
*
      INCLUDE 'commonp.h'
*
*
*       Set square perturber distance.
      R2 = R(1)**2 + R(2)**2 + R(3)**2
      RCRIT2 = R2*CMSEP2
*
*       Select new perturbers.
      NNB1 = 1
      DO 10 J = IFIRST,NMASS
          W1 = X(1,J) - X(1,NTOT)
          W2 = X(2,J) - X(2,NTOT)
          W3 = X(3,J) - X(3,NTOT)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
          IF (RSEP2.LT.RCRIT2) THEN
              NNB1 = NNB1 + 1
              LIST(NNB1) = J
          END IF
   10 CONTINUE
*
*       Save perturber membership.
      LIST(1) = NNB1 - 1
*
      RETURN
*
      END
