      SUBROUTINE TRANSMEANS(                                            &
!                     DIMENSIONS
     &                  NSPEC ,  NANG , NFRE ,                          &
!                     INPUT
     &                  MFRCUT, FRCUT,                                  &
     &                  SPECW ,                                         &
!                     OUTPUT
     &                  MEANTE , MEANTANG , MEANTFRE)
!
!-------------------------------------------------------------------
!
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    , &
                 COSTH    ,SINTH     ,FRTAIL   ,WETAIL   ,FRATIO
      USE YOWPCONS , ONLY : EPSMIN

!-------------------------------------------------------------------

      IMPLICIT NONE
!
!     PURPOSE:
!     --------
!
!     MEAN VALUES OF TRANSFORMED SPECTRA
!     COMPUTES MEAN ENERGY,MEAN ANGLE AND MEAN FREQUENCY
!
!     INTERFACE:
!     ----------

      INTEGER(KIND=JWIM) :: NSPEC, NANG, NFRE
!                     DIMENSIONS.
      INTEGER(KIND=JWIM) :: MFRCUT
!                 INDEX OF THE LAST FREQUENCY FOR THE INTEGRATION
      REAL(KIND=JWRB) :: FRCUT
!                 CUT-OFF FREQUENCY 
      REAL(KIND=JWRB) :: SPECW(0:NSPEC,NANG,NFRE)
!                     2D SPECTRUM.
      REAL(KIND=JWRB) :: MEANTE(NSPEC), MEANTANG(NSPEC), MEANTFRE(NSPEC)
!                     MEAN PARAMETERS.
!     local variables:
!     ----------------

      INTEGER(KIND=JWIM) :: ISPEC, IANG, IFRE
      REAL(KIND=JWRB) :: FACTOR
      REAL(KIND=JWRB) :: SUM0(NSPEC), SUMX0(NSPEC), SUMY0(NSPEC)

!
!     ADDITIONAL (MF) VARIABLES:
!     ----------
!
!  JB: I don't see the point of doing this complicated calculation, instead we should be using a different DFIM
!      that is defined up to MFRCUT such that we perform the trapzoidalintegration method correctly !!!!
      INTEGER(KIND=JWIM) :: M10P
      REAL(KIND=JWRB) :: DELTH2
      REAL(KIND=JWRB) :: C10, D10, C11, C12, DF, SIH, COH, D11, CO, CO1
!
!-------------------------------------------------------------------
!

!     INITIALIZE:
!     -----------
      CO = FRATIO

!     COMPUTE MEAN ENERGIES:
!     ----------------------

      IF ( MFRCUT < NFRE ) THEN

        DELTH2 = DELTH / 2.0_JWRB

        DO ISPEC = 1,NSPEC
          MEANTE(ISPEC) = 0.0_JWRB
        END DO

        M10P = MIN(MFRCUT + 1_JWIM, NFRE)
        CO1 = CO - 1.0_JWRB
        C10  = (FRCUT - FR(MFRCUT))*DELTH - CO1*FR(MFRCUT)*DELTH2
        D10  = (FRCUT - FR(MFRCUT))**2*DELTH2/(CO1*FR(MFRCUT))
        DO ISPEC =1, NSPEC
          DO IANG = 1, NANG
            DO IFRE = 1, MFRCUT
              MEANTE(ISPEC) = MEANTE(ISPEC)+DFIM(IFRE)*SPECW(ISPEC,IANG,IFRE)
            ENDDO
            D11 = SPECW(ISPEC,IANG,M10P)-SPECW(ISPEC,IANG,MFRCUT)
            MEANTE(ISPEC) = MEANTE(ISPEC)+C10*SPECW(ISPEC,IANG,MFRCUT) + D10*D11
          ENDDO
        ENDDO

      ELSE
        DO ISPEC = 1,NSPEC
          MEANTE(ISPEC) = EPSMIN 
        END DO
        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC) = 0.0_JWRB
          END DO
          DO IANG = 1,NANG
            DO ISPEC = 1,NSPEC
              SUM0(ISPEC) = SUM0(ISPEC) + SPECW(ISPEC,IANG,IFRE)
            END DO
          END DO
          DO ISPEC = 1,NSPEC
            MEANTE(ISPEC) = MEANTE(ISPEC) + DFIM(IFRE) * SUM0(ISPEC)
          END DO
        END DO

!       ADD TAIL CONTRIBUTION
        FACTOR = WETAIL*FR(NFRE)*DELTH
        DO ISPEC = 1,NSPEC
          MEANTE(ISPEC) = MEANTE(ISPEC) + FACTOR * SUM0(ISPEC)
        END DO

      ENDIF

!     COMPUTE MEAN ANGLES:
!     --------------------

      DO ISPEC = 1,NSPEC
        SUMX0(ISPEC) = 0.0_JWRB
        SUMY0(ISPEC) = 0.0_JWRB
      END DO
      M10P = MIN(MFRCUT + 1_JWIM, NFRE)
      C10 = (FRCUT-FR(MFRCUT))*DELTH-(FR(M10P)-FR(MFRCUT))*DELTH2
      D10 = (FRCUT - FR(MFRCUT))**2*DELTH2/(FR(M10P)-FR(MFRCUT))

      DO ISPEC = 1, NSPEC
        SUMX0(ISPEC) = 0.0_JWRB
        SUMY0(ISPEC) = 0.0_JWRB
        DO IANG = 1, NANG
          SIH = SINTH(IANG)
          COH = COSTH(IANG)
          DO IFRE = 1, MFRCUT
            IF (IFRE > 1) THEN
              IF (IFRE.NE.NFRE) THEN
                DF=(FR(IFRE+1)+FR(IFRE-1))*DELTH2
              ELSE
                DF=(FR(IFRE)+FR(IFRE-1))*DELTH2
              ENDIF
            ELSE
              DF=(FR(IFRE+1)+FR(IFRE))*DELTH2
            ENDIF
            SUMX0(ISPEC) = SUMX0(ISPEC)+SPECW(ISPEC,IANG,IFRE)*SIH*DF
            SUMY0(ISPEC) = SUMY0(ISPEC)+SPECW(ISPEC,IANG,IFRE)*COH*DF
          ENDDO
          C11 = SPECW(ISPEC,IANG,M10P)-SPECW(ISPEC,IANG,MFRCUT)
          SUMX0(ISPEC) = SUMX0(ISPEC)+C10*SPECW(ISPEC,IANG,MFRCUT)*SIH
          SUMX0(ISPEC) = SUMX0(ISPEC)+D10*C11*SIH
          SUMY0(ISPEC) = SUMY0(ISPEC)+C10*SPECW(ISPEC,IANG,MFRCUT)*COH
          SUMY0(ISPEC) = SUMY0(ISPEC)+D10*C11*COH
        ENDDO

        IF ( (SUMY0(ISPEC) /= 0.0_JWRB) .OR. (SUMX0(ISPEC) /= 0.0_JWRB) ) THEN
          MEANTANG(ISPEC) = ATAN2(SUMX0(ISPEC),SUMY0(ISPEC))
        ELSE
          MEANTANG(ISPEC) = 0.0_JWRB
        END IF

      END DO

!     COMPUTE MEAN FREQUENCIES:
!     -------------------------

      IF ( MFRCUT < NFRE ) THEN
        C10  = (FRCUT - FR(MFRCUT))*DELTH
        C10  = C10/FRCUT
        C11  = CO1*FR(MFRCUT)*DELTH2
        C11  = C11/FR(MFRCUT)
        D10  = (FRCUT - FR(MFRCUT))**2*DELTH2/(CO1*FR(MFRCUT))
        D10  = D10/FRCUT

        DO ISPEC = 1,NSPEC
          MEANTFRE(ISPEC) = 0.0_JWRB
        END DO
        DO ISPEC = 1,NSPEC
          DO IANG=1,NANG
            DO IFRE = 1,MFRCUT
              FACTOR = DFIM(IFRE) / FR(IFRE)
              MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + FACTOR*SPECW(ISPEC,IANG,IFRE)
            ENDDO
            D11 = SPECW(ISPEC,IANG,M10P)-SPECW(ISPEC,IANG,MFRCUT)
            MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + (C10-C11)*SPECW(ISPEC,IANG,MFRCUT)+D10*D11
          ENDDO
        ENDDO

        DO ISPEC = 1, NSPEC
          IF (MEANTFRE(ISPEC) /= 0.0_JWRB) THEN
            MEANTFRE(ISPEC) = MEANTE(ISPEC) / MEANTFRE(ISPEC)
          ELSE
            MEANTFRE(ISPEC) = 0.0_JWRB
          END IF
        END DO

      ELSE
        DO ISPEC = 1,NSPEC
          MEANTFRE(ISPEC) = EPSMIN
        END DO

        DO IFRE = 1,NFRE
          DO ISPEC = 1,NSPEC
            SUM0(ISPEC) = 0.0_JWRB
          END DO
          DO IANG = 1,NANG
            DO ISPEC = 1,NSPEC
              SUM0(ISPEC) = SUM0(ISPEC) + SPECW(ISPEC,IANG,IFRE)
            END DO
          END DO
          DO ISPEC = 1,NSPEC
            MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + DFIMOFR(IFRE)*SUM0(ISPEC)
          END DO
        END DO

        FACTOR = FRTAIL*DELTH
        DO ISPEC = 1,NSPEC
          MEANTFRE(ISPEC) = MEANTFRE(ISPEC) + FACTOR * SUM0(ISPEC)
          IF (MEANTFRE(ISPEC) /= 0.0_JWRB) THEN
            MEANTFRE(ISPEC) = MEANTE(ISPEC) / MEANTFRE(ISPEC)
          ELSE
            MEANTFRE(ISPEC) = 0.0_JWRB
          END IF
        END DO

      ENDIF

      END SUBROUTINE TRANSMEANS
