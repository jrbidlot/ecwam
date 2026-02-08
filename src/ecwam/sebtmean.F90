! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SEBTMEAN (KIJS, KIJL, FL1, TB, TT, EBT)

! ----------------------------------------------------------------------

!**** *SEBTMEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                  WITH  TB <= PERIODS <= TT
!                  IF NEEDED A HIGH FREQUENCY f**-5 TAIL IS ASSUMED !!!!

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SEBTMEAN(KIJS, KIJL, FL1, TB, TT, EBT)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1*     - SPECTRUM.
!          *TB*      - BOTTOM PERIOD
!          *TT*      - TOP PERIOD
!          *EBT*     - MEAN VARIANCE

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR, FR5, DELTH
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), INTENT(IN) :: TB, TT
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: EBT


      INTEGER(KIND=JWIM) :: IJ, M, K, MCUTB, MCUTT

      REAL(KIND=JWRB) :: FCUTB, FCUTT, FBOT, FTOP, ZW
      REAL(KIND=JWRB) :: WL, WR, DF
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: FRLOC 
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: F1D 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',0,ZHOOK_HANDLE)

      FBOT = 1.0_JWRB/MAX(TT,EPSMIN)
      FCUTB = MAX(FR(1),MIN(FBOT,FR(NFRE)))
      FBOT = MAX(FBOT,FR(NFRE))   !! FBOT is used if a tail contribution is needed

      MCUTB=1
      DO WHILE (FR(MCUTB) < FCUTB .AND. MCUTB < NFRE )
        MCUTB = MCUTB+1
      ENDDO

      FTOP = 1.0_JWRB/MAX(TB,EPSMIN)
      FCUTT = MAX(FR(1),MIN(FTOP,FR(NFRE)))
      FTOP = MAX(FTOP,FR(NFRE))   !! FTOP is used if a tail contribution is needed

      MCUTT=NFRE
      DO WHILE (FR(MCUTT) > FCUTB .AND. MCUTT > 1 )
        MCUTT = MCUTT-1
      ENDDO

      IF(FCUTB == FCUTT) MCUTT=MCUTB-1

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=KIJS,KIJL
        EBT(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

!     FREQUENCY SPECTRUM
      IF (MCUTB > 1) THEN
        FRLOC(MCUTB-1) = FCUTB
        WL = (FR(MCUTB) - FCUTB) / (FR(MCUTB) - FR(MCUTB-1))
        WR = 1.0_JWRB - WL
        K=1
        ! Linear interpolation
        DO IJ = KIJS, KIJL
          F1D(IJ,MCUTB-1) = ( WL*FL1(IJ,K,MCUTB-1) + WR*FL1(IJ,K,MCUTB) )*DELTH
        ENDDO
        DO K = 2, NANG
          DO IJ = KIJS, KIJL
            F1D(IJ,MCUTB-1) = F1D(IJ,MCUTB-1) + ( WL*FL1(IJ,K,MCUTB-1) + WR*FL1(IJ,K,MCUTB) )*DELTH
          ENDDO
        ENDDO
      ENDIF

      DO M = MCUTB, MCUTT
        FRLOC(M) = FR(M)
        K=1
        DO IJ = KIJS, KIJL
          F1D(IJ,M) = FL1(IJ,K,M)*DELTH
        ENDDO
        DO K = 2, NANG
          DO IJ = KIJS, KIJL
            F1D(IJ,N) = F1D(IJ,M)+FL1(IJ,K,M)*DELTH
          ENDDO
        ENDDO
      ENDDO

      IF (MCUTT < NFRE) THEN
        FRLOC(MCUTT+1) = FCUTT
        WL = (FR(MCUTT+1) - FCUTT) / (FR(MCUTT+1) - FR(MCUTT))
        WR = 1.0_JWRB - WL
        K=1
        ! Linear interpolation
        DO IJ = KIJS, KIJL
          F1D(IJ,MCUTT+1) = ( WL*FL1(IJ,K,MCUTT) + WR*FL1(IJ,K,MCUTT+1) )*DELTH
        ENDDO
        DO K = 2, NANG
          DO IJ = KIJS, KIJL
            F1D(IJ,MCUTT+1) = F1D(IJ,MCUTT+1) + ( WL*FL1(IJ,K,MCUTT) + WR*FL1(IJ,K,MCUTT+1) )*DELTH
          ENDDO
        ENDDO
      ENDIF


!     TRAPEZOIDAL INTERPOLATION of F1D between FRLOC(MCUTB) and FRLOC(MCUTT)
      DO M = MAX(MCUTB-1,1), MIN(MCUTT,NFRE-1)
        DF = 0.5_JWRB*(FRLOC(M+1)-FRLOC(M))
        DO IJ = KIJS, KIJL
          EBT(IJ) = EBT(IJ) + DF*(F1D(IJ,M+1)+F1D(IJ,M))
        ENDDO
      ENDDO

!     CHECK IF A THE REQUESTED FREQUENCIES ARE ABOVE FR(NFRE)
      IF ( FBOT < FTOP ) THEN

        ZW = 0.25_JWRB * DELTH * FR5(NFRE) * ( 1.0_JWRB/FBOT**4 - 1.0_JWRB/FTOP**4 )
        K=1
        DO IJ = KIJS, KIJL
          TEMP(IJ) = FL1(IJ,K,NFRE)
        ENDDO
        DO K = 2, NANG
          DO IJ = KIJS, KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,NFRE)
          ENDDO
        ENDDO
        DO IJ = KIJS, KIJL
          EBT(IJ) = EBT(IJ) + ZW*TEMP(IJ)
        ENDDO

      ENDIF 

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SEBTMEAN
