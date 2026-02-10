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

      USE YOWFRED  , ONLY : FR, FR5, DFIM, DELTH
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

      REAL(KIND=JWRB) :: FCUTB, FCUTB_FT, FCUTT, FBOT, FTOP, ZW
      REAL(KIND=JWRB) :: DELF, WL, WR
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFIMLOC
      REAL(KIND=JWRB), DIMENSION(KIJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',0,ZHOOK_HANDLE)

      DFIMLOC(:)=DFIM(:)

!     CORRECT DFIMLOC for both side of the integral

!     Left side:
      FBOT = 1.0_JWRB/MAX(TT,EPSMIN)
      FCUTB_FT = MIN(FBOT,FR(NFRE))
      FCUTB = MAX(FR(1),FCUTB_FT)
      FBOT = MAX(FBOT,FR(NFRE))   !! FBOT is used if a tail contribution is needed

      MCUTB=1
!!!!  MCUTB is defined such that the integration will need to be from index MAX(MCUTB-1,1)
      DO WHILE (FR(MCUTB) < FCUTB .AND. MCUTB < NFRE )
        MCUTB = MCUTB+1
      ENDDO

      IF (MCUTB > 1) THEN
        DELF = FR(MCUTB) - FR(MCUTB-1)
        WL = (FR(MCUTB) - FCUTB)/DELF
        WR = 1.0_JWRB - WL

        DFIMLOC(MCUTB-1) = DELTH * 0.5_JWRB * (FR(MCUTB)-FCUTB) * WL
        DFIMLOC(MCUTB) = DELTH * 0.5_JWRB * ( (FR(MCUTB)-FCUTB) * (1.0_JWRB + WR) + (FR(MCUTB+1)-FR(MCUTB)) )
      ENDIF


!     Right side:
      FTOP = 1.0_JWRB/MAX(TB,EPSMIN)
      FCUTT = MAX(FR(1),MIN(FTOP,FR(NFRE)))
      FTOP = MAX(FTOP,FR(NFRE))   !! FTOP is used if a tail contribution is needed

      MCUTT=1
!!!!  MCUTT is defined such that the integration will need to be up to index MIN(MCUTT,NFRE)
      DO WHILE ( FCUTT > FR(MCUTT) .AND. MCUTT < NFRE )
        MCUTT = MCUTT+1
      ENDDO


      DELF = FR(MCUTT) - FR(MCUTT-1)
      WL = (FR(MCUTT)-FCUTT)/DELF
      WR = 1.0_JWRB - WL

      DFIMLOC(MCUTT-1) = DELTH * 0.5_JWRB * ( (FR(MCUTT-1)-FR(MCUTT-2)) + (FRFCUTT-FR(MCUTT-1)) * (1.0_JWRB + WL) )
      DFIMLOC(MCUTT) = DELTH * 0.5_JWRB * (FCUTT-FR(MCUTT-1)) * WR 


      IF(FCUTB == FCUTT) MCUTT=MCUTB-1

!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=KIJS,KIJL
        EBT(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=MAX(MCUTB-1,1),MCUTT
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EBT(IJ) = EBT(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

!     ADD SMALL CONTRIBUTION FOR A LINEAR FRONT TAIL IF NEEDED
      IF (FCUTB_FT < FCUTB .AND. FCUTB == FR(1)) THEN
        WL = (FR(1) - FCUTB_FT) / FR(1)
        WR = 1.0_JWRB - WL
        DELF = 0.5_JWRB * (FR(1) - FCUTB_FT) * (1.0_JWRB + WR)
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,1)*DELTH
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,1)*DELTH
          ENDDO
        ENDDO
        DO IJ = KIJS, KIJL
          EBT(IJ) = EBT(IJ) + DELF*TEMP(IJ)
        ENDDO
      ENDIF


!     CHECK IF A THE REQUESTED FREQUENCIES ARE ABOVE FR(NFRE)
      IF ( FBOT < FTOP ) THEN

        ZW = 0.25_JWRB * DELTH * FR5(NFRE) * ( 1.0_JWRB/FBOT**4 - 1.0_JWRB/FTOP**4 )
        K=1
        DO IJ=KIJS,KIJL
          TEMP(IJ) = FL1(IJ,K,NFRE)
        ENDDO
        DO K=2,NANG
          DO IJ=KIJS,KIJL
            TEMP(IJ) = TEMP(IJ)+FL1(IJ,K,NFRE)
          ENDDO
        ENDDO
        DO IJ=KIJS,KIJL
          EBT(IJ) = EBT(IJ) + ZW*TEMP(IJ)
        ENDDO

      ENDIF 

      IF (LHOOK) CALL DR_HOOK('SEBTMEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SEBTMEAN
