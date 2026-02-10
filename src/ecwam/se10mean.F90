! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SE10MEAN (KIJS, KIJL, FL1, E10)

! ----------------------------------------------------------------------

!**** *SE10MEAN* - COMPUTATION OF SPECTRAL VARIANCE FOR ALL WAVES
!                 WITH PERIOD LARGER THAN 10s

!*    PURPOSE.
!     --------

!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *SE10MEAN(KIJS, KIJL, FL1, E10)*
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *FL1      - SPECTRUM.
!          *E10*     - MEAN ENERGY

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

      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL, NANG, NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: E10


      INTEGER(KIND=JWIM) :: IJ, M, K, MCUT

      REAL(KIND=JWRB), PARAMETER :: FCUT=0.1_JWRB
      REAL(KIND=JWRB) :: DELF, WL, WR
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NFRE) :: DFIMLOC
      REAL(KIND=JWRB), DIMENSION(KIJL) :: TEMP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',0,ZHOOK_HANDLE)

      MCUT=1
      DO WHILE ( FCUT > FR(MCUT) .AND. MCUT < NFRE )
        MCUT = MCUT+1
      ENDDO
      

      DFIMLOC(:)=DFIM(:)

      DELF = FR(MCUT) - FR(MCUT-1)
      WL = (FR(MCUT)-FCUT)/DELF
      WR = 1.0_JWRB - WL

      DFIMLOC(MCUT-1) = DELTH * 0.5_JWRB * ( (FR(MCUT-1)-FR(MCUT-2)) + (FRFCUT-FR(MCUT-1)) * (1.0_JWRB + WL) )
      DFIMLOC(MCUT) = DELTH * 0.5_JWRB * (FCUT-FR(MCUT-1)) * WR 





!*    1. INITIALISE ENERGY ARRAY.
!        ------------------------

      DO IJ=KIJS,KIJL
        E10(IJ) = EPSMIN 
      ENDDO

! ----------------------------------------------------------------------

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTION.
!        -----------------------------------------

      DO M=1,MCUT
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
          E10(IJ) = E10(IJ)+DFIMLOC(M)*TEMP(IJ)
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('SE10MEAN',1,ZHOOK_HANDLE)

      END SUBROUTINE SE10MEAN
