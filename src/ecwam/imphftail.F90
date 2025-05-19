! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE IMPHFTAIL (KIJS, KIJL, MIJ, FCUT, FLM, WAVNUM, XK2CG, FL1) 
! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM


!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ


!**   INTERFACE.
!     ----------

!       *CALL* *IMPHFTAIL (KIJS, KIJL, MIJ, FCUT, FLM, WAVNUM, XK2CG, FL1)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*     - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *FCUT*    - ACTUAL FREQUENCY OF THE PROGNOSTIC RANGE,
!          *FLM*     - SPECTAL DENSITY MINIMUM VALUE
!          *WAVNUM*  - WAVENUMBER
!          *XK2CG*   - (WAVNUM)**2 * GROUP SPEED
!          *FL1*     - SPECTRUM (INPUT AND OUTPUT).

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWFRED  , ONLY : FR       ,FRM5  ,FRATIO,  DELTH
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : EPSMIN
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: FCUT 
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(IN) :: FLM 
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, XK2CG
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: TEWHMIN, TEWHMAX
      REAL(KIND=JWRB), DIMENSION(KIJL) :: TEMP1, TEMP2
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ZW1, ZSCL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',0,ZHOOK_HANDLE)

!*    DIAGNOSTIC TAIL.
!     ----------------

!     APPLY F**-5 TAIL FROM FCUT WHENEVER FCUT < FR(MIJ)

      DO IJ=KIJS,KIJL
        ZSCL(IJ) =  FCUT(IJ)**5 * FRM5(MIJ(IJ))
        ZW1(IJ) = (FR(MIJ(IJ))-FCUT(IJ))/(FR(MIJ(IJ)) - FR(MIJ(IJ)-1))
      ENDDO
      DO K=1,NANG
        DO IJ=KIJS,KIJL
          FL1(IJ,K,MIJ(IJ)) = (ZW1(IJ)*FL1(IJ,K,MIJ(IJ)-1) + (1.0_JWRB-ZW1(IJ))*FL1(IJ,K,MIJ(IJ))) * ZSCL(IJ)
        ENDDO
      ENDDO



!*    MERGE TAIL INTO SPECTRA.
!     ------------------------
      DO IJ=KIJS,KIJL
        TEMP1(IJ) = XK2CG(IJ,MIJ(IJ))*WAVNUM(IJ,MIJ(IJ))
      ENDDO
      DO IJ=KIJS,KIJL
        DO M=MIJ(IJ)+1,NFRE
          TEMP2(IJ) = TEMP1(IJ)/(XK2CG(IJ,M)*WAVNUM(IJ,M))
          DO K=1,NANG
            FL1(IJ,K,M) = MAX(TEMP2(IJ)*FL1(IJ,K,MIJ(IJ)),FLM(IJ,K))
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('IMPHFTAIL',1,ZHOOK_HANDLE)

      END SUBROUTINE IMPHFTAIL
