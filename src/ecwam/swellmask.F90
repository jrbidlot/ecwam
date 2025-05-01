! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SWELLMASK (KIJS, KIJL, FL1, XLLWS, CINV, UFRIC, COSWDIF, SWM)

! ----------------------------------------------------------------------

!**** *SWELLMASK* - DETERNINES A MASK FOR SEPARATING WINDSEA AND SWELL.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : FRIC     ,OLDWSFC
USE YOWPARAM , ONLY : NANG     ,NFRE

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL                      !! INDEX OF FIRST GRIDPOINT AND LAST GRID POINT
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1     !! SPECTRA
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: XLLWS   !! WINDSEA MASK FROM INPUT SOURCE TERM
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: CINV         !! INVERSE PHASE SPEED
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC             !! FRICTION VELOCITY IN M/S.
REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(IN) :: COSWDIF      !! COSINE (WDWAVE - WAVES DIRECTIONS)
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: SWM    !! MASK = 1 IF SWELL = 0 IF WINDSEA


INTEGER(KIND=JWIM) :: IJ, K, M

REAL(KIND=JWRB) :: COEF, CHECKTA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: XINVWVAGE
REAL(KIND=JWRB), DIMENSION(KIJL,NANG) :: DIRCOEF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SWELLMASK',0,ZHOOK_HANDLE)


DO M=1,NFRE
  DO IJ=KIJS,KIJL
    XINVWVAGE(IJ,M)=UFRIC(IJ)*CINV(IJ,M)
  ENDDO
ENDDO

COEF = OLDWSFC*FRIC
DO K=1,NANG
  DO IJ=KIJS,KIJL
    DIRCOEF(IJ,K)=COEF*COSWDIF(IJ,K)
  ENDDO
ENDDO

DO M=1,NFRE
  DO K=1,NANG
    DO IJ=KIJS,KIJL
      IF (XLLWS(IJ,K,M) /= 0.0_JWRB) THEN
        ! this is windsea 
        SWM(IJ,K,M)=0.0_JWRB
      ELSE
        CHECKTA=XINVWVAGE(IJ,M)*DIRCOEF(IJ,K)
        IF (CHECKTA >= 1.0_JWRB) THEN
          ! this is extra windsea 
          SWM(IJ,K,M)=0.0_JWRB
        ELSE
          ! this is swell
          SWM(IJ,K,M)=1.0_JWRB
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SWELLMASK',1,ZHOOK_HANDLE)

END SUBROUTINE SWELLMASK
