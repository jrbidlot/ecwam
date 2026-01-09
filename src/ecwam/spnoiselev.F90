! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SPNOISELEV (KIJS, KIJL, WSWAVE, WDWAVE, CICOVER, FLM) 
! ----------------------------------------------------------------------

!**** *SPNOISELEV* - SETS THE SPECTRAL NOISE LEVEL

! ----------------------------------------------------------------------
USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOWFRED  , ONLY : TH
USE YOWICE   , ONLY : FLMIN
USE YOWPARAM , ONLY : NANG
! ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WSWAVE     !! WIND SPEED
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE     !! WIND DIRECTION
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICOVER    !! SEA ICE COVER 
REAL(KIND=JWRB), DIMENSION(KIJL,NANG), INTENT(OUT) :: FLM  !! SPECTRAL NOISE LEVEL 


INTEGER(KIND=JWIM) :: IJ, K

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB), DIMENSION(KIJL) :: ZRDC

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPNOISELEV',0,ZHOOK_HANDLE)

DO IJ=KIJS,KIJL
  IF (CICOVER(IJ) > 0.0_JWRB) THEN
    ! still allow noise in full sea ice cover, but only ten percent
    ZRDC(IJ) = (1._JWRB - 0.9_JWRB*MIN(CICOVER(IJ),0.99_JWRB))*FLMIN
  ELSE
    ! Reduce it for low winds (not over sea ice for now)
    ZRDC(IJ) = (MIN(WSWAVE(IJ),3._JWRB)/3._JWRB)*FLMIN
  ENDIF
ENDDO

DO K=1,NANG
  DO IJ=KIJS,KIJL
    FLM(IJ,K) = ZRDC(IJ)*MAX(0.0_JWRB, COS(TH(K)-WDWAVE(IJ)) )**2
  ENDDO
ENDDO

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPNOISELEV',1,ZHOOK_HANDLE)

END SUBROUTINE SPNOISELEV
