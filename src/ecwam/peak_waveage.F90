! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PEAK_WAVEAGE (KIJS, KIJL, FL1, DEPTH, UFRIC, ZPWAVEAGE)

! ----------------------------------------------------------------------

!**** *PEAK_WAVEAGE* - DETERNINES THE WAVEAGE AT THE PEAK FREQUENCY

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWPARAM , ONLY : NANG     ,NFRE
USE YOWPCONS , ONLY : ZPI      ,EPSUS

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------

IMPLICIT NONE

#include "aki.intfb.h"
#include "dominant_period.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL                      !! INDEX OF FIRST GRIDPOINT AND LAST GRID POINT
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1     !! SPECTRA
REAL(KIND=JWRB), DIMENSION(KIJS), INTENT(IN) :: DEPTH             !! WATER DEPTH
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC             !! FRICTION VELOCITY IN M/S.
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: ZPWAVEAGE        !! WAVEAGE AT THE PEAK FREQUENCY

INTEGER(KIND=JWIM) :: IJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB) ::  ZOMEGAP, ZKP
REAL(KIND=JWRB), DIMENSION(KIJL) :: DP

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PEAK_WAVEAGE',0,ZHOOK_HANDLE)

CALL DOMINANT_PERIOD (KIJS, KIJL, FL1, DP)

DO IJ=KIJS,KIJL
  IF(DP(IJ) > 0.0_JWRB ) THEN
    ZOMEGAP = ZPI/DP(IJ)
    ZKP = AKI(ZOMEGAP,DEPTH(IJ)) 
    ZPWAVEAGE(IJ) = ZOMEGAP/(ZKP*MAX(UFRIC(IJ),EPSUS))
  ELSE
    ZPWAVEAGE(IJ) = 0.0_JWRB 
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('PEAK_WAVEAGE',1,ZHOOK_HANDLE)

END SUBROUTINE PEAK_WAVEAGE
