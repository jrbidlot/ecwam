! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WVFLD2IFS (KIJS, KIJL, FL1, DEPTH, WAVNUM, WSWAVE, WDWAVE, UFRIC, &
 &                    ZPWVAGE_U10, ZPWVAGE_US, ZMEANSS)

! ----------------------------------------------------------------------

!**** *WVFLD2IFS* - COMPUTES EXTRA FIELDS POTENTIALLY PASSED BACK TO IFS.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWFRED  , ONLY : FR       ,TH       , XKMSS_CUTOFF_IFS
USE YOWPARAM , ONLY : NANG     ,NFRE
USE YOWPCONS , ONLY : ZPI      ,EPSUS

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -----------------------------------------------------------------------

IMPLICIT NONE

#include "meansqs.intfb.h"
#include "peak_waveage.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL                      !! INDEX OF FIRST GRID POINT AND LAST GRID POINT
REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1     !! SPECTRA
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH             !! WATER DEPTH
REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM       !! WAVE NUMBER
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WSWAVE            !! 10m WIND SPEED in M/S
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: WDWAVE            !! 10m WIND DIRECTION 
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: UFRIC             !! FRICTION VELOCITY IN M/S
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: ZPWVAGE_U10      !! WAVEAGE AT THE PEAK FREQUENCY BASED ON U10
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: ZPWVAGE_US       !! WAVEAGE AT THE PEAK FREQUENCY BASED ON USTAR
REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: ZMEANSS          !! MEAN SQUARE SLOPE WITH CUTOFF WAVE NUMBER XKMSS_CUTOFF_IFS 
 

INTEGER(KIND=JWIM) :: K

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JWRB), DIMENSION(KIJL, NANG) :: COSWDIF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WVFLD2IFS',0,ZHOOK_HANDLE)

DO K=1,NANG
  COSWDIF(:,K) = COS(TH(K)-WDWAVE(:))
ENDDO

!! WAVE AGE(s)
CALL PEAK_WAVEAGE (KIJS, KIJL, FL1, DEPTH, WSWAVE, UFRIC, ZPWVAGE_U10, ZPWVAGE_US)

!! MEAN SQUARE SLOPE
CALL MEANSQS (XKMSS_CUTOFF_IFS, KIJS, KIJL, FL1, WAVNUM, UFRIC, WSWAVE, COSWDIF, ZMEANSS)


IF (LHOOK) CALL DR_HOOK('WVFLD2IFS',1,ZHOOK_HANDLE)

END SUBROUTINE WVFLD2IFS
