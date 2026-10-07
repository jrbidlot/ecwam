! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MCURP (KIJS, KIJL, IFROMIJ, JFROMIJ, &
&                 NXS, NXE, NYS, NYE, FIELDG,   &
&                 CICOVER, NEMOUCUR, NEMOVCUR,  &
&                 UCUR, VCUR)

!****  *MCURP* - TRANSFER CURRENTS FROM NEMO2WAM TO WVENVI AND CHECK FOR LAKES, SEA ICE AND MAXIMUM VALUE ALLOWED.

! ------------------------------------------------------------------- 

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

USE YOWCOUP  , ONLY : LWCOU
USE YOWCURR  , ONLY : CURRENT_MAX

USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! --------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ, JFROMIJ
INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(IN) :: CICOVER
REAL(KIND=JWRO), DIMENSION (KIJS:KIJL), INTENT(IN) :: NEMOUCUR, NEMOVCUR
REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(OUT) :: UCUR, VCUR

INTEGER(KIND=JWIM) :: IJ, IX, JY

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! --------------------------------------------------------------------- 

IF (LHOOK) CALL DR_HOOK('MCURP',0,ZHOOK_HANDLE)

! TRANSFER CURRENTS FROM NEMO
! ---------------------------
IF (LWCOU) THEN
  DO IJ = KIJS, KIJL 
    IX = IFROMIJ(IJ)
    JY = JFROMIJ(IJ)
    IF (FIELDG%LKFR(IX,JY) <=  0.0_JWRB ) THEN
!     if lake cover = 0, we assume open ocean point, then get currents directly from NEMO 
!     In sea ice, the currents will be reduced as it not clear how wave - sea ice -current interaction works
      UCUR(IJ) = SIGN(MIN(ABS(NEMOUCUR(IJ)),REAL(CURRENT_MAX,JWRO)), NEMOUCUR(IJ))
      UCUR(IJ) = (1.0_JWRB - CICOVER(IJ)) * UCUR(IJ) 

      VCUR(IJ) = SIGN(MIN(ABS(NEMOVCUR(IJ)),REAL(CURRENT_MAX,JWRO)), NEMOVCUR(IJ)) 
      VCUR(IJ) = (1.0_JWRB - CICOVER(IJ)) * VCUR(IJ)

    ELSE
!     no currents over lakes and land
      UCUR(IJ) = 0.0_JWRB
      VCUR(IJ) = 0.0_JWRB
    ENDIF
  ENDDO
ELSE
  UCUR(:) = NEMOUCUR(:)
  VCUR(:) = NEMOVCUR(:)
ENDIF

IF (LHOOK) CALL DR_HOOK('MCURP',1,ZHOOK_HANDLE)

END SUBROUTINE MCURP
