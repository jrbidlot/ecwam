! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

      SUBROUTINE SWRT_ENH (KIJS, KIJL, FL1, WAVNUM, DEPTH, AKMEAN, ENH)

! ----------------------------------------------------------------------

!     SHALLOW WATER SCALING ENH(IJ,M), M = 1..NFRE, FOR NONLINEAR SOURCE TERMS

      USE PARKIND_WAVE, ONLY : JWIM, JWRB

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : ISNONLIN

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peak_ang.intfb.h"
#include "transf.intfb.h"
#include "transf_snl.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: DEPTH, AKMEAN
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(OUT) :: ENH

      INTEGER(KIND=JWIM) :: IJ, M
      REAL(KIND=JWRB), PARAMETER :: ENH_MAX=10.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: ENH_MIN=0.1_JWRB
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ENHFR, XNU, SIG_TH

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

     IF (LHOOK) CALL DR_HOOK('SNL_ENH',0,ZHOOK_HANDLE)

      SELECT CASE (ISNONLIN)
      CASE(0)
        DO IJ=KIJS,KIJL
          ENHFR(IJ) = MAX(0.75_JWRB*DEPTH(IJ)*AKMEAN(IJ), 0.5_JWRB)
          ENHFR(IJ) = 1.0_JWRB+(5.5_JWRB/ENHFR(IJ)) * (1.0_JWRB-.833_JWRB*ENHFR(IJ)) * EXP(-1.25_JWRB*ENHFR(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=KIJS,KIJL
            ENH(IJ,M) = ENHFR(IJ)
          ENDDO
        ENDDO

      CASE(1)
        DO M=1,NFRE
          DO IJ = KIJS, KIJL
            ENH(IJ,M) = MAX(MIN(ENH_MAX,TRANSF(WAVNUM(IJ,M),DEPTH(IJ))),ENH_MIN)
          ENDDO
        ENDDO

      CASE(2)
        CALL PEAK_ANG(KIJS, KIJL, FL1, XNU, SIG_TH)
        DO M=1,NFRE
          DO IJ = KIJS, KIJL
            ENH(IJ,M) = TRANSF_SNL(WAVNUM(IJ,M), DEPTH(IJ), XNU(IJ), SIG_TH(IJ))
          ENDDO
        ENDDO
      END SELECT

      IF (LHOOK) CALL DR_HOOK('SNL_ENH',1,ZHOOK_HANDLE)

      END SUBROUTINE SWRT_ENH
