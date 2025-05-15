! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE FRCUTINDEX (KIJS, KIJL, FM, FMWS, UFRIC, CICOVER,     &
     &                       MIJ, FCUT, RHOWGDFTH)

! ----------------------------------------------------------------------

!**** *FRCUTINDEX* - RETURNS THE LAST FREQUENCY INDEX OF
!                    PROGNOSTIC PART OF SPECTRUM.

!**   INTERFACE.
!     ----------

!       *CALL* *FRCUTINDEX (KIJS, KIJL, FM, FMWS, CICOVER, MIJ, RHOWGDFTH)
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FM*     - MEAN FREQUENCY
!          *FMWS*   - MEAN FREQUENCY OF WINDSEA
!          *UFRIC*  - FRICTION VELOCITY IN M/S
!          *CICOVER*- CICOVER 
!          *MIJ*    - LAST FREQUENCY INDEX for imposing high frequency tail
!          *FCUT*   - THE ACTUAL FREQUENCY OF THE PROGNOSTIC PART OF SPECTRUM.
!          *RHOWGDFTH - WATER DENSITY * G * DF * DTHETA
!                       FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ) 
!                       !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)


!     METHOD.
!     -------

!*    COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE 2.5*MAX(FMWS,FM).


!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated (see module yowphys) !!!


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM       ,FRATIO   ,FLOGSPRDM1, &
     &                DELTH          ,RHOWG_DFIM ,FRIC
      USE YOWICE   , ONLY : CITHRSH_TAIL
      USE YOWPARAM , ONLY : NFRE
      USE YOWPCONS , ONLY : G        ,ZPI        ,EPSMIN   ,ROWATER
      USE YOWPHYS  , ONLY : TAILFACTOR, TAILFACTOR_PM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), INTENT(OUT) :: MIJ(KIJL)
      REAL(KIND=JWRB),DIMENSION(KIJL), INTENT(IN) :: FM, FMWS, UFRIC, CICOVER
      REAL(KIND=JWRB),DIMENSION(KIJL), INTENT(OUT) :: FCUT
      REAL(KIND=JWRB),DIMENSION(KIJL,NFRE), INTENT(OUT) :: RHOWGDFTH 


      INTEGER(KIND=JWIM) :: IJ, M

      REAL(KIND=JWRB) :: FPMH, FPPM, FM2, FPM, FPM4, ZLOG10FR1
      REAL(KIND=JWRB) :: ZHRWGDELTH, ZW1, ZRCUT, ZLOGFRATIOM1
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',0,ZHOOK_HANDLE)

!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

      FPMH = TAILFACTOR
      FPPM = TAILFACTOR_PM*G/(FRIC*ZPI)
      ZLOG10FR1 = LOG10(FR(1))

      DO IJ=KIJS,KIJL
        IF (CICOVER(IJ) <= CITHRSH_TAIL) THEN
          FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
          FPM = FPPM/MAX(UFRIC(IJ),EPSMIN)
          FPM4 = MAX(FM2,FPM)
          MIJ(IJ) = INT((LOG10(FPM4)-ZLOG10FR1)*FLOGSPRDM1)+2
          MIJ(IJ) = MIN(MAX(2,MIJ(IJ)),NFRE)
          FCUT(IJ) = FPM4
          FCUT(IJ) = MAX(MIN(FCUT(IJ),FR(NFRE)),FR(1))
        ELSE
          MIJ(IJ) = NFRE
          FCUT(IJ) = FR(NFRE)
        ENDIF
!!!silly test
       if ( FCUT(IJ) < FR(MIJ(IJ)-1) .OR. FCUT(IJ) > FR(MIJ(IJ)) ) then
         write (*,*) 'debile we have a problem ',FCUT(IJ), FR(MIJ(IJ)-1),FR(MIJ(IJ))
       endif
      ENDDO

!     SET RHOWGDFTH
      ZHRWGDELTH=0.5_JWRB*ROWATER*G*DELTH
      ZLOGFRATIOM1=1.0_JWRB/LOG(FRATIO)
      DO IJ=KIJS,KIJL
        DO M=1,MIJ(IJ)-2
          RHOWGDFTH(IJ,M) = RHOWG_DFIM(M)
        ENDDO
        ZW1 = LOG(FR(MIJ(IJ))/FCUT(IJ))*ZLOGFRATIOM1
        ZRCUT = LOG(FCUT(IJ)/FR(MIJ(IJ)-1))
        RHOWGDFTH(IJ,MIJ(IJ)-1)=0.5_JWRB*RHOWG_DFIM(MIJ(IJ)-1) + ZHRWGDELTH*FR(MIJ(IJ)-1)*ZRCUT*(1.0_JWRB+ZW1)
        RHOWGDFTH(IJ,MIJ(IJ))=ZHRWGDELTH*FR(MIJ(IJ))*ZRCUT*(1.0_JWRB-ZW1)

        DO M=MIJ(IJ)+1,NFRE
          RHOWGDFTH(IJ,M) = 0.0_JWRB
        ENDDO
      ENDDO

      IF (LHOOK) CALL DR_HOOK('FRCUTINDEX',1,ZHOOK_HANDLE)

      END SUBROUTINE FRCUTINDEX
