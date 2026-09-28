! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPAGS2 (F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)

! ----------------------------------------------------------------------

!**** *PROPAGS2* -  ADVECTION USING THE CORNER TRANSPORT SCHEME IN SPACE

!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS2(F1, F3, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)*
!          *F1*          - SPECTRUM AT TIME T (with exchange halo).
!          *F3*          - SPECTRUM AT TIME T+DELT
!          *NINF:NSUP+1* - 1st DIMENSION OF F1 and F3
!          *KIJS*        - ACTIVE INDEX OF FIRST POINT
!          *KIJL*        - ACTIVE INDEX OF LAST POINT
!          *NANG*        - NUMBER OF DIRECTIONS
!          *ND3SF1*      - LOWER 3rd DIMENSION OF F1 
!          *ND3EF1*      - UPPER 3d DIMENSION OF F1 
!          *ND3S*        - FREQUENCY INDEX SOLVED BY THIS CALL ND3S:ND3E
!          *ND3E*        - FREQUENCY INDEX SOLVED BY THIS CALL ND3S:ND3E

!     METHOD.
!     -------
!
!       CORNER TRANSPORT PROPAGATION SCHEME.
!
!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!      See IFS Documentation, part VII 

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWSTAT  , ONLY : ICASE    ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR      ,             &
     &            WLATN    ,WLONN    ,WCORN    ,WKPMN    ,WMPMN     ,   &
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN   ,   &
     &            SUMWN    ,                                            &
     &            JXO      ,JYO      ,KCR      ,KPM      ,MPM

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT , ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

! ----------------------------------------------------------------------
!     ARGUMENTS.
! ----------------------------------------------------------------------

      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, ND3SF1:ND3EF1), INTENT(IN) :: F1
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, ND3S:ND3E), INTENT(OUT) :: F3 

      INTEGER(KIND=JWIM), INTENT(IN) :: NINF
      INTEGER(KIND=JWIM), INTENT(IN) :: NSUP
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS
      INTEGER(KIND=JWIM), INTENT(IN) :: KIJL
      INTEGER(KIND=JWIM), INTENT(IN) :: NANG
      INTEGER(KIND=JWIM), INTENT(IN) :: ND3SF1
      INTEGER(KIND=JWIM), INTENT(IN) :: ND3EF1
      INTEGER(KIND=JWIM), INTENT(IN) :: ND3S
      INTEGER(KIND=JWIM), INTENT(IN) :: ND3E

! ----------------------------------------------------------------------
!     LOOP INDICES.
! ----------------------------------------------------------------------

      INTEGER(KIND=JWIM) :: K
      INTEGER(KIND=JWIM) :: M
      INTEGER(KIND=JWIM) :: IJ
      INTEGER(KIND=JWIM) :: IC
      INTEGER(KIND=JWIM) :: ICL
      INTEGER(KIND=JWIM) :: ICR

! ----------------------------------------------------------------------
!     PRIVATE INDICES, CONSTANT OVER THE IJ VECTOR LOOP.
! ----------------------------------------------------------------------

      INTEGER(KIND=JWIM) :: KCR1
      INTEGER(KIND=JWIM) :: KCR2
      INTEGER(KIND=JWIM) :: KCR3
      INTEGER(KIND=JWIM) :: KCR4

      INTEGER(KIND=JWIM) :: KPMM
      INTEGER(KIND=JWIM) :: KPMP

! ----------------------------------------------------------------------
!     LOCAL INTEGER VARIABLES.
! ----------------------------------------------------------------------

     INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,2)    :: LOC_KLON 
     INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,2,2)  :: LOC_KLAT 
     INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,4,2)  :: LOC_KCOR 
     INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,-1:1) :: LOC_KPM 
     INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL,-1:1) :: LOC_MPM 

! ----------------------------------------------------------------------
!     LOCAL REAL VARIABLES.
! ----------------------------------------------------------------------

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAGS2',0,ZHOOK_HANDLE)

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE == 1) THEN

!*      SPHERICAL GRID.
!       ---------------

        IF (IREFRA /= 2 .AND. IREFRA /= 3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------

          !$acc parallel loop independent collapse(3) &
          !$acc & present(F1,F3,KLON,KLAT,KCOR,SUMWN,WLONN,WLATN,WCORN,JXO,JYO,KCR,WKPMN,KPM)
          DO K = 1, NANG
            DO M = ND3S, ND3E

!DIR$ IVDEP
!DIR$ PREFERVECTOR
              DO IJ = KIJS, KIJL
                F3(IJ,K,M) =                                            &
     &                (1.0_JWRB-SUMWN(IJ,K,M))* F1(IJ           ,K  ,M) &

     &         + WLONN(IJ,K,M,JXO(K,1))   * F1(KLON(IJ,JXO(K,1))  ,K  ,M) &
     &         + WLATN(IJ,K,M,JYO(K,1),1) * F1(KLAT(IJ,JYO(K,1),1),K  ,M) &
     &         + WLATN(IJ,K,M,JYO(K,1),2) * F1(KLAT(IJ,JYO(K,1),2),K  ,M) &
     &         + WCORN(IJ,K,M,1,1)        * F1(KCOR(IJ,KCR(K,1),1),K  ,M) &
     &         + WCORN(IJ,K,M,1,2)        * F1(KCOR(IJ,KCR(K,1),2),K  ,M) &
     &         + WKPMN(IJ,K,M,-1)         * F1(IJ,KPM(K,-1),M)            &
     &         + WKPMN(IJ,K,M, 1)         * F1(IJ,KPM(K, 1),M)
              ENDDO

            ENDDO
          ENDDO
          !$acc end parallel loop 

        ELSE

!*      DEPTH AND/OR CURRENT REFRACTION.
!       -----------------------------

!*        DEPTH AND/OR CURRENT REFRACTION.
!         --------------------------------
!
          !$acc parallel                                                &
          !$acc & present(F1,F3,SUMWN,WLONN,KLON,WLATN,KLAT,WCORN,KCOR, &
          !$acc &         WKPMN,KPM,WMPMN,MPM)


          !$acc loop gang collapse(2)                                   &
          !$acc& private(LOC_KLON,LOC_KLAT,LOC_KCOR,LOC_KPM,LOC_MPM,    &
          !$acc&         KCR1,KCR2,KCR3,KCR4,KPMM,KPMP)

          DO M = ND3S, ND3E
            DO K = 1, NANG

!*            Local pointers to halo points (if not use point to the central point).
              DO IC=1,2
                IF (LLWLONN(K,M,IC)) THEN
                  LOC_KLON(KIJS:KIJL,IC) = KLON(KIJS:KIJL,IC)
                ELSE
                  LOC_KLON(KIJS:KIJL,IC) = [(IJ, IJ=KIJS,KIJL)]
                ENDIF
              ENDDO

              DO ICL=1,2

                DO IC=1,2
                  IF (LLWLATN(K,M,IC,ICL)) THEN
                    LOC_KLAT(KIJS:KIJL,IC,ICL) = KLAT(KIJS:KIJL,IC,ICL)
                  ELSE
                    LOC_KLAT(KIJS:KIJL,IC,ICL) = [(IJ, IJ=KIJS,KIJL)]
                  ENDIF
                ENDDO

                DO ICR=1,4
                  IF (LLWCORN(K,M,ICR,ICL)) THEN
                    LOC_KCOR(KIJS:KIJL,KCR(K,ICR),ICL) = KCOR(KIJS:KIJL,KCR(K,ICR),ICL)
                  ELSE
                    LOC_KCOR(KIJS:KIJL,KCR(K,ICR),ICL) = [(IJ, IJ=KIJS,KIJL)]
                  ENDIF
                ENDDO
              ENDDO

              DO IC=-1,1,2
                IF (LLWKPMN(K,M,IC)) THEN
                  LOC_KPM(KIJS:KIJL,IC) = KPM(K,IC)
                ELSE
                  LOC_KPM(KIJS:KIJL,IC) = K
                ENDIF

                IF (LLWMPMN(K,M,IC)) THEN
                  LOC_MPM(KIJS:KIJL,IC) = MPM(M,IC)
                ELSE
                  LOC_MPM(KIJS:KIJL,IC) = M
                ENDIF
              ENDDO

!*            K-dependent neighbour indices.

              KCR1 = KCR(K,1)
              KCR2 = KCR(K,2)
              KCR3 = KCR(K,3)
              KCR4 = KCR(K,4)

              KPMM = KPM(K,-1)
              KPMP = KPM(K, 1)


!DIR$ IVDEP
!DIR$ PREFERVECTOR
              DO IJ = KIJS, KIJL

                F3(IJ,K,M) = (1.0_JWRB-SUMWN(IJ,K,M))*F1(IJ,K,M) &

                  &         + WLONN(IJ,K,M,1)*F1(LOC_KLON(IJ,1),K,M) &
                  &         + WLONN(IJ,K,M,2)*F1(LOC_KLON(IJ,2),K,M) &

                  &         + WLATN(IJ,K,M,1,1)*F1(LOC_KLAT(IJ,1,1),K,M) &
                  &         + WLATN(IJ,K,M,2,1)*F1(LOC_KLAT(IJ,2,1),K,M) &
                  &         + WCORN(IJ,K,M,1,1)*F1(LOC_KCOR(IJ,KCR1,1),K,M) &
                  &         + WCORN(IJ,K,M,2,1)*F1(LOC_KCOR(IJ,KCR2,1),K,M) &
                  &         + WCORN(IJ,K,M,3,1)*F1(LOC_KCOR(IJ,KCR3,1),K,M) &
                  &         + WCORN(IJ,K,M,4,1)*F1(LOC_KCOR(IJ,KCR4,1),K,M) &

                  &         + WLATN(IJ,K,M,1,2)*F1(LOC_KLAT(IJ,1,2),K,M) & 
                  &         + WLATN(IJ,K,M,2,2)*F1(LOC_KLAT(IJ,2,2),K,M) &
                  &         + WCORN(IJ,K,M,1,2)*F1(LOC_KCOR(IJ,KCR1,2),K,M) &
                  &         + WCORN(IJ,K,M,2,2)*F1(LOC_KCOR(IJ,KCR2,2),K,M) &
                  &         + WCORN(IJ,K,M,3,2)*F1(LOC_KCOR(IJ,KCR3,2),K,M) &
                  &         + WCORN(IJ,K,M,4,2)*F1(LOC_KCOR(IJ,KCR4,2),K,M) &

                  &         + WKPMN(IJ,K,M,-1)*F1(IJ,LOC_KPM(IJ,-1),M) &
                  &         + WMPMN(IJ,K,M,-1)*F1(IJ,K,LOC_MPM(IJ,-1)) &

                  &         + WKPMN(IJ,K,M, 1)*F1(IJ,LOC_KPM(IJ, 1),M) &
                  &         + WMPMN(IJ,K,M, 1)*F1(IJ,K,LOC_MPM(IJ, 1))

              ENDDO

            ENDDO
          ENDDO

          !$acc end parallel

        ENDIF

      ELSE

!*      CARTESIAN GRID.
!       ---------------

        IF (IREFRA == 2 .OR. IREFRA == 3) THEN

!*        DEPTH AND/OR CURRENT REFRACTION REQUESTED.
!         ------------------------------------------

          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *'
          WRITE (IU06,*) '* FOR CARTESIAN GRID !                   *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'

          CALL ABORT1

        ELSE

!*        NO DEPTH OR CURRENT REFRACTION.
!         -------------------------------

          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *'
          WRITE (IU06,*) '* FOR CARTESIAN GRID !                   *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'

          CALL ABORT1

        ENDIF

      ENDIF

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAGS2',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAGS2
