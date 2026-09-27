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

! ----------------------------------------------------------------------
!     PRIVATE INDICES, CONSTANT OVER THE IJ VECTOR LOOP.
! ----------------------------------------------------------------------

      INTEGER(KIND=JWIM) :: KCR1
      INTEGER(KIND=JWIM) :: KCR2
      INTEGER(KIND=JWIM) :: KCR3
      INTEGER(KIND=JWIM) :: KCR4

      INTEGER(KIND=JWIM) :: KPMM
      INTEGER(KIND=JWIM) :: KPMP

      INTEGER(KIND=JWIM) :: MPMM
      INTEGER(KIND=JWIM) :: MPMP

! ----------------------------------------------------------------------
!     PRIVATE LOGICAL FLAGS, CONSTANT OVER THE IJ VECTOR LOOP.
! ----------------------------------------------------------------------

      LOGICAL :: LWLON1
      LOGICAL :: LWLON2

      LOGICAL :: LWLAT11
      LOGICAL :: LWLAT21
      LOGICAL :: LWLAT12
      LOGICAL :: LWLAT22

      LOGICAL :: LWC11
      LOGICAL :: LWC21
      LOGICAL :: LWC31
      LOGICAL :: LWC41

      LOGICAL :: LWC12
      LOGICAL :: LWC22
      LOGICAL :: LWC32
      LOGICAL :: LWC42

      LOGICAL :: LWKM
      LOGICAL :: LWKP
      LOGICAL :: LWMM
      LOGICAL :: LWMP

! ----------------------------------------------------------------------
!     LOCAL REAL VARIABLES.
! ----------------------------------------------------------------------

      REAL(KIND=JWRB)   :: ZF3
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
!*        The M and K loops are mapped to gangs.
!*        The contiguous IJ dimension is mapped to the vector level.
!*        Logical flags and neighbour indices are loaded once for each
!*        (K,M) pair and reused by all IJ vector iterations.
!*        ZF3 is private to each IJ vector iteration.

          !$acc parallel                                                      &
          !$acc& present(F1,F3,SUMWN,WLONN,KLON,LLWLONN,WLATN,KLAT,           &
          !$acc&         LLWLATN,WCORN,KCOR,KCR,LLWCORN,WKPMN,KPM,            &
          !$acc&         LLWKPMN,WMPMN,MPM,LLWMPMN)

          !$acc loop gang collapse(2)                                         &
          !$acc& private(LWLON1,LWLON2,                                       &
          !$acc&         LWLAT11,LWLAT21,LWLAT12,LWLAT22,                     &
          !$acc&         LWC11,LWC21,LWC31,LWC41,                             &
          !$acc&         LWC12,LWC22,LWC32,LWC42,                             &
          !$acc&         LWKM,LWKP,LWMM,LWMP,                                 &
          !$acc&         KCR1,KCR2,KCR3,KCR4,KPMM,KPMP,MPMM,MPMP)

          DO M = ND3S, ND3E
            DO K = 1, NANG

!*            Logical switches constant over the IJ vector loop.

              LWLON1  = LLWLONN(K,M,1)
              LWLON2  = LLWLONN(K,M,2)

              LWLAT11 = LLWLATN(K,M,1,1)
              LWLAT21 = LLWLATN(K,M,2,1)
              LWLAT12 = LLWLATN(K,M,1,2)
              LWLAT22 = LLWLATN(K,M,2,2)

              LWC11   = LLWCORN(K,M,1,1)
              LWC21   = LLWCORN(K,M,2,1)
              LWC31   = LLWCORN(K,M,3,1)
              LWC41   = LLWCORN(K,M,4,1)

              LWC12   = LLWCORN(K,M,1,2)
              LWC22   = LLWCORN(K,M,2,2)
              LWC32   = LLWCORN(K,M,3,2)
              LWC42   = LLWCORN(K,M,4,2)

              LWKM    = LLWKPMN(K,M,-1)
              LWKP    = LLWKPMN(K,M, 1)

              LWMM    = LLWMPMN(K,M,-1)
              LWMP    = LLWMPMN(K,M, 1)

!*            K-dependent neighbour indices.

              KCR1 = KCR(K,1)
              KCR2 = KCR(K,2)
              KCR3 = KCR(K,3)
              KCR4 = KCR(K,4)

              KPMM = KPM(K,-1)
              KPMP = KPM(K, 1)

!*            M-dependent neighbour indices.

              MPMM = MPM(M,-1)
              MPMP = MPM(M, 1)

              !$acc loop vector private(ZF3)
!DIR$ IVDEP
!DIR$ PREFERVECTOR
              DO IJ = KIJS, KIJL

!*              Central contribution.

                ZF3 = (1.0_JWRB-SUMWN(IJ,K,M))*F1(IJ,K,M)

!*              Longitude contributions.

                IF (LWLON1) THEN
                  ZF3 = ZF3 + WLONN(IJ,K,M,1)*F1(KLON(IJ,1),K,M)
                ENDIF

                IF (LWLON2) THEN
                  ZF3 = ZF3 + WLONN(IJ,K,M,2)*F1(KLON(IJ,2),K,M)
                ENDIF

!*              ICL=1 latitude contributions.

                IF (LWLAT11) THEN
                  ZF3 = ZF3 + WLATN(IJ,K,M,1,1)*F1(KLAT(IJ,1,1),K,M)
                ENDIF

                IF (LWLAT21) THEN
                  ZF3 = ZF3 + WLATN(IJ,K,M,2,1)*F1(KLAT(IJ,2,1),K,M)
                ENDIF

!*              ICL=1 corner contributions.

                IF (LWC11) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,1,1)*F1(KCOR(IJ,KCR1,1),K,M)
                ENDIF

                IF (LWC21) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,2,1)*F1(KCOR(IJ,KCR2,1),K,M)
                ENDIF

                IF (LWC31) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,3,1)*F1(KCOR(IJ,KCR3,1),K,M)
                ENDIF

                IF (LWC41) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,4,1)*F1(KCOR(IJ,KCR4,1),K,M)
                ENDIF

!*              ICL=2 latitude contributions.

                IF (LWLAT12) THEN
                  ZF3 = ZF3 + WLATN(IJ,K,M,1,2)*F1(KLAT(IJ,1,2),K,M)
                ENDIF

                IF (LWLAT22) THEN
                  ZF3 = ZF3 + WLATN(IJ,K,M,2,2)*F1(KLAT(IJ,2,2),K,M)
                ENDIF

!*              ICL=2 corner contributions.

                IF (LWC12) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,1,2)*F1(KCOR(IJ,KCR1,2),K,M)
                ENDIF

                IF (LWC22) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,2,2)*F1(KCOR(IJ,KCR2,2),K,M)
                ENDIF

                IF (LWC32) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,3,2)*F1(KCOR(IJ,KCR3,2),K,M)
                ENDIF

                IF (LWC42) THEN
                  ZF3 = ZF3 + WCORN(IJ,K,M,4,2)*F1(KCOR(IJ,KCR4,2),K,M)
                ENDIF

!*              IC=-1 directional refraction.

                IF (LWKM) THEN
                  ZF3 = ZF3 + WKPMN(IJ,K,M,-1)*F1(IJ,KPMM,M)
                ENDIF

!*              IC=-1 frequency refraction.

                IF (LWMM) THEN
                  ZF3 = ZF3 + WMPMN(IJ,K,M,-1)*F1(IJ,K,MPMM)
                ENDIF

!*              IC=+1 directional refraction.

                IF (LWKP) THEN
                  ZF3 = ZF3 + WKPMN(IJ,K,M,1)*F1(IJ,KPMP,M)
                ENDIF

!*              IC=+1 frequency refraction.

                IF (LWMP) THEN
                  ZF3 = ZF3 + WMPMN(IJ,K,M,1)*F1(IJ,K,MPMP)
                ENDIF

!*              Single write to F3.

                F3(IJ,K,M) = ZF3

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
