! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE GAPINTERPOL(                                           &
!              DIMENSIONS
     &             NSPEC, ISPECS ,ISPECE, MGAP ,NGAP, NANG, NFRE,       &
!              INPUT ( SPEC ALSO OUTPUT)
     &             SPEC ,                                               &
     &             FX , FY , LANG, HANG , LFRE , HFRE ,                 &
     &             IU06 )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     INTERPOLATE GAPS.

!      AUTHOR.
!      -------
!      S.HASSELMANN, MPI HAMBURG, 1993.
!      J. WASZKEWITZ,MPI HAMBURG, 1993.

!     method.
!     -------
!     compute the 6 coefficients of a 2d parabola fit by minimizing
!     a cost function over all solutions gained at spectral points
!     with values > 0.
!     derivatives of the cost function with respect to the parabola
!     coefficients deliver a set of 6 equations with 6 unknowns
!     which are solved with nag routine f04jgf.

!-------------------------------------------------------------------
!     MODULE :
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWABORT,  ONLY : WAM_ABORT
      USE YOWGAP   , ONLY : GAPVAL
      USE YOWLSQ,    ONLY : LSQ

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

!     INTERFACE:
!     ----------


      INTEGER(KIND=JWIM) :: NGAP(ISPECS:ISPECE)
      INTEGER(KIND=JWIM) :: NSPEC, ISPECS, ISPECE, MGAP, NANG, NFRE
!                 DIMENSIONS.
      REAL(KIND=JWRB) :: SPEC(0:NSPEC,NANG,NFRE)
!                 2D SPECTRUM.
      REAL(KIND=JWRB) :: FX(NANG,NFRE), FY(NANG,NFRE)
!                 POLAR COORDINATES.
      INTEGER(KIND=JWIM) :: LANG(ISPECS:ISPECE,MGAP),                   &
     &                      HANG(ISPECS:ISPECE,MGAP),                   &
     &                      LFRE(ISPECS:ISPECE,MGAP),                   &
     &                      HFRE(ISPECS:ISPECE,MGAP)
!                 LOWER LEFT AND UPPER RIGHT CORNER OF GAP.
      INTEGER(KIND=JWIM) :: IU06
!                 STANDARD OUTPUT /UNIT

!     LOCAL VARIABLES:
!     ----------------

      INTEGER(KIND=JWIM) :: KERR 

      INTEGER(KIND=JWIM) :: ISPEC, IGAP, IANG, IFRE, I12, I, J
!                 LOOP INDICES.
      INTEGER(KIND=JWIM) :: LANG0(2),HANG0(2),LFRE0,HFRE0
!                 CORNERS OF FRAMES FOR TEMPORARILY USE
      REAL(KIND=JWRB) :: X, Y
!                 FREQUENCY VALUES IN CARTHESIAN COORDINATES

      REAL(KIND=JWRB), PARAMETER :: RCOND = 5.E-4_JWRB
      REAL(KIND=JWRB) :: A(6,6), B(6), T(6), ZB(6)

!-------------------------------------------------------------------

!     1.1 CHECK POSITION OF CORNERS WITH RESPECT TO 0 DEGREE.
!     ------------------------------------------------------

!     LOOP OVER SPECTRA.
      DO ISPEC = ISPECS,ISPECE 
        DO IGAP = 1,NGAP(ISPEC)
!     LOOP OVER GAPS.
          IF (LANG(ISPEC,IGAP).EQ.1    .AND.                            &
     &        HANG(ISPEC,IGAP).EQ.NANG .OR.                             &
     &        LANG(ISPEC,IGAP)-HANG(ISPEC,IGAP).EQ.2 .OR.               &
     &        LANG(ISPEC,IGAP)-HANG(ISPEC,IGAP)+NANG.EQ.2) THEN
            LANG0(1) = 1
            HANG0(1) = NANG
            LANG0(2) = 1
            HANG0(2) = 0
          ELSE
            IF (LANG(ISPEC,IGAP).GT.HANG(ISPEC,IGAP)) THEN
              LANG0(1) = 1
              HANG0(1) = HANG(ISPEC,IGAP) + 1
              LANG0(2) = LANG(ISPEC,IGAP) - 1
              HANG0(2) = NANG
            ELSE
              LANG0(1) = MAX( LANG(ISPEC,IGAP)-1 , 1)
              HANG0(1) = MIN( HANG(ISPEC,IGAP)+1 , NANG)
              IF (LANG(ISPEC,IGAP).EQ.1) THEN
                LANG0(2) = NANG
                HANG0(2) = NANG
              ELSEIF (HANG(ISPEC,IGAP).EQ.NANG) THEN
                LANG0(2) = 1
                HANG0(2) = 1
              ELSE
                LANG0(2) = 1
                HANG0(2) = 0
              END IF
            END IF
          END IF
          LFRE0 = MAX(1,LFRE(ISPEC,IGAP)-1)
          HFRE0 = MIN(NFRE,HFRE(ISPEC,IGAP)+1)

!     1.2 COMPUTE COEFFICIENT MATRIX FOR LINEAR SYSTEM OF 6 
!         EQUATIONS.
!     ---------------------------------------------------------

          DO I = 1,6
            DO J = 1,6
              A(I,J) = 0.0_JWRB
            END DO
            B(I) = 0.0_JWRB
          END DO
          DO I12 = 1,2
            DO IANG = LANG0(I12),HANG0(I12)
              DO IFRE = LFRE0,HFRE0
                IF (SPEC(ISPEC,IANG,IFRE).NE.0.0_JWRB)THEN
                  X = FX(IANG,IFRE)
                  Y = FY(IANG,IFRE)
                  T(1) = 1.0_JWRB
                  T(2) = X
                  T(3) = Y
                  T(4) = X*X
                  T(5) = Y*Y
                  T(6) = X*Y
                  DO I = 1,6
                    DO J= 1,6
                      A(I,J) = A(I,J) + T(I) * T(J)
                    END DO
                    B(I) = B(I) + T(I) * SPEC(ISPEC,IANG,IFRE)
                  END DO
                END IF
              END DO
            END DO
          END DO

!      1.3 SOLVE SYSTEM OF EQUATIONS WITH RESPECT TO COEFFICIENTS 
!          FOR 2D PARABOLA FIT.
!      ------------------------------------------------------------


          CALL LSQ(A, B, RCOND, ZB, KERR)

          IF (KERR == 0) THEN
            ZB(1) = MAX(ZB(1),0.0_JWRB)
            DO I12 = 1,2
              DO IANG = LANG0(I12),HANG0(I12)
                DO IFRE = LFRE0,HFRE0
                  IF (SPEC(ISPEC,IANG,IFRE).EQ.GAPVAL) THEN
                    X = FX(IANG,IFRE)
                    Y = FY(IANG,IFRE)
                    SPEC(ISPEC,IANG,IFRE) = MAX(0.0_JWRB ,              &
     &                ZB(1)+ZB(2)*X+ZB(3)*Y+ZB(4)*X*X+ZB(5)*Y*Y+ZB(6)*X*Y)
                  END IF
                END DO
              END DO
            END DO
          END IF

!       END OF LOOPS OVER NUMBER OF GAPS.

        END DO
!    
!     END OF LOOP OVER NUMBER OF SPECTRA.

      END DO


      END SUBROUTINE GAPINTERPOL
