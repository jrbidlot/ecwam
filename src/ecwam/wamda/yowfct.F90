! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE YOWFCT

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED,  ONLY : FR
      USE YOWPARAM, ONLY : NFRE
      USE YOWSTAT , ONLY : LSARAS   , LSWIMAS
      USE YOWTEST , ONLY : IU06

      IMPLICIT NONE

#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: M10CONST
      INTEGER(KIND=JWIM), ALLOCATABLE :: M10(:,:)

      REAL(KIND=JWRB) :: FCUTCONST
      REAL(KIND=JWRB), ALLOCATABLE    :: FCUT(:,:)

!     VARIABLE   TYPE     PURPOSE
!     --------   ----     -------
!     *M10CONST* INT.     INDEX OF THE CONSTANT CUT OF FREQUENCY
!     *M10*      INT.     MODEL INDEX OF THE CUT OF FREQUENCY
!     *FCUTCONST*REAL     CONSTANT CUT-OFF FREQUENCY
!     *FCUT*     REAL     CUT-OFF FREQUENCY INTERPOLATED ONTO MODEL GRID
! ----------------------------------------------------------------------

      ! LOOP THROUGH FR AND FIND THE INDEX WHICH CORRESPONDS TO FCUTCONST

      CONTAINS

      SUBROUTINE INIT_YOWFCT()
        INTEGER(KIND=JWIM) :: M

        IF (LSARAS) THEN
          FCUTCONST = 0.121_JWRB ! 0.121 is the MF value
        ELSEIF (LSWIMAS) THEN
          FCUTCONST = 0.151_JWRB ! 0.151 is the MF value
        ELSE
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' + ERROR IN --YOWFCT-- +'
          WRITE(IU06,*) ' + UKNOWN CORRELATION LENGTH  +'
          WRITE(IU06,*) ' + PROGRAM WILL ABORT         +'
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

        M10CONST = NFRE

        DO M = NFRE, 1, -1
          IF (FR(M) < FCUTCONST) THEN
            M10CONST = M
            EXIT
          END IF
        END DO

      END SUBROUTINE INIT_YOWFCT

END MODULE YOWFCT
