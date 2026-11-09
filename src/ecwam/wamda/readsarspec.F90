! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE READSARSPEC(LNODATA)

! WRAPPER READING OF SAR DATA INPUT

! -------------------------------------------------------------------

USE YOWSARAS , ONLY : LECSARDA
USE YOWSTAT  , ONLY : LSWIMAS
USE YOWTEST  , ONLY : IU06

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------

IMPLICIT NONE

#include "abort1.intfb.h"
#include "fcut_interpolate.intfb.h"
#include "readecsarspec.intfb.h"
#include "readmfsarspec.intfb.h"
#include "readsarfct.intfb.h"


LOGICAL, INTENT(OUT) :: LNODATA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE ! Dr HOOK

REAL, ALLOCATABLE :: FCUTALL(:,:), FCUTOBS(:,:)

!--------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('READSARSPEC',0,ZHOOK_HANDLE)

LNODATA = .TRUE.

IF (LECSARDA) THEN
  ! ECMWF DATA ACCESS
  IF (LSWIMAS) THEN
    WRITE(IU06,*) ' ++++++++++++++++++++++++++++'
    WRITE(IU06,*) ' + ERROR IN --READSARSPEC-- +'
    WRITE(IU06,*) ' + NOT READY FOR LSWIMAS    +'
    WRITE(IU06,*) ' + PROGRAM WILL ABORT       +'
    WRITE(IU06,*) ' ++++++++++++++++++++++++++++'
    CALL ABORT1
  ENDIF

  CALL READSARFCT(FCUTALL)
  CALL READECSARSPEC(LNODATA, FCUTALL, FCUTOBS)
  CALL FCUT_INTERPOLATE(FCUTOBS)
ELSE
  ! METEO FRANCE DATA ACCESS
  CALL READMFSARSPEC(LNODATA)
ENDIF

IF (ALLOCATED(FCUTOBS)) DEALLOCATE(FCUTOBS)

IF (LHOOK) CALL DR_HOOK('READSARSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE READSARSPEC
