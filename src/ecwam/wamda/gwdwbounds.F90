! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GWDWBOUNDS(IDATAWL, CDATEAN, CBEGINDT, CENDDT)

!--------------------------------------------------------------------

!**** *GWDWBOUNDS* -  SETS THE DATE BOUNDS ON THE ANALYSIS WINDOW.

!--------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK

! -------------------------------------------------------------------

IMPLICIT NONE

#include "incdate.intfb.h"

INTEGER(KIND=JWIM), INTENT(IN) :: IDATAWL  ! ANALYSIS WINDOW LENGTH IN SECONDS.
CHARACTER(LEN=14), INTENT(IN)  :: CDATEAN  ! ANALYSIS TIME.

CHARACTER(LEN=14), INTENT(OUT) :: CBEGINDT ! BEGINNING DATE OF THE TIME WINDOW DURING WHICH DATA ARE USED FOR ANALYSIS.
CHARACTER(LEN=14), INTENT(OUT) :: CENDDT   ! END DATE OF THE TIME WINDOW DURING WHICH DATA ARE SAVE FOR ANALYSIS.

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GWDWBOUNDS',0,ZHOOK_HANDLE)

IF (IDATAWL <= 0 ) THEN
  CBEGINDT = CDATEAN
  CENDDT = CBEGINDT
ELSEIF (IDATAWL <= 10800) THEN
! FOR WINDOW LENGTH OF 3 HOURS OR LESS, THE ANALYSIS TIME IS AT THE END OF THE WINDOW.
  CBEGINDT = CDATEAN
  CALL INCDATE (CBEGINDT, -IDATAWL)
  CENDDT = CBEGINDT
  CALL INCDATE (CENDDT, IDATAWL)
  ! ADD 1 SECOND TO AVOID OVERLAPPING WINDOWS
  CALL INCDATE (CBEGINDT, 1)
ELSE
! OTHERWISE, THE ANALYSIS TIME IS IN THE MIDDLE OF THE WINDOW.
  CBEGINDT = CDATEAN
  CALL INCDATE (CBEGINDT, -IDATAWL/2)
  CENDDT = CBEGINDT
  CALL INCDATE (CENDDT, IDATAWL)
ENDIF

IF (LHOOK) CALL DR_HOOK('GWDWBOUNDS',1,ZHOOK_HANDLE)

END SUBROUTINE GWDWBOUNDS
