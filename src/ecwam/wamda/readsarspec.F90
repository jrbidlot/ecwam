SUBROUTINE READSARSPEC(LNODATA)

! WRAPPER READING OF SAR DATA INPUT

! -------------------------------------------------------------------

USE YOWSARAS , ONLY : LECSARDA

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------

IMPLICIT NONE

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
  CALL READSARFCT(FCUTALL)
  CALL READECSARSPEC(LNODATA, FCUTALL, FCUTOBS)
  CALL FCUT_INTERPOLATE(FCUTOBS)
ELSE
  CALL READMFSARSPEC(LNODATA)
ENDIF

IF (ALLOCATED(FCUTOBS)) DEALLOCATE(FCUTOBS)

IF (LHOOK) CALL DR_HOOK('READSARSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE READSARSPEC
