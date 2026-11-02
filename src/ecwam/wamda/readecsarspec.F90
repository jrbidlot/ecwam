SUBROUTINE READECSARSPEC(LNODATA)

! WRAPPER FOR ECMWF READING OF SAR DATA INPUT

! -------------------------------------------------------------------

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------

IMPLICIT NONE

#include "fcut_interpolate.intfb.h"
#include "readsarfct.intfb.h"
#include "readsarspec.intfb.h"

LOGICAL, INTENT(OUT) :: LNODATA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE ! Dr HOOK

REAL, ALLOCATABLE :: FCUTALL(:,:), FCUTOBS(:,:)

!--------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('READECSARSPEC',0,ZHOOK_HANDLE)

LNODATA = .TRUE.

CALL READSARFCT(FCUTALL)

CALL READSARSPEC(LNODATA, FCUTALL, FCUTOBS)

CALL FCUT_INTERPOLATE(FCUTOBS)

IF (ALLOCATED(FCUTOBS)) DEALLOCATE(FCUTOBS)

IF (LHOOK) CALL DR_HOOK('READECSARSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE READECSARSPEC
