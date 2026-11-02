SUBROUTINE READSARFCT(FCUTALL)
  !-----------------------------------------------------------------------
  ! This subroutine imports FCUTALL spectra from an open FCT file.
  !-----------------------------------------------------------------------

  USE PARKIND_WAVE  , ONLY : JWIM, JWRB

  USE YOWSARAS , ONLY : NSAR              ! Number of SAR observations
  USE YOWSTAT  , ONLY : CDTPRO
  USE YOWTEST  , ONLY : IU06              ! Output log unit
  USE YOWDES   , ONLY : JPTSAR            ! Number of directions

! ----------------------------------------------------------------------
  IMPLICIT NONE

#include "abort1.intfb.h"
#include "iwam_get_unit.intfb.h"

  REAL, ALLOCATABLE, INTENT(OUT) :: FCUTALL(:,:)

  INTEGER :: IUFCT
  INTEGER :: ISAR, IDIR
  CHARACTER(LEN=3)  :: CFCT
  CHARACTER(LEN=12) :: CL_DATE
  CHARACTER(LEN=15) :: FCTNAME
  LOGICAL :: LLFCTEXIST

! ----------------------------------------------------------------------

  CFCT='FCT' ! FCT file prefix

  ! Check for FCT file (it should have been produced by the pre-processing) 
  WRITE(FCTNAME,'(A3,A12)') CFCT, CDTPRO
  INQUIRE (FILE=FCTNAME, EXIST=LLFCTEXIST)
  IF (LLFCTEXIST) THEN
    IUFCT=IWAM_GET_UNIT(IU06, FCTNAME, 's', 'u', 0, 'READ')
    WRITE(IU06,*) '  READSARFCT: ',FCTNAME, ' FOUND'
  ELSE
    WRITE(IU06,*) '  READSARFCT: ',FCTNAME, ' NOT FOUND'
    CALL FLUSH(IU06)
    CALL ABORT1
  ENDIF
  WRITE(IU06,*) ''

  ! Read header line with date string
  READ(IUFCT, END=2100, ERR=2000) CL_DATE
  WRITE(IU06,*) '  READSARFCT: CL_DATE FROM FCT = ', CL_DATE

  ! Read FCUT spectra for each SAR observation
  ALLOCATE(FCUTALL(JPTSAR, NSAR))
  FCUTALL(:,:) = -999._JWRB
  DO ISAR = 1, NSAR
    READ(IUFCT, ERR=2000, END=2100) (FCUTALL(IDIR, ISAR), IDIR = 1, JPTSAR)
  END DO

  ! Close file after reading
  CLOSE(IUFCT)
  RETURN

! Error handling
2000 CONTINUE
  WRITE(IU06,*) '  READSARFCT: Error reading FCT file.'
  CALL ABORT1

2100 CONTINUE
  WRITE(IU06,*) '  READSARFCT: Unexpected end of FCT file.'
  WRITE(IU06,*) '  READSARFCT: NSAR = ', NSAR
  WRITE(IU06,*) '  READSARFCT: ISAR = ', ISAR
  WRITE(IU06,*) '  READSARFCT: JPTSAR = ', JPTSAR
  CALL ABORT1

END SUBROUTINE READSARFCT
