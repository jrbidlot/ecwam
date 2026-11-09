! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!------------------------------------------------------------------------------
!  WRITEPRTOUT
!
!  Writes all of DEV()'s INPUT arguments to an unformatted file
!  named "PRT"//CDTPRO, following the conventions used in WRITESASOUT().
!
!  File structure (unformatted records):
!    1) CHARACTER*12: CDTPRO(1:12)
!    2) INTS: NSPEC, MPART
!    3) ARRAYS (shapes as passed to DEV):
!       - NPARTW(NSPEC), NPARTS(NSPEC)
!       - XKLOW(NSPEC), XKUP(NSPEC), XKCUTOFF(NSPEC)
!       - XMEANEW(NSPEC,MPART), XMEANANGW(NSPEC,MPART), XMEANFREW(NSPEC,MPART)
!       - PKFREW(NSPEC,MPART), PARTINFOW(NSPEC,MPART)
!       - XMEANES(NSPEC,MPART), XMEANANGS(NSPEC,MPART), XMEANFRES(NSPEC,MPART)
!       - PKFRES(NSPEC,MPART), PARTINFOS(NSPEC,MPART)
!       - NCORRELTAB(NSPEC,MPART)
!       - CONTROL (LOGICAL)
!
!  Notes:
!   - Uses same open pattern as WRITESASOUT: IWAM_GET_UNIT(IU06, name, 's','u',0)
!   - Writes full-sized arrays (NSPEC, MPART) even if some rows are partially
!     used (e.g. only NPARTW(JSPEC) partitions are meaningful). This keeps the
!     file layout simple and deterministic for later reading.
!------------------------------------------------------------------------------

SUBROUTINE WRITEPRTOUT( &
!   DIMENSIONS
  &  NSPEC , MPART,        &
!   INPUTS (all of DEV's)
  &  NCORRELTAB, XKLOW, XKUP, XKCUTOFF,               &
  &  XMEANEW , XMEANANGW, XMEANFREW , PKFREW, PARTINFOW, &
  &  XMEANES , XMEANANGS, XMEANFRES , PKFRES, PARTINFOS, &
  &  NPARTW, NPARTS,                                  &
  &  CONTROL, COUT, IREAD)

!------------------------------------------------------------------------------

  USE PARKIND_WAVE  , ONLY : JWIM, JWRB

  USE YOWSTAT       , ONLY : CDTPRO
  USE YOWTEST       , ONLY : IU06
  USE YOWMPP        , ONLY : IRANK

!------------------------------------------------------------------------------

  IMPLICIT NONE

#include "iwam_get_unit.intfb.h"

  !-------------------------
  ! Arguments
  !-------------------------
  INTEGER(KIND=JWIM), INTENT(IN) :: NSPEC, MPART

  INTEGER(KIND=JWIM), INTENT(IN) :: NPARTW(NSPEC), NPARTS(NSPEC)

  REAL   (KIND=JWRB), INTENT(IN) :: XKLOW   (NSPEC)
  REAL   (KIND=JWRB), INTENT(IN) :: XKUP    (NSPEC)
  REAL   (KIND=JWRB), INTENT(IN) :: XKCUTOFF(NSPEC)

  REAL   (KIND=JWRB), INTENT(IN) :: XMEANEW   (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: XMEANANGW (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: XMEANFREW (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: PKFREW    (NSPEC,MPART)
  INTEGER(KIND=JWIM), INTENT(IN) :: PARTINFOW (NSPEC,MPART)

  REAL   (KIND=JWRB), INTENT(IN) :: XMEANES   (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: XMEANANGS (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: XMEANFRES (NSPEC,MPART)
  REAL   (KIND=JWRB), INTENT(IN) :: PKFRES    (NSPEC,MPART)
  INTEGER(KIND=JWIM), INTENT(IN) :: PARTINFOS (NSPEC,MPART)

  INTEGER(KIND=JWIM), INTENT(IN) :: NCORRELTAB(NSPEC,MPART)

  LOGICAL, INTENT(IN) :: CONTROL
  CHARACTER(2), INTENT(IN)  :: COUT
  INTEGER(KIND=JWIM), INTENT(IN) :: IREAD

  !-------------------------
  ! Locals
  !-------------------------
  CHARACTER(15)      :: FNAME
  INTEGER(KIND=JWIM) :: IUPRT
  CHARACTER(3)       :: PCOUT

  !-------------------------
  ! Begin
  !-------------------------

  IF (IRANK == IREAD) THEN

    PCOUT = 'P' // COUT

    WRITE(IU06,*) 'ATTEMPTING TO WRITE PRT FILE'
    WRITE(IU06,*) 'PREFIX = ', PCOUT
    WRITE(IU06,*) 'CDTPRO = ', CDTPRO
    WRITE(IU06,*) 'NSPEC  = ', NSPEC
    WRITE(IU06,*) 'MPART  = ', MPART

    ! Build filename and open (same convention as WRITESASOUT)
    WRITE(FNAME,'(A3,A12)') PCOUT, CDTPRO
    IUPRT = IWAM_GET_UNIT(IU06, FNAME, 's', 'u', 0, 'READWRITE')

    WRITE(IU06,*) 'PRT FILENAME = ', FNAME

    !-------------------------
    ! Write header / sizes
    !-------------------------
    WRITE(IUPRT) CDTPRO(1:12)
    WRITE(IUPRT) NSPEC, MPART

    !-------------------------
    ! Write all DEV inputs
    !-------------------------
    ! Dimensions per spectrum
    WRITE(IUPRT) NPARTW
    WRITE(IUPRT) NPARTS

    ! SAR wavenumber limits
    WRITE(IUPRT) XKLOW
    WRITE(IUPRT) XKUP
    WRITE(IUPRT) XKCUTOFF

    ! WAM partition means / peaks / info
    WRITE(IUPRT) XMEANEW
    WRITE(IUPRT) XMEANANGW
    WRITE(IUPRT) XMEANFREW
    WRITE(IUPRT) PKFREW
    WRITE(IUPRT) PARTINFOW

    ! SAR partition means / peaks / info
    WRITE(IUPRT) XMEANES
    WRITE(IUPRT) XMEANANGS
    WRITE(IUPRT) XMEANFRES
    WRITE(IUPRT) PKFRES
    WRITE(IUPRT) PARTINFOS

    ! Cross-assignment table
    WRITE(IUPRT) NCORRELTAB

    ! Control flag (for traceability)
    WRITE(IUPRT) CONTROL

    CLOSE(IUPRT)

    WRITE(IU06,*) 'PRT FILE WRITE COMPLETE'

  ENDIF

END SUBROUTINE WRITEPRTOUT
