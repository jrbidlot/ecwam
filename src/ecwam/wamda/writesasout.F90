! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WRITESASOUT(COUT)

  USE PARKIND_WAVE , ONLY : JWIM, JWRB

  USE YOWSTAT       , ONLY : CDTPRO
  USE YOWTEST       , ONLY : IU06
  USE YOWSARAS      , ONLY : NSPEC, LONG, LAT, SPEC
  USE YOWFCT        , ONLY : FCUT, M10
  USE YOWPARAM      , ONLY : NANG, NFRE

!------------------------------------------------------------------

  IMPLICIT NONE

#include "iwam_get_unit.intfb.h"

!------------------------------------------------------------------

  CHARACTER(2), INTENT(IN) :: COUT

  CHARACTER(3)  :: CSAS
  CHARACTER(15) :: FSASNAME
  INTEGER(JWIM) :: IUSAS, ISPEC

  !-------------------------------------------------------
  ! 1) Write SAS* output (diagnostic)
  !-------------------------------------------------------

  CSAS = 'S' // COUT

  WRITE(IU06,*) 'ATTEMPTING TO WRITE SAS FILE'
  WRITE(IU06,*) 'CSAS   = ', CSAS
  WRITE(IU06,*) 'CDTPRO = ', CDTPRO
  WRITE(IU06,*) 'NSPEC  = ', NSPEC

  ! Build filename and open
  WRITE(FSASNAME,'(A3,A12)') CSAS, CDTPRO
  IUSAS = IWAM_GET_UNIT(IU06, FSASNAME, 's', 'u', 0, 'READWRITE')
  WRITE(IU06,*) 'SAS FILENAME = ', FSASNAME

  ! Observations & co-located spectra
  WRITE(IU06,*) 'WRITING OBS AND COLLOCATED SPECTRA'

  WRITE(IUSAS) CDTPRO(1:12)
  WRITE(IUSAS) NANG, NFRE
  WRITE(IUSAS) NSPEC
  DO ISPEC = 1, NSPEC
    WRITE(IUSAS) ISPEC
    WRITE(IUSAS) LAT (ISPEC,1), LONG(ISPEC,1)
    WRITE(IUSAS) LAT (ISPEC,2), LONG(ISPEC,2)
    WRITE(IUSAS) FCUT(:,ISPEC)
    WRITE(IUSAS) M10 (:,ISPEC)
    WRITE(IUSAS) SPEC(ISPEC,:,:,1), SPEC(ISPEC,:,:,2)
  END DO

  CLOSE(IUSAS)

END SUBROUTINE WRITESASOUT
