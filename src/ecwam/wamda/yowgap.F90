! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWGAP

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      REAL(KIND=JWRB), PARAMETER :: GAPVAL=0.1E-19

      INTEGER(KIND=JWIM) :: MGAP

      INTEGER(KIND=JWIM), ALLOCATABLE :: LANG(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: HANG(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: LFRE(:,:) 
      INTEGER(KIND=JWIM), ALLOCATABLE :: HFRE(:,:) 

      LOGICAL, ALLOCATABLE :: EQUIGAP(:,:,:)

!     VARIABLE   TYPE     PURPOSE                                  
!     --------   ----     -------                                  
!     *GAPVAL*   REAL     VALUE GIVEN TO SPECTRAL TO TAG THEM AS GAPS
!     *MGAP*     INTEGER  MAXIMUM NUMBER OF GAPS IN THE TRANSFORMED
!                         SPECTRA.
!     *LANG*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.  
!     *HANG*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *LFRE*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *HFRE*     INTEGER  LOWER LEFT AND UPPER RIGHT CORNER OF A GAP.
!     *EQUIGAP*  LOGICAL  FLAG INDICATING WHETHER GAPS ARE CONTAINED
!                         IN ONE ANOTHER OR OVERLAPPING.
!                         EQUIVALENCE OF GAPS, WHERE
!                         EQUIGAP(ISPEC,IGAP,GAP) = .TRUE. MEANS,
!                         THAT IN SPECTRUM NUMBER ISPEC THE GAP
!                        x  NUMBERS IGAP AND JGAP ARE EQUIVALENT,
!                         ONLY USED IN SUBROUTINE MAKEFRAMES

! ----------------------------------------------------------------------
      END MODULE YOWGAP
