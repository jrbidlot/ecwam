! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITEDIAOUT( ispec, iparts, njs, jpsiw,                 &
     &                        obs_e,  mod_e,                             &
     &                        obs_kx, mod_kx, obs_ky, mod_ky,            &
     &                        cond_e, cond_kx, cond_ky,                  &
     &                        w_e, w_kx, w_ky,                           &
     &                        idx_used,                                  &
     &                        COUT )

!------------------------------------------------------------------
! PURPOSE:
!   Write per-partition OI diagnostics: identifiers, error settings,
!   condition numbers, index list, and reduced weight vectors for
!   Energy and for Kx/Ky separately.
!------------------------------------------------------------------

      USE PARKIND_WAVE , ONLY : JWIM, JWRB

      USE YOWSTAT      , ONLY : CDTPRO
      USE YOWTEST      , ONLY : IU06

!------------------------------------------------------------------

      IMPLICIT NONE

#include "iwam_get_unit.intfb.h"

!-------------------- args --------------------
      INTEGER(JWIM), INTENT(IN) :: ispec, iparts, njs, jpsiw
      REAL(JWRB)   , INTENT(IN) :: obs_e,  mod_e
      REAL(JWRB)   , INTENT(IN) :: obs_kx, mod_kx, obs_ky, mod_ky
      REAL(JWRB)   , INTENT(IN) :: cond_e, cond_kx, cond_ky
      REAL(JWRB)   , INTENT(IN) :: w_e (jpsiw)     ! weights for E
      REAL(JWRB)   , INTENT(IN) :: w_kx(jpsiw)     ! weights for Kx
      REAL(JWRB)   , INTENT(IN) :: w_ky(jpsiw)     ! weights for Ky
      INTEGER(JWIM), INTENT(IN) :: idx_used(njs)   ! positions used (1..jpsiw)
      CHARACTER(2) , INTENT(IN) :: COUT

!-------------------- locals ------------------
      CHARACTER(15) :: FNAME
      CHARACTER(3)  :: DCOUT
      INTEGER       :: IUDIA, j
      REAL(JWRB)    :: w_e_used(njs), w_kx_used(njs), w_ky_used(njs)

!-------------------- body --------------------
      DCOUT = 'D' // COUT

! Reduce to used entries (preserve order given in idx_used)
      DO j = 1, njs
        IF (idx_used(j) < 1 .OR. idx_used(j) > jpsiw) THEN
          WRITE(IU06,*) 'WRITEDIAOUT: idx_used out of range:', idx_used(j), ' jpsiw=', jpsiw
          CYCLE
        END IF
        w_e_used (j) = w_e ( idx_used(j) )
        w_kx_used(j) = w_kx( idx_used(j) )
        w_ky_used(j) = w_ky( idx_used(j) )
      END DO

! Name like DCOUT//CDTPRO to mirror your existing convention
      WRITE(FNAME,'(A3,A12)') DCOUT, CDTPRO
      IUDIA = IWAM_GET_UNIT(IU06, FNAME, 's', 'u', 0, 'READWRITE')

! --- file layout (compact, per-call record group) ---
! 1) identifiers + sizes
      WRITE(IUDIA) ispec, iparts, njs, jpsiw
! 2) scalar settings (errors + condition numbers)
      WRITE(IUDIA) obs_e, mod_e, obs_kx, mod_kx, obs_ky, mod_ky,        &
     &              cond_e, cond_kx, cond_ky
! 3) index list (as written) and the reduced weights
      WRITE(IUDIA) idx_used
      WRITE(IUDIA) w_e_used
      WRITE(IUDIA) w_kx_used
      WRITE(IUDIA) w_ky_used

      CLOSE(IUDIA)

      END SUBROUTINE WRITEDIAOUT
