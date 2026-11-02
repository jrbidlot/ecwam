SUBROUTINE FCUT_INTERPOLATE(FCUTOBS)
!-----------------------------------------------------------------------
! LINEARLY INTERPOLATE FCUTOBS(JPTSAR,NSPEC) FROM THE SAR ANGULAR GRID
! (JPTSAR BINS) TO THE MODEL ANGULAR GRID (NANG BINS), WITHOUT ROTATION.
! RESULTS GO INTO MODULE ARRAYS FCUT(NANG,NSPEC) AND M10(NANG,NSPEC).
! M10 is an interval index; defaults to 1 and is clamped to [1, NFRE-1].
!-----------------------------------------------------------------------

  USE YOWSARAS , ONLY : NSPEC               ! # spectra
  USE YOWDES   , ONLY : JPTSAR              ! # SAR angle bins
  USE YOWFCT   , ONLY : FCUT, M10
  USE YOWFRED  , ONLY : FR                  ! frequency grid (size NFRE)
  USE YOWPARAM , ONLY : NFRE, NANG          ! # freqs, # model angle bins
  IMPLICIT NONE

  ! Arguments
  REAL   , INTENT(IN) :: FCUTOBS(JPTSAR,NSPEC)

  ! Locals
  INTEGER :: ISPEC, IANG
  INTEGER :: I0, I1, IINT, ILOW
  REAL    :: IDX, FRAC, FVAL

  ! (Re)allocate module arrays safely
  IF (ALLOCATED(FCUT)) THEN
     IF (SIZE(FCUT,1) /= NANG .OR. SIZE(FCUT,2) /= NSPEC) DEALLOCATE(FCUT)
  END IF
  IF (ALLOCATED(M10)) THEN
     IF (SIZE(M10,1) /= NANG .OR. SIZE(M10,2) /= NSPEC) DEALLOCATE(M10)
  END IF
  IF (.NOT. ALLOCATED(FCUT)) ALLOCATE(FCUT(NANG,NSPEC))
  IF (.NOT. ALLOCATED(M10))  ALLOCATE(M10(NANG,NSPEC))

  ! Default M10 = 1 everywhere
  M10(:,:) = 1

  DO ISPEC = 1, NSPEC
    DO IANG = 1, NANG

      ! Fractional SAR-bin index corresponding to this model angle
      IDX  = REAL(IANG-1) * REAL(JPTSAR) / REAL(NANG)

      ! Neighbour SAR bins with periodic wrap
      IINT = INT(IDX)                    ! floor since IDX >= 0
      I0   = IINT + 1
      IF (I0 > JPTSAR) I0 = JPTSAR
      I1   = I0 + 1
      IF (I1 > JPTSAR) I1 = 1

      FRAC = IDX - REAL(IINT)

      ! Linear interpolation between SAR bins I0 and I1
      FVAL = (1.0 - FRAC) * FCUTOBS(I0,ISPEC) + FRAC * FCUTOBS(I1,ISPEC)
      FCUT(IANG,ISPEC) = FVAL

      ! Find largest i with FR(i) <= FVAL; clamp to an **interval** [1, NFRE-1]
      ILOW = 1
      IF (NFRE >= 2) THEN
        DO WHILE (ILOW < NFRE .AND. FR(ILOW+1) <= FVAL)
          ILOW = ILOW + 1
        END DO
        ! Force interval index so M10+1 is always valid downstream
        ILOW = MAX(1, MIN(ILOW, NFRE-1))
      ELSE
        ILOW = 1
      END IF

      M10(IANG,ISPEC) = ILOW
    END DO
  END DO

END SUBROUTINE FCUT_INTERPOLATE
