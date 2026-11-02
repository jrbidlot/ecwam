SUBROUTINE READMFSARSPEC (LNODATA)

!     READ OBSERVED SPECTRA FILE AS INPUT BY METEOFRANCE
!
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : FLMINFR,FLMIN
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWSARAS , ONLY : NSPEC, SPEC, LAT, LONG, XKLOW, XKUP, &
     &                      XKCUTOFF, USAR, RA, CUT, NVI
      USE YOWSTAT  , ONLY : CDTPRO  ,LSARINV
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: LNODATA


      INTEGER(KIND=JWIM) :: INUNIT
      INTEGER(KIND=JWIM) :: ISPEC, ISS, IC
      INTEGER(KIND=JWIM) :: IANG, IFRE

      REAL(KIND=JWRB), ALLOCATABLE :: ANG_TR(:)

      CHARACTER(LEN=12) :: DATESAR
      CHARACTER(LEN=15) :: YONAME
      CHARACTER(LEN=50) :: FORM_SPEC


      NSPEC = 15000


      ALLOCATE(SPEC(0:NSPEC,NANG,NFRE,2))
      ALLOCATE(LONG(NSPEC,2))
      ALLOCATE(LAT(NSPEC,2))
      ALLOCATE(XKLOW(NSPEC))
      ALLOCATE(XKUP(NSPEC))
      ALLOCATE(XKCUTOFF(NSPEC))
      IF (.NOT.ALLOCATED(USAR)) ALLOCATE(USAR(NSPEC))
      IF (.NOT.ALLOCATED(RA)) ALLOCATE(RA(NSPEC))
      IF (.NOT.ALLOCATED(CUT)) ALLOCATE(CUT(NSPEC))
      IF (.NOT.ALLOCATED(NVI)) ALLOCATE(NVI(NSPEC))
      IF (.NOT.ALLOCATED(ANG_TR)) ALLOCATE(ANG_TR(NSPEC))

      DO IC=1,2
        DO IANG=1,NANG
          DO IFRE=1,NFRE
            DO ISS=0,NSPEC
              SPEC(ISS,IANG,IFRE,IC) = 0.25_JWRB*FLMINFR(1,IFRE)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      INUNIT = 60
      YONAME(1:3) = 'SAR'
      YONAME(4:15) = CDTPRO
      OPEN (UNIT=INUNIT, FILE=YONAME)

      ISPEC = 1

      LONG(:,:) = 0.0_JWRB
      LAT(:,:) = 0.0_JWRB
      XKLOW(:) = 0.0_JWRB
      XKUP(:) = 0.0_JWRB
      XKCUTOFF(:) = 0.0_JWRB

      LNODATA=.FALSE.

      WRITE(FORM_SPEC,'("(",i,"F8.3 )")' ) NANG

      DO WHILE(.TRUE.)
        READ(INUNIT,1500,END=20,ERR=20) DATESAR, LONG(ISPEC,2), LAT(ISPEC,2),     &
     &    USAR(ISPEC), RA(ISPEC), CUT(ISPEC), NVI(ISPEC), ANG_TR(ISPEC)

        DO IFRE=1,NFRE
          READ(INUNIT,FORM_SPEC,END=20) (SPEC(ISPEC,IANG,IFRE,2),IANG=1,NANG)
        END DO

        ISPEC=ISPEC+1
      ENDDO

   20 CONTINUE

 1500 FORMAT(A12,7F10.3)

!! JB : this is very bad, all your arrays have been declare with NSPEC=15000 as in line 28
!!! might want to adapt to using resize_sarinput !!!!
      NSPEC=ISPEC-1
      WRITE(*,*) 'test ispec ', NSPEC,NFRE

      IF (NSPEC == 0) LNODATA=.TRUE.

      XKCUTOFF(:)=CUT(:)

!     CALCULATE MEAN PARAMETERS (not any more but still check SPEC for negative value (needed??)

      DO ISPEC=1,NSPEC
        DO IFRE=1,NFRE
          DO IANG=1,NANG
            IF (SPEC(ISPEC,IANG,IFRE,2) < 0.0_JWRB) SPEC(ISPEC,IANG,IFRE,2)=FLMIN
          END DO
        END DO
      END DO


END SUBROUTINE READMFSARSPEC
