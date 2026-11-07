      SUBROUTINE READMFSARSPEC (LNODATA)

!     READ OBSERVED SPECTRA FILE AS INPUT BY METEOFRANCE
!
!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWPCONS , ONLY : EPSMIN
      USE YOWSARAS , ONLY : NSPEC, SPEC, LAT, LONG, XKLOW, XKUP, XKCUTOFF
      USE YOWSTAT  , ONLY : CDTPRO  ,LSARINV
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

      LOGICAL, INTENT(OUT) :: LNODATA


      INTEGER(KIND=JWIM) :: INUNIT
      INTEGER(KIND=JWIM) :: ISPEC, ISS, IC
      INTEGER(KIND=JWIM) :: IANG, IFRE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE ! Dr HOOK

!    These following variables are read in but appears to not not be used !
      REAL(KIND=JWRB), ALLOCATABLE :: USAR(:)   !! WIND SPEED FROM SAR
      REAL(KIND=JWRB), ALLOCATABLE :: RA(:)     !! RATIO OF SIGNAL TO NOISE
      REAL(KIND=JWRB), ALLOCATABLE :: NVI(:)    !! NORMALIZED VARIANCE OF IMAGETTES
      REAL(KIND=JWRB), ALLOCATABLE :: ANG_TR(:) !! ??

      CHARACTER(LEN=12) :: DATESAR
      CHARACTER(LEN=15) :: YONAME
      CHARACTER(LEN=50) :: FORM_SPEC

      LOGICAL :: LLOPENED, LLTEST

!------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READMFSARSPEC',0,ZHOOK_HANDLE)


      NSPEC = 15000

!     note that in parallel run spec, long and lat will be allocated
!     for the other PE's in mpbcastsarin.
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
              SPEC(ISS,IANG,IFRE,IC) = 0.25_JWRB*FLMIN
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      LONG(:,:) = 0.0_JWRB
      LAT(:,:) = 0.0_JWRB
      XKLOW(:) = 0.0_JWRB
      XKUP(:) = 0.0_JWRB
      XKCUTOFF(:) = 0.0_JWRB

      LNODATA=.FALSE.

      UNIT: DO INUNIT=34,99,1
        INQUIRE ( UNIT=INUNIT, OPENED=LLOPENED)
        IF ( .NOT. LLOPENED) THEN
          WRITE(IU06,*) '  READMLSARSPEC: TAKES UNIT ',INUNIT,' AS SAR DATA INPUT UNIT'
          EXIT UNIT
        ENDIF
      ENDDO UNIT
      IF(INUNIT.GT.99) THEN
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) ' +                                        +'
        WRITE(IU06,*) ' + ERROR IN --READMLSARSPEC--             +'
        WRITE(IU06,*) ' + =========================              +'
        WRITE(IU06,*) ' + RUN OUT OF UNIT TO OPEN SAR DATA INPUT +'
        WRITE(IU06,*) ' + PROGRAM WILL ABORT                     +'
        WRITE(IU06,*) ' +                                        +'
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF


      YONAME(1:3) = 'SAR'
      YONAME(4:15) = CDTPRO
      OPEN (UNIT=INUNIT, FILE=YONAME)

      ISPEC = 1

      WRITE(FORM_SPEC,'("(",i,"F8.3 )")' ) NANG

      DO WHILE(.TRUE.)
        READ(INUNIT,1500,END=20,ERR=20) DATESAR, LONG(ISPEC,2), LAT(ISPEC,2),     &
     &    USAR(ISPEC), RA(ISPEC),XKCUTOFF(ISPEC), NVI(ISPEC), ANG_TR(ISPEC)

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

to finish

      IF (LHOOK) CALL DR_HOOK('READMFSARSPEC',1,ZHOOK_HANDLE)

      END SUBROUTINE READMFSARSPEC
