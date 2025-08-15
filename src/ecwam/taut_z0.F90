! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TAUT_Z0(KIJS, KIJL, IUSFG,          &
&                  HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC, CICOVER, &
&                  USTAR, Z0, Z0B, CHRNCK)

! ----------------------------------------------------------------------

!**** *TAUT_Z0* - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.


!**   INTERFACE.
!     ----------

!       *CALL* *TAUT_Z0(KIJS, KIJL, IUSFG, FL1, WAVNUM,
!                       UTOP, UDIR, TAUW, TAUWDIR, RNFAC, CICOVER
!                       USTAR, Z0, Z0B, CHRNCK)
!          *KIJS*    - INDEX OF FIRST GRIDPOINT
!          *KIJL*    - INDEX OF LAST GRIDPOINT
!          *IUSFG*   - IF = 1 THEN USE THE FRICTION VELOCITY (US) AS FIRST GUESS in TAUT_Z0
!                           0 DO NOT USE THE FIELD USTAR
!          *FL1*     - 2D-SPECTRA
!          *WAVNUM*  - WAVE NUMBER
!          *HALP*    - 1/2 PHILLIPS PARAMETER
!          *UTOP*    - WIND SPEED AT REFERENCE LEVEL XNLEV
!          *UDIR*    - WIND SPEED DIRECTION AT REFERENCE LEVEL XNLEV
!          *TAUW*    - WAVE STRESS.
!          *TAUWDIR* - WAVE STRESS DIRECTION.
!          *RNFAC*   - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
!          *CICOVER* - SEA ICE COVER

!          *USTAR*   - FRICTION VELOCITY
!          *Z0*      - ROUGHNESS LENGTH
!          *Z0B*     - BACKGROUND ROUGHNESS LENGTH
!          *CHRNCK*  - CHARNOCK COEFFICIENT

!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESS LENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!       IT WAS EXTENDED TO INCLUDE THE GRAVITY-CAPILLARY MODEL FOR THE CALCULATION
!       OF THE BACKGROUND ROUGHNESS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LLCAPCHNK, LLGCBZ0
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G, GM1, EPSUS, EPSMIN, ACD, BCD, ACDLIN, BCDLIN, CDMAX
      USE YOWPHYS  , ONLY : XKAPPA, XNLEV, RNU, RNUM, ALPHA, ALPHAMIN, ALPHAMAX, &
     &                      ANG_GC_A, ANG_GC_B, ANG_GC_C, ANG_GC_MIN
      USE YOWTABL  , ONLY : EPS1 
      USE YOWWIND  , ONLY : WSPMIN_WAVE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "chnkmin.intfb.h"
#include "stress_gc.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL, IUSFG
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: HALP, UTOP, UDIR, TAUW, TAUWDIR, RNFAC, CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(INOUT) :: USTAR
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(OUT) :: Z0, Z0B, CHRNCK


      INTEGER(KIND=JWIM), PARAMETER :: NITER=18

      REAL(KIND=JWRB), PARAMETER :: TWOXMP1=3.0_JWRB
      REAL(KIND=JWRB), PARAMETER :: PMAX=0.99_JWRB

      INTEGER(KIND=JWIM) :: IJ, ITER
      INTEGER(KIND=JWIM) :: IFRPH

      REAL(KIND=JWRB) :: ALPHAGM1

      REAL(KIND=JWRB), PARAMETER :: Z0MIN = 0.000001_JWRB
      REAL(KIND=JWRB), PARAMETER :: PCHARMAX_UTHRS=1.61_JWRB
!          WHEN FINDING Z0 WHICH SATISTIES THE NEUTRAL WIND PROFILE
!          UTOP = (USTAR/XKAPPA)*LOG(1.+XNLEV/Z0)
!          WHERE  Z0=RNUM/USTAR+PCHAR*USTAR**2/G, 
!          IT CAN BE SHOWN THAT THE SCHEME CONVERGES FOR PCHAR VALUES <= PCHARMAX_UTHRS*G*XNLEV/XKAPPA/UTOP**2

      REAL(KIND=JWRB) :: PCE_GC
      REAL(KIND=JWRB) :: Z0MINRST
      REAL(KIND=JWRB) :: CHARNOCK_MIN, ZNCHARMAX
      REAL(KIND=JWRB) :: COSDIFF 
      REAL(KIND=JWRB) :: ZCHAR
      REAL(KIND=JWRB) :: US2TOTAUW, USMAX
      REAL(KIND=JWRB) :: XLOGXL, XKUTOP, XOLOGZ0
      REAL(KIND=JWRB) :: USTOLD, USTNEW, TAUOLD, TAUNEW, X, F, DELF, CDFG
      REAL(KIND=JWRB) :: USNRF, Z0NRF, Z0BNRF, ALPOG
      REAL(KIND=JWRB) :: USTM1, Z0TOT, Z0CH, Z0VIS, HZ0VISO1MX, ZZ
      REAL(KIND=JWRB) :: CONST, TAUV, DEL
      REAL(KIND=JWRB) :: ZNURDC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ALPHAMINOG, XMIN, ZNURDCKAPM1, ZALPHAMAX 
      REAL(KIND=JWRB), DIMENSION(KIJL) :: W1
      REAL(KIND=JWRB), DIMENSION(KIJL) :: TAUWACT, TAUWEFF
      REAL(KIND=JWRB), DIMENSION(KIJL) :: ANG_GC, TAUUNR

      LOGICAL,  DIMENSION(KIJL) :: LLCOSDIFF

! ----------------------------------------------------------------------

#include "cdm.func.h"

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',0,ZHOOK_HANDLE)

      XLOGXL=LOG(XNLEV)
      US2TOTAUW=1.0_JWRB+EPS1

      DO IJ = KIJS, KIJL
!       Below WSPMIN_WAVE we assume very little waves would exist and the viscous stress should not be reduced
        ZNURDC =  0.04_JWRB + 0.48_JWRB*(1.0_JWRB-TANH(5.0_JWRB*(UTOP(IJ)-WSPMIN_WAVE)))
        ZNURDCKAPM1(IJ) = ZNURDC*RNU/XKAPPA
      ENDDO

      PCE_GC = 0.001_JWRB * IUSFG + (1-IUSFG) * 0.005_JWRB

!     ONLY take the contribution of TAUW that is in the wind direction
      DO IJ = KIJS, KIJL
        COSDIFF = COS(UDIR(IJ)-TAUWDIR(IJ))
        TAUWACT(IJ) = MAX(TAUW(IJ)*COSDIFF, EPSMIN )
        LLCOSDIFF(IJ) = (COSDIFF > 0.9_JWRB )
      ENDDO

!  USING THE CG MODEL:
IF (LLGCBZ0) THEN

      IF (LLCAPCHNK) THEN
        DO IJ=KIJS,KIJL
          CHARNOCK_MIN = CHNKMIN(UTOP(IJ), CICOVER(IJ))
          ALPHAMINOG(IJ) = CHARNOCK_MIN*GM1
        ENDDO
      ELSE
        ALPHAMINOG(KIJS:KIJL)= 0.0_JWRB
      ENDIF

      ZNCHARMAX = PCHARMAX_UTHRS*G*XNLEV/XKAPPA
      DO IJ = KIJS, KIJL
         ZALPHAMAX(IJ) = MIN(ZNCHARMAX/UTOP(IJ)**2, ALPHAMAX) 
      ENDDO

      DO IJ = KIJS, KIJL
        USMAX = MAX(-0.15_JWRB + 0.093698_JWRB*UTOP(IJ) -0.0020944_JWRB*UTOP(IJ)**2 + 5.5091E-5_JWRB*UTOP(IJ)**3, 0.03_JWRB)
        TAUWEFF(IJ) = MIN(TAUWACT(IJ)*US2TOTAUW, USMAX**2 )
      ENDDO

      IF (IUSFG == 0 ) THEN
        ALPHAGM1 = ALPHA*GM1
        DO IJ = KIJS, KIJL
          IF ( UTOP(IJ) < WSPMIN_WAVE ) THEN
            CDFG = 0.0011_JWRB
          ELSEIF ( LLCOSDIFF(IJ) ) THEN
            X = MIN(TAUWACT(IJ)/MAX(USTAR(IJ),EPSUS)**2,PMAX)
            ZCHAR = MIN( ALPHAGM1 * USTAR(IJ)**2 / SQRT(1.0_JWRB - X), 0.05_JWRB*EXP(-0.05_JWRB*(UTOP(IJ)-35._JWRB)) )
            ZCHAR = MIN(ZCHAR,ZALPHAMAX(IJ))
            CDFG = ACDLIN + BCDLIN*SQRT(ZCHAR) * UTOP(IJ)
          ELSE
            CDFG = CDM(UTOP(IJ))
          ENDIF
          USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
        ENDDO
      ENDIF

      DO IJ = KIJS, KIJL
        W1(IJ) = 0.85_JWRB - 0.05_JWRB*( TANH(10.0_JWRB*(UTOP(IJ)-5.0_JWRB)) + 1.0_JWRB )
      ENDDO

      DO IJ = KIJS, KIJL
        XKUTOP = XKAPPA * UTOP(IJ)

        USTOLD = USTAR(IJ)
        TAUOLD = USTOLD**2

        DO ITER=1,NITER
!         Z0 IS DERIVED FROM THE NEUTRAL LOG PROFILE: UTOP = (USTAR/XKAPPA)*LOG((XNLEV+Z0)/Z0)
          Z0(IJ) = MAX(XNLEV/(EXP(MIN(XKUTOP/USTOLD, 50.0_JWRB))-1.0_JWRB), Z0MIN)
          ! Viscous kinematic stress nu_air * dU/dz at z=0 of the neutral log profile reduced by factor 25 (0.04)
          TAUV = ZNURDCKAPM1(IJ)*USTOLD/Z0(IJ)

          ANG_GC(IJ) = MAX( ANG_GC_A + ANG_GC_B * TANH(ANG_GC_C * TAUOLD), ANG_GC_MIN)

          TAUUNR(IJ) = STRESS_GC(ANG_GC(IJ), USTAR(IJ), Z0(IJ), Z0MIN, HALP(IJ), RNFAC(IJ))

!         TOTAL kinematic STRESS:
          TAUNEW = TAUWEFF(IJ) + TAUV + TAUUNR(IJ)
          USTNEW = SQRT(TAUNEW)
          USTAR(IJ) = W1(IJ)*USTOLD+(1.0_JWRB-W1(IJ))*USTNEW

!         CONVERGENCE ?
          DEL = USTAR(IJ)-USTOLD
          IF (ABS(DEL) < PCE_GC*USTAR(IJ)) EXIT
          TAUOLD = USTAR(IJ)**2
          USTOLD = USTAR(IJ)
        ENDDO

        X = TAUWEFF(IJ)/TAUOLD

        ! protection just in case there is no convergence
        IF (ITER > NITER .AND. X >= PMAX ) THEN
          CDFG = CDM(UTOP(IJ))
          USTAR(IJ) = UTOP(IJ)*SQRT(CDFG)
          Z0MINRST = USTAR(IJ)**2 * ALPHA*GM1
          Z0(IJ) = MAX(XNLEV/(EXP(XKUTOP/USTAR(IJ))-1.0_JWRB), Z0MINRST)
          Z0B(IJ) = Z0MINRST
        ELSE
          Z0(IJ) = MAX(XNLEV/(EXP(XKUTOP/USTAR(IJ))-1.0_JWRB), Z0MIN)
          Z0B(IJ) = Z0(IJ)*SQRT(TAUUNR(IJ)/TAUOLD)
        ENDIF

!       Refine solution

        IF (X < PMAX) THEN

          USNRF = USTAR(IJ)
          Z0NRF = Z0(IJ)
          Z0BNRF = Z0B(IJ)

          USTOLD = USTAR(IJ)
          TAUOLD = MAX(USTOLD**2,TAUWEFF(IJ))

          ALPOG = MAX(MIN(Z0B(IJ)/TAUOLD,GM1*ZALPHAMAX(IJ)), ALPHAMINOG(IJ))

          DO ITER=1,NITER
!           FIND USTAR WHICH SATISTIES THE NEUTRAL WIND PROFILE
!           UTOP = (USTAR/XKAPPA)*LOG(1.+XNLEV/Z0)
!           WHERE Z0 = HZ0VISO1MX+SQRT(HZ0VISO1MX**2+Z0B**2/(1.0-X))
!           WITH HZ0VISO1MX = 0.5*Z0VIS/(1.0-X), Z0VIS = RNUM/USTAR, Z0B=ALPHOG*USTAR**2, X = TAUEFF/USTAR**2
!
!           A NEWTON METHOD IS USED TO IETRATIVELY FIND UST AND HENCE Z0:
!           GIVEN F=UTOP-(UST/PKAP)*LOG(1.+Z/Z0)
!           UST(n+1) = UST(n) - F/(dF/dUST)

            X = MIN(TAUWEFF(IJ)/TAUOLD, PMAX)
            USTM1 = 1.0_JWRB/MAX(USTOLD,EPSUS)
            Z0VIS = RNUM*USTM1
            HZ0VISO1MX = 0.5_JWRB*Z0VIS/(1.0_JWRB-X)
            Z0B(IJ) = ALPOG*TAUOLD
            Z0(IJ) = HZ0VISO1MX+SQRT(HZ0VISO1MX**2+Z0B(IJ)**2/(1.0_JWRB-X))

            XOLOGZ0= 1.0_JWRB/LOG(XNLEV/Z0(IJ)+1.0_JWRB)
            F = USTOLD-XKUTOP*XOLOGZ0
            ZZ = 2.0_JWRB*USTM1*(3.0_JWRB*Z0B(IJ)**2+0.5_JWRB*Z0VIS*Z0(IJ)-Z0(IJ)**2) &
&                / (2.0_JWRB*Z0(IJ)**2*(1.0_JWRB-X)-Z0VIS*Z0(IJ))

            DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ

            IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTOLD-F/DELF
!           CONVERGENCE ?
            TAUNEW = MAX(USTAR(IJ)**2,TAUWEFF(IJ))
            USTAR(IJ) = SQRT(TAUNEW)
            DEL = TAUNEW-TAUOLD
            IF (ABS(DEL) < PCE_GC*TAUOLD) EXIT
            TAUOLD = TAUNEW
            USTOLD = USTAR(IJ)

          ENDDO
          ! protection just in case there is no convergence
          IF (ITER > NITER ) THEN
            USTAR(IJ) = USNRF
            Z0(IJ) = Z0NRF
            Z0B(IJ) = Z0BNRF
            USTM1 = 1.0_JWRB/MAX(USTAR(IJ), EPSUS)
            Z0VIS = RNUM*USTM1
            CHRNCK(IJ) = MAX(G*(Z0(IJ)-Z0VIS) * USTM1**2, ALPHAMIN)

          ELSE
            CHRNCK(IJ) = MIN(MAX( G*(Z0B(IJ)/SQRT(1.0_JWRB-X))/MAX(USTAR(IJ),EPSUS)**2, ALPHAMIN), ZALPHAMAX(IJ))
          ENDIF

        ELSE
          USTM1 = 1.0_JWRB/MAX(USTAR(IJ), EPSUS)
          Z0VIS = RNUM*USTM1
          CHRNCK(IJ) = MAX(G*(Z0(IJ)-Z0VIS) * USTM1**2, ALPHAMIN)
        ENDIF

      ENDDO


ELSE

      DO IJ = KIJS, KIJL
        TAUWEFF(IJ) = TAUWACT(IJ)*US2TOTAUW
      ENDDO

      IF (LLCAPCHNK) THEN
        DO IJ=KIJS,KIJL
          CHARNOCK_MIN = CHNKMIN(UTOP(IJ), CICOVER(IJ))
          XMIN(IJ) = 0.15_JWRB*(ALPHA-CHARNOCK_MIN)
          ALPHAMINOG(IJ) = CHARNOCK_MIN*GM1
        ENDDO
      ELSE
        DO IJ=KIJS,KIJL
          XMIN(IJ)= 0.0_JWRB
          ALPHAMINOG(IJ)= ALPHA*GM1
        ENDDO
      ENDIF

      DO IJ=KIJS,KIJL
        XKUTOP = XKAPPA * UTOP(IJ)

        USTOLD = (1-IUSFG)*UTOP(IJ)*SQRT(MIN(ACD+BCD*UTOP(IJ),CDMAX)) + IUSFG*USTAR(IJ)
        TAUOLD = MAX(USTOLD**2,TAUWEFF(IJ))
        USTAR(IJ) = SQRT(TAUOLD)
        USTM1 = 1.0_JWRB/MAX(USTAR(IJ),EPSUS) 

        DO ITER=1,NITER
          X = MAX(TAUWACT(IJ)/TAUOLD,XMIN(IJ))
          Z0CH = ALPHAMINOG(IJ)*TAUOLD/SQRT(1.0_JWRB-X)
          Z0VIS = RNUM*USTM1
          Z0TOT = Z0CH+Z0VIS

          XOLOGZ0= 1.0_JWRB/(XLOGXL-LOG(Z0TOT))
          F = USTAR(IJ)-XKUTOP*XOLOGZ0
          ZZ = USTM1*(Z0CH*(2.0_JWRB-TWOXMP1*X)/(1.0_JWRB-X)-Z0VIS)/Z0TOT
          DELF= 1.0_JWRB-XKUTOP*XOLOGZ0**2*ZZ

          IF (DELF /= 0.0_JWRB) USTAR(IJ) = USTAR(IJ)-F/DELF
          TAUNEW = MAX(USTAR(IJ)**2,TAUWEFF(IJ))
          USTAR(IJ) = SQRT(TAUNEW)
          IF (TAUNEW == TAUOLD) EXIT
          USTM1 = 1.0_JWRB/MAX(USTAR(IJ),EPSUS)
          TAUOLD = TAUNEW
        ENDDO

        Z0(IJ) = Z0CH
        Z0B(IJ) = ALPHAMINOG(IJ)*TAUOLD
        CHRNCK(IJ) = MAX(G*Z0(IJ)*USTM1**2, ALPHAMIN)

      ENDDO

ENDIF

IF (LHOOK) CALL DR_HOOK('TAUT_Z0',1,ZHOOK_HANDLE)

END SUBROUTINE TAUT_Z0
