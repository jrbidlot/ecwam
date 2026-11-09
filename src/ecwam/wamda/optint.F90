! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

       SUBROUTINE OPTINT(                                               &
!           INPUT
     &        NPARTW,NPASC2W,                                           &
     &        XMEANEW,                                                  &
     &        EW,XKW,YKW,                                               &
     &        ED,XKD,YKD,                                               &
     &        SPECSA,CORRELTAS,                                         &
     &        PARTS,                                                    &
!            OUTPUT
     &        EA,XANGA,XFREA)

!------------------------------------------------------------------

!     PURPOSE
!     -------
!     OPTIMAL INTERPOLATION SCHEME TO COMPUTE CORRECTION FACTORS FOR
!     ALL WAVE SYSTEMS AT WAM POINTS INFLUENCED BY SAR DATA POINTS.

!     AUTHOR
!     ------
!     S. HASSELMANN, MAX PLANCK INSTITUT FUER METEOROLOGIE,
!                    HAMBURG, GER.

!     EXTERNALS
!     ---------
!     WAM_SYMINV

!     METHOD.
!     -------
!     THE SAR DATA ERROR CORRELATION MATRIX IS COMPUTED FROM AN INVERSE
!     EXPONENTIAL DISTANCE BETWEEN THE CORRELATED POINTS.
!     THE MODEL ERROR CORRELATION MATRIX IS ASSUMED TO BE THE IDENTITY
!     MATRIX. 

!**********************************************************************

!     MODULES:
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,R        ,RAD      ,ZPI
      USE YOWSARAS , ONLY : NSPEC    ,NSPECW   ,JPSIW    ,MPARTSW  ,    &
     &            MAXMPART ,LONG     ,LAT      ,COSSARLAT,SINSARLAT,    &
     &            NTOSIW1D ,NSIW1D   ,SPECW    ,DIST     ,SARCORDIA,    &
     &            COEFPHI
      USE YOWSTAT  , ONLY : LSARAS   , LSWIMAS
      USE YOWTEST  , ONLY : IU06
      USE EC_LUN   , ONLY : NULERR

   
!     INTERFACE
!     ---------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "wam_syminv.intfb.h"

!     INTEGERS WAM - SAR INDICES
      INTEGER(KIND=JWIM) ::                                             &
     &      NPARTW(NSPECW),                                             &
!              NUMBER OF PARTITIONINGS IN WAM FIRST GUESS SPECTRUM.
     &      NPASC2W(NSPECW,MPARTSW,JPSIW) 
!              INDEX OF PARTITIONING IN SAR SPECTRUM NSIW1D(MSPECW,JPSIW)
!              CROSSASSIGNED TO PARTITIONING MPARTSW OF WAM SPECTRUM MSPEC.

      REAL(KIND=JWRB) ::                                                &
     &      SARCORLEN,                                                  &
!              CORRELATION LENGTHS (in radian) !!!
     &      WEIGHT(NSPECW,MPARTSW,JPSIW), WEIGHTX(NSPECW,MPARTSW,JPSIW),&
!              INTERPOLATION WEIGHTS,
     &      EW(NSPECW,MPARTSW),                                         &
     &      XKW(NSPECW,MPARTSW),                                        &
     &      YKW(NSPECW,MPARTSW),                                        &
!              CHARACTERISTIC PARAMETERS OF WAM FIRST GUESS WAVE SYSTEMS.
     &      ED(NSPEC,0:MPARTSW),                                        &
     &      XKD(NSPEC,0:MPARTSW),                                       &
     &      YKD(NSPEC,0:MPARTSW)
!              ERRORS (difference) BETWEEN CHARACTERISTIC PARAMETERS
!              OF SAR WAVE SYSTEMS AND WAM WAVE SYSTEMS AT SAR DATA POINTS
!              (ENERGY AND MEAN WAVE NUMBER IN LAT AND LON DIRECTION).
!              MEAN FREQUENCY OF SAR WAVE SYSTEMS.

      REAL(KIND=JWRB) :: XMEANEW(NSPEC,MPARTSW)
!              MODEL FG  MEAN ENERGIES OF PARTITIONINGS AT SAR LOCATIONS
      REAL(KIND=JWRB) :: EA(NSPECW,MPARTSW), XKA(NSPECW,MPARTSW),       &
     &                   YKA(NSPECW,MPARTSW), XANGA(NSPECW,MPARTSW),    &
     &                   XFREA(NSPECW,MPARTSW)
!              CORRECTION PARAMETERS FOR WAM WAVE SYSTEMS.
      REAL(KIND=JWRB) :: PHI(NSPECW,MPARTSW,JPSIW), PHIX(NSPECW,MPARTSW,JPSIW)
!              CROSS CORRELATION BETWEEN  MODEL AND DATA ERROR.

!!!! in the long run more arrays will have to be allocatable
      REAL(KIND=JWRB), ALLOCATABLE :: SIGMA(:,:)
!              ERROR COVARIANCE MATRIX
!     
      INTEGER(KIND=JWIM) :: IFRE, IANG, IDSAR, IPARTS
      REAL(KIND=JWRB) :: SPECSA(0:NSPEC,NANG,NFRE)
      INTEGER(KIND=JWIM) :: CORRELTAS(NSPECW,MPARTSW,JPSIW),            &
     &                      PARTS(NSPEC,NANG,NFRE)

      REAL(KIND=JWRB) :: WSUME, WSUMEX, COND
!            SUMS OF ERRORS.
      REAL(KIND=JWRB) :: SIGMASAR(NSPEC,NSPEC)
!     CORRELATION FUNCTION BETWEEN SAR OBS.

      REAL(KIND=JWRB) :: CORDIST,STRETCHCOR,YY

      INTEGER(KIND=JWIM) :: ISPEC, JSAR, KSAR, IPARTSAR, JPARTSAR,      &
     &                      NJS, K, J, KK, JJ
!            LOOP INDICES.
      INTEGER(KIND=JWIM) :: NJSAR
!             NUMBER OF SAR DATA POINTS INFLUENCING WAM POINT.

      INTEGER(KIND=JWIM) :: MAXMPARTPE
!     MAXIMIUM OF PARTITIONS ON THAT PE.

      REAL(KIND=JWRB) :: ZPHIEAMP, ZPHIEEXP , ZPHIXAMP, ZPHIXEXP
!     CROSS CORRELATION BETWEEN MODEL AND DATA ERROR ARE SPECIFIED AS
!     FOR ENERGY :      ZPHIEAMP*EXP(-(DIST(ISPEC,JSAR)/CORDIST)**ZPHIEEXP)
!     FOR WAVE NUMBER:  ZPHIXAMP*EXP(-(DIST(ISPEC,JSAR)/CORDIST)**ZPHIXEXP)

!------------------------------------------------------------------

!     1. 

      NJSAR=0
      DO ISPEC=1,NSPECW
        NJSAR=MAX(NJSAR,NTOSIW1D(ISPEC))
      ENDDO
      IF (NJSAR > JPSIW) THEN
        WRITE(NULERR,*)'OPTINT ERROR: NUMBER OF SAR POINTS INFLUENCING WAM'
        WRITE(NULERR,*)'POINTS, NJSAR: ',NJSAR
        WRITE(NULERR,*)'LARGER THEN DIMENSION OF INVERSE MATRIX.'
        WRITE(NULERR,*)'JPSIW: ',JPSIW
        CALL ABORT1
      ENDIF

!------------------------------------------------------------------

!      2. COMPUTE CROSS CORRELATION BETWEEN MODEL AND DATA ERROR.
!      -----------------------------------------------------------

      SARCORLEN=SARCORDIA/R

      IF (LSARAS) THEN
        ZPHIEAMP = 0.9_JWRB
        ZPHIEEXP = 0.9_JWRB
        ZPHIXAMP = 0.8_JWRB
        ZPHIXEXP = 0.9_JWRB
      ELSEIF (LSWIMAS) THEN
!     JB: this was kept from the MF implememtation but a better control of the error function PHI in OPTINT shoult be introduced !!
        ZPHIEAMP = COEFPHI
        ZPHIEEXP = 0.9_JWRB
        ZPHIXAMP = 0.8_JWRB
        ZPHIXEXP = 0.9_JWRB
      ELSE
        WRITE(IU06,*) ' +++++++++++++++++++++++++++++'
        WRITE(IU06,*) ' + ERROR IN --OPTINT--       +'
        WRITE(IU06,*) ' + UKNOWN PARAMETERS FOR PHI +'
        WRITE(IU06,*) ' + PROGRAM WILL ABORT        +'
        WRITE(IU06,*) ' +++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      MAXMPARTPE=0
      DO ISPEC=1,NSPECW
        MAXMPARTPE=MAX(MAXMPARTPE,NPARTW(ISPEC))
      ENDDO

      DO JSAR=1,JPSIW
        DO IPARTS=1,MAXMPARTPE
          DO ISPEC=1,NSPECW
            IF (IPARTS.LE.NPARTW(ISPEC).AND.JSAR.LE.NTOSIW1D(ISPEC))THEN
              STRETCHCOR=1.0_JWRB
              IPARTSAR = NPASC2W(ISPEC,IPARTS,JSAR)
              IF (IPARTSAR.GT.0 ) THEN
                YY=XMEANEW(NSIW1D(ISPEC,JSAR),IPARTSAR)
                IF (YY.GT.0.0_JWRB) THEN
                  STRETCHCOR = ABS(EW(ISPEC,IPARTS)-YY)/YY
                  STRETCHCOR = MAX(EXP(-STRETCHCOR),0.1_JWRB) 
                ENDIF
              ENDIF
              CORDIST=STRETCHCOR*SARCORLEN
              PHI(ISPEC,IPARTS,JSAR) = ZPHIEAMP*EXP(-(DIST(ISPEC,JSAR)/CORDIST)**ZPHIEEXP)
              PHIX(ISPEC,IPARTS,JSAR)= ZPHIXAMP*EXP(-(DIST(ISPEC,JSAR)/CORDIST)**ZPHIXEXP)
            ENDIF
          ENDDO
        ENDDO
      ENDDO


!-------------------------------------------------------------------

!      3. COMPUTE INTERPOLATION WEIGHTS.
!      ---------------------------------

      ALLOCATE(SIGMA(NJSAR,NJSAR))


      DO K = 1,NSPEC
        DO J = 1,K-1
          SIGMASAR(K,J) = COS(RAD*(LONG(K,2)-LONG(J,2))) *              &
     &                      COSSARLAT(J)*COSSARLAT(K) +                 &
     &                      SINSARLAT(J)*SINSARLAT(K)
          SIGMASAR(K,J) = MAX(MIN(SIGMASAR(K,J),1.0_JWRB),-1.0_JWRB) 

          ! MF ADDITION: GEOGRAPHIC CORELLATION LENGTH DEPENDENCE
          IF (ABS(LAT(K,2)) <= 20._JWRB) THEN
            SARCORLEN = 3.0_JWRB*SARCORDIA/R
          ELSEIF (ABS(LAT(K,2)) <= 50._JWRB) THEN
            SARCORLEN = 1.6_JWRB*SARCORDIA/R
          ELSE
            SARCORLEN = 0.8_JWRB*SARCORDIA/R
          ENDIF
          SIGMASAR(K,J) = EXP(-(ACOS(SIGMASAR(K,J))/SARCORLEN))

          SIGMASAR(J,K) = SIGMASAR(K,J)
        ENDDO
        SIGMASAR(K,K) = 1.0_JWRB
      ENDDO

      WEIGHT(:,:,:) = 0.0_JWRB
      WEIGHTX(:,:,:) = 0.0_JWRB

      DO ISPEC=1,NSPECW

        DO IPARTS=1,NPARTW(ISPEC)
!         create symmetric matrix SIGMA (WAM_SYMINV only requires the lower
!         triangular part.)
          NJS=0
          KSAR=0
          DO KK = 1,NTOSIW1D(ISPEC)
            K=NSIW1D(ISPEC,KK)
            IPARTSAR = NPASC2W(ISPEC,IPARTS,KK)
            IF (IPARTSAR.GT.0) THEN
              NJS=NJS+1
              KSAR=KSAR+1
              JSAR=0
              DO JJ = 1,KK-1
                J=NSIW1D(ISPEC,JJ)
                JPARTSAR = NPASC2W(ISPEC,IPARTS,JJ)
                IF (JPARTSAR.GT.0) THEN
                  JSAR=JSAR+1
                  SIGMA(KSAR,JSAR) = SIGMASAR(K,J) 
                ENDIF
              ENDDO
              SIGMA(KSAR,KSAR) = 1.0_JWRB + SIGMASAR(K,K) 
            ENDIF
          ENDDO

!         INVERT MATRIX SIGMA.

          IF (NJS.GT.1) THEN
            COND=0.0_JWRB
            CALL WAM_SYMINV(SIGMA,NJSAR,NJS,COND)
          ELSE IF (NJS.EQ.1) THEN
            SIGMA(1,1)=1.0_JWRB/SIGMA(1,1)
          ENDIF

          DO KSAR=1,NJS
            DO JSAR=1,KSAR-1
              SIGMA(JSAR,KSAR)=SIGMA(KSAR,JSAR)
            ENDDO
          ENDDO

          JSAR=0
          DO JJ = 1,NTOSIW1D(ISPEC)
            JPARTSAR = NPASC2W(ISPEC,IPARTS,JJ)
            IF (JPARTSAR.GT.0) THEN
              JSAR=JSAR+1
              WSUME=0.0_JWRB
              WSUMEX=0.0_JWRB
              KSAR=0
              DO KK = 1,NTOSIW1D(ISPEC)
                K=NSIW1D(ISPEC,KK)
                IPARTSAR = NPASC2W(ISPEC,IPARTS,KK)
                IF (IPARTSAR.GT.0) THEN
                  KSAR=KSAR+1
                  WSUME = WSUME + SIGMA(JSAR,KSAR)*PHI(ISPEC,IPARTS,KK)
                  WSUMEX=WSUMEX+SIGMA(JSAR,KSAR)*PHIX(ISPEC,IPARTS,KK)
                ENDIF
              ENDDO
              WEIGHT(ISPEC,IPARTS,JJ) = WSUME
              WEIGHTX(ISPEC,IPARTS,JJ)=WSUMEX
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(SIGMA)

!------------------------------------------------------------------

!     4.  COMPUTE VALUES TO CORRECT WAM FIRST GUESS WAVE SYSTEMS
!         INFLUENCING WAM POINT.
!     -----------------------------------------------------------

      EA(:,:) = 0.0_JWRB
      XKA(:,:) = 0.0_JWRB
      YKA(:,:) = 0.0_JWRB

      DO JSAR=1,JPSIW
        DO IPARTS=1,MAXMPARTPE
          DO ISPEC=1,NSPECW
            IF (IPARTS.LE.NPARTW(ISPEC).AND.JSAR.LE.NTOSIW1D(ISPEC))THEN
              IPARTSAR = NPASC2W(ISPEC,IPARTS,JSAR)

              IF (IPARTSAR.GT.0) THEN
                EA(ISPEC,IPARTS) = EA(ISPEC,IPARTS)+                    &
     &                             WEIGHT(ISPEC,IPARTS,JSAR)*           &
     &                             ED(NSIW1D(ISPEC,JSAR),IPARTSAR)
                XKA(ISPEC,IPARTS) = XKA(ISPEC,IPARTS)+                  &
     &                              WEIGHTX(ISPEC,IPARTS,JSAR)*         &
     &                              XKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
                YKA(ISPEC,IPARTS) = YKA(ISPEC,IPARTS)+                  &
     &                              WEIGHTX(ISPEC,IPARTS,JSAR)*         &
     &                              YKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!!    insure that the increment does no produce negative energy
      DO IPARTS=1,MAXMPART
        DO ISPEC=1,NSPECW
          EA(ISPEC,IPARTS)=MAX(EW(ISPEC,IPARTS)+EA(ISPEC,IPARTS),0.0_JWRB)
!!!! test
!!!! always produce an analysis for wave number even when the
!!!! system is reduced to zero energy. It still needs to be shifted and
!!!! rotated !!!!!
!!!!          IF (EA(ISPEC,IPARTS).GT.0.) THEN
            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS) + XKA(ISPEC,IPARTS)
            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS) + YKA(ISPEC,IPARTS)
!!!!          ELSE
!!!!            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS)
!!!!            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS)
!!!!          ENDIF
        ENDDO
      ENDDO

!     COMPUTE MEAN ANALYSED DIRECTIONS AND FREQUENCIES FROM 
!     WAVE NUMBERS.

      DO IPARTS=1,MAXMPARTPE
        DO ISPEC=1,NSPECW
          IF (IPARTS.LE.NPARTW(ISPEC)) THEN
            XANGA(ISPEC,IPARTS) =                                       &
     &               ATAN2(YKA(ISPEC,IPARTS),XKA(ISPEC,IPARTS))
            XFREA(ISPEC,IPARTS) = SQRT(SQRT(XKA(ISPEC,IPARTS)**2+       &
     &               YKA(ISPEC,IPARTS)**2)*G)/ZPI
          ENDIF
        ENDDO
      ENDDO

!------------------------------------------------------------------

!     5. SUPERIMPOSE EXTRA WAVE SYSTEMS FROM SAR RETRIEVAL ON 
!        FIRST GUESS SPECTRUM.
!     ----------------------------------------------------------

!!!   this is no longer used !!!!

      END SUBROUTINE OPTINT
