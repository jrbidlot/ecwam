! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDICE3 (KIJS, KIJL, FL1, FLD, SL, SLICE,     &
     &                   WAVNUM, CGROUP,                      &
     &                   CICV, CITH, ALPFAC)
! ----------------------------------------------------------------------

!**** *SDICE3* - COMPUTATION OF SEA ICE ATTENUATION DUE TO VISCOUS FRICTION


!     LOTFI AOUF       METEO FRANCE 2023
!     JOSH KOUSAL      ECMWF 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *SDICE3 (KIJS, KIJL, FL1, FLD,SL,
!                       WAVNUM, CGROUP,                      
!                       CICV, CITH, ALPFAC)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *SLICE*  - TOTAL ICE SOURCE FUNCTION ARRAY
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *CICV*   - SEA ICE COVER
!          *CITH*   - SEA ICE THICKNESS
!          *ALPFAC* - FACTOR TO REDUCE ATTENUATION IN ICE. EQUAL TO 1 IN SOLID ICE

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     Jie Yu * , W. Erick Rogers, David W. Wang, 2022 
!     DOI:10.1016/j.coldregions.2022.103582

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWNEMOCOUWRS
      USE YOWFRED  , ONLY : FR      
      USE YOWICE   , ONLY : ZALPFACB
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPCONS , ONLY : G       ,ZPI
      USE YOWSTAT  , ONLY : IDELT   ,XIMP

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(OUT) :: SLICE
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICV
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CITH
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: ALPFAC


      INTEGER(KIND=JWIM) :: IMODEL       !! DAMPING MODEL: 1=FIT TO TEMPELFJORD DATA, 2=Jie Yu 2022
      INTEGER(KIND=JWIM) :: IJ, K, M

      REAL(KIND=JWRB)    :: DELTM, DELT5, DELT, GTEMP1
      REAL(KIND=JWRB)    :: CDICE, CDICE2
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(KIJL) :: ALPTMP
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: ALP  !! ALP=SPATIAL ATTENUATION RATE OF ENERGY
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE) :: FLDICE 
      REAL(KIND=JWRB), DIMENSION(NFRE)      :: FRP

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDICE3',0,ZHOOK_HANDLE)

      DELT  = IDELT
      DELTM = 1.0_JWRB/DELT
      DELT5 = XIMP*DELT

      IMODEL = 2

      SELECT CASE (IMODEL)
         CASE (1)
!          Best fit w Tempelfjorde obs from Lotfi Aouf
           CDICE=0.0656_JWRB

           DO M = 1,NFRE
             DO IJ = KIJS,KIJL
                ALP(IJ,M) = (CDICE*CITH(IJ)*WAVNUM(IJ,M)**2) * ALPFAC(IJ)
             ENDDO
           ENDDO

         CASE (2)
!          Jie Yu, W. Erik Rogers, David W. Wang 2022
           CDICE=0.1274_JWRB*( ZPI/SQRT(G) )**(4.5_JWRB)

           DO M = 1,NFRE
             FRP(M) = FR(M)**4.5_JWRB
           ENDDO
           CDICE2 = 2._JWRB*CDICE*ZALPFACB
           DO IJ = KIJS,KIJL
             ALPTMP(IJ) = (CDICE2*(CITH(IJ)**1.25_JWRB)) * ALPFAC(IJ) 
           ENDDO
           DO M = 1,NFRE
             DO IJ = KIJS,KIJL
               ALP(IJ,M) = ALPTMP(IJ)*FRP(M)
             ENDDO
           ENDDO
         
      END SELECT

      DO M = 1,NFRE
         DO IJ = KIJS,KIJL
           FLDICE(IJ,M) = -ALP(IJ,M)*CGROUP(IJ,M)
         ENDDO
      ENDDO

      DO M = 1,NFRE
        DO K = 1,NANG
          DO IJ = KIJS,KIJL
!           apply the source term
            SLICE(IJ,K,M) =  FL1(IJ,K,M) * FLDICE(IJ,M)
            SL(IJ,K,M)    =  SL(IJ,K,M)  + CICV(IJ)*SLICE(IJ,K,M)
            FLD(IJ,K,M)   =  FLD(IJ,K,M) + CICV(IJ)*FLDICE(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      IF (LWNEMOCOUWRS) THEN
        DO M = 1,NFRE
          DO K = 1,NANG
            DO IJ = KIJS,KIJL
!             to be used for wave radiative stress calculation
              GTEMP1        =  MAX((1.0_JWRB-DELT5*FLDICE(IJ,M)),1.0_JWRB)
              SLICE(IJ,K,M) =  SLICE(IJ,K,M)/GTEMP1
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      
      IF (LHOOK) CALL DR_HOOK('SDICE3',1,ZHOOK_HANDLE)

      END SUBROUTINE SDICE3
