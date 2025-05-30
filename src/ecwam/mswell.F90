! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MSWELL (XLON, YLAT, FL1)
! ----------------------------------------------------------------------

!**** *MSWELL* - MAKES START SWELL FIELDS FOR WAMODEL.

!     J. BIDLOT     ECMWF    

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSWELL (XLON, YLAT, FL1)
!      *XLON*     REAL      GRID POINT LONGITUDES
!      *XLAT*     REAL      GRID POINT LATITUDES
!      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 

!     METHOD.
!     -------

!       NLOC SWELL SYSTEMS ARE SPECIFIED (SEE BELOW), CENTERED AT
!       DIFFERENT (XLON0, YLAT0) COORDINATES.
!
!       EACH SWELL SPECTRUM AT (XLON0, YLAT0) IS SPECIFIED AS
!       F(OMEGA,THETA)=F1(OMEGA)*Q(THETA)
!       WHERE
!
!       F1(OMEGA)=(N+1)*m0*(OMEGA_p**N/OMEGA**(N+1))*EXP[-((N+1)/N)*(OMEGA_p/OMEGA)**N]
!       WITH
!          OMEGA=2PI*f ; f: wave frequency
!          m0 the zero spectral moment = Hs**2/(2PI) ; Hs the wave height (to be specified) 
!          OMEGA_p=2PI*f_p  ; f_p the peak frequency 
!          N=5
!
!       Q(THETA)=(8/3PI)*COS(THETA-THETA0)**4 for ABS(THETA-THATA0) >  PI/2 
!       Q(THETA)=0                            for ABS(THETA-THATA0) <= PI/2 
!       WITH
!           THETA the wave direction
!           THETA0 the mean wave direction of the swell system (to be specified)
!
!       THE SPECTRUM AT EACH (XLON0,YLAT0) is SPREAD ON THE SPHERE BY MULTIPLYING
!       F(OMEGA,THETA) BY SPEAD(XLON,YLAT)
!       WHERE
!       SPREAD(XLON,YLAT)=EXP[-2*DIST/L}]
!       WITH
!          DIST THE DISTANCE BETWEEN (XLON0,YLAT0)  and (XLON,YLAT)
!          L the correlation distance of the initial disturbance to be specified (in m)
!
!       NOTE THAT XLON,YLAT,XLON0,YLAT0 WILL SPECIFIED IN DEGREE BUT ARE CONVERTED 
!       IN RADIAN FOR CALCULATIONS.
!

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       I. V. LAVRENOV, 2003: WIND-WAVES IN OCEANS, DYNAMICS AND NUMERICAl SIMILATIONS
!       SPRINGER-VERLAG, 376pp.
!       see pages 52-53

!       NONE.

! ----------------------------------------------------------------------

!
!       F1(OMEGA)=(N+1)*m0*(OMEGA_p**N/OMEGA**(N+1))*EXP[-((N+1)/N)*(OMEGA_p/OMEGA)**N]
!       WITH
!          OMEGA=2PI*f ; f: wave frequency
!          m0 the zero spectral moment = Hs**2/(2PI) ; Hs the wave height (to be specified) 
!          OMEGA_p=2PI*f_p  ; f_p the peak frequency 
!          N=5
!
!       Q(THETA)=(8/3PI)*COS(THETA-THETA0)**4 for ABS(THETA-THATA0) >  PI/2 
!       Q(THETA)=0                            for ABS(THETA-THATA0) <= PI/2 
!       WITH
!           THETA the wave direction
!           THETA0 the mean wave direction of the swell system (to be specified)
!
!       THE SPECTRUM AT EACH (XLON0,YLAT0) is SPREAD ON THE SPHERE BY MULTIPLYING
!       F(OMEGA,THETA) BY SPEAD(XLON,YLAT)
!       WHERE
!       SPREAD(XLON,YLAT)=EXP[-2*DIST/L}]
!       WITH
!          DIST THE DISTANCE BETWEEN (XLON0,YLAT0)  and (XLON,YLAT)
!          L the correlation distance of the initial disturbance to be specified (in m)
!
!       NOTE THAT XLON,YLAT,XLON0,YLAT0 WILL SPECIFIED IN DEGREE BUT ARE CONVERTED 
!       IN RADIAN FOR CALCULATIONS.
!

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       I. V. LAVRENOV, 2003: WIND-WAVES IN OCEANS, DYNAMICS AND NUMERICAl SIMILATIONS
!       SPRINGER-VERLAG, 376pp.
!       see pages 52-53

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM,NCHNK), INTENT(IN) :: XLON, YLAT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(OUT) :: FL1

      FL1=0._JWRB

      END SUBROUTINE MSWELL
