MODULE YOWLSQ

  USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: LSQ

CONTAINS

!-----------------------------------------------------------------------
! Solve the least-squares problem
!
!     A X = B
!
! using singular value decomposition.
!
! Singular values smaller than
!
!     TOL * MAXVAL(W)
!
! are set to zero before back substitution.
!
! KERR = 0 : successful completion
! KERR > 0 : singular value KERR failed to converge
!-----------------------------------------------------------------------

SUBROUTINE LSQ(A, B, TOL, X, KERR)

  IMPLICIT NONE

  REAL(KIND=JWRB), INTENT(IN)  :: A(:,:)
  REAL(KIND=JWRB), INTENT(IN)  :: B(:)
  REAL(KIND=JWRB), INTENT(IN)  :: TOL

  REAL(KIND=JWRB), INTENT(OUT) :: X(:)
  INTEGER(KIND=JWIM), INTENT(OUT) :: KERR

  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM) :: N

  REAL(KIND=JWRB) :: WMIN
  REAL(KIND=JWRB) :: WMAX
  REAL(KIND=JWRB) :: U(SIZE(A, DIM=1), SIZE(A,DIM=2))
  REAL(KIND=JWRB) :: V(SIZE(A, DIM=1), SIZE(A,DIM=2))
  REAL(KIND=JWRB) :: W(SIZE(A, DIM=2))

  KERR = 0

  M = SIZE(A, DIM=1)
  N = SIZE(A, DIM=2)

!-----------------------------------------------------------------------
! Basic consistency checks.
!
! These retain STOP because they indicate programming/interface errors,
! rather than a numerical convergence failure.
!-----------------------------------------------------------------------

  IF (N > M) THEN
    WRITE(*,*) 'LSQ: number of columns must not exceed number of rows'
    STOP
  END IF

  IF (SIZE(B) < M) THEN
    WRITE(*,*) 'LSQ: vector B is too small'
    STOP
  END IF

  IF (SIZE(X) < N) THEN
    WRITE(*,*) 'LSQ: vector X is too small'
    STOP
  END IF

  IF (SIZE(U, DIM=1) < M .OR. SIZE(U, DIM=2) < N) THEN
    WRITE(*,*) 'LSQ: array U has incorrect dimensions'
    STOP
  END IF

  IF (SIZE(V, DIM=1) < N .OR. SIZE(V, DIM=2) < N) THEN
    WRITE(*,*) 'LSQ: array V has incorrect dimensions'
    STOP
  END IF

  IF (SIZE(W) < N) THEN
    WRITE(*,*) 'LSQ: vector W is too small'
    STOP
  END IF

  IF (TOL < 0.0_JWRB) THEN
    WRITE(*,*) 'LSQ: tolerance must not be negative'
    STOP
  END IF

!-----------------------------------------------------------------------
! Copy A because SVDCMP overwrites its matrix argument.
!-----------------------------------------------------------------------

  U(1:M,1:N) = A(1:M,1:N)
  X(1:N)     = 0.0_JWRB

!-----------------------------------------------------------------------
! Singular value decomposition.
!-----------------------------------------------------------------------

  CALL SVDCMP(U(1:M,1:N), W(1:N), V(1:N,1:N), KERR)

!-----------------------------------------------------------------------
! If the SVD has not converged, return control to the caller.
! Do not attempt the back substitution with a partial decomposition.
!-----------------------------------------------------------------------

  IF (KERR /= 0) THEN
    RETURN
  END IF

!-----------------------------------------------------------------------
! Remove sufficiently small singular values.
!-----------------------------------------------------------------------

  WMAX = MAXVAL(W(1:N))
  WMIN = TOL * WMAX

  WHERE (W(1:N) < WMIN)
    W(1:N) = 0.0_JWRB
  END WHERE

!-----------------------------------------------------------------------
! Solve the least-squares problem.
!-----------------------------------------------------------------------

  CALL SVBKSB(U(1:M,1:N), W(1:N), V(1:N,1:N), B(1:M), X(1:N))

END SUBROUTINE LSQ


!=======================================================================
! SVD back substitution.
!
! Given
!
!     A = U W V^T
!
! solve A X = B using
!
!     X = V W^-1 U^T B
!=======================================================================

SUBROUTINE SVBKSB(U, W, V, B, X)

  IMPLICIT NONE

  REAL(KIND=JWRB), INTENT(IN)  :: U(:,:)
  REAL(KIND=JWRB), INTENT(IN)  :: W(:)
  REAL(KIND=JWRB), INTENT(IN)  :: V(:,:)
  REAL(KIND=JWRB), INTENT(IN)  :: B(:)
  REAL(KIND=JWRB), INTENT(OUT) :: X(:)

  INTEGER(KIND=JWIM) :: I
  INTEGER(KIND=JWIM) :: J
  INTEGER(KIND=JWIM) :: JJ
  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM) :: N

  REAL(KIND=JWRB) :: S
  REAL(KIND=JWRB) :: TMP(SIZE(W))

  M = SIZE(U, DIM=1)
  N = SIZE(U, DIM=2)

!-----------------------------------------------------------------------
! Compute W^-1 U^T B.
!-----------------------------------------------------------------------

  DO J = 1, N

    S = 0.0_JWRB

    IF (W(J) /= 0.0_JWRB) THEN

      DO I = 1, M
        S = S + U(I,J) * B(I)
      END DO

      S = S / W(J)

    END IF

    TMP(J) = S

  END DO

!-----------------------------------------------------------------------
! Compute X = V TMP.
!-----------------------------------------------------------------------

  DO J = 1, N

    S = 0.0_JWRB

    DO JJ = 1, N
      S = S + V(J,JJ) * TMP(JJ)
    END DO

    X(J) = S

  END DO

END SUBROUTINE SVBKSB


!=======================================================================
! Singular value decomposition.
!
! On entry:
!
!     A = input M x N matrix, N <= M
!
! On return:
!
!     A = U
!     W = singular values
!     V = V
!
! such that
!
!     A_original = U W V^T
!
! KERR = 0 : successful decomposition
! KERR > 0 : singular value KERR failed to converge
!=======================================================================

SUBROUTINE SVDCMP(A, W, V, KERR)

  IMPLICIT NONE

  REAL(KIND=JWRB), INTENT(INOUT) :: A(:,:)
  REAL(KIND=JWRB), INTENT(OUT)   :: W(:)
  REAL(KIND=JWRB), INTENT(OUT)   :: V(:,:)

  INTEGER(KIND=JWIM), INTENT(OUT) :: KERR

  INTEGER(KIND=JWIM), PARAMETER :: MAX_ITERATIONS = 30

  INTEGER(KIND=JWIM) :: I
  INTEGER(KIND=JWIM) :: ITS
  INTEGER(KIND=JWIM) :: J
  INTEGER(KIND=JWIM) :: JJ
  INTEGER(KIND=JWIM) :: K
  INTEGER(KIND=JWIM) :: L
  INTEGER(KIND=JWIM) :: M
  INTEGER(KIND=JWIM) :: N
  INTEGER(KIND=JWIM) :: NM

  REAL(KIND=JWRB) :: ANORM
  REAL(KIND=JWRB) :: C
  REAL(KIND=JWRB) :: F
  REAL(KIND=JWRB) :: G
  REAL(KIND=JWRB) :: H
  REAL(KIND=JWRB) :: S
  REAL(KIND=JWRB) :: ZSCALE
  REAL(KIND=JWRB) :: X
  REAL(KIND=JWRB) :: Y
  REAL(KIND=JWRB) :: Z

  REAL(KIND=JWRB) :: RV1(SIZE(A,DIM=2))

  LOGICAL :: CANCELLATION_REQUIRED

  KERR = 0

  M = SIZE(A, DIM=1)
  N = SIZE(A, DIM=2)

!-----------------------------------------------------------------------
! Householder reduction to bidiagonal form.
!-----------------------------------------------------------------------

  G      = 0.0_JWRB
  ZSCALE = 0.0_JWRB
  ANORM  = 0.0_JWRB

  DO I = 1, N

    L      = I + 1
    RV1(I) = ZSCALE * G
    G      = 0.0_JWRB
    S      = 0.0_JWRB
    ZSCALE = 0.0_JWRB

    IF (I <= M) THEN

      DO K = I, M
        ZSCALE = ZSCALE + ABS(A(K,I))
      END DO

      IF (ZSCALE /= 0.0_JWRB) THEN

        DO K = I, M
          A(K,I) = A(K,I) / ZSCALE
          S      = S + A(K,I) * A(K,I)
        END DO

        F      = A(I,I)
        G      = -SIGN(SQRT(S),F)
        H      = F * G - S
        A(I,I) = F - G

        IF (I < N) THEN

          DO J = L, N

            S = 0.0_JWRB

            DO K = I, M
              S = S + A(K,I) * A(K,J)
            END DO

            F = S / H

            DO K = I, M
              A(K,J) = A(K,J) + F * A(K,I)
            END DO

          END DO

        END IF

        DO K = I, M
          A(K,I) = ZSCALE * A(K,I)
        END DO

      END IF

    END IF

    W(I)  = ZSCALE * G
    G     = 0.0_JWRB
    S     = 0.0_JWRB
    ZSCALE= 0.0_JWRB

    IF (I <= M .AND. I /= N) THEN

      DO K = L, N
        ZSCALE = ZSCALE + ABS(A(I,K))
      END DO

      IF (ZSCALE /= 0.0_JWRB) THEN

        DO K = L, N
          A(I,K) = A(I,K) / ZSCALE
          S      = S + A(I,K) * A(I,K)
        END DO

        F      = A(I,L)
        G      = -SIGN(SQRT(S),F)
        H      = F * G - S
        A(I,L) = F - G

        DO K = L, N
          RV1(K) = A(I,K) / H
        END DO

        DO J = L, M

          S = 0.0_JWRB

          DO K = L, N
            S = S + A(J,K) * A(I,K)
          END DO

          DO K = L, N
            A(J,K) = A(J,K) + S * RV1(K)
          END DO

        END DO

        DO K = L, N
          A(I,K) = ZSCALE * A(I,K)
        END DO

      END IF

    END IF

    ANORM = MAX(ANORM, ABS(W(I)) + ABS(RV1(I)))

  END DO

!-----------------------------------------------------------------------
! Accumulate right-hand transformation matrix V.
!-----------------------------------------------------------------------

  G = 0.0_JWRB
  L = N + 1

  DO I = N, 1, -1

    IF (I < N) THEN

      IF (G /= 0.0_JWRB) THEN

        DO J = L, N
          V(J,I) = (A(I,J) / A(I,L)) / G
        END DO

        DO J = L, N

          S = 0.0_JWRB

          DO K = L, N
            S = S + A(I,K) * V(K,J)
          END DO

          DO K = L, N
            V(K,J) = V(K,J) + S * V(K,I)
          END DO

        END DO

      END IF

      DO J = L, N
        V(I,J) = 0.0_JWRB
        V(J,I) = 0.0_JWRB
      END DO

    END IF

    V(I,I) = 1.0_JWRB
    G      = RV1(I)
    L      = I

  END DO

!-----------------------------------------------------------------------
! Accumulate left-hand transformation matrix U.
! U is returned in A.
!-----------------------------------------------------------------------

  DO I = MIN(M,N), 1, -1

    L = I + 1
    G = W(I)

    IF (I < N) THEN
      DO J = L, N
        A(I,J) = 0.0_JWRB
      END DO
    END IF

    IF (G /= 0.0_JWRB) THEN

      G = 1.0_JWRB / G

      IF (I < N) THEN

        DO J = L, N

          S = 0.0_JWRB

          DO K = L, M
            S = S + A(K,I) * A(K,J)
          END DO

          F = (S / A(I,I)) * G

          DO K = I, M
            A(K,J) = A(K,J) + F * A(K,I)
          END DO

        END DO

      END IF

      DO J = I, M
        A(J,I) = A(J,I) * G
      END DO

    ELSE

      DO J = I, M
        A(J,I) = 0.0_JWRB
      END DO

    END IF

    A(I,I) = A(I,I) + 1.0_JWRB

  END DO

!-----------------------------------------------------------------------
! Diagonalise the bidiagonal matrix.
!-----------------------------------------------------------------------

  DO K = N, 1, -1

    DO ITS = 1, MAX_ITERATIONS

      CANCELLATION_REQUIRED = .FALSE.
      L = K

!-----------------------------------------------------------------------
! Locate the splitting point.
!-----------------------------------------------------------------------

      DO

        NM = L - 1

        IF (ABS(RV1(L)) + ANORM == ANORM) THEN
          EXIT
        END IF

        IF (NM >= 1) THEN
          IF (ABS(W(NM)) + ANORM == ANORM) THEN
            CANCELLATION_REQUIRED = .TRUE.
            EXIT
          END IF
        END IF

        L = L - 1

      END DO

!-----------------------------------------------------------------------
! Cancel RV1(L), if required.
!-----------------------------------------------------------------------

      IF (CANCELLATION_REQUIRED) THEN

        C = 0.0_JWRB
        S = 1.0_JWRB

        DO I = L, K

          F      = S * RV1(I)
          RV1(I) = C * RV1(I)

          IF (ABS(F) + ANORM == ANORM) EXIT

          G    = W(I)
          H    = PYTHAG(F,G)
          W(I) = H

          IF (H /= 0.0_JWRB) THEN
            H = 1.0_JWRB / H
          END IF

          C =  G * H
          S = -F * H

          DO J = 1, M
            Y       = A(J,NM)
            Z       = A(J,I)
            A(J,NM) = Y * C + Z * S
            A(J,I)  = -Y * S + Z * C
          END DO

        END DO

      END IF

      Z = W(K)

!-----------------------------------------------------------------------
! Singular value K has converged.
!-----------------------------------------------------------------------

      IF (L == K) THEN

        IF (Z < 0.0_JWRB) THEN

          W(K) = -Z

          DO J = 1, N
            V(J,K) = -V(J,K)
          END DO

        END IF

        EXIT

      END IF

!-----------------------------------------------------------------------
! No convergence.
!
! The old implementation used PAUSE here. Instead, return the index of
! the singular value that failed to converge. The caller can then decide
! how the failure should be handled without terminating ecWAM.
!-----------------------------------------------------------------------

      IF (ITS == MAX_ITERATIONS) THEN
        KERR = K
        RETURN
      END IF

!-----------------------------------------------------------------------
! Shift from the bottom 2 x 2 minor.
!-----------------------------------------------------------------------

      X  = W(L)
      NM = K - 1
      Y  = W(NM)
      G  = RV1(NM)
      H  = RV1(K)

      F = ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / &
          (2.0_JWRB * H * Y)

      G = PYTHAG(F,1.0_JWRB)

      F = ((X-Z)*(X+Z) + H*(Y/(F+SIGN(G,F))-H)) / X

!-----------------------------------------------------------------------
! QR transformation.
!-----------------------------------------------------------------------

      C = 1.0_JWRB
      S = 1.0_JWRB

      DO J = L, NM

        I = J + 1
        G = RV1(I)
        Y = W(I)
        H = S * G
        G = C * G

        Z = PYTHAG(F,H)

        RV1(J) = Z

        IF (Z /= 0.0_JWRB) THEN
          C = F / Z
          S = H / Z
        END IF

        F = X * C + G * S
        G = -X * S + G * C
        H = Y * S
        Y = Y * C

        DO JJ = 1, N
          X       = V(JJ,J)
          Z       = V(JJ,I)
          V(JJ,J) = X * C + Z * S
          V(JJ,I) = -X * S + Z * C
        END DO

        Z    = PYTHAG(F,H)
        W(J) = Z

        IF (Z /= 0.0_JWRB) THEN
          Z = 1.0_JWRB / Z
          C = F * Z
          S = H * Z
        END IF

        F = C * G + S * Y
        X = -S * G + C * Y

        DO JJ = 1, M
          Y       = A(JJ,J)
          Z       = A(JJ,I)
          A(JJ,J) = Y * C + Z * S
          A(JJ,I) = -Y * S + Z * C
        END DO

      END DO

      RV1(L) = 0.0_JWRB
      RV1(K) = F
      W(K)   = X

    END DO

  END DO

END SUBROUTINE SVDCMP


!=======================================================================
! Compute SQRT(A**2+B**2) while reducing the possibility of numerical
! overflow or underflow.
!=======================================================================

FUNCTION PYTHAG(A,B) RESULT(HYPOT)

  IMPLICIT NONE

  REAL(KIND=JWRB), INTENT(IN) :: A
  REAL(KIND=JWRB), INTENT(IN) :: B

  REAL(KIND=JWRB) :: HYPOT
  REAL(KIND=JWRB) :: ABSA
  REAL(KIND=JWRB) :: ABSB

  ABSA = ABS(A)
  ABSB = ABS(B)

  IF (ABSA > ABSB) THEN

    HYPOT = ABSA * SQRT(1.0_JWRB + (ABSB/ABSA)**2)

  ELSE IF (ABSB > 0.0_JWRB) THEN

    HYPOT = ABSB * SQRT(1.0_JWRB + (ABSA/ABSB)**2)

  ELSE

    HYPOT = 0.0_JWRB

  END IF

END FUNCTION PYTHAG

END MODULE YOWLSQ
