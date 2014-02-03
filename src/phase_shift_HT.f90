SUBROUTINE Phase_shift_HT(dim_X, X_grid, avec_X, lambda, Index, ee, ek, eta_1, eta_2)
  !
  !     EQ(13) Hazi & Taylor PRA 1 ( 1970) 1109
  !
  !
  ! $Id: phase_shift_HT.f90,v 1.3 2013/05/05 22:27:25 curro Exp $
  !
  ! by Currix TM.
  !
  USE nrtype
  USE constants
  USE pot_param
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: dim_X, Index
  REAL(KIND = DP), INTENT(IN) :: lambda, ee, ek
  REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid
  REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
  !
  REAL(KIND = DP), INTENT(OUT) :: eta_1, eta_2
  !
  !
  INTERFACE F_FUNCTION
     ELEMENTAL FUNCTION F(lambda,X)
       USE nrtype
       IMPLICIT NONE
       REAL(KIND = DP), INTENT(IN) :: lambda, X
       REAL(KIND = DP) ::  F
       !
     END FUNCTION F
  END INTERFACE F_FUNCTION
  !
  INTERFACE FP_FUNCTION
     ELEMENTAL FUNCTION FP(lambda,X)
       USE nrtype
       IMPLICIT NONE
       REAL(KIND = DP), INTENT(IN) :: lambda, X
       REAL(KIND = DP) ::  FP
     END FUNCTION FP
  END INTERFACE FP_FUNCTION
  !
  INTERFACE FPP_FUNCTION
     ELEMENTAL FUNCTION FPP(lambda,X)
       USE nrtype
       IMPLICIT NONE
       REAL(KIND = DP), INTENT(IN) :: lambda, X
       REAL(KIND = DP) ::  FPP
     END FUNCTION FPP
  END INTERFACE FPP_FUNCTION
  !     
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE ::  aux
  REAL(KIND = DP) ::  ainteg1, ainteg2, error
  !
  INTEGER(KIND = I4B) :: I, Ifail, IERR
  !
  !
  INTERFACE Potf
     !
     ELEMENTAL FUNCTION Potf(X)
       !
       !     WOODS-SAXON 1D POTENTIAL
       !
       USE nrtype
       USE constants
       USE pot_param
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL(KIND = DP), INTENT(IN) :: X
       REAL(KIND = DP) ::  Potf
     END FUNCTION Potf
     !
  END INTERFACE Potf
  !
  !     ZERO TO INFINITY
  !     NUMERATOR
  ALLOCATE(aux(1:dim_X), STAT = IERR)    
  IF (IERR /= 0) THEN
     PRINT*, "aux allocation request denied."
     STOP
  ENDIF
  !
  !
  DO I = dim_X/2 + 1, dim_X
     aux(I) = avec_X(I,Index)*ee*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) + &
          H_SQ_OVER_M*0.5_DP*avec_X(I,Index)*( &
          FPP(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) + &
          2.0_DP*ek*FP(lambda,X_Grid(I))*COS(ek*X_Grid(I))- &
          ek*ek*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) &
          ) - &
          avec_X(I,Index)*Potf(X_Grid(I))*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I))
  ENDDO
  !
  Ifail = 0
  CALL D01GAF(X_Grid(dim_X/2+1), aux(dim_X/2+1), dim_X/2, ainteg1, error, Ifail)
  !
  !     DENOMINATOR
  DO I = dim_X/2 + 1, dim_X
     aux(I) = avec_X(I,Index)*ee*F(lambda,X_Grid(I))*COS(ek*X_Grid(I)) + &
          H_SQ_OVER_M*0.5_DP*avec_X(I,Index)*( &
          FPP(lambda,X_Grid(I))*COS(ek*X_Grid(I)) - &
          2.0_DP*ek*FP(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) - &
          ek*ek*F(lambda,X_Grid(I))*COS(ek*X_Grid(I)) &
          ) - &
          avec_X(I,Index)*Potf(X_Grid(I))*F(lambda,X_Grid(I))*COS(ek*X_Grid(I))
  ENDDO
  !
  Ifail = 0
  CALL D01GAF(X_Grid(dim_X/2+1), aux(dim_X/2+1), dim_X/2, ainteg2, error, Ifail)
  !
  eta_1 = ATAN(-ainteg1/ainteg2)
  !     -INFINITY TO ZERO
  !     NUMERATOR
  DO I = 1, dim_X/2 + 1
     aux(I) = avec_X(I,Index)*ee*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) + &
          H_SQ_OVER_M*0.5_DP*avec_X(I,Index)*( &
          FPP(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) + &
          2.0_DP*ek*FP(lambda,X_Grid(I))*COS(ek*X_Grid(I)) - &
          ek*ek*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) &
          ) - &
          avec_X(I,Index)*Potf(X_Grid(I))*F(lambda,X_Grid(I))*SIN(ek*X_Grid(I))
  ENDDO
  !
  Ifail = 0
  CALL D01GAF(X_Grid, aux, dim_X/2 + 1, ainteg1, error, Ifail)
  !
  !     DENOMINATOR
  DO I = 1, dim_X/2 + 1
     aux(I) = avec_X(I,Index)*ee*F(lambda,X_Grid(I))*COS(ek*X_Grid(I)) + &
          H_SQ_OVER_M*0.5_DP*avec_X(I,Index)*( &
          FPP(lambda,X_Grid(I))*COS(ek*X_Grid(I)) - &
          2.0_DP*ek*FP(lambda,X_Grid(I))*SIN(ek*X_Grid(I)) - &
          ek*ek*F(lambda,X_Grid(I))*COS(ek*X_Grid(I)) & 
          ) - &
          avec_X(I,Index)*Potf(X_Grid(I))*F(lambda,X_Grid(I))*COS(ek*X_Grid(I))
  ENDDO
  !
  Ifail = 0
  CALL D01GAF(X_Grid, aux, dim_X/2 + 1, ainteg2, error, Ifail)
  !
  eta_2 = ATAN(ainteg1/ainteg2)
  !
  DEALLOCATE(aux, STAT = IERR)    
  IF (IERR /= 0) THEN
     PRINT*, "aux deallocation request denied."
     STOP
  ENDIF
  !
  RETURN
  !
END SUBROUTINE PHASE_SHIFT_HT
!
ELEMENTAL FUNCTION F(lambda,X)
  USE nrtype
  IMPLICIT NONE
  REAL(KIND = DP), INTENT(IN) :: lambda, X
  REAL(KIND = DP) :: F
  F = 1.0D0 - EXP(-lambda*X*X)
  RETURN
END FUNCTION F
!
ELEMENTAL FUNCTION FP(lambda,X)
  USE nrtype
  IMPLICIT NONE
  REAL(KIND = DP), INTENT(IN) :: lambda, X
  REAL(KIND = DP) :: FP
  FP = 2.0D0*lambda*X*EXP(-lambda*X*X)
  RETURN
END FUNCTION FP
!
ELEMENTAL FUNCTION FPP(lambda,X)
  USE nrtype
  IMPLICIT NONE
  REAL(KIND = DP), INTENT(IN) :: lambda, X
  REAL(KIND = DP) :: FPP
  FPP = 2.0D0*lambda*(1.0D0-2.0D0*lambda*X*X)*EXP(-lambda*X*X)
  RETURN
END FUNCTION FPP
