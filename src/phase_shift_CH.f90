SUBROUTINE Phase_shift_Ch(Index,ek,eta)
  !     
  !     EQ(20) Chadan et al JMP 42 (2001) 4031
  !     
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_ho_f
  !
  IMPLICIT NONE
  !
  !      
  INTEGER(KIND = I4B), INTENT(IN) :: Index
  REAL(KIND = DP), INTENT(IN)  :: ek
  !
  REAL(KIND = DP), INTENT(OUT) :: eta
  !    
  !     
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: AUX
  !REAL(KIND = DP) :: X,VPOT,WF2,WFP2
  REAL(KIND = DP) :: temp, Error
  INTEGER(KIND = DP) :: Ifail, Ierr
  !  INTEGER(KIND = DP) :: I
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
       REAL(KIND = DP) :: Potf
     END FUNCTION Potf
     !
  END INTERFACE Potf
  !
  ALLOCATE(AUX(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "AUX allocation request denied."
     STOP
  ENDIF
  !
  ! NUMERATOR
  !
  ! ZERO TO INFINITY 
  AUX = Potf(X_GRID)*AVEC_HAR_X(:,Index)**2/(AVEC_HAR_DER_X(:,Index)**2 + ek*ek*AVEC_HAR_X(:,Index)**2)
  AUX(1:dim_X/2) = 0.0_DP
!  DO I = dim_X/2 + 1, dim_X
!     X = X_GRID(I)
!     VPOT = Potf(X)
!     WF2 = AVEC_HAR_X(I,Index)**2
!     WFP2 = AVEC_HAR_DER_X(I,Index)**2
!     AUX(I) = VPOT*WF2/(WFP2+ek*ek*WF2)
!  ENDDO
  !     
  Ifail = 0
  CALL D01GAF(X_Grid, AUX, dim_X, temp, Error, Ifail)
  !
  eta = -ek*temp
  !
  RETURN
  !
END SUBROUTINE Phase_shift_Ch
