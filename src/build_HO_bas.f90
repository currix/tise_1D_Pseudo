SUBROUTINE HO_1D_BASIS(apar, NdimH, Iprint)
  !     
  !     COMPUTES A 1D HARMONIC BASIS
  !
  !     INPUT  :: X_GRID  --> VECTOR WITH X GRID
  !               apar    --> LENGTH SCALE OF THE PROBLEM
  !               NdimH   --> DIMENSION + 1 OF THE HARMONIC BASIS
  !
  !
  !     OUTPUT :: HARBAS  --> MATRIX WITH HARMONIC BASIS
  !
  !     FORMAT :: Iprint  --> VERBOSITY CONTROL
  !
  !     $Id: build_HO_bas.f90,v 1.4 2013/06/13 10:53:11 curro Exp $
  !
  !     by Currix TM.
  !
  USE nrtype
  USE egs_ho_f
  !
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: NdimH, Iprint
  REAL(KIND = DP), INTENT(IN) :: apar 
  !
  !
  ! OTHER VARIABLES
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: HO_norm_test
  REAL(KIND = DP) :: PI14, apar2, ainteg, Error
  INTEGER(KIND = I4B) :: kx, Ifail, Ierr
  !
  IF (Iprint > 2) PRINT*, "BUILDING HARMONIC BASIS"
  !
  PI14 = SQRT(SQRT(PI_D))
  apar2 = apar*apar
  !
  !
  !     HO n = 0
  HAR_BAS(:,1) = SQRT(apar)/(PI14)*EXP(-apar2*X_GRID*X_GRID/2.0_DP)
  !
  !     HO n = 1
  IF (NdimH > 1) HAR_BAS(:,2) = (SQRT(apar)/(PI14))*SQRT(2.0_DP)*apar*X_GRID*EXP(-apar2*X_GRID*X_GRID/2.0_DP)
  !
  !    RECURRENCE RELATION (WATCH OUT THE INDEXES :: MATRIX START AT 1, NOT 0)
  DO kx = 2, NdimH-1
     HAR_BAS(:,kx+1) = &
          SQRT(2.0_DP/(1.0_DP*kx))*apar*X_GRID*HAR_BAS(:,kx) - &
          SQRT((1.0_DP*(kx-1))/(1.0_DP*kx))*HAR_BAS(:,kx-1)
  ENDDO
  !
  !     TESTING NORMALIZATION 
  IF (Iprint > 1) THEN
     !
     ! DEFINE HO_norm_test
     ALLOCATE(HO_norm_test(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "HO_norm_test allocation request denied."
        STOP
     ENDIF
     !
     DO kx = 1, NdimH
        HO_norm_test = HAR_BAS(:,kx)*HAR_BAS(:,kx)
        Ifail = 0
        !
        CALL D01GAF(X_GRID, HO_norm_test, dim_X, ainteg, Error, Ifail)
        PRINT*, "HARMONIC FUNCTION ", kx, " NORMALIZATION", ainteg
     ENDDO
     DEALLOCATE(HO_norm_test, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "HO_norm_test deallocation request denied."
        STOP
     ENDIF
     !
     IF (Iprint > 2) WRITE(*,*) "DONE"
     !
  ENDIF
  !     
  !     
END SUBROUTINE HO_1D_BASIS
