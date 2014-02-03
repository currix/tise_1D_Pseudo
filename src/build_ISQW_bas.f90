SUBROUTINE ISQW_1D_BASIS(XMpar, Iprint)
  !     
  !     COMPUTES THE 1D BOX BASIS
  !
  !     INPUT     :: XMPAR  --> PI/(2 X_M)
  !    
  !     VARIABLES :: DIM_X  --> XGRID DIMENSION
  !                  X_GRID --> VECTOR WITH X GRID
  !                  X_MAX  --> BOX RADIUS
  !                  NDIM   --> DIMENSION OF THE BOX BASIS
  !
  !     OUTPUT    :: BOXBAS --> MATRIX WITH BOX BASIS
  !
  !     FORMAT    :: IPRINT --> VERBOSITY CONTROL
  !
  !     by Currix TM.
  !
  USE nrtype
  USE egs_mod_isqw
  !
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint
  REAL(KIND = DP), INTENT(IN) ::  XMpar 
  !
  !
  ! OTHER VARIABLES
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: ISQW_norm_test
  REAL(KIND = DP) :: xmi, error, ainteg
  INTEGER(KIND = I4B) :: kx, ind, Ifail, Ierr
  !
  IF (Iprint > 2) PRINT*, "BUILDING BOX BASIS"
  !
  xmi = 1.0_dp/SQRT(X_max)
  !
  DO kx = 1, dim_X
     DO ind = 1, dim_BOX, 2
        Box_Bas(kx,ind) = xmi*COS(ind*XMpar*X_grid(kx))
     ENDDO
     DO ind = 2, dim_BOX, 2
        Box_Bas(kx,ind) = xmi*SIN(ind*XMpar*X_grid(kx))
     ENDDO
  ENDDO
  !
  !    TESTING NORMALIZATION 
  IF (Iprint > 1) THEN
     !
     ! DEFINE ISQW_norm_test
     ALLOCATE(ISQW_norm_test(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "ISQW_norm_test allocation request denied."
        STOP
     ENDIF
     !
     DO kx = 1, dim_BOX
        ISQW_norm_test = Box_Bas(:,kx)*Box_Bas(:,kx)
        Ifail = 0
        !
        CALL D01GAF(X_grid, ISQW_norm_test, dim_X, ainteg, error, Ifail)
        PRINT*, "BOX FUNCTION ", kx, " NORMALIZATION", ainteg
     ENDDO
     !
     DEALLOCATE(ISQW_norm_test, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "ISQW_norm_test deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !     
  IF (Iprint > 2) WRITE(*,*) "DONE"
  !     
END SUBROUTINE ISQW_1D_BASIS
