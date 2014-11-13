MODULE egs_mod_isqw
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !
  ! XGRID DIMENSION 
  INTEGER(KIND = I4B) :: dim_X
  ! MIN AND MAX VALUE OF THE XGRID, INTERVAL OF THE GRID
  REAL(KIND = DP) :: X_min, X_max, Delta_X    
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_grid    
  !
  ! 
  ! BOX BASIS DIMENSION
  INTEGER(KIND = I4B) :: dim_BOX, dim_BOX_diag
  ! BOX BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_Box
  REAL(KIND = DP)  :: max_aval_BOX ! Maximum value of energies considered
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  Box_Bas, Avec_Box, Avec_Box_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_Box_Der_X
  !
END MODULE egs_mod_isqw
