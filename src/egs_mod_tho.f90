MODULE egs_tho_f
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !
  ! NDX -> XGRID DIMENSION
  INTEGER(KIND = I4B) :: dim_X
  ! HARDIM -> DIMENSION OF THE HO BASIS
  INTEGER(KIND = I4B) :: dim_THO, dim_THO_diag
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_Grid 
  !
  ! Xmin Xmax and delta_X step
  REAL(KIND = DP) :: X_min, X_max, Delta_X
  !     
  ! TRANSFORMED HARMONIC OSCILLATOR BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_THO
  REAL(KIND = DP)  :: max_aval_tho ! Maximum value of energies considered
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE ::  THO_Bas, Avec_THO, Avec_THO_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_THO_Der_X
  !
  !
  ! Last bound state (for optimization of the oscillator length)
  INTEGER(KIND = I4B) :: last_bound_state
  !
END MODULE egs_tho_f
