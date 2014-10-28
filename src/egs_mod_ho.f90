MODULE egs_ho_f
  !
  ! $Id: egs_mod_ho.f90,v 1.3 2013/05/05 09:29:14 curro Exp $
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !
  ! NDX -> XGRID DIMENSION
  INTEGER(KIND = I4B) :: dim_X
  ! HARDIM -> DIMENSION OF THE HO BASIS
  INTEGER(KIND = I4B) :: dim_HO
  !
  ! Xmin Xmax and delta_X step
  REAL(KIND = DP) :: X_min, X_max, Delta_X
  ! XGRID
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: X_Grid      
  ! HARMONIC BASIS, EIGENVALUES AND EIGENVECTORS
  REAL(KIND = DP),  DIMENSION(:), ALLOCATABLE :: Aval_Har
  REAL(KIND = DP),  DIMENSION(:,:), ALLOCATABLE :: Har_Bas, Avec_Har, Avec_Har_X
  ! EIGENVECTOR DERIVATIVES
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Avec_Har_Der_X
  !
  ! Last bound state (for optimization of the oscillator length)
  INTEGER(KIND = I4B) :: last_bound_state
  !
END MODULE egs_ho_f
