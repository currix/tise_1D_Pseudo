PROGRAM ANALYTICAL_GINOCCHIO
  !
  ! PROGRAM THAT CALCULATES ANALYTICAL GINOCCHIO EIGENVALUES AND EIGENVECTORS
  !
  ! $Id:$
  !
  ! by Prototipo 1.0
  !
  USE nrtype
  USE constants
  USE pot_param
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP) :: lambda, nu, reduced_mass
  INTEGER(KIND = I4B) :: n, dim, Ierr
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Energy
  !
  ! DATA INPUT
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_MASS/  reduced_mass
  NAMELIST/INP_DIM/   dim
  NAMELIST/INP_POT/   Param_pot
  !
  !POTENTIAL PARAMETERS
  nu = PARAM_POT(1) 
  lambda = PARAM_POT(2)
  !
  ALLOCATE(Energy(1:dim), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Energy allocation request denied."
     STOP
  ENDIF
  !
DO n = 1, dim
   !
   Energy(n) = -(lambda*(nu + 0.5_DP))**2 - (2.0_DP - lambda**2)*((n + 0.5_DP)**2)&
               + (2.0_DP*n + 1.0_DP)*SQRT( (lambda*(nu + 0.5_DP))**2 + (1.0_DP - lambda**2)*((n + 0.5_DP)**2) )
   !
   PRINT*, "Energy ", n, " -> ", Energy(n)
ENDDO
  !
  DEALLOCATE(Energy, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Energy deallocation request denied."
     STOP
  ENDIF
  !
END PROGRAM ANALYTICAL_GINOCCHIO
