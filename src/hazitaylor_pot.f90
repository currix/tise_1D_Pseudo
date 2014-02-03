ELEMENTAL FUNCTION Potf(x)
  !
  !     HAZI&TAYLOR 1D POTENTIAL
  !
  USE nrtype
  USE constants
  USE pot_param
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  REAL(KIND = DP), INTENT(IN) :: x
  !
  REAL(KIND = DP) ::  Potf
  !
  !     POTENTIAL PARAMETERS
  !     PARAM_POT(1) --> V0 Potential depth
  !     PARAM_POT(2) --> LAMBDA diffusivity
  !
  IF (x > 0.0_DP) THEN
     POTf = 0.5_dp*PARAM_POT(1)*(x**2)*exp(-PARAM_POT(2)*x**2)
  ELSE
     POTf = 0.5_dp*PARAM_POT(1)*x**2
  ENDIF
  !
  !
END FUNCTION Potf
