ELEMENTAL FUNCTION Potf(x)
  !
  !     MORSE 1D POTENTIAL
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
  !     PARAM_POT(1) --> VP0 DEPTH (MeV)
  !     PARAM_POT(2) --> a inverse of potential range (fm⁻1)
  !
  Potf = Param_pot(1) * ( EXP(-2.0_DP*PARAM_POT(2)*X) - 2.0_DP*EXP(-PARAM_POT(2)*X)  )
  !
  !
END FUNCTION Potf
