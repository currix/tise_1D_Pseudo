ELEMENTAL FUNCTION Potf(x)
  !
  !     MEXICAN HAT 1D POTENTIAL
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
  !     PARAM_POT(1) --> VP0 DEPTH (MeV fm^-2)
  !     PARAM_POT(2) --> to regulate the range (MeVfm^-4)
  !
  !	THE TRUNCATED MEXICAN HAT
  !
  !IF ( ABS(X) < 2.0_DP) THEN
  !	Potf = Param_pot(1) * X**2 + PARAM_POT(2)*X**4
  !ELSE
  !	Potf = 0.0_DP
  !ENDIF
  !
  !	OR TWO SYMMETRICAL OSCILLATOR
  IF ( 1.0_DP <= X .AND. X <= 5.0_DP) THEN
        Potf = Param_pot(1) - PARAM_POT(2)*X + PARAM_POT(3)*X**2
  ELSE IF ( -5.0_DP <= X .AND. X <= -1.0_DP) THEN
        Potf = Param_pot(1) + PARAM_POT(2)*X + PARAM_POT(3)*X**2
  ELSE
        Potf = 0.0_DP
  ENDIF
  !
END FUNCTION Potf
