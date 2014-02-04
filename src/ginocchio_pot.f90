FUNCTION Potf(x)
  !
  !     GINOCCHIO 1D POTENTIAL
  !
  USE nrtype
  USE constants
  USE pot_param
  !
  IMPLICIT NONE
  !
  INTERFACE r
     ELEMENTAL FUNCTION r(y,l)
       !
       USE nrtype
       USE constants
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL(KIND = DP), INTENT(IN) :: y, l
       !
       REAL(KIND = DP) ::  r
     END FUNCTION r
  END INTERFACE r
  !
  ! ARGUMENTS
  REAL(KIND = DP), INTENT(IN) :: x
  !
  REAL(KIND = DP) ::  Potf
  REAL(KIND = DP) :: y, y0, y1, ym, tol, nu, lambda, f0, f1, fm, s0, s1, sm, V0
  !
  !POTENTIAL PARAMETERS
  nu = PARAM_POT(1) 
  lambda = PARAM_POT(2)
  V0 = PARAM_POT(3) 
  !
  !tolerance
  tol = 0.001_DP
  !
  !Initial values
  y0 = -1.0_DP
  f0 = r(y0,lambda) - X
  check_y0: IF (f0 == 0.0_DP) THEN !checks if y0 is a solution
     Potf = 0.0_DP
     RETURN
  ENDIF check_y0
  S0 = abs(F0)/F0
  !
  y1 = 1.0_DP
  f1 =  r(y1,lambda) - X
  check_y1 : IF(F1 == 0.0_DP) THEN !checks if y1 is a solution
     Potf = 0.0_DP
     RETURN
  ENDIF check_y1
  S1 = abs(F1)/F1
  !
  bisection_method : DO
     !middle point
     ym = (y0 + y1)/2.0_DP
     fm =  r(ym,lambda) - X
     check_ym : IF(Fm == 0.0) THEN !checks if ym is  a solution
        y = ym
        Potf = V0*(-(lambda**2)*nu*(nu+1.0_DP)*(1.0_DP-y**2) &
             + ((1.0_DP-y**2)*(1.0_DP-lambda**2)/4.0_DP)* &
             (2.0_DP - (7.0_DP-lambda**2)*(y**2) + 5.0_DP*(1.0_DP-lambda**2)*(y**4) ) )
        RETURN
     ENDIF check_ym
     Sm = abs(fm)/fm
     !
     IF( S0 /= Sm) THEN
        y1 = ym
        F1 =  r(y1,lambda) - X
        S1 = abs(F1)/F1
     ELSE
        y0 = ym
        F0 =  r(y0,lambda) - X
        S0 =  abs(F0)/F0
     ENDIF
     !
     IF (ABS(y0-y1) < tol )THEN
        y = y0
        Potf = V0*(-(lambda**2)*nu*(nu+1.0_DP)*(1.0_DP-y**2) &
             + ((1.0_DP-y**2)*(1.0_DP-lambda**2)/4.0_DP)* &
             (2.0_DP - (7.0_DP-lambda**2)*(y**2) + 5.0_DP*(1.0_DP-lambda**2)*(y**4) ))
        RETURN
     ENDIF
  ENDDO bisection_method
  !
END FUNCTION Potf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELEMENTAL FUNCTION r(y,l)
   !
  USE nrtype
  USE constants
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  REAL(KIND = DP), INTENT(IN) :: y, l
  !
  REAL(KIND = DP) ::  r
  REAL(KIND = DP) :: aux 
  aux = SQRT(l**2 -1.0_DP)
  !
  IF( y >= 1.0_DP) THEN
     r = (l**(-2))* (100.0_DP + aux*ATAN(aux*y) )
     RETURN
  ELSE IF (y <= -1.0_DP) THEN
     r = (l**(-2))* (  -100.0_DP + aux*ATAN(aux*y) )
     RETURN
  ELSE
     r = (l**(-2))* ( ATANH(y) + aux*ATAN(aux*y) )
     RETURN
  ENDIF
  !
END FUNCTION r

