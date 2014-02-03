MODULE lst_param
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !THO CONSTANTS
  REAL(KIND = DP) ::  m = 4.0_dp, gamma, b, K_eff, ratio
  !b = oscillator length
  !ratio = gamma/b
  !
  !SCALING FUNCTION
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: S_x, der_S_x, der2_S_x, S_over_x
  !
  !
CONTAINS
  !
  ELEMENTAL FUNCTION lst(x)
    !
    !     ANALYTICAL 1D THO LST
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  lst
    !
    REAL(KIND = DP) :: aux
    !
    aux = ABS(x)**(-m)
    lst = ( aux + SQRT(aux) * (gamma**(-m)) )**(-1.0_DP/REAL(m,DP))
    IF (x < 0.0_DP) lst = - lst
    !
  END FUNCTION lst
  !
!  ELEMENTAL FUNCTION inv(s)
    !
    !     INVERSE FUNCTION
    !
!    USE nrtype
    !
!    IMPLICIT NONE
    !
    ! ARGUMENTS
!    REAL(KIND = DP), INTENT(IN) :: 
    !
!    REAL(KIND = DP) ::  inv
    !
!    REAL(KIND = DP) :: aux
    !
!    aux = ABS(x)**(-m)
!    lst = ( aux + SQRT(aux) * (gamma**(-m)) )**(-1.0_DP/REAL(m,DP))
!    IF (x < 0.0_DP) lst = - lst
    !
!  END FUNCTION inv
  !
  ELEMENTAL FUNCTION sox(x)
    !
    !     ANALYTICAL 1D THO LST s(x)/x
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  sox
    !
    sox = ( 1.0_DP + (gamma**(-m))*SQRT(ABS(x))**m )**(-1.0_DP/REAL(m,DP))
    !
  END FUNCTION sox
  !
  ELEMENTAL FUNCTION d_lst(x)
    !
    !     ANALYTICAL 1D THO LST (1st derivative)
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  d_lst
    !
    d_lst = (sox(x)/2.0_DP)*( 1.0_DP + sox(x)**m )
    !
  END FUNCTION d_lst
  !
  ELEMENTAL FUNCTION d2_lst(x)
    !
    !     ANALYTICAL 1D THO LST (2nd derivative)
    !
    USE nrtype
    !
    IMPLICIT NONE
    !
    ! ARGUMENTS
    REAL(KIND = DP), INTENT(IN) :: x
    !
    REAL(KIND = DP) ::  d2_lst
    !
    REAL(KIND = DP) :: aux1, aux2, vx
    !
    vx = ABS(x)
    aux1 = SQRT(vx)**(m-2)
    aux2 = SQRT(vx)**m
    !
    d2_lst = -(sox(x)/4.0_DP)*(aux1)/(gamma**m + aux2)*( 1.0_DP + (m+1.0_DP)*sox(x)**m )
    IF (x < 0.0_DP) d2_lst = - d2_lst
    !
  END FUNCTION d2_lst
  !
END MODULE lst_param
