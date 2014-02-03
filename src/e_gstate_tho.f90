SUBROUTINE EGS(AVAL,FVAL)
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_tho_f
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP), INTENT(IN)  :: AVAL
  REAL(KIND = DP), INTENT(OUT) :: FVAL
  !
  !
  INTERFACE THO_1D_BASIS
     !
     SUBROUTINE THO_1D_BASIS(APAR, NDIMH, IPRINT)
       !
       USE nrtype
       USE egs_tho_f
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: NDIMH, IPRINT
       REAL(KIND = DP), INTENT(IN) :: APAR 
       !
     END SUBROUTINE THO_1D_BASIS
     !
  END INTERFACE THO_1D_BASIS
  !
  INTERFACE THARDIAG
     !
     SUBROUTINE THARDIAG(APT, IPRINT, IFLAG)
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_tho_f
       !
       ! Lapack 95
       USE LA_PRECISION, ONLY: WP => DP
       USE F95_LAPACK, ONLY: LA_SYEVR
       !
       IMPLICIT NONE
       !
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: IPRINT, IFLAG
       REAL(KIND = DP), INTENT(IN) ::  APT
       !
       !
     END SUBROUTINE THARDIAG
     !
  END INTERFACE THARDIAG
  !
  !
  CALL THO_1D_BASIS(AVAL, 2, 1)
  !
  CALL THARDIAG(AVAL, 0, 1)
  !
  FVAL = AVAL_THO(1)
  !
  RETURN
END SUBROUTINE EGS
