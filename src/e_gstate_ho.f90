SUBROUTINE EGS(KVAL,FVAL)
  !
  ! $Id: e_gstate_ho.f90,v 1.5 2013/05/22 16:12:43 curro Exp $
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_ho_f
  !
  IMPLICIT NONE
  !
  REAL(KIND = DP), INTENT(IN)  :: KVAL
  REAL(KIND = DP), INTENT(OUT) :: FVAL
  !
  REAL(KIND = DP) :: APAR
  !
  INTERFACE HO_1D_BASIS
     !
     SUBROUTINE HO_1D_BASIS(APAR, NDIMH, IPRINT)
       !
       USE nrtype
       USE egs_ho_f
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: NDIMH, IPRINT
       REAL(KIND = DP), INTENT(IN) :: APAR 
       !
     END SUBROUTINE HO_1D_BASIS
     !
  END INTERFACE HO_1D_BASIS
  !
  INTERFACE HARDIAG
     !
     SUBROUTINE HARDIAG(APT, IPRINT, IFLAG)
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_ho_f
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
     END SUBROUTINE HARDIAG
     !
  END INTERFACE HARDIAG
  !
  APAR = SQRT(SQRT(KVAL/H_SQ_OVER_M))
  !
  CALL HO_1D_BASIS(APAR, 1, 1)
  !
  CALL HARDIAG(APAR, 0, 1)
  !
  FVAL = AVAL_HAR(1)
  !
  RETURN
END SUBROUTINE EGS
