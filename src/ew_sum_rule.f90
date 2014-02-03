SUBROUTINE Ew_Strength(ndim, dim_X, X_grid, aval, avec_X, iprint)
  !
  ! $Id: ew_sum_rule.f90,v 1.3 2013/05/05 20:31:56 curro Exp $
  !
  !     by Currix & Laura M.
  !
  USE nrtype
  USE constants
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, Iprint
  REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, aval
  REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
  !
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail
  INTEGER(KIND = I4B) :: KI
  REAL(KIND = DP) ::  Ews_temp, Ews_total, Error
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: x2_gstate
   !
  IF (Iprint > 1) WRITE(*,*) "CALCULATING E-W STRENGTH"
  !    
  ALLOCATE(x2_gstate(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "X2_GSTATE allocation request denied."
     STOP
  ENDIF
  !
  Ews_total = 0.0_DP
  !
  DO KI = 1, ndim
     x2_gstate = avec_X(:,1)*X_Grid(:)*avec_X(:,KI)
     Ifail = 0
     !
     Ews_temp = 0.0_DP
     !
     CALL D01GAF(X_Grid, x2_gstate, dim_X, Ews_temp, Error, Ifail)
     !     only to avoid an optimization error
     if (ki.eq.ndim+1) write(*,*) ki,Ews_temp,Ews_total
     Ews_total = Ews_total + Ews_temp*Ews_temp*(aval(KI) - aval(1))
  ENDDO
  !
  WRITE(*,*) "E-W SUM RULE E_W = ", Ews_total/h_sq_over_m
  !
  RETURN
  !
END SUBROUTINE EW_STRENGTH
