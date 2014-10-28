SUBROUTINE Total_Strength(ndim, dim_X, X_grid, avec_X, Aval, iprint)
  !
  !  ndim   :: pseudostate basis dimension
  !  dim_X  :: X grid dimension
  !  X_grid :: X grid
  !  avec_X   :: eigenstates array (X dependence)
  !  iprint :: verbosity control
  !
  !  $Id: total_strength.f90,v 1.5 2013/05/23 08:07:32 laura Exp $
  !
  !     by Currix TM.
  !
  USE nrtype
  USE constants
  !
  IMPLICIT NONE
  !
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, Iprint
  REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, Aval
  REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
  !
  INTEGER(KIND = I4B) :: Ifail, Ierr
  INTEGER(KIND = I4B) :: ki, i
  REAL(KIND = DP) ::  Total_Str_0, Total_Str_Temp, Total_Str, error
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: x2_vec
  !
  IF (Iprint > 1) PRINT*, "CALCULATING TOTAL STRENGTH"
  !
  ALLOCATE(x2_vec(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "x2_vec allocation request denied."
     STOP
  ENDIF
  !
  DO i = 1, ndim
     IF (Aval(i) >= 0.0_DP) EXIT
     !
     X2_vec = (avec_X(:,i)*X_grid(:))**2
     !
     !
     Ifail = 0
     Total_Str_0 = 0.0_DP
     !
     CALL D01GAF(X_Grid, x2_vec, dim_X, Total_Str_0, error, Ifail)
     !
     WRITE(*,11) "TOTAL STRENGTH ", i, "-th bound state: <",i, "| X^2 |",i,"> = ", Total_Str_0
  ENDDO
11 FORMAT (A,i2,A,i2,A,i2,A,E14.6)
  !
  Total_Str = 0.0_DP
  DO ki = 1, ndim
     !
     x2_vec = avec_X(:,1)*X_Grid(:)*avec_X(:,ki)
     !
     Ifail = 0
     Total_Str_Temp = 0.0_DP
     !
     CALL D01GAF(X_Grid, x2_vec, dim_X, Total_Str_Temp, error, Ifail)
     !     only to avoid an optimization error
     if (ki == ndim+1) write(123,*) ki, Total_Str_Temp, Total_Str
     Total_Str = Total_Str + Total_Str_Temp*Total_Str_Temp
  ENDDO
  !
  WRITE(*,*) "SUM RULE S_T = ", Total_Str
  !
  DEALLOCATE(x2_vec, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "x2_vec deallocation request denied."
     STOP
  ENDIF
  !
  RETURN
  !
END SUBROUTINE TOTAL_STRENGTH
