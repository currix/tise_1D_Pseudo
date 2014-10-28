SUBROUTINE E2(ndim, dim_X, X_max, X_grid, aval, avec_X, nstates)
  !
  ! SUM_i/=b |<AVEC(b)|X^2|AVEC(i)>|^2 = <AVEC(b)|Y^2|AVEC(b)> - <AVEC(b)|Y|AVEC(b)>^2 
  ! where Y = X^2 - <AVEC(b)|X^2|AVEC(b)>
  ! by LauPK
  !
  USE constants
  USE nrtype
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  REAL(KIND = DP), INTENT(IN) :: X_max
  INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, nstates
  REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, aval
  REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, j, n
  REAL(KIND = DP) :: error, E2_numerical, Total_E2, Test_value, res_scale
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, TEST_element, Rescaled_X_Grid, scale_vector
  !
  CHARACTER(LEN=65) :: filename, filename_tot
  CHARACTER(LEN=65) :: file
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(TEST_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Rescaled_X_Grid(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Rescaled_X_Grid allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(scale_vector(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "scale_vector allocation request denied."
     STOP
  ENDIF
  !
  PRINT*, "  "
  PRINT*, "CALCULATING E2"
  !
  DO n = 1, nstates
     !
     ! TO PLOT E1 AS A FUNCTION OF ENERGY
     file = 'E2_from_'
     WRITE(filename, '(A,I1,"_N",I2, ".dat")')   TRIM(file), n-1, ndim
     IF ( ndim < 10) THEN !to avoid spaces
        WRITE(filename, '(A,I1,"_N",I1, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename, '(A,I1,"_N",I3, ".dat")')    TRIM(file), n-1, ndim
     ENDIF
     OPEN(UNIT = 79, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(79,*) "# E2 from ", n-1, "-th state  dim_bas = ", ndim, " Box radius = ", X_max, " fm"
     WRITE(79,*) "# Aval(i), (E2_analytical)**2"
     !
     ! TO PLOT TOT_E1 AS A FUNCTION OF ENERGY
     WRITE(filename_tot, '(A,I1,"_TOT_N",I2, ".dat")')  TRIM(file), n-1,  ndim
     IF ( ndim < 10) THEN 
        WRITE(filename_tot, '(A,I1,"_TOT_N",I1, ".dat")')  TRIM(file), n-1,  ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename_tot, '(A,I1,"_TOT_N",I3, ".dat")')   TRIM(file), n-1,  ndim
     ENDIF
     OPEN(UNIT = 88, FILE = filename_tot, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(88,*) "# TOT E2 from ", n-1, "-th state  dim_bas = ", ndim, " Box radius = ", X_max, " fm"
     WRITE(88,*) "#Aval(i)     Tot_E2"
     !
     Total_E2 = 0.0_DP
     DO i = 1, ndim -2 ! to avoid ISQW divergence
        !     
        IF(i == n) CYCLE 
        !
        matrix_element = Avec_x(:,n)*X_grid(:)*X_grid(:)*Avec_x(:,i)  
        !
        E2_numerical = 0.0_DP
        Ifail = 0
        CALL D01GAF(X_Grid, matrix_element, dim_X, E2_numerical, error, Ifail)
        !
        Total_E2 = Total_E2 + (E2_numerical)**2
        !
        ! SAVING E1 and Total_E1
        WRITE(79,11)  Aval(i), (E2_numerical)**2
        WRITE(88,11)  Aval(i), Total_E2
     ENDDO
     CLOSE(UNIT = 78)
     CLOSE(UNIT = 88)
     !
     !
     ! OPERATOR SCALING
     scale_vector = Avec_x(:,n)*(X_grid(:)**2)*Avec_x(:,n)
     res_scale = 0.0_DP
     Ifail = 0
     CALL D01GAF(X_Grid, scale_vector, dim_X, res_scale, error, Ifail)
     !
     DO j = 1, dim_X
        Rescaled_X_Grid(j) = X_grid(j)**2 - res_scale
     ENDDO
     !
     ! EVALUATE m_0
     !
     TEST_element = Avec_x(:,n)*(Rescaled_X_Grid(:)**2)*Avec_x(:,n)  
     !
     Test_value = 0.0_DP
     !
     Ifail = 0
     CALL D01GAF(X_Grid, TEST_element, dim_X, Test_value, error, Ifail)
     !
     ! PRINT RESULTS
     PRINT*, "Total E2 from ", n-1, "-th state: ", Total_E2, " <-------> Test value: ", Test_value
     !
  ENDDO
  !
11 FORMAT (1X,E16.8,1X,E17.8)
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(TEST_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_1 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Rescaled_X_Grid, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Rescaled_X_Grid deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(scale_vector, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "scale_vector deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE E2
