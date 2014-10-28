SUBROUTINE E1(ndim, dim_X, X_max, X_grid, aval, avec_X, nstates)
  !
  ! SUM_i |<AVEC(b)|X|AVEC(i)>|^2 = <AVEC(b)|X^2|AVEC(2)> - <AVEC(b)|X|AVEC(2)>^2
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
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, n
  REAL(KIND = DP) :: error, E1_numerical, Total_E1, Test_value, res_1, res_2
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, TEST_element_1, TEST_element_2 
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
  ALLOCATE(TEST_element_1(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_1 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(TEST_element_2(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_2 allocation request denied."
     STOP
  ENDIF
  !
  PRINT*, "  "
  PRINT*, "CALCULATING E1"
  !
  DO n = 1, nstates
     !
     Total_E1 = 0.0_DP
     !
     ! TO PLOT E1 AS A FUNCTION OF ENERGY
     file = 'E1_from_'
     WRITE(filename, '(A,I1,"_N",I2, ".dat")') TRIM(file), n-1, ndim
     IF ( ndim < 10) THEN !to avoid spaces
        WRITE(filename, '(A,I1,"_N",I1, ".dat")') TRIM(file), n-1, ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename, '(A,I1,"_N",I3, ".dat")') TRIM(file), n-1, ndim
     ENDIF
     OPEN(UNIT = 78, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(78,*) "# E1 from ", n-1,"-th state  dim_bas = ", ndim, " Box radius = ", X_max, " fm"
     WRITE(78,*) "#Aval(i)     (E1_numerical)**2"
     !
     ! TO PLOT TOT_E1 AS A FUNCTION OF ENERGY
     WRITE(filename_tot, '(A,I1,"_TOT_N",I2, ".dat")') TRIM(file), n-1, ndim
     IF ( ndim < 10) THEN 
        WRITE(filename_tot, '(A,I1,"_TOT_N",I1, ".dat")') TRIM(file), n-1, ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename_tot, '(A,I1,"_TOT_N",I3, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     OPEN(UNIT = 88, FILE = filename_tot, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(88,*) "# TOT E1 from ", n-1,"-th state   dim_bas = ", ndim, " Box radius = ", X_max, " fm"
     WRITE(88,*) "#Aval(i)     Tot_E1"
     !
     DO i = 1, ndim - 2 !to avoid ISQW divergence
        !
        matrix_element = Avec_x(:,n)*X_grid(:)*Avec_x(:,i)  
        !
        E1_numerical = 0.0_DP
        Ifail = 0
        CALL D01GAF(X_Grid, matrix_element, dim_X, E1_numerical, error, Ifail)
        !
        Total_E1 = Total_E1 + (E1_numerical)**2
        !
        ! SAVING E1 and Total_E1
        WRITE(78,11)  Aval(i), (E1_numerical)**2
        WRITE(88,11)  Aval(i), Total_E1
     ENDDO
     CLOSE(UNIT = 78)
     CLOSE(UNIT = 88)
     !
     ! EVALUATE m_0
     TEST_element_1 = Avec_x(:,n)*X_grid(:)*X_grid(:)*Avec_x(:,n)  
     TEST_element_2 = Avec_x(:,n)*X_grid(:)*Avec_x(:,n)  
     !
     res_1 = 0.0_DP
     res_2 = 0.0_DP
     !
     Ifail = 0
     CALL D01GAF(X_Grid, TEST_element_1, dim_X, res_1, error, Ifail)
     Ifail = 0
     CALL D01GAF(X_Grid, TEST_element_2, dim_X, res_2, error, Ifail)
     !
     Test_value = res_1 - res_2**2
     !
     ! PRINT RESULTS
     PRINT*, "Total E1 from ", n-1, "-th state: ", Total_E1, " <-------> Test value: ", Test_value
     !
  ENDDO
  !
11 FORMAT (1X,E16.8,1X,E17.8)
  !
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(TEST_element_1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_1 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(TEST_element_2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "TEST_element_2 deallocation request denied."
     STOP
  ENDIF
  !
  RETURN
  !
END SUBROUTINE E1
