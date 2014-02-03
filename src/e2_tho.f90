SUBROUTINE E2_THO(Iprint, I_toten)
  !
  ! <GS|X^2|AVEC(i)>
  !
  ! $Id: e2_tho.f90,v 1.5 2013/10/05 17:48:34 curro Exp $
  !
  USE constants
  USE nrtype
  USE egs_tho_f
  USE lst_param
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, k, j, i_state
  REAL(KIND = DP) :: error, E2_numerical,  E2_num2, integ, Total_E2
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, matel
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_X_THO_BASIS
  !
  CHARACTER(LEN=65) :: filename_TM, filename_E2
  CHARACTER(LEN=65) :: file
  CHARACTER(LEN=56) :: prog = 'tho'
  !
  IF (Iprint > 2) PRINT*, "CALCULATING E2"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matrix_X_THO_BASIS(1:dim_THO,1:dim_THO), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_X_THO_BASIS allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matel(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matel allocation request denied."
     STOP
  ENDIF
  !
  ! X**2 matrix
  matrix_X_THO_BASIS = 0.0_DP
  DO j = 1, dim_THO
     DO k = 1, dim_THO
        !
        integ = 0.0_DP
        Ifail = 0
        !
        !        IF (j /= k) THEN
        matel =(X_Grid(:)**2)*THO_Bas(:,k)*THO_Bas(:,j)*der_S_x(:)
        !
        CALL D01GAF(X_Grid, matel, dim_X, integ, error, Ifail)
        !        ENDIF
        matrix_X_THO_BASIS(k,j) = integ
     ENDDO
  ENDDO
  !
  DO i_state = 1, I_toten
     !
     IF (Iprint > 0) PRINT*, "E2 ", i_state,"-th state :: Numerical method - 1"
     !
     DO i = 1, dim_THO
        !     
        matrix_element = Avec_THO_x(:,i_state)*X_grid(:)*X_grid(:)*Avec_THO_x(:,i)  
        !
        E2_numerical = 0.0_DP
        Ifail = 0
        !
        CALL D01GAF(X_Grid, matrix_element, dim_X, E2_numerical, error, Ifail)
        !
        PRINT 12, i, "-th state energy: ", Aval_THO(i), " <state | X^2 |Avec(",i,")> = ", E2_numerical
     ENDDO
     !
     !
     PRINT*, "Numerical method - 2"
     !
     !
     Total_E2 = 0.0_DP
     file = 'E2'
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     WRITE(filename_E2, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
        WRITE(filename_E2, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ELSE IF ( dim_THO > 99) THEN 
        WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
        WRITE(filename_E2, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ENDIF
     !
     OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
     OPEN(UNIT = 78, FILE = filename_E2, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(76,*) "# THO  E2 dim_THO = ", dim_THO, " Integ radius = ", X_max, " fm"
     WRITE(76,*) "#Aval_THO(i)    E2_analytical"
     WRITE(78,*) "# THO  dim_THO = ", dim_THO, " Integ radius = ", X_max, " fm"
     WRITE(78,*) "#Aval_THO(i)    E2_analytical**2"
     !
     DO i = 1, dim_THO
        E2_num2 = DOT_PRODUCT(Avec_THO(:,i_state),MATMUL(matrix_X_THO_BASIS, Avec_THO(:,i)))
        PRINT 12, i, "-th state energy: ", Aval_THO(i), " <state | X^2 |Avec(",i,")> = ", E2_num2
        ! SAVING E2
        WRITE(76,11)  Aval_THO(i), E2_num2
        WRITE(78,11)  Aval_THO(i), E2_num2**2
        Total_E2 = Total_E2 + (E2_num2**2)
     ENDDO
     !
     CLOSE(UNIT = 76)
     CLOSE(UNIT = 78)
     !
     PRINT*, "Total E2: ", Total_E2
     !
  ENDDO
  !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (2X,I2,A,E16.8,2X,A,I2,A,E17.8)
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(matel, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matel deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(matrix_X_THO_BASIS, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE E2_THO
