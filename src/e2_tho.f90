SUBROUTINE B2_THO(Iprint, I_toten, B_met1, B_met2)
  !
  !    < bound | X^2 | AVEC(i) > THO pseudostate basis
  !
  !
  USE constants
  USE nrtype
  USE egs_tho_f
  USE lst_param
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  LOGICAL, INTENT(IN) :: B_met1, B_met2
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, k, j, i_state
  REAL(KIND = DP) :: error, B2_num1, integ
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, matel, Total_B2
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_X2_THO_BASIS, B2_matrix
  !
  CHARACTER(LEN=65) :: filename_TM, filename_B2
  CHARACTER(LEN=65) :: file = "B2"
  CHARACTER(LEN=56) :: prog = 'tho'
  !
  IF (Iprint > 1) WRITE(*,*) "CALCULATING B2"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matel(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "matel allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_B2(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "Total_E1  allocation request denied."
     STOP
  ENDIF
  !
  Total_B2 = 0.0_DP
  !
  !
  ALLOCATE(B2_matrix(1:dim_THO,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "B2_matrix  allocation request denied."
     STOP
  ENDIF
  !
  ! Build the < x^2 > (matrix_x2_THO_BASIS) in the THO basis
  !
  IF (B_met2) THEN
     !
     ALLOCATE(matrix_X2_THO_BASIS(1:dim_THO,1:dim_THO), STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*) "matrix_X2_THO_BASIS allocation request denied."
        STOP
     ENDIF
     !
     matrix_X2_THO_BASIS = 0.0_DP
     !
     ! To do : Matrix is symmetric, do not compute each element.
     DO j = 1, dim_THO
        DO k = 1, dim_THO
           !
           integ = 0.0_DP
           Ifail = 0
           !
           matel =(X_Grid(:)**2)*THO_Bas(:,k)*THO_Bas(:,j)*der_S_x(:)
           !
           CALL D01GAF(X_Grid, matel, dim_X, integ, error, Ifail)
           !
           matrix_X2_THO_BASIS(k,j) = integ
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! Main loop
  !
  DO i_state = 1, I_toten
     !
     ! Define output filename
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename_B2, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ELSE IF ( dim_THO < 100) THEN 
        WRITE(filename_B2, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ELSE
        WRITE(filename_B2, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_B2, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     !
     !
     !
     IF (B_met1) THEN
        !
        IF (Iprint > 1) WRITE(*,*) "B2 ", i_state,"-th state :: Numerical method - 1"
        !
        !
        ! File header
        WRITE(78,*) "# THO  dim_THO = ", dim_THO, "dim_THO_diag = ", dim_THO_diag, " Integ radius = ", X_max, " fm"
        WRITE(78,*) "# Aval_THO(i)    B2_met1**2"
        !
        DO i = 1, dim_THO_diag
           !     
           matrix_element = Avec_THO_x(:,i_state)*X_grid(:)*X_grid(:)*Avec_THO_x(:,i)  
           !
           B2_num1 = 0.0_DP
           Ifail = 0
           !
           CALL D01GAF(X_Grid, matrix_element, dim_X, B2_num1, error, Ifail)
           !
           IF (Iprint > 1) WRITE(*,15), i, "-th state energy: ", Aval_THO(i), " <",i_state," | X^2 |Avec(",i,")> = ", B2_num1
           !
           !
           B2_matrix(i, i_state) =  B2_num1
           !
           !
           ! SAVING B2
           WRITE(78,11)  Aval_THO(i), B2_num1**2
           !
           Total_B2(i_state) = Total_B2(i_state) + (B2_num1**2)
           !
        ENDDO
        !
        WRITE(*,*) "Total B2: ", i_state, Total_B2(i_state)
        !
     ENDIF
     !
     !
     IF (B_met2) THEN
        !
        IF (Iprint > 1) WRITE(*,*), "Numerical method - 2"
        !
        ! File header
        WRITE(78,*) "# THO  dim_THO = ", dim_THO, "dim_THO_diag = ", dim_THO_diag, " Integ radius = ", X_max, " fm"
        WRITE(78,*) "# Aval_THO(i)    B1_met2**2"
        !
        DO i = 1, dim_THO_diag
           !
           B2_matrix(i, i_state) = DOT_PRODUCT(Avec_THO(:,i_state),MATMUL(matrix_X2_THO_BASIS, Avec_THO(:,i)))
           !
           IF (Iprint > 1) &
                WRITE(*,15), i, "-th state energy: ", Aval_THO(i), " <",I_state," | X^2 |Avec(",i,")> = ", B2_matrix(i, i_state)
           !
           ! SAVING B2
           WRITE(78,11)  Aval_THO(i), B2_matrix(i, i_state)**2
           !
           Total_B2(i_state) = Total_B2(i_state) + B2_matrix(i, i_state)**2
           !
        ENDDO
        !
        WRITE(*,*) "Total B2: ", Total_B2(i_state)
        !
     ENDIF
     !
     CLOSE(UNIT = 78)
     !
     !
  ENDDO
  !
  !
 IF (dim_THO < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,".dat")') TRIM(prog), TRIM(file), dim_THO
  ELSE IF ( dim_THO < 100) THEN 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,".dat")') TRIM(prog), TRIM(file), dim_THO
  ELSE 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,".dat")') TRIM(prog), TRIM(file), dim_THO
  ENDIF
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# THO  dim_THO = ", dim_THO, "dim_THO_diag = ", dim_THO_diag, " Integ. radius = ", X_max, " fm"
  WRITE(76,*) "# Aval_tho(i)    B2_num(1:n_states)"
  !
  DO i = 1, dim_THO_diag
     !
     ! SAVING B2
     WRITE(76,12)  Aval_THO(i), B2_matrix(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)!
    !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (1X,E16.8,1X,10E17.8)
15 FORMAT (2X,I2,A,E16.8,2X,A,I2,A,I2,A,E17.8)
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(matel, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "matel deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Total_B2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "Total_E2  deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(B2_matrix, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*) "B2_matrix deallocation request denied."
     STOP
  ENDIF
  !
  IF (B_met2) THEN
     DEALLOCATE(matrix_X2_THO_BASIS, STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*) "matrix_X2_THO_BASIS deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE B2_THO
