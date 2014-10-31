SUBROUTINE B1_THO(Iprint, I_toten, B_met1, B_met2)
  !
  ! <GS|X|AVEC(i)> THO pseudostate basis
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
  REAL(KIND = DP) :: error, B1_num1,  B1_num2, integ
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, matel, Total_B1
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_X_THO_BASIS, B1_matrix
  !
  CHARACTER(LEN=65) :: filename_TM, filename_E1
  CHARACTER(LEN=65) :: file = 'B1'
  CHARACTER(LEN=56) :: prog = 'tho'
  !
  IF (Iprint > 1) PRINT*, "CALCULATING B1"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matel(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matel allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_B1(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Total_B1  allocation request denied."
     STOP
  ENDIF
  !
  Total_B1 = 0.0_DP
  !
  !
  ALLOCATE(B1_matrix(1:dim_THO,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "B1_matrix  allocation request denied."
     STOP
  ENDIF
  !
  !
  ! X matrix
  !
  IF (B_met2) THEN
     !
     ALLOCATE(matrix_X_THO_BASIS(1:dim_THO,1:dim_THO), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "matrix_X_THO_BASIS allocation request denied."
        STOP
     ENDIF
     !
     DO j = 1, dim_THO
        DO k = 1, dim_THO
           !
           integ = 0.0_DP
           Ifail = 0
           !
           IF (j /= k) THEN
              matel = X_Grid(:)*THO_Bas(:,k)*THO_Bas(:,j)*der_S_x(:)
              !
              CALL D01GAF(X_Grid, matel, dim_X, integ, error, Ifail)
           ENDIF
           matrix_X_THO_BASIS(k,j) = integ
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  !Main loop
  !
  DO i_state = 1, I_toten
     !
     ! Define output filename
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename_E1, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ELSE IF ( dim_THO < 100) THEN 
        WRITE(filename_E1, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ELSE
        WRITE(filename_E1, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_THO, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_E1, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     !
     IF (B_met1) THEN
        !
        IF (Iprint > 1) PRINT*, "B1 ", i_state,"-th state :: Numerical method - 1"
        !
        ! File header
        WRITE(78,*) "# THO  dim_THO = ", dim_THO, " Integ radius = ", X_max, " fm"
        WRITE(78,*) "#Aval_THO(i)    B1_met1**2"
        !
        DO i = 1, dim_THO
           !
           matrix_element = 0.0_DP
           !     
           matrix_element = Avec_THO_x(:,i_state)*X_grid(:)*Avec_THO_x(:,i)  
           !
           B1_num1 = 0.0_DP
           Ifail = 0
           !
           CALL D01GAF(X_Grid, matrix_element, dim_X, B1_num1, error, Ifail)
           !
           IF (Iprint > 1) WRITE(*,15), i, "-th state energy: ", Aval_THO(i), " <",I_state," | X |Avec(",i,")> = ", B1_num1
           !
           B1_matrix(i, i_state) =  B1_num1
           !
           !
           ! SAVING B1
           WRITE(78,11)  Aval_THO(i), B1_num1**2
           !
           Total_B1(i_state) = Total_B1(i_state) + (B1_num1**2)
           !
         ENDDO
        !
         PRINT*, "Total B1: ", i_state, Total_B1(i_state)
        !
     ENDIF
     !
     IF (B_met2) THEN
        !
        IF (Iprint > 1) PRINT*, "B1 ", i_state,"-th state :: Numerical method - 2"
        !
        ! File header
        WRITE(78,*) "# THO  dim_THO = ", dim_THO, " Integ radius = ", X_max, " fm"
        WRITE(78,*) "#Aval_THO(i)    B1_met2**2"
        !
        DO i = 1, dim_THO
           B1_num2 = DOT_PRODUCT(Avec_THO(:,i_state),MATMUL(matrix_X_THO_BASIS, Avec_THO(:,i)))
           IF (Iprint > 1) WRITE(*,15), i, "-th state energy: ", Aval_THO(i), " <",I_state," | X |Avec(",i,")> = ", B1_num2
           !
           B1_matrix(i, i_state) =  B1_num2
           !
           ! SAVING B1
           WRITE(78,11)  Aval_THO(i), B1_num2**2
           !
           Total_B1(i_state) = Total_B1(i_state) + (B1_num2**2)
           !
        ENDDO
        !
        PRINT*, "Total B1: ", i_state, Total_B1(i_state)
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
  WRITE(76,*) "# THO  dim_THO = ", dim_THO, " Integ. radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_tho(i)    B1(i,1:n_states)"
  !
  DO i = 1, dim_THO
     !
     ! SAVING B1
     WRITE(76,12)  Aval_THO(i), B1_matrix(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !
  !
  !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (1X,E16.8,1X,10E17.8)
15 FORMAT (2X,I2,A,E16.8,2X,A,I2,A,I2,A,E17.8)
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
  DEALLOCATE(Total_B1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Total_B1  deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(B1_matrix, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "B1_matrix deallocation request denied."
     STOP
  ENDIF
  !
  ! 
  IF (B_met2) THEN
     DEALLOCATE(matrix_X_THO_BASIS, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "matrix_element deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE B1_THO
