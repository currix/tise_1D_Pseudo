SUBROUTINE B1_HO(Iprint, I_toten, apar, B_analytical, B_numerical)
  !
  ! < bound |X| AVEC(i) >
  !
  USE constants
  USE nrtype
  USE egs_ho_f
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  REAL(KIND = DP), INTENT(IN) :: apar
  LOGICAL, INTENT(IN) :: B_numerical, B_analytical
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, j, i_state
  REAL(KIND = DP) :: error, rj, B1_numerical
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_B1
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: matrix_x,  B1_matrix
  !
  CHARACTER(LEN=65) :: filename_TM, filename_B1
  CHARACTER(LEN=65) :: file = 'B1'
  CHARACTER(LEN=56) :: prog = 'ho'
  !
  !
  IF (Iprint > 2) PRINT*, "CALCULATING B1"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element allocation request denied."
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
  ALLOCATE(B1_matrix(1:dim_HO,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "B1_matrix allocation request denied."
     STOP
  ENDIF
  !
  ! Build x matrix
  !
  IF (B_analytical) THEN
     !
     ALLOCATE(matrix_x(1:dim_HO,1:dim_HO), STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "matrix_x allocation request denied."
        STOP
     ENDIF
     !
     matrix_x = 0.0_DP
     !
     DO j = 1, dim_HO - 1
        rj = j*1.0_dp
        matrix_x(j, j + 1) = SQRT(rj/2.0_DP)/apar
     ENDDO
     !
     DO j = 0, dim_HO-2
        rj = j*1.0_dp
        matrix_x(j + 2, j + 1) = SQRT((rj + 1.0_DP)/2.0_DP)/apar
     ENDDO
     !
  ENDIF
  !
  !
  ! Main Loop
  DO i_state = 1, I_toten
     !
     ! Define output filename
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename_B1, '(A, "_",A,"_N",I1,"_",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     ELSE IF ( dim_HO < 100) THEN !to avoid spaces
        WRITE(filename_B1, '(A, "_",A,"_N",I2,"_",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
        WRITE(filename_B1, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_B1, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     IF (B_numerical) THEN
        !
        IF (Iprint > 0) WRITE(*,*) "B1 :: Numerical method, state ", i_state
        !
        WRITE(78,*) "#   dim_BOX = ", dim_HO, "dim_HO_diag = ", dim_HO_diag, " Integ. radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i   k_i     B1_i,numerical**2"
        !
        DO i = 1, dim_HO_diag
           !     
           matrix_element = Avec_Har_x(:,i_state)*X_grid(:)*Avec_Har_x(:,i)  
           !
           B1_numerical = 0.0_DP
           Ifail = 0
           !
           CALL D01GAF(X_Grid, matrix_element, dim_X, B1_numerical, error, Ifail)
           !
           IF (Iprint > 0) &
                WRITE(*,15), i, "-th state energy: ", Aval_Har(i), " <",i_state,"| X |Avec(",i,")> = ", B1_numerical
           !
           B1_matrix(i, i_state) =  B1_numerical
           !
           ! SAVING B1
           WRITE(78,11)  Aval_Har(i), (SIGN(Aval_Har(i),Aval_Har(i))/ABS(Aval_Har(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_Har(i))/h_sq_over_m), B1_numerical**2
           !
           Total_B1(i_state) = Total_B1(i_state) + B1_numerical**2
           !
        ENDDO
        !
        !
        IF (Iprint > 0) WRITE(*,*) "Total B1: ", i_state, Total_B1(i_state)
        !
     ENDIF
     !
     !
     IF (B_analytical) THEN 
        !
        IF (Iprint > 0) WRITE(*,*) "B1 :: Analytical method, state ", i_state
        !
        ! File header
        WRITE(78,*) "#   dim_BOX = ", dim_HO, "dim_HO_diag = ", dim_HO_diag, " Integ. radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i    k_i    B1_i_analytical**2"
        !
        DO i = 1, dim_HO_diag
           !
           B1_matrix(i, i_state) =  DOT_PRODUCT(Avec_Har(:,i_state), MATMUL(matrix_x,Avec_Har(:,i)))
           !
           IF (Iprint > 0) &
                WRITE(*,15), i, "-th state energy: ", Aval_Har(i), " <",i_state,"| X |Avec(",i,")> = ", B1_matrix(i, i_state)
           !
           ! SAVING B1
           WRITE(78,11)  Aval_Har(i), (SIGN(Aval_Har(i),Aval_Har(i))/ABS(Aval_Har(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_Har(i))/h_sq_over_m), B1_matrix(i, i_state)**2
           !
           Total_B1(i_state) = Total_B1(i_state) + B1_matrix(i, i_state)**2
           !
        ENDDO
        !
        IF (Iprint > 0) WRITE(*,*) "Total B1: ", i_state, Total_B1(i_state)
        !
     ENDIF
     !
     CLOSE(UNIT = 78)
     !
  ENDDO
  !
  IF ( dim_HO < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
  ELSE IF ( dim_HO < 100) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
  ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
  ENDIF
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# HO  dim_HO = ", dim_HO, " Integ. radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_har(i)    B1(i,1:n_states)"
  !
  DO i = 1, dim_HO_diag
     !
     ! SAVING B1
     WRITE(76,12)  Aval_Har(i), B1_matrix(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !
11 FORMAT (1X,E16.8,1X,E16.8,1X,E18.9)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I3,A,E16.8,2X,A,I3,A,I3,A,E17.8)
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Total_B1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Total_E1  deallocation request denied."
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
  IF (B_analytical) THEN
     DEALLOCATE(matrix_x, STAT = Ierr)
     IF (Ierr /= 0) THEN
        PRINT*, "matrix_x deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE B1_HO
