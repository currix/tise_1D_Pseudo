SUBROUTINE E1_HO(Iprint, I_toten, apar)
  !
  ! <GS|X|AVEC(i)>
  !
  USE constants
  USE nrtype
  USE egs_ho_f
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  REAL(KIND = DP), INTENT(IN) :: apar
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, j, i_state
  REAL(KIND = DP) :: error, rj, E1_numerical
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_E1
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: matrix_x,  E1_analytical
  !
  CHARACTER(LEN=65) :: filename_TM, filename_E1
  CHARACTER(LEN=65) :: file = 'E1'
  CHARACTER(LEN=56) :: prog = 'ho'
  !
  IF (Iprint > 2) PRINT*, "CALCULATING E1"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matrix_x(1:dim_HO,1:dim_HO), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_E1(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Total_E1  allocation request denied."
     STOP
  ENDIF
  !
  Total_E1 = 0.0_DP
  !
  ALLOCATE(E1_analytical(1:dim_HO,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "E1_analytical  allocation request denied."
     STOP
  ENDIF
  !
  ! Build x matrix
  !
  matrix_x = 0.0_DP
  !
  DO j = 1, dim_HO - 1
     rj = j*1.0_dp
     matrix_x(j, j + 1) = SQRT(rj/2.0_DP)
  ENDDO
  !
  DO j = 0, dim_HO-2
     rj = j*1.0_dp
     matrix_x(j + 2, j + 1) = SQRT((rj + 1.0_DP)/2.0_DP)
  ENDDO
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Main Loop
  DO i_state = 1, I_toten
     !
     PRINT*, "Numerical method"
     !
     DO i = 1, dim_HO
        !     
        matrix_element = Avec_Har_x(:,i_state)*X_grid(:)*Avec_Har_x(:,i)  
        !
        E1_numerical = 0.0_DP
        Ifail = 0
        !
        CALL D01GAF(X_Grid, matrix_element, dim_X, E1_numerical, error, Ifail)
        !
        WRITE(*, 15), i, "-th state energy: ", Aval_Har(i), " <GS | X |Avec(",i,")> = ", E1_numerical
     ENDDO
     !
     !
     PRINT*, "Analytical method"
     !
     !
     WRITE(filename_E1, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename_E1, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     ELSE IF ( dim_HO > 99) THEN 
        WRITE(filename_E1, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_HO, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_E1, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(78,*) "#   dim_BOX = ", dim_HO, " Integ. radius = ", X_max, " fm"
     WRITE(78,*) "#Aval_har(i)    E1_analytical**2"
     !
     DO i = 1, dim_HO
        E1_analytical(i, i_state) =  DOT_PRODUCT(Avec_Har(:,i_state), MATMUL(matrix_x,Avec_Har(:,i)))
        WRITE(*, 15),  i, "-th state energy: ", Aval_Har(i), " <GS | X |Avec(",i,")> = ", E1_analytical(i, i_state)/apar
        !
        ! SAVING E1
        WRITE(78,11)  Aval_Har(i), (E1_analytical(i, i_state)/apar)**2
        !
        Total_E1(i_state) = Total_E1(i_state) + (E1_analytical(i, i_state)/apar)**2
        !
     ENDDO
     !
     CLOSE(UNIT = 78)
     !
     PRINT*, "Total E1: ", i_state, Total_E1(i_state)
     !
  ENDDO
  !
  WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,".dat")') TRIM(prog), TRIM(file), dim_HO
  IF ( dim_HO < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,".dat")') TRIM(prog), TRIM(file), dim_HO
  ELSE IF ( dim_HO > 99) THEN 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,".dat")') TRIM(prog), TRIM(file), dim_HO
  ENDIF
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# HO  dim_HO = ", dim_HO, " Integ. radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_har(i)    E1_analytical(1:n_states)"
  !
  DO i = 1, dim_HO
     !
     ! SAVING E1
     WRITE(76,12)  Aval_Har(i), E1_analytical(i,1:I_toten)/apar
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I2,A,E16.8,2X,A,I2,A,E17.8)
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(matrix_x, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "matrix_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Total_E1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Total_E1  deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(E1_analytical, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "E1_analytical deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE E1_HO