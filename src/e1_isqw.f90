SUBROUTINE B1_ISQW(Iprint, I_toten, B_numerical, B_analytical)
  !
  ! < bound | X | AVEC(i) >
  !
  !
  USE constants
  USE nrtype
  USE egs_mod_isqw
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  LOGICAL, INTENT(IN) :: B_numerical, B_analytical
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, k, j, i_state
  REAL(KIND = DP) :: error, B1_analytical, aj, ak
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_B1
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_x,  B1_matrix
  !
  CHARACTER(LEN=65) :: filename_TM, filename_B1
  CHARACTER(LEN=65) :: file = 'B1'
  CHARACTER(LEN=56) :: prog = 'isqw'
  !
  IF (Iprint > 2) WRITE(*,*), "CALCULATING B1"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_B1(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_B1  allocation request denied."
     STOP
  ENDIF
  !
  !
  Total_B1 = 0.0_DP
  !
  ALLOCATE(B1_matrix(1:dim_BOX,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "B1_matrix allocation request denied."
     STOP
  ENDIF
  !
  !
  ! Compute X matrix in the ISQW basis
  !
  IF (B_analytical) THEN
     !
     ALLOCATE(matrix_x(1:DIM_BOX,1:DIM_BOX), STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*), "matrix_x allocation request denied."
        STOP
     ENDIF
     !
     matrix_x = 0.0_DP
     !
     DO j = 1, DIM_BOX
        aj = 1.0_dp*j
        DO k = 1, DIM_BOX
           ak = 1.0_dp*k
           ! Different parity case
           IF ( (MOD(j,2)==0.0_dp .and. MOD(k,2)/=0.0_dp) .or. (MOD(j,2)/=0.0_dp .and. MOD(k,2)==0.0_dp)) THEN
              matrix_x(k,j) = ((4.0_dp*X_max)/(PI_D**2))*(((-1)**((j-k-1)/2))/((aj-ak)**2)  + ((-1)**((j+k-1)/2))/((aj+ak)**2))
           ENDIF
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  ! Main Loop
  DO i_state = 1, I_toten
     !
     ! Define output filename
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename_B1, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename_B1, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE 
        WRITE(filename_B1, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_B1, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     !
     IF (B_numerical) THEN
        !
        IF (Iprint > 0) WRITE(*,*), "B1 :: ", i_state,"-th state :: Numerical method"
        !
        ! File header
        WRITE(78,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i    k_i    B1_i_numerical**2"
        !
        DO i = 1, dim_BOX_diag
           !     
           matrix_element = Avec_box_x(:,i_state)*X_grid(:)*Avec_box_x(:,i)  
           !
           B1_matrix(i, i_state) = 0.0_DP
           Ifail = 0
           !
           CALL D01GAF(X_Grid, matrix_element, dim_X, B1_matrix(i, i_state), error, Ifail)
           !
           IF (Iprint > 0) &
                WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <",i_state,"| X |Avec(",i,")> = ", B1_matrix(i, i_state)
           !
           Total_B1(i_state) = Total_B1(i_state) + (B1_matrix(i, i_state)**2)
           !
           ! SAVING B1
           WRITE(78,11)  Aval_box(i), (SIGN(Aval_box(i),Aval_box(i))/ABS(Aval_box(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_box(i))/h_sq_over_m), B1_matrix(i, i_state)**2
           !
        ENDDO
        !
        IF (Iprint > 0) WRITE(*,*) "Total B1: ", i_state, Total_B1(i_state)
        !
     ENDIF
     !
     !
     IF (B_analytical) THEN 
        !
        IF (Iprint > 0) WRITE(*,*), "B1 :: ", i_state,"-th state :: Analytical method"
        !
        ! File header
        WRITE(78,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i    k_i    B1_i_analytical**2"
        !
        DO i = 1, DIM_BOX_diag
           !
           B1_analytical = DOT_PRODUCT(Avec_box(:,i_state),MATMUL(matrix_x,Avec_box(:,i)))
           !
           IF (Iprint > 0) WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <",i_state,"| X |Avec(",i,")> = ", B1_analytical
           !
           B1_matrix(i, i_state) = B1_analytical
           !
           ! SAVING B1
           WRITE(78,11)  Aval_box(i), (SIGN(Aval_box(i),Aval_box(i))/ABS(Aval_box(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_box(i))/h_sq_over_m), B1_matrix(i, i_state)**2
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
     !
  ENDDO
  !
  !
  IF ( dim_BOX < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ELSE IF ( dim_BOX < 100) THEN 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ELSE
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ENDIF
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_box(i)    B1(i,1:n_states)"
  !
  DO i = 1, dim_BOX     
     !
     ! SAVING B1
     WRITE(76,12)  Aval_box(i), B1_matrix(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !  
  !
11 FORMAT (1X,E16.8,1X,E16.8,1X,E18.9)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I3,A,E16.8,2X,A,I3,A,I3,A,E17.8)
  !
  IF (B_analytical) THEN
     DEALLOCATE(matrix_x, STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*), "matrix_x deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  !
  DEALLOCATE(Total_B1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_B1  deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(B1_matrix, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "B1_matrix deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE B1_ISQW
