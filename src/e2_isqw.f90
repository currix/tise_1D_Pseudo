SUBROUTINE B2_ISQW(Iprint, I_toten, B_numerical, B_analytical)
  !
  ! < bound | X**2 | AVEC(i) >
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
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, k, j, i_state
  REAL(KIND = DP) :: error, B2_analytical
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_B2
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_x2,  B2_matrix
  !
  CHARACTER(LEN=65) :: filename_B2, filename_TM
  CHARACTER(LEN=65) :: file = 'B2'
  CHARACTER(LEN=56) :: prog = 'isqw'
  !
  !
  IF (Iprint > 2) WRITE(*,*), "CALCULATING B2"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_B2(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_B2  allocation request denied."
     STOP
  ENDIF
  !
  Total_B2 = 0.0_DP
  !
  !
  ALLOCATE(B2_matrix(1:dim_BOX,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "B2_matrix  allocation request denied."
     STOP
  ENDIF
  !
  ! Compute X**2 matrix in the ISQW basis
  IF (B_analytical) THEN
     !
     ALLOCATE(matrix_x2(1:DIM_BOX,1:DIM_BOX), STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*), "matrix_x2 allocation request denied."
        STOP
     ENDIF
     !
     matrix_x2 = 0.0_DP
     !
     DO j = 1, DIM_BOX
        DO k = 1, DIM_BOX  ! can make just half matrix and consider it symmetric!!!
           ! Equal parity case
           IF (MOD(j,2) == MOD(k,2)) THEN
              !Same parity case (odd-odd)
              IF ( MOD(j,2) /= 0.0_dp ) THEN  
                 IF ( j == k ) THEN
                    matrix_x2(k,j) = 2.0_DP * &
                         ( &
                         1.0_DP/6.0_DP - 1.0_DP/((PI_D*REAL(j,DP))**2) &
                         )    
                 ELSE 
                    matrix_x2(k,j) = ( 8.0_dp /(PI_D**2)) * &
                         ( &
                         ((-1.0_DP)**((k-j)/2)) / (REAL(k-j,DP)**2)  &
                         + &
                         ((-1.0_DP)**((j+k)/2)) / (REAL(j+k,DP)**2) &
                         )
                 ENDIF
                 !Same parity case (even-even)
              ELSE 
                 IF ( j == k ) THEN
                    matrix_x2(k,j) = 2.0_DP * &
                         ( &
                         1.0_DP/6.0_DP - 1.0_DP/((PI_D*REAL(j,DP))**2) &
                         )            
                 ELSE
                    matrix_x2(k,j) = ( 8.0_DP/(PI_D**2)) * &
                         ( &
                         ((-1.0_DP)**((k-j)/2)) / (REAL(k-j,DP)**2) &
                         - &
                         ((-1.0_DP)**((j+k)/2)) / (REAL(j+k,DP)**2) &
                         )
                 ENDIF
              ENDIF
           !
           ENDIF
        ENDDO
     ENDDO
     !
     matrix_x2 = matrix_x2*(X_max**2)
     !
  ENDIF
  !
  ! Main Loop
  DO i_state = 1, I_toten
     !
     ! Define output filename
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename_B2, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename_B2, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE 
        WRITE(filename_B2, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ENDIF
     !
     !
     OPEN(UNIT = 78, FILE = filename_B2, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     !
     IF (B_numerical) THEN
        !        
        IF (Iprint > 0) WRITE(*,*), "B2 :: ", i_state,"-th state :: Numerical method"
        !
        ! File header
        WRITE(78,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i    k_i    B1_i_analytical**2"
        !
        DO i = 1, DIM_BOX_diag
           !     
           matrix_element = Avec_box_x(:,i_state)*X_grid(:)*X_grid(:)*Avec_box_x(:,i)  
           !
           B2_matrix(i, i_state) = 0.0_DP
           Ifail = 0
           !
           CALL D01GAF(X_Grid, matrix_element, dim_X, B2_matrix(i, i_state), error, Ifail)
           !
           IF (Iprint > 0) &
                WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <",i_state,"| X^2 |Avec(",i,")> = ", B2_matrix(i, i_state)
           !
           Total_B2(i_state) = Total_B2(i_state) + (B2_matrix(i, i_state)**2)
           !
           ! SAVING B2
           WRITE(78,11)  Aval_box(i), (SIGN(Aval_box(i),Aval_box(i))/ABS(Aval_box(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_box(i))/h_sq_over_m), B2_matrix(i, i_state)**2
           !
        ENDDO
        !
        IF (Iprint > 0) WRITE(*,*) "Total B2: ", i_state, Total_B2(i_state)
        !
     ENDIF
     !
     !
     IF (B_analytical) THEN 
        !
        !
        IF (Iprint > 0) WRITE(*,*) "B2 :: Analytical method, state ", i_state
        !
        ! File header
        WRITE(78,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
        WRITE(78,*) "#   E_i    k_i    B1_i_analytical**2"
        !
        DO i = 1, DIM_BOX_diag
           !
           B2_analytical = DOT_PRODUCT(Avec_box(:,i_state), MATMUL(matrix_x2,Avec_box(:,i)))
           !
           B2_matrix(i, i_state) = B2_analytical
           !
           IF (Iprint > 0) &
                WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <", i_state,"| X^2 |Avec(",i,")> = ",  B2_analytical
           !
           ! SAVING B2
           WRITE(78,11)  Aval_box(i), (SIGN(Aval_box(i),Aval_box(i))/ABS(Aval_box(i)))* & ! To consider bound states
                SQRT(2.0_DP*ABS(Aval_box(i))/h_sq_over_m), B2_matrix(i, i_state)**2
           !
           Total_B2(i_state) = Total_B2(i_state) + B2_matrix(i, i_state)**2
           !
        ENDDO
        !
        !
        IF (Iprint > 0) WRITE(*,*) "Total B2: ", i_state, Total_B2(i_state)
        !
     ENDIF
     !
     CLOSE(UNIT = 78)
     !
     !
  ENDDO
  !
  IF ( dim_BOX < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ELSE IF ( dim_BOX < 100) THEN 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ELSE 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ENDIF
  !
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_box(i)    B2_numerical"
  !
  DO i = 1, dim_BOX     
     ! SAVING B2
     WRITE(76,12)  Aval_box(i), B2_matrix(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !
11 FORMAT (1X,E16.8,1X,E16.8,1X,E18.9)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I3,A,E16.8,2X,A,I3,A,I3,A,E17.8)
  !
  IF (B_analytical) THEN
     !
     DEALLOCATE(matrix_x2, STAT = Ierr)
     IF (Ierr /= 0) THEN
        WRITE(*,*), "matrix_x2 deallocation request denied."
        STOP
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  !
  DEALLOCATE(Total_B2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_E2 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(B2_matrix, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "E2_numerical deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE B2_ISQW
