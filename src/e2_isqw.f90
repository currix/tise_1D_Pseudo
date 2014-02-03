SUBROUTINE E2_ISQW(Iprint, I_toten)
  !
  ! <GS|X**2|AVEC(i)>
  !
  ! $Id: e2_isqw.f90,v 1.6 2013/10/05 17:56:24 curro Exp curro $
  !
  USE constants
  USE nrtype
  USE egs_mod_isqw
  !
  IMPLICIT NONE
  !
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail, i, k, j, i_state
  REAL(KIND = DP) :: error, E2_numerical
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_E2
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_x2,  E2_analytical
  !
  CHARACTER(LEN=65) :: filename_E2, filename_TM
  CHARACTER(LEN=65) :: file = 'E2'
  CHARACTER(LEN=56) :: prog = 'isqw'
  !
  !
  IF (Iprint > 2) WRITE(*,*), "CALCULATING E2"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matrix_x2(1:DIM_BOX,1:DIM_BOX), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_x2 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_E2(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_E2  allocation request denied."
     STOP
  ENDIF
  !
  Total_E2 = 0.0_DP
  !
  ALLOCATE(E2_analytical(1:dim_BOX,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "E2_analytical  allocation request denied."
     STOP
  ENDIF
  !
  !
  ! Compute X**2 matrix in the ISQW basis
  matrix_x2 = 0.0_DP
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
  DO i_state = 1, I_toten
     !
     WRITE(*,*), " ", i_state,"-th state :: Numerical method"
     !
     DO i = 1, DIM_BOX
        !     
        matrix_element = Avec_box_x(:,i_state)*X_grid(:)*X_grid(:)*Avec_box_x(:,i)  
        !
        E2_numerical = 0.0_DP
        Ifail = 0
        !
        CALL D01GAF(X_Grid, matrix_element, dim_X, E2_numerical, error, Ifail)
        !
        WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <GS | X^2 |Avec(",i,")> = ", E2_numerical
     ENDDO
     !
     !
     WRITE(*,*), "Analytical method"
     !
     Total_E2 = 0.0_DP
     !
     WRITE(filename_E2, '(A, "_",A,"_N",I2,"_",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename_E2, '(A, "_",A,"_N",I1,"_",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE IF ( dim_BOX > 99) THEN 
        WRITE(filename_E2, '(A, "_",A,"_N",I3,"_",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ENDIF
     !
     OPEN(UNIT = 79, FILE = filename_E2, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     WRITE(79,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
     WRITE(79,*) "# Aval_box(i)     E2_analytical**2"
     !
     DO i = 1, DIM_BOX
        E2_analytical(i, i_state) = DOT_PRODUCT(Avec_box(:,i_state), MATMUL(matrix_x2,Avec_box(:,i)))
        WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <bnd | X^2 |Avec(",i,")> = ",  E2_analytical(i, i_state)
        ! SAVING E2
        WRITE(79,11)  Aval_box(i), (E2_analytical(i, i_state))**2
        !
        Total_E2(i_state) = Total_E2(i_state) + (E2_analytical(i, i_state))**2
     ENDDO
     !
     CLOSE(UNIT = 79)
     !
     WRITE(*,*), "Total E2: ", i_state, Total_E2(i_state)
     !
  ENDDO
  !
  WRITE(filename_TM, '(A, "_",A,"_TM_N",I2,".dat")') TRIM(prog), TRIM(file), dim_BOX
  IF ( dim_BOX < 10) THEN !to avoid spaces
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ELSE IF ( dim_BOX > 99) THEN 
     WRITE(filename_TM, '(A, "_",A,"_TM_N",I3,".dat")') TRIM(prog), TRIM(file), dim_BOX
  ENDIF
  !
  !
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_box(i)    E2_analytical"
  !
  DO i = 1, dim_BOX     
     ! SAVING E2
     WRITE(76,12)  Aval_box(i), E2_analytical(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I3,A,E16.8,2X,A,I2,A,E17.8)
  !
  DEALLOCATE(matrix_x2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_x2 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(matrix_element, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element deallocation request denied."
     STOP
  ENDIF
  !
  !
  DEALLOCATE(Total_E2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_E2 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(E2_analytical, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "E2_analytical deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE E2_ISQW
