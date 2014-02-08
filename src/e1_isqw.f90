SUBROUTINE E1_ISQW(Iprint, I_toten)
  !
  ! <GS|X|AVEC(i)>
  !
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
  REAL(KIND = DP) :: error, E1_numerical, aj, ak
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: matrix_element, Total_E1
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  matrix_x,  E1_analytical
  !
  CHARACTER(LEN=65) :: filename_TM, filename_E1
  CHARACTER(LEN=65) :: file = 'E1'
  CHARACTER(LEN=56) :: prog = 'isqw'
  !
  IF (Iprint > 2) WRITE(*,*), "CALCULATING E1"
  !
  ALLOCATE(matrix_element(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_element allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(matrix_x(1:DIM_BOX,1:DIM_BOX), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Total_E1(1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_E1  allocation request denied."
     STOP
  ENDIF
  !
  Total_E1 = 0.0_DP
  !
  ALLOCATE(E1_analytical(1:dim_BOX,1:I_toten), STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "E1_analytical  allocation request denied."
     STOP
  ENDIF
  !
  !
  ! Compute X matrix in the ISQW basis
  Total_E1 = 0.0_DP
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
  !
  ! Numerical method
  DO i_state = 1, I_toten
     !
     WRITE(*,*), " ", i_state,"-th state :: Numerical method"
     !
     DO i = 1, DIM_BOX
        !     
        matrix_element = Avec_box_x(:,i_state)*X_grid(:)*Avec_box_x(:,i)  
        !
        E1_numerical = 0.0_DP
        Ifail = 0
        !
        CALL D01GAF(X_Grid, matrix_element, dim_X, E1_numerical, error, Ifail)
        !
     WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <",i_state,"| X |Avec(",i,")> = ", E1_numerical
     ENDDO
     !
     !
     WRITE(*,*), "Analytical method"
     !
     !
     WRITE(filename_E1, '(A, "_",A,"_N",I2,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename_E1, '(A, "_",A,"_N",I1,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ELSE IF ( dim_BOX > 99) THEN 
        WRITE(filename_E1, '(A, "_",A,"_N",I3,"_",I1,".dat")') TRIM(prog), TRIM(file), dim_BOX, i_state
     ENDIF
     !
     OPEN(UNIT = 78, FILE = filename_E1, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     WRITE(78,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
     WRITE(78,*) "#Aval_box(i)    E1_analytical**2"
     !
     !
     !
     DO i = 1, DIM_BOX
        E1_analytical(i, i_state) = DOT_PRODUCT(Avec_box(:,i_state),MATMUL(matrix_x,Avec_box(:,i)))
        !
        WRITE(*,15), i, "-th state energy: ", Aval_box(i), " <",i_state,"| X |Avec(",i,")> = ", E1_analytical(i, i_state)
        ! SAVING E1
        WRITE(78,11)  Aval_box(i), E1_analytical(i, i_state)**2
        !
        Total_E1(i_state) = Total_E1(i_state) + (E1_analytical(i, i_state)**2)
        !
     ENDDO
     !
     CLOSE(UNIT = 78)
     !
     WRITE(*,*), "Total E1: ", i_state, Total_E1(i_state)
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
  OPEN(UNIT = 76, FILE = filename_TM, STATUS = "UNKNOWN", ACTION = "WRITE")
  !
  WRITE(76,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
  WRITE(76,*) "#Aval_box(i)    E1_analytical"
  !
  DO i = 1, dim_BOX     !
     ! SAVING E1
     WRITE(76,12)  Aval_box(i), E1_analytical(i,1:I_toten)
     !
  ENDDO
  !
  CLOSE(UNIT = 76)
  !  
  !
11 FORMAT (1X,E16.8,1X,E17.8)
12 FORMAT (1X,E16.8,1X,10E17.8) !!!! Take care of the number of bound states I_toten
15 FORMAT (2X,I3,A,E16.8,2X,A,I3,A,I3,A,E17.8)
  !
  DEALLOCATE(matrix_x, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "matrix_x deallocation request denied."
     STOP
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
  DEALLOCATE(Total_E1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "Total_E1  deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(E1_analytical, STAT = Ierr)
  IF (Ierr /= 0) THEN
     WRITE(*,*), "E1_analytical deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE E1_ISQW
