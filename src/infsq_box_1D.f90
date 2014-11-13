PROGRAM SQBOX_1BODY_1D
  !
  !
  !     PROGRAM THAT SOLVES NUMERICALLY THE 1D SCHROEDINGER EQUATION 
  !     USING AN INF SQUARE WELL BASIS FOR A GENERAL POTENTIAL.
  !
  !     OUTPUT :: See documentation
  !
  !       BASIS FUNCTIONS ::  X_i F_1(X_i) ... F_dim(X_i)
  !
  !       EIGENFUNCTIONS  ::  X_i EF_1(X_i) ... EF_dim(X_i)
  !
  !       EIGENVALUES     ::  i E_i
  !
  !       EIGENVALUES + EIGENFUNCTIONS AND POTENTIAL
  !                         ::  X_i POT(X_i) 10*EF_1(X_i) + E_1 ... X_dim 10*EF_dim(X_i) + E_dim
  !       EIGENVALUES, EIGENFUNCTIONS SQUARED AND POTENTIAL
  !                         ::  X_i 10*EF_1(X_i)^2 + E_1 ... X_dim 10*EF_dim(X_i)^2 + E_dim
  !     96  Energy (positive values) - phase shift :: E_i \delta_i
  !
  ! by Currix TM and Lau
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_mod_isqw
  !
  IMPLICIT NONE
  !
  !     FLAGS FOR SAVING PHASE SHIFTS, BASIS AND EIGENVECTORS, sum rules, B1 and B2
  INTEGER(KIND = I4B) :: Iphase, IsaveEN, IsaveBAS, IsaveWF, I_sumr, Igs, I_toten
  LOGICAL :: B_analytical, B_numerical
  !    
  !     ADIMENSIONAL CASE
  INTEGER(KIND = I4B) :: Iad
  !
  !     OTHER VARIABLES
  INTEGER(KIND = I4B) :: i, kx, kkx, Ierr
  REAL(KIND = DP) :: XMpar , ee, ek, eta1, eta2, lambda
  !
  CHARACTER(LEN=65) :: filename
  CHARACTER(LEN=65) :: file
  CHARACTER(LEN=56) :: prog 
  !
  !     VERBOSITY CONTROL 
  INTEGER(KIND = I4B) :: Iprint
  !
 ! REAL(KIND = DP), EXTERNAL :: POT
  REAL(KIND = DP), EXTERNAL :: X01AAF
  !
  !
  ! INTERFACES
  INTERFACE Potf
     !
     ELEMENTAL FUNCTION Potf(x)
       !
       !     WOODS-SAXON 1D POTENTIAL
       !
       USE nrtype
       USE constants
       USE pot_param
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL(KIND = DP), INTENT(IN) :: x
       REAL(KIND = DP) :: Potf
     END FUNCTION Potf
     !
  END INTERFACE Potf
  !
  INTERFACE ISQW_1D_BASIS
     SUBROUTINE ISQW_1D_BASIS(XMpar, Iprint)
       !
       USE nrtype
       USE egs_mod_isqw
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint
       REAL(KIND = DP), INTENT(IN) :: XMpar
     END SUBROUTINE ISQW_1D_BASIS
  END INTERFACE ISQW_1D_BASIS
  !
  INTERFACE BOXDIAG
     SUBROUTINE BOXDIAG(XMpar, Iprint)
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_mod_isqw
       !
       !
       IMPLICIT NONE
       !
       !ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint
       REAL(KIND = DP), INTENT(IN) :: XMpar
     END SUBROUTINE BOXDIAG
  END INTERFACE BOXDIAG
  !
  !
  INTERFACE Phase_Shift_HT
     SUBROUTINE Phase_shift_HT(dim_X, X_grid, avec_X, lambda, Index, ee, ek, eta_1, eta_2)
       !
       !     EQ(13) Hazi & Taylor PRA 1 ( 1970) 1109
       !
       !
       ! by Currix TM.
       !
       USE nrtype
       USE constants
       USE pot_param
       !
       IMPLICIT NONE
       !
       INTEGER(KIND = I4B), INTENT(IN) :: dim_X, Index
       REAL(KIND = DP), INTENT(IN) :: lambda, ee, ek
       REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
       !
       REAL(KIND = DP), INTENT(OUT) :: eta_1, eta_2
       !
     END SUBROUTINE Phase_shift_HT
     !
  END INTERFACE Phase_Shift_HT
  !
  INTERFACE B1_ISQW
     !
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
     END SUBROUTINE B1_ISQW
     !
  END INTERFACE B1_ISQW
  !
  INTERFACE B2_ISQW
     !
     SUBROUTINE B2_ISQW(Iprint, I_toten, B_numerical, B_analytical)
       !
       ! < bound | X^2 | AVEC(i) >
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
     END SUBROUTINE B2_ISQW
     !
  END INTERFACE B2_ISQW
  !
  INTERFACE Total_Strength
     !
     SUBROUTINE Total_Strength(nstates, ndim, dim_X, X_grid, avec_X, iprint)
       !
       USE nrtype
       USE constants
       !
       IMPLICIT NONE
       !
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: nstates, ndim, dim_X, Iprint
       REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
       !
       !
     END SUBROUTINE Total_Strength
     !
  END INTERFACE Total_Strength
  !
  INTERFACE Ew_Strength
     !
     SUBROUTINE Ew_Strength(ndim, dim_X, X_grid, aval, avec_X, nstates)
       !
       USE nrtype
       USE constants
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, nstates
       REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, aval
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
       !
     END SUBROUTINE Ew_Strength
     !
  END INTERFACE Ew_Strength
  !
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_X/     X_min, X_max
  NAMELIST/INP_DIM/   dim_X, dim_BOX, max_aval_BOX
  NAMELIST/INP_MASS/  Iad, reduced_mass
  NAMELIST/INP_POT/   Param_pot
  NAMELIST/INP_SHIFT/ Iphase, lambda
  NAMELIST/INP_AUX/   Igs, IsaveEN, IsaveWF, IsaveBAS, I_sumr, I_toten, B_analytical, B_numerical, Iprint 
  !
  !
  ! READING INPUT
  !
  READ(UNIT=*,NML=INP_AUX)
  IF (Iprint>3) WRITE(*,*) ' Reading IGS, ISAVEEN, ISAVEWF, ISAVEBAS, I_sumr, IPRINT'
  !
  READ(UNIT=*,NML=INP_X)
  IF (Iprint>3) WRITE(*,*) ' Reading xmin, xmax'
  !
  READ(UNIT=*,NML=INP_DIM) 
  IF (Iprint>3) WRITE(*,*) ' Reading box and xgrid basis dimension'
  dim_X = dim_X + 2 ! TO ACCOMODATE X_min AND X_max
  !
  READ(UNIT=*,NML=INP_MASS)
  IF (Iprint>3) WRITE(*,*) ' Reading mass'
  IF (Iprint>3) WRITE(*,*) ' Reading IAD'
  !
  READ(UNIT=*,NML=INP_POT)
  IF (Iprint>3) WRITE(*,*) ' Reading potential parameters'
  !
  READ(UNIT=*,NML=INP_SHIFT)
  IF (Iprint>3) WRITE(*,*) ' Reading IPHASE, lambda'
  !
  !
  !     PROGRAM VERSION
  prog = 'isqw' !to set output file's names
  !
  !     DEFINE PROBLEM UNITS
  IF (IAD == 1) THEN 
     h_sq_over_m = 1.0_DP/reduced_mass
  ELSE
     h_sq_over_m = H2OM/reduced_mass
  ENDIF
  !
  Mass_mev = UM0*reduced_mass
  !
  !
  !     DEFINE X GRID
  ALLOCATE(X_grid(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid allocation request denied."
     STOP
  ENDIF
  !
  Delta_X = (X_max-X_min)/REAL(dim_X - 1, DP)  
  X_grid = X_min + Delta_X*REAL((/ (I, I = 0, dim_X - 1) /),DP)
  !
  IF (Iprint > 1) PRINT*, ' X grid step = ', Delta_X, 'fm'
  IF (Iprint > 5) PRINT*, ' X grid = ', X_grid, 'fm'
  !    
  !     BUILDING BOX BASIS     
  !
  !     LENGTH SCALE PARAMETER \pi/(2x_m)
  XMpar = PI_D/(2.0_dp*X_max)
  !
  !
  !     BOX BASIS CALCULATION
  IF (Iprint>2) & 
       WRITE(*,*) "BOX BASIS CALCULATION XMPAR = ", XMpar, &
       " DIMENSION ", dim_BOX     
  !
  ALLOCATE(BOX_BAS(1:dim_X, 1:dim_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Har_Bas allocation request denied."
     STOP
  ENDIF
  !
  CALL ISQW_1D_BASIS(XMpar, Iprint)
  !
  !
  ALLOCATE(Aval_Box(1:dim_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_Box allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Box(1:dim_BOX, 1:dim_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Box allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Box_X(1:dim_X, 1:dim_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Box_X allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Box_Der_X(1:dim_X, 1:dim_BOX), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "AVEC_BOX_DER_X allocation request denied."
     STOP
  ENDIF
  !
  !
  !     HAMILTONIAN DIAGONALIZATION
  !
  !     PROBLEM DIAGONALIZATION IN BOX BASIS
  CALL BOXDIAG(XMpar, Iprint)
  !
  !     COMPUTING EIGENVECTOR DERIVATIVES (no)
  !
  !     SAVING GROUND STATE, XGRID, AND POTENTIAL INFORMATION
  IF (Igs>1) THEN
     !
     OPEN(UNIT=9,FILE='gs_wavefunction.dat')
     !         
     WRITE(9,*) reduced_mass, "         # Reduced Mass"
     WRITE(9,*)(Param_pot(i), i = 1,5),"        # Potential parameters"
     WRITE(9,*) X_min, X_max, Delta_X, "      # xmin xmax dx"
     WRITE(9,*) Aval_Box(1), "     # G.S. energy"
     DO kx = 1, dim_X
        WRITE(9,*) X_grid(KX), Avec_Box(kx,1)
     ENDDO
     !     
     CLOSE(9)
     !
  ENDIF
  !     SAVING BOX BASIS
  IF (IsaveBAS==1) THEN
     file = 'basis'
     !
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(70,*) "# INFSQ  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
     WRITE(70,*) "#Grid     Infinite squared well basis"
     DO kx = 1, dim_X
        WRITE(70,11) X_grid(kx), (Box_bas(kx,kkx), kkx=1, dim_BOX)
     ENDDO
     CLOSE(UNIT = 70)
  ENDIF
  !
  !     SAVING EIGENVECTORS
  IF (IsaveWF==1) THEN
     file = 'eigenvectors'
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(71,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
     WRITE(71,*) "#Grid    Eigenvectors"
     !
     DO kx = 1, dim_X
        WRITE(71,11) X_grid(kx), (Avec_Box_X(kx,kkx), kkx=1, dim_BOX_diag)
     ENDDO
     !
     CLOSE(UNIT = 71)
  ENDIF
  !
  !
  !     SAVING ENERGIES
  IF (IsaveEN==1) THEN
     !
     file = 'eigenvalues'
     !
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 100) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(73,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
     WRITE(73,*) "# Eigenvalues  Max_aval_BOX = ", Max_aval_box
     !
     DO i = 1, dim_BOX_diag
        WRITE(73,10) i, Aval_Box(i)
     ENDDO
     !
     CLOSE(UNIT = 73)
     !
     file = 'pot_eigvec'
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     WRITE(74,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
     WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
     !
     file = 'pot_eigvec2'
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 75, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     !
     WRITE(75,*) "# BOX  dim_BOX = ", dim_BOX, "dim_BOX_diag = ", dim_BOX_diag, " Box radius = ", X_max, " fm"
     WRITE(75,*) "#Grid    Potential    10*Eigenfunctions^2+eigenvalue"
     DO kx = 1, dim_X
        WRITE(74,11) X_grid(kx), Potf(X_grid(kx)), &
             10.0_dp*Avec_Box_X(kx,1:dim_BOX_diag)+Aval_Box(1:dim_BOX_diag)
     ENDDO
     DO kx = 1, dim_X
        WRITE(75,11) X_grid(kx), Potf(X_grid(kx)), &
             10.0_dp*Avec_Box_X(kx,1:dim_BOX_diag)**2+Aval_Box(1:dim_BOX_diag)
     ENDDO
     !
     CLOSE(UNIT = 74)
     CLOSE(UNIT = 75)
  ENDIF
  !
10 FORMAT (1X,I6,1X,E16.8)
11 FORMAT (1X,E14.6,1X,300E17.8E3)
  !
  !
  IF (i_SUMR /= 0) THEN
     ! COMPUTE TOTAL SUM RULE STRENGTH
     CALL TOTAL_STRENGTH(i_SUMR, dim_BOX_diag, dim_X, X_grid, Avec_Box_X, Iprint)
     ! COMPUTE ENERGY WEIGHTED SUM RULE STRENGTH
     CALL EW_STRENGTH(dim_BOX_diag, dim_X, X_grid, Aval_Box, Avec_Box_X, I_sumr)
  ENDIF
  !
  !
  ! CALCULATING B1 and B2
  IF(I_toten /= 0) THEN
     !
     CALL B1_ISQW(Iprint, I_toten, B_numerical, B_analytical)
     !
     !
     CALL B2_ISQW(Iprint, I_toten, B_numerical, B_analytical)
     !
  ENDIF
  !
  ! PHASE SHIFT CALCULATION
  IF (Iphase == 1) THEN
     IF (Iprint > 1) PRINT*, "PHASE SHIFT CALCULATION"
     file = 'phase_shift_ht'
     IF ( dim_BOX < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE IF ( dim_BOX < 100) THEN 
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ELSE
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_BOX
     ENDIF
     !
     OPEN(UNIT = 76, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(76,*) "# BOX  dim_BOX = ", dim_BOX, " Box radius = ", X_max, " fm"
     WRITE(76,*) "#2ek/Pi    Phase shift (Hazi and Taylor calculation)"
     !
     DO i = 1, dim_BOX_diag
        IF (Aval_Box(i) > 0.0_DP) THEN
           ee = Aval_Box(i)
           ek = SQRT((2.0D0*ee)/h_sq_over_m)
           !     Hazi and Taylor formula
           CALL PHASE_SHIFT_HT(dim_X, X_grid, Avec_Box_X, lambda, I, ee, ek, eta1, eta2)
           !     
           IF (Iprint > 1) WRITE(*,*) i,"-TH STATE E = ", ee, &
                "   K = ", ek, " :: eta1 = ", eta1, ", eta2 = ", eta2
           !
           WRITE(76,*) ek*2.0_DP/PI_D, eta2
           !
        ENDIF
     ENDDO
     CLOSE(UNIT = 76)
  ENDIF
  !
  ! DEALLOCATE ARRAYS
  !
  DEALLOCATE(X_grid, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(BOX_BAS, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Box_Bas allocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Aval_Box, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_Box deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Box, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Box deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Box_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Box_X deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Box_Der_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "AVEC_BOX_DER_X deallocation request denied."
     STOP
  ENDIF
  !
  !
  STOP "SAYONARA BABY..."
  !
END PROGRAM SQBOX_1BODY_1D
