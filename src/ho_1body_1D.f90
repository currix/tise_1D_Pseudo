PROGRAM HO_1BODY_1D
  !
  ! PROGRAM THAT SOLVES THE 1D TISE USING A TRUNCATED HO BASIS
  !
  ! $Id: ho_1body_1D.f90,v 1.22 2014/01/28 15:05:37 curro Exp $
  !
  ! by Currix TM
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_ho_f
  !
  IMPLICIT NONE
  !
  !
  ! FLAGS FOR SAVING BASIS AND EIGENVECTORS
  INTEGER(KIND = I4B) :: isave_EN, isave_BAS, isave_WF, i_GS
  ! FLAG TO DISPLAY SUM RULES, E1 and E2
  INTEGER(KIND = I4B) ::  i_SUMR, I_toten
  !    
  ! FLAG FOR NATURAL UNITS (H2OM = 1)
  INTEGER(KIND = I4B) ::  iad  
  !
  ! AUXILIARY VARIABLES
  INTEGER(KIND = I4B) :: Ierr, I, J, kx, Ifail, dim_HO_temp, nev
  REAL(KIND = DP) :: tol, fmult, k0, k1, k10, kmin
  REAL(KIND = DP) :: apar, egs0, ee, ek, eta1, eta2
  REAL(KIND = DP) :: eta
  !
  CHARACTER(LEN=65) :: filename
  CHARACTER(LEN=65) :: file
  CHARACTER(LEN=56) :: prog 
  !
  ! PHASE SHIFTS
  INTEGER(KIND = I4B) :: I_phase
  REAL(KIND = DP) :: lambda
  !
  ! VERBOSITY CONTROL 
  INTEGER(KIND = I4B) :: Iprint
  !
  !
  EXTERNAL EGS ! Declaring EGS as external
  !
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
  INTERFACE HO_1D_BASIS
     !
     SUBROUTINE HO_1D_BASIS(apar, ndimH, Iprint)
       !
       USE nrtype
       USE egs_ho_f
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: ndimH, Iprint
       REAL(KIND = DP), INTENT(IN) :: apar 
       !
     END SUBROUTINE HO_1D_BASIS
     !
  END INTERFACE HO_1D_BASIS
  !
  INTERFACE HARDIAG
     !
     SUBROUTINE HARDIAG(apt, Iprint, Iflag)
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_ho_f
       !
       ! Lapack 95
       USE LA_PRECISION, ONLY: WP => DP
       USE F95_LAPACK, ONLY: LA_SYEVR
       !
       IMPLICIT NONE
       !
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint, Iflag
       REAL(KIND = DP), INTENT(IN) ::  apt
       !
     END SUBROUTINE HARDIAG
     !
  END INTERFACE HARDIAG
  !
  INTERFACE WFP
     !
     SUBROUTINE WFP(indho, indx, apar, ndim)
       !
       USE nrtype
       USE egs_ho_f
       !
       IMPLICIT NONE
       !
       INTEGER(KIND = I4B), INTENT(IN) :: indho, indx, ndim
       REAL(KIND = DP), INTENT(IN) :: apar
       !
     END SUBROUTINE WFP
     !
  END INTERFACE WFP
  !
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
     END SUBROUTINE Total_Strength
     !
  END INTERFACE Total_Strength
  !
  INTERFACE Ew_Strength
     !
     SUBROUTINE Ew_Strength(ndim, dim_X, X_grid, aval, avec_X, iprint)
       !
       USE nrtype
       USE constants
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, Iprint
       REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, aval
       REAL(KIND = DP), DIMENSION(:,:), INTENT(IN) :: avec_X
       !
     END SUBROUTINE Ew_Strength
     !
  END INTERFACE Ew_Strength
  !
  !
  INTERFACE Phase_shift_Ch
     !
     SUBROUTINE Phase_shift_Ch(Index,ek,eta)
       !     
       !     EQ(20) Chadan et al JMP 42 (2001) 4031
       !     
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_ho_f
       !
       IMPLICIT NONE
       !
       !      
       INTEGER(KIND = I4B), INTENT(IN) :: Index
       REAL(KIND = DP), INTENT(IN)  :: ek
       !
       REAL(KIND = DP), INTENT(OUT) :: eta
       !
     END SUBROUTINE Phase_shift_Ch
     !
  END INTERFACE Phase_shift_Ch
  !
  !
  INTERFACE Phase_shift_HT 
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
  END INTERFACE Phase_shift_HT
  !
  !
  INTERFACE E1_HO
     !
     SUBROUTINE E1_HO(Iprint, I_toten, apar)
       !
       USE nrtype
       USE constants
       USE egs_ho_f
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
       REAL(KIND = DP), INTENT(IN) :: apar
       !
     END SUBROUTINE E1_HO
     !
  END INTERFACE E1_HO
  !
  INTERFACE E2_HO
     !
     SUBROUTINE E2_HO(Iprint, I_toten, apar)
       !
       USE nrtype
       USE constants
       USE egs_ho_f
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint, I_toten
       REAL(KIND = DP), INTENT(IN) :: apar
       !
     END SUBROUTINE E2_HO
     !
  END INTERFACE E2_HO
  !
  ! DATA INPUT
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_X/     X_min, X_max
  NAMELIST/INP_DIM/   dim_X, dim_HO
  NAMELIST/INP_MASS/  iad, reduced_mass
  NAMELIST/INP_POT/   Param_pot
  NAMELIST/INP_SHIFT/ I_phase, lambda
  NAMELIST/INP_AUX/   i_GS, i_SUMR, I_toten, isave_EN, isave_WF, isave_BAS, Iprint
  !
  ! PROGRAM VERSION
  IF (Iprint > 1) PRINT*, "$Id: ho_1body_1D.f90,v 1.22 2014/01/28 15:05:37 curro Exp $"
  !
  prog = 'ho' !to set output file's names
  !
  ! NAMELIST FILE
  ! OPEN(UNIT=10,FILE='sec_order.inp',STATUS='OLD')
  !
  ! READING INPUT
  !
  READ(UNIT=*,NML=INP_X)
  !
  READ(UNIT=*,NML=INP_DIM)
  dim_X = dim_X + 2 ! TO ACCOMODATE X_min AND X_max
  !
  READ(UNIT=*,NML=INP_MASS)
  !
  READ(UNIT=*,NML=INP_POT)
  !
  READ(UNIT=*,NML=INP_SHIFT)
  !
  READ(UNIT=*,NML=INP_AUX)
  !  
  !
  ! DEFINE PROBLEM UNITS
  IF (iad == 1) THEN 
     h_sq_over_m = 1.0_DP/reduced_mass
  ELSE
     h_sq_over_m = H2OM/reduced_mass
  ENDIF
  !
  MASS_MEV = UM0*reduced_mass
  !
  !
  ! DEFINE X GRID
  ALLOCATE(X_grid(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid allocation request denied."
     STOP
  ENDIF
  !
  Delta_X = (X_max-X_min)/REAL(dim_X - 1, DP)  
  X_grid = X_min + Delta_X*REAL((/ (I, I = 0, dim_X-1) /),DP)
  !
  IF (Iprint > 1) PRINT*, ' X grid step = ', Delta_X, 'fm'
  IF (Iprint > 5) PRINT*, ' X grid = ', X_grid, 'fm'
  !
  IF (IAD == 1) THEN
     !
     kmin = 1.0_DP
     !
!  ELSE IF (IAD == 3) THEN 
!     PRINT*, "ENTER KMIN"
!     READ*, kmin
  ELSE
     !
     ! Compute optimum apar value
     dim_HO_temp = dim_HO
     dim_HO = 1
     !
     fmult = 5.0_DP/1.5_DP
     !
     !
     ALLOCATE(Har_Bas(1:dim_X, 1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Har_Bas allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Aval_Har(1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Aval_Har allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Avec_Har(1:1, 1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_Har allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Avec_Har_X(1:dim_X, 1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_Har_X allocation request denied."
        STOP
     ENDIF
     !
     !
     DO
        fmult = fmult*1.5_DP
        ! LENGTH SCALE PARAMETER INITIAL EVALUATION
        k0 = 0.0_DP
        !     PARAM(1) -> Well depth  &  PARAM(2)  -> Well length unit 
        IF (Param_pot(2) /= 0.0_DP) THEN
           k1 = ABS(fmult*Param_pot(1)/(Param_pot(2)*Param_pot(2)))
        ELSE
           k1 = ABS(fmult*Param_pot(1))
        ENDIF
        k10 = k1
        !     
        !
        !
        Ifail = 0
        tol = 0.0_DP
        nev = 50
        !
        CALL E04ABF(EGS, tol, tol, k0, k1, nev, kmin, egs0, Ifail)
        !
        IF (Iprint > 2) PRINT*, "kmin = ", kmin, " kmin - k10 = ", kmin-k10, & 
             " a = ", SQRT(sqrt(kmin/h_sq_over_m)), &
             " EGS = ", egs0
        !     
        IF (ABS(kmin-k10) > 1.0D-6) EXIT
        !
     ENDDO
     !
     dim_HO = dim_HO_temp
     !
     DEALLOCATE(Har_Bas, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Har_Bas deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Aval_Har, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Aval_Har deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Avec_Har, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_Har deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Avec_Har_X, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_Har_X deallocation request denied."
        STOP
     ENDIF
     !
     IF (Iprint > 2) PRINT*, "MINIMUM:: kmin = ", kmin, & 
          " a = ", SQRT(sqrt(kmin/h_sq_over_m)), &
          " EGS = ", egs0
     !
  ENDIF
  !    INVERSE OSCILLATOR LENGTH a = (\nu K/\hbar^2)^(1/4)  (fm^{-1})
  apar = SQRT(SQRT(kmin/h_sq_over_m)) 
  !
  IF (Iprint > 2) PRINT*,  "HARMONIC BASIS CALCULATION apar = ", apar, " DIMENSION ", dim_HO     
  !
  !  Add one for the calculation of derivatives in the wfp subroutine
  ALLOCATE(Har_Bas(1:dim_X, 1:dim_HO + 1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Har_Bas allocation request denied."
     STOP
  ENDIF
  !
  CALL HO_1D_BASIS(apar, dim_HO + 1, Iprint)
  !
  ALLOCATE(Aval_Har(1:dim_HO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_Har allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Har(1:dim_HO, 1:dim_HO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Har_X(1:dim_X, 1:dim_HO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har_X allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_Har_Der_X(1:dim_X, 1:dim_HO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har_Der_X allocation request denied."
     STOP
  ENDIF
  !
  !     HAMILTONIAN DIAGONALIZATION
  CALL HARDIAG(apar, Iprint, 1)
  !
  !
  !     COMPUTING EIGENVECTOR DERIVATIVES
  !     EIGENSTATE DERIVATIVE
  DO J = 1, dim_X
     DO I = 1, dim_HO
        CALL WFP(I, J, apar, dim_HO)
     ENDDO
  ENDDO
  !
  ! SAVING GROUND STATE, XGRID, AND POTENTIAL INFORMATION
  IF (i_GS == 1) THEN
     !
     OPEN(UNIT=9,FILE='gs_wavefunction.dat')
     !         
     WRITE(9,*) "# ", reduced_mass, "         # Reduced Mass"
     WRITE(9,*) "# ", Param_pot, "        # Potential parameters"
     WRITE(9,*) "# ", X_min, X_max, Delta_X, "      # xmin xmax dx"
     WRITE(9,*) "# ", Aval_Har(1), "     # G.S. energy"
     !
     DO kx = 1, dim_X
        WRITE(9,*) X_grid(kx), Avec_Har(kx,1)
     ENDDO
     !  
     CLOSE(9)
     !
  ENDIF
  ! SAVING HARMONIC BASIS
  IF (isave_BAS == 1) THEN
     file = 'basis'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     IF ( dim_HO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     !
     OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(70,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(70,*) "#Grid    Harmonic basis"
     DO kx = 1, dim_X
        WRITE(70,11) X_grid(kx), Har_Bas(kx,:)
     ENDDO
     CLOSE(UNIT = 70)
  ENDIF
  !
  ! SAVING EIGENVECTORS
  IF (isave_WF == 1) THEN
     file = 'eigenvectors'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     IF ( dim_HO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     !
     OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(71,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(71,*) "#Grid     Eigenvectors"
     DO kx = 1, dim_X
        WRITE(71,11) X_grid(kx), Avec_Har_X(kx,1:dim_HO)
     ENDDO
     CLOSE(UNIT = 71)
  ENDIF
  !
  ! SAVING EIGENVECTOR DERIVATIVES
  IF (isave_WF == 1) THEN
     file = 'eigvec_der'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE IF ( dim_HO < 100) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     !
     OPEN(UNIT = 72, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(72,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(72,*) "#Grid    Eigenvector derivatives"
     DO kx = 1, dim_X
        WRITE(72,11) X_grid(kx), Avec_Har_Der_X(kx,1:dim_HO)
     ENDDO
     CLOSE(UNIT = 72)
  ENDIF
  !
  !
  ! SAVING ENERGIES
  IF (isave_EN == 1) THEN
     file = 'eigenvalues'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE IF ( dim_HO < 100) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     !
     OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(73,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(73,*) "# Eigenvalues"
     DO I = 1, dim_HO
        WRITE(73,10) I, Aval_Har(I)
     ENDDO
     CLOSE(UNIT = 73)
     !
     file = 'pot_eigvec'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE IF ( dim_HO < 100) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(74,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
     !
     file = 'pot_eigvec2'
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE IF ( dim_HO < 100) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ELSE ! Max dim 999 (dim > 999 -> asterisks will appear)
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     OPEN(UNIT = 75, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(75,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(75,*) "#Grid    Potential    10*Eigenfunctions^2+eigenvalue"
     DO kx = 1, dim_X
        WRITE(74,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_Har_X(kx,1:dim_HO) + Aval_Har(1:dim_HO)
        WRITE(75,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_Har_X(kx,1:dim_HO)**2 + Aval_Har(1:dim_HO)
     ENDDO
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
     CALL Total_Strength(i_SUMR, dim_HO, dim_X, X_grid, Avec_Har_X, Iprint)
     ! COMPUTE ENERGY WEIGHTED SUM RULE STRENGTH
     CALL EW_STRENGTH(dim_HO, dim_X, X_grid, Aval_Har, Avec_Har_X, iprint)
  ENDIF
  !
  IF(I_toten /= 0) THEN
     !CALCULATING E1
     CALL E1_HO(Iprint, I_toten, apar)
     !
     !CALCULATING E2
     CALL E2_HO(Iprint, I_toten, apar)
  ENDIF
  !
  ! PHASE SHIFT CALCULATION
  IF (I_phase == 1) THEN
     IF (Iprint > 1) PRINT*, "PHASE SHIFT CALCULATION"
     file = 'phase_shift_ht'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     IF ( dim_HO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     OPEN(UNIT = 76, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(76,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(76,*) "#2ek/Pi    Phase shift (Hazi and Taylor calculation)"
     !
     file = 'phase_shift_ch'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_HO
     IF ( dim_HO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     IF ( dim_HO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_HO
     ENDIF
     OPEN(UNIT = 77, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(77,*) "# HO  dim_HO = ", dim_HO, " Box radius = ", X_max, " fm"
     WRITE(77,*) "#2ek/Pi      Phase shift (Chadan calculation)"
     !
     DO I = 1, dim_HO
        IF (Aval_Har(I) > 0.0_DP) THEN
           ee = Aval_Har(I)
           ek = SQRT((2.0D0*ee)/h_sq_over_m)
           !     Hazi and Taylor formula
           CALL PHASE_SHIFT_HT(dim_X, X_grid, Avec_Har_X, lambda, I, ee, ek, eta1, eta2)
           !     
           IF (Iprint > 1) WRITE(*,*) I,"-TH STATE E = ", ee, &
                "   K = ", ek, " :: eta1 = ", eta1, ", eta2 = ", eta2
           !
           WRITE(76,*) ek*2.0_DP/PI_D, eta2
           !
           ! Chadan et al. formula
           CALL PHASE_SHIFT_CH(I,ek,eta)
           !
           IF (Iprint > 1) WRITE(*,*) I,"-TH STATE E = ", ee, &
                "   K = ", ek, " :: eta = ", eta
           !
           WRITE(77,*) ek*2.0_DP/PI_D, eta
           !     
        ENDIF
     ENDDO
     CLOSE(UNIT = 76)
     CLOSE(UNIT = 77)
  ENDIF
  !
  DEALLOCATE(Har_Bas, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Har_Bas deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Aval_Har, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_Har deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Har, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Har_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har_X deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_Har_Der_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_Har_Der_X deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(X_grid, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid deallocation request denied."
     STOP
  ENDIF
  !
  !
  STOP "SAYONARA BABY..."
  !
END PROGRAM HO_1BODY_1D
    
