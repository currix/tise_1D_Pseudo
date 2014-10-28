PROGRAM THO_1BODY_1D
  !
  ! PROGRAM THAT SOLVES THE 1D TISE USING A TRUNCATED HO BASIS
  !
  ! $Id: tho_1body_1D.f90,v 1.10 2013/05/22 18:04:42 laura Exp laura $
  !
  ! by Currix TM
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_tho_f
  USE lst_param
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=65) :: filename
  CHARACTER(LEN=65) :: file
  CHARACTER(LEN=56) :: prog 
  !
  ! FLAGS FOR SAVING BASIS AND EIGENVECTORS
  INTEGER(KIND = I4B) :: Isave_EN, Isave_BAS, Isave_WF, I_GS
  ! FLAG TO DISPLAY SUM RULES, E1 and E2
  INTEGER(KIND = I4B) ::  I_sumr, I_toten
  !    
  ! FLAG FOR NATURAL UNITS (H2OM = 1)
  INTEGER(KIND = I4B) ::  iad  
  !
  ! AUXILIARY VARIABLES
  INTEGER(KIND = I4B) :: Ierr, I, J, kx, Ifail, nev, dim_THO_temp
  REAL(KIND = DP) ::   ee,    ek,  eta1, eta2
  REAL(KIND = DP) :: amin,  apar,  egs0
  REAL(KIND = DP) ::  tol, fmult,  a0, a1, a10
  !  REAL(KIND = DP) :: eta
  !
  ! PHASE SHIFTS
  INTEGER(KIND = I4B) :: I_phase
  REAL(KIND = DP) :: lambda
  !
  ! VERBOSITY CONTROL 
  INTEGER(KIND = I4B) :: Iprint
  !
  !
  EXTERNAL EGS ! Declaring EGS as external (E04ABF's request)
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
  INTERFACE THO_1D_BASIS
     !
     SUBROUTINE THO_1D_BASIS(apar, ndimH, Iprint)
       !
       USE nrtype
       USE egs_tho_f
       USE constants
       USE lst_param
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: ndimH, Iprint
       REAL(KIND = DP), INTENT(IN) :: apar 
       !
     END SUBROUTINE THO_1D_BASIS
     !
  END INTERFACE THO_1D_BASIS
  !
  INTERFACE THARDIAG
     !
     SUBROUTINE THARDIAG(apt, Iprint, Iflag)
       !
       USE nrtype
       USE constants
       USE pot_param
       USE egs_tho_f
       USE lst_param
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
     END SUBROUTINE THARDIAG
     !
  END INTERFACE THARDIAG
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
  INTERFACE Total_Strength
     !
     SUBROUTINE Total_Strength(ndim, dim_X, X_grid, avec_X, Aval, iprint)
       !
       USE nrtype
       USE constants
       !
       IMPLICIT NONE
       !
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: ndim, dim_X, Iprint
       REAL(KIND = DP), DIMENSION(:), INTENT(IN) :: X_grid, Aval
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
  INTERFACE E1
     !
     SUBROUTINE E1(Iprint)
       !
       USE nrtype
       USE constants
       USE egs_tho_f
       USE lst_param
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint
       !
     END SUBROUTINE E1
     !
  END INTERFACE E1
  !
  INTERFACE E2
     !
     SUBROUTINE E2(Iprint)
       !
       USE nrtype
       USE constants
       USE egs_tho_f
       USE lst_param
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       INTEGER(KIND = I4B), INTENT(IN) :: Iprint
       !
     END SUBROUTINE E2
     !
  END INTERFACE E2
  !
  INTERFACE WFP
     SUBROUTINE WFP(indho, indx, apar, ndim)
       !
       USE nrtype
       USE egs_tho_f
       USE lst_param
       !
       IMPLICIT NONE
       !
       INTEGER(KIND = I4B), INTENT(IN) :: indho, indx, ndim
       REAL(KIND = DP), INTENT(IN) :: apar
       !
     END SUBROUTINE WFP
  END INTERFACE WFP
  !
  !
  !
  ! DATA INPUT
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_X/     X_min, X_max, ratio
  NAMELIST/INP_DIM/   dim_X, dim_THO
  NAMELIST/INP_MASS/  iad, reduced_mass
  NAMELIST/INP_POT/   Param_pot
  NAMELIST/INP_SHIFT/ I_phase, lambda
  NAMELIST/INP_AUX/   i_GS, i_SUMR, I_toten, isave_EN, isave_WF, isave_BAS, Iprint
  !
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
  ! PROGRAM VERSION
  IF (Iprint > 1) PRINT*, "$Id: tho_1body_1D.f90,v 1.10 2013/05/22 18:04:42 laura Exp laura $"
  !
  prog = 'tho' !to set output file's names
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
  ! K_eff = SQRT(2.0_dp*reduced_mass*70.0_dp)/HBAR !70MeV is the maximum excitation energy populated in the reaction
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
  !     OSCILLATOR LENGTH b = (\nu K/\hbar^2)^(-1/4)  (fm)
  IF (IAD == 1) THEN
     !
     b = (SQRT(SQRT(h_sq_over_m)))
     apar = 1.0_DP/b
!  ELSE IF (IAD == 3) THEN
!     PRINT*, "ENTER b"
!     READ*, b
  ELSE 
     ! Compute optimum apar value
     dim_THO_temp = dim_THO
     dim_THO = 1
     !
     !fmult = 5.0_DP/1.5_DP
     fmult = 1.0_DP/1.5_DP
     !
     !
     ! Basis dimension = dim + 1 = 2 for Hamiltonian diagonalization
     ALLOCATE(THO_Bas(1:dim_X, 1:2), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "THO_Bas allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(S_x(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "S_x allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(der_S_x(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "der_S_x allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(der2_S_x(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "der2_S_x allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(S_over_x(1:dim_X), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "S_over_x allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Aval_THO(1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Aval_THO allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Avec_THO(1:1, 1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_THO allocation request denied."
        STOP
     ENDIF
     !
     ALLOCATE(Avec_THO_X(1:dim_X, 1:1), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_THO_X allocation request denied."
        STOP
     ENDIF
     !
     !
     !
     DO
        fmult = fmult*1.5_DP
        ! LENGTH SCALE PARAMETER INITIAL EVALUATION
        a0 = 0.0_DP
        ! DEFINE THO CONSTANTS
        !
        IF (Param_pot(2) /= 0.0_DP) THEN
           a1 = 1.0_DP/SQRT(SQRT(ABS(fmult*Param_pot(1)/(Param_pot(2)*Param_pot(2)))/h_sq_over_m))
        ELSE
           a1 = 1.0_DP/SQRT(SQRT(ABS(fmult*Param_pot(1))/h_sq_over_m))
        ENDIF
        a1 = fmult*0.5_DP
        a10 = a1
        !     
        !
        !
        Ifail = 0
        tol = 0.0_DP ! se la fissi a zero come defoult è la sqrt(epsilon)
        nev = 50
        !
        IF (Iprint > 2) print*, "a1 = ", a1
        CALL EGS(a1,egs0)
        IF (Iprint > 2) print*, "energy = ", egs0
        a1 = 2.0_DP*a1
        IF (Iprint > 2) print*, "a1 = ", a1
        CALL EGS(a1,egs0)
        IF (Iprint > 2) print*, "energy = ", egs0
        CALL E04ABF(EGS, tol, tol, a0, a1, nev, amin, egs0, Ifail)
        !
        IF (Iprint > 2) PRINT*, "amin = ", amin, " amin - a10 = ", amin-a10, & 
             " bmin = ", 1/amin," fm; EGS = ", egs0
        !     
        IF (ABS(amin-a10) > 1.0D-6) EXIT
        EXIT
        !
     ENDDO
     !
     dim_THO = dim_THO_temp
     !
     DEALLOCATE(THO_Bas, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "THO_Bas deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(S_x, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "S_x deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(der_S_x, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "der_S_x deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(der2_S_x, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "der2_S_x deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(S_over_x, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "S_over_x deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Aval_THO, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Aval_THO deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Avec_THO, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_THO deallocation request denied."
        STOP
     ENDIF
     !
     DEALLOCATE(Avec_THO_X, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "Avec_THO_X deallocation request denied."
        STOP
     ENDIF
     !
     apar = amin
     IF (Iprint > 2) PRINT*, "MINIMUM:: bmin = ", 1/amin, & 
          " a = ", apar, &
          " EGS = ", egs0
     !
  ENDIF
  !
  !
  !
  IF (Iprint > 1) PRINT*,  "TRANSFORMED HARMONIC BASIS DIMENSION ", dim_THO     
  !
  !  Add one for the calculation of derivatives in the wfp subroutine and the kinetic energy calculation
  ALLOCATE(THO_Bas(1:dim_X, 1:dim_THO + 1), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "THO_Bas allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(der_S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der_S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(der2_S_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der2_S_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(S_over_x(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_over_x allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Aval_THO(1:dim_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_THO allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_THO(1:dim_THO, 1:dim_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_THO_X(1:dim_X, 1:dim_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO_X allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Avec_THO_Der_X(1:dim_X, 1:dim_THO), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO_Der_X allocation request denied."
     STOP
  ENDIF
  !
  !     BASIS CONSTRUCTION
  CALL THO_1D_BASIS(apar, dim_THO + 1, Iprint)
  !
  !     HAMILTONIAN DIAGONALIZATION
  CALL THARDIAG(apar, Iprint, 1)
  !
  !
  !     COMPUTING EIGENVECTOR DERIVATIVES
  !     EIGENSTATE DERIVATIVE
  DO J = 1, dim_X
     DO I = 1, dim_THO
        CALL WFP(I, J, apar, dim_THO)
     ENDDO
  ENDDO
  !
  !
  !
  IF (i_SUMR == 1) THEN
     ! COMPUTE TOTAL SUM RULE STRENGTH
     CALL Total_Strength(dim_THO, dim_X, X_grid, Avec_THO_X, Aval_THO, Iprint)
     ! COMPUTE ENERGY WEIGHTED SUM RULE STRENGTH
     CALL EW_STRENGTH(dim_THO, dim_X, X_grid, Aval_THO, Avec_THO_X, Iprint)
  ENDIF
  !
  !
  IF( I_toten == 1) THEN
     !CALCULATING E1
     CALL E1(Iprint)
     !
     !CALCULATING E2
     CALL E2(Iprint)
  ENDIF
  !
  !
  ! PHASE SHIFT
  IF (I_phase == 1) THEN
     IF (Iprint > 1) PRINT*, "PHASE SHIFT CALCULATION"
     file = 'phase_shift_ht'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 76, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(76,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(76,*) "# Phase shift (Hazi and Taylor calculation)"
     DO I = 1, dim_THO
        IF (Aval_THO(I) > 0.0_DP) THEN
           ee = Aval_THO(I)
           ek = SQRT((2.0D0*ee)/h_sq_over_m)
           !     Hazi and Taylor formula
           CALL PHASE_SHIFT_HT(dim_X, X_grid, Avec_THO_X, lambda, I, ee, ek, eta1, eta2)
           !     
           IF (Iprint > 1) WRITE(*,*) I,"-TH STATE E = ", ee, &
                "   K = ", ek, " :: eta1 = ", eta1, ", eta2 = ", eta2
           !
           WRITE(76,*) ek*2.0_DP/PI_D, eta2
           !
           !           ! Chadan et al. formula
           !           CALL PHASE_SHIFT_CH(I,ek,eta)
           !
           !           IF (Iprint > 0) WRITE(*,*) I,"-TH STATE E = ", ee, &
           !                "   K = ", ek, " :: eta = ", eta
           !
           !           WRITE(97,*) ek*2.0_DP/PI_D, eta
           !     
        ENDIF
     ENDDO
     CLOSE(UNIT = 76)
  ENDIF
  !
  ! SAVING GROUND STATE, XGRID, AND POTENTIAL INFORMATION
  IF (i_GS == 1) THEN
     !
     OPEN(UNIT=9,FILE='gs_wavefunction.dat')
     !         
     WRITE(9,*) "# ", reduced_mass, "         # Reduced Mass"
     WRITE(9,*) "# ", Param_pot, "        # Potential parameters"
     WRITE(9,*) "# ", X_min, X_max, Delta_X, "      # xmin xmax dx"
     WRITE(9,*) "# ", Aval_THO(1), "     # G.S. energy"
     !
     DO kx = 1, dim_X
        WRITE(9,*) X_grid(kx), Avec_THO(kx,1)
     ENDDO
     !  
     CLOSE(9)
     !
  ENDIF
  !
  ! SAVING T-HARMONIC BASIS
  IF (isave_BAS == 1) THEN
     file = 'basis'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 70, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(70,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(70,*) "#Grid    T-Harmonic basis"
     DO kx = 1, dim_X
        WRITE(70,11) X_grid(kx), THO_Bas(kx,1:dim_THO)*SQRT(der_S_x(kx))
     ENDDO
     CLOSE(UNIT = 70)
  ENDIF
  !
  ! SAVING EIGENVECTORS
  IF (isave_WF == 1) THEN
     file = 'eigenvectors'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 71, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(71,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(71,*) "#Grid  Eigenvectors"
     DO kx = 1, dim_X
        WRITE(71,11) X_grid(kx), Avec_THO_X(kx,1:dim_THO)
     ENDDO
     CLOSE(UNIT = 71)
  ENDIF
  !
  ! SAVING EIGENVECTOR DERIVATIVES
  IF (isave_WF == 1) THEN
     file = 'eigvec_der'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 72, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(72,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(72,*) "#Grid  Eigenvector derivatives"
     DO kx = 1, dim_X
        WRITE(72,11) X_grid(kx), Avec_THO_Der_X(kx,1:dim_THO)
     ENDDO
     CLOSE(UNIT = 72)
  ENDIF
  !
  !
  ! SAVING ENERGIES
  IF (isave_EN == 1) THEN
     file = 'eigenvalues'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 73, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(73,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(73,*) "# Eigenvalues"
     DO I = 1, dim_THO
        WRITE(73,10) I, Aval_THO(I)
     ENDDO
     CLOSE(UNIT = 73)
     !
     file = 'pot_eigvec'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 74, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(74,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(74,*) "#Grid    Potential    10*Eigenfunctions+eigenvalue"
     !
     file = 'pot_eigvec2'
     WRITE(filename, '(A, "_",A,"_N",I2, ".dat")') TRIM(prog), TRIM(file), dim_THO
     IF ( dim_THO < 10) THEN !to avoid spaces
        WRITE(filename, '(A, "_",A,"_N",I1, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     IF ( dim_THO > 99) THEN 
        WRITE(filename, '(A, "_",A,"_N",I3, ".dat")') TRIM(prog), TRIM(file), dim_THO
     ENDIF
     OPEN(UNIT = 75, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(75,*) "# THO  dim_THO = ", dim_THO, " Box radius = ", X_max, " fm"
     WRITE(75,*) "#Grid    Potential    10*Eigenfunctions^2+eigenvalue"
     DO kx = 1, dim_X
        WRITE(74,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_THO_X(kx,1:dim_THO) + Aval_THO(1:dim_THO)
        WRITE(75,11) X_grid(kx), Potf(X_grid(kx)), 10.0_DP*Avec_THO_X(kx,1:dim_THO)**2 + Aval_THO(1:dim_THO)
     ENDDO
     CLOSE(UNIT = 74)
     CLOSE(UNIT = 75)
  ENDIF
  !
  !
10 FORMAT (1X,I6,1X,E16.8)
11 FORMAT (1X,E14.6,1X,300E16.8)
  !
  !
  DEALLOCATE(THO_Bas, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "THO_Bas deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Aval_THO, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Aval_THO deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_THO, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_THO_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO_X deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Avec_THO_Der_X, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Avec_THO_Der_X deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(X_grid, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(der_S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der_S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(der2_S_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "der2_S_x deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(S_over_x, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "S_over_x deallocation request denied."
     STOP
  ENDIF
  !
  STOP "SAYONARA BABY..."
  !
END PROGRAM THO_1BODY_1D
    
