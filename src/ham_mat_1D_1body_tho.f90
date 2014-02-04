SUBROUTINE THARDIAG(apt, Iprint, Iflag)
  !     
  !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A HARMONIC BASIS
  !
  !     INPUT  :: 
  !
  !     FORMAT :: IPRINT  --> VERBOSITY CONTROL
  !
  !  $Id: ham_mat_1D_1body_tho.f90,v 1.2 2013/05/22 18:03:25 laura Exp laura $
  !
  !     by Currix TM.
  !
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
  ! POTENTIAL FUNCTION
  INTERFACE Potf
     !
     ELEMENTAL FUNCTION Potf(X)
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
       REAL(KIND = DP), INTENT(IN) :: X
       REAL(KIND = DP) :: Potf
     END FUNCTION Potf
     !
  END INTERFACE Potf
  !
  !     HAMILTONIAN MATRIX
  REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE ::  Ham_Mat
  !     POTENTIAL VECTOR
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE ::  Pot_vec, Kin_vec_1, Kin_vec_2, Kin_vec_3, Kin_vec_4
  !
  !
  INTEGER(KIND = I4B) ::  I, J, KX, IFAIL, Ierr
  REAL(KIND = DP) ::  AI, AJ, AHAM, ERROR, K1, K2, K3, K4, Kinetic_En
  !
  !
  IF (Iprint > 2) PRINT*, "BUILDING HAMILTONIAN MATRIX"
  !     
  !   
  ALLOCATE(Ham_Mat(1:dim_THO,1:dim_THO), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Ham_Mat allocation request denied."
     STOP
  ENDIF
  ! 
  ALLOCATE(Pot_vec(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Pot_vec allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Kin_vec_1(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_1 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Kin_vec_2(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_2 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Kin_vec_3(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_3 allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Kin_vec_4(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_4 allocation request denied."
     STOP
  ENDIF
  ! 
  Ham_Mat = 0.0_DP
  !
  !
  !     HAMILTONIAN MATRIX
  ! 
  DO J = 1, dim_THO
     AJ = (J-1)*1.0_dp
     !    
     DO I = 1, dim_THO
        AI = (I-1)*1.0_dp
        !     KINETIC ENERGY (ADIMENSIONAL UNITS)
        Kin_vec_1 = 0.0_dp
        Kin_vec_2 = 0.0_dp
        Kin_vec_3 = 0.0_dp
        Kin_vec_4 = 0.0_dp
        !
        !same parity case
        parity_if : IF (mod(I,2) == mod(J,2)) THEN
           !
           ! Four Terms
           !
           ! Term 1
           Kin_vec_1 = ((der_S_x(:))**(-1))*((der2_S_x(:))**2)*THO_Bas(:,I)*THO_Bas(:,J)
           !
           ! Term 2
           IF ( J == 1 ) THEN
              Kin_vec_2 = der_S_x(:)*der2_S_x(:)*THO_Bas(:,I)*(-SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
           ELSE
              Kin_vec_2 = der_S_x(:)*der2_S_x(:)*THO_Bas(:,I)*(SQRT(AJ/2.0_dp)*THO_Bas(:,J-1) -&
                   SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
           ENDIF
           !
           ! Term 3
           IF ( I == 1) THEN
              Kin_vec_3 = der_S_x(:)*der2_S_x(:) * THO_Bas(:,J) * (-SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1))
           ELSE 
              Kin_vec_3 = der_S_x(:)*der2_S_x(:) * THO_Bas(:,J) * (SQRT(AI/2.0_dp)*THO_Bas(:,I-1) -&
                   SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1))
           ENDIF
           !
           ! Term 4
           IF (I==1) THEN
              IF(J==1) THEN
                 Kin_vec_4 = ((der_S_x(:))**3) * &
                      ( - SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1)) * &
                      ( - SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
              ELSE
                 Kin_vec_4 = ((der_S_x(:))**3) * &
                      ( - SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1)) * &
                      (SQRT(AJ/2.0_dp)*THO_Bas(:,J-1) - SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
              ENDIF
           ELSE
              IF(J==1) THEN
                 Kin_vec_4 = ((der_S_x(:))**3) * &
                      (SQRT(AI/2.0_dp)*THO_Bas(:,I-1) - SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1)) * &
                      ( - SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
              ELSE
                 Kin_vec_4 = ((der_S_x(:))**3) * &
                      (SQRT(AI/2.0_dp)*THO_Bas(:,I-1) - SQRT((AI+1)/2.0_dp)*THO_Bas(:,I+1)) * &
                      (SQRT(AJ/2.0_dp)*THO_Bas(:,J-1) - SQRT((AJ+1)/2.0_dp)*THO_Bas(:,J+1))
              ENDIF
           ENDIF
           !
        ENDIF parity_if
        !
        !     INTEGRATION 
        IFAIL = 0
        CALL D01GAF(X_Grid,Kin_vec_1,dim_X,K1,ERROR,IFAIL)
        IFAIL = 0
        CALL D01GAF(X_Grid,Kin_vec_2,dim_X,K2,ERROR,IFAIL)
        IFAIL = 0
        CALL D01GAF(X_Grid,Kin_vec_3,dim_X,K3,ERROR,IFAIL)
        IFAIL = 0
        CALL D01GAF(X_Grid,Kin_vec_4,dim_X,K4,ERROR,IFAIL)
        !
        Kinetic_en =((1.0_DP/8.0_DP)*h_sq_over_m)*K1 +((apt/4.0_DP)*h_sq_over_m)*K2 &
             +((apt/4.0_DP)*h_sq_over_m)*K3 + ( ((apt*apt)/2.0_DP )*h_sq_over_m)*K4
        !
        !     POTENTIAL ENERGY
        Pot_vec = Potf(X_Grid)*THO_BAS(:,I)*THO_BAS(:,J)*der_S_x(:)
        !
        !     INTEGRATION
        IFAIL = 0
        CALL D01GAF(X_Grid,Pot_vec,dim_X,AHAM,ERROR,IFAIL)
        !
        Ham_Mat(I,J) = Kinetic_en + AHAM
        !
     ENDDO
     !
  ENDDO
  !
  !
  DEALLOCATE(Pot_vec, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Pot_vec deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Kin_vec_1, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_1 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Kin_vec_2, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_2 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Kin_vec_3, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_3 deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Kin_vec_4, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Kin_vec_4 deallocation request denied."
     STOP
  ENDIF
  !     DIAGONALIZATION USING LAPACK
  IF (Iflag == 1) THEN
     CALL LA_SYEVR(A=Ham_Mat, W=AVAL_THO, JOBZ='V', UPLO='L')
  ELSE
     CALL LA_SYEVR(A=Ham_Mat, W=AVAL_THO, JOBZ='N', UPLO='L')
  ENDIF
  !
  IF (Iprint >= 1) THEN
     PRINT*, "EIGENVALUES IN A ", dim_THO, " DIM HARMONIC BASIS"
     DO I = 1, dim_THO
        WRITE(*,*) I, AVAL_THO(I)
     ENDDO
  ENDIF
  !
!print*, " AVAL_THO(1) ",  AVAL_THO(1)
  !     EIGENVECTOR MATRIX
  !
  IF (Iflag == 1) THEN
     !
     AVEC_THO = Ham_Mat
     !
     AVEC_THO_X = 0.0_DP
     !
     DO KX = 1, dim_X
        DO I = 1, dim_THO
           DO J = 1, dim_THO
              AVEC_THO_X(KX,I) = AVEC_THO_X(KX,I) + Ham_Mat(J,I)*THO_BAS(KX,J)*SQRT(der_S_x(KX))
           ENDDO
        ENDDO
     ENDDO
     !    CHECK THIS
     !        DO I = 1, dim_THO
     !           DO J = 1, dim_THO
     !              AVEC_THO_X(:,I) = AVEC_THO_X(:,I) + Ham_Mat(J,I)*HAR_BAS(:,J)
     !           ENDDO
     !        ENDDO
     !
  ENDIF
  !
  DEALLOCATE(Ham_Mat, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Ham_Mat deallocation request denied."
     STOP
  ENDIF
  !    
END SUBROUTINE THARDIAG
