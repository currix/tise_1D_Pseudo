SUBROUTINE HARDIAG(apt, Iprint, Iflag)
  !     
  !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A HARMONIC BASIS
  !
  !     INPUT  :: 
  !               apt     --> PROBLEM LENGTH SCALE (fm^-1)
  !               NDIM    --> HARMONIC BASIS DIMENSION
  !               IFLAG   --> INTEGER :: 0 EIGENVALUES 
  !                                      1 EIGENVALUES AND EIGENVECTORS
  !
  !     FORMAT :: IPRINT  --> VERBOSITY CONTROL
  !
  !     by Currix TM.
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
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE ::  Pot_vec
  !
  !
  INTEGER(KIND = I4B) ::  I, J, KX, IFAIL, Ierr
  REAL(KIND = DP) ::  apt2, AJ, AHAM, ERROR
  !
  apt2 = apt*apt*H_SQ_OVER_M/2.0_DP
  !
  IF (Iprint > 2) PRINT*, "BUILDING HAMILTONIAN MATRIX"
  !     
  !   
  ALLOCATE(Ham_Mat(1:dim_HO,1:dim_HO), STAT = Ierr)
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
  Ham_Mat = 0.0_DP
  !
  !
  !     HAMILTONIAN MATRIX
  !     
  DO I = 1, dim_HO
     !
     DO J = I, dim_HO
        AJ = REAL(J - 1, DP)
        !     KINETIC ENERGY (ADIMENSIONAL UNITS)
        IF (I.EQ.J-2) Ham_Mat(J,I) = Ham_Mat(J,I) -  apt2*SQRT(AJ*(AJ-1.0_DP))/2.0_DP
        IF (I.EQ.J) Ham_Mat(J,I)   = Ham_Mat(J,I) +  apt2*(2.0_DP*AJ+1.0_DP)/2.0_DP
        !
        !     POTENTIAL ENERGY
        Pot_vec = Potf(X_Grid)*HAR_BAS(:,I)*HAR_BAS(:,J)
        !
        !     INTEGRATION
        IFAIL = 0
        CALL D01GAF(X_Grid,Pot_vec,dim_X,AHAM,ERROR,IFAIL)
        Ham_Mat(J,I) = Ham_Mat(J,I) + AHAM
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(Pot_vec, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Pot_vec deallocation request denied."
     STOP
  ENDIF
  !     DIAGONALIZATION USING LAPACK 95
  IF (max_aval_har == 0.0) THEN
     IF (Iflag == 1) THEN
        CALL LA_SYEVR(A=Ham_Mat, W=AVAL_HAR, JOBZ='V', UPLO='L')
     ELSE
        CALL LA_SYEVR(A=Ham_Mat, W=AVAL_HAR, JOBZ='N', UPLO='L')
     ENDIF
     !
     dim_HO_diag = dim_HO
  ELSE
     IF (Iflag == 1) THEN
        CALL LA_SYEVR(A=Ham_Mat, W=AVAL_HAR, JOBZ='V', VU = max_aval_har, M = dim_HO_diag, UPLO='L')
     ELSE
        CALL LA_SYEVR(A=Ham_Mat, W=AVAL_HAR, JOBZ='N', VU = max_aval_har, M = dim_HO_diag, UPLO='L')
     ENDIF
  ENDIF
  !
  IF (Iprint >= 1) THEN
     PRINT*, dim_HO_diag," EIGENVALUES COMPUTED IN A ", dim_HO, " DIM HARMONIC BASIS"
     DO I = 1, dim_HO_diag
        WRITE(*,*) I, AVAL_HAR(I)
     ENDDO
  ENDIF
  !
  !     EIGENVECTOR MATRIX
  !
  IF (Iflag == 1) THEN
     !
     AVEC_HAR = Ham_Mat
     !
     AVEC_HAR_X = 0.0_DP
     !
     DO KX = 1, dim_X
        DO I = 1, dim_HO_diag
           DO J = 1, dim_HO
              AVEC_HAR_X(KX,I) = AVEC_HAR_X(KX,I) + Ham_Mat(J,I)*HAR_BAS(KX,J)
           ENDDO
        ENDDO
     ENDDO
     !    CHECK THIS
     !        DO I = 1, dim_HO
     !           DO J = 1, dim_HO
     !              AVEC_HAR_X(:,I) = AVEC_HAR_X(:,I) + Ham_Mat(J,I)*HAR_BAS(:,J)
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
END SUBROUTINE HARDIAG
