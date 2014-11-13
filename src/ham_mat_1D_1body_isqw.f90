SUBROUTINE BOXDIAG(XMpar, Iprint)
  !
  !     COMPUTES AND DIAGONALIZES 1D HAMILTONIAN USING A BOXBASIS
  !
  !     INPUT      :: XMPAR    --> PI/(2 X_MAX)
  !
  !     FORMAT     :: IPRINT  --> VERBOSITY CONTROL
  !
  !     by Currix TM.
  !
  USE nrtype
  USE constants
  USE pot_param
  USE egs_mod_isqw
  !
  !
  ! Lapack 95
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR
  !
  IMPLICIT NONE
  !
  !
  !
  !ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: Iprint
  REAL(KIND = DP), INTENT(IN) :: XMpar
  !
  INTEGER(KIND = I4B) :: Ierr
  REAL(KIND = DP) :: XMpar2
  !
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
  REAL(KIND = DP) :: ai, aham, error
  INTEGER(KIND = I4B) :: i, j, kx, Ifail
  !
  XMpar2 = XMpar*XMpar*h_sq_over_m/2.0_dp
  !
  IF (IPRINT>2) WRITE(*,*) "BUILDING BOX HAMILTONIAN MATRIX"
  !
  !   
  ALLOCATE(Ham_Mat(1:dim_BOX,1:dim_BOX), STAT = Ierr)
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
  !    HAMILTONIAN MATRIX
  !
  DO i = 1, dim_BOX
     ai = REAL(i, DP)
     !     KINETIC ENERGY 
     Ham_Mat(i,i) = Ham_Mat(i,i) + XMpar2*ai*ai
     !     
     DO j = i, dim_BOX
        !     POTENTIAL ENERGY
        Pot_vec = Potf(X_Grid)*Box_Bas(:,i)*Box_Bas(:,j)
        !
        !     INTEGRATION
        Ifail = 0
        CALL D01GAF(X_grid, Pot_vec,dim_X,aham,error,Ifail)
        Ham_Mat(j,i) = Ham_Mat(j,i) + aham
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
  !
  !
  !     DIAGONALIZATION USING LAPACK 95
  IF (max_aval_box == 0.0) THEN
     CALL LA_SYEVR(A=Ham_Mat, W=Aval_Box, JOBZ='V', UPLO='L')
     !
     dim_BOX_diag = dim_BOX
  ELSE
     CALL LA_SYEVR(A=Ham_Mat, W=AVAL_Box, JOBZ='V', VU = max_aval_box, M = dim_BOX_diag, UPLO='L')
  ENDIF
  !     
  IF (Iprint>=1) THEN
     WRITE(*,*) "EIGENVALUES IN A ", dim_BOX, " DIM BOX BASIS"
     DO i = 1, dim_BOX_diag
        WRITE(*,*) i, Aval_Box(i)
     ENDDO
  ENDIF
  !
  IF (Iprint >= 1) THEN
     PRINT*, dim_BOX_diag," EIGENVALUES COMPUTED IN A ", dim_BOX, " DIM ISQW BASIS"
     DO I = 1, dim_BOX_diag
        WRITE(*,*) I, AVAL_BOX(I)
     ENDDO
  ENDIF
  !
  !     EIGENVECTOR MATRIX
  !
  Avec_Box = Ham_Mat
  Avec_Box_x = 0.0_dp
  !
  DO kx = 1, dim_X
     DO i = 1, dim_BOX
        DO j = 1, dim_BOX
           Avec_Box_x(kx,i) = Avec_Box_x(kx,i) + Ham_Mat(j,i)*Box_Bas(kx,j)
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(Ham_Mat, STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "Ham_Mat deallocation request denied."
     STOP
  ENDIF
  !
END SUBROUTINE BOXDIAG
