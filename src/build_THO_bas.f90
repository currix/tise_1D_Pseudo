SUBROUTINE THO_1D_BASIS(apar, NdimH, Iprint)
  !     
  !     COMPUTES A 1D ANALYTIC THO BASIS
  !
  !     INPUT  :: 
  !               apar    --> LENGTH SCALE OF THE PROBLEM
  !               NdimH   --> DIMENSION + 1 OF THE TRANSFORMED HARMONIC BASIS
  !               THO constants (INCLUDE THEM AS ARGUMENTS)
  !
  !     OUTPUT :: HARBAS  --> MATRIX WITH HARMONIC BASIS
  !
  !     FORMAT :: Iprint  --> VERBOSITY CONTROL
  !
  !     $Id: build_THO_bas.f90,v 1.3 2013/06/13 10:57:54 curro Exp $
  !
  !     by Currix TM.
  !
  USE nrtype
  USE egs_tho_f
  USE constants
  USE lst_param
  USE pot_param
  !
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  INTEGER(KIND = I4B), INTENT(IN) :: NdimH, Iprint
  REAL(KIND = DP), INTENT(IN) :: apar 
  !
  !
  ! OTHER VARIABLES
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: THO_norm_test
  REAL(KIND = DP) :: PI14, apar2, ainteg, Error
  INTEGER(KIND = I4B) :: kx, Ifail, Ierr, X_dim
  !
  !
  !fisso il rapporto (automaticamente gamma)
  gamma = ratio/apar ! ratio articolo = 1.65
  !print*, "gamma = ", gamma
  !print*, "b = ", 1.0_DP/apar
  !
  IF (Iprint > 2) PRINT*, "BUILDING THO BASIS"
  !
  ! DEFINING SCALING FUNCTION 
  S_x = 0.0_dp
  der_S_x = 0.0_dp
  der2_S_x = 0.0_dp
  !
  S_x = lst(X_grid)
  !
  S_over_x = sox(X_grid)
  !
  der_S_x = d_lst(X_grid)
  !
  der2_S_x = d2_lst(X_grid)
  !
  PI14 = SQRT(SQRT(PI_D))
  apar2 = apar*apar
  !
  !     HO n = 0
  THO_BAS(:,1) = SQRT(apar)/(PI14)*EXP(-apar2*S_x*S_x/2.0_DP)
  !
  !     HO n = 1
  IF (NdimH > 1) THO_BAS(:,2) = (SQRT(apar)/(PI14))*SQRT(2.0_DP)*apar*S_x*EXP(-apar2*S_x*S_x/2.0_DP)
  !
  !    RECURRENCE RELATION (WATCH OUT THE INDEXES :: MATRIX START AT 1, NOT 0)
  DO kx = 2, NdimH-1
     THO_BAS(:,kx+1) = &
          SQRT(2.0_DP/(1.0_DP*kx))*apar*S_x*THO_BAS(:,kx) - &
          SQRT((1.0_DP*(kx-1))/(1.0_DP*kx))*THO_BAS(:,kx-1)
  ENDDO
  !
  !
  !BUILDING THO BASIS 
!  DO n = 1, NdimH
!        THO_bas(:,n)= SQRT(der_S_x(:))*THO_Bas(:,n) 
!  ENDDO
  !
  !     TESTING NORMALIZATION (DESTROYS HARBAS(X_GRID,DIMENSION + 2 )!!!)
  IF (Iprint > 1) THEN
     !
     ! DEFINE THO_norm_test
     X_dim = SIZE(X_GRID)
     ALLOCATE(THO_norm_test(1:X_dim), STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "THO_norm_test allocation request denied."
        STOP
     ENDIF
     !
     DO kx = 1, NdimH
        THO_norm_test = THO_BAS(:,kx)*THO_BAS(:,kx)*der_S_x(:)
        !        Print*, "der2_S_x(1+dim_X/2) ",  der2_S_x(1+dim_X/2)
        !        Ifail = 0
        !
        !        CALL D01GAF(X_GRID, THO_norm_test, X_dim, ainteg, Error, Ifail)
        !
        !        THO_Bas(:,kx)=THO_Bas(:,kx)/SQRT(ainteg)
        !
        !        THO_norm_test = THO_BAS(:,kx)*THO_BAS(:,kx)
        !        Print*, "der2_S_x(1+dim_X/2) ",  der2_S_x(1+dim_X/2)
        Ifail = 0
        !
        CALL D01GAF(X_GRID, THO_norm_test, X_dim, ainteg, Error, Ifail)
        PRINT*, "THO FUNCTION ", kx, " NORMALIZATION", ainteg
     ENDDO
     DEALLOCATE(THO_norm_test, STAT = Ierr)    
     IF (Ierr /= 0) THEN
        PRINT*, "THO_norm_test deallocation request denied."
        STOP
     ENDIF
  ENDIF
  !     
  IF (Iprint > 2) WRITE(*,*) "DONE"
  !     
END SUBROUTINE THO_1D_BASIS
