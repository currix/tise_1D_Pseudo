SUBROUTINE Ew_Strength(ndim, dim_X, X_grid, aval, avec_X, nstates)
  !
  ! $Id: ew_sum_rule.f90,v 1.3 2013/05/05 20:31:56 curro Exp $
  !
  !     by Currix & LauPK
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
  !
  INTEGER(KIND = I4B) :: Ierr, Ifail
  INTEGER(KIND = I4B) :: KI, n
  REAL(KIND = DP) ::  Ews_temp, Ews_total, Error
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: x2_gstate
  !
  CHARACTER(LEN=65) :: filename
  CHARACTER(LEN=65) :: file
  !
  ALLOCATE(x2_gstate(1:dim_X), STAT = Ierr)
  IF (Ierr /= 0) THEN
     PRINT*, "X2_GSTATE allocation request denied."
     STOP
  ENDIF
  !
  WRITE(*,*) "   "
  WRITE(*,*) "CALCULATING Energy Weighted SUM RULES"
  !    
  DO n = 1, nstates
     !
     ! TO PLOT EW_sum_rule (x operator) AS A FUNCTION OF ENERGY
     file = 'EW_sumrule_X_from_'
     WRITE(filename, '(A,I1,"_N",I2, ".dat")')  TRIM(file), n-1, ndim
     IF ( ndim < 10) THEN 
        WRITE(filename, '(A,I1,"_N",I1, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename, '(A,I1,"_N",I3, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     OPEN(UNIT = 80, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(80,*) "# EW sum rule from ", n-1, "-th state  dim_bas = ", ndim, " Box radius = ", X_grid(dim_X), " fm"
     WRITE(80,*) "# Aval(i), EW_sum_rule (x operator)"
     !
     Ews_total = 0.0_DP
     DO KI = 1, ndim - 2 !to avoid ISQW divergence
        x2_gstate = avec_X(:,n)*X_Grid(:)*avec_X(:,KI)
        Ifail = 0
        !
        Ews_temp = 0.0_DP
        !
        CALL D01GAF(X_Grid, x2_gstate, dim_X, Ews_temp, Error, Ifail)
        !     only to avoid an optimization error
        if (ki.eq.ndim+1) write(*,*) ki,Ews_temp,Ews_total
        Ews_total = Ews_total + Ews_temp*Ews_temp*(aval(KI) - aval(n))
        WRITE(80,11)  Aval(ki), Ews_total/h_sq_over_m
     ENDDO
     !
     WRITE(*,*) "E-W SUM RULE (x operator)   from ", n-1, "-th state  = ", Ews_total/h_sq_over_m
     !
     ! TO PLOT EW_sum_rule (x^2 operator) AS A FUNCTION OF ENERGY
     file = 'EW_sumrule_X2_from_'
     WRITE(filename, '(A,I1,"_N",I2, ".dat")')  TRIM(file), n-1, ndim
     IF ( ndim < 10) THEN 
        WRITE(filename, '(A,I1,"_N",I1, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     IF ( ndim > 99) THEN 
        WRITE(filename, '(A,I1,"_N",I3, ".dat")')  TRIM(file), n-1, ndim
     ENDIF
     OPEN(UNIT = 90, FILE = filename, STATUS = "UNKNOWN", ACTION = "WRITE")
     WRITE(90,*) "# EW sum rule from ", n-1, "-th state   dim_bas = ", ndim, " Box radius = ", X_grid(dim_X), " fm"
     WRITE(90,*) "# Aval(i), EW_sum_rule (x^2 operator)"
     !
     Ews_total = 0.0_DP
     x2_gstate = 0.0_DP
     !
     ! EW sum rule
     DO KI = 1, ndim - 6 !to avoid ISQW divergence
        x2_gstate = avec_X(:,n)*X_Grid(:)*X_Grid(:)*avec_X(:,KI)
        Ifail = 0
        !
        Ews_temp = 0.0_DP
        !
        CALL D01GAF(X_Grid, x2_gstate, dim_X, Ews_temp, Error, Ifail)
        !     only to avoid an optimization error
        if (ki.eq.ndim+1) write(*,*) ki,Ews_temp,Ews_total
        Ews_total = Ews_total + Ews_temp*Ews_temp*(aval(KI) - aval(n))
        WRITE(90,11)  Aval(ki), Ews_total/h_sq_over_m
     ENDDO
     !
     ! EW test calculation
     x2_gstate = 0.0_DP
     x2_gstate = avec_X(:,n)*X_Grid(:)*X_Grid(:)*avec_X(:,n)
     Ifail = 0
     !
     Ews_temp = 0.0_DP
     !
     CALL D01GAF(X_Grid, x2_gstate, dim_X, Ews_temp, Error, Ifail)
     !
     WRITE(*,*) "E-W SUM RULE (x^2 operator) from ", n-1, "-th state  = ", Ews_total/h_sq_over_m, &
          " <-------> test value ", (Ews_temp*2.0_DP)
     !
     CLOSE(UNIT = 80)
     CLOSE(UNIT = 90)
     !
  ENDDO
  !
11 FORMAT (1X,E16.8,1X,E17.8)
  !
  DEALLOCATE(x2_gstate, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "x2_gstate deallocation request denied."
     STOP
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE EW_STRENGTH
