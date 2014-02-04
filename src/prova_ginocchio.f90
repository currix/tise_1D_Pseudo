PROGRAM NEWTON_RAPHSON

  !
  !     GINOCCHIO 1D POTENTIAL
  !
  IMPLICIT NONE
  !
  REAL :: nu, lambda, Delta_X, X_min, X_max
  REAL :: y, y0, y1,ym,s0, s1, sm, tol, func, func_der, f0, f1, fm
  INTEGER :: count, MaxIt, dim_X, i, kx, Ierr, Iter
  REAL, DIMENSION(:), ALLOCATABLE :: X_grid
  REAL, DIMENSION(:), ALLOCATABLE :: Potf
  !
  !
  INTERFACE r
     ELEMENTAL FUNCTION r(y,l)
       !
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL, INTENT(IN) :: y, l
       !
       REAL ::  r
     END FUNCTION r
  END INTERFACE r
  !
  INTERFACE r_der
     ELEMENTAL FUNCTION r_der(y,l)
       !
       IMPLICIT NONE
       !
       ! ARGUMENTS
       REAL, INTENT(IN) :: y, l
       !
       REAL ::  r_der
     END FUNCTION r_der
  END INTERFACE r_der
  !
  !
  Dim_X = 501
  !
  ALLOCATE(X_grid(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "X_grid allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Potf(1:dim_X), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Potf allocation request denied."
     STOP
  ENDIF
  !
  ! DEFINE X GRID  
  !
  X_min = -10.0
  X_max = 10.0
  !
  Delta_X = (X_max-X_min)/REAL(dim_X - 1)  
  X_grid = X_min + Delta_X*REAL((/ (I, I = 0, dim_X-1) /))
  !
!print*, X_grid
  !
  PRINT*, "ENTER NU"
  READ*, nu
  PRINT*, "ENTER LAMBDA"
  READ*, lambda
  !
  grid : DO i = 1, dim_X
     count = 0
     Iter = 0
     MaxIt = 20
     tol = 0.01
     y0 = 0.5 
     !y1 = 1.0 
     !
!BISECTION METHOD
!     f0 = r(y0,lambda) - X_grid(i)
!     check_Ea: IF (f0 == 0.0) THEN !checks if y0 is a solution
!        Potf(i) = 0.0
!        CYCLE grid
!     ENDIF check_Ea
!     S0 = abs(F0)/F0
!     !
!     !
!     f1 =  r(y1,lambda) - X_grid(i) !Fb
!     check_Eb : IF(F1 == 0.0) THEN !checks if y1 is a solution
!        Potf(i) = 0.0
!        CYCLE grid
!     ENDIF check_Eb
!     S1 = abs(F1)/F1
!     !
!     !loop che dimezza l'intervallo
!     bisection_method : DO
!        ym = (y0 + y1)/2.0
!        fm =  r(ym,lambda) - X_grid(i)
!        !
!        check_m : IF(Fm == 0.0) THEN !checks if ym is  a solution
!           y = ym
!           Potf(i) = -(lambda**2)*nu*(nu+1.0)*(1.0-y**2) &
!                + ((1.0-y**2)*(1.0-lambda**2)/4.0)* &
!                (2.0 - (7.0-lambda**2)*(y**2) + 5.0*(1.0-lambda**2)*(y**4) )
!           CYCLE grid
!        ENDIF check_m
!        Sm = abs(fm)/fm
!        !print*,"f1 ", f1, "fm ", fm, "f0", f0
!        !
!        IF( S0 /= Sm) THEN
!           y1 = ym
!           F1 =  r(y1,lambda) - X_grid(i)
!           S1 = abs(F1)/F1
!        ELSE
!           y0 = ym
!           F0 =  r(y0,lambda) - X_grid(i)
!           S0 =  abs(F0)/F0
!        ENDIF
!        !
!        IF (ABS(y0-y1) < tol )THEN
!           y = y0
!           Potf(i) = -(lambda**2)*nu*(nu+1.0)*(1.0-y**2) &
!                + ((1.0-y**2)*(1.0-lambda**2)/4.0)* &
!                (2.0 - (7.0-lambda**2)*(y**2) + 5.0*(1.0-lambda**2)*(y**4) )
!           CYCLE grid
!        ENDIF
!     ENDDO bisection_method
 !
       NewtonRaphson : DO
          func = r(y0,lambda) - X_grid(i)
          func_der = r_der(y0,lambda)
          y1 = y0 - func/func_der
          print*, "y0 ", y0, "func/func_der", func/func_der, "func_der", func_der,  "func", func
          !
          IF( ABS(y1 - y0) < tol) THEN
             y = y1
          ENDIF
          !
          y0 = y1
          Iter = Iter +1
          !
          IF( Iter > MaxIt) THEN
             PRINT*, "Too much iterations."
             !y = 1.0
             EXIT NewtonRaphson
          ENDIF
       ENDDO NewtonRaphson
  !
           Potf(i) = -(lambda**2)*nu*(nu+1.0)*(1.0-y**2) &
                + ((1.0-y**2)*(1.0-lambda**2)/4.0)* &
                (2.0 - (7.0-lambda**2)*(y**2) + 5.0*(1.0-lambda**2)*(y**4) )
ENDDO grid
!
!print*, X_grid
!
OPEN(UNIT = 74, FILE = 'test_ginocchio.dat', STATUS = "UNKNOWN", ACTION = "WRITE")
DO kx = 1, dim_X
  WRITE(74,11) X_grid(kx), Potf(kx)
ENDDO
CLOSE(UNIT = 74)
!
11 FORMAT (1X,E14.6,1X,E14.6)
!
DEALLOCATE(X_grid, STAT = Ierr)    
IF (Ierr /= 0) THEN
  PRINT*, "X_grid deallocation request denied."
  STOP
ENDIF
!
DEALLOCATE(Potf, STAT = Ierr)    
IF (Ierr /= 0) THEN
  PRINT*, "Potf deallocation request denied."
  STOP
ENDIF
!
!
PRINT*, "SAYONARA BABY..."
END PROGRAM
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELEMENTAL FUNCTION r(y,l)
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  REAL, INTENT(IN) :: y, l
  !
  REAL ::  r
  REAL :: aux 
  aux = SQRT(l**2 -1.0)
  !
  IF( y>=1) THEN
     r = (l**(-2))* (100 + aux*ATAN(aux*y) )
     RETURN
  ELSE IF (y<=-1) THEN
     r = (l**(-2))* (  -100 + aux*ATAN(aux*y) )
     RETURN
  ELSE
     r = (l**(-2))* ( ATANH(y) + aux*ATAN(aux*y) )
     RETURN
  ENDIF
  !
  RETURN
END FUNCTION r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELEMENTAL FUNCTION r_der(y,l)
  !
  IMPLICIT NONE
  !
  ! ARGUMENTS
  REAL, INTENT(IN) :: y, l
  !
  REAL ::  r_der
  REAL :: aux 
  aux = SQRT(l**2 -1.0)
  !
  IF( y>=1) THEN
   r_der = (l**(-2))* ( (aux**2)/(1.0 + (aux*y)**2 ) )
     RETURN
  ELSE IF (y<=-1) THEN
   r_der = (l**(-2))* ( (aux**2)/(1.0 + (aux*y)**2 ) )
     RETURN
  ELSE
   r_der = (l**(-2))* ( 1.0/(1.0-y**2)  + (aux**2)/(1.0 + (aux*y)**2 ) )
     RETURN
  ENDIF
  !
  RETURN
END FUNCTION r_der
