PROGRAM B1_B2_INTEGRAL
  !
  !	PLOT OF THE ANALYTICAL B1 AND B2 AS A FUNCTION OF E => F(E) = |<Psi_b|O(x)|Psi_k>|^2 where O(x) = x or x^2
  !
  !	E_max: maximum energy (related to the continuum wave)
  !	E_b:  bound state energy
  !
  !
  ! by LauPK
  !
  USE nrtype
  !
  IMPLICIT NONE
  !
  !VARIABLES and FLAGS
  INTEGER(KIND = I4B) :: NdE, Ierr, NE, I
  REAL(KIND = DP) :: E_b, dE, h_sq, h_sq_o_mu, h_sq_o_mu_div, red_mass, E_max
  !
  REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Func_of_E, Func_2_of_E, E
  !
  !
  ! NAMELIST DEFINITIONS
  NAMELIST/INP_nrg/   E_b, E_max, NdE
  NAMELIST/INP_mass/  red_mass
  !
  !
  ! READING INPUT
  !
   READ(UNIT=*,NML=INP_nrg)
   READ(UNIT=*,NML=INP_mass)
  !
  ALLOCATE(Func_of_E(1:NdE), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Func_of_E allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(Func_2_of_E(1:NdE), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Func_2_of_E allocation request denied."
     STOP
  ENDIF
  !
  ALLOCATE(E(1:NdE), STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "E allocation request denied."
     STOP
  ENDIF
  !
  ! Energy discretization
  dE = E_max/REAL(NdE, DP)
  E = dE*REAL((/ (I, I = 0, NdE) /),DP)
  !
  ! Costants
  h_sq = 41.471291_DP
  h_sq_o_mu = h_sq/red_mass
  h_sq_o_mu_div = 2.0_DP/h_sq_o_mu
  !
  OPEN(UNIT = 70, FILE = 'B1_analytical.dat', STATUS = "UNKNOWN", ACTION = "WRITE")
  WRITE(70,*) "# E       F(E)     for B1"
  OPEN(UNIT = 80, FILE = 'B2_analytical.dat', STATUS = "UNKNOWN", ACTION = "WRITE")
  WRITE(80,*) "# E       F(E)     for B2"
  !
  DO NE = 1, NdE
     Func_of_E(NE)= (2.0_DP*h_sq_o_mu*sqrt(E_b*E(NE)))/(E_b + E(NE))**2
     Func_2_of_E(NE)= ( sqrt(16.0_DP*E_b*h_sq_o_mu_div)*E(NE) )/( (h_sq_o_mu_div**2)*(E_b + E(NE))**3 )
     !
     WRITE(70,*) E(NE), Func_of_E(NE)**2
     WRITE(80,*) E(NE), Func_2_of_E(NE)**2
     !
  ENDDO
  !
  CLOSE(UNIT = 70)
  CLOSE(UNIT = 80)
  !
  PRINT*, "E grid size: ", SIZE(E), " F(E) size: ", SIZE(Func_of_E)
  !
  DEALLOCATE(E, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "E deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Func_of_E, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Func_of_E deallocation request denied."
     STOP
  ENDIF
  !
  DEALLOCATE(Func_2_of_E, STAT = Ierr)    
  IF (Ierr /= 0) THEN
     PRINT*, "Func_of_E deallocation request denied."
     STOP
  ENDIF
  !
  STOP "DON'T PANIC BABY..."
  !
END PROGRAM B1_B2_INTEGRAL
