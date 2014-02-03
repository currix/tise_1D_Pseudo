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
  !
  INTEGER(KIND = I4B) :: n
  REAL(KIND = DP) ::  an, t1, t2
  !
  t1 = 0.0_DP; t2 = 0.0_DP
  !
  DO n = 1, ndim-1
     an = SQRT(REAL(n,DP))
     t1 = t1 + an*Avec_THO(N+1,indho)*SQRT(der_S_x(indx))*THO_Bas(indx,n)
  ENDDO
  !
  DO n = 0, ndim-1
     an = SQRT(REAL(n+1,DP))
     t2 = t2 + an*Avec_THO(n+1,indho)*SQRT(der_S_x(indx))*THO_Bas(indx,n+2)
  ENDDO
  !
  Avec_THO_Der_x(indx,indho) = apar*(t1 - t2)/SQRT(2.0_DP)
  !
  RETURN
END SUBROUTINE WFP
