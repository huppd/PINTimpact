!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE pseudocall(phi)
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  END SUBROUTINE pseudocall
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE pseudocall_int(phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  END SUBROUTINE pseudocall_int
  
  
  
  
