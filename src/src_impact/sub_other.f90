!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine pseudocall(phi)
  
  implicit none
  
  real   , intent(in)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  end subroutine pseudocall
  
  
  
  
  
  
  
  
  
  
  
  subroutine pseudocall_int(phi)
  
  implicit none
  
  integer, intent(in)    ::  phi
  
  !---------------------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Pseudo-CALL fuer mod_exchange. Darf auf keinen Fall vom Compiler geinlined werden!                   !
  !---------------------------------------------------------------------------------------------------------------------!
  
  
  end subroutine pseudocall_int
  
  
  
  
