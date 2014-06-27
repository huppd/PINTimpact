!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief user defined variables
MODULE usr_vars
  
  USE mod_dims
  USE mod_vars
  
  IMPLICIT NONE
  
  
  !--- specification of additional, user-defined variables ---
  !
  REAL                   ::  pi
  INTEGER                ::  wallnorm_dir
  INTEGER                ::  amax, bmax
  REAL                   ::  L1c, L2c, L3c
  
#ifdef ALLOC
  REAL                   ::  energy_visc
  REAL                   ::  diss_viscInt_old
  
  REAL   , ALLOCATABLE   ::  energy_stokes(:)
  REAL   , ALLOCATABLE   ::  mass_old     (:)
  
  REAL   , ALLOCATABLE   ::  Epot_diff    (:)
  REAL   , ALLOCATABLE   ::  diss_conc_old(:)



#else
  REAL                   ::  energy_visc
  REAL                   ::  diss_viscInt_old
  
  REAL                   ::  energy_stokes(1:n_conc)
  REAL                   ::  mass_old     (1:n_conc)
  
  REAL                   ::  Epot_diff    (1:n_conc)
  REAL                   ::  diss_conc_old(1:n_conc)

!  REAL                   ::  velp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) !< \brief stored velocity, used as referece to previous period
#endif
  
  !> \{ \brief origin of coordinate system
  REAL                   ::  y1_origin
  REAL                   ::  y2_origin
  REAL                   ::  y3_origin
  !> \}
  
  
END MODULE usr_vars
