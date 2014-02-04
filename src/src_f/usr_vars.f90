!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief user defined variables
module usr_vars
  
  use mod_dims
  use mod_vars
  
  implicit none
  
  
  !--- specification of additional, user-defined variables ---
  !
  real                   ::  pi
  integer                ::  wallnorm_dir
  integer                ::  amax, bmax
  real                   ::  L1c, L2c, L3c
  
#ifdef ALLOC
  real                   ::  energy_visc
  real                   ::  diss_viscInt_old
  
  real   , allocatable   ::  energy_stokes(:)
  real   , allocatable   ::  mass_old     (:)
  
  real   , allocatable   ::  Epot_diff    (:)
  real   , allocatable   ::  diss_conc_old(:)



#else
  real                   ::  energy_visc
  real                   ::  diss_viscInt_old
  
  real                   ::  energy_stokes(1:n_conc)
  real                   ::  mass_old     (1:n_conc)
  
  real                   ::  Epot_diff    (1:n_conc)
  real                   ::  diss_conc_old(1:n_conc)

!  REAL                   ::  velp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3) !< \brief stored velocity, used as referece to previous period
#endif
  
  !> \{ \brief origin of coordinate system
  real                   ::  y1_origin
  real                   ::  y2_origin
  real                   ::  y3_origin
  !> \}
  
  
end module usr_vars
