!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  

  
  
  
    !> \brief user defines initializes vector field
  !! \todo could be moved to \c cmod_VectorField
  subroutine cinitial_conditions_sca( phi ) bind ( c, name='SF_init_field' )
  ! (basic subroutine)
  
  use iso_c_binding

  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                ::  i, j, k
  
  
  !--- initial conditions for pressure ---
  ! note: - initialization is generally not necessary
  !       - specification provides only an initial guess for the first pressure iteration
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !        grid points in the domain and on the boundary
  !        |       |       |
  ! pre(S1p:N1p,S2p:N2p,S3p:N3p)
  !
  phi = 0.
  
  end subroutine cinitial_conditions_sca
  
  
