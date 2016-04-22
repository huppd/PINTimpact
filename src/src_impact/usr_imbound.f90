!> \brief module introducing a discrete direct immersed boundaries
module usr_imbound

  use mod_dims
  use mod_vars
  use usr_vars
!  USE usr_func

  implicit none

  real,parameter :: amp = 1.0/8.0

  contains



  !> \brief calculate distance to immersed boundary
  !! \param[in] x
  !! \param[in] y
  !! \param[in] z
  function distance2ib(x,y,z) result(dis)

  implicit none

  real    ,   intent(in)  ::  x,y,z
  real                    ::  dis
  real                    ::  x0
  real                    ::  y0
  real                    ::  z0
  real                    ::  R
  real                    ::  dr
  integer                 ::  m

  ! geometric properties of disc
  x0 = L1/2. + L1*amp*SIN(2.*pi*freq*subtime)
  y0 = L2/4.
  z0 = L3/2.

  R = L1/10.
  dr = SQRT( ( L1/REAL(M1-1) )**2 + ( L2/REAL(M2-1) )**2 )
!  write(*,*) 'dr=',dr

  dis = SQRT( (x0-x)**2 + (y0-y)**2 )

  if( dis <= R) then
    dis = 1.
  else if( dis-R <= dr) then
    dis = (1./( 1. + EXP(1./(-(dis-R)/dr) + 1./(1.-(dis-R)/dr) ) ));
!    dis = (1-(dis-R)/dr)
  else
    dis = 0.
  end if
!    dis = 0.

  end function distance2ib



  !> \brief returns velocity of immersed boundaries at time
  !! \param[in] m velocity direction
  function vel_ib(x,y,z,m) result(velocity_ib)

  implicit none

  real    ,   intent(in)  ::  x,y,z
  real                    ::  velocity_ib
  integer ,   intent(in)  ::  m

  if( m==1 ) velocity_ib = L1*amp*2.*pi*freq*COS(2.*pi*freq*subtime)
  if( m==2 ) velocity_ib = 0.
  if( m==3 ) velocity_ib = 0.

  end function vel_ib


end module usr_imbound
