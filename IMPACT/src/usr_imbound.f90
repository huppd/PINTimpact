!> \brief module introducing a discrete direct immersed boundaries
!! \author huppd
MODULE usr_imbound

  USE mod_dims
  USE mod_vars
  USE usr_vars
!  USE usr_func

  IMPLICIT NONE

  real,parameter :: amp = 1.0/8.0

  CONTAINS


  !> \todo add obstacle characteristic function

  !> \brief calculate distance to immersed boundary
  !! \author huppd
  !! \param[in] x
  !! \param[in] y
  !! \param[in] z
  FUNCTION distance2ib(x,y,z) RESULT(dis)

  IMPLICIT NONE

  REAL    ,   INTENT(IN)  ::  x,y,z
  REAL                    ::  dis
  REAL                    ::  x0
  REAL                    ::  y0
  REAL                    ::  z0
  REAL                    ::  R
  REAL                    ::  dr
  INTEGER                 ::  m

  ! geometric properties of disc
  x0 = L1/2. + L1*amp*SIN(2.*pi*freq*subtime)
  y0 = L2/4.
  z0 = L3/2.

  R = L1/10.
  dr = SQRT( ( L1/REAL(M1-1) )**2 + ( L2/REAL(M2-1) )**2 )
!  write(*,*) 'dr=',dr

  dis = SQRT( (x0-x)**2 + (y0-y)**2 )

  IF( dis <= R) THEN
    dis = 1.
  ELSE IF( dis-R <= dr) THEN
    dis = (1./( 1. + EXP(1./(-(dis-R)/dr) + 1./(1.-(dis-R)/dr) ) ));
!    dis = (1-(dis-R)/dr)
  ELSE
    dis = 0.
  END IF
!    dis = 0.

  END FUNCTION distance2ib



  !> \brief returns velocity of immersed boundaries at time
  !! \author huppd
  !! \param[in] m velocity direction
  FUNCTION vel_ib(x,y,z,m) RESULT(velocity_ib)

  IMPLICIT NONE

  REAL    ,   INTENT(IN)  ::  x,y,z
  REAL                    ::  velocity_ib
  INTEGER ,   INTENT(IN)  ::  m

  IF( m==1 ) velocity_ib = L1*amp*2.*pi*freq*COS(2.*pi*freq*subtime)
  IF( m==2 ) velocity_ib = 0.
  IF( m==3 ) velocity_ib = 0.

  END FUNCTION vel_ib


END MODULE usr_imbound
