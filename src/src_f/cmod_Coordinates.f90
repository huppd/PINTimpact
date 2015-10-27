!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************


!> \brief module providing routine, that initializes the coordinates
!!  by using get_coords
module cmod_Coordinates

  use iso_c_binding

  implicit none

contains


  !> \brief sets xx and dx equidistant 
  !!
  !! \param[in] Lmax Li
  !! \param[in] iimax Mi
  !! \param[in] ii0U index
  !! \param[in] ii0L index
  !! \param[out] xx
  !! \param[out] dx
  subroutine coord_equi(  &
      Lmax,               &
      iimax,              &
      ii,                 &
      xx,                 &
      dx )

    implicit none

    real(c_double), intent(in)    ::  Lmax
    real(c_double), intent(in)    ::  iimax

    real(c_double), intent(in)    ::  ii
    real(c_double), intent(out)   ::  xx
    real(c_double), intent(out)   ::  dx

    xx = (ii-1.)*Lmax/(iimax-1.)
    dx =         Lmax/(iimax-1.)

  end subroutine coord_equi


  !> \brief sets xx and dx so that one gets a nice grid streching
  !!
  !! \note if ii0U==ii0L equidastant grid is used
  !! \param[in] Lmax Li
  !! \param[in] iimax Mi
  !! \param[in] ii0U index
  !! \param[in] ii0L index
  !! \param[out] xx
  !! \param[out] dx
  subroutine coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
    ! (sample subroutine)

    implicit none

    real(c_double), intent(in)    ::  Lmax
    real(c_double), intent(in)    ::  iimax

    real(c_double), intent(in)    ::  ii0L
    real(c_double), intent(in)    ::  ii0U

    real(c_double), intent(in)    ::  ii
    real(c_double), intent(out)   ::  xx
    real(c_double), intent(out)   ::  dx

    real(c_double) ::  yy, cmin, cmax


    if( ii0U == ii0L ) then
      ! equidistant grid:
      xx = (ii-1.)*Lmax/(iimax-1.)
      dx =         Lmax/(iimax-1.)
    else
      cmax =  TANH(ii0U)
      cmin = -TANH(ii0L)

      ! coordinate transformation (index i=1. is the origin):
      ! y = (i-1.)/(imax-1.)*(cmax-cmin) + cmin
      yy = (ii-1.)/(iimax-1.)*(cmax-cmin) + cmin

      ! mapping funktion f(yy)
      ! x = L * (f(y)-f(cmin)) / (f(cmax)-f(cmin))
      xx = Lmax * (atanh(yy)-atanh(cmin)) / (ii0U+ii0L)

      ! dx/di = L / (f(cmax)-f(cmin)) * dy/di                 * df(y(i))/dy
      !       = L / (f(cmax)-f(cmin)) * (cmax-cmin)/(imax-1.) * df(y(i))/dy
      dx = Lmax / (atanh(cmax)-atanh(cmin)) * (cmax-cmin)/(iimax-1.) * 1./(1.-yy**2)
    end if


  end subroutine coord_tanh



  !> \brief here can the user sets the coordinates according to his favorite stretching
  !! \realtes GridCoordinatesGlobal
  !! \note: - local processor-block coordinates and grid spacings are
  !!          automatically derived from global grid
  !!        - dy3p = dy3w = 1. for 2D (may simplify programming)
  !!        - ensure that for all i
  !!             y1p(i) < y1p(i+1)
  !!             y1u(i) < y1u(i+1)
  !!             y1p(i) < y1u(i) < y1p(i+1)
  !!                etc.
  !!        - code is tested only for
  !!             y1p(1 ) = 0.
  !!             y1p(M1) = L1
  !!                etc.
  !! \deprecated
  subroutine PI_getGlobalCoordinates( &
      stretchType,                    &
      L,                              &
      M,                              &
      y_origin,                       &
      ys,                             &
      yv,                             &
      dys,                            &
      dyv                             &
      ) bind ( c, name='PI_getGlobalCoordinates' )


    implicit none

    integer(c_int),intent(in)   :: stretchType

    real(c_double),intent(in)   :: L
    integer(c_int),intent(in)   :: M

    real(c_double),intent(in)   :: y_origin

    real(c_double),intent(inout):: ys(1:M)

    real(c_double),intent(inout):: yv(0:M)

    real(c_double),intent(inout):: dys(1:M)

    real(c_double),intent(inout):: dyv(0:M)

    integer(c_int) ::  i


    !--- specify global coordinates and grid spacings ---
    do i = 1, M

      if( 0==stretchType ) then
        call coord_equi( L, REAL(M), REAL(i), ys(i), dys(i) )
      end if

    end do

    do i = 0, M

      if( 0==stretchType ) then
        call coord_equi( L, REAL(M), REAL(i)+0.5, yv(i), dyv(i) )
      end if

    end do

    ys = ys - y_origin
    yv = yv - y_origin

  end subroutine PI_getGlobalCoordinates



  !pgi$g unroll = n:8
  !!pgi$r unroll = n:8
  !!pgi$l unroll = n:8
  !> \realtes GridCoordinatesLocal
  subroutine PI_getLocalCoordinates(  &
      L,                              &
      M,                              &
      N,                              &
      bL,bU,                          &
      BC_L_global,                    &
      BC_U_global,                    &
      BC_L,                           &
      BC_U,                           &
      iB,                             &
      ys,                             &
      yv,                             &
      xs,                             &
      xv,                             &
      dxs,                            &
      dxv ) bind(c,name='PI_getLocalCoordinates')

    implicit none


    real(c_double), intent(in)   :: L

    integer(c_int), intent(in)   :: M

    integer(c_int), intent(in)   :: N

    integer(c_int), intent(in)   :: bL
    integer(c_int), intent(in)   :: bU

    integer(c_int), intent(in)   :: BC_L_global
    integer(c_int), intent(in)   :: BC_U_global

    integer(c_int), intent(in)   :: BC_L
    integer(c_int), intent(in)   :: BC_U

    integer(c_int), intent(in)   :: iB

    real(c_double),intent(in)    ::  ys(1:M)
    real(c_double),intent(in)    ::  yv(0:M)

    real(c_double),intent(inout) ::  xs(bL:(N+bU))
    real(c_double),intent(inout) ::  xv(bL:(N+bU))

    real(c_double),intent(inout) ::  dxs(1:N)
    real(c_double),intent(inout) ::  dxv(0:N)

    integer(c_int)               ::  i, iiShift

    real(c_double)               ::  ysR(bL:(M+bU)), yvR(bL:(M+bU))


    xs  = 0.
    xv  = 0.

    dxs = 0.
    dxv = 0.

    ysR = 0.
    yvR = 0.

    ysR(1:M) = ys(1:M)
    yvR(0:M) = yv(0:M)

    !--- periodic BC ----------------------------------------------------------------------------
    if (BC_L_global == -1) then
      ysR(bL:0) = ysR((M-1+bL):(M-1)) - L
      yvR(bL:0) = yvR((M-1+bL):(M-1)) - L

      ysR(M:(M+bU)) = ysR(1:(1+bU)) + L
      yvR(M:(M+bU)) = yvR(1:(1+bU)) + L
    end if


    !--- symmetic BC or fixed walls -------------------------------------------------------------
    if( BC_L_global == -2 .or. BC_L_global > 0 ) then
      if( BC_L_global > 0 ) then
        ysR(bL: 0) = 2.*ysR(1) - yvR((2-bL):2:-1)
        yvR(bL:-1) = 2.*ysR(1) - yvR((1-bL):2:-1)
      else
        ysR(bL: 0) = 2.*ysR(1) - ysR((2-bL):2:-1)
        yvR(bL: 0) = 2.*ysR(1) - yvR((1-bL):1:-1)
      end if
    end if

    if (BC_U_global == -2 .or. BC_U_global > 0) then
      if (BC_L_global > 0 ) then
        ysR((M+1):(M+bU)) = 2.*ysR(M) - ysR((M-1):(M  -bU):-1)
        yvR((M+1):(M+bU)) = 2.*ysR(M) - yvR((M-2):(M-1-bU):-1)
      else
        ysR((M+1):(M+bU)) = 2.*ysR(M) - ysR((M-1):(M  -bU):-1)
        yvR( M   :(M+bU)) = 2.*ysR(M) - yvR((M-1):(M-1-bU):-1)
      end if
    end if


    !--- distribution to processor blocks -------------------------------------------------------
    iiShift = (iB-1)*(N-1)

    xs(bL:(N+bU)) = ysR((iiShift+bL):(iiShift+N+bU))
    xv(bL:(N+bU)) = yvR((iiShift+bL):(iiShift+N+bU))


    ! is better for evalution of quadrature weights
    do i = 1, N
      dxs(i) = xv(i)-xv(i-1)
    end do

    if( BC_L > 0 .or. BC_L == -2 ) dxs(1) = xv(1)-xs(1  )
    if( BC_U > 0 .or. BC_U == -2 ) dxs(N) = xs(N)-xv(N-1)

    do i = 1, N-1
      dxv(i) = xs(i+1)-xs(i)
    end do

    dxv(0 ) = dxv(1   ) ! necessary for particle feedback loop
    dxv(N) = dxv(N-1)

    if( BC_L == 0 .or. BC_L == -1 ) dxv(0) = xs(1  )-xs(0)
    if( BC_U == 0 .or. BC_U == -1 ) dxv(N) = xs(N+1)-xs(N)


  end subroutine PI_getLocalCoordinates


end module cmod_Coordinates
