!> IMPACT
!! \author Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)
!! \date Mai 2005 - Dec 2011


!> \brief module providing routine, that initializes the coordinates
module cmod_Coordinates

  use iso_c_binding

  implicit none

contains

  !> \brief extracts local coordinates from global
  !!
  !! \param L domain length
  !! \param M global grid size
  !! \param N local grid size
  !! \param bl left max stencil width
  !! \param bu right max stencil width
  !! \param BC_L_global global lower boundary conditions
  !! \param BC_U_global global upper boundary conditions
  !! \param BC_L local lower boundary conditions
  !! \param BC_U local upper boundary conditiosn
  !! \param iB block coordinate
  !! \param ys global scalar coordinates
  !! \param yv global vector coordinates
  !! \param xs local scalar coordinates
  !! \param xv local vector coordinates
  !! \param dxs local scalar deltas 
  !! \param dxv local vector deltas
  !! \relates GridCoordinatesLocal
  subroutine PI_getLocalCoordinates(  &
      time,                           &
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

    integer(c_int), intent(in)   :: time

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

    real(c_double),intent(  out) ::  xs(bL:(N+bU))
    real(c_double),intent(  out) ::  xv(bL:(N+bU))

    real(c_double),intent(  out) ::  dxs(1:N)
    real(c_double),intent(  out) ::  dxv(0:N)

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


    !--- symmetric BC or fixed walls -------------------------------------------------------------
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
    if( time==1 ) then
      iiShift = (iB-1)*N
    else
      iiShift = (iB-1)*(N-1)
    end if

    xs(bL:(N+bU)) = ysR((iiShift+bL):(iiShift+N+bU))
    xv(bL:(N+bU)) = yvR((iiShift+bL):(iiShift+N+bU))


    ! Is better for evaluation of quadrature weights
    do i = 1, N
      dxs(i) = xv(i)-xv(i-1)
    end do

    if( BC_L > 0 .or. BC_L == -2 ) dxs(1) = xv(1)-xs(1  )
    if( BC_U > 0 .or. BC_U == -2 ) dxs(N) = xs(N)-xv(N-1)

    do i = 1, N-1
      dxv(i) = xs(i+1)-xs(i)
    end do

    dxv(0) = dxv(1   ) ! Necessary for particle feedback loop
    dxv(N) = dxv(N-1)

    if( BC_L == 0 .or. BC_L == -1 ) dxv(0) = xs(1  )-xs(0)
    if( BC_U == 0 .or. BC_U == -1 ) dxv(N) = xs(N+1)-xs(N)


  end subroutine PI_getLocalCoordinates


end module cmod_Coordinates
