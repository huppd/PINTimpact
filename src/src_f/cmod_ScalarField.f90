!>  Modul: cmod_ScalarVector
!!
!! Mature
!! Impact functions \c Pimpact::ScalarVector e.g. scales, norms ...
!! \author huppd
module cmod_ScalarField

  use iso_c_binding
  use mpi

  implicit none

contains



  !> \brief init vector field with 2d Poiseuille flow in x-direction
  !!
  !! \f[ u(x) = x*( L - x )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleX(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleX' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(1):(N(1)+bU(1)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(i)*( L - x(i) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleX



  !> \brief init vector field with 2d Poiseuille flow in y-direction
  !!
  !! \f[ u(y) = y*( L - y )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleY(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleY' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(2):(N(2)+bU(2)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(j)*( L - x(j) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleY



  !> \brief init vector field with 2d Poiseuille flow in y-direction
  !!
  !! \f[ u(z) = y*( L - z )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleZ(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleZ' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(3):(N(3)+bU(3)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(k)*( L - x(k) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleZ



  !> \brief init vector field with constant gradient in x-direction
  !!
  !! \f[ u(x) = \frac{x}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradX( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradX' )

    implicit none

    integer(c_int), intent(in)  ::  N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: L

    real(c_double), intent(in)  :: x( bL(1):(N(1)+bU(1)) )

    real(c_double), intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)  :: alpha

    integer(c_int)              ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) =  alpha*( x(i) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradX



  !> \brief init vector field with constant gradient in y-direction
  !!
  !! \f[ u(y) = \frac{y}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradY( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradY' )

    implicit none

    integer(c_int), intent(in)  ::  N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: L

    real(c_double), intent(in)  :: x( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)  :: alpha

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = alpha*( x(j) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradY



  !> \brief init vector field with constant gradient in z-direction
  !!
  !! \f[ u(z) = \frac{z}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradZ( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradZ' )

    implicit none

    integer(c_int), intent(in)   ::  N(3)

    integer(c_int), intent(in)   :: bL(3)
    integer(c_int), intent(in)   :: bU(3)

    integer(c_int), intent(in)   :: SS(3)
    integer(c_int), intent(in)   :: NN(3)

    real(c_double), intent(in)   :: L

    real(c_double), intent(in)   :: x( bL(3):(N(3)+bU(3)) )

    real(c_double),  intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)   :: alpha

    integer(c_int)               ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = alpha*( x(k) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradZ



  !> \brief init forcing point
  !!
  !! \f[ phi(x) = amp* \exp( - \frac{||x-xc||_2^2}{sig^2} ) \f]
  subroutine SF_init_Vpoint(  &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      x1,                     &
      x2,                     &
      x3,                     &
      xc,                     &
      amp,                    &
      sig,                    &
      phi ) bind ( c, name='SF_init_Vpoint' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)


    real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )
    real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )
    real(c_double), intent(in)     :: x3( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(in)     :: xc(3)

    real(c_double), intent(in)     :: sig(3)
    real(c_double), intent(in)     :: amp

    real(c_double),  intent(inout) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    integer(c_int)               ::  i, j, k


    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        do i = SS(1), NN(1)
          phi(i,j,k) = amp*exp( -( (x1(i)-xc(1) )/sig(1) )**2 -( ( x2(j)-xc(2) )/sig(2) )**2 -( ( x3(k)-xc(3) )/sig(3) )**2 )
        end do
      end do
    end do

  end subroutine SF_init_Vpoint



end module cmod_ScalarField
