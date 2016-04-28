!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing Helmholtz
module cmod_GradOp

  use iso_c_binding

  implicit none

contains

  !> \brief computes \f$ \mathrm{grad = \nabla\phi } \f$
  !! is used for inner field
  !! \deprecated
  subroutine OP_S2VOp( &
      dir,             &
      N,               &
      bL,bU,           &
      gL,gU,           &
      SS,NN,           &
      c,               &
      phi,             &
      grad ) bind (c,name='OP_S2VOp')


    implicit none

    integer(c_int), intent(in)  :: dir

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: gL(3)
    integer(c_int), intent(in)  :: gU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: c(gL(dir):gU(dir),0:N(dir))

    real(c_double), intent(in)  :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(out) :: grad( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )


    integer(c_int)               ::  i, ii
    integer(c_int)               ::  j, jj
    integer(c_int)               ::  k, kk


    !===========================================================================================================
    if( 1==dir ) then

      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            grad(i,j,k) = c(gL(1),i)*phi(i+gL(1),j,k)
            !pgi$ unroll = n:8
            do ii = gL(1)+1, gU(1)
              grad(i,j,k) = grad(i,j,k) + c(ii,i)*phi(i+ii,j,k)
            end do
          end do
        end do
      end do

    end if
    !===========================================================================================================
    if( 2==dir ) then

      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            grad(i,j,k) = c(gL(2),j)*phi(i,j+gL(2),k)
            !pgi$ unroll = n:8
            do jj = gL(2)+1, gU(2)
              grad(i,j,k) = grad(i,j,k) + c(jj,j)*phi(i,j+jj,k)
            end do
          end do
        end do
      end do

    end if

    !===========================================================================================================
    if( 3==dir ) then

      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            grad(i,j,k) = c(gL(3),k)*phi(i,j,k+gL(3))
            !pgi$ unroll = n:8
            do kk = gL(3)+1, gU(3)
              grad(i,j,k) = grad(i,j,k) + c(kk,k)*phi(i,j,k+kk)
            end do
          end do
        end do
      end do

    end if
    !===========================================================================================================

  end subroutine OP_S2VOp



  !>  \brief sets non block boundary conditions to zero(not joust corners)
  !!
  !! should be obsolete by now
  !! \param[in] N local grid size
  !! \param[in] bl lower storage offset
  !! \param[in] bu upper storage offset
  !! \param[in] BC_L lower boundary conditions
  !! \param[in] BC_U upper boundary condtions
  !! \param[in] SB start index including boundaries
  !! \param[in] NB end index including boundaries
  !! \param[out] grad field
  subroutine OP_SetBCZero(  &
      N,                    &
      bL,                   &
      bU,                   &
      BC_L,                 &
      BC_U,                 &
      SB,                   &
      NB,                   &
      grad ) bind ( c, name='OP_SetBCZero' )

    implicit none

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: BC_L(3)
    integer(c_int), intent(in)  :: BC_U(3)

    integer(c_int), intent(in)  :: SB(3)
    integer(c_int), intent(in)  :: NB(3)

    real(c_double), intent(inout) :: grad(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    !--- Boundary conditions --------------------------------------------------------------------
    if (BC_L(1) > 0) grad(SB(1)      ,SB(2):NB(2),SB(3):NB(3)) = 0.
    if (BC_U(1) > 0) grad(      NB(1),SB(2):NB(2),SB(3):NB(3)) = 0.
    if (BC_L(2) > 0) grad(SB(1):NB(1),SB(2)      ,SB(3):NB(3)) = 0.
    if (BC_U(2) > 0) grad(SB(1):NB(1),      NB(2),SB(3):NB(3)) = 0.
    if (BC_L(3) > 0) grad(SB(1):NB(1),SB(2):NB(2),SB(3)      ) = 0.
    if (BC_U(3) > 0) grad(SB(1):NB(1),SB(2):NB(2),      NB(3)) = 0.

  end subroutine OP_SetBCZero



  !>  \brief sets non block boundary conditions to zero
  subroutine OP_extrapolateBC(  &
      m,                        &
      N,                        &
      bL,                       &
      bU,                       &
      dL,                       &
      dU,                       &
      BC_L,                     &
      BC_U,                     &
      SB,                       &
      NB,                       &
      c,                        &
      phi ) bind (c,name='OP_extrapolateBC')


    implicit none

    integer(c_int), intent(in)  :: m

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: dL
    integer(c_int), intent(in)  :: dU

    integer(c_int), intent(in)  :: BC_L(3)
    integer(c_int), intent(in)  :: BC_U(3)

    integer(c_int), intent(in)  :: SB(3)
    integer(c_int), intent(in)  :: NB(3)

    real(c_double), intent(in)  :: c(dL:dU,0:N(m))

    real(c_double), intent(inout):: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )

    integer(c_int) :: i, j, k
    integer(c_int) :: ii,jj,kk

    !--------------------------------------------------------------------------------------------
    IF (m == 1) THEN
      IF (BC_L(1) > 0) THEN
        i = SB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            !pgi$ unroll = n:8
            phi(i,j,k) = 0.
            DO ii = 0, dU
              phi(i,j,k) = phi(i,j,k) - c(ii,1)*phi(1+ii,j,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
            !print *, c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(1) > 0) THEN
        i = NB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            phi(i,j,k) = 0.
            !pgi$ unroll = n:8
            DO ii = dL, -1
              phi(i,j,k) = phi(i,j,k) - c(ii,i)*phi(i+ii,j,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(0,i)
          END DO
        END DO
      END IF
    END IF
    !--------------------------------------------------------------------------------------------
    IF (m == 2) THEN
      IF (BC_L(2) > 0) THEN
        j = SB(2)
        DO k = SB(3), NB(3)
          DO i = SB(1), NB(1)
            phi(i,j,k) = 0.
            !pgi$ unroll = n:8
            DO jj = 0, dU
              phi(i,j,k) = phi(i,j,k) - c(jj,1)*phi(i,1+jj,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(2) > 0) THEN
        j = NB(2)
        DO k = SB(3), NB(3)
          DO i = SB(1), NB(1)
            phi(i,j,k) = 0.
            !pgi$ unroll = n:8
            DO jj = dL, -1
              phi(i,j,k) = phi(i,j,k) - c(jj,j)*phi(i,j+jj,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(0,j)
          END DO
        END DO
      END IF
    END IF
    !--------------------------------------------------------------------------------------------
    IF (m == 3) THEN
      IF (BC_L(3) > 0) THEN
        k = SB(3)
        DO j = SB(2), NB(2)
          DO i = SB(1), NB(1)
            !pgi$ unroll = n:8
            phi(i,j,k) = 0.
            DO kk = 0, dU
              phi(i,j,k) = phi(i,j,k) - c(kk,1)*phi(i,j,1+kk)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(3) > 0) THEN
        k = NB(3)
        DO j = SB(2), NB(2)
          DO i = SB(1), NB(1)
            phi(i,j,k) = 0.
            !pgi$ unroll = n:8
            DO kk = dL, -1
              phi(i,j,k) = phi(i,j,k) - c(kk,k)*phi(i,j,k+kk)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(0,k)
          END DO
        END DO
      END IF
    END IF
    !--------------------------------------------------------------------------------------------
  end subroutine OP_extrapolateBC



  !> \brief computes \f$ \mathrm{grad = \nabla\phi } \f$
  !! is used for inner field
  !! \depcrecated
  subroutine OP_grad( &
      dimens,         &
      N,              &
      bL,bU,          &
      gL,gU,          &
      SU,NU,          &
      SV,NV,          &
      SW,NW,          &
      c1,             &
      c2,             &
      c3,             &
      phi,            &
      grad ) bind (c,name='OP_grad')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: gL(3)
    integer(c_int), intent(in)  :: gU(3)

    integer(c_int), intent(in)  :: SU(3)
    integer(c_int), intent(in)  :: NU(3)

    integer(c_int), intent(in)  :: SV(3)
    integer(c_int), intent(in)  :: NV(3)

    integer(c_int), intent(in)  :: SW(3)
    integer(c_int), intent(in)  :: NW(3)

    real(c_double), intent(in)  :: c1(gL(1):gU(1),0:N(1))
    real(c_double), intent(in)  :: c2(gL(2):gU(2),0:N(2))
    real(c_double), intent(in)  :: c3(gL(3):gU(3),0:N(3))

    real(c_double), intent(in)  :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(out) :: grad( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3 )


    integer(c_int)               ::  i, ii
    integer(c_int)               ::  j, jj
    integer(c_int)               ::  k, kk


    !============================================================================================

    do k = SU(3), NU(3)
      do j = SU(2), NU(2)
        do i = SU(1), NU(1)
          grad(i,j,k,1) = c1(gL(1),i)*phi(i+gL(1),j,k)
          !pgi$ unroll = n:8
          do ii = gL(1)+1, gU(1)
            grad(i,j,k,1) = grad(i,j,k,1) + c1(ii,i)*phi(i+ii,j,k)
          end do
        end do
      end do
    end do

    !============================================================================================

    do k = SV(3), NV(3)
      do j = SV(2), NV(2)
        do i = SV(1), NV(1)
          grad(i,j,k,2) = c2(gL(2),j)*phi(i,j+gL(2),k)
          !pgi$ unroll = n:8
          do jj = gL(2)+1, gU(2)
            grad(i,j,k,2) = grad(i,j,k,2) + c2(jj,j)*phi(i,j+jj,k)
          end do
        end do
      end do
    end do


    !============================================================================================
    if( 3==dimens ) then

      do k = SW(3), NW(3)
        do j = SW(2), NW(2)
          do i = SW(1), NW(1)
            grad(i,j,k,3) = c3(gL(3),k)*phi(i,j,k+gL(3))
            !pgi$ unroll = n:8
            do kk = gL(3)+1, gU(3)
              grad(i,j,k,3) = grad(i,j,k,3) + c3(kk,k)*phi(i,j,k+kk)
            end do
          end do
        end do
      end do

    end if
    !============================================================================================

  end subroutine OP_grad


end module cmod_GradOp
