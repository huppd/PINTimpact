!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing Helmholtz
module cmod_DivOp

  use iso_c_binding

  implicit none

contains

  !>  \brief computes \f$ \mathrm{div = \nabla\cdot\phi } \f$
  !! \deprecated
  subroutine OP_div(  &
      dimens,         &
      N,              &
      bL,bU,          &
      dL,dU,          &
      SS,NN,          &
      c1,c2,c3,       &
      phiU,           &
      div ) bind (c,name='OP_div')

    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: dL(3)
    integer(c_int), intent(in)  :: dU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: c1(dL(1):dU(1),0:N(1))
    real(c_double), intent(in)  :: c2(dL(2):dU(2),0:N(2))
    real(c_double), intent(in)  :: c3(dL(3):dU(3),0:N(3))

    real(c_double), intent(in)  :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)), 1:3 )

    real(c_double), intent(out) :: div (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    integer(c_int)                ::  i, ii
    integer(c_int)                ::  j, jj
    integer(c_int)                ::  k, kk


    !===========================================================================================================
    if( 3==dimens) then

      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            div(i,j,k) = c1(dL(1),i)*phiU(i+dL(1),j,k,1)
            !pgi$ unroll = n:8
            do ii = dL(1)+1, dU(1)
              div(i,j,k) = div(i,j,k) + c1(ii,i)*phiU(i+ii,j,k,1)
            end do
            !pgi$ unroll = n:8
            do jj = dL(2), dU(2)
              div(i,j,k) = div(i,j,k) + c2(jj,j)*phiU(i,j+jj,k,2)
            end do
            !pgi$ unroll = n:8
            do kk = dL(3), dU(3)
              div(i,j,k) = div(i,j,k) + c3(kk,k)*phiU(i,j,k+kk,3)
            end do
          end do
        end do
      end do

    else

      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            div(i,j,k) = c1(dL(1),i)*phiU(i+dL(1),j,k,1)
            !pgi$ unroll = n:8
            do ii = dL(1)+1, dU(1)
              div(i,j,k) = div(i,j,k) + c1(ii,i)*phiU(i+ii,j,k,1)
            end do
            !pgi$ unroll = n:8
            do jj = dL(2), dU(2)
              div(i,j,k) = div(i,j,k) + c2(jj,j)*phiU(i,j+jj,k,2)
            end do
          end do
        end do
      end do

    end if
    !===========================================================================================================


  end subroutine OP_div


  !>  \brief sets non block boundary conditions to zero
  subroutine OP_extrapolateBC_transp(  &
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
      phi ) bind (c,name='OP_extrapolateBC_transp')


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
            DO ii = 0, dU
              phi(i+ii+1,j,k) = phi(i+ii+1,j,k) - phi(i,j,k)*c(ii,1)/c(-1,1)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(1) > 0) THEN
        i = NB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            !pgi$ unroll = n:8
            DO ii = dL, -1
              phi(i+ii,j,k) = phi(i+ii,j,k) - phi(i,j,k)*c(ii,i)/c(0,i)
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
            !pgi$ unroll = n:8
            DO jj = 0, dU
              phi(i,j+jj+1,k) = phi(i,j+jj+1,k) - phi(i,j,k)*c(jj,1)/c(-1,1)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(2) > 0) THEN
        j = NB(2)
        DO k = SB(3), NB(3)
          DO i = SB(1), NB(1)
            !pgi$ unroll = n:8
            DO jj = dL, -1
              phi(i,j+jj,k) = phi(i,j+jj,k) - phi(i,j,k)*c(jj,j)/c(0,j)
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
            DO kk = 0, dU
              phi(i,j,k+kk+1) = phi(i,j,k+kk+1) - phi(i,j,k)*c(kk,1)/c(-1,1)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1)
          END DO
        END DO
      END IF
      IF (BC_U(3) > 0) THEN
        k = NB(3)
        DO j = SB(2), NB(2)
          DO i = SB(1), NB(1)
            !pgi$ unroll = n:8
            DO kk = dL, -1
              phi(i,j,k+kk) = phi(i,j,k+kk) - phi(i,j,k)*c(kk,k)/c(0,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(0,k)
          END DO
        END DO
      END IF
    END IF
    !--------------------------------------------------------------------------------------------
  end subroutine OP_extrapolateBC_transp


end module cmod_DivOp
