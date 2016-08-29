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



  !>  \brief extrapolates for non-block BC outside value
  subroutine OP_extrapolateBC(  &
      m,                        &
      N,                        &
      bL,                       &
      bU,                       &
      dL,                       &
      dU,                       &
      BCL,                      &
      BCU,                      &
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

    integer(c_int), intent(in)  :: BCL
    integer(c_int), intent(in)  :: BCU

    integer(c_int), intent(in)  :: SB(3)
    integer(c_int), intent(in)  :: NB(3)

    real(c_double), intent(in)  :: c(dL:dU,0:N(m))

    real(c_double), intent(inout):: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )

    integer(c_int) :: i, j, k
    integer(c_int) :: ii,jj,kk

    !--------------------------------------------------------------------------------------------
    IF (m == 1) THEN
      IF( BCL > 0 ) THEN
        i = SB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            !pgi$ unroll = n:8
            phi(i,j,k) = 0.
            DO ii = 0, dU
              phi(i,j,k) = phi(i,j,k) - c(ii,1)*phi(1+ii,j,k)
            END DO
            phi(i,j,k) = phi(i,j,k) / c(-1,1) 
          END DO
        END DO
      END IF
      IF( BCU > 0 ) THEN
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
      IF( BCL > 0 ) THEN
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
      IF( BCU > 0 ) THEN
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
      IF( BCL > 0 ) THEN
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
      IF( BCU > 0 ) THEN
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



  !>  \brief extrapolates for non-block BC outside value
  !!
  !! using Neville-Aitken scheme
  subroutine OP_extrapolateBC2(  &
      m,                         &
      N,                         &
      bL,                        &
      bU,                        &
      dL,                        &
      dU,                        &
      BCL,                       &
      BCU,                       &
      SB,                        &
      NB,                        &
      xu,                        &
      phi ) bind (c,name='OP_extrapolateBC2')


    implicit none

    integer(c_int), intent(in)    :: m

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: dL
    integer(c_int), intent(in)    :: dU

    integer(c_int), intent(in)    :: BCL
    integer(c_int), intent(in)    :: BCU

    integer(c_int), intent(in)    :: SB(3)
    integer(c_int), intent(in)    :: NB(3)

    !real(c_double), intent(in)    :: c(dL:dU,0:N(m))

    real(c_double), intent(in)    :: xu(bl(m):N(m)+bu(m))

    real(c_double), intent(inout) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )


    real(c_double)                :: y(1:dU)
    real(c_double)                :: t(1:dU)
    real(c_double)                :: x

    integer(c_int) :: i, j, k
    integer(c_int) :: ii,jj,kk

    !--------------------------------------------------------------------------------------------
    !write(*,*) dl
    !write(*,*) dU
    if( m == 1 ) then
      if( BCL > 0 ) then
        i = SB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            ! load data
            y(1:dU) = phi(i+1:i+dU,j,k)
            t(1:dU) = xu (i+1:i+dU)
            x = xu(i)
            ! Lagrange extrapolation
            do ii = 1, dU
              do kk = ii-1, 1, -1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
      IF( BCU > 0 ) THEN
        i = NB(1)
        DO k = SB(3), NB(3)
          DO j = SB(2), NB(2)
            ! load data
            y(1:dL) = phi(i-dL:i-1,j,k)
            t(1:dL) = xu (i-dL:i-1)
            x = xu(i)
            ! Lagrange extrapolation
            do ii = 1, dL
              do kk=ii-1,1,-1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
    end if
    !--------------------------------------------------------------------------------------------
    if( m == 2 ) then
      if( BCL > 0 ) then
        j = SB(2)
        DO k = SB(3), NB(3)
          DO i = SB(1), NB(1)
            ! load data
            y(1:dU) = phi(i,j+1:j+dU,k)
            t(1:dU) = xu (j+1:j+dU)
            x = xu(j)
            ! Lagrange extrapolation
            do ii = 1, dU
              do kk = ii-1, 1, -1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
      if( BCU > 0 ) then
        j = NB(2)
        DO k = SB(3), NB(3)
          DO i = SB(1), NB(1)
            ! load data
            y(1:dL) = phi(i,j-dL:j-1,k)
            t(1:dL) = xu (j-dL:j-1)
            x = xu(j)
            ! Lagrange extrapolation
            do ii = 1, dL
              do kk=ii-1,1,-1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
    end if
    !--------------------------------------------------------------------------------------------
    if( m == 3 ) then
      if( BCL > 0 ) THEN
        k = SB(3)
        DO j = SB(2), NB(2)
          DO i = SB(1), NB(1)
            ! load data
            y(1:dU) = phi(i,j,k+1:k+dU)
            t(1:dU) = xu (k+1:k+dU)
            x = xu(k)
            ! Lagrange extrapolation
            do ii = 1, dU
              do kk = ii-1, 1, -1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
      if( BCU > 0 ) then
        k = NB(3)
        DO j = SB(2), NB(2)
          DO i = SB(1), NB(1)
            ! load data
            y(1:dL) = phi(i,j,k-dL:k-1)
            t(1:dL) = xu (k-dL:k-1)
            x = xu(k)
            ! Lagrange extrapolation
            do ii = 1, dL
              do kk=ii-1,1,-1
                y(kk) = y(kk+1)+(y(kk+1)-y(kk)) * (x-t(ii))/(t(ii)-t(kk));
              end do
            end do
            phi(i,j,k) = y(1);
          END DO
        END DO
      end if
    end if
    !--------------------------------------------------------------------------------------------
  end subroutine OP_extrapolateBC2



end module cmod_GradOp
