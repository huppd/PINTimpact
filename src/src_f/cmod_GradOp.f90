!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing Helmholtz
module cmod_GradOp
  
  
    use iso_c_binding
  
  
contains

    !>  \brief computes \f$ \mathrm{grad = \nabla\phi } \f$
    subroutine OP_grad( &
        dir,            &
        N,              &
        bL,bU,          &
        gL,gU,          &
        SS,NN,          &
        c,              &
        phi,            &
        grad ) bind (c,name='OP_grad')
  

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

        real(c_double), intent(out) :: grad(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  

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
        if( 3==dir )then

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

    end subroutine OP_grad


    !>  \brief sets non block boundary conditions to zero
    subroutine OP_SetBCZero(    &
        N,                      &
        bL,bU,                  &
        BC_L,BC_U,              &
        SB,NB,                  &
        grad ) bind (c,name='OP_SetBCZero')


        implicit none


        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: BC_L(3)
        integer(c_int), intent(in)  :: BC_U(3)

        integer(c_int), intent(in)  :: SB(3)
        integer(c_int), intent(in)  :: NB(3)

        real(c_double), intent(inout) :: grad(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        !--- Boundary conditions ------------------------------------------------------------------------------------
        if (BC_L(1) > 0) grad(SB(1)      ,SB(2):NB(2),SB(3):NB(3)) = 0.
        if (BC_U(1) > 0) grad(      NB(1),SB(2):NB(2),SB(3):NB(3)) = 0.
        if (BC_L(2) > 0) grad(SB(1):NB(1),SB(2)      ,SB(3):NB(3)) = 0.
        if (BC_U(2) > 0) grad(SB(1):NB(1),      NB(2),SB(3):NB(3)) = 0.
        if (BC_L(3) > 0) grad(SB(1):NB(1),SB(2):NB(2),SB(3)      ) = 0.
        if (BC_U(3) > 0) grad(SB(1):NB(1),SB(2):NB(2),      NB(3)) = 0.

  
    end subroutine OP_SetBCZero
  
end module cmod_GradOp
