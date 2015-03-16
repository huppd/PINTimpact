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
    subroutine OP_div(  &
        dimens,         &
        N,              &
        bL,bU,          &
        dL,dU,          &
        SS,NN,          &
        c1,c2,c3,       &
        phiU,           &
        phiV,           &
        phiW,           &
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

        real(c_double), intent(in)  :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in)  :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in)  :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(out) :: div (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        real                   ::  dd1
  

        !===========================================================================================================
        if( 3==dimens) then

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        div(i,j,k) = c1(dL(1),i)*phiU(i+dL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = dL(1)+1, dU(1)
                            div(i,j,k) = div(i,j,k) + c1(ii,i)*phiU(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = dL(2), dU(2)
                            div(i,j,k) = div(i,j,k) + c2(jj,j)*phiV(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = dL(3), dU(3)
                            div(i,j,k) = div(i,j,k) + c3(kk,k)*phiV(i,j,k+kk)
                        end do
                    end do
                end do
            end do

        else

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        div(i,j,k) = c1(dL(1),i)*phiU(i+dL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = dL(1)+1, dU(1)
                            div(i,j,k) = div(i,j,k) + c1(ii,i)*phiU(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = dL(2), dU(2)
                            div(i,j,k) = div(i,j,k) + c2(jj,j)*phiV(i,j+jj,k)
                        end do
                    end do
                end do
            end do

        end if
    !===========================================================================================================

  
    end subroutine OP_div

  
end module cmod_DivOp
