!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing Helmholtz
module cmod_InterpolateV2SOp
  
    use iso_c_binding
  
    implicit none

contains

    !>  \brief interpolates vel 2 scalar
    subroutine OP_interpolateV2S(   &
        m,                          &
        N,                          &
        bL,bU,                      &
        dL,dU,                      &
        SS,NN,                      &
        c,                          &
        phi,                        &
        inter ) bind (c,name='OP_interpolateV2S')
  
        implicit none
  
        integer(c_int), intent(in)  :: m

        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: dL
        integer(c_int), intent(in)  :: dU

        integer(c_int), intent(in)  :: SS(3)
        integer(c_int), intent(in)  :: NN(3)

        real(c_double), intent(in)  :: c(dL:dU,0:N(m))

        real(c_double), intent(in)  :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(out) :: inter (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  

        !===========================================================================================================
        if( m == 1 )then

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        inter(i,j,k) = c(dL,i)*phi(i+dL,j,k)
                        do ii = dL+1, dU
                            inter(i,j,k) = inter(i,j,k) + c(ii,i)*phi(i+ii,j,k)
                        end do
                    end do
                end do
            end do

        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 2) then

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        inter(i,j,k) = c(dL,j)*phi(i,j+dL,k)
                        do jj = dL+1, dU
                            inter(i,j,k) = inter(i,j,k) + c(jj,j)*phi(i,j+jj,k)
                        end do
                    end do
                end do
            end do

        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 3) then
            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        inter(i,j,k) = c(dL,k)*phi(i,j,k+dL)
                        do kk = dL+1, dU
                            inter(i,j,k) = inter(i,j,k) + c(kk,k)*phi(i,j,k+kk)
                        end do
                    end do
                end do
            end do
        end if
        !===========================================================================================================

  
    end subroutine OP_interpolateV2S

  
end module cmod_InterpolateV2SOp
