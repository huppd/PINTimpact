!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing DtHelmholtz
module cmod_DtHelmholtzOp
  
  
    use iso_c_binding
  
  
contains


    !>  \brief computes \f$ \mathrm{lap_m = mulI phiI_m - mulL \Delta phiL_m} \f$
    !!
    !! \param[in] dimens dimension 2 or 3
    !! \param[in] N local amount of grid points
    !! \param[in] m dimension from one to three
    !! \param[in] mulI factor which coresponds to the factor of the identity part
    !! \param[in] mulL factor which coresponds to the factor of the laplace part
    !! \param[in] phiI
    !! \param[in] phiL
    !! \param[out] Lap
    subroutine OP_DtHelmholtz(  &
        dimens,                 &
        N,                      &
        bL,bU,                  &
        SS,NN,                  &
        c11,c22,c33,            &
        mulI,                   &
        mulL,                   &
        phiI,                   &
        phiL,                   &
        lap ) bind (c,name='OP_DtHelmholtz')
  
        implicit none
  
        integer(c_int), intent(in)  :: dimens

        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: SS(3)
        integer(c_int), intent(in)  :: NN(3)

        real(c_double), intent(in)  :: c11(bL(1):bU(1),0:N(1))
        real(c_double), intent(in)  :: c22(bL(2):bU(2),0:N(2))
        real(c_double), intent(in)  :: c33(bL(3):bU(3),0:N(3))

        real(c_double), intent(in)  :: mulI
        real(c_double), intent(in)  :: mulL

        real(c_double), intent(in)   :: phiI (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in)   :: phiL (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(inout):: lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        real                   ::  dd1
  

        !===========================================================================================================
        if( 3==dimens) then

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        dd1 = c11(bL(1),i)*phiL(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phiL(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phiL(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = bL(3), bU(3)
                            dd1 = dd1 + c33(kk,k)*phiL(i,j,k+kk)
                        end do
                        Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                    end do
                end do
            end do

        else

            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        dd1 = c11(bL(1),i)*phiL(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phiL(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phiL(i,j+jj,k)
                        end do
                        Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                    end do
                end do
            end do

        end if
        !===========================================================================================================

  
    end subroutine OP_DtHelmholtz

  
end module cmod_DtHelmholtzOp