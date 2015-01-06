!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing convection
module cmod_ConvectionOp
  
    use iso_c_binding
  
contains


    subroutine OP_convection(   &
        dimens,                 &
        N,                      &
        bL,bU,                  &
        nL,nU,                  &
        SS,NN,                  &
        c1D,c2D,c3D,            &
        c1U,c2U,c3U,            &
        phiU,phiV,phiW,         &
        phi,                    &
        nlu,                    &
        mul ) bind (c,name='OP_convection')

        implicit none

        integer(c_int), intent(in)  :: dimens

        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: nL(3)
        integer(c_int), intent(in)  :: nU(3)

        integer(c_int), intent(in)  :: SS(3)
        integer(c_int), intent(in)  :: NN(3)

        real(c_double), intent(in)  :: c1D(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2D(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3D(nL(3):nU(3),0:N(3))

        real(c_double), intent(in)  :: c1U(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2U(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3U(nL(3):nU(3),0:N(3))

        real(c_double), intent(in   ) ::  phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in   ) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(inout) ::  nlu(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in)    :: mul

        real                   ::  ddU, ddV, ddW

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk



        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                do i = SS(1), NN(1)
                    !===========================================================================================================
                    !=== u*d /dx ===============================================================================================
                    !===========================================================================================================
                    if( phiU(i,j,k) >= 0. ) then
                        ddU = c1U(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            ddU = ddU + c1U(ii,i)*phi(i+ii,j,k)
                        end do
                    else
                        ddU = c1D(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            ddU = ddU + c1D(ii,i)*phi(i+ii,j,k)
                        end do
                    end if

                    !===========================================================================================================
                    !=== v*d /dy ===============================================================================================
                    !===========================================================================================================
                    if( phiV(i,j,k) >= 0. ) then
                        ddV = c2U(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            ddV = ddV + c2U(jj,j)*phi(i,j+jj,k)
                        end do
                    else
                        ddV = c2D(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            ddV = ddV + c2D(jj,j)*phi(i,j+jj,k)
                        end do
                    end if

                    if (dimens == 3) then

                        !===========================================================================================================
                        !=== w*d /dz ===============================================================================================
                        !===========================================================================================================
                        if( phiW(i,j,k) >= 0. ) then
                            ddW = c3U(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                ddW = ddW + c3U(kk,k)*phi(i,j,k+kk)
                            end do
                        else
                            ddW = c3D(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                ddW = ddW + c3D(kk,k)*phi(i,j,k+kk)
                            end do
                        end if

                        nlu(i,j,k) = nlu(i,j,k) + mul*(phiU(i,j,k)*ddU + phiV(i,j,k)*ddV + mul*phiW(i,j,k)*ddW )

                    else

                        nlu(i,j,k) = nlu(i,j,k) + mul*( phiU(i,j,k)*ddU + phiV(i,j,k)*ddV  )

                    end if

                end do
            end do
        end do

    end subroutine OP_convection
  
end module cmod_ConvectionOp
