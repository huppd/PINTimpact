!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing convection
module cmod_ConvectionDiffusionOp
  
  
    use iso_c_binding
  
  
contains


    subroutine OP_convectionDiffusion(  &
        dimens,                         &
        N,                              &
        bL,bU,                          &
        nL,nU,                          &
        SS,NN,                          &
        c1D,c2D,c3D,                    &
        c1U,c2U,c3U,                    &
        c11,c22,c33,                    &
        phiU,phiV,phiW,                 &
        phi,                            &
        nlu,                            &
        mulI,                           &
        mulL,                           &
        mul ) bind (c,name='OP_convectionDiffusion')

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

        real(c_double), intent(in)  :: c11(bL(1):bU(1),0:N(1))
        real(c_double), intent(in)  :: c22(bL(2):bU(2),0:N(2))
        real(c_double), intent(in)  :: c33(bL(3):bU(3),0:N(3))


        real(c_double), intent(in   ) ::  phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in   ) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(inout) ::  nlu(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in)  :: mul
        real(c_double), intent(in)  :: mulI
        real(c_double), intent(in)  :: mulL

        real                   ::  ddU, ddV, ddW, dd1

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

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = bL(3), bU(3)
                            dd1 = dd1 + c33(kk,k)*phi(i,j,k+kk)
                        end do

                        nlu(i,j,k) = mulI*nlu(i,j,k) + mul*(phiU(i,j,k)*ddU + phiV(i,j,k)*ddV + mul*phiW(i,j,k)*ddW ) - mulL*dd1

                    else

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do

                        nlu(i,j,k) = mulI*nlu(i,j,k) + mul*(phiU(i,j,k)*ddU + phiV(i,j,k)*ddV  ) - mulL*dd1

                    end if

                end do
            end do
        end do


    end subroutine OP_convectionDiffusion



    subroutine OP_convectionDiffusionSOR(   &
        dimens,                             &
        N,                                  &
        bL,bU,                              &
        nL,nU,                              &
        SS,NN,                              &
        dir,                                &
        c1D,c2D,c3D,                        &
        c1U,c2U,c3U,                        &
        c11,c22,c33,                        &
        phiU,phiV,phiW,                     &
        b,                                  &
        phi,                                &
        mulI,                               &
        mulL,                               &
        mul,                                &
        om) bind (c,name='OP_convectionDiffusionSOR')

        implicit none

        integer(c_int), intent(in)  :: dimens

        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: nL(3)
        integer(c_int), intent(in)  :: nU(3)

        integer(c_int), intent(in)  :: SS(3)
        integer(c_int), intent(in)  :: NN(3)

        integer(c_int), intent(in)  :: dir(3)

        real(c_double), intent(in)  :: c1D(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2D(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3D(nL(3):nU(3),0:N(3))

        real(c_double), intent(in)  :: c1U(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2U(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3U(nL(3):nU(3),0:N(3))

        real(c_double), intent(in)  :: c11(bL(1):bU(1),0:N(1))
        real(c_double), intent(in)  :: c22(bL(2):bU(2),0:N(2))
        real(c_double), intent(in)  :: c33(bL(3):bU(3),0:N(3))


        real(c_double), intent(in   ) ::  phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in   ) ::  b(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        real(c_double), intent(in)  :: mul
        real(c_double), intent(in)  :: mulI
        real(c_double), intent(in)  :: mulL
        real(c_double), intent(in)  :: om


        real                   ::  ddU, ddV, ddW, dd1, dd

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        integer(c_int)         :: SSS(3)
        integer(c_int)         :: NNN(3)

        do i = 1,3
            if( dir(i)>0 )then
                SSS(i) = SS(i)
                NNN(i) = NN(i)
            else
                SSS(i) = NN(i)
                NNN(i) = SS(i)
            end if
        end do



        do k = SSS(3), NNN(3), dir(3)
            do j = SSS(2), NNN(2), dir(2)
                do i = SSS(1), NNN(1), dir(1)
                    dd = 0
                    !===========================================================================================================
                    !=== u*d /dx ===============================================================================================
                    !===========================================================================================================
                    if( phiU(i,j,k) >= 0. ) then
                        ddU = c1U(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            if( ii/=0 ) ddU = ddU + c1U(ii,i)*phi(i+ii,j,k)
                        end do
                        dd = dd + mul*c1U(0,i)*phiU(i,j,k)
                    else
                        ddU = c1D(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            if( ii/=0 ) ddU = ddU + c1D(ii,i)*phi(i+ii,j,k)
                        end do
                        dd = dd + mul*c1D(0,i)*phiU(i,j,k)
                    end if

                    !===========================================================================================================
                    !=== v*d /dy ===============================================================================================
                    !===========================================================================================================
                    if( phiV(i,j,k) >= 0. ) then
                        ddV = c2U(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            if( jj/=0 )ddV = ddV + c2U(jj,j)*phi(i,j+jj,k)
                        end do
                        dd = dd + mul*c2U(0,j)*phiV(i,j,k)
                    else
                        ddV = c2D(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            if( jj/=0 ) ddV = ddV + c2D(jj,j)*phi(i,j+jj,k)
                        end do
                        dd = dd + mul*c2D(0,j)*phiV(i,j,k)
                    end if

                    if (dimens == 3) then

                        !===========================================================================================================
                        !=== w*d /dz ===============================================================================================
                        !===========================================================================================================
                        if( phiW(i,j,k) >= 0. ) then
                            ddW = c3U(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                if( kk/=0 )ddW = ddW + c3U(kk,k)*phi(i,j,k+kk)
                            end do
                            dd = dd + mul*c3U(0,k)*phiW(i,j,k)
                        else
                            ddW = c3D(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                if( kk/=0 ) ddW = ddW + c3D(kk,k)*phi(i,j,k+kk)
                            end do
                            dd = dd + mul*c3D(0,k)*phiW(i,j,k)
                        end if

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            if( ii/=0 ) dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            if( jj/=0 ) dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = bL(3), bU(3)
                            if( kk/=0 )dd1 = dd1 + c33(kk,k)*phi(i,j,k+kk)
                        end do
                        dd = dd - mulL*( c11(0,i) + c22(0,j) +c33(0,k) ) + mulI

                        phi(i,j,k) = (1-om)*phi(i,j,k) + om/dd*( b(i,j,k) - mul*(phiU(i,j,k)*ddU + phiV(i,j,k)*ddV + phiW(i,j,k)*ddW ) + mulL*dd1 )

                    else

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            if( ii/=0 )dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            if( jj/=0 )dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do

                        dd = dd - mulL*( c11(0,i) + c22(0,j) ) + mulI

                        !! check
                        phi(i,j,k) = (1-om)*phi(i,j,k) + om/dd*( b(i,j,k) - mul*phiU(i,j,k)*ddU - mul*phiV(i,j,k)*ddV   + mulL*dd1 )

                    end if

                end do
            end do
        end do


    end subroutine OP_convectionDiffusionSOR


       subroutine OP_convectionDiffusionJSmoother(  &
        dimens,                                     &
        N,                                          &
        bL,bU,                                      &
        nL,nU,                                      &
        SS,NN,                                      &
!        dir,                                        &
        c1D,c2D,c3D,                                &
        c1U,c2U,c3U,                                &
        c11,c22,c33,                                &
        phiU,phiV,phiW,                             &
        b,                                          &
        phi,                                        &
        phiout,                                     &
        mulI,                                       &
        mulL,                                       &
        mul,                                        &
        om) bind (c,name='OP_convectionDiffusionJSmoother')

        implicit none

        integer(c_int), intent(in)  :: dimens

        integer(c_int), intent(in)  :: N(3)

        integer(c_int), intent(in)  :: bL(3)
        integer(c_int), intent(in)  :: bU(3)

        integer(c_int), intent(in)  :: nL(3)
        integer(c_int), intent(in)  :: nU(3)

        integer(c_int), intent(in)  :: SS(3)
        integer(c_int), intent(in)  :: NN(3)

!        integer(c_int), intent(in)  :: dir(3)

        real(c_double), intent(in)  :: c1D(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2D(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3D(nL(3):nU(3),0:N(3))

        real(c_double), intent(in)  :: c1U(nL(1):nU(1),0:N(1))
        real(c_double), intent(in)  :: c2U(nL(2):nU(2),0:N(2))
        real(c_double), intent(in)  :: c3U(nL(3):nU(3),0:N(3))

        real(c_double), intent(in)  :: c11(bL(1):bU(1),0:N(1))
        real(c_double), intent(in)  :: c22(bL(2):bU(2),0:N(2))
        real(c_double), intent(in)  :: c33(bL(3):bU(3),0:N(3))


        real(c_double), intent(in   ) ::  phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(in   ) ::  phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in   ) ::  b(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in   ) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(  out) ::  phiout(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        real(c_double), intent(in)  :: mul
        real(c_double), intent(in)  :: mulI
        real(c_double), intent(in)  :: mulL
        real(c_double), intent(in)  :: om

        real                   ::  ddU, ddV, ddW, dd1, dd

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk




        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                do i = SS(1), NN(1)
                    dd = 0
                    !===========================================================================================================
                    !=== u*d /dx ===============================================================================================
                    !===========================================================================================================
                    if( phiU(i,j,k) >= 0. ) then
                        ddU = c1U(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            if( ii/=0 ) ddU = ddU + c1U(ii,i)*phi(i+ii,j,k)
                        end do
                        dd = dd + mul*c1U(0,i)*phiU(i,j,k)
                    else
                        ddU = c1D(nL(1),i)*phi(i+nL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = nL(1)+1, nU(1)
                            if( ii/=0 ) ddU = ddU + c1D(ii,i)*phi(i+ii,j,k)
                        end do
                        dd = dd + mul*c1D(0,i)*phiU(i,j,k)
                    end if

                    !===========================================================================================================
                    !=== v*d /dy ===============================================================================================
                    !===========================================================================================================
                    if( phiV(i,j,k) >= 0. ) then
                        ddV = c2U(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            if( jj/=0 )ddV = ddV + c2U(jj,j)*phi(i,j+jj,k)
                        end do
                        dd = dd + mul*c2U(0,j)*phiV(i,j,k)
                    else
                        ddV = c2D(nL(2),j)*phi(i,j+nL(2),k)
                        !pgi$ unroll = n:8
                        do jj = nL(2)+1, nU(2)
                            if( jj/=0 ) ddV = ddV + c2D(jj,j)*phi(i,j+jj,k)
                        end do
                        dd = dd + mul*c2D(0,j)*phiV(i,j,k)
                    end if

                    if (dimens == 3) then

                        !===========================================================================================================
                        !=== w*d /dz ===============================================================================================
                        !===========================================================================================================
                        if( phiW(i,j,k) >= 0. ) then
                            ddW = c3U(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                if( kk/=0 )ddW = ddW + c3U(kk,k)*phi(i,j,k+kk)
                            end do
                            dd = dd + mul*c3U(0,k)*phiW(i,j,k)
                        else
                            ddW = c3D(nL(3),k)*phi(i,j,k+nL(3))
                            !pgi$ unroll = n:8
                            do kk = nL(3)+1, nU(3)
                                if( kk/=0 ) ddW = ddW + c3D(kk,k)*phi(i,j,k+kk)
                            end do
                            dd = dd + mul*c3D(0,k)*phiW(i,j,k)
                        end if

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            if( ii/=0 ) dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            if( jj/=0 ) dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = bL(3), bU(3)
                            if( kk/=0 )dd1 = dd1 + c33(kk,k)*phi(i,j,k+kk)
                        end do
                        dd = dd - mulL*( c11(0,i) + c22(0,j) +c33(0,k) ) + mulI

                        phiout(i,j,k) = (1-om)*phi(i,j,k) + om/dd*( b(i,j,k) - mul*(phiU(i,j,k)*ddU + phiV(i,j,k)*ddV + phiW(i,j,k)*ddW ) + mulL*dd1 )

                    else

                        !===========================================================================================================
                        !=== d^2 /dx^2 + d^2 /dy^2 ===============================================================================================
                        !===========================================================================================================
                        dd1 = c11(bL(1),i)*phi(i+bL(1),j,k)
                        !pgi$ unroll = n:8
                        do ii = bL(1)+1, bU(1)
                            if( ii/=0 )dd1 = dd1 + c11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = bL(2), bU(2)
                            if( jj/=0 )dd1 = dd1 + c22(jj,j)*phi(i,j+jj,k)
                        end do

                        dd = dd - mulL*( c11(0,i) + c22(0,j) ) + mulI

                        !! check
                        phiout(i,j,k) = (1-om)*phi(i,j,k) + om/dd*( b(i,j,k) - mul*phiU(i,j,k)*ddU - mul*phiV(i,j,k)*ddV   + mulL*dd1 )

                    end if

                end do
            end do
        end do


    end subroutine OP_convectionDiffusionJSmoother
  
end module cmod_ConvectionDiffusionOp
