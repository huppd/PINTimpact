!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing convection
module cmod_ConvectionDiffusionOp

  use iso_c_binding

  implicit none

contains


  !> \brief convection diffusion 
  !! \f[ nlu = mul nlu + mulI phi + mulC ( U\cdot\nabla) phi - mulL \Delta phi \f]
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
      mul,                            &
      mulI,                           &
      mulC,                           &
      mulL ) bind (c,name='OP_convectionDiffusion')

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
    real(c_double), intent(in)  :: mulC
    real(c_double), intent(in)  :: mulL

    real(c_double)              ::  ddU, ddV, ddW, dd1

    integer(c_int)              ::  i, ii
    integer(c_int)              ::  j, jj
    integer(c_int)              ::  k, kk


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
            !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 =====================================================================
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

            nlu(i,j,k) = mul*nlu(i,j,k) + mulI*phi(i,j,k) + mulC*( phiU(i,j,k)*ddU + phiV(i,j,k)*ddV + phiW(i,j,k)*ddW ) - mulL*dd1

          else

            !===========================================================================================================
            !=== d^2 /dx^2 + d^2 /dy^2 =================================================================================
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

            nlu(i,j,k) = mul*nlu(i,j,k) + mulI*phi(i,j,k) + mulC*( phiU(i,j,k)*ddU + phiV(i,j,k)*ddV  ) - mulL*dd1

          end if

        end do
      end do
    end do


  end subroutine OP_convectionDiffusion



  !> \brief convection diffusion SOR
  subroutine OP_convectionDiffusionSOR(   &
      dimens,                             &
      N,                                  &
      bL,bU,                              &
      nL,nU,                              &
      SS,NN,                              &
      dir,                                &
      loopOrder,                          &
      c1D,c2D,c3D,                        &
      c1U,c2U,c3U,                        &
      c11,c22,c33,                        &
      phiU,phiV,phiW,                     &
      b,                                  &
      phi,                                &
      mulI,                               &
      mulC,                               &
      mulL,                               &
      om ) bind (c,name='OP_convectionDiffusionSOR')

    implicit none

    integer(c_int),   intent(in)  :: dimens

    integer(c_int),   intent(in)    :: N(3)

    integer(c_int),   intent(in)    :: bL(3)
    integer(c_int),   intent(in)    :: bU(3)

    integer(c_int),   intent(in)    :: nL(3)
    integer(c_int),   intent(in)    :: nU(3)

    integer(c_int),   intent(in)    :: SS(3)
    integer(c_int),   intent(in)    :: NN(3)

    integer(c_short), intent(in)    :: dir(3)
    integer(c_short), intent(in)    :: loopOrder(3)

    real(c_double),   intent(in)    :: c1D(nL(1):nU(1),0:N(1))
    real(c_double),   intent(in)    :: c2D(nL(2):nU(2),0:N(2))
    real(c_double),   intent(in)    :: c3D(nL(3):nU(3),0:N(3))

    real(c_double),   intent(in)    :: c1U(nL(1):nU(1),0:N(1))
    real(c_double),   intent(in)    :: c2U(nL(2):nU(2),0:N(2))
    real(c_double),   intent(in)    :: c3U(nL(3):nU(3),0:N(3))

    real(c_double),   intent(in)    :: c11(bL(1):bU(1),0:N(1))
    real(c_double),   intent(in)    :: c22(bL(2):bU(2),0:N(2))
    real(c_double),   intent(in)    :: c33(bL(3):bU(3),0:N(3))


    real(c_double),   intent(in)    ::  phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),   intent(in)    ::  phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),   intent(in)    ::  phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double),   intent(in)    ::  b(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double),   intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    real(c_double), intent(in)  :: mulI
    real(c_double), intent(in)  :: mulC
    real(c_double), intent(in)  :: mulL
    real(c_double), intent(in)  :: om


    real(c_double) ::  ddU, ddV, ddW, dd1, dd

    integer(c_int) ::  i, ii
    integer(c_int) ::  j, jj
    integer(c_int) ::  k, kk

    integer(c_int) ::  ind(3)

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

    do k = SSS(loopOrder(3)), NNN(loopOrder(3)), dir(loopOrder(3))
      do j = SSS(loopOrder(2)), NNN(loopOrder(2)), dir(loopOrder(2))
        do i = SSS(loopOrder(1)), NNN(loopOrder(1)), dir(loopOrder(1))
          !        do k = SSS(3), NNN(3), dir(3)
          !            do j = SSS(2), NNN(2), dir(2)
          !                do i = SSS(1), NNN(1), dir(1)

          ind(loopOrder(1)) = i
          ind(loopOrder(2)) = j
          ind(loopOrder(3)) = k

          dd = 0
          !===========================================================================================================
          !=== u*d /dx ===============================================================================================
          !===========================================================================================================
          if( phiU( ind(1),ind(2),ind(3) ) >= 0. ) then
            ddU = c1U( nL(1),ind(1) )*phi( ind(1)+nL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = nL(1)+1, nU(1)
              if( ii/=0 ) ddU = ddU + c1U( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            dd = dd + mulC*c1U( 0,ind(1) )*phiU( ind(1),ind(2),ind(3) )
          else
            ddU = c1D( nL(1),ind(1) )*phi( ind(1)+nL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = nL(1)+1, nU(1)
              if( ii/=0 ) ddU = ddU + c1D( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            dd = dd + mulC*c1D( 0,ind(1) )*phiU( ind(1),ind(2),ind(3) )
          end if
          !===========================================================================================================
          !=== v*d /dy ===============================================================================================
          !===========================================================================================================
          if( phiV( ind(1),ind(2),ind(3) ) >= 0. ) then
            ddV = c2U( nL(2),ind(2) )*phi( ind(1),ind(2)+nL(2),ind(3) )
            !pgi$ unroll = n:8
            do jj = nL(2)+1, nU(2)
              if( jj/=0 ) ddV = ddV + c2U( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do
            dd = dd + mulC*c2U( 0,ind(2) )*phiV( ind(1),ind(2),ind(3) )
          else
            ddV = c2D(nL(2),ind(2))*phi( ind(1),ind(2)+nL(2),ind(3) )
            !pgi$ unroll = n:8
            do jj = nL(2)+1, nU(2)
              if( jj/=0 ) ddV = ddV + c2D( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do
            dd = dd + mulC*c2D( 0,ind(2) )*phiV( ind(1),ind(2),ind(3) )
          end if


          if (dimens == 3) then

            !===========================================================================================================
            !=== w*d /dz ===============================================================================================
            !===========================================================================================================
            if( phiW( ind(1),ind(2),ind(3) ) >= 0. ) then
              ddW = c3U(nL(3),ind(3))*phi(ind(1),ind(2),ind(3)+nL(3))
              !pgi$ unroll = n:8
              do kk = nL(3)+1, nU(3)
                if( kk/=0 )ddW = ddW + c3U(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
              end do
              dd = dd + mulC*c3U(0,ind(3))*phiW(ind(1),ind(2),ind(3))
            else
              ddW = c3D(nL(3),ind(3))*phi(ind(1),ind(2),ind(3)+nL(3))
              !pgi$ unroll = n:8
              do kk = nL(3)+1, nU(3)
                if( kk/=0 ) ddW = ddW + c3D(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
              end do
              dd = dd + mulC*c3D(0,ind(3))*phiW(ind(1),ind(2),ind(3))
            end if

            !===========================================================================================================
            !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 ===============================================================================================
            !===========================================================================================================
            dd1 = c11(bL(1),ind(1))*phi(ind(1)+bL(1),ind(2),ind(3))
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
              if( ii/=0 ) dd1 = dd1 + c11(ii,ind(1))*phi(ind(1)+ii,ind(2),ind(3))
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
              if( jj/=0 ) dd1 = dd1 + c22(jj,ind(2))*phi(i,ind(2)+jj,ind(3))
            end do
            !pgi$ unroll = n:8
            do kk = bL(3), bU(3)
              if( kk/=0 ) dd1 = dd1 + c33(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
            end do
            dd = dd - mulL*( c11(0,ind(1)) + c22(0,ind(2)) +c33(0,ind(3)) ) + mulI

            phi(ind(1),ind(2),ind(3)) = (1-om)*phi(ind(1),ind(2),ind(3)) + om/dd*( b(ind(1),ind(2),ind(3)) - mulC*(phiU(ind(1),ind(2),ind(3))*ddU + phiV(ind(1),ind(2),ind(3))*ddV + phiW(ind(1),ind(2),ind(3))*ddW ) + mulL*dd1 )

          else

            !===========================================================================================================
            !=== d^2 /dx^2 + d^2 /dy^2 ===============================================================================================
            !===========================================================================================================
            dd1 = c11( bL(1),ind(1) )*phi( ind(1)+bL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
              if( ii/=0 ) dd1 = dd1 + c11( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
              if( jj/=0 ) dd1 = dd1 + c22( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do

            dd = dd - mulL*( c11( 0,ind(1) ) + c22( 0,ind(2) ) ) + mulI

            phi( ind(1),ind(2),ind(3) ) =                       &
              (1-om)*phi(ind(1),ind(2),ind(3))               &
              + om/dd*( b(ind(1),ind(2),ind(3))               &
              - mulC*phiU(ind(1),ind(2),ind(3))*ddU    &
              - mulC*phiV(ind(1),ind(2),ind(3))*ddV    &
              + mulL*dd1 )

          end if

        end do
      end do
    end do


  end subroutine OP_convectionDiffusionSOR



  subroutine OP_convectionDiffusionJSmoother( &
      dimens,                                 &
      N,                                      &
      bL,bU,                                  &
      nL,nU,                                  &
      SS,NN,                                  &
      c1D,c2D,c3D,                            &
      c1U,c2U,c3U,                            &
      c11,c22,c33,                            &
      phiU,phiV,phiW,                         &
      b,                                      &
      phi,                                    &
      phiout,                                 &
      mulI,                                   &
      mulC,                                   &
      mulL,                                   &
      om ) bind (c,name='OP_convectionDiffusionJSmoother')

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


    real(c_double), intent(in)  :: mulI
    real(c_double), intent(in)  :: mulC
    real(c_double), intent(in)  :: mulL
    real(c_double), intent(in)  :: om

    real(c_double) :: ddU, ddV, ddW, dd1, dd

    integer(c_int) :: ind(1:3)
    integer(c_int) :: i, ii
    integer(c_int) :: j, jj
    integer(c_int) :: k, kk

    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        do i = SS(1), NN(1)

          ind(1) = i
          ind(2) = j
          ind(3) = k

          dd = 0
          !===========================================================================================================
          !=== u*d /dx ===============================================================================================
          !===========================================================================================================
          if( phiU( ind(1),ind(2),ind(3) ) >= 0. ) then
            ddU = c1U( nL(1),ind(1) )*phi( ind(1)+nL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = nL(1)+1, nU(1)
              if( ii/=0 ) ddU = ddU + c1U( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            dd = dd + mulC*c1U( 0,ind(1) )*phiU( ind(1),ind(2),ind(3) )
          else
            ddU = c1D( nL(1),ind(1) )*phi( ind(1)+nL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = nL(1)+1, nU(1)
              if( ii/=0 ) ddU = ddU + c1D( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            dd = dd + mulC*c1D( 0,ind(1) )*phiU( ind(1),ind(2),ind(3) )
          end if
          !===========================================================================================================
          !=== v*d /dy ===============================================================================================
          !===========================================================================================================
          if( phiV( ind(1),ind(2),ind(3) ) >= 0. ) then
            ddV = c2U( nL(2),ind(2) )*phi( ind(1),ind(2)+nL(2),ind(3) )
            !pgi$ unroll = n:8
            do jj = nL(2)+1, nU(2)
              if( jj/=0 ) ddV = ddV + c2U( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do
            dd = dd + mulC*c2U( 0,ind(2) )*phiV( ind(1),ind(2),ind(3) )
          else
            ddV = c2D(nL(2),ind(2))*phi( ind(1),ind(2)+nL(2),ind(3) )
            !pgi$ unroll = n:8
            do jj = nL(2)+1, nU(2)
              if( jj/=0 ) ddV = ddV + c2D( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do
            dd = dd + mulC*c2D( 0,ind(2) )*phiV( ind(1),ind(2),ind(3) )
          end if


          if (dimens == 3) then

            !===========================================================================================================
            !=== w*d /dz ===============================================================================================
            !===========================================================================================================
            if( phiW( ind(1),ind(2),ind(3) ) >= 0. ) then
              ddW = c3U(nL(3),ind(3))*phi(ind(1),ind(2),ind(3)+nL(3))
              !pgi$ unroll = n:8
              do kk = nL(3)+1, nU(3)
                if( kk/=0 )ddW = ddW + c3U(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
              end do
              dd = dd + mulC*c3U(0,ind(3))*phiW(ind(1),ind(2),ind(3))
            else
              ddW = c3D(nL(3),ind(3))*phi(ind(1),ind(2),ind(3)+nL(3))
              !pgi$ unroll = n:8
              do kk = nL(3)+1, nU(3)
                if( kk/=0 ) ddW = ddW + c3D(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
              end do
              dd = dd + mulC*c3D(0,ind(3))*phiW(ind(1),ind(2),ind(3))
            end if

            !===========================================================================================================
            !=== d^2 /dx^2 + d^2 /dy^2 + d^2 /dz^2 ===============================================================================================
            !===========================================================================================================
            dd1 = c11(bL(1),ind(1))*phi(ind(1)+bL(1),ind(2),ind(3))
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
              if( ii/=0 ) dd1 = dd1 + c11(ii,ind(1))*phi(ind(1)+ii,ind(2),ind(3))
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
              if( jj/=0 ) dd1 = dd1 + c22(jj,ind(2))*phi(i,ind(2)+jj,ind(3))
            end do
            !pgi$ unroll = n:8
            do kk = bL(3), bU(3)
              if( kk/=0 ) dd1 = dd1 + c33(kk,ind(3))*phi(ind(1),ind(2),ind(3)+kk)
            end do
            dd = dd - mulL*( c11(0,ind(1)) + c22(0,ind(2)) +c33(0,ind(3)) ) + mulI

            phiout(ind(1),ind(2),ind(3)) = (1-om)*phi(ind(1),ind(2),ind(3)) + om/dd*( b(ind(1),ind(2),ind(3)) - mulC*(phiU(ind(1),ind(2),ind(3))*ddU + phiV(ind(1),ind(2),ind(3))*ddV + phiW(ind(1),ind(2),ind(3))*ddW ) + mulL*dd1 )

          else

            !===========================================================================================================
            !=== d^2 /dx^2 + d^2 /dy^2 ===============================================================================================
            !===========================================================================================================
            dd1 = c11( bL(1),ind(1) )*phi( ind(1)+bL(1),ind(2),ind(3) )
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
              if( ii/=0 ) dd1 = dd1 + c11( ii,ind(1) )*phi( ind(1)+ii,ind(2),ind(3) )
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
              if( jj/=0 ) dd1 = dd1 + c22( jj,ind(2) )*phi( ind(1),ind(2)+jj,ind(3) )
            end do

            dd = dd - mulL*( c11( 0,ind(1) ) + c22( 0,ind(2) ) ) + mulI

            phiout( ind(1),ind(2),ind(3) ) =                       &
              (1-om)*phi(ind(1),ind(2),ind(3))               &
              + om/dd*( b(ind(1),ind(2),ind(3))               &
              - mulC*phiU(ind(1),ind(2),ind(3))*ddU    &
              - mulC*phiV(ind(1),ind(2),ind(3))*ddV    &
              + mulL*dd1 )

          end if

        end do
      end do
    end do

  end subroutine OP_convectionDiffusionJSmoother

end module cmod_ConvectionDiffusionOp
