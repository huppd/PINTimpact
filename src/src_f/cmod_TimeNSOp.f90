!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



module cmod_TimeNSOp

  use iso_c_binding

  implicit none

contains


  !> \brief computes \f$ \begin{bmatrix} r\_vel \\ r\_p \end{bmatrix} = \begin{pmatrix} mulI\partial_t -mulL\Delta & \nabla \\ \nabla \cdot \end{pmatrix}\begin{bmatrix} r\_vel \\ r\_p \end{bmatrix}  \f$
  !!
  !! \param dimens denotes spatial dimension can be either two or three is
  !!        defined in \c Space::dim
  !! \param N is the spatial local grid size in the three spatial direction can be get from \c Space::nLoc()
  !! \param bl is the lower stencil width used for ghost layer and boundary
  !!        condtions can be get from Space::bl()
  !! \param bu is the upper stencil width used for ghost layer and boundary
  !!        condtions can be get from \c Space::bu()
  !! \param dl is the lower stencil width for the divergence stencil used for ghost layer and boundary
  !!        condtions can be get from \c Space::dl()
  !! \param bu is the upper stencil width for the divergence stencil used for ghost layer and boundary
  !!        condtions can be get from \c Space::du()
  !! \param gl is the lower stencil width for the gradiant stencil used for ghost layer and boundary
  !!        condtions can be get from \c Space::gl()
  !! \param gu is the upper stencil width for the gradiant stencil used for ghost layer and boundary
  !!        condtions can be get from \c Space::gu()
  !! \param SS starting index for pressure fields, can be get from \c
  !!        Space::sInd(EFieldType)
  !! \param NN end index for pressure fields, can be get from \c
  !!        Space::eInd(EFieldType)
  !! \param SU starting index for velocity fields in x-direction, can be get from \c
  !!        Space::sInd(EFieldType::U)
  !! \param NU end index for velocity fields in x-direction, can be get from \c
  !!        Space::eInd(EFieldType::U)
  !! \param SV starting index for velocity fields in y-direction, can be get from \c
  !!        Space::sInd(EFieldType::V)
  !! \param NV end index for velocity fields in y-direction, can be get from \c
  !!        Space::eInd(EFieldType::V)
  !! \param SW starting index for velocity fields in z-direction, can be get from \c
  !!        Space::sInd(EFieldType::W)
  !! \param NW end index for velocity fields in z-direction, can be get from \c
  !!        Space::eInd(EFieldType::W)
  !! \param c11p stencil coefficients for laplace operator in x-direction on
  !!        pressure coordinates, can be get from \c HelmholtzOp::getC(X,S)
  !! \param c22p stencil coefficients for laplace operator in y-direction on
  !!        pressure coordinates, can be get from \c HelmholtzOp::getC(Y,S)
  !! \param c33p stencil coefficients for laplace operator in z-direction on
  !!        pressure coordinates, can be get from \c HelmholtzOp::getC(Z,S)
  !! \param c11u stencil coefficients for laplace operator in x-direction on
  !!        velocity coordinates, can be get from \c HelmholtzOp::getC(X,U)
  !! \param c22v stencil coefficients for laplace operator in y-direction on
  !!        velocity coordinates, can be get from \c HelmholtzOp::getC(Y,V)
  !! \param c33w stencil coefficients for laplace operator in z-direction on
  !!        velocity coordinates, can be get from \c HelmholtzOp::getC(Z,W)
  !! \param cD1 stencil coefficients for divergence operator in x-direction, can be get from \c DivOp::getC(X)
  !! \param cD2 stencil coefficients for divergence operator in y-direction, can be get from \c DivOp::getC(Y)
  !! \param cD3 stencil coefficients for divergence operator in z-direction, can be get from \c DivOp::getC(Z)
  !! \param cG1 stencil coefficients for divergence operator in x-direction, can be get from \c GradOp::getC(X)
  !! \param cG2 stencil coefficients for divergence operator in y-direction, can be get from \c GradOp::getC(Y)
  !! \param cG3 stencil coefficients for divergence operator in z-direction, can be get from \c GradOp::getC(Z)
  !! \param mulI correspond to multiplicator of time derivative should be \f$ \frac{\alpha^2}{\mathrm{Re} \Delta t} \f$
  !! \param mulL correspond to multiplicator of laplace operator should be \f$ \frac{1}{\mathrm{Re}} \f$
  !! \param veln velocity field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param pn pressure field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param[out] r_vel residual velocity field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param[out] r_p residual pressure field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  subroutine OP_TimeNS(     &
      dimens,               &
      N,                    &
      bL,bU,                &
      cL,cU,                &
      dL,dU,                &
      gL,gU,                &
      SS,NN,                &
      SU,NU,                &
      SV,NV,                &
      SW,NW,                &
      c1uD,c2vD,c3wD,       &
      c1uU,c2vU,c3wU,       &
      c1pD,c2pD,c3pD,       &
      c1pU,c2pU,c3pU,       &
      c11p,c22p,c33p,       &
      c11u,c22v,c33w,       &
      cD1,                  &
      cD2,                  &
      cD3,                  &
      cG1,                  &
      cG2,                  &
      cG3,                  &
      mulI,                 &
      mulL,                 &
      windU,                &
      windV,                &
      windW,                &
      veln,                 &
      pn,                   &
      r_vel,                &
      r_p ) bind (c,name='OP_TimeNS')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(4)

    integer(c_int), intent(in)  :: bL(4)
    integer(c_int), intent(in)  :: bU(4)

    integer(c_int), intent(in)  :: dL(4)
    integer(c_int), intent(in)  :: dU(4)

    integer(c_int), intent(in)  :: cL(3)
    integer(c_int), intent(in)  :: cU(3)

    integer(c_int), intent(in)  :: gL(4)
    integer(c_int), intent(in)  :: gU(4)

    integer(c_int), intent(in)  :: SS(4)
    integer(c_int), intent(in)  :: NN(4)

    integer(c_int), intent(in)  :: SU(4)
    integer(c_int), intent(in)  :: NU(4)

    integer(c_int), intent(in)  :: SV(4)
    integer(c_int), intent(in)  :: NV(4)

    integer(c_int), intent(in)  :: SW(4)
    integer(c_int), intent(in)  :: NW(4)

    real(c_double), intent(in)  :: c1uD(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2vD(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3wD(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1uU(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2vU(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3wU(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1pD(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2pD(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3pD(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1pU(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2pU(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3pU(cL(3):cU(3),0:N(3))


    real(c_double), intent(in)  :: c11p( bL(1):bU(1), 0:N(1) )
    real(c_double), intent(in)  :: c22p( bL(2):bU(2), 0:N(2) )
    real(c_double), intent(in)  :: c33p( bL(3):bU(3), 0:N(3) )

    real(c_double), intent(in)  :: c11u( bL(1):bU(1), 0:N(1) )
    real(c_double), intent(in)  :: c22v( bL(2):bU(2), 0:N(2) )
    real(c_double), intent(in)  :: c33w( bL(3):bU(3), 0:N(3) )

    real(c_double), intent(in)  :: cD1( dL(1):dU(1), 0:N(1) )
    real(c_double), intent(in)  :: cD2( dL(2):dU(2), 0:N(2) )
    real(c_double), intent(in)  :: cD3( dL(3):dU(3), 0:N(3) )

    real(c_double), intent(in)  :: cG1( gL(1):gU(1), 0:N(1) )
    real(c_double), intent(in)  :: cG2( gL(2):gU(2), 0:N(2) )
    real(c_double), intent(in)  :: cG3( gL(3):gU(3), 0:N(3) )

    real(c_double), intent(in)  :: mulI
    real(c_double), intent(in)  :: mulL

    real(c_double), intent(in)  :: windU ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )
    real(c_double), intent(in)  :: windV ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )
    real(c_double), intent(in)  :: windW ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in)  :: veln ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in)  :: pn   ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(out) :: r_vel( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(out) :: r_p  ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )


    real(c_double)              :: dd1

    integer(c_int)              ::  i, ii
    integer(c_int)              ::  j, jj
    integer(c_int)              ::  k, kk
    integer(c_int)              ::  t 


    !===========================================================================================================

    do t = SS(4), N(4)
      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            !===========================================================================================================
            !=== computing velocity residual in x-direction ============================================================
            !===========================================================================================================
            if( SU(1)<=i .and. i<=NU(1) .and. SU(2)<=j .and. j<=NU(2) .and. SU(3)<=k .and. k<=NU(3) ) then
              !--- compute time derivative ----------------------------------------------------------------------------- 
              r_vel(i,j,k,1,t) =  mulI*( veln(i,j,k,1,t) - veln(i,j,k,1,t-1) )
              !--- compute convection -----------------------------------------------------------------------------------
              !--- u*d/dx ---
              if( windU(i,j,k,1,t) >= 0. ) then
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windU(i,j,k,1,t)*c1uU(ii,i)*veln(i+ii,j,k,1,t)
                end do
              else
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windU(i,j,k,1,t)*c1uD(ii,i)*veln(i+ii,j,k,1,t)
                end do
              end if
              !--- v*d/dy ---
              if( windV(i,j,k,1,t) >= 0. ) then
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windV(i,j,k,1,t)*c2pU(jj,j)*veln(i,j+jj,k,1,t)
                end do
              else
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windV(i,j,k,1,t)*c2pD(jj,j)*veln(i,j+jj,k,1,t)
                end do
              end if
              !--- w*d/dz ---
              if (dimens == 3) then

                if( windW(i,j,k,1,t) >= 0. ) then
                  !pgi$ unroll = n:8
                  do kk = cL(3), cU(3)
                    r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windW(i,j,k,1,t)*c3pU(kk,k)*veln(i,j,k+kk,1,t)
                  end do                                                                        
                else                                                                            
                  !pgi$ unroll = n:8                                                            
                  do jj = cL(3), cU(3)                                                          
                    r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + windW(i,j,k,1,t)*c3pD(kk,k)*veln(i,j,k+kk,1,t)
                  end do
                end if

              endif
              !--- compute diffusion -----------------------------------------------------------------------------------
              dd1 = c11u(bL(1),i)*veln(i+bL(1),j,k,1,t)
              !pgi$ unroll = n:8
              do ii = bL(1)+1, bU(1)
                dd1 = dd1 + c11u(ii,i)*veln(i+ii,j,k,1,t)
              end do
              !pgi$ unroll = n:8
              do jj = bL(2), bU(2)
                dd1 = dd1 + c22p(jj,j)*veln(i,j+jj,k,1,t)
              end do
              if( 3==dimens) then
                !pgi$ unroll = n:8
                do kk = bL(3), bU(3)
                  dd1 = dd1 + c33p(kk,k)*veln(i,j,k+kk,1,t)
                end do
              endif
              r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) - mulL*dd1

              !--- compute gradient ------------------------------------------------------------------------------------ 
              !pgi$ unroll = n:8
              do ii = gL(1), gU(1)
                r_vel(i,j,k,1,t) = r_vel(i,j,k,1,t) + cG1(ii,i)*pn(i+ii,j,k,t)
              end do

            endif

            !===========================================================================================================
            !=== computing velocity residual in y-direction ============================================================
            !===========================================================================================================
            if( SV(1)<=i .and. i<=NV(1) .and. SV(2)<=j .and. j<=NV(2) .and. SV(3)<=k .and. k<=NV(3) ) then
              !--- compute time derivative ----------------------------------------------------------------------------- 
              r_vel(i,j,k,2,t) =  mulI*( veln(i,j,k,2,t) - veln(i,j,k,2,t-1) )
              !--- compute convection -----------------------------------------------------------------------------------
              !--- u*d/dx ---
              if( windU(i,j,k,2,t) >= 0. ) then
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windU(i,j,k,2,t)*c1pU(ii,i)*veln(i+ii,j,k,2,t)
                end do
              else
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windU(i,j,k,2,t)*c1pD(ii,i)*veln(i+ii,j,k,2,t)
                end do
              end if
              !--- v*d/dy ---
              if( windV(i,j,k,2,t) >= 0. ) then
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windV(i,j,k,2,t)*c2vU(jj,j)*veln(i,j+jj,k,2,t)
                end do
              else
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windV(i,j,k,2,t)*c2vD(jj,j)*veln(i,j+jj,k,2,t)
                end do
              end if
              !--- w*d/dz ---
              if (dimens == 3) then

                if( windW(i,j,k,2,t) >= 0. ) then
                  !pgi$ unroll = n:8
                  do kk = cL(3), cU(3)
                    r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windW(i,j,k,2,t)*c3pU(kk,k)*veln(i,j,k+kk,2,t)
                  end do                                                                        
                else                                                                            
                  !pgi$ unroll = n:8                                                            
                  do jj = cL(3), cU(3)                                                          
                    r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + windW(i,j,k,2,t)*c3pD(kk,k)*veln(i,j,k+kk,2,t)
                  end do
                end if

              endif
              !--- compute diffusion ----------------------------------------------------------------------------------- 
              dd1 = c11p(bL(1),i)*veln(i+bL(1),j,k,2,t)
              !pgi$ unroll = n:8
              do ii = bL(1)+1, bU(1)
                dd1 = dd1 + c11p(ii,i)*veln(i+ii,j,k,2,t)
              end do
              !pgi$ unroll = n:8
              do jj = bL(2), bU(2)
                dd1 = dd1 + c22v(jj,j)*veln(i,j+jj,k,2,t)
              end do
              if( 3==dimens) then
                !pgi$ unroll = n:8
                do kk = bL(3), bU(3)
                  dd1 = dd1 + c33p(kk,k)*veln(i,j,k+kk,2,t)
                end do
              endif
              r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) - mulL*dd1

              !--- compute gradient ------------------------------------------------------------------------------------ 
              !pgi$ unroll = n:8
              do jj = gL(2), gU(2)
                r_vel(i,j,k,2,t) = r_vel(i,j,k,2,t) + cG2(jj,j)*pn(i,j+jj,k,t)
              end do
            endif

            !===========================================================================================================
            !=== computing velocity residual in z-direction ============================================================
            !===========================================================================================================
            if( 3==dimens .and. SW(1)<=i .and. i<=NW(1) .and. SW(2)<=j .and. j<=NW(2) .and. SW(3)<=k .and. k<=NW(3) ) then
              !--- compute time derivative ----------------------------------------------------------------------------- 
              r_vel(i,j,k,3,t) =  mulI*( veln(i,j,k,3,t) - veln(i,j,k,3,t-1) )
              !--- compute convection -----------------------------------------------------------------------------------
              !--- u*d/dx ---
              if( windU(i,j,k,3,t) >= 0. ) then
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windU(i,j,k,3,t)*c1pU(ii,i)*veln(i+ii,j,k,3,t)
                end do
              else
                !pgi$ unroll = n:8
                do ii = cL(1), cU(1)
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windU(i,j,k,3,t)*c1pD(ii,i)*veln(i+ii,j,k,3,t)
                end do
              end if
              !--- v*d/dy ---
              if( windV(i,j,k,3,t) >= 0. ) then
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windV(i,j,k,3,t)*c2pU(jj,j)*veln(i,j+jj,k,3,t)
                end do
              else
                !pgi$ unroll = n:8
                do jj = cL(2), cU(2)
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windV(i,j,k,3,t)*c2pD(jj,j)*veln(i,j+jj,k,3,t)
                end do
              end if
              !--- w*d/dz ---
              if( windW(i,j,k,3,t) >= 0. ) then
                !pgi$ unroll = n:8
                do kk = cL(3), cU(3)
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windW(i,j,k,3,t)*c3wU(kk,k)*veln(i,j,k+kk,3,t)
                end do                                                                       
              else                                                                           
                !pgi$ unroll = n:8                                                           
                do jj = cL(3), cU(3)                                                         
                  r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + windW(i,j,k,3,t)*c3wD(kk,k)*veln(i,j,k+kk,3,t)
                end do
              end if
              !--- compute diffusion ----------------------------------------------------------------------------------- 
              dd1 = c11p(bL(1),i)*veln(i+bL(1),j,k,3,t)
              !pgi$ unroll = n:8
              do ii = bL(1)+1, bU(1)
                dd1 = dd1 + c11p(ii,i)*veln(i+ii,j,k,3,t)
              end do
              !pgi$ unroll = n:8
              do jj = bL(2), bU(2)
                dd1 = dd1 + c22p(jj,j)*veln(i,j+jj,k,3,t)
              end do
              if( 3==dimens) then
                !pgi$ unroll = n:8
                do kk = bL(3), bU(3)
                  dd1 = dd1 + c33w(kk,k)*veln(i,j,k+kk,3,t)
                end do
              endif
              r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) - mulL*dd1

              !--- compute gradient ------------------------------------------------------------------------------------ 
              !pgi$ unroll = n:8
              do kk = gL(3), gU(3)
                r_vel(i,j,k,3,t) = r_vel(i,j,k,3,t) + cG3(kk,k)*pn(i,j,k+kk,t)
              end do
            endif

            !===========================================================================================================
            !=== computing pressure residual ===========================================================================
            !===========================================================================================================
            if( 3==dimens ) then

              r_p(i,j,k,t) = cD1(dL(1),i)*veln(i+dL(1),j,k,1,t)
              !pgi$ unroll = n:8
              do ii = dL(1)+1, dU(1)
                r_p(i,j,k,t) = r_p(i,j,k,t) + cD1(ii,i)*veln(i+ii,j,k,1,t)
              end do
              !pgi$ unroll = n:8
              do jj = dL(2), dU(2)
                r_p(i,j,k,t) = r_p(i,j,k,t) + cD2(jj,j)*veln(i,j+jj,k,2,t)
              end do
              !pgi$ unroll = n:8
              do kk = dL(3), dU(3)
                r_p(i,j,k,t) = r_p(i,j,k,t) + cD3(kk,k)*veln(i,j,k+kk,3,t)
              end do

            else

              r_p(i,j,k,t) = cD1(dL(1),i)*veln(i+dL(1),j,k,1,t)
              !pgi$ unroll = n:8
              do ii = dL(1)+1, dU(1)
                r_p(i,j,k,t) = r_p(i,j,k,t) + cD1(ii,i)*veln(i+ii,j,k,1,t)
              end do
              !pgi$ unroll = n:8
              do jj = dL(2), dU(2)
                r_p(i,j,k,t) = r_p(i,j,k,t) + cD2(jj,j)*veln(i,j+jj,k,2,t)
              end do

            end if

          end do
        end do
      end do
    end do

    !===========================================================================================================

  end subroutine OP_TimeNS



  !> \brief computes Time dependent NS operator
  !! is used for inner field
  !! \todo implement: generat small systems solve them, update solution
subroutine OP_TimeNSBSmoother(  &
        dimens,               &
        N,                    &
        bL,bU,                &
        BCL,BCU,              &
        cL,cU,                &
        dL,dU,                &
        gL,gU,                &
        SS,NN,                &
        SU,NU,                &
        SV,NV,                &
        SW,NW,                &
        c1uD,c2vD,c3wD,       &
        c1uU,c2vU,c3wU,       &
        c1pD,c2pD,c3pD,       &
        c1pU,c2pU,c3pU,       &
        c11p,c22p,c33p,       &
        c11u,c22v,c33w,       &
        cD1,                  &
        cD2,                  &
        cD3,                  &
        cG1,                  &
        cG2,                  &
        cG3,                  &
        mulI,                 &
        mulL,                 &
        windU,                &
        windV,                &
        windW,                &
        rhs_vel,              &
        rhs_p,                &
        vel,                  &
        p ) bind (c,name='OP_TimeNSBSmoother')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(4)

    integer(c_int), intent(in)  :: bL(4)
    integer(c_int), intent(in)  :: bU(4)

    integer(c_int), intent(in)  :: BCL(3)
    integer(c_int), intent(in)  :: BCU(3)

    integer(c_int), intent(in)  :: cL(3)
    integer(c_int), intent(in)  :: cU(3)

    integer(c_int), intent(in)  :: dL(4)
    integer(c_int), intent(in)  :: dU(4)

    integer(c_int), intent(in)  :: gL(4)
    integer(c_int), intent(in)  :: gU(4)

    integer(c_int), intent(in)  :: SS(4)
    integer(c_int), intent(in)  :: NN(4)

    integer(c_int), intent(in)  :: SU(4)
    integer(c_int), intent(in)  :: NU(4)

    integer(c_int), intent(in)  :: SV(4)
    integer(c_int), intent(in)  :: NV(4)

    integer(c_int), intent(in)  :: SW(4)
    integer(c_int), intent(in)  :: NW(4)

    real(c_double), intent(in)  :: c1uD(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2vD(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3wD(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1uU(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2vU(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3wU(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1pD(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2pD(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3pD(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c1pU(cL(1):cU(1),0:N(1))
    real(c_double), intent(in)  :: c2pU(cL(2):cU(2),0:N(2))
    real(c_double), intent(in)  :: c3pU(cL(3):cU(3),0:N(3))

    real(c_double), intent(in)  :: c11p(bL(1):bU(1),0:N(1))
    real(c_double), intent(in)  :: c22p(bL(2):bU(2),0:N(2))
    real(c_double), intent(in)  :: c33p(bL(3):bU(3),0:N(3))

    real(c_double), intent(in)  :: c11u(bL(1):bU(1),0:N(1))
    real(c_double), intent(in)  :: c22v(bL(2):bU(2),0:N(2))
    real(c_double), intent(in)  :: c33w(bL(3):bU(3),0:N(3))

    real(c_double), intent(in)  :: cD1(dL(1):dU(1),0:N(1))
    real(c_double), intent(in)  :: cD2(dL(2):dU(2),0:N(2))
    real(c_double), intent(in)  :: cD3(dL(3):dU(3),0:N(3))

    real(c_double), intent(in)  :: cG1(gL(1):gU(1),0:N(1))
    real(c_double), intent(in)  :: cG2(gL(2):gU(2),0:N(2))
    real(c_double), intent(in)  :: cG3(gL(3):gU(3),0:N(3))

    real(c_double), intent(in)  :: windU ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )
    real(c_double), intent(in)  :: windV ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )
    real(c_double), intent(in)  :: windW ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in)  :: mulI
    real(c_double), intent(in)  :: mulL
        
    real(c_double), intent(inout)  :: vel  ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(inout)  :: p    ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(out) :: rhs_vel( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(out) :: rhs_p  ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    integer(c_int)              ::  i
    integer(c_int)              ::  j
    integer(c_int)              ::  k
    integer(c_int)              ::  t

    !===========================================================================================================


    integer (c_int), parameter :: block_size = 7
    integer(c_int)             :: l
    integer(c_int)             :: ipiv(block_size)

    real (c_double) omega

    real (c_double) A(block_size,block_size)
    real (c_double) b(block_size)

     ! ----- for the test -----!
     integer(c_int) :: info
     !real(c_double) :: vt(block_size,block_size)
     !real(c_double) :: u(block_size,block_size)
     !real(c_double) :: s(block_size)
     !integer(c_int), parameter :: lwork = 5*block_size
     !real(c_double) :: work(lwork)

     ! ------------------------!

    omega = 0.7

do t = SS(4), N(4)
    do k = SW(3), NW(3)
        do j = SV(2), NV(2)
            do i = SU(1), NU(1)

                A(1:block_size,1:block_size) = 0.0
                b(1:block_size) = 0.0

                !==============================================================
                !========== assembling the matrix A ===========================
                !==============================================================

                ! diagonal: time derivative + diffusion

                A(1,1) = mulI - mulL*( c11u(0,i  ) + c22p(0,j  ) + c33p(0,k  ) )
                A(2,2) = mulI - mulL*( c11u(0,i-1) + c22p(0,j-1) + c33p(0,k-1) )
                A(3,3) = mulI - mulL*( c11p(0,i  ) + c22v(0,j  ) + c33p(0,k  ) )
                A(4,4) = mulI - mulL*( c11p(0,i-1) + c22v(0,j-1) + c33p(0,k-1) )
                A(5,5) = mulI - mulL*( c11p(0,i  ) + c22p(0,j  ) + c33w(0,k  ) )
                A(6,6) = mulI - mulL*( c11p(0,i-1) + c22p(0,j-1) + c33w(0,k-1) )

                ! sub/super-diagonal: diffusion

                A(1,2) = - mulL*c11u(bL(1),i)
                A(2,1) = - mulL*c11u(bU(1),i-1)

                A(3,4) = - mulL*c22v(bL(2),j)
                A(4,3) = - mulL*c22v(bU(2),j-1)

                A(5,6) = - mulL*c33w(bL(3),k)
                A(6,5) = - mulL*c33w(bU(3),k-1)


                ! convection
                ! ------------------ u -------------------
                ! u*d/dx
                if( windU(i,j,k,1,t) >= 0. ) then
                    A(1,1) = A(1,1) + windU(i,j,k,1,t)*c1uU(0,i)
                    A(1,2) = A(1,2) + windU(i,j,k,1,t)*c1uU(cL(1),i)
                else
                    A(1,1) = A(1,1) + windU(i,j,k,1,t)*c1uD(0,i)
                    b(1) = b(1) + windU(i,j,k,1,t)*c1uD(cU(1),i)*vel(i+1,j,k,1,t)
                end if
                ! v*d/dy
                if( windV(i,j,k,1,t) >= 0. ) then
                    A(1,1) = A(1,1) + windV(i,j,k,1,t)*c2pU(0,j)
                    b(1) = b(1) + windV(i,j,k,1,t)*c2pU(cL(2),j)*vel(i,j-1,k,1,t)
                else
                    A(1,1) = A(1,1) + windV(i,j,k,1,t)*c2pD(0,j)
                    b(1) = b(1) + windV(i,j,k,1,t)*c2pD(cU(2),j)*vel(i,j+1,k,1,t)
                end if
                ! w*d/dz
                if( windW(i,j,k,1,t) >= 0. ) then
                    A(1,1) = A(1,1) + windW(i,j,k,1,t)*c3pU(0,k)
                    b(1) = b(1) + windW(i,j,k,1,t)*c3pU(cL(3),k)*vel(i,j,k-1,1,t)
                else
                    A(1,1) = A(1,1) + windW(i,j,k,1,t)*c3pD(0,k)
                    b(1) = b(1) + windW(i,j,k,1,t)*c3pD(cU(3),k)*vel(i,j,k+1,1,t)
                end if
                ! u*d/dx
                if( windU(i-1,j,k,1,t) >= 0. ) then
                    A(2,2) = A(2,2) + windU(i-1,j,k,1,t)*c1uU(0,i-1)
                    b(2) = b(2) + windU(i-1,j,k,1,t)*c1uU(cL(1),i-1)*vel(i-1,j,k,1,t)
                else
                    A(2,2) = A(2,2) + windU(i-1,j,k,1,t)*c1uD(0,i-1)
                    A(2,1) = A(2,1) + windU(i-1,j,k,1,t)*c1uD(cU(1),i-1)
                end if
                ! v*d/dy
                if( windV(i-1,j,k,1,t) >= 0. ) then
                    A(2,2) = A(2,2) + windV(i-1,j,k,1,t)*c2pU(0,j-1)
                    b(2) = b(2) + windV(i-1,j,k,1,t)*c2pU(cL(2),j-1)*vel(i-1,j-1,k,1,t)
                else
                    A(2,2) = A(2,2) + windV(i-1,j,k,1,t)*c2pD(0,j-1)
                    b(2) = b(2) + windV(i-1,j,k,1,t)*c2pD(cU(2),j-1)*vel(i-1,j+1,k,1,t)
                end if
                ! w*d/dz
                if( windW(i-1,j,k,1,t) >= 0. ) then
                    A(2,2) = A(2,2) + windW(i-1,j,k,1,t)*c3pU(0,j-1)
                    b(2) = b(2) + windW(i-1,j,k,1,t)*c3pU(cL(3),j-1)*vel(i-1,j,k-1,1,t)
                else
                    A(2,2) = A(2,2) + windW(i-1,j,k,1,t)*c3pD(0,j-1)
                    b(2) = b(2) + windW(i-1,j,k,1,t)*c3pD(cU(3),j-1)*vel(i-1,j,k+1,1,t)
                end if


                ! ------------------ v -------------------
                ! u*d/dx
                if( windU(i,j,k,2,t) >= 0. ) then
                    A(3,3) = A(3,3) + windU(i,j,k,2,t)*c1pU(0,i)
                    b(3) = b(3) + windU(i,j,k,2,t)*c1pU(cL(1),i)*vel(i-1,j,k,2,t)
                else
                    A(3,3) = A(3,3) + windU(i,j,k,2,t)*c1pD(0,i)
                    b(3) = b(3) + windU(i,j,k,2,t)*c1pD(cU(1),i)*vel(i+1,j,k,2,t)
                end if
                ! v*d/dy
                if( windV(i,j,k,2,t) >= 0. ) then
                    A(3,3) = A(3,3) + windV(i,j,k,2,t)*c2vU(0,j)
                    A(3,4) = A(3,4) + windV(i,j,k,2,t)*c2vU(cL(2),j)
                else
                    A(3,3) = A(3,3) + windV(i,j,k,2,t)*c2vD(0,j)
                    b(3) = b(3) + windV(i,j,k,2,t)*c2vD(cU(2),j)*vel(i,j+1,k,2,t)
                end if
                ! w*d/dz
                if( windW(i,j,k,2,t) >= 0. ) then
                    A(3,3) = A(3,3) + windW(i,j,k,2,t)*c3pU(0,k)
                    b(3) = b(3) + windW(i,j,k,2,t)*c3pU(cL(3),k)*vel(i,j,k-1,2,t)
                else
                    A(3,3) = A(3,3) + windW(i,j,k,2,t)*c3pD(0,k)
                    b(3) = b(3) + windW(i,j,k,2,t)*c3pD(cU(3),k)*vel(i,j,k+1,2,t)
                end if


                ! u*d/dx
                if( windU(i,j-1,k,2,t) >= 0. ) then
                    A(4,4) = A(4,4) + windU(i,j-1,k,2,t)*c1pU(0,i-1)
                    b(4) = b(4) + windU(i,j-1,k,2,t)*c1pU(cL(1),i-1)*vel(i-1,j-1,k,2,t)
                else
                    A(4,4) = A(4,4) + windU(i,j-1,k,2,t)*c1pD(0,i-1)
                    b(4) = b(4) + windU(i,j-1,k,2,t)*c1pD(cU(1),i-1)*vel(i+1,j-1,k,2,t)
                end if
                ! v*d/dy
                if( windV(i,j-1,k,2,t) >= 0. ) then
                    A(4,4) = A(4,4) + windV(i,j-1,k,2,t)*c2vU(0,j-1)
                    b(4) = b(4) + windV(i,j-1,k,2,t)*c2vU(cL(2),j-1)*vel(i,j-2,k,2,t)
                else
                    A(4,4) = A(4,4) + windV(i,j-1,k,2,t)*c2vD(0,j-1)
                    A(4,3) = A(4,3) + windV(i,j-1,k,2,t)*c2vD(cU(2),j-1)
                end if
                ! w*d/dz
                if( windW(i,j-1,k,2,t) >= 0. ) then
                    A(4,4) = A(4,4) + windW(i,j-1,k,2,t)*c3pU(0,k-1)
                    b(4) = b(4) + windW(i,j-1,k,2,t)*c3pU(cL(3),k-1)*vel(i,j-1,k-1,2,t)
                else
                    A(4,4) = A(4,4) + windW(i,j-1,k,2,t)*c3pD(0,k-1)
                    b(4) = b(4) + windW(i,j-1,k,2,t)*c3pD(cU(3),k-1)*vel(i,j-1,k+1,2,t)
                end if

                ! ------------------ w -------------------
                ! u*d/dx
                if( windU(i,j,k,3,t) >= 0. ) then
                    A(5,5) = A(5,5) + windU(i,j,k,3,t)*c1pU(0,i)
                    b(5) = b(5) + windU(i,j,k,3,t)*c1pU(cL(1),i)*vel(i-1,j,k,3,t)
                else
                    A(5,5) = A(5,5) + windU(i,j,k,3,t)*c1pD(0,i)
                    b(5) = b(5) + windU(i,j,k,3,t)*c1pD(cU(1),i)*vel(i+1,j,k,3,t)
                end if
                    ! v*d/dy
                if( windV(i,j,k,3,t) >= 0. ) then
                    A(5,5) = A(5,5) + windV(i,j,k,3,t)*c2pU(0,j)
                    b(5) = b(5) + windV(i,j,k,3,t)*c2pU(cL(2),j)*vel(i,j-1,k,3,t)
                else
                    A(5,5) = A(5,5) + windV(i,j,k,3,t)*c2pD(0,j)
                    b(5) = b(5) + windV(i,j,k,3,t)*c2pD(cU(2),j)*vel(i,j+1,k,3,t)
                end if
                ! w*d/dz
                if( windV(i,j,k,3,t) >= 0. ) then
                    A(5,5) = A(5,5) + windW(i,j,k,3,t)*c3wU(0,k)
                    A(4,5) = A(4,5) + windW(i,j,k,3,t)*c3wU(cL(3),k)
                else
                    A(5,5) = A(5,5) + windW(i,j,k,3,t)*c3wD(0,k)
                    b(5) = b(5) + windW(i,j,k,3,t)*c3wD(cU(3),k)*vel(i,j,k+1,3,t)
                end if

                ! u*d/dx
                if( windU(i,j,k-1,3,t) >= 0. ) then
                    A(6,6) = A(6,6) + windU(i,j,k-1,3,t)*c1pU(0,i-1)
                    b(6) = b(6) + windU(i,j,k-1,3,t)*c1pU(cL(1),i-1)*vel(i-1,j,k-1,3,t)
                else
                    A(6,6) = A(6,6) + windU(i,j,k-1,3,t)*c1pD(0,i-1)
                    b(6) = b(6) + windU(i,j,k-1,3,t)*c1pD(cU(1),i-1)*vel(i+1,j,k-1,3,t)
                end if
                ! v*d/dy
                if( windV(i,j,k-1,3,t) >= 0. ) then
                    A(6,6) = A(6,6) + windV(i,j,k-1,3,t)*c2pU(0,j-1)
                    b(6) = b(6) + windV(i,j,k-1,3,t)*c2pU(cL(2),j-1)*vel(i,j-1,k-1,3,t)
                else
                    A(6,6) = A(6,6) + windV(i,j,k-1,3,t)*c2pD(0,j-1)
                    b(6) = b(6) + windV(i,j,k-1,3,t)*c2pD(cU(2),j-1)*vel(i,j+1,k-1,3,t)
                end if
                ! w*d/dz
                if( windV(i,j,k-1,3,t) >= 0. ) then
                    A(6,6) = A(6,6) + windW(i,j,k-1,3,t)*c3wU(0,k-1)
                    b(6) = b(6) + windW(i,j,k-1,3,t)*c3wU(cL(3),k-1)*vel(i,j,k-2,3,t)
                else
                    A(6,6) = A(6,6) + windW(i,j,k-1,3,t)*c3wD(0,k-1)
                    A(6,5) = A(6,5) + windW(i,j,k-1,3,t)*c3wD(cU(3),k-1)
                end if
                
                !---------- end of convection ---------------
                ! pressure gradient
                A(1:6,7) = (/ cG1(gL(1),i), cG1(gU(1),i-1), cG2(gL(2),j),cG2(gU(2),j-1), cG3(gL(3),k), cG3(gU(3),k-1) /)

                ! divergence on pressure points indeces are correct?
                A(7,1:6) = (/ cD1(dU(1),i), cD1(dL(1),i), cD2(dU(2),j), cD2(dL(2),j),cD3(dU(3),k), cD3(dL(3),k)/)

                !===================================================
                !========== assembling the RHS =====================
                !===================================================

                ! diffusion
                do l = 0,1

                        b(l+1) = b(l+1) +  mulL*(c22p(bU(1),j-l)*vel(i-l,j+1,k,1,t) + c22p(bL(1),j-l)*vel(i-l,j-1,k,1,t) + &
                        c33p(bU(1),k-l)*vel(i-l,j,k+1,1,t) + c33p(bL(1),k-l)*vel(i-l,j,k-1,1,t))

                        b(l+3) = b(l+3) + mulL*(c11p(bU(2),i-l)*vel(i+1,j-l,k,2,t) + c11p(bL(2),i-l)*vel(i-1,j-l,k,2,t) + &
                        c33p(bU(2),k-l)*vel(i,j-l,k+1,2,t) + c33p(bL(2),k-l)*vel(i,j-l,k-1,2,t))

                        b(l+5) = b(l+5) + mulL*(c22p(bU(3),j-l)*vel(i,j+1,k-l,3,t) + c22p(bL(3),j-l)*vel(i,j-1,k-l,3,t) + &
                        c11p(bU(3),i-l)*vel(i+1,j,k-l,3,t) + c11p(bL(3),i-l)*vel(i-1,j,k-l,3,t))
                end do

                b = b + mulL*(/ c11u(bU(1),i)*vel(i+1,j,k,1,t),c11u(bL(1),i-1)*vel(i-2,j,k,1,t),&
                                c22v(bU(2),j)*vel(i,j+1,k,2,t),c22v(bL(2),j-1)*vel(i,j-2,k,2,t),&
                                c33w(bU(3),k)*vel(i,j,k+1,3,t),c33w(bL(3),k-1)*vel(i,j,k-2,3,t), 0. /)

                ! pressure gradient
                b = b + (/ cG1(gU(1),i)*p(i+1,j,k,t),cG1(gL(1),i-1)*p(i-1,j,k,t),&
                            cG2(gU(2),j)*p(i,j+1,k,t),cG2(gL(2),j-1)*p(i,j-1,k,t),&
                            cG3(gU(3),k)*p(i,j,k+1,t),cG3(gL(3),k-1)*p(i,j,k-1,t), 0. /)

                ! time stencil
                b(1:6) = b(1:6) - mulI*(/ vel(i,j,k,1,t-1), vel(i-1,j,k,1,t-1),&
                                            vel(i,j,k,2,t-1), vel(i,j-1,k,2,t-1),&
                                            vel(i,j,k,3,t-1), vel(i,j,k-1,3,t-1) /)

                ! move on the rhs
                b = - b

                b = b + (/ rhs_vel(i,j,k,1,t),rhs_vel(i-1,j,k,1,t),&
                            rhs_vel(i,j,k,2,t),rhs_vel(i,j-1,k,2,t),&
                            rhs_vel(i,j,k,3,t),rhs_vel(i,j,k-1,3,t),&
                            rhs_p(i,j,k,t) /)


                ! ---------- test the matrix A --------------!
                !if (i == 1 .and. j==1 .and. k==1 .and. t==1 ) then
                !write(*,*) i,j,k,t
                !print *,'-------- the matrix A ----------'

                !do l = 1, block_size
                !write(*,'(14F7.2)') (A(l,c), c = 1, block_size)
                !end do
                !write(*,*)

                !print*, b
                !end if

                !call dgesvd ( 'A', 'A', block_size, block_size, A, block_size, s, u, block_size, vt, block_size, work, lwork , info )

                !if ( info /= 0 ) then
                !write ( *, '(a)' ) ' '
                !write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
                !return
                !end if

                !print*,'------- singular values ---------'
                !print*, s
                !-------------------------------------------!

                !print*, '! -------- Solve the matrix A --------------!'
                ! subroutine     dgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
                call dgesv( block_size, 1, A, block_size, ipiv, b, block_size, info )

                if ( info /= 0 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
                write ( *, '(a)' ) '  The matrix is numerically singular.'
                return
                end if

                !assign new solution

                vel(i,j,k,1,t) =  b(1)*omega + (1-omega)*vel(i,j,k,1,t)
                vel(i-1,j,k,1,t) =  b(2)*omega + (1-omega)*vel(i-1,j,k,1,t)
                vel(i,j,k,2,t) =  b(3)*omega + (1-omega)*vel(i,j,k,2,t)
                vel(i,j-1,k,2,t) =  b(4)*omega + (1-omega)*vel(i,j-1,k,2,t)
                vel(i,j,k,3,t) =  b(5)*omega + (1-omega)*vel(i,j,k,3,t)
                vel(i,j,k-1,3,t) =  b(6)*omega + (1-omega)*vel(i,j,k-1,3,t)

                p(i,j,k,t) =  b(7)*omega + (1-omega)*p(i,j,k,t)


            end do
        end do
    end do

    ! boundary pressure points
    go to 100 !!! 
    ! in X direction
        if (BCL(1) > 0) then
                i = SS(1)
                do k = SS(3), NN(3)
                        do j = SS(2), NN(2)
                                p(i,j,k,t) = p(i+1,j,k,t)
                        end do
                end do
         end if

        if (BCU(1) > 0) then
                i = NN(1)
                do k = SS(3), NN(3)
                        do j = SS(2), NN(2)
                                p(i,j,k,t) = p(i-1,j,k,t)
                        end do
                end do
         end if

     ! in Y direction
        if (BCL(2) > 0) then
                j = SS(2)
                do k = SS(3), NN(3)
                        do i = SS(1), NN(1)
                                p(i,j,k,t) = p(i,j+1,k,t)
                        end do
                end do
        end if

        if (BCU(2) > 0) then
                j = NN(2)
                do k = SS(3), NN(3)
                        do i = SS(1), NN(1)
                                p(i,j,k,t) = p(i,j-1,k,t)
                        end do
                end do
        end if

    ! in Z direction
        if (BCL(3) > 0) then
                k = SS(3)
                do i = SS(1), NN(1)
                        do j = SS(2), NN(2)
                                p(i,j,k,t) = p(i,j,k+1,t)
                        end do
                end do
        end if

        if (BCL(3) > 0) then
                k = NN(3)
                do i = SS(1), NN(1)
                        do j = SS(2), NN(2)
                                p(i,j,k,t) = p(i,j,k-1,t)
                        end do
                end do
        end if
100 continue !!!!!!!!!!
end do ! loop over time ends



  end subroutine OP_TimeNSBSmoother

end module cmod_TimeNSOp
