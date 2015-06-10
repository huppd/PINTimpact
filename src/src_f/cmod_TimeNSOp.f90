!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



!> \brief module providing Helmholtz
module cmod_TimeStokesOp

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
  !! \param velp velocity field array of previous time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param veln velocity field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param pn pressure field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param[out] r_vel residual velocity field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  !! \param[out] r_p residual pressure field array of next time step, can be get from \c VelocityField::getConstRawPtr    
  subroutine OP_TimeStokes( &
      dimens,               &
      N,                    &
      bL,bU,                &
      dL,dU,                &
      gL,gU,                &
      SS,NN,                &
      SU,NU,                &
      SV,NV,                &
      SW,NW,                &
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
      velp,                 &
      veln,                 &
      pn,                   &
      r_vel,                &
      r_p ) bind (c,name='OP_TimeStokes')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: dL(3)
    integer(c_int), intent(in)  :: dU(3)

    integer(c_int), intent(in)  :: gL(3)
    integer(c_int), intent(in)  :: gU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    integer(c_int), intent(in)  :: SU(3)
    integer(c_int), intent(in)  :: NU(3)

    integer(c_int), intent(in)  :: SV(3)
    integer(c_int), intent(in)  :: NV(3)

    integer(c_int), intent(in)  :: SW(3)
    integer(c_int), intent(in)  :: NW(3)

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

    real(c_double), intent(in)  :: mulI
    real(c_double), intent(in)  :: mulL

    real(c_double), intent(in)  :: velp ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3 )

    real(c_double), intent(in)  :: veln ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3 )

    real(c_double), intent(out) :: pn   ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )

    real(c_double), intent(out) :: r_vel( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3 )

    real(c_double), intent(out) :: r_p   (bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )


    real(c_double)              :: dd1

    integer(c_int)              ::  i, ii
    integer(c_int)              ::  j, jj
    integer(c_int)              ::  k, kk


    !===========================================================================================================

    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        do i = SS(1), NN(1)
          !===========================================================================================================
          !=== computing velocity residual in x-direction ============================================================
          !===========================================================================================================
          if( SU(1)<=i .and. i<=NU(1) .and. SU(2)<=j .and. NU(2)<=j .and. SU(3)<=k .and. k<=NU(3) ) then
            !--- compute time derivative ----------------------------------------------------------------------------- 
            r_vel(i,j,k,1) =  mulI*( veln(i,j,k,1) - velp(i,j,k,1) )
            !--- compute diffusion -----------------------------------------------------------------------------------
            dd1 = c11u(bL(1),i)*veln(i+bL(1),j,k,1)
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
            dd1 = dd1 + c11u(ii,i)*veln(i+ii,j,k,1)
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
            dd1 = dd1 + c22p(jj,j)*veln(i,j+jj,k,1)
            end do
            if( 3==dimens) then
              !pgi$ unroll = n:8
              do kk = bL(3), bU(3)
              dd1 = dd1 + c33p(kk,k)*veln(i,j,k+kk,1)
              end do
            endif
            r_vel(i,j,k,1) = r_vel(i,j,k,1) - mulL*dd1

            !--- compute gradient ------------------------------------------------------------------------------------ 
            !pgi$ unroll = n:8
            do ii = gL(1), gU(1)
            r_vel(i,j,k,1) = r_vel(i,j,k,1) + cG1(ii,i)*pn(i+ii,j,k)
            end do

          endif

          !===========================================================================================================
          !=== computing velocity residual in y-direction ============================================================
          !===========================================================================================================
          if( SV(1)<=i .and. i<=NV(1) .and. SV(2)<=j .and. NV(2)<=j .and. SV(3)<=k .and. k<=NV(3) ) then
            !--- compute time derivative ----------------------------------------------------------------------------- 
            r_vel(i,j,k,2) =  mulI*( veln(i,j,k,2) - velp(i,j,k,2) )
            !--- compute diffusion ----------------------------------------------------------------------------------- 
            dd1 = c11p(bL(1),i)*veln(i+bL(1),j,k,2)
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
            dd1 = dd1 + c11p(ii,i)*veln(i+ii,j,k,2)
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
            dd1 = dd1 + c22v(jj,j)*veln(i,j+jj,k,2)
            end do
            if( 3==dimens) then
              !pgi$ unroll = n:8
              do kk = bL(3), bU(3)
              dd1 = dd1 + c33p(kk,k)*veln(i,j,k+kk,2)
              end do
            endif
            r_vel(i,j,k,2) = r_vel(i,j,k,2) - mulL*dd1

            !--- compute gradient ------------------------------------------------------------------------------------ 
            !pgi$ unroll = n:8
            do jj = gL(2), gU(2)
              r_vel(i,j,k,2) = r_vel(i,j,k,2) + cG2(jj,j)*pn(i,j+jj,k)
            end do
          endif

          !===========================================================================================================
          !=== computing velocity residual in z-direction ============================================================
          !===========================================================================================================
          if( 3==dimens .and. SW(1)<=i .and. i<=NW(1) .and. SW(2)<=j .and. NW(2)<=j .and. SW(3)<=k .and. k<=NW(3) ) then
            !--- compute time derivative ----------------------------------------------------------------------------- 
            r_vel(i,j,k,3) =  mulI*( veln(i,j,k,3) - velp(i,j,k,3) )
            !--- compute diffusion ----------------------------------------------------------------------------------- 
            dd1 = c11p(bL(1),i)*veln(i+bL(1),j,k,3)
            !pgi$ unroll = n:8
            do ii = bL(1)+1, bU(1)
            dd1 = dd1 + c11p(ii,i)*veln(i+ii,j,k,3)
            end do
            !pgi$ unroll = n:8
            do jj = bL(2), bU(2)
            dd1 = dd1 + c22p(jj,j)*veln(i,j+jj,k,3)
            end do
            if( 3==dimens) then
              !pgi$ unroll = n:8
              do kk = bL(3), bU(3)
              dd1 = dd1 + c33w(kk,k)*veln(i,j,k+kk,3)
              end do
            endif
            r_vel(i,j,k,3) = r_vel(i,j,k,3) - mulL*dd1

            !--- compute gradient ------------------------------------------------------------------------------------ 
            !pgi$ unroll = n:8
            do kk = gL(3), gU(3)
              r_vel(i,j,k,3) = r_vel(i,j,k,3) + cG3(kk,k)*pn(i,j,k+kk)
            end do
          endif

          !===========================================================================================================
          !=== computing pressure residual ===========================================================================
          !===========================================================================================================
          if( 3==dimens ) then

            r_p(i,j,k) = cD1(dL(1),i)*veln(i+dL(1),j,k,1)
            !pgi$ unroll = n:8
            do ii = dL(1)+1, dU(1)
              r_p(i,j,k) = r_p(i,j,k) + cD1(ii,i)*veln(i+ii,j,k,1)
            end do
            !pgi$ unroll = n:8
            do jj = dL(2), dU(2)
              r_p(i,j,k) = r_p(i,j,k) + cD2(jj,j)*veln(i,j+jj,k,2)
            end do
            !pgi$ unroll = n:8
            do kk = dL(3), dU(3)
              r_p(i,j,k) = r_p(i,j,k) + cD3(kk,k)*veln(i,j,k+kk,3)
            end do

          else

            r_p(i,j,k) = cD1(dL(1),i)*veln(i+dL(1),j,k,1)
            !pgi$ unroll = n:8
            do ii = dL(1)+1, dU(1)
              r_p(i,j,k) = r_p(i,j,k) + cD1(ii,i)*veln(i+ii,j,k,1)
            end do
            !pgi$ unroll = n:8
            do jj = dL(2), dU(2)
              r_p(i,j,k) = r_p(i,j,k) + cD2(jj,j)*veln(i,j+jj,k,2)
            end do

          end if

        end do
      end do
    end do

    !===========================================================================================================

  end subroutine OP_TimeStokes


end module cmod_TimeStokesOp
