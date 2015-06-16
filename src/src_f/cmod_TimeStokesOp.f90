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
      !velp,                 &
      veln,                 &
      pn,                   &
      r_vel,                &
      r_p ) bind (c,name='OP_TimeStokes')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(4)

    integer(c_int), intent(in)  :: bL(4)
    integer(c_int), intent(in)  :: bU(4)

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

    real(c_double), intent(in)  :: c11p( bL(1):bU(1), 0:N(1) )
    real(c_double), intent(in)  :: c22p( bL(2):bU(2), 0:N(2) )
    real(c_double), intent(in)  :: c33p( bL(3):bU(3), 0:N(3) )

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

    !real(c_double), intent(in)  :: velp ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3 )

    real(c_double), intent(in)  :: veln ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, (bL(4)+1):(N(4)+bU(4)) )

    real(c_double), intent(in)  :: pn   ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), (bL(4)+1):(N(4)+bU(4)))

    real(c_double), intent(out) :: r_vel( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, (bL(4)+1):(N(4)+bU(4)) )

    real(c_double), intent(out) :: r_p   (bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), (bL(4)+1):(N(4)+bU(4)) )


    real(c_double)              :: dd1

    integer(c_int)              ::  i, ii
    integer(c_int)              ::  j, jj
    integer(c_int)              ::  k, kk
    integer(c_int)              ::  t 


    !===========================================================================================================

    do t = SS(4), NN(4)
      do k = SS(3), NN(3)
        do j = SS(2), NN(2)
          do i = SS(1), NN(1)
            !===========================================================================================================
            !=== computing velocity residual in x-direction ============================================================
            !===========================================================================================================
            if( SU(1)<=i .and. i<=NU(1) .and. SU(2)<=j .and. j<=NU(2) .and. SU(3)<=k .and. k<=NU(3) ) then
              !--- compute time derivative ----------------------------------------------------------------------------- 
              r_vel(i,j,k,1,t) =  mulI*( veln(i,j,k,1,t) - veln(i,j,k,1,t-1) )
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

  end subroutine OP_TimeStokes


  !> \brief computes Time dependent Stokes operator
  !! is used for inner field
  !! \todo implement: generat small systems solve them, update solution
  subroutine OP_TimeStokesBSmoother( &
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
      vel,                  &
      p,                    &
      r_vel,                &
      r_p ) bind (c,name='OP_TimeStokesBSmoother')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(4)

    integer(c_int), intent(in)  :: bL(4)
    integer(c_int), intent(in)  :: bU(4)

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
        
    real(c_double), intent(in)  :: vel ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, (bL(4)+1):(N(4)+bU(4)) )

    real(c_double), intent(in)  :: p   ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), (bL(4)+1):(N(4)+bU(4)))

    real(c_double), intent(out) :: r_vel( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, (bL(4)+1):(N(4)+bU(4)))

    real(c_double), intent(out) :: r_p   (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), (bL(4)+1):(N(4)+bU(4)))

    integer(c_int)              ::  i, ii
    integer(c_int)              ::  j, jj
    integer(c_int)              ::  k, kk
    integer(c_int)              ::  t

    !===========================================================================================================


    integer (c_int), parameter :: block_size = 14 ! decompose this in time points + pressure...
    integer(c_int)             :: l
    integer(c_int)             :: ll
    integer(c_int)             :: ipiv(block_size)
     

    real (c_double) A(block_size,block_size) ! fill with zeros?
    real (c_double) b(block_size)

     ! ----- for the test -----!
     integer(c_int) :: info
     !real(c_double) :: vt(block_size,block_size)
     !real(c_double) :: u(block_size,block_size)
     !real(c_double) :: s(block_size)
     !integer(c_int), parameter :: lwork = 5*block_size
     !real(c_double) :: work(lwork)

     ! ------------------------!

     A(1:block_size,1:block_size) = 0.0

    do t = SS(4), NN(4) - 2 ! here should be -1
    do k = SW(3), NW(3)
      do j = SV(2), NV(2)
        do i = SU(1), NU(1)

            ! from dL to DU loop? no for b.s. just small stencil. time
            ! periodcity?
            ! if no strched grids you can compute it just one time

        !==============================================================
        !========== assembling the matrix A ===========================
        !==============================================================

            ! diagonal: time derivative + diffusion (sign of dt?)

            A(1,1) = mulI -  mulL*c11u(bL(1)+1,i)
            A(2,2) = mulI -  mulL*c11u(bL(1)+1,i-1)

            A(3,3) = mulI -  mulL*c22v(bL(2)+1,j)
            A(4,4) = mulI -  mulL*c22v(bL(2)+1,j-1)

            A(5,5) = mulI -  mulL*c33w(bL(3)+1,k)
            A(6,6) = mulI -  mulL*c33w(bL(3)+1,k-1)

            ! sub/super-diagonal: diffusion

            A(1,2) = - mulL*c11u(bL(1),i)
            A(2,1) = - mulL*c11u(bU(1),i-1)

            A(3,4) = - mulL*c22v(bL(2),j)
            A(4,3) = - mulL*c22v(bU(2),j-1)

            A(5,6) = - mulL*c33w(bL(3),k)
            A(6,5) = - mulL*c33w(bU(3),k-1)

            ! copy in the lower block for the next time point

            A(7:12,7:12) = A(1:6,1:6)

            ! pressure gradient 
            A(1:6,13) = (/ cG1(gL(1),i), cG1(gU(1),i-1), cG2(gL(2),j),cG2(gU(2),j-1), cG3(gL(3),k), cG3(gU(3),k-1) /)
            A(7:12,14) = A(1:6,13)

            ! divergence on pressure points indeces are correct?
            A(13,1:6) = (/ cD1(dU(1),i), cD1(dL(1),i), cD2(dU(2),j), cD2(dL(2),j),cD3(dU(3),k), cD3(dL(3),k)/)
            A(14,7:12) = A(13,1:6)

            ! time derivative 

            do l = 1,6
            A(l+6,l) = - mulI
            end do

            !===================================================
            !========== assembling the RHS =====================
            !===================================================

            ! diffusion
                
            b = mulL*(/ c11u(bU(1),i)*vel(i+1,j,k,1,t),c11u(bL(1),i-1)*vel(i-2,j,k,1,t),&
                       c22v(bU(2),j)*vel(i,j+1,k,2,t),c22v(bL(2),j-1)*vel(i,j-2,k,2,t),&
                       c33w(bU(3),k)*vel(i,j,k+1,3,t),c33w(bL(3),j-1)*vel(i,j,k-2,3,t),&
                       c11u(bU(1),i)*vel(i+1,j,k,1,t+1),c11u(bL(1),i-1)*vel(i-2,j,k,1,t+1),&
                       c22v(bU(2),j)*vel(i,j+1,k,2,t+1),c22v(bL(2),j-1)*vel(i,j-2,k,2,t+1),&
                       c33w(bU(3),k)*vel(i,j,k+1,3,t+1),c33w(bL(3),j-1)*vel(i,j,k-2,3,t+1),&
                       0, 0 /) ! -> divergence

            ! pressure gradient
            b = b + (/ cG1(gU(1),i)*p(i+1,j,k,t),cG1(gL(1),i-1)*p(i-1,j,k,t),&
                       cG2(gU(2),j)*p(i,j+1,k,t),cG2(gL(2),j-1)*p(i,j-1,k,t),&
                       cG3(gU(3),k)*p(i,j,k+1,t),cG3(gL(3),k-1)*p(i,j,k-1,t),&
                       cG1(gU(1),i)*p(i+1,j,k,t+1),cG1(gL(1),i-1)*p(i-1,j,k,t+1),&
                       cG2(gU(2),j)*p(i,j+1,k,t+1),cG2(gL(2),j-1)*p(i,j-1,k,t+1),&
                       cG3(gU(3),k)*p(i,j,k+1,t+1),cG3(gL(3),k-1)*p(i,j,k-1,t+1),&
                       0, 0 /) ! -> divergence

            ! time stencil (just in the second time slice)
            b(1:6) = b(1:6) - mulI*(/ vel(i,j,k,1,t-1), vel(i-1,j,k,1,t-1),&
                                        vel(i,j,k,2,t-1), vel(i,j-1,k,2,t-1),&
                                        vel(i,j,k,3,t-1), vel(i,j,k-1,3,t-1) /)

! ---------- test the matrix A --------------!

!write(*,*) i,j,k

!print *,'-------- the matrix A ----------'

!do l = 1, block_size
!write(*,'(14F11.3)') (A(l,ll), ll = 1, block_size) 
!end do
!write(*,*)


!call dgesvd ( 'A', 'A', block_size, block_size, A, block_size, s, u, block_size, vt, block_size, work, lwork , info )

!if ( info /= 0 ) then
!write ( *, '(a)' ) ' '
!write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
!return
!end if

!print*,'------- singular values ---------'
!print*, s

! -------------------------------------------!


!print*, '! -------- Factor the matrix A --------------!'

  call dgetrf ( block_size, block_size, A, block_size, ipiv, info )
  
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
    write ( *, '(a)' ) '  The matrix is numerically singular.'
    return
  end if
  
! print*, '------ matrix A -----'
!do l = 1, block_size
!write(*,'(14F8.1)') (A(l,ll), ll = 1, block_size) 
!end do
!write(*,*)
  !print*, '---------------------'


  call dgetrs ( 'N', block_size, 1, A, block_size, ipiv, b, block_size, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution procedure failed!'
    write ( *, '(a,i8)' ) '  INFO = ', info
    return
  end if

!print*, '------------ result -------------'
!print*, b
!print*, '---------------------------------'

        ! assign the solution (more general)     
        
        !print*, '------ max -------'
        !print*, NU(1),NV(2),NW(3),NN(4)
        !print*, '------- i j k t+1 -------'
        !print*, i,j,k,t+1
       
        r_vel(i,j,k,1,t) = b(1)
        r_vel(i,j,k,1,t+1) = b(7)
        !r_vel(i-1,j,k,1,t:t+1) = (/ b(2), b(8) /)
        !r_vel(i,j,k,2,t:t+1) = (/ b(3), b(9) /)
        !r_vel(i,j-1,k,2,t:t+1) = (/ b(4), b(10) /)
        !r_vel(i,j,k,3,t:t+1) = (/ b(5), b(11) /)
        !r_vel(i,j,k-1,3,t:t+1) = (/ b(6), b(12) /) 

        !r_p(i,j,k,t:t+1) = b(13:14)

        end do
      end do
    end do
    end do

  end subroutine OP_TimeStokesBSmoother

end module cmod_TimeStokesOp
