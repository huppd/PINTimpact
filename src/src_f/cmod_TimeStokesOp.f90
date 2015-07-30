!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************



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
      BCL,BCU,              &
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
      rhs_vel,              &
      rhs_p,                &
      vel,                  &
      p ) bind (c,name='OP_TimeStokesBSmoother')


    implicit none

    integer(c_int), intent(in)  :: dimens

    integer(c_int), intent(in)  :: N(4)

    integer(c_int), intent(in)  :: bL(4)
    integer(c_int), intent(in)  :: bU(4)

    integer(c_int), intent(in)  :: BCL(3)
    integer(c_int), intent(in)  :: BCU(3)

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
        
    real(c_double), intent(inout)  :: vel  ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(inout)  :: p    ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in) :: rhs_vel( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in) :: rhs_p  ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    integer(c_int)              ::  i
    integer(c_int)              ::  j
    integer(c_int)              ::  k
    integer(c_int)              ::  t

    !===========================================================================================================


    integer(c_int), parameter :: block_size = 14 ! decompose this in time points + pressure...
    integer(c_int)             :: l, ll, d, c
    integer(c_int)             :: ipiv(block_size)
     

    real (c_double) A(block_size,block_size) ! fill with zeros?
    real (c_double) b(block_size)

    integer(c_int) :: info

     ! ----- for the test -----!
     !real(c_double) :: vt(block_size,block_size)
     !real(c_double) :: u(block_size,block_size)
     !real(c_double) :: s(block_size)
     !integer(c_int), parameter :: lwork = 5*block_size
     !real(c_double) :: work(lwork)

     ! ------------------------!
! this is bad
do t = SS(4), N(4) - 1
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        do i = SS(1), NN(1)


          
    !    if (    (i==SU(1) .and. j==SV(2)) .or. &
     !           (i==SU(1) .and. j==NV(2)) .or. & 
      !          (i==SU(1) .and. k==SW(3)) .or. & 
       !         (i==SU(1) .and. k==NW(3)) .or. & 
        !        (i==NU(1) .and. j==SV(2)) .or. & 
         !       (i==NU(1) .and. j==NV(2)) .or. & 
          !      (i==NU(1) .and. k==SW(3)) .or. & 
           !     (i==NU(1) .and. k==SW(3)) .or. & 
            !    (j==SV(2) .and. k==SW(3)) .or. & 
             !   (j==SV(2) .and. k==SW(3)) .or. &
              !  (j==SV(2) .and. k==NW(3)) .or. &
               ! (j==SV(2) .and. k==NW(3))  ) then
        !if (i==SS(1) .or. i==NN(1) .or. j==SS(2) .or. j==NN(2) .or. k==SS(3) .or. k==NN(3)) then         
                
                p(i,j,k,t) = 0
                ! this is not working. why?
        !end if       

   end do
  end do
 end do
end do


        
    do t = SS(4), N(4) - 1 
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

            ! copy in the lower block for the next time point

            A(7:12,7:12) = A(1:6,1:6)

            ! pressure gradient 
            A(1:6,13) = (/ cG1(gL(1),i), cG1(gU(1),i-1), cG2(gL(2),j),cG2(gU(2),j-1), cG3(gL(3),k), cG3(gU(3),k-1) /)
            A(7:12,14) = A(1:6,13)

            ! divergence on pressure points indeces are correct?
            A(13,1:6) = (/ cD1(dU(1),i), cD1(dL(1),i), cD2(dU(2),j), cD2(dL(2),j),cD3(dU(3),k), cD3(dL(3),k)/)
            A(14,7:12) = A(13,1:6)
                
            ! time derivative 

            do d = 1,6
            A(d+6,d) = - mulI
            end do

            !===================================================
            !========== assembling the RHS =====================
            !===================================================
                
                ! boundary term !

            ! diffusion 
            do l = 0,1
               do ll =0,1
                
                b(l+1+6*ll) = c22p(bU(1),j-l)*vel(i-l,j+1,k,1,t+ll) + c22p(bL(1),j-l)*vel(i-l,j-1,k,1,t+ll) + &
                           c33p(bU(1),k-l)*vel(i-l,j,k+1,1,t+ll) + c33p(bL(1),k-l)*vel(i-l,j,k-1,1,t+ll) 
                        
                b(l+3+6*ll) = c11p(bU(2),i-l)*vel(i+1,j-l,k,2,t+ll) + c11p(bL(2),i-l)*vel(i-1,j-l,k,2,t+ll) + &
                           c33p(bU(2),k-l)*vel(i,j-l,k+1,2,t+ll) + c33p(bL(2),k-l)*vel(i,j-l,k-1,2,t+ll) 

                b(l+5+6*ll) = c22p(bU(3),j-l)*vel(i,j+1,k-l,3,t+ll) + c22p(bL(3),j-l)*vel(i,j-1,k-l,3,t+ll) + &
                           c11p(bU(3),i-l)*vel(i+1,j,k-l,3,t+ll) + c11p(bL(3),i-l)*vel(i-1,j,k-l,3,t+ll) 
                end do
           end do
            
        b = b + (/ c11u(bU(1),i)*vel(i+1,j,k,1,t),c11u(bL(1),i-1)*vel(i-2,j,k,1,t),&
                       c22v(bU(2),j)*vel(i,j+1,k,2,t),c22v(bL(2),j-1)*vel(i,j-2,k,2,t),&
                       c33w(bU(3),k)*vel(i,j,k+1,3,t),c33w(bL(3),k-1)*vel(i,j,k-2,3,t),&
                       c11u(bU(1),i)*vel(i+1,j,k,1,t+1),c11u(bL(1),i-1)*vel(i-2,j,k,1,t+1),&
                       c22v(bU(2),j)*vel(i,j+1,k,2,t+1),c22v(bL(2),j-1)*vel(i,j-2,k,2,t+1),&
                       c33w(bU(3),k)*vel(i,j,k+1,3,t+1),c33w(bL(3),k-1)*vel(i,j,k-2,3,t+1),&
                       0., 0. /) ! -> divergence
                
            b = mulL*b
            
            ! pressure gradient
            b = b + (/ cG1(gU(1),i)*p(i+1,j,k,t),cG1(gL(1),i-1)*p(i-1,j,k,t),&
                       cG2(gU(2),j)*p(i,j+1,k,t),cG2(gL(2),j-1)*p(i,j-1,k,t),&
                       cG3(gU(3),k)*p(i,j,k+1,t),cG3(gL(3),k-1)*p(i,j,k-1,t),&
                       cG1(gU(1),i)*p(i+1,j,k,t+1),cG1(gL(1),i-1)*p(i-1,j,k,t+1),&
                       cG2(gU(2),j)*p(i,j+1,k,t+1),cG2(gL(2),j-1)*p(i,j-1,k,t+1),&
                       cG3(gU(3),k)*p(i,j,k+1,t+1),cG3(gL(3),k-1)*p(i,j,k-1,t+1),&
                       0., 0. /) ! -> divergence
            
            ! time stencil (just in the first time slice)
            b(1:6) = b(1:6) - mulI*(/ vel(i,j,k,1,t-1), vel(i-1,j,k,1,t-1),&
                                        vel(i,j,k,2,t-1), vel(i,j-1,k,2,t-1),&
                                        vel(i,j,k,3,t-1), vel(i,j,k-1,3,t-1) /)
            ! new
            b = - b
            
            b = b + (/ rhs_vel(i,j,k,1,t),rhs_vel(i-1,j,k,1,t),&
                       rhs_vel(i,j,k,2,t),rhs_vel(i,j-1,k,2,t),&
                       rhs_vel(i,j,k,3,t),rhs_vel(i,j,k-1,3,t),&
                    rhs_vel(i,j,k,1,t+1),rhs_vel(i-1,j,k,1,t+1),&
                    rhs_vel(i,j,k,2,t+1),rhs_vel(i,j-1,k,2,t+1),&
                    rhs_vel(i,j,k,3,t+1),rhs_vel(i,j,k-1,3,t+1),&
                    rhs_p(i,j,k,t),rhs_p(i,j,k,t+1) /)


! ---------- test the matrix A --------------!
!if (i==1 .and. k==1 .and. j==1 .and. t==1 ) then

!write(*,*) i,j,k
!print *,'-------- the matrix A ----------'

!do l = 1, block_size
!write(*,'(14F7.2)') (A(l,c), c = 1, block_size) 
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
       
        vel(i,j,k,1,t) =  b(1)
        vel(i-1,j,k,1,t) =  b(2)
        vel(i,j,k,2,t) =  b(3)
        vel(i,j-1,k,2,t) =  b(4)
        vel(i,j,k,3,t) =  b(5)
        vel(i,j,k-1,3,t) =  b(6)

        vel(i,j,k,1,t+1) =  b(7)
        vel(i-1,j,k,1,t+1) =  b(8)
        vel(i,j,k,2,t+1) =  b(9)
        vel(i,j-1,k,2,t+1) =  b(10)
        vel(i,j,k,3,t+1) =  b(11)
        vel(i,j,k-1,3,t+1) =  b(12)

        !p(i,j,k,t) =  b(13)
        !p(i,j,k,t+1) =  b(14)
        
        end do
      end do
    end do

go to 100
    ! boundary pressure points
       
    ! in X direction
        if (BCL(1) > 0) then
                i = SS(1)
                do k = SW(3), NW(3)
                        do j = SV(2), NV(2)
                            do ll = 0,1

                                p(i,j,k,t+ll) =  ( rhs_vel(i,j,k,1,t+ll) - (   &
                                                  cG1(gU(1),i)*p(i+1,j,k,t+ll) + mulI*(vel(i,j,k,1,t+ll) - vel(i,j,k,1,t-1+ll)) &
                                                  - mulL*(vel(i,j,k,1,t+ll)*(c11u(0,i)+c22p(0,j)+c33p(0,k)) +            &
                                                        c22p(bU(1),j)*vel(i,j+1,k,1,t+ll) + c22p(bL(1),j)*vel(i,j-1,k,1,t+ll) + &
                                                        c33p(bU(1),k)*vel(i,j,k+1,1,t+ll) + c33p(bL(1),k)*vel(i,j,k-1,1,t+ll) + &
                                                        c11u(bU(1),i)*vel(i+1,j,k,1,t+ll)))) /cG1(gL(1),i)
                                end do
                        end do
                end do
        end if

        if (BCU(1) > 0) then
                i = NN(1)
                do k = SW(3), NW(3)
                        do j = SV(2), NV(2)
                            do ll = 0,1


                                p(i,j,k,t+ll) =  ( rhs_vel(i-1,j,k,1,t+ll) - (  &
                                                   cG1(gL(1),i-1)*p(i-1,j,k,t+ll) + mulI*(vel(i-1,j,k,1,t+ll) - vel(i-1,j,k,1,t-1+ll)) &
                                                   - mulL*(vel(i-1,j,k,1,t+ll)*(c11u(0,i-1)+c22p(0,j) +c33p(0,k)) +            &
                                                        c22p(bU(1),j)*vel(i-1,j+1,k,1,t+ll) + c22p(bL(1),j)*vel(i-1,j-1,k,1,t+ll) + &
                                                        c33p(bU(1),k)*vel(i-1,j,k+1,1,t+ll) + c33p(bL(1),k)*vel(i-1,j,k-1,1,t+ll) + &
                                                        c11u(bL(1),i-1)*vel(i-2,j,k,1,t+ll)))) /cG1(gU(1),i-1)
                           end do
                        end do
                end do
        end if
    
        ! in Y direction
        if (BCL(2) > 0) then
                j = SS(2)
                do k = SW(3), NW(3)
                        do i = SU(1), NU(1)
                            do ll = 0,1

                                p(i,j,k,t+ll) =  ( rhs_vel(i,j,k,2,t+ll) - ( &
                                                   cG2(gU(2),j)*p(i,j+1,k,t+ll) + mulI*(vel(i,j,k,2,t+ll) - vel(i,j,k,2,t-1+ll)) &
                                                   -mulL*(vel(i,j,k,2,t+ll)*(c11p(0,i)+c22v(0,j)+c33p(0,k)) +            &
                                                        c11p(bU(2),i)*vel(i+1,j,k,2,t+ll) + c11p(bL(2),i)*vel(i-1,j,k,2,t+ll) + &
                                                        c33p(bU(2),k)*vel(i,j,k+1,2,t+ll) + c33p(bL(2),k)*vel(i,j,k-1,2,t+ll) + &
                                                        c22v(bU(2),j)*vel(i,j+1,k,2,t+ll)))) /cG2(gL(2),j)
                            end do
                        end do
                end do
        end if

        if (BCU(2) > 0) then
                j = NN(2)
                do k = SW(3), NW(3)
                        do i = SU(1), NU(1)
                            do ll = 0,1

                                p(i,j,k,t+ll) =  ( rhs_vel(i,j-1,k,2,t+ll) - ( &
                                                   cG2(gL(2),j-1)*p(i,j-1,k,t+ll) + mulI*(vel(i,j-1,k,2,t+ll) - vel(i,j-1,k,2,t-1+ll)) &
                                                   -mulL*(vel(i,j-1,k,2,t+ll)*(c11p(0,i)+c22v(0,j-1)+c33p(0,k))+            &
                                                        c11p(bU(2),i)*vel(i+1,j-1,k,2,t+ll) + c11p(bL(2),i)*vel(i-1,j-1,k,2,t+ll) + &
                                                        c33p(bU(2),k)*vel(i,j-1,k+1,2,t+ll) + c33p(bL(2),k)*vel(i,j-1,k-1,2,t+ll) + &
                                                        c22v(bL(2),j-1)*vel(i,j-2,k,2,t+ll)))) /cG2(gU(2),j-1)
                            end do
                        end do
                end do
        end if
   
    ! in Z direction
        if (BCL(3) > 0) then
                k = SS(3)
                do i = SU(1), NU(1)
                        do j = SV(2), NV(2)
                            do ll = 0,1

                                p(i,j,k,t+ll) =  (rhs_vel(i,j,k,3,t+ll) - ( &
                                                  cG3(gU(3),k)*p(i,j,k+1,t+ll) + mulI*(vel(i,j,k,3,t+ll) - vel(i,j,k,3,t-1+ll)) &
                                                  -mulL*(vel(i,j,k,3,t+ll)*(c11p(0,i)+c22p(0,j)+c33w(0,k)) +            &
                                                        c22p(bU(3),j)*vel(i,j+1,k,3,t+ll) + c22p(bL(3),j)*vel(i,j-1,k,3,t+ll) + &
                                                        c11p(bU(3),i)*vel(i+1,j,k,3,t+ll) + c11p(bL(3),i)*vel(i-1,j,k,3,t+ll) + &
                                                        c33w(bU(3),k)*vel(i,j,k+1,3,t+ll)))) /cG3(gL(3),k)
                                end do
                        end do
                end do
        end if
        
        if (BCU(3) > 0) then
                k = NN(3)
                do i = SU(1), NU(1)
                        do j = SV(2), NV(2)
                            do ll = 0,1

                                
                                p(i,j,k,t+ll) =  (rhs_vel(i,j,k-1,3,t+ll) - ( &
                                                  cG3(gL(3),k-1)*p(i,j,k-1,t+ll) + mulI*(vel(i,j,k-1,3,t+ll) - vel(i,j,k-1,3,t-1+ll)) &
                                                  -mulL*(vel(i,j,k-1,3,t+ll)*(c11p(0,i)+c22p(0,j)+c33w(0,k-1)) +            &
                                                        c22p(bU(3),j)*vel(i,j+1,k-1,3,t+ll) + c22p(bL(3),j)*vel(i,j-1,k-1,3,t+ll) + &
                                                        c11p(bU(3),i)*vel(i+1,j,k-1,3,t+ll)+ c11p(bL(3),i)*vel(i-1,j,k-1,3,t+ll) + &
                                                        c33w(bL(3),k-1)*vel(i,j,k-2,3,t+ll)))) / cG3(gU(3),k-1)
                            end do
                        end do
                end do
        end if
100 continue
    end do

  end subroutine OP_TimeStokesBSmoother

! new stuff

  subroutine OP_TimeStokesLSmoother( &
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
      rhs_vel,              &
      rhs_p,                &
      vel,                  &
      p,                    &
      t_size) bind (c,name='OP_TimeStokesLSmoother')


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
        
    real(c_double), intent(inout)  :: vel  ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(inout)  :: p    ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),   0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in) :: rhs_vel  ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)), 1:3, 0:(N(4)+bU(4)-bL(4)) )

    real(c_double), intent(in) :: rhs_p    ( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)),      0:(N(4)+bU(4)-bL(4)) )

    integer(c_int), intent(in)  :: t_size

    integer(c_int)              ::  i
    integer(c_int)              ::  j
    integer(c_int)              ::  k
    integer(c_int)              ::  t

    !===========================================================================================================

    integer(c_int)             :: block_size 
    integer(c_int)             :: l, ll, d, c
    integer(c_int), allocatable :: ipiv(:)
        
    real (c_double), allocatable :: A(:,:) 
    real (c_double), allocatable :: b(:)

    integer(c_int) :: info

    block_size = 7*t_size

    allocate(A(block_size,block_size),b(block_size),ipiv(block_size))
 
    do t = SS(4),N(4) - (t_size-1),t_size-1
    do k = SW(3), NW(3)
      do j = SV(2), NV(2)
        do i = SU(1), NU(1)
            
         A(1:block_size,1:block_size) = 0.0
         b(1:block_size) = 0.0

        !==============================================================
        !========== assembling the matrix A ===========================
        !==============================================================
                
            ! diagonal: time derivative + diffusion 
            
            A(1,1) = mulI - mulL*(c11u(0,i) + c22p(0,i) + c33p(0,i))
            A(2,2) = mulI - mulL*(c11u(0,i-1) + c22p(0,i-1) + c33p(0,i-1))
            A(3,3) = mulI - mulL*(c22v(0,j) + c11p(0,j) + c33p(0,j))
            A(4,4) = mulI - mulL*(c22v(0,j-1) + c11p(0,j-1) + c33p(0,j-1))
            A(5,5) = mulI - mulL*(c33w(0,k) + c11p(0,k) + c22p(0,k))
            A(6,6) = mulI - mulL*(c33w(0,k-1) + c11p(0,k-1) + c22p(0,k-1))
            
            ! sub/super-diagonal: diffusion
                
            A(1,2) = - mulL*c11u(bL(1),i)
            A(2,1) = - mulL*c11u(bU(1),i-1)

            A(3,4) = - mulL*c22v(bL(2),j)
            A(4,3) = - mulL*c22v(bU(2),j-1)

            A(5,6) = - mulL*c33w(bL(3),k)
            A(6,5) = - mulL*c33w(bU(3),k-1)

            ! pressure gradient
            A(1:6,block_size-t_size+1) = (/ cG1(gL(1),i), cG1(gU(1),i-1),cG2(gL(2),j),cG2(gU(2),j-1), cG3(gL(3),k), cG3(gU(3),k-1) /)

            ! divergence on pressure points indeces are correct?
            A(block_size-t_size+1,1:6) = (/ cD1(dU(1),i), cD1(dL(1),i), cD2(dU(2),j),cD2(dL(2),j),cD3(dU(3),k), cD3(dL(3),k)/)

            ! copy for others time slices
            if ( t_size > 1) then
                do l = 1, t_size - 1
                    ! diffusion
                    A(l*6+1:(l+1)*6,l*6+1:(l+1)*6) = A(1:6,1:6)

                    ! p gradient
                    A(l*6+1:(l+1)*6,6*t_size + 1 + l) = A(1:6,6*t_size + 1)

                    ! divergence
                    A(6*t_size + 1 + l,l*6+1:(l+1)*6) = A(6*t_size + 1,1:6)

                    ! time derivative
                    do d = l*6 - 5, l*6
                        A(d+6,d) = - mulI
                    end do

                end do
            end if

            !===================================================
            !========== assembling the RHS =====================
            !===================================================

            ! diffusion 
            do ll = 0,t_size-1

                do l = 0,1

                    b(l+1+6*ll) = c22p(bU(1),i-l)*vel(i-l,j+1,k,1,t+ll) + c22p(bL(1),i-l)*vel(i-l,j-1,k,1,t+ll) + &
                           c33p(bU(1),i-l)*vel(i-l,j,k+1,1,t+ll) + c33p(bL(1),i-l)*vel(i-l,j,k-1,1,t+ll) 
                        
                    b(l+3+6*ll) = c11p(bU(2),j-l)*vel(i+1,j-l,k,2,t+ll) + c11p(bL(2),j-l)*vel(i-1,j-l,k,2,t+ll) + &
                           c33p(bU(2),j-l)*vel(i,j-l,k+1,2,t+ll) + c33p(bL(2),j-l)*vel(i,j-l,k-1,2,t+ll) 

                    b(l+5+6*ll) = c22p(bU(3),k-l)*vel(i,j+1,k-l,3,t+ll) + c22p(bL(3),k-l)*vel(i,j-1,k-l,3,t+ll) + &
                           c11p(bU(3),k-l)*vel(i+1,j,k-l,3,t+ll) + c11p(bL(3),k-l)*vel(i-1,j,k-l,3,t+ll) 
                end do

                b(6*ll + 1 : 6*(ll+1)) = b(6*ll + 1 : 6*(ll+1)) + &
                                             (/ c11u(bU(1),i)*vel(i+1,j,k,1,t+ll),c11u(bL(1),i-1)*vel(i-2,j,k,1,t+ll),&
                                                c22v(bU(2),j)*vel(i,j+1,k,2,t+ll),c22v(bL(2),j-1)*vel(i,j-2,k,2,t+ll),&
                                                c33w(bU(3),k)*vel(i,j,k+1,3,t+ll),c33w(bL(3),k-1)*vel(i,j,k-2,3,t+ll) /)
            end do

            b = mulL*b
            
            do ll = 0,t_size-1

            ! pressure gradient
            b(6*ll + 1 : 6*(ll+1)) = b(6*ll + 1 : 6*(ll+1)) + &
                    (/ cG1(gU(1),i)*p(i+1,j,k,t+ll),cG1(gL(1),i-1)*p(i-1,j,k,t+ll),&
                       cG2(gU(2),j)*p(i,j+1,k,t+ll),cG2(gL(2),j-1)*p(i,j-1,k,t+ll),&
                       cG3(gU(3),k)*p(i,j,k+1,t+ll),cG3(gL(3),k-1)*p(i,j,k-1,t+ll)/)

            end do
            
            ! time stencil (just in the first time slice)
            b(1:6) = b(1:6) - mulI*(/ vel(i,j,k,1,t-1), vel(i-1,j,k,1,t-1),&
                                        vel(i,j,k,2,t-1), vel(i,j-1,k,2,t-1),&
                                        vel(i,j,k,3,t-1), vel(i,j,k-1,3,t-1) /)
            ! move it to the rhs
            b = - b

            ! add the RHS to b (boundary)
            do ll = 0,t_size-1

                b(6*ll + 1 : 6*(ll+1)) = b(6*ll + 1 : 6*(ll+1)) + &
                        (/ rhs_vel(i,j,k,1,t+ll),rhs_vel(i-1,j,k,1,t+ll),&
                           rhs_vel(i,j,k,2,t+ll),rhs_vel(i,j-1,k,2,t+ll),&
                           rhs_vel(i,j,k,3,t+ll),rhs_vel(i,j,k-1,3,t+ll) /)

                b(block_size - t_size + ll + 1) = b(block_size - t_size + ll + 1) + rhs_p(i,j,k,t+ll)

            end do

! ---------- test the matrix A --------------!


!write(*,*) i,j,k,t
!print *,'-------- the matrix A ----------'

!do l = 1, block_size
!write(*,'(14F7.2)') (A(l,c), c = 1, block_size) 
!end do
!write(*,*)

! -------- Solve the matrix A --------------!'
! subroutine     dgesv (N, NRHS, A, LDA, IPIV, B, LDB, INFO)
 call dgesv( block_size, 1, A, block_size, ipiv, b, block_size, info )
  
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
    write ( *, '(a)' ) '  The matrix is numerically singular.'
    return
  end if

       !assign new solution  
        do ll = 0,t_size-1
            vel(i,j,k,1,t+ll) = b(6*ll+1)
            vel(i-1,j,k,1,t+ll) = b(6*ll+2)
            vel(i,j,k,2,t+ll) = b(6*ll+3)
            vel(i,j-1,k,2,t+ll) = b(6*ll+4)
            vel(i,j,k,3,t+ll) = b(6*ll+5)
            vel(i,j,k-1,3,t+ll) = b(6*ll+6)

            p(i,j,k,t+ll) = b(block_size-t_size+1+ll)
        end do

        end do
      end do
    end do
    end do

    deallocate(A,b,ipiv)

  end subroutine OP_TimeStokesLSmoother

end module cmod_TimeStokesOp
