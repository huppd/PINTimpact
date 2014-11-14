!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing many routines to compute derivatives ...
!! \deprecated Operators get there own module
module cmod_operator
  
  
    !  use mod_dims
    use mod_vars
    use mod_exchange, only: exchange, exchange2

    use iso_c_binding
  
    use mpi
  
    private
  
!    public divergence2
!    public gradient
    public Helmholtz,   Helmholtz_conc_explicit
    public nonlinear
    public interpolate_vel, interpolate_conc
    public outflow_bc, sediment_bc
    public filter
    public bc_extrapolation, bc_extrapolation_transp
    public interpolate_pre_vel , interpolate_vel_pre
    public interpolate2_pre_vel, interpolate2_vel_pre
    public first_pre_vel, first_vel_pre
    public first_adv_pre, first_adv_vel
    public Helmholtz_pre_explicit ! TEST!!!
  
  
contains
  
  
  
    ! TEST!!! Generell bei Randbehandlung bei ii=0,1 beginnen, bzw. bis ii=N1 rechnen!
  


    !  subroutine divergence_transp(m,phi,div)
    !
    !  implicit none
    !
    !  integer, intent(in   ) ::  m
    !
    !  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
    !  real   , intent(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
    !
    !  integer                ::  i, ii
    !  integer                ::  j, jj
    !  integer                ::  k, kk
    !
    !
    !  !===========================================================================================================
    !  if (m == 1) then
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_div_yes) then
    !
    !        do k = S3p, N3p ! Intervall ist ok!
    !           do j = S2p, N2p
    !!pgi$ unroll = n:8
    !              do i = S1p, N1p
    !                 com(i,j,k) = phi(i,j,k)*dx1DM(i)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(1,0,S1p,S2p,S3p,N1p,N2p,N3p,S11B,S21B,S31B,N11B,N21B,N31B,N1,ndL,ndR,dimS1,cDu1CLT,cDu1CLT_LU,cDu1CRT,WDu1T,SDu1T,com,div)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,0,phi)
    !
    !        do k = S31B, N31B ! "B"-Grenzen etwas konsistenter, aber nicht wirklich notwendig, da ohnehin mit bc_extrapolation nachmultipliziert wird ...
    !           do j = S21B, N21B
    !              do i = S11B, N11B
    !                 div(i,j,k) = cDu1T(g1L,i)*phi(i+g1L,j,k)
    !!pgi$ unroll = n:8
    !                 do ii = g1L+1, g1U
    !                    div(i,j,k) = div(i,j,k) + cDu1T(ii,i)*phi(i+ii,j,k)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !  end if
    !  !===========================================================================================================
    !  if (m == 2) then
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_div_yes) then
    !
    !        do k = S3p, N3p
    !           do j = S2p, N2p
    !!pgi$ unroll = n:8
    !              do i = S1p, N1p
    !                 com(i,j,k) = phi(i,j,k)*dx2DM(j)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(2,0,S1p,S2p,S3p,N1p,N2p,N3p,S12B,S22B,S32B,N12B,N22B,N32B,N2,ndL,ndR,dimS2,cDv2CLT,cDv2CLT_LU,cDv2CRT,WDv2T,SDv2T,com,div)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,0,phi)
    !
    !        do k = S32B, N32B
    !           do j = S22B, N22B
    !              do i = S12B, N12B
    !                 div(i,j,k) = cDv2T(g2L,j)*phi(i,j+g2L,k)
    !!pgi$ unroll = n:8
    !                 do jj = g2L+1, g2U
    !                    div(i,j,k) = div(i,j,k) + cDv2T(jj,j)*phi(i,j+jj,k)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !  end if
    !  !===========================================================================================================
    !  if (m == 3) then
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_div_yes) then
    !
    !        do k = S3p, N3p
    !           do j = S2p, N2p
    !!pgi$ unroll = n:8
    !              do i = S1p, N1p
    !                 com(i,j,k) = phi(i,j,k)*dx3DM(k)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(3,0,S1p,S2p,S3p,N1p,N2p,N3p,S13B,S23B,S33B,N13B,N23B,N33B,N3,ndL,ndR,dimS3,cDw3CLT,cDw3CLT_LU,cDw3CRT,WDw3T,SDw3T,com,div)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,0,phi)
    !
    !        do k = S33B, N33B
    !           do j = S23B, N23B
    !              do i = S13B, N13B
    !                 div(i,j,k) = cDw3T(g3L,k)*phi(i,j,k+g3L)
    !!pgi$ unroll = n:8
    !                 do kk = g3L+1, g3U
    !                    div(i,j,k) = div(i,j,k) + cDw3T(kk,k)*phi(i,j,k+kk)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !  end if
    !  !===========================================================================================================
    !
    !
    !  end subroutine divergence_transp





!    !> \brief applies gradient
!    !! to phi, which is a scalar field, and stores the grad in a vector field, the
!    !! boundary condtions for Dirichlet and Neumann are set to zero
!    !! \param[in] m direction in which gradient si calculated
!    !! \param[inout] phi ScalarField from which the gradient is taken
!    !! \param[out] grad gradient in m direction
!    !! \todo add parameterers:
!    !!    - N1,N2,N3,
!    !!    - b1L, b2L, b3L, b1U,b2U,b3U
!    !!    - bC1L, bC2L, bC3L, bC1U,bC2U,bC3U
!    !!    - S11,S21,S31,N11,N21,N31
!    !!    - S12,S22,S32,N12,N22,N32
!    !!    - S13,S23,S33,N13,N23,N33
!    !!    - S11B,S21B,S31B,N11B,N21B,N31B
!    !!    - S12B,S22B,S32B,N12B,N22B,N32B
!    !!    - S13B,S23B,S33B,N13B,N23B,N33B
!    !!    - dimens
!    !!    - cGp1, cGp2, cGp3
!    !!    - g1L,g2L,g3L
!    !!    - g1U,g2U,g3U
!    !! maybe extract setting BC to zero
!    subroutine OP_grad(m,phi,grad) bind(c,name='OP_grad')
!
!        implicit none
!
!        integer(c_int), intent(in   ) ::  m
!
!        real(c_double), intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!        real(c_double), intent(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!
!        integer                ::  i, ii
!        integer                ::  j, jj
!        integer                ::  k, kk
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
!        !                vermutlich nicht wirklich lohnt.                                                          !
!        !----------------------------------------------------------------------------------------------------------!
!        !===========================================================================================================
!        if (m == 1) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S31, N31
!                do j = S21, N21
!                    do i = S11, N11
!                        grad(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
!                        !pgi$ unroll = n:8
!                        do ii = g1L+1, g1U
!                            grad(i,j,k) = grad(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
!                        end do
!                    end do
!                end do
!            end do
!            !--------------------------------------------------------------------------------------------------------
!
!            !--- Randbedingungen ------------------------------------------------------------------------------------
!            if (BC_1L > 0) grad(0 ,S21B:N21B,S31B:N31B) = 0.
!            if (BC_1U > 0) grad(N1,S21B:N21B,S31B:N31B) = 0.
!            if (BC_2L > 0) grad(S11B:N11B,1 ,S31B:N31B) = 0.
!            if (BC_2U > 0) grad(S11B:N11B,N2,S31B:N31B) = 0.
!            if (BC_3L > 0) grad(S11B:N11B,S21B:N21B,1 ) = 0.
!            if (BC_3U > 0) grad(S11B:N11B,S21B:N21B,N3) = 0.
!
!        end if
!        !===========================================================================================================
!        if (m == 2) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S32, N32
!                do j = S22, N22
!                    do i = S12, N12
!                        grad(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
!                        !pgi$ unroll = n:8
!                        do jj = g2L+1, g2U
!                            grad(i,j,k) = grad(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
!                        end do
!                    end do
!                end do
!            end do
!            !--------------------------------------------------------------------------------------------------------
!
!            !--- Randbedingungen ------------------------------------------------------------------------------------
!            if (BC_1L > 0) grad(1 ,S22B:N22B,S32B:N32B) = 0.
!            if (BC_1U > 0) grad(N1,S22B:N22B,S32B:N32B) = 0.
!            if (BC_2L > 0) grad(S12B:N12B,0 ,S32B:N32B) = 0.
!            if (BC_2U > 0) grad(S12B:N12B,N2,S32B:N32B) = 0.
!            if (BC_3L > 0) grad(S12B:N12B,S22B:N22B,1 ) = 0.
!            if (BC_3U > 0) grad(S12B:N12B,S22B:N22B,N3) = 0.
!
!        end if
!        !===========================================================================================================
!        if (m == 3) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S33, N33
!                do j = S23, N23
!                    do i = S13, N13
!                        grad(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
!                        !pgi$ unroll = n:8
!                        do kk = g3L+1, g3U
!                            grad(i,j,k) = grad(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
!                        end do
!                    end do
!                end do
!            end do
!            !--------------------------------------------------------------------------------------------------------
!
!            !--- Randbedingungen ------------------------------------------------------------------------------------
!            if (BC_1L > 0) grad(1 ,S23B:N23B,S33B:N33B) = 0.
!            if (BC_1U > 0) grad(N1,S23B:N23B,S33B:N33B) = 0.
!            if (BC_2L > 0) grad(S13B:N13B,1 ,S33B:N33B) = 0.
!            if (BC_2U > 0) grad(S13B:N13B,N2,S33B:N33B) = 0.
!            if (BC_3L > 0) grad(S13B:N13B,S23B:N23B,0 ) = 0.
!            if (BC_3U > 0) grad(S13B:N13B,S23B:N23B,N3) = 0.
!
!        end if
!    !===========================================================================================================
!
!
!    end subroutine OP_grad


    !  subroutine gradient_transp(m,phi,grad)
    !
    !  implicit none
    !
    !  integer, intent(in   ) ::  m
    !
    !  real   , intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
    !  real   , intent(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
    !
    !  integer                ::  i, ii
    !  integer                ::  j, jj
    !  integer                ::  k, kk
    !
    !
    !  !----------------------------------------------------------------------------------------------------------!
    !  ! Anmerkungen: - Umgekehrte Reihenfolge im Vergleich zu Subroutine "gradient".                             !
    !  !              - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
    !  !                vermutlich nicht wirklich lohnt.                                                          !
    !  !----------------------------------------------------------------------------------------------------------!
    !
    !
    !  !===========================================================================================================
    !  if (m == 1) then
    !
    !     !--- Randbedingungen ------------------------------------------------------------------------------------
    !     if (BC_1L > 0) phi(0 ,S21B:N21B,S31B:N31B) = 0.
    !     if (BC_1U > 0) phi(N1,S21B:N21B,S31B:N31B) = 0.
    !     if (BC_2L > 0) phi(S11B:N11B,1 ,S31B:N31B) = 0.
    !     if (BC_2U > 0) phi(S11B:N11B,N2,S31B:N31B) = 0.
    !     if (BC_3L > 0) phi(S11B:N11B,S21B:N21B,1 ) = 0.
    !     if (BC_3U > 0) phi(S11B:N11B,S21B:N21B,N3) = 0.
    !
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_grad_yes) then
    !
    !        if (BC_1L > 0) com(0 ,S21B:N21B,S31B:N31B) = 0.
    !        if (BC_1U > 0) com(N1,S21B:N21B,S31B:N31B) = 0.
    !        if (BC_2L > 0) com(S11B:N11B,1 ,S31B:N31B) = 0.
    !        if (BC_2U > 0) com(S11B:N11B,N2,S31B:N31B) = 0.
    !        if (BC_3L > 0) com(S11B:N11B,S21B:N21B,1 ) = 0.
    !        if (BC_3U > 0) com(S11B:N11B,S21B:N21B,N3) = 0.
    !
    !        do k = S31, N31 ! Intervall ist ok!
    !           do j = S21, N21
    !!pgi$ unroll = n:8
    !              do i = S11, N11
    !                 com(i,j,k) = phi(i,j,k)*dx1GM(i)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(1,1,S11,S21,S31,N11,N21,N31,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cGp1CLT,cGp1CLT_LU,cGp1CRT,WGp1T,SGp1T,com,grad)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
    !
    !        do k = S3p, N3p
    !           do j = S2p, N2p
    !              do i = S1p, N1p
    !                 grad(i,j,k) = cGp1T(d1L,i)*phi(i+d1L,j,k)
    !!pgi$ unroll = n:8
    !                 do ii = d1L+1, d1U
    !                    grad(i,j,k) = grad(i,j,k) + cGp1T(ii,i)*phi(i+ii,j,k)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !
    !  end if
    !  !===========================================================================================================
    !  if (m == 2) then
    !
    !     !--- Randbedingungen ------------------------------------------------------------------------------------
    !     if (BC_1L > 0) phi(1 ,S22B:N22B,S32B:N32B) = 0.
    !     if (BC_1U > 0) phi(N1,S22B:N22B,S32B:N32B) = 0.
    !     if (BC_2L > 0) phi(S12B:N12B,0 ,S32B:N32B) = 0.
    !     if (BC_2U > 0) phi(S12B:N12B,N2,S32B:N32B) = 0.
    !     if (BC_3L > 0) phi(S12B:N12B,S22B:N22B,1 ) = 0.
    !     if (BC_3U > 0) phi(S12B:N12B,S22B:N22B,N3) = 0.
    !
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_grad_yes) then
    !
    !        if (BC_1L > 0) com(1 ,S22B:N22B,S32B:N32B) = 0.
    !        if (BC_1U > 0) com(N1,S22B:N22B,S32B:N32B) = 0.
    !        if (BC_2L > 0) com(S12B:N12B,0 ,S32B:N32B) = 0.
    !        if (BC_2U > 0) com(S12B:N12B,N2,S32B:N32B) = 0.
    !        if (BC_3L > 0) com(S12B:N12B,S22B:N22B,1 ) = 0.
    !        if (BC_3U > 0) com(S12B:N12B,S22B:N22B,N3) = 0.
    !
    !        do k = S32, N32
    !           do j = S22, N22
    !!pgi$ unroll = n:8
    !              do i = S12, N12
    !                 com(i,j,k) = phi(i,j,k)*dx2GM(j)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(2,2,S12,S22,S32,N12,N22,N32,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cGp2CLT,cGp2CLT_LU,cGp2CRT,WGp2T,SGp2T,com,grad)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
    !
    !        do k = S3p, N3p
    !           do j = S2p, N2p
    !              do i = S1p, N1p
    !!pgi$ unroll = n:8
    !                 do jj = d2L, d2U
    !                    grad(i,j,k) = grad(i,j,k) + cGp2T(jj,j)*phi(i,j+jj,k)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !
    !  end if
    !  !===========================================================================================================
    !  if (m == 3) then
    !
    !     !--- Randbedingungen ------------------------------------------------------------------------------------
    !     if (BC_1L > 0) phi(1 ,S23B:N23B,S33B:N33B) = 0.
    !     if (BC_1U > 0) phi(N1,S23B:N23B,S33B:N33B) = 0.
    !     if (BC_2L > 0) phi(S13B:N13B,1 ,S33B:N33B) = 0.
    !     if (BC_2U > 0) phi(S13B:N13B,N2,S33B:N33B) = 0.
    !     if (BC_3L > 0) phi(S13B:N13B,S23B:N23B,0 ) = 0.
    !     if (BC_3U > 0) phi(S13B:N13B,S23B:N23B,N3) = 0.
    !
    !     !--------------------------------------------------------------------------------------------------------
    !     if (comp_grad_yes) then
    !
    !        if (BC_1L > 0) com(1 ,S23B:N23B,S33B:N33B) = 0.
    !        if (BC_1U > 0) com(N1,S23B:N23B,S33B:N33B) = 0.
    !        if (BC_2L > 0) com(S13B:N13B,1 ,S33B:N33B) = 0.
    !        if (BC_2U > 0) com(S13B:N13B,N2,S33B:N33B) = 0.
    !        if (BC_3L > 0) com(S13B:N13B,S23B:N23B,0 ) = 0.
    !        if (BC_3U > 0) com(S13B:N13B,S23B:N23B,N3) = 0.
    !
    !        do k = S33, N33
    !           do j = S23, N23
    !!pgi$ unroll = n:8
    !              do i = S13, N13
    !                 com(i,j,k) = phi(i,j,k)*dx3GM(k)
    !              end do
    !           end do
    !        end do
    !
    !        call apply_compact_transp(3,3,S13,S23,S33,N13,N23,N33,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cGp3CLT,cGp3CLT_LU,cGp3CRT,WGp3T,SGp3T,com,grad)
    !     !--------------------------------------------------------------------------------------------------------
    !     else
    !        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
    !
    !        do k = S3p, N3p
    !           do j = S2p, N2p
    !              do i = S1p, N1p
    !!pgi$ unroll = n:8
    !                 do kk = d3L, d3U
    !                    grad(i,j,k) = grad(i,j,k) + cGp3T(kk,k)*phi(i,j,k+kk)
    !                 end do
    !              end do
    !           end do
    !        end do
    !     end if
    !     !--------------------------------------------------------------------------------------------------------
    !
    !  end if
    !  !===========================================================================================================
    !
    !
    !  end subroutine gradient_transp



  
  
  
  
!    !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
!    !!
!    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
!    !! \param[in] m dimension from one to three
!    !! \param[in] mulI factor which coresponds to the factor of the identity part
!    !! \param[in] mulL factor which coresponds to the factor of the laplace part
!    !! \param[inout] phi
!    !! \param[out] Lap
!    !! \todo parameters: bl,bu,N, SS,NN,dimens, cu11,cp11,cv22,cp22, cw33, cp33 ...
!    subroutine OP_helmholtz( m, mulI, mulL, phi, Lap ) bind (c,name='OP_helmholtz')
!
!        implicit none
!
!        integer(c_int)  , intent(in   ) ::  m
!        real(c_double)  , intent(in   ) ::  mulI
!        real(c_double)  , intent(in   ) ::  mulL
!
!        real(c_double)  , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!        real(c_double)  , intent(  out) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!
!        integer                ::  i, ii
!        integer                ::  j, jj
!        integer                ::  k, kk
!
!        real                   ::  dd1
!
!
!
!        !===========================================================================================================
!        if (m == 1) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S31, N31
!                    do j = S21, N21
!                        do i = S11, N11
!                            dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 2) then
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do kk = b3L, b3U
!                                dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            else
!                do k = S32, N32
!                    do j = S22, N22
!                        do i = S12, N12
!                            dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                            !pgi$ unroll = n:8
!                            do ii = b1L+1, b1U
!                                dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                            end do
!                            !pgi$ unroll = n:8
!                            do jj = b2L, b2U
!                                dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
!                            end do
!                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                        end do
!                    end do
!                end do
!            end if
!        end if
!           !--------------------------------------------------------------------------------------------------------
!        !===========================================================================================================
!        if (m == 3 .and. dimens == 3) then
!            !--------------------------------------------------------------------------------------------------------
!            do k = S33, N33
!                do j = S23, N23
!                    do i = S13, N13
!                        dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!                        !pgi$ unroll = n:8
!                        do ii = b1L+1, b1U
!                            dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do jj = b2L, b2U
!                            dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
!                        end do
!                        !pgi$ unroll = n:8
!                        do kk = b3L, b3U
!                            dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
!                        end do
!                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
!                    end do
!                end do
!            end do
!        end if
!           !--------------------------------------------------------------------------------------------------------
!    !===========================================================================================================
!
!
!    end subroutine OP_helmholtz
!

  
    !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
    !!
    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
    !! \param[in] m dimension from one to three
    !! \param[in] mulI factor which coresponds to the factor of the identity part
    !! \param[in] mulL factor which coresponds to the factor of the laplace part
    !! \param[inout] phi
    !! \param[out] Lap
    !! \deprecated
    subroutine OP_innerhelmholtz(    &
        m,    &
        SS,   &
        NN,   &
        mulI, &
        mulL, &
        phi,  &
        Lap ) bind (c,name='OP_innerhelmholtz')

        implicit none

        integer(c_int), intent(in ) ::  m

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in   ) ::  mulI
        real(c_double), intent(in   ) ::  mulL

        real(c_double), intent(inout) ::  phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))
        real(c_double), intent(inout) ::  Lap(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        real                   ::  dd1


        !===========================================================================================================
        if (m == 1) then
            if (dimens == 3) then
                do k = SS(3), NN(3)
                    do j = SS(2), NN(2)
                        do i = SS(1), NN(1)
                            dd1 = 0.
                            !pgi$ unroll = n:8
                            do ii = b1L, b1U
                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                            end do
                            !pgi$ unroll = n:8
                            do kk = b3L, b3U
                                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                            end do
                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            else
                do k = SS(3), NN(3)
                    do j = SS(2), NN(2)
                        do i = SS(1), NN(1)
                            dd1 = 0.
                            !pgi$ unroll = n:8
                            do ii = b1L, b1U
                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                            end do
                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            end if
        end if
        !===========================================================================================================
        if (m == 2) then
            if (dimens == 3) then
                do k = SS(3), NN(3)
                    do j = SS(2), NN(2)
                        do i = SS(1), NN(1)
                            dd1 = 0.
                            !pgi$ unroll = n:8
                            do ii = b1L, b1U
                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                            end do
                            !pgi$ unroll = n:8
                            do kk = b3L, b3U
                                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                            end do
                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            else
                do k = SS(3), NN(3)
                    do j = SS(2), NN(2)
                        do i = SS(1), NN(1)
                            dd1 = 0.
                            !pgi$ unroll = n:8
                            do ii = b1L, b1U
                                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                            end do
                            Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            end if
        end if
        !===========================================================================================================
        if (m == 3 .and. dimens == 3) then
            do k = SS(3), NN(3)
                do j = SS(2), NN(2)
                    do i = SS(1), NN(1)
                        dd1 = 0.
                        !pgi$ unroll = n:8
                        do ii = b1L, b1U
                            if( SS(1) <= i+ii .and. i+ii <= NN(1) ) dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = b2L, b2U
                            if( SS(2) <= j+jj .and. j+jj <= NN(2) ) dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = b3L, b3U
                            if( SS(3) <= k+kk .and. k+kk <= NN(3) ) dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
                        end do
                        Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                    end do
                end do
            end do
        end if
        !===========================================================================================================


    end subroutine OP_innerhelmholtz




    !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
    !!
    !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
    !! \deprecated
    subroutine OP_HelmholtzGetRowEntries(    &
        m,             &
        SS,            &
        NN,            &
        mulI,          &
        mulL,          &
        i,             &
        j,             &
        k,             &
        maxRowEntries, &
        ic,            &
        jc,            &
        kc,            &
        val,           &
        rowEntries     &
        ) bind (c,name='OP_HelmholtzGetRowEntries')

        implicit none

        integer(c_int), intent(in ) ::  m

        integer(c_int), intent(in ) ::  i
        integer(c_int), intent(in ) ::  j
        integer(c_int), intent(in ) ::  k

        integer(c_int), intent(in ) ::  SS(3)
        integer(c_int), intent(in ) ::  NN(3)

        integer(c_int), intent(in ) ::  maxRowEntries

        integer(c_int), intent(out) ::  ic(maxRowEntries)
        integer(c_int), intent(out) ::  jc(maxRowEntries)
        integer(c_int), intent(out) ::  kc(maxRowEntries)

        real(c_double), intent(out) ::  val(maxRowEntries)

        integer(c_int), intent(out) ::  rowEntries

        real(c_double), intent(in   ) ::  mulI
        real(c_double), intent(in   ) ::  mulL


        integer                ::  ii
        integer                ::  jj
        integer                ::  kk
        !        integer                ::  counter


        rowEntries = 1;
        !===========================================================================================================
        if (m == 1) then
            if (dimens == 3) then
                val(rowEntries) = mulI - mulL*( cu11(0,i)+cp22(0,j)+cp33(0,k) )
                ic(rowEntries) = i
                jc(rowEntries) = j
                kc(rowEntries) = k
                rowEntries = rowEntries + 1
                !pgi$ unroll = n:8
                do ii = b1L, -1
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cu11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do ii = 1, b1U
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cu11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = b2L, -1
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cp22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = 1, b2U
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cp22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do kk = b3L, -1
                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                        val(rowEntries) = - mulL*( cp33(kk,k) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j
                        kc(rowEntries) = k+kk
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do kk = 1, b3U
                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                        val(rowEntries) = - mulL*( cp33(kk,k) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j
                        kc(rowEntries) = k+kk
                        rowEntries = rowEntries + 1
                    end if
                end do
            else
                val(rowEntries) = mulI - mulL*( cu11(0,i)+cp22(0,j) )
                ic(rowEntries) = i
                jc(rowEntries) = j
                kc(rowEntries) = k
                rowEntries = rowEntries + 1
                !pgi$ unroll = n:8
                do ii = b1L, -1
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cu11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do ii = 1, b1U
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cu11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = b2L, -1
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cp22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = 1, b2U
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cp22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
            end if
        end if
        !===========================================================================================================
        if (m == 2) then
            if (dimens == 3) then
                val(rowEntries) = mulI - mulL*( cp11(0,i)+cv22(0,j)+cp33(0,k) )
                ic(rowEntries) = i
                jc(rowEntries) = j
                kc(rowEntries) = k
                rowEntries = rowEntries + 1
                !pgi$ unroll = n:8
                do ii = b1L, -1
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cp11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do ii = 1, b1U
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cp11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = b2L, -1
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cv22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = 1, b2U
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cv22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do kk = b3L, -1
                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                        val(rowEntries) = - mulL*( cp33(kk,k) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j
                        kc(rowEntries) = k+kk
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do kk = 1, b3U
                    if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                        val(rowEntries) = - mulL*( cp33(kk,k) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j
                        kc(rowEntries) = k+kk
                        rowEntries = rowEntries + 1
                    end if
                end do
            else
                val(rowEntries) = mulI - mulL*( cp11(0,i)+cv22(0,j) )
                ic(rowEntries) = i
                jc(rowEntries) = j
                kc(rowEntries) = k
                rowEntries = rowEntries + 1
                !pgi$ unroll = n:8
                do ii = b1L, -1
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cp11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do ii = 1, b1U
                    if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                        val(rowEntries) = - mulL*( cp11(ii,i) )
                        ic(rowEntries) = i+ii
                        jc(rowEntries) = j
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = b2L, -1
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cv22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
                !pgi$ unroll = n:8
                do jj = 1, b2U
                    if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                        val(rowEntries) = - mulL*( cv22(jj,j) )
                        ic(rowEntries) = i
                        jc(rowEntries) = j+jj
                        kc(rowEntries) = k
                        rowEntries = rowEntries + 1
                    end if
                end do
            end if
        end if
        !===========================================================================================================
        if (m == 3 .and. dimens == 3) then
            val(rowEntries) = mulI - mulL*( cp11(0,i)+cp22(0,j)+cw33(0,k) )
            ic(rowEntries) = i
            jc(rowEntries) = j
            kc(rowEntries) = k
            rowEntries = rowEntries + 1
            !pgi$ unroll = n:8
            do ii = b1L, -1
                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                    val(rowEntries) = - mulL*( cp11(ii,i) )
                    ic(rowEntries) = i+ii
                    jc(rowEntries) = j
                    kc(rowEntries) = k
                    rowEntries = rowEntries + 1
                end if
            end do
            !pgi$ unroll = n:8
            do ii = 1, b1U
                if( SS(1) <= i+ii .and. i+ii <= NN(1) ) then
                    val(rowEntries) = - mulL*( cp11(ii,i) )
                    ic(rowEntries) = i+ii
                    jc(rowEntries) = j
                    kc(rowEntries) = k
                    rowEntries = rowEntries + 1
                end if
            end do
            !pgi$ unroll = n:8
            do jj = b2L, -1
                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                    val(rowEntries) = - mulL*( cp22(jj,j) )
                    ic(rowEntries) = i
                    jc(rowEntries) = j+jj
                    kc(rowEntries) = k
                    rowEntries = rowEntries + 1
                end if
            end do
            !pgi$ unroll = n:8
            do jj = 1, b2U
                if( SS(2) <= j+jj .and. j+jj <= NN(2) ) then
                    val(rowEntries) = - mulL*( cp22(jj,j) )
                    ic(rowEntries) = i
                    jc(rowEntries) = j+jj
                    kc(rowEntries) = k
                    rowEntries = rowEntries + 1
                end if
            end do
            !pgi$ unroll = n:8
            do kk = b3L, -1
                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                    val(rowEntries) = - mulL*( cw33(kk,k) )
                    ic(rowEntries) = i
                    jc(rowEntries) = j
                    kc(rowEntries) = k+kk
                    rowEntries = rowEntries + 1
                end if
            end do
            !pgi$ unroll = n:8
            do kk = 1, b3U
                if( SS(3) <= k+kk .and. k+kk <= NN(3) ) then
                    val(rowEntries) = - mulL*( cw33(kk,k) )
                    ic(rowEntries) = i
                    jc(rowEntries) = j
                    kc(rowEntries) = k+kk
                    rowEntries = rowEntries + 1
                end if
            end do
        end if
        !===========================================================================================================


    end subroutine OP_HelmholtzGetRowEntries



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
    !! \todo shoudl include stencil from PIMP, move to folder src_f
    subroutine OP_DtHelmholtz(   &
        dimens,                 &
        N,                      &
        bL,bU,                  &
        m,                      &
        mulI,                   &
        mulL,                   &
        phiI,                   &
        phiL,                   &
        Lap ) bind (c,name='OP_DtHelmholtz')

        implicit none

        integer(c_int), intent(in)   ::  dimens

        integer(c_int), intent(in)   ::  N(3)

        integer(c_int), intent(in)   ::  bL(3)
        integer(c_int), intent(in)   ::  bU(3)

        integer(c_int), intent(in   )::  m
        real(c_double), intent(in   )::  mulI
        real(c_double), intent(in   )::  mulL

        real(c_double), intent(inout)::  phiI(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double), intent(inout)::  phiL(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(out)  ::  Lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        real                   ::  dd1

        !        write(*,*) cu11(:,4)

        !===========================================================================================================
        if (m == 1) then
            !--------------------------------------------------------------------------------------------------------
            if (dimens == 3) then
                do k = S31, N31
                    do j = S21, N21
                        do i = S11, N11
                            dd1 = cu11(b1L,i)*phiL(i+b1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = b1L+1, b1U
                                dd1 = dd1 + cu11(ii,i)*phiL(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
                            end do
                            !pgi$ unroll = n:8
                            do kk = b3L, b3U
                                dd1 = dd1 + cp33(kk,k)*phiL(i,j,k+kk)
                            end do
                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            else
                do k = S31, N31
                    do j = S21, N21
                        do i = S11, N11
                            dd1 = cu11(b1L,i)*phiL(i+b1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = b1L+1, b1U
                                dd1 = dd1 + cu11(ii,i)*phiL(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
                            end do
                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            end if
        end if
           !--------------------------------------------------------------------------------------------------------
        !===========================================================================================================
        if (m == 2) then
            !--------------------------------------------------------------------------------------------------------
            if (dimens == 3) then
                do k = S32, N32
                    do j = S22, N22
                        do i = S12, N12
                            dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = b1L+1, b1U
                                dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                dd1 = dd1 + cv22(jj,j)*phiL(i,j+jj,k)
                            end do
                            !pgi$ unroll = n:8
                            do kk = b3L, b3U
                                dd1 = dd1 + cp33(kk,k)*phiL(i,j,k+kk)
                            end do
                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            else
                do k = S32, N32
                    do j = S22, N22
                        do i = S12, N12
                            dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = b1L+1, b1U
                                dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
                            end do
                            !pgi$ unroll = n:8
                            do jj = b2L, b2U
                                dd1 = dd1 + cv22(jj,j)*phiL(i,j+jj,k)
                            end do
                            Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                        end do
                    end do
                end do
            end if
        end if
           !--------------------------------------------------------------------------------------------------------
        !===========================================================================================================
        if (m == 3 .and. dimens == 3) then
            !--------------------------------------------------------------------------------------------------------
            do k = S33, N33
                do j = S23, N23
                    do i = S13, N13
                        dd1 = cp11(b1L,i)*phiL(i+b1L,j,k)
                        !pgi$ unroll = n:8
                        do ii = b1L+1, b1U
                            dd1 = dd1 + cp11(ii,i)*phiL(i+ii,j,k)
                        end do
                        !pgi$ unroll = n:8
                        do jj = b2L, b2U
                            dd1 = dd1 + cp22(jj,j)*phiL(i,j+jj,k)
                        end do
                        !pgi$ unroll = n:8
                        do kk = b3L, b3U
                            dd1 = dd1 + cw33(kk,k)*phiL(i,j,k+kk)
                        end do
                        Lap(i,j,k) = mulI*phiI(i,j,k) - mulL*dd1
                    end do
                end do
            end do
        end if
           !--------------------------------------------------------------------------------------------------------
    !===========================================================================================================


    end subroutine OP_DtHelmholtz












    !> \brief computes nonlinear terms.
    !!
    !! first interpolates worki to pp, then \f$\partial_i\f$ \c vel then
    !! \f[ \mathrm{nl = pp*\partial_i vel}\f]
    !! \test Teile davon (nur zentrale Operationen!) koennten jeweils durch interpolate2_pre_vel/interpolate2_vel_pre, first_adv_pre/first_adv_vel
    !!         ersetzt werden (beachte aber Addition von nl!)
    !! \test umbenennen in convection ... (?)
    !! \relates Pimpact::Nonlinear
    !! \deprecated
    subroutine OP_nonlinear(   &
        phi1U,phi1V,phi1W,  &
        phi2U,phi2V,phi2W,  &
        nlU,nlV,nlW,        &
        mul ) bind (c,name='OP_nonlinear')

        implicit none


        real(c_double),  intent(inout) ::  phi1U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  phi1V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  phi1W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double),  intent(inout) ::  phi2U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  phi2V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  phi2W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double),  intent(inout) ::  nlU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  nlV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double),  intent(inout) ::  nlW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double),  intent(in)    :: mul

        real                   ::  dd1

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk


        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]                                         !
        !              - Feld "res" wird mehrfach ausgetauscht, um mit einem statt drei skalaren Feldern arbeiten  !
        !                zu können. Im Prinzip könnte auch rhs(:,:,:,1:3) für die Zwischenspeicherung verwendet    !
        !                werden, jedoch sind dazu einige Umbaumassnahmen notwendig bei geringem Effizienzgewinn.   !
        !----------------------------------------------------------------------------------------------------------!

        call interpolate2_vel_pre(1,phi1U(b1L,b2L,b3L),work1)
        call interpolate2_vel_pre(2,phi1V(b1L,b2L,b3L),work2)
        if (dimens == 3) call interpolate2_vel_pre(3,phi1W(b1L,b2L,b3L),work3)


        call exchange(1,0,work1)
        call exchange(2,0,work1)
        call exchange(3,0,work1)

        call exchange(1,0,work2)
        call exchange(2,0,work2)
        call exchange(3,0,work2)

        if (dimens == 3) then
            call exchange(1,0,work3)
            call exchange(2,0,work3)
            call exchange(3,0,work3)
        end if


        !===========================================================================================================
        !=== u*du/dx ===============================================================================================
        !===========================================================================================================
        call interpolate2_pre_vel(1,work1,pp) ! work1 -> pp, why?
        !-----------------------------------------------------------------------------------------------------------
        if (upwind_yes) then
            do k = S31, N31
                do j = S21, N21
                    do i = S11, N11
                        if (pp(i,j,k) >= 0.) then
                            dd1 = cNu1U(n1L,i)*phi2U(i+n1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = n1L+1, n1U
                                dd1 = dd1 + cNu1U(ii,i)*phi2U(i+ii,j,k)
                            end do
                        else
                            dd1 = cNu1D(n1L,i)*phi2U(i+n1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = n1L+1, n1U
                                dd1 = dd1 + cNu1D(ii,i)*phi2U(i+ii,j,k)
                            end do
                        end if

                        nlU(i,j,k) = nlU(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        else
            do k = S31, N31
                do j = S21, N21
                    do i = S11, N11
                        dd1 = cu1(b1L,i)*phi2U(i+b1L,j,k)
                        !pgi$ unroll = n:8
                        do ii = b1L+1, b1U
                            dd1 = dd1 + cu1(ii,i)*phi2U(i+ii,j,k)
                        end do

                        nlU(i,j,k) = nlU(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        end if
        !===========================================================================================================


        !===========================================================================================================
        !=== v*du/dy ===============================================================================================
        !===========================================================================================================
        call interpolate2_pre_vel(1,work2,pp)
        !-----------------------------------------------------------------------------------------------------------
        if (upwind_yes) then
            do k = S31, N31
                do j = S21, N21
                    do i = S11, N11
                        if (pp(i,j,k) >= 0.) then
                            dd1 = cNp2U(n2L,j)*phi2U(i,j+n2L,k)
                            !pgi$ unroll = n:8
                            do jj = n2L+1, n2U
                                dd1 = dd1 + cNp2U(jj,j)*phi2U(i,j+jj,k)
                            end do
                        else
                            dd1 = cNp2D(n2L,j)*phi2U(i,j+n2L,k)
                            !pgi$ unroll = n:8
                            do jj = n2L+1, n2U
                                dd1 = dd1 + cNp2D(jj,j)*phi2U(i,j+jj,k)
                            end do
                        end if

                        nlU(i,j,k) =  nlU(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        else
            do k = S31, N31
                do j = S21, N21
                    do i = S11, N11
                        dd1 = cp2(b2L,j)*phi2U(i,j+b2L,k)
                        !pgi$ unroll = n:8
                        do jj = b2L+1, b2U
                            dd1 = dd1 + cp2(jj,j)*phi2U(i,j+jj,k)
                        end do

                        nlU(i,j,k) = nlU(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        end if
        !===========================================================================================================


        if (dimens == 3) then
            !===========================================================================================================
            !=== w*du/dz ===============================================================================================
            !===========================================================================================================
            call interpolate2_pre_vel(1,work3,pp)
            !-----------------------------------------------------------------------------------------------------------
            if (upwind_yes) then
                do k = S31, N31
                    do j = S21, N21
                        do i = S11, N11
                            if (pp(i,j,k) >= 0.) then
                                dd1 = cNp3U(n3L,k)*phi2U(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNp3U(kk,k)*phi2U(i,j,k+kk)
                                end do
                            else
                                dd1 = cNp3D(n3L,k)*phi2U(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNp3D(kk,k)*phi2U(i,j,k+kk)
                                end do
                            end if

                            nlU(i,j,k) =  nlU(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            else
                do k = S31, N31
                    do j = S21, N21
                        do i = S11, N11
                            dd1 = cp3(b3L,k)*phi2U(i,j,k+b3L)
                            !pgi$ unroll = n:8
                            do kk = b3L+1, b3U
                                dd1 = dd1 + cp3(kk,k)*phi2U(i,j,k+kk)
                            end do

                            nlU(i,j,k) =  nlU(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            end if
        !===========================================================================================================
        end if




        !===========================================================================================================
        !=== u*dv/dx ===============================================================================================
        !===========================================================================================================
        call interpolate2_pre_vel(2,work1,pp)
        !-----------------------------------------------------------------------------------------------------------
        if (upwind_yes) then
            do k = S32, N32
                do j = S22, N22
                    do i = S12, N12
                        if (pp(i,j,k) >= 0.) then
                            dd1 = cNp1U(n1L,i)*phi2V(i+n1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = n1L+1, n1U
                                dd1 = dd1 + cNp1U(ii,i)*phi2V(i+ii,j,k)
                            end do
                        else
                            dd1 = cNp1D(n1L,i)*phi2V(i+n1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = n1L+1, n1U
                                dd1 = dd1 + cNp1D(ii,i)*phi2V(i+ii,j,k)
                            end do
                        end if

                        nlV(i,j,k) = nlV(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        else
            do k = S32, N32
                do j = S22, N22
                    do i = S12, N12
                        dd1 = cp1(b1L,i)*phi2V(i+b1L,j,k)
                        !pgi$ unroll = n:8
                        do ii = b1L+1, b1U
                            dd1 = dd1 + cp1(ii,i)*phi2V(i+ii,j,k)
                        end do

                        nlV(i,j,k) = nlV(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        end if
        !===========================================================================================================


        !===========================================================================================================
        !=== v*dv/dy ===============================================================================================
        !===========================================================================================================
        call interpolate2_pre_vel(2,work2,pp)
        !-----------------------------------------------------------------------------------------------------------
        if (upwind_yes) then
            do k = S32, N32
                do j = S22, N22
                    do i = S12, N12
                        if (pp(i,j,k) >= 0.) then
                            dd1 = cNv2U(n2L,j)*phi2V(i,j+n2L,k)
                            !pgi$ unroll = n:8
                            do jj = n2L+1, n2U
                                dd1 = dd1 + cNv2U(jj,j)*phi2V(i,j+jj,k)
                            end do
                        else
                            dd1 = cNv2D(n2L,j)*phi2V(i,j+n2L,k)
                            !pgi$ unroll = n:8
                            do jj = n2L+1, n2U
                                dd1 = dd1 + cNv2D(jj,j)*phi2V(i,j+jj,k)
                            end do
                        end if

                        nlV(i,j,k) =  nlV(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        else
            do k = S32, N32
                do j = S22, N22
                    do i = S12, N12
                        dd1 = cv2(b2L,j)*phi2V(i,j+b2L,k)
                        !pgi$ unroll = n:8
                        do jj = b2L+1, b2U
                            dd1 = dd1 + cv2(jj,j)*phi2V(i,j+jj,k)
                        end do

                        nlV(i,j,k) =  nlV(i,j,k) + mul*dd1*pp(i,j,k)
                    end do
                end do
            end do
        end if
        !===========================================================================================================


        if (dimens == 3) then
            !===========================================================================================================
            !=== w*dv/dz ===============================================================================================
            !===========================================================================================================
            call interpolate2_pre_vel(2,work3,pp)
            !-----------------------------------------------------------------------------------------------------------
            if (upwind_yes) then
                do k = S32, N32
                    do j = S22, N22
                        do i = S12, N12
                            if (pp(i,j,k) >= 0.) then
                                dd1 = cNp3U(n3L,k)*phi2V(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNp3U(kk,k)*phi2V(i,j,k+kk)
                                end do
                            else
                                dd1 = cNp3D(n3L,k)*phi2V(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNp3D(kk,k)*phi2V(i,j,k+kk)
                                end do
                            end if

                            nlV(i,j,k) =  nlV(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            else
                do k = S32, N32
                    do j = S22, N22
                        do i = S12, N12
                            dd1 = cp3(b3L,k)*phi2V(i,j,k+b3L)
                            !pgi$ unroll = n:8
                            do kk = b3L+1, b3U
                                dd1 = dd1 + cp3(kk,k)*phi2V(i,j,k+kk)
                            end do

                            nlV(i,j,k) =  nlV(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            end if
        !===========================================================================================================
        end if





        if (dimens == 3) then
            !===========================================================================================================
            !=== u*dw/dx ===============================================================================================
            !===========================================================================================================
            call interpolate2_pre_vel(3,work1,pp)
            !-----------------------------------------------------------------------------------------------------------
            if (upwind_yes) then
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            if (pp(i,j,k) >= 0.) then
                                dd1 = cNp1U(n1L,i)*phi2W(i+n1L,j,k)
                                !pgi$ unroll = n:8
                                do ii = n1L+1, n1U
                                    dd1 = dd1 + cNp1U(ii,i)*phi2W(i+ii,j,k)
                                end do
                            else
                                dd1 = cNp1D(n1L,i)*phi2W(i+n1L,j,k)
                                !pgi$ unroll = n:8
                                do ii = n1L+1, n1U
                                    dd1 = dd1 + cNp1D(ii,i)*phi2W(i+ii,j,k)
                                end do
                            end if

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            else
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            dd1 = cp1(b1L,i)*phi2W(i+b1L,j,k)
                            !pgi$ unroll = n:8
                            do ii = b1L+1, b1U
                                dd1 = dd1 + cp1(ii,i)*phi2W(i+ii,j,k)
                            end do

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            end if
            !===========================================================================================================


            !===========================================================================================================
            !=== v*dw/dy ===============================================================================================
            !===========================================================================================================
            call interpolate2_pre_vel(3,work2,pp)
            !-----------------------------------------------------------------------------------------------------------
            if (upwind_yes) then
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            if (pp(i,j,k) >= 0.) then
                                dd1 = cNp2U(n2L,j)*phi2W(i,j+n2L,k)
                                !pgi$ unroll = n:8
                                do jj = n2L+1, n2U
                                    dd1 = dd1 + cNp2U(jj,j)*phi2W(i,j+jj,k)
                                end do
                            else
                                dd1 = cNp2D(n2L,j)*phi2W(i,j+n2L,k)
                                !pgi$ unroll = n:8
                                do jj = n2L+1, n2U
                                    dd1 = dd1 + cNp2D(jj,j)*phi2W(i,j+jj,k)
                                end do
                            end if

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            else
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            dd1 = cp2(b2L,j)*phi2W(i,j+b2L,k)
                            !pgi$ unroll = n:8
                            do jj = b2L+1, b2U
                                dd1 = dd1 + cp2(jj,j)*phi2W(i,j+jj,k)
                            end do

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            end if
            !===========================================================================================================
  
  
            !===========================================================================================================
            !=== w*dw/dz ===============================================================================================
            !===========================================================================================================
            call interpolate2_pre_vel(3,work3,pp)
            !-----------------------------------------------------------------------------------------------------------
            if (upwind_yes) then
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            if (pp(i,j,k) >= 0.) then
                                dd1 = cNw3U(n3L,k)*phi2W(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNw3U(kk,k)*phi2W(i,j,k+kk)
                                end do
                            else
                                dd1 = cNw3D(n3L,k)*phi2W(i,j,k+n3L)
                                !pgi$ unroll = n:8
                                do kk = n3L+1, n3U
                                    dd1 = dd1 + cNw3D(kk,k)*phi2W(i,j,k+kk)
                                end do
                            end if

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            else
                do k = S33, N33
                    do j = S23, N23
                        do i = S13, N13
                            dd1 = cw3(b3L,k)*phi2W(i,j,k+b3L)
                            !pgi$ unroll = n:8
                            do kk = b3L+1, b3U
                                dd1 = dd1 + cw3(kk,k)*phi2W(i,j,k+kk)
                            end do

                            nlW(i,j,k) = nlW(i,j,k) + mul*dd1*pp(i,j,k)
                        end do
                    end do
                end do
            end if
        !===========================================================================================================
        end if


    end subroutine OP_nonlinear








    !> \brief interpolates pressure to vel nodes
    !!
    !! \f[ inter = inter + interpolated(phi) \f]
    !! Wie interpolate_pre_vel, allerdings mit fixen Index-Limiten (ohne Rand)
    !! \param[in] m dimension
    !! \param[inout] phi input field
    !! \param[inout] inter output field
    !! \deprecated
    subroutine interpolate2_pre_vel(m,phi,inter)

        implicit none

        !  logical, intent(in   ) ::  exch_yes
        integer, intent(in   ) ::  m

        real   , intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
  
        !  if (exch_yes) call exchange(m,0,phi)
  
  
        !===========================================================================================================
        if (m == 1) then
            do k = S31, N31
                do j = S21, N21
                    do i = S11, N11
                        inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                        do ii = g1L+1, g1U
                            inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                        end do
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 2) then
            do k = S32, N32
                do j = S22, N22
                    do i = S12, N12
                        inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                        do jj = g2L+1, g2U
                            inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                        end do
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 3) then
            do k = S33, N33
                do j = S23, N23
                    do i = S13, N13
                        inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                        do kk = g3L+1, g3U
                            inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                        end do
                    end do
                end do
            end do
        end if
    !===========================================================================================================


    end subroutine interpolate2_pre_vel

  
  
  
  
  


    !> \brief interpolates velocity to pressure grid
    !!
    !! Wie interpolate_vel_pre, allerdings mit fixen Index-Limiten (ohne Rand)
    !! \deprecated
    subroutine interpolate2_vel_pre(m,phi,inter)

        implicit none

        !  logical, intent(in   ) ::  exch_yes
        integer, intent(in)    ::  m

        real   , intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk


        !  if (exch_yes) call exchange(m,m,phi)


        !===========================================================================================================
        if (m == 1) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                        do ii = d1L+1, d1U
                            inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                        end do
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 2) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                        do jj = d2L+1, d2U
                            inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                        end do
                    end do
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 3) then
            do k = S3p, N3p
                do j = S2p, N2p
                    do i = S1p, N1p
                        inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                        do kk = d3L+1, d3U
                            inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                        end do
                    end do
                end do
            end do
        end if
    !===========================================================================================================


    end subroutine interpolate2_vel_pre




    !> \brief extrapolates Direchlet-BC to pressure
    !! \param[in] m direction
    !! \param[inout] phi pressure
    !! \attention Wird in "product_div_grad" und "explicit" verwendet almoste ver after grad
    subroutine bc_extrapolation( m, phi ) bind (c,name='OP_bc_extrapolation')

        implicit none

        integer(c_int), intent(in   ) ::  m
        real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk




        ! TEST!!! bislang nur Dirichlet-RB eingebaut!
        !-----------------------------------------------------------------------------------------------------------
        if (m == 1) then
            if (BC_1L > 0) then
                i = 0
                do k = S31B, N31B
                    do j = S21B, N21B
                        !pgi$ unroll = n:8
                        do ii = 0, d1U
                            phi(i,j,k) = phi(i,j,k) - cIup(ii,1)*phi(1+ii,j,k)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
                    end do
                end do
            end if
            if (BC_1U > 0) then
                i = N1
                do k = S31B, N31B
                    do j = S21B, N21B
                        !pgi$ unroll = n:8
                        do ii = d1L, -1
                            phi(i,j,k) = phi(i,j,k) - cIup(ii,i)*phi(i+ii,j,k)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIup(0,i)
                    end do
                end do
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 2) then
            if (BC_2L > 0) then
                j = 0
                do k = S32B, N32B
                    do i = S12B, N12B
                        !pgi$ unroll = n:8
                        do jj = 0, d2U
                            phi(i,j,k) = phi(i,j,k) - cIvp(jj,1)*phi(i,1+jj,k)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
                    end do
                end do
            end if
            if (BC_2U > 0) then
                j = N2
                do k = S32B, N32B
                    do i = S12B, N12B
                        !pgi$ unroll = n:8
                        do jj = d2L, -1
                            phi(i,j,k) = phi(i,j,k) - cIvp(jj,j)*phi(i,j+jj,k)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
                    end do
                end do
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 3) then
            if (BC_3L > 0) then
                k = 0
                do j = S23B, N23B
                    do i = S13B, N13B
                        !pgi$ unroll = n:8
                        do kk = 0, d3U
                            phi(i,j,k) = phi(i,j,k) - cIwp(kk,1)*phi(i,j,1+kk)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
                    end do
                end do
            end if
            if (BC_3U > 0) then
                k = N3
                do j = S23B, N23B
                    do i = S13B, N13B
                        !pgi$ unroll = n:8
                        do kk = d3L, -1
                            phi(i,j,k) = phi(i,j,k) - cIwp(kk,k)*phi(i,j,k+kk)
                        end do
                        phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
                    end do
                end do
            end if
        end if
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine bc_extrapolation
  
  
  
  
  
  
  
  
  
  
  
!  subroutine bc_extrapolation_transp(m,phi)
!
!  implicit none
!
!  integer, intent(in   ) ::  m
!  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!
!  integer                ::  i, ii
!  integer                ::  j, jj
!  integer                ::  k, kk
!
!
!  !----------------------------------------------------------------------------------------------------------!
!  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
!  !----------------------------------------------------------------------------------------------------------!
!
!
!  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
!  !-----------------------------------------------------------------------------------------------------------
!  if (m == 1) then
!     if (BC_1L > 0) then
!        i = 0
!        do k = S31B, N31B
!           do j = S21B, N21B
!!pgi$ unroll = n:8
!              do ii = 0, d1U
!                 phi(i+ii+1,j,k) = phi(i+ii+1,j,k) - phi(i,j,k)*cIup(ii,1)/cIup(-1,1)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
!           end do
!        end do
!     end if
!     if (BC_1U > 0) then
!        i = N1
!        do k = S31B, N31B
!           do j = S21B, N21B
!!pgi$ unroll = n:8
!              do ii = d1L, -1
!                 phi(i+ii,j,k) = phi(i+ii,j,k) - phi(i,j,k)*cIup(ii,i)/cIup(0,i)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
!           end do
!        end do
!     end if
!  end if
!  !-----------------------------------------------------------------------------------------------------------
!  if (m == 2) then
!     if (BC_2L > 0) then
!        j = 0
!        do k = S32B, N32B
!           do i = S12B, N12B
!!pgi$ unroll = n:8
!              do jj = 0, d2U
!                 phi(i,j+jj+1,k) = phi(i,j+jj+1,k) - phi(i,j,k)*cIvp(jj,1)/cIvp(-1,1)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
!           end do
!        end do
!     end if
!     if (BC_2U > 0) then
!        j = N2
!        do k = S32B, N32B
!           do i = S12B, N12B
!!pgi$ unroll = n:8
!              do jj = d2L, -1
!                 phi(i,j+jj,k) = phi(i,j+jj,k) - phi(i,j,k)*cIvp(jj,j)/cIvp(0,j)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
!           end do
!        end do
!     end if
!  end if
!  !-----------------------------------------------------------------------------------------------------------
!  if (m == 3) then
!     if (BC_3L > 0) then
!        k = 0
!        do j = S23B, N23B
!           do i = S13B, N13B
!!pgi$ unroll = n:8
!              do kk = 0, d3U
!                 phi(i,j,k+kk+1) = phi(i,j,k+kk+1) - phi(i,j,k)*cIwp(kk,1)/cIwp(-1,1)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
!           end do
!        end do
!     end if
!     if (BC_3U > 0) then
!        k = N3
!        do j = S23B, N23B
!           do i = S13B, N13B
!!pgi$ unroll = n:8
!              do kk = d3L, -1
!                 phi(i,j,k+kk) = phi(i,j,k+kk) - phi(i,j,k)*cIwp(kk,k)/cIwp(0,k)
!              end do
!              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
!           end do
!        end do
!     end if
!  end if
!  !-----------------------------------------------------------------------------------------------------------
!
!
!  end subroutine bc_extrapolation_transp
  
  
  
end module cmod_operator
