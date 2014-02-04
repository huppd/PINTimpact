!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> \brief allows user to define starting? boundary conditions for the velocity
  !! Randbedingungen (instantant)
  subroutine boundary_vel_stat
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  integer                ::  i, j, k
  
  
  !--- boundary conditions for velocity ---
  ! note: - advective boundary conditions are automatically initialized from the initial condition
  !       - if constant, bc11 etc. can also be specified in "initial_conditions_vel" in file
  !         "usr_initcond.f90"
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !   velocity component
  !   |orientation of the boundary normal
  !   ||     grid points on the boundary
  !   ||     |         |      lower/upper boundary
  !   ||     |         |      |
  ! bc11(S21B:N21B,S31B:N31B,1:2)
  ! bc12(S11B:N11B,S31B:N31B,1:2)
  ! bc13(S11B:N11B,S21B:N21B,1:2)
  ! bc21(S22B:N22B,S32B:N32B,1:2)
  ! bc22(S12B:N12B,S32B:N32B,1:2)
  ! bc23(S12B:N12B,S22B:N22B,1:2)
  ! bc31(S23B:N23B,S33B:N33B,1:2)
  ! bc32(S13B:N13B,S33B:N33B,1:2)
  ! bc33(S13B:N13B,S23B:N23B,1:2)
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L
  !
  bc11(S21B:N21B,S31B:N31B,1) = 0.
  bc21(S22B:N22B,S32B:N32B,1) = 0.
  bc31(S23B:N23B,S33B:N33B,1) = 0.
  
  bc12(S11B:N11B,S31B:N31B,1) = 0.
  bc32(S13B:N13B,S33B:N33B,1) = 0.
  
  if ((.not. outlet(2,1,2)) .and. BC_2L == 1) then ! Note: this case differentiation may save computational time but is not stringently necessary
     do k = S32B, N32B
        do i = S12B, N12B
!           bc22(i,k,1) = interface((x1p(i)-0.75*L1)/(0.01*L1))
!           bc22(i,k,1) = x2v(i)*( L2 - x2v(i)/L2 )*4/L2
           bc22(i,k,1) = x1p(i)*( L1 - x1p(i) )*4/L1/L1
        end do
     end do
  end if
  
  if ((.not. outlet(2,1,2)) .and. BC_2U == 1) then ! Note: this case differentiation may save computational time but is not stringently necessary
     do k = S32B, N32B
        do i = S12B, N12B
!           bc22(i,k,1) = interface((x1p(i)-0.75*L1)/(0.01*L1))
!           bc22(i,k,1) = x2v(i)*( L2 - x2v(i)/L2 )*4/L2
           bc22(i,k,2) = x1p(i)*( L1 - x1p(i) )*4/L1/L1
        end do
     end do
  end if
  
  end subroutine boundary_vel_stat
  
  
  
  
  
  
  
  
  
  
  
  !> \brief allows user to define time boundary conditions for the velocity
  !! instationaer + Zeitintegration
  subroutine boundary_vel_tint
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  
  !--- additional terms on RHS of time-integrated velocity boundary conditions ---
  ! note: - du/dt = RHS-nlbc
  !       - if not specified, advective boundary conditions are used
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !     velocity component
  !     |orientation of the boundary normal
  !     ||     grid points on the boundary
  !     ||     |         |      lower/upper boundary
  !     ||     |         |      |
  ! nlbc11(S21B:N21B,S31B:N31B,1:2)
  ! nlbc12(S11B:N11B,S31B:N31B,1:2)
  ! nlbc13(S11B:N11B,S21B:N21B,1:2)
  ! nlbc21(S22B:N22B,S32B:N32B,1:2)
  ! nlbc22(S12B:N12B,S32B:N32B,1:2)
  ! nlbc23(S12B:N12B,S22B:N22B,1:2)
  ! nlbc31(S23B:N23B,S33B:N33B,1:2)
  ! nlbc32(S13B:N13B,S33B:N33B,1:2)
  ! nlbc33(S13B:N13B,S23B:N23B,1:2)
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L
  !
  
  
  end subroutine boundary_vel_tint
  
  
  
  
  
  
  
  
  
  
  !> \brief allows user to define boundary conditions for the concentration
  subroutine boundary_conc_stat
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  integer                ::  i, j, k
  
  
  !--- boundary conditions for concentrations ---
  ! note: - boundary conditions for concentration index "conc_nu"
  !       - "res" is overwritten by other routines after each sub-time step
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !     lower/upper boundary
  !     |          grid points on the boundary
  !     |          |       |
  ! res(S1p    ,S2p:N2p,S3p:N3p)
  ! res(N1p    ,S2p:N2p,S3p:N3p)
  ! res(S1p:N1p,S1p    ,S3p:N3p)
  ! res(S1p:N1p,N2p    ,S3p:N3p)
  ! res(S1p:N1p,S2p:N2p,S3p    )
  ! res(S1p:N1p,S2p:N2p,N3p    )
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC:       BC =  3   _|
  ! *not yet implemented for concentrations
  !
  !     orientation of boundary normal
  !     |lower/upper boundary
  !     ||  concenctration index
  !     ||  |
  ! BCc_1L(1:n_conc)
  !
  if (conc_nu == 1) then
     if (BCc_2L(conc_nu) == 1 .and. isopycnal(2,1,conc_nu)) then
        j = 1
        do k = S3p, N3p
           do i = S1p, N1p
              res(i,j,k) = 1.-interface((x1p(i)-0.75*L1)/(0.01*L1))
           end do
        end do
     end if
     
     if (BCc_2U(conc_nu) == 1 .and. (.not. isopycnal(2,2,conc_nu))) then
        j = N2
        do k = S3p, N3p
           do i = S1p, N1p
              ! work is the velocity interpolated on the pressure grid
              if (work2(i,j,k) <= 0.) res(i,j,k) = 1.
           end do
        end do
     end if
  end if
  
  
  if (BCc_1L(conc_nu) == 3) res(1 ,S2p:N2p,S3p:N3p) = 0.
  if (BCc_1U(conc_nu) == 3) res(N1,S2p:N2p,S3p:N3p) = 0.
  
  if (BCc_2L(conc_nu) == 3) res(S1p:N1p,1 ,S3p:N3p) = 0.
  if (BCc_2U(conc_nu) == 3) res(S1p:N1p,N2,S3p:N3p) = 0.
  
  if (BCc_3L(conc_nu) == 3) res(S1p:N1p,S2p:N2p,1 ) = 0.
  if (BCc_3U(conc_nu) == 3) res(S1p:N1p,S2p:N2p,N3) = 0.
  
  
  end subroutine boundary_conc_stat
  
  
  
  
  
  
  
  
  
  
  !> \brief allows user to define boundary conditions for the concentration
  subroutine boundary_conc_tint
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  
  !--- additional terms on RHS of time-integrated concentration boundary conditions ---
  ! note: - dc/dt = RHS-sed1L etc.
  !       - if not specified, Robin or advective boundary conditions are used
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !    orientation of the boundary normal
  !    |lower/upper boundary
  !    ||    grid points on the boundary
  !    ||    |       |     concenctration index
  !    ||    |       |     |
  ! sed1L(S2p:N2p,S3p:N3p,1:n_conc)
  ! sed1U(S2p:N2p,S3p:N3p,1:n_conc)
  ! sed2L(S1p:N1p,S3p:N3p,1:n_conc)
  ! sed2U(S1p:N1p,S3p:N3p,1:n_conc)
  ! sed3L(S1p:N1p,S2p:N2p,1:n_conc)
  ! sed3U(S1p:N1p,S2p:N2p,1:n_conc)
  !
  !
  ! processor-block boundary condition types:
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1    |- symmetric, central FD stencils
  ! neighbor block: BC =  0   _|
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC:       BC =  3   _|
  ! *not yet implemented for concentrations
  !
  !     orientation of boundary normal
  !     |lower/upper boundary
  !     ||  concenctration index
  !     ||  |
  ! BCc_1L(1:n_conc)
  !
  
  
  end subroutine boundary_conc_tint
  
  
  
  
  
  
  
  
  
  
  
  subroutine add_particles
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use mod_particles
  use usr_vars
  use usr_func
  
  implicit none
  
  
  !--- under construction ---
  
  
  end subroutine add_particles
