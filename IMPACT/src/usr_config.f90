!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> \brief configurates all specification by the user first the mandatory ones, then the user-defined at the end
  SUBROUTINE configuration
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  
  !--- all specifications of the configuration ---
  ! note: - mandatory specifications
  !       - additional, user-defined specifications (cf. end of this subroutine)
  !
  
  
  !===========================================================================================================
  !=== job control ===========================================================================================
  !===========================================================================================================
  !--- tasks ---
  ! task = 1 (time integration)
  ! task = 2 (post-processing)
  ! task = 3 (Eigenvalue problems, under construction)
  ! task = 4 (other matrix analysis, under construction)
  !
  task = 1
  
  !--- restart number ---
  restart = 0
  
  
  
  !===========================================================================================================
  !=== resolution ============================================================================================
  !===========================================================================================================
  !--- max. number of time steps (independend of restart count) ---
  ! note: - time integration is terminated when "n_timesteps" or "time_end" is reached
  !
  n_timesteps = 1000000000
  
#ifdef ALLOC
  !--- number of grid points (global) ---
  ! note: - multigrid works most efficiently if Mi = a*2**q+1 and a is a "small" integer
  !       - Mi must be sufficiently large such that at least one central finite difference stencil fits in
  !         (will be improved in future releases)
  !       - for 2D simulations set M3 = 2 (all terms in the third direction are switched off)
  !
  M1 = 257 
  M2 = 65 
!  M3 = 2*2**2+1
  M3 = 129 
  
  !--- numbers of processor blocks ---
  ! note: - multigrid requires that MOD(Mi-1,NBi) = 0
  !       - 2D simulations (M3 = 2) require NB3 = 1
  !
  NB1 = 4
  NB2 = 1
  NB3 = 1
#endif
  
  
  
  !===========================================================================================================
  !=== concentrations / particles ============================================================================
  !===========================================================================================================
  !--- time integration of Euler concentration(s) ---
  concentration_yes = .FALSE.
  
  !--- time integration of Lagrange particles ---
  particles_yes = .FALSE.
  
#ifdef ALLOC
  !--- number of Euler concentrations ---
  n_conc = 1
#endif
  
  !--- total number of Lagrangian particles ---
  n_part_tot = 400000
  
#ifdef ALLOC
  !--- maximum number of Lagrangian particles allowed on each subdomain ---
  n_part_max = 100000
  
  !--- number particle species ---
  n_spec = 1
  
  !--- particle properties ---
  n_args = 14
#endif
  
  !--- particle groups ---
  n_groups = 4
  
  
  
  !===========================================================================================================
  !=== non-dimensional numbers ===============================================================================
  !===========================================================================================================
#ifdef ALLOC
  ALLOCATE(Sc (1:n_conc))
  ALLOCATE(Ric(1:n_conc))
  ALLOCATE(Rip(1:n_spec))
  ALLOCATE(usc(1:n_conc))
  ALLOCATE(usp(1:n_spec))
  ALLOCATE(Stp(1:n_spec))
#endif
  !--- Reynolds number ---
  Re = 200
  !--- frequenzy ---
  freq = 0.2
  
  !--- Froude number ---
#ifdef NONBOUSSINESQ
  nonBoussinesq_yes = .TRUE.
  Fro = 0.5
#endif
  
  !--- Schmidt numbers ---
  Sc(1:n_conc) = 1.
  
  !--- Richardson numbers (concentrations) ---
  ! note: - concentrations are limted to 0 <= conc <= 1
  !
  Ric(1:n_conc) = 1.
  
  !--- Richardson numbers (particles) ---
  Rip(1:n_spec) = 0.25
  
  !--- Stokes settling velocities (concentrations) ---
  ! note: - absolute velocity (scalar), velocity vector is obtained by multiplication with "gravity"
  !
  usc(1:n_conc) = 0.
  
  !--- Stokes settling velocities (particles) ---
  ! note: - absolute velocity (scalar), velocity vector is obtained by multiplication with "gravity"
  !
  usp(1:n_spec) = 0.01
  
  !--- Stokes numbers (particles) ---
  Stp(1:n_spec) = 1.
  
  
  
  !===========================================================================================================
  !=== extents ===============================================================================================
  !===========================================================================================================
  !--- start time ---
  time_start = 0.
  
  !--- end time ---
  ! note: - time integration is terminated when "n_timesteps" or "time_end" is reached
  !
  time_end = 10*1./freq
  periodic_tol = 1.e-11
  
  !--- extents of the physical domain ---
  L1 = 8
  L2 = 2
  L3 = 4
  

  
  
  !===========================================================================================================
  !=== configuration =========================================================================================
  !===========================================================================================================
  !--- velocity boundary conditions (BC) ---
  ! note: - concentration boundary conditions are controlled via "isopycnal" and may depend on the velocity
  !       - for advective boundary conditions you have to specify Dirichlet boundary conditions
  !                           _
  ! symmetry BC:    BC = -2    |
  ! periodicity BC: BC = -1   _|- symmetric, central FD stencils
  ! Dirichlet BC:   BC =  1    |
  ! Neumann BC*:    BC =  2    |- non-symmetric, skew FD stencils
  ! Robin BC*:      BC =  3   _|
  ! *not yet implemented for velocity
  !
  !    orientation of boundary normal
  !    |lower/upper boundary
  !    ||
  ! BC_1L_global
  !
  BC_1L_global =  1
  BC_1U_global =  1
  
  BC_2L_global =  1
  BC_2U_global =  1
  
  BC_3L_global =  1
  BC_3U_global =  1
  
  !--- advective velocity boundaries ---
  ! note: - requires BC_global = 1
  !
  !         orientation of the boundary normal
  !         |   lower/upper boundary
  !         |   |   velocity component
  !         |   |   |
  ! outlet(1:3,1:2,1:3)
  !
  outlet          = .FALSE.
!  outlet(2,2,1:3) = .TRUE.
  
  !--- Dirichlet boundaries for concentrations ---
  ! note: - parameter enforces Dirichlet boundary conditions
  !       - if set to .FALSE.: boundary conditions for the concentrations are set to Robin boundary conditions
  !         if the advection velocity (i.e. fluid + Stokes particle settling velocity) points inside the
  !         domain. Boundary conditions for the concentrations are set to advective boundary conditions if the
  !         advection velocity (i.e. fluid + Stokes particle settling velocity) points outwards the domain
  !
  !            orientation of the boundary normal
  !            |   lower/upper boundary
  !            |   |   concenctration index
  !            |   |   |
  ! isopycnal(1:3,1:2,1:n_conc)
  !
#ifdef ALLOC
  ALLOCATE(isopycnal(1:3,1:2,1:n_conc))
#endif
  isopycnal               = .FALSE.
  isopycnal(2,1,1:n_conc) = .TRUE.
  
  !--- gravity vector ---
  gravity(1:3) = (/-1.,0.,0./)
  
  !--- bulk flow forcing mode ---
  ! huppd: meaning mass-flow control
  ! (under construction, currently limited to channel flows)
  ! note: - fixes bulk flow velocity to 2./3. in periodic channels
  !
  ! forcing_mode = 0 (no forcing)
  ! forcing_mode = 1 (explicit forcing in time)
  ! forcing_mode = 2 (implicit forcing in time)
  !
  forcing_mode = 0
  
  !--- direction of bulk flow (relevant only for bulk flow forcing) ---
  ! (under construction, currently limited to channel flows)
  !
  bulkflow_dir = 2
  
  !--- bulk flow velocity ---
  vel_bulk = 2./3.
  
  
  !===========================================================================================================
  !=== time discretization ===================================================================================
  !===========================================================================================================
  !--- time integration scheme ---
  ! timeint_mode = 0 (CN-RK3)
  ! timeint_mode = 1 (RK3-O3)
  !
  timeint_mode = 0
  
  !--- Crank-Nicolson factor ---
  ! note: - must be 0. < thetaL < 1.
  !
  ! theta = 0. (fully explicit time integration of viscous terms)
  ! theta = 1. (fully implicit time integration of viscous terms)
  !
!  thetaL = 1.0
  thetaL = 0.5
!  thetaL = 0.0
  
  !--- Euler flow (no physical viscosity) ---
  Euler_yes = .FALSE.
  
  !--- Stokes flow (no fluid advection in the domain) ---
  Stokes_yes = .FALSE.
  
  !--- two-step time integration (relevant only for CN-RK3) ---
  twostep_yes = .FALSE.
  
  !--- CFL number ---
  ! note: - normalized with maximum stable number, i.e. should be 0 < CFL < 1
  !       - applies only to advective and viscous terms of velocity and concentrations but not on other user-
  !         defined volume forces
  !       - full RK3 and CN-RK3 stability domains are taken into account
  !       - symbols of the spatial discretization are taken into account 
  !
  CFL = 0.9
  
  !--- max. time step size ---
  dtime_max = 1./freq/320
  
  !--- max. time step size at beginning of time integration ---
  ! note: - "dtime" is limited by "dtime0" for the first "Int_dtime" time steps
  !
  dtime0 = 1./freq/320
  
  !--- number of time steps after which time step size is recomputed ---
  Int_dtime = 900000
  
  
  
  !===========================================================================================================
  !=== spatial discretization ================================================================================
  !===========================================================================================================
  !--- upwinding for advective terms of velocity momentum equation ---
  ! note: - option not available for compact finite differences (cf. comp_conv_yes)
  !
  upwind_yes = .TRUE.
  
  !--- upwinding for advective terms of concentration equations ---
  ! note: - option not available for compact finite differences (cf. comp_conv_yes)
  !
  upwind_conc_yes = .TRUE.
  
  !--- computation of finite difference coefficients (relevant only for explicit finite differences) ---
  ! note: - affects the accuracy and temporal stability of the discretization
  !
  ! mapping_yes = .TRUE.  (computation on computational, equidistant grid and subsequent mapping on
  !                        physical grid)
  ! mapping_yes = .FALSE. (direct computation on physical grid)
  !
  mapping_yes = .TRUE.
  
  !--- compact finite differences ---
  ! note: - currently limited to 10th-order in the inner field
  !
  ! from top: viscous terms, advective terms, interpolations, divergence, gradient
  !
  comp_visc_yes  = .FALSE.
  comp_conv_yes  = .FALSE.
  comp_inter_yes = .FALSE.
  comp_div_yes   = .FALSE.
  comp_grad_yes  = .FALSE.
  
  
  
  !===========================================================================================================
  !=== LES ===================================================================================================
  !===========================================================================================================
#ifdef ALLOC
  ALLOCATE(n_lp_conc(1:n_conc))
  ALLOCATE(n_hp_conc(1:n_conc))
  ALLOCATE(chi_conc (1:n_conc))
#endif
  !--- LES model ---
  ! LES_mode = 0 (no explicit LES model)
  ! LES_mode = 1 (high-order filtering, "ADM-RT")
  ! LES_mode = 2 (Smagorinsky, high-pass filtered Smagorinsky model)
  !
  LES_mode = 0
  
  !--- number of low-pass filter applications (relevant only for LES) ---
  ! from top: velocity, concentrations
  !
  n_lp_vel            = 1
  n_lp_conc(1:n_conc) = 1
  
  !--- number of high-pass filter applications (relevant only for LES) ---
  ! from top: velocity, concentrations
  !
  n_hp_vel            = 6
  n_hp_conc(1:n_conc) = 2
  
  !--- relaxation factor (relevant only for LES) ---
  ! from top: velocity, concentrations
  !
  chi_vel             = 10.
  chi_conc (1:n_conc) = 10.
  
  
  
  !===========================================================================================================
  !=== iterative solver ======================================================================================
  !===========================================================================================================
  !--- infinity norm of velocity divergence ---
  ! note: - cf. comments on "weighting_yes"
  !
  epsU = 10.**(-6)
  
  !--- infinity norm of velocity divergence ---
  ! weighting_yes = .TRUE.  (div(vel) is normalized with the local grid spacing)
  ! weighting_yes = .FALSE. (div(vel) is not normalized)
  !
  weighting_yes = .FALSE.
  
  !--- number of time steps after which the mean pressure is set to zero ---
  ! note: - the mean pressure level is not implicitly fixed in the pressure solver
  !
  Int_lev_pre = 1
  
  !--- preconditioners ---
  ! precond            = 0: no preconditioner, i.e. only primary solvers (Richardson iteration or BiCGstab)
  ! precond_outer      = 1: Richardson iteration with Laplace preconditioner
  ! precond_outer      = 2: Richardson iteration with commuation-based preconditioner
  ! precond_Poisson    = 1: BiCGstab with V-cycle multigrid
  ! precond_Poisson    = 2: BiCGstab with F-cycle multigrid
  ! precond_Helmh_vel  = 1: BiCGstab with V-cycle multigrid
  ! precond_Helmh_conc = 1: BiCGstab with V-cycle multigrid
  !
  ! from top: outer pressure iteration, Helmholtz equation (velocity), Poisson equation,
  !           Helmholtz equations (concentrations)
  !
  precond_outer      = 2
  precond_Poisson    = 1
  precond_Helmh_vel  = 0
  precond_Helmh_conc = 0
  
  !--- max. numbers of iterations ---
  ! from top: outer pressure iteration, Helmholtz equation (velocity), Poisson equation,
  !           Helmholtz equations (concentrations)
  !
  n_it_outer      = 10
  n_it_Poisson    = 10
  n_it_Helmh_vel  = 10
  n_it_Helmh_conc = 10
  
  !--- initialization before implicit problems are solved in a Runge-Kutta sub-step ---
  ! note: - initialization of velocity and concentration applies only to CN-RK3 time integration
  !
  ! from top: pressure, velocity, concentrations
  !
  init_pre (1:RK_steps) = .FALSE.
  init_vel (1:RK_steps) = .FALSE.
  init_conc(1:RK_steps) = .FALSE.
  
  !--- accuracy of preconditioner problems in outer pressure iteration (relevant only for CN-RK3) ---
  ! note: - sets the solution accuracy of the inner Poisson problems relative to the expected accuracy of the
  !         outer pressure iteration cycle
  !
  precOffset0(1:RK_steps) = 0.5
  
  !--- expected convervengence ratios of the first time step of a restart (relevant only for CN-RK3) ---
  precRatio0 (1:RK_steps) = 10.**(-4)
  
  !--- use flux corrections on fine/coarse grids ---
  ! note: - all pressure-Poisson problems are singular (including all multigrid levels)
  !       - to guarantee at least one solution, each right-hand side of a pressure-Poisson problem must be
  !         fully covered by the column space of the corresponding pressure matrix, i.e. each right-hand side
  !         must be orthogonal to the left null space of the corresponding pressure matrix
  !       - setting these variables enforces corrections of the fluxes over the boundaries in each sub-time
  !         step such that the corresponding pressure-Poisson problems have at least one solution
  !       - "nullspace_coarse_yes" (for the multigrid problems) is normally not necessary
  !
  nullspace_yes        = .FALSE.
  nullspace_coarse_yes = .FALSE.
  
  !--- read null space base vectors instead of recomputing them ---
  read_nullspace_yes = .FALSE.
  
  !--- specify if the fluxes are corrected by means of orthogonal projection ---
  ! note: - a vector of the left null space is used as correction vector
  !       - leads to minimal L2-norms of the flux corrections
  !       - fluxes are corrected on the entire boundary (excluding symmetry boundaries)
  !       - if orthogonal projection is not used you need to specify the correction vector "th" in
  !         "usr_initcond.f90"
  !
  nullspace_ortho_yes = .TRUE.
  
  
  
  !===========================================================================================================
  !=== multigrid =============================================================================================
  !===========================================================================================================
  !--- numbers of relaxation sweeps ---
  ! from top: restriction, prolongation, coarsest grid
  !
  n_relax_down   = 4
  n_relax_up     = 4
  n_relax_bottom = 4
  
  !--- directions of (alternating) line-relaxations ---
  ! note: - line-relaxation is necessary if the grid is strongly unisotropic, i.e. when the grid spacings
  !         differ strongly in the different spatial directions
  !
  ! impl_dir = 0 ((partially) red-black-ordered point-block relaxation)
  ! impl_dir = 1 (lexicographically-ordered line-relaxation)
  !
  impl_dir(1:3) = 0
  
  !--- Jacobi smoothing instead of Gauss-Seidel ---
  Jacobi_yes = .FALSE.
  
  !--- max. number of coarse grid levels ---
  ! note: - must be 1 <= n_grids_limit <= n_grids_max
  !
  n_grids_limit = n_grids_max
  
  
  
  !===========================================================================================================
  !=== output ================================================================================================
  !===========================================================================================================
  !--- write standard output ---
  write_stout_yes = .TRUE.
  
  !--- write log file for iterations ---
  log_iteration_yes = .FALSE.
  
  !--- write restart files ---
  write_restart_yes = .TRUE.
  
  !--- time interval of field output ---
  dtime_out_vect = 1./freq

  !--- time interval of other ouput (statistics) ---
  dtime_out_scal = 1.
  
  !--- compute and write lambda2 fields from velocity ---
  write_lambda2_yes = .FALSE.
  
  !--- coarsening factors for 3D fields ---
  ! note: - helpful to save disc space
  !       - 2D data are written without strides
  !       - stride = 0 indicates no output
  !
  stride_large(1:3) = 1
  stride_med  (1:3) = 2
  stride_small(1:3) = 4
  
  !--- write debugging files ---
  ! (deprecated, will vanish in future releases)
  !
  write_test_yes = .FALSE.
  
  
  
  !===========================================================================================================
  !=== additional, user-defined input ========================================================================
  !===========================================================================================================
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  !--- direction of wall normal ---
  wallnorm_dir = 1
  
  !--- max. number of wavenumbers of discrete Fourier transform ---
  amax = 4
  bmax = 4
  
  !--- origin of coordinate system ---
  y1_origin = 0.
  y2_origin = 0.
  y3_origin = 0.
  
  !--- sub-volume for initial particle distribution ---
  L1c = L1
  L2c = 0.25*L1
  L3c = L3
  
  !--- tp ---
  tp = 0
!  velp(:,:,:,:) = vel(:,:,:,:)
!  ALLOCATE(velp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)) !< added


  
  END SUBROUTINE configuration
  
