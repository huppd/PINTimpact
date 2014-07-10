!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
!> \brief user defines initial conditions for the velocity, the connection to the boundary conditions is
!! the commenter not yet clear
SUBROUTINE initial_conditions_vel
    ! (basic subroutine)
  
    USE mod_dims
    USE mod_vars
    USE usr_vars
    USE usr_func
  
    IMPLICIT NONE
  
    INTEGER                ::  i, j, k
  
  
    !--- initial conditions for pressure ---
    ! note: - initialization is generally not necessary
    !       - specification provides only an initial guess for the first pressure iteration
    !       - cf. sketch in file "usr_geometry.f90"
    !
    !        grid points in the domain and on the boundary
    !        |       |       |
    ! pre(S1p:N1p,S2p:N2p,S3p:N3p)
    !
    pre = 0.
  
  
    !--- initial conditions for velocity ---
    ! note: - cf. sketch in file "usr_geometry.f90"
    !
    !         grid points in the domain and on the boundary
    !         |         |         |     velocity component
    !         |         |         |     |
    ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
    ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
    ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
    !
    vel = 0.
!    DO k = S32B, N32B
!        DO j = S22B, N22B
!            DO i = S12B, N12B
!                !vel(i,j,k,2) = vel(i,j,k,2) + interface((x1p(i)-0.75*L1)/(0.01*L1))
!                vel(i,j,k,2) = x1p(i)*( L1 - x1p(i) )*4/L1/L1
!            !           vel(i,j,k,2) = 1.0
!            END DO
!        END DO
!    END DO
    DO k = S31B, N31B
        DO j = S21B, N21B
            DO i = S11B, N11B
                vel(i,j,k,1) = x2p(j)*( L2 - x2p(j) )*4/L2/L2
            END DO
        END DO
    END DO
  
  
    !--- advection velocity ---
    ! note: - only boundary-normal velocity component is specified, i.e. velocity
    !         components tangential to the boundary are not considered
    !       - must be specified for the pressure grid points (values for advection
    !         velocity of boundary-tangential velocity components are interpolated)
    !       - cf. sketch in file "usr_geometry.f90"
    !
    !      orientation of the boundary normal
    !      |    grid points on the boundary
    !      |    |       |     lower/upper boundary
    !      |    |       |     |
    ! drift1(S2p:N2p,S3p:N3p,1:2) = 1.
    ! drift2(S1p:N1p,S3p:N3p,1:2) = 1.
    ! drift3(S1p:N1p,S2p:N2p,1:2) = 1.
    !
    drift1 = 1.
    drift2 = 1.
    drift3 = 1.
  
  
    !--- specify vector for flux corrections ---
    ! note: - needs to be specified only if orthogonal projection is NOT used, i.e. for
    !         "nullspace_ortho_yes = .FALSE."
    !       - correction can be applied only at boundaries (not symmetry)
    !       - correction vector must not be orthogonal to the null space vector "psi_vel"
    !       - only the shape matters, not the amplitude or the sign
    !       - the correction vector will be normalized subsequently
    !       - the correction vector is applied only in sub-domains located at the corresponding boundaries
    !       - cf. sketch in file "usr_geometry.f90"
    !
    !   velocity component
    !   |orientation of the boundary normal
    !   ||     grid points on the boundary
    !   ||     |         |      lower/upper boundary
    !   ||     |         |      |
    ! th11(S21B:N21B,S31B:N31B,1:2)
    ! th12(S11B:N11B,S31B:N31B,1:2)
    ! th13(S11B:N11B,S21B:N21B,1:2)
    ! th21(S22B:N22B,S32B:N32B,1:2)
    ! th22(S12B:N12B,S32B:N32B,1:2)
    ! th23(S12B:N12B,S22B:N22B,1:2)
    ! th31(S23B:N23B,S33B:N33B,1:2)
    ! th32(S13B:N13B,S33B:N33B,1:2)
    ! th33(S13B:N13B,S23B:N23B,1:2)
    !
    th11 = 1.
    th12 = 0.
    th13 = 0.
  
    th21 = 0.
    th22 = 1.
    th23 = 0.
  
    th31 = 0.
    th32 = 0.
    th33 = 1.
  
  
    !--- additional, user-defined initializations ---
    !
    ! for statistics:
    energy_visc       = 0.
    diss_viscInt_old  = 0.
  
  
END SUBROUTINE initial_conditions_vel
  
  
  
  
  
  
  
  
  
  
  
SUBROUTINE initial_conditions_conc
    ! (basic subroutine)
  
    USE mod_dims
    USE mod_vars
    USE usr_vars
    USE usr_func
  
    IMPLICIT NONE
  
    INTEGER                ::  i, j, k
  
  
    !--- initial conditions for concentrations ---
    ! note: - cf. sketch in file "usr_geometry.f90"
    !
    !         grid points in the domain and on the boundary
    !         |       |       |     concentration index
    !         |       |       |     |
    ! conc(S1p:N1p,S2p:N2p,S3p:N3p,1:n_conc)
    !
    conc = 0.
    DO k = S3p, N3p
        DO j = S2p, N2p
            DO i = S1p, N1p
                conc(i,j,k,1) = 1. - interface((x1p(i)-0.75*L1)/(0.01*L1))
            END DO
        END DO
    END DO
  
  
    !--- initial conditions for concentration deposit ---
    ! note: - cf. sketch in file "usr_geometry.f90"
    !       - I/O is implemented currently only for "dep1L_conc"
    !
    !    orientation of the boundary normal
    !    |lower/upper boundary
    !    ||         grid points on the boundary
    !    ||         |       |     concenctration index
    !    ||         |       |     |
    ! dep1L_conc(S2p:N2p,S3p:N3p,1:n_conc)
    ! dep1U_conc(S2p:N2p,S3p:N3p,1:n_conc)
    ! dep2L_conc(S1p:N1p,S3p:N3p,1:n_conc)
    ! dep2U_conc(S1p:N1p,S3p:N3p,1:n_conc)
    ! dep3L_conc(S1p:N1p,S2p:N2p,1:n_conc)
    ! dep3U_conc(S1p:N1p,S2p:N2p,1:n_conc)
    !
    dep1L_conc = 0.
    dep1U_conc = 0.
    dep2L_conc = 0.
    dep2U_conc = 0.
    dep3L_conc = 0.
    dep3U_conc = 0.
  
  
    !--- additional, user-defined initializations ---
    !
    ! for statistics:
#ifdef ALLOC
    ALLOCATE(energy_stokes(1:n_conc))
    ALLOCATE(mass_old     (1:n_conc))
    ALLOCATE(Epot_diff    (1:n_conc))
    ALLOCATE(diss_conc_old(1:n_conc))
#endif
  
    energy_stokes = 0.
    mass_old      = 0.
    Epot_diff     = 0.
    diss_conc_old = 0.
  
  
END SUBROUTINE initial_conditions_conc
  
  
  
  
  
  
  
  
  
  
  
SUBROUTINE initial_conditions_part
    ! (basic subroutine)
  
    USE mod_dims
    USE mod_vars
    USE mod_particles
    USE usr_vars
    USE usr_func
  
    IMPLICIT NONE
  
    !--- under construction ---
  
  
    !--- initial conditions for particle deposit ---
    ! note: - cf. sketch in file "usr_geometry.f90"
    !       - I/O is implemented currently only for "dep1L_part"
    !
    !    orientation of the boundary normal
    !    |lower/upper boundary
    !    ||         grid points on the boundary
    !    ||         |       |
    ! dep1L_part(S2p:N2p,S3p:N3p)
    ! dep1U_part(S2p:N2p,S3p:N3p)
    ! dep2L_part(S1p:N1p,S3p:N3p)
    ! dep2U_part(S1p:N1p,S3p:N3p)
    ! dep3L_part(S1p:N1p,S2p:N2p)
    ! dep3U_part(S1p:N1p,S2p:N2p)
    !
    dep1L_part = 0.
    dep1U_part = 0.
    dep2L_part = 0.
    dep2U_part = 0.
    dep3L_part = 0.
    dep3U_part = 0.
  
  
    !--- kinetic energy of the particle outflow ---
    !
    Ekin_part  = 0.
  
  
    CALL init_particles(n_part_tot,particles)
  
  
END SUBROUTINE initial_conditions_part
  
