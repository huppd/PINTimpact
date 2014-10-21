!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module cmod_solvers
  

    use iso_c_binding
    use mpi
  
    use mod_dims
    use mod_vars
    use mod_exchange
    use mod_diff
    use mod_laplace
    use mod_helmholtz
    use mod_inout
    use mod_solvers
  
    private
  

    !    public outer_iteration, explicit, twostep
    !#ifdef NONBOUSSINESQ
    !    public non_Boussinesq ! TEST!!!
    !#endif
    !    public force_massflow, apply_nullspace, get_nullspace, solve_nullspace
    public csolve_Helmholtz
    !    public solve_conc, solve_conc_explicit
    public cmultigridV, cmultigridF !, restrict, interpolate
!    public restrict_Helmholtz, interpolate_Helmholtz
!    public BiCGstab, Richardson
!    public get_norms, product_scalar, multadd1, multadd2
!    public status_iteration
!    public handle_corner_rhs, handle_corner_rhs_conc
!
!    public BiCGstab2, get_norms2, product_scalar2, apply_nullspace2 ! TEST!!!
!    public relax_restrict, interpolate_relax, plain_restrict, interpolate_mg, relax_bottom
!
!    public get_norms_vel
!
contains
  
    !> \brief solves Helmholtz probelm with IMPACT bicg+multigrid
    !!
    !! \param m direction
    !! \param epsU tolerance
    !! \param n_it_max number of iterations
    !! \param init_yes init starting vector with zero
    !! \param bb right hand side
    !! \param phi solution
    !! \param quiet_yes1
    !! \param quiet_yes2
    subroutine csolve_Helmholtz( &
        m,                      &
        mulL,                   &
        epsU,                   &
        n_it_max,               &
        bb,                     &
        phi,                    &
        cquiet_yes1,             &
        cquiet_yes2 ) bind ( c, name='OP_SolveHelmholtz' )
  
        implicit none
  
        integer(c_int), intent(in   ) ::  m
  
        real(c_double), intent(in   ) ::  mulL

        real(c_double), intent(in   ) ::  epsU
        integer(c_int), intent(in   ) ::  n_it_max
  
        real(c_double), intent(in   ) ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical(c_bool), intent(in  ) ::  cquiet_yes1
        logical(c_bool), intent(in  ) ::  cquiet_yes2
  
        logical ::  quiet_yes1
        logical ::  quiet_yes2

        multL = mulL

        quiet_yes1 = cquiet_yes1
        quiet_yes2 = cquiet_yes2
  
        if (m == 1) then
            direction = 1
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(*,'(a)') 'velocity component 1:'
            call BiCGstab(epsU,n_it_max,.false.,S11B,S21B,S31B,N11B,N21B,N31B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
        end if
  
        if (m == 2) then
            direction = 2
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(*,'(a)') 'velocity component 2:'
            call BiCGstab(epsU,n_it_max,.false.,S12B,S22B,S32B,N12B,N22B,N32B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
        end if
  
        if (m == 3) then
            direction = 3
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(*,'(a)') 'velocity component 3:'
            call BiCGstab(epsU,n_it_max,.false.,S13B,S23B,S33B,N13B,N23B,N33B,bb,phi,1,quiet_yes1,quiet_yes2,precond_Helmh_vel)
        end if
  
  
    end subroutine csolve_Helmholtz
  
  
  
  
  
    !> \brief heavily overloaded fortran interface for GeometricMultiGrid v-solver
    !!
    !! \param[in] problem_type 1: Helmholtz, 3: conc Helmholtz, 2: divGrad, 4: divGrad_trans, 5: divgrad?
    subroutine cmultigridV( &
        m,                  &
        mulL,               &
        cinit_yes,          &
        gstart,             &
        bb,                 &
        phi,                &
        problem_type )  bind ( c, name='OP_MGV' )
  
        implicit none
  
        integer(c_int), intent(in   ) ::  m
        real(c_double), intent(in   ) :: mulL

        logical(c_bool), intent(in  ) ::  cinit_yes

        integer(c_int), intent(in   ) ::  gstart
        integer(c_int), intent(in   ) ::  problem_type
  
        real(c_double), intent(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
        real(c_double), intent(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
        logical ::  init_yes
        integer                ::  g
  
  
         !----------------------------------------------------------------------------------------------------------!
         ! Anmerkungen: - übergebene Felder sind generell zu gross (b1L:(N1+b1U) vs. 0:(N1+1)), andererseits stört  !
         !                der Umstand i.A. auch nicht!                                                              !
         !              - vec1C bzw. vecC können eingespart werden, wenn Produkt und Restriktion zusammengelegt     !
         !                werden. Darauf wird hier der übersicht halber verzichtet, zumal vecC nicht auf dem        !
         !                feinsten Gitter existiert und somit der Speicher nicht überstrapaziert wird. Diese Felder !
         !                werden zudem an die Interpolationsroutinen als Arbeitsfelder übergeben, wären aber auch   !
         !                dort mit entsprechendem Programmieraufwand eliminierbar.                                  !
         !              - INTENT(inout) für bb ist wegen restrict_Helmholtz notwendig.                              !
         !              - Es wird vorausgesetzt, dass das feinstes Gitterlevel keine Nullraum-Korrektur benötigt    !
         !                (sollte normalerweise auch erfüllt sein).                                                 !
         !----------------------------------------------------------------------------------------------------------!

        direction = m

        multL = mulL

        if(cinit_yes) then
            init_yes = .true.
        else
            init_yes = .false.
        end if

  
        !===========================================================================================================
        do g = gstart, n_grids-1
     
            if (participate_yes(g)) then
                if (g == 1 )    call relax_restrict(init_yes,.false.             ,g,psi_rel1 ,bb    ,phi   ,vec1C ,vec2A ,problem_type)
                if (g == gstart) then
                    if (g == 2 ) call relax_restrict(init_yes,.false.             ,g,psi_rel2 ,bb    ,phi   ,vec2C ,vec3A ,problem_type)
                    if (g == 3 ) call relax_restrict(init_yes,.false.             ,g,psi_rel3 ,bb    ,phi   ,vec3C ,vec4A ,problem_type)
                    if (g == 4 ) call relax_restrict(init_yes,.false.             ,g,psi_rel4 ,bb    ,phi   ,vec4C ,vec5A ,problem_type)
                    if (g == 5 ) call relax_restrict(init_yes,.false.             ,g,psi_rel5 ,bb    ,phi   ,vec5C ,vec6A ,problem_type)
                    if (g == 6 ) call relax_restrict(init_yes,.false.             ,g,psi_rel6 ,bb    ,phi   ,vec6C ,vec7A ,problem_type)
                    if (g == 7 ) call relax_restrict(init_yes,.false.             ,g,psi_rel7 ,bb    ,phi   ,vec7C ,vec8A ,problem_type)
                    if (g == 8 ) call relax_restrict(init_yes,.false.             ,g,psi_rel8 ,bb    ,phi   ,vec8C ,vec9A ,problem_type)
                    if (g == 9 ) call relax_restrict(init_yes,.false.             ,g,psi_rel9 ,bb    ,phi   ,vec9C ,vec10A,problem_type)
                    if (g == 10) call relax_restrict(init_yes,.false.             ,g,psi_rel10,bb    ,phi   ,vec10C,vec11A,problem_type)
                    if (g == 11) call relax_restrict(init_yes,.false.             ,g,psi_rel11,bb    ,phi   ,vec11C,vec12A,problem_type)
                    if (g == 12) call relax_restrict(init_yes,.false.             ,g,psi_rel12,bb    ,phi   ,vec12C,vec13A,problem_type)
                    if (g == 13) call relax_restrict(init_yes,.false.             ,g,psi_rel13,bb    ,phi   ,vec13C,vec14A,problem_type)
                    if (g == 14) call relax_restrict(init_yes,.false.             ,g,psi_rel14,bb    ,phi   ,vec14C,vec15A,problem_type)
                else
                    if (g == 2 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel2 ,vec2A ,vec2B ,vec2C ,vec3A ,problem_type)
                    if (g == 3 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel3 ,vec3A ,vec3B ,vec3C ,vec4A ,problem_type)
                    if (g == 4 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel4 ,vec4A ,vec4B ,vec4C ,vec5A ,problem_type)
                    if (g == 5 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel5 ,vec5A ,vec5B ,vec5C ,vec6A ,problem_type)
                    if (g == 6 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel6 ,vec6A ,vec6B ,vec6C ,vec7A ,problem_type)
                    if (g == 7 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel7 ,vec7A ,vec7B ,vec7C ,vec8A ,problem_type)
                    if (g == 8 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel8 ,vec8A ,vec8B ,vec8C ,vec9A ,problem_type)
                    if (g == 9 ) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel9 ,vec9A ,vec9B ,vec9C ,vec10A,problem_type)
                    if (g == 10) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel10,vec10A,vec10B,vec10C,vec11A,problem_type)
                    if (g == 11) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel11,vec11A,vec11B,vec11C,vec12A,problem_type)
                    if (g == 12) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel12,vec12A,vec12B,vec12C,vec13A,problem_type)
                    if (g == 13) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel13,vec13A,vec13B,vec13C,vec14A,problem_type)
                    if (g == 14) call relax_restrict(.true.  ,nullspace_coarse_yes,g,psi_rel14,vec14A,vec14B,vec14C,vec15A,problem_type)
                end if
            end if
     
        end do
        !===========================================================================================================
  
        !--- Grob-Gitter Lösung -------------------------------------
        if (participate_yes(n_grids)) then
            if (gstart == n_grids) then
                call relax_bottom(init_yes,.false.             ,n_grids,psi_rel1 ,bb    ,phi   ,problem_type) ! Achtung: psi_rel1 ist i.A. zu gross! (nullspace_coarse == .FALSE. <-- unschön!)
            else
                if (n_grids == 2 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel2 ,vec2A ,vec2B ,problem_type)
                if (n_grids == 3 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel3 ,vec3A ,vec3B ,problem_type)
                if (n_grids == 4 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel4 ,vec4A ,vec4B ,problem_type)
                if (n_grids == 5 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel5 ,vec5A ,vec5B ,problem_type)
                if (n_grids == 6 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel6 ,vec6A ,vec6B ,problem_type)
                if (n_grids == 7 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel7 ,vec7A ,vec7B ,problem_type)
                if (n_grids == 8 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel8 ,vec8A ,vec8B ,problem_type)
                if (n_grids == 9 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel9 ,vec9A ,vec9B ,problem_type)
                if (n_grids == 10) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel10,vec10A,vec10B,problem_type)
                if (n_grids == 11) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel11,vec11A,vec11B,problem_type)
                if (n_grids == 12) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel12,vec12A,vec12B,problem_type)
                if (n_grids == 13) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel13,vec13A,vec13B,problem_type)
                if (n_grids == 14) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel14,vec14A,vec14B,problem_type)
                if (n_grids == 15) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel15,vec15A,vec15B,problem_type)
            end if
        end if
  
        !===========================================================================================================
        do g = n_grids-1, gstart, -1
     
            if (participate_yes(g)) then
                if (g == gstart) then
                    if (g == 14) call interpolate_relax(g,bb    ,phi   ,vec14C,vec15B,problem_type)
                    if (g == 13) call interpolate_relax(g,bb    ,phi   ,vec13C,vec14B,problem_type)
                    if (g == 12) call interpolate_relax(g,bb    ,phi   ,vec12C,vec13B,problem_type)
                    if (g == 11) call interpolate_relax(g,bb    ,phi   ,vec11C,vec12B,problem_type)
                    if (g == 10) call interpolate_relax(g,bb    ,phi   ,vec10C,vec11B,problem_type)
                    if (g == 9 ) call interpolate_relax(g,bb    ,phi   ,vec9C ,vec10B,problem_type)
                    if (g == 8 ) call interpolate_relax(g,bb    ,phi   ,vec8C ,vec9B ,problem_type)
                    if (g == 7 ) call interpolate_relax(g,bb    ,phi   ,vec7C ,vec8B ,problem_type)
                    if (g == 6 ) call interpolate_relax(g,bb    ,phi   ,vec6C ,vec7B ,problem_type)
                    if (g == 5 ) call interpolate_relax(g,bb    ,phi   ,vec5C ,vec6B ,problem_type)
                    if (g == 4 ) call interpolate_relax(g,bb    ,phi   ,vec4C ,vec5B ,problem_type)
                    if (g == 3 ) call interpolate_relax(g,bb    ,phi   ,vec3C ,vec4B ,problem_type)
                    if (g == 2 ) call interpolate_relax(g,bb    ,phi   ,vec2C ,vec3B ,problem_type)
                else
                    if (g == 14) call interpolate_relax(g,vec14A,vec14B,vec14C,vec15B,problem_type)
                    if (g == 13) call interpolate_relax(g,vec13A,vec13B,vec13C,vec14B,problem_type)
                    if (g == 12) call interpolate_relax(g,vec12A,vec12B,vec12C,vec13B,problem_type)
                    if (g == 11) call interpolate_relax(g,vec11A,vec11B,vec11C,vec12B,problem_type)
                    if (g == 10) call interpolate_relax(g,vec10A,vec10B,vec10C,vec11B,problem_type)
                    if (g == 9 ) call interpolate_relax(g,vec9A ,vec9B ,vec9C ,vec10B,problem_type)
                    if (g == 8 ) call interpolate_relax(g,vec8A ,vec8B ,vec8C ,vec9B ,problem_type)
                    if (g == 7 ) call interpolate_relax(g,vec7A ,vec7B ,vec7C ,vec8B ,problem_type)
                    if (g == 6 ) call interpolate_relax(g,vec6A ,vec6B ,vec6C ,vec7B ,problem_type)
                    if (g == 5 ) call interpolate_relax(g,vec5A ,vec5B ,vec5C ,vec6B ,problem_type)
                    if (g == 4 ) call interpolate_relax(g,vec4A ,vec4B ,vec4C ,vec5B ,problem_type)
                    if (g == 3 ) call interpolate_relax(g,vec3A ,vec3B ,vec3C ,vec4B ,problem_type)
                    if (g == 2 ) call interpolate_relax(g,vec2A ,vec2B ,vec2C ,vec3B ,problem_type)
                end if
                if (g == 1 )    call interpolate_relax(g,bb    ,phi   ,vec1C ,vec2B ,problem_type)
            end if
     
        end do
    !===========================================================================================================
  
  
    end subroutine cmultigridV
  
  
  
  
  
  
  
  
  
  
  
    subroutine cmultigridF(init_yes,gstart,bb,phi,problem_type)
  
        implicit none
  
        logical, intent(in   ) ::  init_yes
        integer, intent(in   ) ::  gstart
        integer, intent(in   ) ::  problem_type
  
        real   , intent(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
        real   , intent(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
        integer                ::  g
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Es wird vorausgesetzt, dass das feinstes Gitterlevel keine Nullraum-Korrektur benötigt    !
        !                (sollte normalerweise auch erfüllt sein)                                                  !
        !              - Es wird immer initialisiert, mit Ausnahme von gstart == n_grids.                          !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        do g = gstart, n_grids-1
     
            if (participate_yes(g)) then
                if (g == 1 )    call plain_restrict(.false.             ,g,psi_rel1 ,bb    ,vec1C ,vec2A ,problem_type)
                if (g == gstart) then
                    if (g == 2 ) call plain_restrict(.false.             ,g,psi_rel2 ,bb    ,vec2C ,vec3A ,problem_type)
                    if (g == 3 ) call plain_restrict(.false.             ,g,psi_rel3 ,bb    ,vec3C ,vec4A ,problem_type)
                    if (g == 4 ) call plain_restrict(.false.             ,g,psi_rel4 ,bb    ,vec4C ,vec5A ,problem_type)
                    if (g == 5 ) call plain_restrict(.false.             ,g,psi_rel5 ,bb    ,vec5C ,vec6A ,problem_type)
                    if (g == 6 ) call plain_restrict(.false.             ,g,psi_rel6 ,bb    ,vec6C ,vec7A ,problem_type)
                    if (g == 7 ) call plain_restrict(.false.             ,g,psi_rel7 ,bb    ,vec7C ,vec8A ,problem_type)
                    if (g == 8 ) call plain_restrict(.false.             ,g,psi_rel8 ,bb    ,vec8C ,vec9A ,problem_type)
                    if (g == 9 ) call plain_restrict(.false.             ,g,psi_rel9 ,bb    ,vec9C ,vec10A,problem_type)
                    if (g == 10) call plain_restrict(.false.             ,g,psi_rel10,bb    ,vec10C,vec11A,problem_type)
                    if (g == 11) call plain_restrict(.false.             ,g,psi_rel11,bb    ,vec11C,vec12A,problem_type)
                    if (g == 12) call plain_restrict(.false.             ,g,psi_rel12,bb    ,vec12C,vec13A,problem_type)
                    if (g == 13) call plain_restrict(.false.             ,g,psi_rel13,bb    ,vec13C,vec14A,problem_type)
                    if (g == 14) call plain_restrict(.false.             ,g,psi_rel14,bb    ,vec14C,vec15A,problem_type)
                else
                    if (g == 2 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel2 ,vec2A ,vec2C ,vec3A ,problem_type)
                    if (g == 3 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel3 ,vec3A ,vec3C ,vec4A ,problem_type)
                    if (g == 4 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel4 ,vec4A ,vec4C ,vec5A ,problem_type)
                    if (g == 5 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel5 ,vec5A ,vec5C ,vec6A ,problem_type)
                    if (g == 6 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel6 ,vec6A ,vec6C ,vec7A ,problem_type)
                    if (g == 7 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel7 ,vec7A ,vec7C ,vec8A ,problem_type)
                    if (g == 8 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel8 ,vec8A ,vec8C ,vec9A ,problem_type)
                    if (g == 9 ) call plain_restrict(nullspace_coarse_yes,g,psi_rel9 ,vec9A ,vec9C ,vec10A,problem_type)
                    if (g == 10) call plain_restrict(nullspace_coarse_yes,g,psi_rel10,vec10A,vec10C,vec11A,problem_type)
                    if (g == 11) call plain_restrict(nullspace_coarse_yes,g,psi_rel11,vec11A,vec11C,vec12A,problem_type)
                    if (g == 12) call plain_restrict(nullspace_coarse_yes,g,psi_rel12,vec12A,vec12C,vec13A,problem_type)
                    if (g == 13) call plain_restrict(nullspace_coarse_yes,g,psi_rel13,vec13A,vec13C,vec14A,problem_type)
                    if (g == 14) call plain_restrict(nullspace_coarse_yes,g,psi_rel14,vec14A,vec14C,vec15A,problem_type)
                end if
            end if
     
        end do
        !===========================================================================================================
  
        !--- Grob-Gitter Lösung -------------------------------------
        if (participate_yes(n_grids)) then
            if (gstart == n_grids) then
                call relax_bottom(init_yes,.false.             ,n_grids,psi_rel1 ,bb    ,phi   ,problem_type) ! Achtung: psi_rel1 ist i.A. zu gross! (nullspace_coarse == .FALSE. <-- unsch�n!)
            else
                if (n_grids == 2 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel2 ,vec2A ,vec2B ,problem_type)
                if (n_grids == 3 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel3 ,vec3A ,vec3B ,problem_type)
                if (n_grids == 4 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel4 ,vec4A ,vec4B ,problem_type)
                if (n_grids == 5 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel5 ,vec5A ,vec5B ,problem_type)
                if (n_grids == 6 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel6 ,vec6A ,vec6B ,problem_type)
                if (n_grids == 7 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel7 ,vec7A ,vec7B ,problem_type)
                if (n_grids == 8 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel8 ,vec8A ,vec8B ,problem_type)
                if (n_grids == 9 ) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel9 ,vec9A ,vec9B ,problem_type)
                if (n_grids == 10) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel10,vec10A,vec10B,problem_type)
                if (n_grids == 11) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel11,vec11A,vec11B,problem_type)
                if (n_grids == 12) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel12,vec12A,vec12B,problem_type)
                if (n_grids == 13) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel13,vec13A,vec13B,problem_type)
                if (n_grids == 14) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel14,vec14A,vec14B,problem_type)
                if (n_grids == 15) call relax_bottom(.true.  ,nullspace_coarse_yes,n_grids,psi_rel15,vec15A,vec15B,problem_type)
            end if
        end if
  
        !===========================================================================================================
        do g = n_grids-1, gstart, -1
     
            if (participate_yes(g)) then
                if (g == gstart) then
                    if (g == 14) call interpolate_mg(g,bb    ,phi   ,vec14C,vec15B,problem_type)
                    if (g == 13) call interpolate_mg(g,bb    ,phi   ,vec13C,vec14B,problem_type)
                    if (g == 12) call interpolate_mg(g,bb    ,phi   ,vec12C,vec13B,problem_type)
                    if (g == 11) call interpolate_mg(g,bb    ,phi   ,vec11C,vec12B,problem_type)
                    if (g == 10) call interpolate_mg(g,bb    ,phi   ,vec10C,vec11B,problem_type)
                    if (g == 9 ) call interpolate_mg(g,bb    ,phi   ,vec9C ,vec10B,problem_type)
                    if (g == 8 ) call interpolate_mg(g,bb    ,phi   ,vec8C ,vec9B ,problem_type)
                    if (g == 7 ) call interpolate_mg(g,bb    ,phi   ,vec7C ,vec8B ,problem_type)
                    if (g == 6 ) call interpolate_mg(g,bb    ,phi   ,vec6C ,vec7B ,problem_type)
                    if (g == 5 ) call interpolate_mg(g,bb    ,phi   ,vec5C ,vec6B ,problem_type)
                    if (g == 4 ) call interpolate_mg(g,bb    ,phi   ,vec4C ,vec5B ,problem_type)
                    if (g == 3 ) call interpolate_mg(g,bb    ,phi   ,vec3C ,vec4B ,problem_type)
                    if (g == 2 ) call interpolate_mg(g,bb    ,phi   ,vec2C ,vec3B ,problem_type)
                else
                    if (g == 14) call interpolate_mg(g,vec14A,vec14B,vec14C,vec15B,problem_type)
                    if (g == 13) call interpolate_mg(g,vec13A,vec13B,vec13C,vec14B,problem_type)
                    if (g == 12) call interpolate_mg(g,vec12A,vec12B,vec12C,vec13B,problem_type)
                    if (g == 11) call interpolate_mg(g,vec11A,vec11B,vec11C,vec12B,problem_type)
                    if (g == 10) call interpolate_mg(g,vec10A,vec10B,vec10C,vec11B,problem_type)
                    if (g == 9 ) call interpolate_mg(g,vec9A ,vec9B ,vec9C ,vec10B,problem_type)
                    if (g == 8 ) call interpolate_mg(g,vec8A ,vec8B ,vec8C ,vec9B ,problem_type)
                    if (g == 7 ) call interpolate_mg(g,vec7A ,vec7B ,vec7C ,vec8B ,problem_type)
                    if (g == 6 ) call interpolate_mg(g,vec6A ,vec6B ,vec6C ,vec7B ,problem_type)
                    if (g == 5 ) call interpolate_mg(g,vec5A ,vec5B ,vec5C ,vec6B ,problem_type)
                    if (g == 4 ) call interpolate_mg(g,vec4A ,vec4B ,vec4C ,vec5B ,problem_type)
                    if (g == 3 ) call interpolate_mg(g,vec3A ,vec3B ,vec3C ,vec4B ,problem_type)
                    if (g == 2 ) call interpolate_mg(g,vec2A ,vec2B ,vec2C ,vec3B ,problem_type)
                end if
                if (g == 1 )    call interpolate_mg(g,bb    ,phi   ,vec1C ,vec2B ,problem_type)
            end if
     
        end do
    !===========================================================================================================
  
  
    end subroutine cmultigridF
  
  


    subroutine clevel_pressure( &
        phi ) bind ( c, name='SF_level' )

        implicit none

        real(c_double),  intent(inout) :: phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        integer                ::  i, j, k
        real                   ::  pre0, pre0_global



        pre0 = 0.

        do k = S3p, N3p
            do j = S2p, N2p
                !pgi$ unroll = n:8
                do i = S1p, N1p
                    pre0 = pre0 + phi(i,j,k)
                end do
            end do
        end do

        call MPI_ALLREDUCE(pre0,pre0_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)

        pre0 = pre0_global/REAL(dim1)/REAL(dim2)/REAL(dim3) ! TEST!!! wegen i4 gefaehrlich!

        phi(S1p:N1p,S2p:N2p,S3p:N3p) = phi(S1p:N1p,S2p:N2p,S3p:N3p) - pre0


    end subroutine clevel_pressure


end module cmod_solvers
