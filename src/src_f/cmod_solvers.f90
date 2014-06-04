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
  
    private
  

    public outer_iteration, explicit, twostep
#ifdef NONBOUSSINESQ
    public non_Boussinesq ! TEST!!!
#endif
    public force_massflow, apply_nullspace, get_nullspace, solve_nullspace
    public solve_Helmholtz
    public solve_conc, solve_conc_explicit
    public multigridV, multigridF, restrict, interpolate
    public restrict_Helmholtz, interpolate_Helmholtz
    public BiCGstab, Richardson
    public get_norms, product_scalar, multadd1, multadd2
    public status_iteration
    public handle_corner_rhs, handle_corner_rhs_conc
  
    public BiCGstab2, get_norms2, product_scalar2, apply_nullspace2 ! TEST!!!
    public relax_restrict, interpolate_relax, plain_restrict, interpolate_mg, relax_bottom
  
    public get_norms_vel
  
contains
  
    !pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
  
  
!    !> brief solver for timeint_mod=0 (CN-RK3).
!    !!
!    !! RHS for Druck-Gleichungssystem:
!    !!   - Helmholtz-Gleichung lösen( solve_Helmholtz )
!    !!   - divergenz bilde( divergence2 )
!    !!   - norm des residuum bilden( get_norms )
!    !! Druckiteration:
!    !!   -
!    subroutine outer_iteration
!
!        implicit none
!
!        real                   ::  norm, norm_prev
!        real                   ::  epsP
!        real                   ::  schur
!
!        integer                ::  i, j, k
!        integer                ::  m, counter
!        logical                ::  exit_yes
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
!        !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
!        !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        if (rank == 0 .and. write_stout_yes .and. concentration_yes .and. timeint_mode == 0) write(*,'(a)') 'flow field:'
!
!        !===========================================================================================================
!        !=== RHS fuer Druck-Gleichungssystem =======================================================================
!        !===========================================================================================================
!        if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Helmholtz equations for RHS ...'
!
!        do m = 1, dimens
!
!            !--- Helmholtz-Gleichung lösen --------------------------------------------------------------------------
!            if (init_pre(substep)) then
!                call solve_Helmholtz(m,epsU,n_it_Helmh_vel,init_vel(substep),rhs(b1L,b2L,b3L,m),vel(b1L,b2L,b3L,m),.true.,.true.)
!            else
!                call gradient(m,pre,gpre)
!
!                if (m == 1) gpre(S11B:N11B,S21B:N21B,S31B:N31B) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
!                if (m == 2) gpre(S12B:N12B,S22B:N22B,S32B:N32B) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
!                if (m == 3) gpre(S13B:N13B,S23B:N23B,S33B:N33B) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
!
!                call solve_Helmholtz(m,epsU,n_it_Helmh_vel,init_vel(substep),gpre,vel(b1L,b2L,b3L,m),.true.,.true.)
!            end if
!
!        end do
!
!        !--- Divergenz bilden --------------------------------------------------------------------------------------
!        call divergence2(vel,res)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Norm des Residuums ====================================================================================
!        !===========================================================================================================
!        if (corner_yes) call handle_corner_rhs(1,res)
!        call get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.true.,.false.,normInf=norm)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Druckiteration ========================================================================================
!        !===========================================================================================================
!        counter = 0
!
!        LOOP: do
!
!            !========================================================================================================
!            !=== Norm an Abbruchkriterium testen ====================================================================
!            !========================================================================================================
!            if (rank == 0 .and. log_iteration_yes) then
!                write(10,'(a)')
!                write(10,'(a)') '-------------------------------------------------------------------------------------'
!                write(10,'(a)') 'current residual:'
!            end if
!
!            call status_iteration(epsU,norm,counter,n_it_outer,exit_yes,.false.,.false.)
!
!            ! ACHTUNG!!! Das findet die xt3 nicht so toll:
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)') '-------------------------------------------------------------------------------------'
!            if (rank == 0 .and. log_iteration_yes) call flush(10)
!
!            if (exit_yes .or. norm == 0.) exit LOOP
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== nächste Iteration ==================================================================================
!            !========================================================================================================
!            counter = counter + 1
!
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)')
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Poisson equations ...'
!
!            epsP = precRatio(substep,counter)*norm
!
!            if (epsP <= epsU) epsP = epsU
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Vorkonditionierer ==================================================================================
!            !========================================================================================================
!            if (precond_outer == 1 .or. precond_outer == 2) then
!                ! RHS korrigieren (sollte eigentlich nicht mehr notwendig sein): ! TEST!!!
!                if (nullspace_yes) call apply_nullspace(res)
!
!                number_poisson = 1
!                call BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.true.,S1p,S2p,S3p,N1p,N2p,N3p,                &
!                    &                       res,work3,2,.true.,.true.,precond_Poisson)
!               !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.TRUE.,S1p,S2p,S3p,N1p,N2p,N3p,               &
!               !           &                       res,work3,2,.TRUE.,.TRUE.,precond_Poisson)
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            if (precond_outer == 2) then
!
!                do m = 1, dimens
!                    direction = m
!                    call gradient(m,work3,gpre)
!                    call bc_extrapolation(m,gpre)
!                    call product_Helmholtz(gpre,work1)
!                    call bc_extrapolation(m,work1)
!                    call divergence(m,work1,work2)
!                end do
!
!                ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
!                if (corner_yes   ) call handle_corner_rhs(1,work2)
!
!                ! RHS korrigieren:
!                if (nullspace_yes) call apply_nullspace(work2)
!
!                number_poisson = 2
!                call BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.false.,S1p,S2p,S3p,N1p,N2p,N3p,               &
!                    &                       work2,work3,2,.true.,.true.,precond_Poisson)
!               !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.FALSE.,S1p,S2p,S3p,N1p,N2p,N3p,              &
!               !           &                       work2,work3,2,.TRUE.,.TRUE.,precond_Poisson)
!            end if
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== aktuelles Residuum =================================================================================
!            !========================================================================================================
!            if (init_pre(substep) .and. counter == 1) then
!                pre(S1p:N1p,S2p:N2p,S3p:N3p) = work3(S1p:N1p,S2p:N2p,S3p:N3p)
!            else
!                pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p)
!               !pre(S1p:N1p,S2p:N2p,S3p:N3p) = 1.25*pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p) ! Relaxation ...
!            end if
!
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)')
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Helmholtz equations for residual ...'
!
!            do m = 1, dimens
!
!                call gradient(m,work3,gpre)
!                call solve_Helmholtz(m,epsU,n_it_Helmh_vel,.true.,gpre,work1,.true.,.true.)
!
!                if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1) - work1(S11B:N11B,S21B:N21B,S31B:N31B)
!                if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2) - work1(S12B:N12B,S22B:N22B,S32B:N32B)
!                if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3) - work1(S13B:N13B,S23B:N23B,S33B:N33B)
!
!            end do
!
!            !--- Divergenz bilden -----------------------------------------------------------------------------------
!            call divergence2(vel,res)
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Norm des Residuums =================================================================================
!            !========================================================================================================
!            norm_prev = norm
!            if (corner_yes) call handle_corner_rhs(1,res)
!            call get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.true.,.false.,normInf=norm)
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Abbruchkriterium für innere Iteration ==============================================================
!            !========================================================================================================
!            precRatio(substep,counter) = norm/norm_prev
!            if (precRatio(substep,counter) > 1.) precRatio(substep,counter) = 1.
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Iterationsstatistik ================================================================================
!            !========================================================================================================
!            ratioO(substep) = ratioO(substep) + LOG10(norm/norm_prev)
!            countO(substep) = countO(substep) + 1
!           !========================================================================================================
!
!        end do LOOP
!
!
!        !===========================================================================================================
!        !=== Massenfluss korrigieren ===============================================================================
!        !===========================================================================================================
!        if (forcing_mode == 2) call force_massflow(.true.)
!
!
!    end subroutine outer_iteration
!
  
  
  
  
  
  
  
  
  
  
!#ifdef NONBOUSSINESQ
!    subroutine non_Boussinesq
!
!        implicit none
!
!        real                   ::  norm, norm_prev
!        real                   ::  epsP
!        real                   ::  schur
!
!        integer                ::  i, j, k
!        integer                ::  m, counter
!        logical                ::  exit_yes
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
!        !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
!        !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        if (rank == 0 .and. write_stout_yes .and. concentration_yes .and. timeint_mode == 0) write(*,'(a)') 'flow field:'
!
!        !===========================================================================================================
!        !=== RHS fuer Druck-Gleichungssystem =======================================================================
!        !===========================================================================================================
!        do m = 1, dimens
!
!            !--- Druckgradient --------------------------------------------------------------------------------------
!            if (init_pre(substep)) then
!                ! TEST!!! weglassen und stattdessen rhs verwenden? Testroutinen lassen sich dann nicht mehr verwenden ...
!                if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1)
!                if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2)
!                if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3)
!            else
!                call gradient(m,pre,gpre)
!
!                if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
!                if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
!                if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
!            end if
!
!            if (m == 1) then
!                do k = S31, N31
!                    do j = S21, N21
!                        !pgi$ unroll = n:8
!                        do i = S11, N11
!                            vel(i,j,k,m) = vel(i,j,k,m)/dens(i,j,k,m)
!                        end do
!                    end do
!                end do
!            end if
!            if (m == 2) then
!                do k = S32, N32
!                    do j = S22, N22
!                        !pgi$ unroll = n:8
!                        do i = S12, N12
!                            vel(i,j,k,m) = vel(i,j,k,m)/dens(i,j,k,m)
!                        end do
!                    end do
!                end do
!            end if
!            if (m == 3) then
!                do k = S33, N33
!                    do j = S23, N23
!                        !pgi$ unroll = n:8
!                        do i = S13, N13
!                            vel(i,j,k,m) = vel(i,j,k,m)/dens(i,j,k,m)
!                        end do
!                    end do
!                end do
!            end if
!
!            !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
!            call bc_extrapolation(m,vel(b1L,b2L,b3L,m))
!
!        end do
!
!        !--- Divergenz bilden --------------------------------------------------------------------------------------
!        call divergence2(vel,res)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Norm des Residuums ====================================================================================
!        !===========================================================================================================
!        if (corner_yes) call handle_corner_rhs(1,res)
!        call get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.true.,.false.,normInf=norm)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Druckiteration ========================================================================================
!        !===========================================================================================================
!        counter = 0
!
!        LOOP: do
!
!            !========================================================================================================
!            !=== Norm an Abbruchkriterium testen ====================================================================
!            !========================================================================================================
!            if (rank == 0 .and. log_iteration_yes) then
!                write(10,'(a)')
!                write(10,'(a)') '-------------------------------------------------------------------------------------'
!                write(10,'(a)') 'current residual:'
!            end if
!
!            call status_iteration(epsU,norm,counter,n_it_outer,exit_yes,.false.,.false.)
!
!            ! ACHTUNG!!! Das findet die xt3 nicht so toll:
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)') '-------------------------------------------------------------------------------------'
!            if (rank == 0 .and. log_iteration_yes) call flush(10)
!
!            if (exit_yes .or. norm == 0.) exit LOOP
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== nächste Iteration ==================================================================================
!            !========================================================================================================
!            counter = counter + 1
!
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)')
!            if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Poisson equations ...'
!
!            epsP = precRatio(substep,counter)*norm
!
!            if (epsP <= epsU) epsP = epsU
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Vorkonditionierer ==================================================================================
!            !========================================================================================================
!            if (precond_outer == 1 .or. precond_outer == 2) then
!
!                number_poisson = 1
!                call BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.true.,S1p,S2p,S3p,N1p,N2p,N3p,                &
!                    &                       res,work3,2,.true.,.true.,precond_Poisson)
!               !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.TRUE.,S1p,S2p,S3p,N1p,N2p,N3p,               &
!               !           &                       res,work3,2,.TRUE.,.TRUE.,precond_Poisson)
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            if (precond_outer == 2) then
!
!                do m = 1, dimens
!                    direction = m
!                    call gradient(m,work3,gpre)
!                    call bc_extrapolation(m,gpre)
!                    call product_Helmholtz(gpre,work1)
!                    call bc_extrapolation(m,work1)
!                    call divergence(m,work1,work2)
!                end do
!
!                ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
!                if (corner_yes   ) call handle_corner_rhs(1,work2)
!
!                number_poisson = 2
!                call BiCGstab  (precOffset(substep,counter)*epsP,n_it_Poisson,.false.,S1p,S2p,S3p,N1p,N2p,N3p,               &
!                    &                       work2,work3,2,.true.,.true.,precond_Poisson)
!               !CALL Richardson(precOffset(substep,counter)*epsP,n_it_Poisson,.FALSE.,S1p,S2p,S3p,N1p,N2p,N3p,              &
!               !           &                       work2,work3,2,.TRUE.,.TRUE.,precond_Poisson)
!            end if
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== aktuelles Residuum =================================================================================
!            !========================================================================================================
!            if (init_pre(substep) .and. counter == 1) then
!                pre(S1p:N1p,S2p:N2p,S3p:N3p) = work3(S1p:N1p,S2p:N2p,S3p:N3p)
!            else
!                pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p)
!               !pre(S1p:N1p,S2p:N2p,S3p:N3p) = 1.25*pre(S1p:N1p,S2p:N2p,S3p:N3p) + work3(S1p:N1p,S2p:N2p,S3p:N3p) ! Relaxation ...
!            end if
!
!            do m = 1, dimens
!
!                call gradient(m,work3,gpre)
!
!                if (m == 1) then
!                    do k = S31, N31
!                        do j = S21, N21
!                            !pgi$ unroll = n:8
!                            do i = S11, N11
!                                gpre(i,j,k) = gpre(i,j,k)/dens(i,j,k,m)
!                            end do
!                        end do
!                    end do
!                end if
!                if (m == 2) then
!                    do k = S32, N32
!                        do j = S22, N22
!                            !pgi$ unroll = n:8
!                            do i = S12, N12
!                                gpre(i,j,k) = gpre(i,j,k)/dens(i,j,k,m)
!                            end do
!                        end do
!                    end do
!                end if
!                if (m == 3) then
!                    do k = S33, N33
!                        do j = S23, N23
!                            !pgi$ unroll = n:8
!                            do i = S13, N13
!                                gpre(i,j,k) = gpre(i,j,k)/dens(i,j,k,m)
!                            end do
!                        end do
!                    end do
!                end if
!
!                !--- Randbedingungen extrapolieren -------------------------------------------------------------------
!                call bc_extrapolation(m,gpre)
!
!                if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
!                if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
!                if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
!
!            end do
!
!            !--- Divergenz bilden -----------------------------------------------------------------------------------
!            call divergence2(vel,res)
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Norm des Residuums =================================================================================
!            !========================================================================================================
!            norm_prev = norm
!            if (corner_yes) call handle_corner_rhs(1,res)
!            call get_norms(S1p,S2p,S3p,N1p,N2p,N3p,res,2,.true.,.false.,normInf=norm)
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Abbruchkriterium für innere Iteration ==============================================================
!            !========================================================================================================
!            precRatio(substep,counter) = norm/norm_prev
!            if (precRatio(substep,counter) > 1.) precRatio(substep,counter) = 1.
!            !========================================================================================================
!
!
!            !========================================================================================================
!            !=== Iterationsstatistik ================================================================================
!            !========================================================================================================
!            ratioO(substep) = ratioO(substep) + LOG10(norm/norm_prev)
!            countO(substep) = countO(substep) + 1
!           !========================================================================================================
!
!        end do LOOP
!
!
!        !===========================================================================================================
!        !=== Massenfluss korrigieren ===============================================================================
!        !===========================================================================================================
!        if (forcing_mode == 2) call force_massflow(.true.)
!
!
!    end subroutine non_Boussinesq
!#endif
  
  
  
  
  
  
  
  
  
    !> \brief "solves" one timestep explicitly( simply advances)
    !!
    !! - RHS fuer Druck-Gleichungssystem
    !! - Poisson Loesung mit bicgstab
    !! - Korektur der Loesung
    !! - Massenfluss korrigieren
!    subroutine explicit
!
!        implicit none
!
!        integer                ::  i, j, k, m
!        real                   ::  schur
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Beim zweiten Poisson-Problem im Vorkonditionierer hat es sich bislang bewährt, die Lösung !
!        !                aus dem ersten Poisson-Problem als Startfeld zu benutzen.                                 !
!        !              - Richardson-Iteration scheint langsamer als BiCGstab zu sein (~30%).                       !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        if (rank == 0 .and. write_stout_yes .and. concentration_yes .and. timeint_mode == 0) write(*,'(a)') 'flow field:'
!
!        !===========================================================================================================
!        !=== RHS fuer Druck-Gleichungssystem =======================================================================
!        !===========================================================================================================
!        do m = 1, dimens
!
!            ! TEST!!! weglassen und stattdessen rhs verwenden? Testroutinen lassen sich dann nicht mehr verwenden ...
!            if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1)
!            if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2)
!            if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3)
!
!            !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
!            call bc_extrapolation(m,vel(b1L,b2L,b3L,m))
!
!        end do
!
!        !--- Divergenz bilden --------------------------------------------------------------------------------------
!        call divergence2(vel,res)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Poisson-Lösung ========================================================================================
!        !===========================================================================================================
!        ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
!        if (corner_yes) call handle_corner_rhs(1,res)
!
!        if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'Poisson problem:'
!
!        number_poisson = 1
!        call BiCGstab  (epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.true.,.true.,precond_Poisson)
!        !CALL Richardson(epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.TRUE.,.TRUE.,precond_Poisson)
!
!        if (rank == 0 .and. log_iteration_yes) call flush(10)
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Korrektur der Lösung ==================================================================================
!        !===========================================================================================================
!        do m = 1, dimens
!
!            call gradient(m,pre,gpre)
!
!            !--- Randbedingungen extrapolieren ----------------------------------------------------------------------
!            call bc_extrapolation(m,gpre)
!
!
!            if (m == 1) then
!                do k = S31B, N31B
!                    do j = S21B, N21B
!                        !pgi$ unroll = n:8
!                        do i = S11B, N11B
!                            vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
!                        end do
!                    end do
!                end do
!            end if
!            if (m == 2) then
!                do k = S32B, N32B
!                    do j = S22B, N22B
!                        !pgi$ unroll = n:8
!                        do i = S12B, N12B
!                            vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
!                        end do
!                    end do
!                end do
!            end if
!            if (m == 3) then
!                do k = S33B, N33B
!                    do j = S23B, N23B
!                        !pgi$ unroll = n:8
!                        do i = S13B, N13B
!                            vel(i,j,k,m) = vel(i,j,k,m) - gpre(i,j,k)
!                        end do
!                    end do
!                end do
!            end if
!
!        end do
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Massenfluss korrigieren ===============================================================================
!        !===========================================================================================================
!        if (forcing_mode == 2) call force_massflow(.false.)
!    !===========================================================================================================
!
!
!    end subroutine explicit
  
  
  
  
  
  
  
  
  
  
  
!    subroutine twostep
!
!        implicit none
!
!        integer                ::  m
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Diese Zeitintegration hat trotz Runge-Kutta-Schema max. Konvergenzordnung eins, da hier   !
!        !                Gp anstelle von H��Gp gerechnet wird, d.h. es liegt ein Fehler in den Impulsgleichungen   !
!        !                vor, abhängig von Delta t/(Re Delta x^2).                                                 !
!        !----------------------------------------------------------------------------------------------------------!
!
!
!        !=== RHS fuer Druck-Gleichungssystem =======================================================================
!        !--- Helmholtz-Gleichung lösen -----------------------------------------------------------------------------
!        if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Helmholtz equation for RHS ...'
!
!        do m = 1, dimens
!            call solve_Helmholtz(m,epsU,n_it_Helmh_vel,.true.,rhs(b1L,b2L,b3L,m),vel(b1L,b2L,b3L,m),.true.,.true.)
!        end do
!
!        !--- Divergenz bilden --------------------------------------------------------------------------------------
!        call divergence2(vel,res)
!        !===========================================================================================================
!
!
!        !=== H-transformierten Druck bestimmen =====================================================================
!        if (rank == 0 .and. log_iteration_yes) write(10,'(a)')
!        if (rank == 0 .and. log_iteration_yes) write(10,'(a)') 'solving Poisson equation ...'
!
!        ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
!        if (corner_yes   ) call handle_corner_rhs(1,res)
!
!        call BiCGstab(epsU,n_it_Poisson,init_pre(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,pre,2,.true.,.true.,precond_Poisson)
!        !===========================================================================================================
!
!
!        !=== Geschwindigkeiten bestimmen ===========================================================================
!        do m = 1, dimens
!            call gradient(m,pre,gpre)
!
!            if (m == 1) vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
!            if (m == 2) vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
!            if (m == 3) vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
!        end do
!        !===========================================================================================================
!
!
!        !===========================================================================================================
!        !=== Massenfluss korrigieren ===============================================================================
!        !===========================================================================================================
!        if (forcing_mode == 2) call force_massflow(.false.)
!    !===========================================================================================================
!
!
!    end subroutine twostep
  
  
  
  
  
  
  
  
  
  
  
    subroutine force_massflow(impl_yes)
  
        implicit none
  
        logical, intent(in)    ::  impl_yes
  
        integer                ::  i, j, k, m
        real                   ::  flux(1:2), flux_global(1:2), dV
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Massenflusskorrektur kann Divergenzfreiheit nicht stoeren, daher ist nachtraegliche       !
        !                Korrektur moeglich:                                                                       !
        !                                                                                                          !
        ! m  = desired massflow                                                                                    !
        !                                                                                                          !
        ! A  = |H  G|                                                                                              !
        !      |D  0|                                                                                              !
        !                                                                                                          !
        ! A^-1 = |H^-1*(I-G*(D*H^-1*G)^-1*D*H^-1)  ... |                                                           !
        !        |          (D*H^-1*G)^-1*D*H^-1   ... |                                                           !
        ! x  = [u ,p]^T                                                                                            !
        ! b  = [f ,0]^T                                                                                            !
        ! j  = [1 ,0]^T                                                                                            !
        ! jt = [1t,0]^T                                                                                            !
        !                                                                                                          !
        ! |A  j|*|x| = |b|  ==>  |I     A^-1*j|*|x| = |   A^-1*b  |                                                !
        ! |jt 0| |k|   |m|       |0  jt*A^-1*j| |k|   |jt*A^-1*b-m|                                                !
        !                                                                                                          !
        ! s := jt*A^-1                                                                                             !
        ! t :=    A^-1*j                                                                                           !
        !                                                                                                          !
        ! ==> x = A^-1*b - t*(s*b-m)/(s*j)                                                                         !
        !                                                                                                          !
        !                                                                                                          !
        ! s = |1t*H^-1 ... |,  because jt*H^-1*G = 0 on equidistant grid in streamwise direction                   !
        !                                                                                                          !
        ! t = |H^-1*1|,        because  D*H^-1*j = 0                                                               !
        !     |  0   |                                                                                             !
        !                                                                                                          !
        ! ==> u = A^-1*f - H^-1*1*(1t*H^-1*f-m)/(1t*H^-1*1)                                                        !
        !----------------------------------------------------------------------------------------------------------!
  
  
        m = bulkflow_dir
        direction = m
  
        !===========================================================================================================
        if (m == 1) then
     
            if (impl_yes) then
                gpre = 1.
        
                if (BC_2L > 0) gpre(S11B:N11B,1 ,S31B:N31B) = 0.
                if (BC_2U > 0) gpre(S11B:N11B,N2,S31B:N31B) = 0.
                if (BC_3L > 0) gpre(S11B:N11B,S21B:N21B,1 ) = 0.
                if (BC_3U > 0) gpre(S11B:N11B,S21B:N21B,N3) = 0.
        
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S11B,S21B,S31B,N11B,N21B,N31B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
            else
                dig = 1.
        
                if (BC_2L > 0) dig (S11B:N11B,1 ,S31B:N31B) = 0.
                if (BC_2U > 0) dig (S11B:N11B,N2,S31B:N31B) = 0.
                if (BC_3L > 0) dig (S11B:N11B,S21B:N21B,1 ) = 0.
                if (BC_3U > 0) dig (S11B:N11B,S21B:N21B,N3) = 0.
            end if
     
            !--------------------------------------------------------------------------------------------------------
            !--- Massenfluss ----------------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
            ! flux(1) = 1t*H^-1*1
            ! flux(2) =    H^-1*f
     
            flux = 0.
     
            do k = S31B, N31B
                do j = S21B, N21B
                    !pgi$ unroll = n:8
                    do i = S11B, N11B ! Anmerkung: dx1p(i) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
                        dV = dx1p(i)*dx2p(j)*dx3p(k)
                        flux(1) = flux(1) + dig(i,j,k  )*dV
                        flux(2) = flux(2) + vel(i,j,k,m)*dV
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
            !--------------------------------------------------------------------------------------------------------
            !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
            gpre(S11B:N11B,S21B:N21B,S31B:N31B) = flux(1)
     
            if (BC_2L > 0) gpre(S11B:N11B,1 ,S31B:N31B) = 0.
            if (BC_2U > 0) gpre(S11B:N11B,N2,S31B:N31B) = 0.
            if (BC_3L > 0) dig (S11B:N11B,S21B:N21B,1 ) = 0.
            if (BC_3U > 0) dig (S11B:N11B,S21B:N21B,N3) = 0.
     
            if (impl_yes) then
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S11B,S21B,S31B,N11B,N21B,N31B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
                vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - dig (S11B:N11B,S21B:N21B,S31B:N31B)
            else
                vel(S11B:N11B,S21B:N21B,S31B:N31B,m) = vel(S11B:N11B,S21B:N21B,S31B:N31B,m) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
            end if
     
        !===========================================================================================================
        else if (m == 2) then
     
            if (impl_yes) then
                gpre = 1.
        
                if (BC_1L > 0) gpre(1 ,S22B:N22B,S32B:N32B) = 0.
                if (BC_1U > 0) gpre(N1,S22B:N22B,S32B:N32B) = 0.
                if (BC_3L > 0) gpre(S12B:N12B,S22B:N22B,1 ) = 0.
                if (BC_3U > 0) gpre(S12B:N12B,S22B:N22B,N3) = 0.
        
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S12B,S22B,S32B,N12B,N22B,N32B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
            else
                dig = 1.
        
                if (BC_1L > 0) dig (1 ,S22B:N22B,S32B:N32B) = 0.
                if (BC_1U > 0) dig (N1,S22B:N22B,S32B:N32B) = 0.
                if (BC_3L > 0) dig (S12B:N12B,S22B:N22B,1 ) = 0.
                if (BC_3U > 0) dig (S12B:N12B,S22B:N22B,N3) = 0.
            end if
     
            !--------------------------------------------------------------------------------------------------------
            !--- Massenfluss ----------------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
            ! flux(1) = 1t*H^-1*1
            ! flux(2) =    H^-1*f
     
            flux = 0.
     
            do k = S32B, N32B
                do j = S22B, N22B
                    !pgi$ unroll = n:8
                    do i = S12B, N12B ! Anmerkung: dx2p(j) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
                        dV = dx1p(i)*dx2p(j)*dx3p(k)
                        flux(1) = flux(1) + dig(i,j,k  )*dV
                        flux(2) = flux(2) + vel(i,j,k,m)*dV
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
            !--------------------------------------------------------------------------------------------------------
            !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
            gpre(S12B:N12B,S22B:N22B,S32B:N32B) = flux(1)
     
            if (BC_1L > 0) gpre(1 ,S22B:N22B,S32B:N32B) = 0.
            if (BC_1U > 0) gpre(N1,S22B:N22B,S32B:N32B) = 0.
            if (BC_3L > 0) dig (S12B:N12B,S22B:N22B,1 ) = 0.
            if (BC_3U > 0) dig (S12B:N12B,S22B:N22B,N3) = 0.
     
            if (impl_yes) then
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S12B,S22B,S32B,N12B,N22B,N32B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
                vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - dig (S12B:N12B,S22B:N22B,S32B:N32B)
            else
                vel(S12B:N12B,S22B:N22B,S32B:N32B,m) = vel(S12B:N12B,S22B:N22B,S32B:N32B,m) - gpre(S12B:N12B,S22B:N22B,S32B:N32B)
            end if
     
        !===========================================================================================================
        else if (m == 3) then
     
            if (impl_yes) then
                gpre = 1.
        
                if (BC_1L > 0) gpre(1 ,S23B:N23B,S33B:N33B) = 0.
                if (BC_1U > 0) gpre(N1,S23B:N23B,S33B:N33B) = 0.
                if (BC_2L > 0) gpre(S13B:N13B,1 ,S33B:N33B) = 0.
                if (BC_2U > 0) gpre(S13B:N13B,N2,S33B:N33B) = 0.
        
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S13B,S23B,S33B,N13B,N23B,N33B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
            else
                dig = 1.
        
                if (BC_1L > 0) dig (1 ,S23B:N23B,S33B:N33B) = 0.
                if (BC_1U > 0) dig (N1,S23B:N23B,S33B:N33B) = 0.
                if (BC_2L > 0) dig (S13B:N13B,1 ,S33B:N33B) = 0.
                if (BC_2U > 0) dig (S13B:N13B,N2,S33B:N33B) = 0.
            end if
     
            !--------------------------------------------------------------------------------------------------------
            !--- Massenfluss ----------------------------------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            ! Anmerkung: Eigentlich wuerde eine einzige Gitterpunktflaeche reichen, so ist es aber genauer.
            ! flux(1) = 1t*H^-1*1
            ! flux(2) =    H^-1*f
     
            flux = 0.
     
            do k = S33B, N33B
                do j = S23B, N23B
                    !pgi$ unroll = n:8
                    do i = S13B, N13B ! Anmerkung: dx3p(k) darf verwendet werden, weil Gitter ohnehin aequidistant (s. Anmerkungen!)
                        dV = dx1p(i)*dx2p(j)*dx3p(k)
                        flux(1) = flux(1) + dig(i,j,k  )*dV
                        flux(2) = flux(2) + vel(i,j,k,m)*dV
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE(flux(1:2),flux_global(1:2),2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
            !--------------------------------------------------------------------------------------------------------
            !--- Forcing (k=(H^-1*f-m)/(1t*H^-1*1)) -----------------------------------------------------------------
            !--------------------------------------------------------------------------------------------------------
            flux(1) = (flux_global(2) - vel_bulk*L1*L2*L3) / flux_global(1) ! TEST!!! nur auf Kanalstroemung beschraenkt!!
     
            gpre(S13B:N13B,S23B:N23B,S33B:N33B) = flux(1)
     
            if (BC_1L > 0) gpre(1 ,S23B:N23B,S33B:N33B) = 0.
            if (BC_1U > 0) gpre(N1,S23B:N23B,S33B:N33B) = 0.
            if (BC_2L > 0) dig (S13B:N13B,1 ,S33B:N33B) = 0.
            if (BC_2U > 0) dig (S13B:N13B,N2,S33B:N33B) = 0.
     
            if (impl_yes) then
                call BiCGstab(epsU,n_it_Helmh_vel,init_vel(substep),S13B,S23B,S33B,N13B,N23B,N33B,gpre,dig,1,.true.,.true.,precond_Helmh_vel)
                vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - dig (S13B:N13B,S23B:N23B,S33B:N33B)
            else
                vel(S13B:N13B,S23B:N23B,S33B:N33B,m) = vel(S13B:N13B,S23B:N23B,S33B:N33B,m) - gpre(S13B:N13B,S23B:N23B,S33B:N33B)
            end if
     
        end if
    !===========================================================================================================
  
  
    end subroutine force_massflow
  
  
  
  
  
  
  
  
  
  
  
    subroutine apply_nullspace(res)
  
        implicit none
  
        real   , intent(inout) ::  res(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  flux, flux_global
        integer                ::  i, j, k
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! note: - null space correction of pressure RHS and velocity RHS is currently not "consistent" (even the   !
        !         orthogonal projection is NOT equivalent!)                                                        !
        !       - null space correction of velocity RHS should make correction of pressure RHS obsolete            !
        !----------------------------------------------------------------------------------------------------------!
  
  
        flux = 0.
  
        do k = S3p, N3p
            do j = S2p, N2p
                do i = S1p, N1p
                    flux = flux + psi(i,j,k)*res(i,j,k)
                end do
            end do
        end do
  
        call MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
        flux = flux_global
  
        if (rank == 0) write(10,'(a,E25.17)') 'flux =', flux
  
  
        ! Orthogonale Projektion:
        res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - flux*psi(S1p:N1p,S2p:N2p,S3p:N3p)
  
  
    end subroutine apply_nullspace
  
  
  
  
  
  
  
  
  
  
  
    subroutine apply_nullspace2(g,psi,res,problem_type) ! TEST!!! Original ersetzen!
  
        implicit none
  
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(inout) ::  res(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
        real                   ::  flux, flux_global
  
        integer                ::  S1R, N1R, i, dim1
        integer                ::  S2R, N2R, j, dim2
        integer                ::  S3R, N3R, k, dim3
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! note: - null space correction of pressure RHS and velocity RHS is currently not "consistent" (even the   !
        !         orthogonal projection is NOT equivalent!)                                                        !
        !       - also RHS correction of multigrid Poisson problem (finest grid) is normally slightly different    !
        !         from actual Poisson problem                                                                      !
        !----------------------------------------------------------------------------------------------------------!
  
  
        S1R  = SNB(1,1,g)
        S2R  = SNB(1,2,g)
        S3R  = SNB(1,3,g)
  
        N1R  = SNB(2,1,g)
        N2R  = SNB(2,2,g)
        N3R  = SNB(2,3,g)
  
        flux = 0.
  
        if (problem_type == 2) then
            do k = S3R, N3R
                do j = S2R, N2R
                    do i = S1R, N1R
                        flux = flux + psi(i,j,k)*res(i,j,k)
                    end do
                end do
            end do
        end if
  
        if (problem_type == 4 .or. problem_type == 5) then
            do k = S3R, N3R
                do j = S2R, N2R
                    do i = S1R, N1R
                        flux = flux + res(i,j,k) ! The left null space of the transposed matrix is a constant ...
                    end do
                end do
            end do
        end if
  
        call MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
        flux = flux_global
  
        !IF (rank == 0) WRITE(10,'(a,i3,a,E25.17)') 'grid #', g, ', flux =', flux
  
  
        ! Orthogonale Projektion:
        if (problem_type == 2) then
            res(S1R:N1R,S2R:N2R,S3R:N3R) = res(S1R:N1R,S2R:N2R,S3R:N3R) - flux*psi(S1R:N1R,S2R:N2R,S3R:N3R)
        end if
  
        if (problem_type == 4 .or. problem_type == 5) then
            dim1 = (NN(1,g)-1)*NB(1,g)+1
            dim2 = (NN(2,g)-1)*NB(2,g)+1
            dim3 = (NN(3,g)-1)*NB(3,g)+1
            if (BC_1L_global == -1) dim1 = dim1-1
            if (BC_2L_global == -1) dim2 = dim2-1
            if (BC_3L_global == -1) dim3 = dim3-1
     
            flux = flux / REAL(dim1*dim2*dim3) ! hier keine Wurzel, weil oben auch nicht mit einem psi durchmultipliziert wird!
     
            res(S1R:N1R,S2R:N2R,S3R:N3R) = res(S1R:N1R,S2R:N2R,S3R:N3R) - flux
        end if
  
        if (corner_yes) call handle_corner_rhs(g,res) ! TEST!!! verifizieren ...
  
  
    end subroutine apply_nullspace2
  
  
  
  
  
  
  
  
  
  
  
    subroutine get_nullspace
  
        implicit none
  
        integer                ::  g, m
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        real   , allocatable   ::  work1(:,:,:)
        real   , allocatable   ::  work2(:,:,:)
  
        real                   ::  eps, eps_global
        real                   ::  norm(1:3), norm_global(1:3)
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: -                                                                                           !
        !----------------------------------------------------------------------------------------------------------!
  
  
        ! Richtig wäre eine Korrektur direkt von rhs, so dass <PSI,rhs> = 0: PSI = H^-1 D^T psi mit D H^-1 G psi^T = 0
        ! Allerdings: Viel zu aufwändig!
  
  
        !===========================================================================================================
        !=== transponierte Stencils bilden =========================================================================
        !===========================================================================================================
        cDu1T = 0.
        cDv2T = 0.
        cDw3T = 0.
  
        cGp1T = 0.
        cGp2T = 0.
        cGp3T = 0.
        !-----------------------------------------------------------------------------------------------------------
        allocate(work1(b1L:b1U,b1L:(M1+b1U),1:2))
        allocate(work2(b1L:b1U,b1L:(M1+b1U),1:2))
  
        work1 = 0.
        do i = S1p, N1p
            work1(d1L:d1U,i+iShift,1) = cDu1(d1L:d1U,i)
        end do
        do i = S11B, N11B
            work1(g1L:g1U,i+iShift,2) = cGp1(g1L:g1U,i)
        end do
  
        call MPI_ALLREDUCE(work1,work2,2*(b1U-b1L+1)*(M1+b1U-b1L+1),MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
  
        if (BC_1L_global == -1) then
            do i = b1L, -1
                work2(b1L:b1U,2 +ls1+i,1) = work2(b1L:b1U,M1+ls1+1+i,1)
                work2(b1L:b1U,2 +ls1+i,2) = work2(b1L:b1U,M1+ls1+1+i,2)
            end do
            do i = 1, b1U
                work2(b1L:b1U,M1+ls1+i,1) = work2(b1L:b1U,1 +ls1  +i,1)
                work2(b1L:b1U,M1+ls1+i,2) = work2(b1L:b1U,1 +ls1  +i,2)
            end do
        end if
  
        do i = S11B, N11B
            do ii = g1L, g1U
                cDu1T(ii,i) = work2(-ii,i+ii+iShift,1)
            end do
        end do
  
        do i = S1p, N1p
            do ii = d1L, d1U
                cGp1T(ii,i) = work2(-ii,i+ii+iShift,2)
            end do
        end do
  
        deallocate(work1)
        deallocate(work2)
        !-----------------------------------------------------------------------------------------------------------
        allocate(work1(b2L:b2U,b2L:(M2+b2U),1:2))
        allocate(work2(b2L:b2U,b2L:(M2+b2U),1:2))
  
        work1 = 0.
        do j = S2p, N2p
            work1(d2L:d2U,j+jShift,1) = cDv2(d2L:d2U,j)
        end do
        do j = S22B, N22B
            work1(g2L:g2U,j+jShift,2) = cGp2(g2L:g2U,j)
        end do
  
        call MPI_ALLREDUCE(work1,work2,2*(b2U-b2L+1)*(M2+b2U-b2L+1),MPI_REAL8,MPI_SUM,COMM_BAR2,merror)
  
        if (BC_2L_global == -1) then
            do j = b2L, -1
                work2(b2L:b2U,2 +ls2+j,1) = work2(b2L:b2U,M2+ls2+1+j,1)
                work2(b2L:b2U,2 +ls2+j,2) = work2(b2L:b2U,M2+ls2+1+j,2)
            end do
            do j = 1, b2U
                work2(b2L:b2U,M2+ls2+j,1) = work2(b2L:b2U,1 +ls2  +j,1)
                work2(b2L:b2U,M2+ls2+j,2) = work2(b2L:b2U,1 +ls2  +j,2)
            end do
        end if
  
        do j = S22B, N22B
            do jj = g2L, g2U
                cDv2T(jj,j) = work2(-jj,j+jj+jShift,1)
            end do
        end do
  
        do j = S2p, N2p
            do jj = d2L, d2U
                cGp2T(jj,j) = work2(-jj,j+jj+jShift,2)
            end do
        end do
  
        deallocate(work1)
        deallocate(work2)
        !-----------------------------------------------------------------------------------------------------------
        if (dimens == 3) then
  
            allocate(work1(b3L:b3U,b3L:(M3+b3U),1:2))
            allocate(work2(b3L:b3U,b3L:(M3+b3U),1:2))
  
            work1 = 0.
            do k = S3p, N3p
                work1(d3L:d3U,k+kShift,1) = cDw3(d3L:d3U,k)
            end do
            do k = S33B, N33B
                work1(g3L:g3U,k+kShift,2) = cGp3(g3L:g3U,k)
            end do
  
            call MPI_ALLREDUCE(work1,work2,2*(b3U-b3L+1)*(M3+b3U-b3L+1),MPI_REAL8,MPI_SUM,COMM_BAR3,merror)
  
            if (BC_3L_global == -1) then
                do k = b3L, -1
                    work2(b3L:b3U,2 +ls3+k,1) = work2(b3L:b3U,M3+ls3+1+k,1)
                    work2(b3L:b3U,2 +ls3+k,2) = work2(b3L:b3U,M3+ls3+1+k,2)
                end do
                do k = 1, b3U
                    work2(b3L:b3U,M3+ls3+k,1) = work2(b3L:b3U,1 +ls3  +k,1)
                    work2(b3L:b3U,M3+ls3+k,2) = work2(b3L:b3U,1 +ls3  +k,2)
                end do
            end if
  
            do k = S33B, N33B
                do kk = g3L, g3U
                    cDw3T(kk,k) = work2(-kk,k+kk+kShift,1)
                end do
            end do
  
            do k = S3p, N3p
                do kk = d3L, d3U
                    cGp3T(kk,k) = work2(-kk,k+kk+kShift,2)
                end do
            end do
  
            deallocate(work1)
            deallocate(work2)
  
        end if
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== Null-Raum lesen / berechnen ===========================================================================
        !===========================================================================================================
        if (rank == 0 .and. write_stout_yes) then
            if (read_nullspace_yes) then
                write(*,'(a)') 'reading left null-space of preconditioner ...'
            else
                write(*,'(a)') 'determining left null-space of preconditioner ...'
            end if
        end if
  
        ! rein empirisch:
        eps = 0.001*epsU
  
        if (.not. read_nullspace_yes) then
            open(10,file='test_log_iterations_nullspace.txt',status='UNKNOWN')
     
            if (rank == 0) write(10,'(a,E27.15)') ' accuracy for null-space solution: eps =', eps
            if (rank == 0) write(10,*)
        end if
  
        call solve_nullspace(eps,1,psi,4)
  
        if (nullspace_coarse_yes) then
            if (n_grids >= 1 ) call solve_nullspace(eps,1 ,psi_rel1 ,5)
            if (n_grids >= 2 ) call solve_nullspace(eps,2 ,psi_rel2 ,5)
            if (n_grids >= 3 ) call solve_nullspace(eps,3 ,psi_rel3 ,5)
            if (n_grids >= 4 ) call solve_nullspace(eps,4 ,psi_rel4 ,5)
            if (n_grids >= 5 ) call solve_nullspace(eps,5 ,psi_rel5 ,5)
            if (n_grids >= 6 ) call solve_nullspace(eps,6 ,psi_rel6 ,5)
            if (n_grids >= 7 ) call solve_nullspace(eps,7 ,psi_rel7 ,5)
            if (n_grids >= 8 ) call solve_nullspace(eps,8 ,psi_rel8 ,5)
            if (n_grids >= 9 ) call solve_nullspace(eps,9 ,psi_rel9 ,5)
            if (n_grids >= 10) call solve_nullspace(eps,10,psi_rel10,5)
            if (n_grids >= 11) call solve_nullspace(eps,11,psi_rel11,5)
            if (n_grids >= 12) call solve_nullspace(eps,12,psi_rel12,5)
            if (n_grids >= 13) call solve_nullspace(eps,13,psi_rel13,5)
            if (n_grids >= 14) call solve_nullspace(eps,14,psi_rel14,5)
            if (n_grids >= 15) call solve_nullspace(eps,15,psi_rel15,5)
        end if
  
        if (.not. read_nullspace_yes) close(10)
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== Null-Raum fuer vel ====================================================================================
        !===========================================================================================================
        ! Richtig wäre eine Korrektur direkt von rhs, so dass <PSI,rhs> = 0: PSI = H^-1 D^T psi mit D H^-1 G psi^T = 0
        ! Allerdings: Viel zu aufwändig!
        psi_vel = 0. ! ACHTUNG!!! Nur für timeint_mode == 1
        do m = 1, dimens
            call divergence_transp  (m,psi,psi_vel(b1L,b2L,b3L,m))
            call bc_extrapolation_transp(m,psi_vel(b1L,b2L,b3L,m))
        end do
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        !=== Normieren =============================================================================================
        !===========================================================================================================
        rhs  = 0.
        if (nullspace_ortho_yes) then
            rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = psi_vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
            rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = psi_vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
            if (dimens == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = psi_vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        else
            ! Umspeichern auf rhs notwendig wegen überlappung an Ecken/Kanten ...
            !--------------------------------------------------------------------------------------------------------
            if (BC_2L > 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,1)
            if (BC_2U > 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,2)
     
            if (BC_3L > 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = th13(S11B:N11B,S21B:N21B,1)
            if (BC_3U > 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = th13(S11B:N11B,S21B:N21B,2)
     
            if (BC_1L > 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,1)
            if (BC_1U > 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,2)
            !--------------------------------------------------------------------------------------------------------
            if (BC_1L > 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,1)
            if (BC_1U > 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,2)
     
            if (BC_3L > 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = th23(S12B:N12B,S22B:N22B,1)
            if (BC_3U > 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = th23(S12B:N12B,S22B:N22B,2)
     
            if (BC_2L > 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,1)
            if (BC_2U > 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,2)
            !--------------------------------------------------------------------------------------------------------
            if (BC_1L > 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = th31(S23B:N23B,S33B:N33B,1)
            if (BC_1U > 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = th31(S23B:N23B,S33B:N33B,2)
     
            if (BC_2L > 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = th32(S13B:N13B,S33B:N33B,1)
            if (BC_2U > 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = th32(S13B:N13B,S33B:N33B,2)
     
            if (BC_3L > 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = th33(S13B:N13B,S23B:N23B,1)
            if (BC_3U > 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = th33(S13B:N13B,S23B:N23B,2)
           !--------------------------------------------------------------------------------------------------------
        end if
        !===========================================================================================================
        norm = 0.
  
        m = 1
        do k = S31B, N31B
            do j = S21B, N21B
                do i = S11B, N11B
                    norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
                    norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
                    norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
                end do
            end do
        end do
  
        m = 2
        do k = S32B, N32B
            do j = S22B, N22B
                do i = S12B, N12B
                    norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
                    norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
                    norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
                end do
            end do
        end do
  
        if (dimens == 3) then
            m = 3
            do k = S33B, N33B
                do j = S23B, N23B
                    do i = S13B, N13B
                        norm(1) = norm(1) + psi_vel(i,j,k,m)*rhs    (i,j,k,m)
                        norm(2) = norm(2) + psi_vel(i,j,k,m)*psi_vel(i,j,k,m)
                        norm(3) = norm(3) + rhs    (i,j,k,m)*rhs    (i,j,k,m)
                    end do
                end do
            end do
        end if
  
        call MPI_ALLREDUCE(norm,norm_global,3,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  
        norm(1) =      norm_global(1)
        norm(2) = SQRT(norm_global(2))
        norm(3) = SQRT(norm_global(3))
  
        if (rank == 0) write(*,*)
        if (rank == 0) write(*,*) '       |psi_vel| =', norm(2)
        if (rank == 0) write(*,*) '         |theta| =', norm(3)
        if (rank == 0) write(*,*) ' <psi_vel,theta> =', norm(1)
        if (rank == 0) write(*,*) ' (note that |<psi_vel,theta>| should not be much smaller than |psi_vel||theta|.)'
        if (rank == 0) write(*,*)
  
        if (norm(1) == 0.) then
            if (rank == 0) write(*,*) 'ERROR! <psi_vel,theta> == 0.'
            call MPI_FINALIZE(merror)
            stop
        end if
  
        norm(1) = 1./SQRT(norm_global(1))
  
        psi_vel = psi_vel * norm(1)
  
        if (.not. nullspace_ortho_yes) then
            !--------------------------------------------------------------------------------------------------------
            if (BC_2L > 0) th12(S11B:N11B,S31B:N31B,1) = th12(S11B:N11B,S31B:N31B,1) * norm(1)
            if (BC_2U > 0) th12(S11B:N11B,S31B:N31B,2) = th12(S11B:N11B,S31B:N31B,2) * norm(1)
     
            if (BC_3L > 0) th13(S11B:N11B,S21B:N21B,1) = th13(S11B:N11B,S21B:N21B,1) * norm(1)
            if (BC_3U > 0) th13(S11B:N11B,S21B:N21B,2) = th13(S11B:N11B,S21B:N21B,2) * norm(1)
     
            if (BC_1L > 0) th11(S21B:N21B,S31B:N31B,1) = th11(S21B:N21B,S31B:N31B,1) * norm(1)
            if (BC_1U > 0) th11(S21B:N21B,S31B:N31B,2) = th11(S21B:N21B,S31B:N31B,2) * norm(1)
            !--------------------------------------------------------------------------------------------------------
            if (BC_1L > 0) th21(S22B:N22B,S32B:N32B,1) = th21(S22B:N22B,S32B:N32B,1) * norm(1)
            if (BC_1U > 0) th21(S22B:N22B,S32B:N32B,2) = th21(S22B:N22B,S32B:N32B,2) * norm(1)
     
            if (BC_3L > 0) th23(S12B:N12B,S22B:N22B,1) = th23(S12B:N12B,S22B:N22B,1) * norm(1)
            if (BC_3U > 0) th23(S12B:N12B,S22B:N22B,2) = th23(S12B:N12B,S22B:N22B,2) * norm(1)
     
            if (BC_2L > 0) th22(S12B:N12B,S32B:N32B,1) = th22(S12B:N12B,S32B:N32B,1) * norm(1)
            if (BC_2U > 0) th22(S12B:N12B,S32B:N32B,2) = th22(S12B:N12B,S32B:N32B,2) * norm(1)
            !--------------------------------------------------------------------------------------------------------
            if (BC_1L > 0) th31(S23B:N23B,S33B:N33B,1) = th31(S23B:N23B,S33B:N33B,1) * norm(1)
            if (BC_1U > 0) th31(S23B:N23B,S33B:N33B,2) = th31(S23B:N23B,S33B:N33B,2) * norm(1)
     
            if (BC_2L > 0) th32(S13B:N13B,S33B:N33B,1) = th32(S13B:N13B,S33B:N33B,1) * norm(1)
            if (BC_2U > 0) th32(S13B:N13B,S33B:N33B,2) = th32(S13B:N13B,S33B:N33B,2) * norm(1)
     
            if (BC_3L > 0) th33(S13B:N13B,S23B:N23B,1) = th33(S13B:N13B,S23B:N23B,1) * norm(1)
            if (BC_3U > 0) th33(S13B:N13B,S23B:N23B,2) = th33(S13B:N13B,S23B:N23B,2) * norm(1)
           !--------------------------------------------------------------------------------------------------------
        end if
    !===========================================================================================================
  
  
    end subroutine get_nullspace
  
  
  
  
  
  
  
  
  
  
  
    subroutine solve_nullspace(eps,g,psi,problem_type)
  
        implicit none
  
        real   , intent(in   ) ::  eps
        integer, intent(in   ) ::  g
        integer, intent(in   ) ::  problem_type
  
        real   , intent(inout) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
        real                   ::  norm2, norm2_global
        integer                ::  i, j, k
        integer                ::  n_grids_orig
  
        real                   ::  line1(b1L:(N1+b1U),1:2)
        real                   ::  line2(b2L:(N2+b2U),1:2)
        real                   ::  line3(b3L:(N3+b3U),1:2)
  
  
        !--- KEIN Multigrid, nur Glätter auf feinem Gitter ---
        n_grids_orig = n_grids ! TEST!!!
        n_grids      = g
  
  
        if (g == 1 .and. problem_type == 4) then
     
            if (1 == 2) then
        
                !=====================================================================================================
                if (read_nullspace_yes) then
                    !--- Lesen ---
                    call read2_hdf('nullspace_fine'  ,'nullspace_fine'  ,S1p,S2p,S3p,N1p,N2p,N3p,0,          psi)
                else
                    !--- Lösen ---
                    res = 0.
                    psi = 1.
                    call BiCGstab(eps,5000,.false.,S1p,S2p,S3p,N1p,N2p,N3p,res,psi,problem_type,.true.,.false.,1)
           
                    !--- Eckenbehandlung ---
                    if (corner_yes) call handle_corner_Lap(g,psi)
           
                    !--- Normieren ---
                    call get_norms2(g,N1,N2,N3,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),psi,problem_type,.false.,.true.,normTwo=norm2)
                    psi = psi / SQRT(norm2)
           
                    !--- Schreiben ---
                    call write_hdf('nullspace_fine'  ,'nullspace_fine'  ,S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),psi)
           
                    if (rank == 0) write(10,*)
                    if (rank == 0) call flush(10)
                end if
               !=====================================================================================================
        
            else
        
                !=====================================================================================================
                if (read_nullspace_yes) then
                    !--- Lesen ---
                    call read2_1D_hdf('nullspace_fine_1','nullspace_fine_1',N1,S1p,N1p,iShift,1,line1(1,2))
                else
                    if (rank == 0) write(*,*) ' direction 1:'
           
                    !--- Lösen ---
                    line1(:,1) = 0.
                    line1(:,2) = 1.
                    call BiCGstab2_1D(eps,10*M1,1,b1L,b1U,N1,line1(b1L,1),line1(b1L,2),.true.,.false.,0) ! Achtung: Bei exakter Arithmetik muesste der Loeser die exakte Loesung nach spaetestens M1 Iterationen finden!
           
                    !-- Normieren (soll under-/overflow vermeiden) ---
                    norm2 = 0.
                    do i = SNB(1,1,g), SNB(2,1,g)
                        norm2 = norm2 + line1(i,2)**2
                    end do
                    call MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
                    line1 = line1 / SQRT(norm2_global)
           
                    !--- Schreiben ---
                    call write_1D_hdf('nullspace_fine_1','nullspace_fine_1',N1,S1p,N1p,iShift,1,line1(1,2))
           
                    if (rank == 0) write(10,*)
                    if (rank == 0) call flush(10)
                end if
                !-----------------------------------------------------------------------------------------------------
                if (read_nullspace_yes) then
                    !--- Lesen ---
                    call read2_1D_hdf('nullspace_fine_2','nullspace_fine_2',N2,S2p,N2p,jShift,2,line2(1,2))
                else
                    if (rank == 0) write(*,*) ' direction 2:'
           
                    !--- Lösen ---
                    line2(:,1) = 0.
                    line2(:,2) = 1.
                    call BiCGstab2_1D(eps,10*M2,2,b2L,b2U,N2,line2(b2L,1),line2(b2L,2),.true.,.false.,0)
           
                    !-- Normieren (soll under-/overflow vermeiden) ---
                    norm2 = 0.
                    do j = SNB(1,2,g), SNB(2,2,g)
                        norm2 = norm2 + line2(j,2)**2
                    end do
                    call MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR2,merror)
                    line2 = line2 / SQRT(norm2_global)
           
                    !--- Schreiben ---
                    call write_1D_hdf('nullspace_fine_2','nullspace_fine_2',N2,S2p,N2p,jShift,2,line2(1,2))
           
                    if (rank == 0) write(10,*)
                    if (rank == 0) call flush(10)
                end if
                !-----------------------------------------------------------------------------------------------------
                if (dimens == 3) then
                    if (read_nullspace_yes) then
                        !--- Lesen ---
                        call read2_1D_hdf('nullspace_fine_3','nullspace_fine_3',N3,S3p,N3p,kShift,3,line3(1,2))
                    else
                        if (rank == 0) write(*,*) ' direction 3:'
           
                        !--- L�sen ---
                        line3(:,1) = 0.
                        line3(:,2) = 1.
                        call BiCGstab2_1D(eps,10*M3,3,b3L,b3U,N3,line3(b3L,1),line3(b3L,2),.true.,.false.,0)
           
                        !-- Normieren (soll under-/overflow vermeiden) ---
                        norm2 = 0.
                        do k = SNB(1,3,g), SNB(2,3,g)
                            norm2 = norm2 + line3(k,2)**2
                        end do
                        call MPI_ALLREDUCE(norm2,norm2_global,1,MPI_REAL8,MPI_SUM,COMM_BAR3,merror)
                        line3 = line3 / SQRT(norm2_global)
           
                        !--- Schreiben ---
                        call write_1D_hdf('nullspace_fine_3','nullspace_fine_3',N3,S3p,N3p,kShift,3,line3(1,2))
           
                        if (rank == 0) write(10,*)
                        if (rank == 0) call flush(10)
                    end if
                end if
                !=====================================================================================================
                !--- "Dyadisches" Produkt bilden ---
                psi = 0.
                if (dimens == 3) then
                    do k = SNB(1,3,g), SNB(2,3,g)
                        do j = SNB(1,2,g), SNB(2,2,g)
                            do i = SNB(1,1,g), SNB(2,1,g)
                                psi(i,j,k) = line1(i,2)*line2(j,2)*line3(k,2)
                            end do
                        end do
                    end do
                else
                    do k = SNB(1,3,g), SNB(2,3,g)
                        do j = SNB(1,2,g), SNB(2,2,g)
                            do i = SNB(1,1,g), SNB(2,1,g)
                                psi(i,j,k) = line1(i,2)*line2(j,2)
                            end do
                        end do
                    end do
                end if
        
                !--- Eckenbehandlung ---
                if (corner_yes) call handle_corner_Lap(g,psi)
               !=====================================================================================================
            end if
     
        else if (n_grids >= g) then
     
            if (read_nullspace_yes .and. g == 1) then
                !--- Lesen ---
                call read2_hdf('nullspace_coarse','nullspace_coarse',S1p,S2p,S3p,N1p,N2p,N3p,0,psi)
            else
                !--- L�sen ---
                res = 0.
                psi = 1.
                call BiCGstab2(eps,200,.false.,NN(1,g),NN(2,g),NN(3,g),g,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),res,psi,problem_type,.false.,.true.,1)
        
                !--- Eckenbehandlung ---
                if (corner_yes) call handle_corner_Lap(g,psi)
        
                !--- Normieren ---
                call get_norms2(g,N1,N2,N3,SNB(1,1,g),SNB(1,2,g),SNB(1,3,g),SNB(2,1,g),SNB(2,2,g),SNB(2,3,g),psi,problem_type,.false.,.true.,normTwo=norm2)
                psi = psi / SQRT(norm2)
        
                !--- Schreiben ---
                if (g == 1) call write_hdf('nullspace_coarse','nullspace_coarse',S1p,S2p,S3p,N1p,N2p,N3p,0,(/1,1,1/),psi)
        
                if (rank == 0) write(10,*)
                if (rank == 0) call flush(10)
            end if
     
        end if
  
  
        n_grids = n_grids_orig ! TEST!!!
  
  
    end subroutine solve_nullspace
  
  
  
  
  
  
  
  
  
  
    ! TEST!!! Verteilen auf andere Files, falls es sich bewaehrt ...
    !***********************************************************************
    !***********************************************************************
    subroutine BiCGstab2_1D(eps,n_it_max,m,bbL,bbU,Nmax,bb,phi,quiet_yes1,quiet_yes2,preconditioner)
  
        implicit none
  
        real   , intent(in)    ::  eps
        integer, intent(in)    ::  n_it_max
  
        integer, intent(in   ) ::  m
        integer, intent(in   ) ::  bbL, bbU, Nmax
  
        integer, intent(in)    ::  preconditioner
  
        real   , intent(in)    ::  bb (bbL:(Nmax+bbU))
        real   , intent(inout) ::  phi(bbL:(Nmax+bbU))
  
        real                   ::  pp (bbL:(Nmax+bbU))
        real                   ::  Ap (bbL:(Nmax+bbU))
        real                   ::  rr (bbL:(Nmax+bbU))
        real                   ::  rh (bbL:(Nmax+bbU))
        real                   ::  Ar (bbL:(Nmax+bbU))
        real                   ::  z1 (bbL:(Nmax+bbU))
        real                   ::  z2 (bbL:(Nmax+bbU))
  
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
        integer                ::  i, counter
        integer                ::  SS, NN
  
        integer                ::  comm
  
        real                   ::  norm_inf, norm_inf_global
        real                   ::  rhr, rhr_prev, ArAr, rAr, rhAp, rhr_global, rhAp_global
        real                   ::  scalar_global(1:2)
        real                   ::  alpha, beta, omega
  
        logical                ::  exit_yes
  
  
        !===========================================================================================================
        if      (m == 1) then
            SS = S1p
            NN = N1p
            comm = COMM_BAR1
        else if (m == 2) then
            SS = S2p
            NN = N2p
            comm = COMM_BAR2
        else if (m == 3) then
            SS = S3p
            NN = N3p
            comm = COMM_BAR3
        else
            if (rank == 0) write(*,*) 'ERROR! m /= 1,2,3!'
            call MPI_FINALIZE(merror)
            stop
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        call product_div_grad_transp_1D(m,bbL,bbU,Nmax,phi,rr)
  
        rr(SS:NN) = bb(SS:NN) - rr(SS:NN)
        rh(SS:NN) = rr(SS:NN)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Residuum ==============================================================================================
        !===========================================================================================================
        norm_inf = 0.
        do i = SS, NN
            norm_inf = MAX(ABS(rr(i)),norm_inf)
        end do
  
        call MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
        norm_inf = norm_inf_global
        !===========================================================================================================
  
        counter = 0
  
  
        ITERATE: do
     
            !========================================================================================================
            !=== überprüfen des Konvergenzkriteriums ================================================================
            !========================================================================================================
            call status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
            if (exit_yes) exit ITERATE
            !========================================================================================================
     
     
            !========================================================================================================
            !=== nächster Durchlauf =================================================================================
            !========================================================================================================
            counter = counter + 1
            !========================================================================================================
     
     
            !========================================================================================================
            rhr_prev = rhr
            rhr = 0.
            do i = SS, NN
                rhr = rhr + rr(i)*rh(i)
            end do
     
            call MPI_ALLREDUCE(rhr,rhr_global,1,MPI_REAL8,MPI_SUM,comm,merror)
            rhr = rhr_global
            !========================================================================================================
            if (ABS(rhr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhr =', rhr
                if (rank == 0) write(10,'(a,E13.5)') 'rhr =', rhr
        
                rh(SS:NN) = rr(SS:NN) ! Neuer Referenzvektor ...
                rhr = 0.
                do i = SS, NN
                    rhr = rhr + rr(i)*rh(i)
                end do
        
                call MPI_ALLREDUCE(rhr,rhr_global,1,MPI_REAL8,MPI_SUM,comm,merror)
                rhr = rhr_global
            end if
            !========================================================================================================
            if (counter >= 2) then
                if (omega == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'omega =', omega
                    if (rank == 0) write(10,'(a,E13.5)') 'omega =', omega
                    exit ITERATE
                end if
                if (rhr_prev == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    if (rank == 0) write(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    exit ITERATE
                end if
                beta = (alpha/omega)*(rhr/rhr_prev)
                omega = -beta*omega
                pp(SS:NN) = rr(SS:NN) + beta*pp(SS:NN) + omega*Ap(SS:NN)
            else
                pp(SS:NN) = rr(SS:NN)
            end if
            !========================================================================================================
            if (preconditioner == 0) then
                call product_div_grad_transp_1D(m,bbL,bbU,Nmax,pp,Ap)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                !IF (preconditioner == 1) CALL multigridV(.TRUE.,g,pp,z1,problem_type)
                !IF (preconditioner == 2) CALL multigridF(.TRUE.,g,pp,z1,problem_type)
                call product_div_grad_transp_1D(m,bbL,bbU,Nmax,z1,Ap)
            else
                if (rank == 0) write(*,'(a)') 'ERROR! Specify valid preconditioner!'
                call MPI_FINALIZE(merror)
                stop
            end if
            !========================================================================================================
            rhAp = 0.
            do i = SS, NN
                rhAp = rhAp + rh(i)*Ap(i)
            end do
     
            call MPI_ALLREDUCE(rhAp,rhAp_global,1,MPI_REAL8,MPI_SUM,comm,merror)
            rhAp = rhAp_global
            !========================================================================================================
            if (ABS(rhAp) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhAp =', rhAp
                if (rank == 0) write(10,'(a,E13.5)') 'rhAp =', rhAp
                exit ITERATE
                alpha = 0.
            else
                alpha = rhr / rhAp
            end if
            !========================================================================================================
            rr(SS:NN) = rr(SS:NN) - alpha*Ap(SS:NN)
            !========================================================================================================
            if (preconditioner == 0) then
                call product_div_grad_transp_1D(m,bbL,bbU,Nmax,rr,Ar)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                !IF (preconditioner == 1) CALL multigridV(.TRUE.,g,rr,z2,problem_type)
                !IF (preconditioner == 2) CALL multigridF(.TRUE.,g,rr,z2,problem_type)
                call product_div_grad_transp_1D(m,bbL,bbU,Nmax,z2,Ar)
            end if
            !========================================================================================================
            rAr  = 0.
            ArAr = 0.
            do i = SS, NN
                rAr  = rAr  + rr(i)*Ar(i)
                ArAr = ArAr + Ar(i)*Ar(i)
            end do
     
            call MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,comm,merror)
            rAr  = scalar_global(1)
            ArAr = scalar_global(2)
            !========================================================================================================
            if (ABS(rAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rAr =', rAr
                if (rank == 0) write(10,'(a,E13.5)') 'rAr =', rAr
            end if
            if (ABS(ArAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'ArAr =', ArAr
                if (rank == 0) write(10,'(a,E13.5)') 'ArAr =', ArAr
                exit ITERATE
                omega = 0.
            else
                omega = rAr / ArAr
            end if
            !========================================================================================================
            norm_inf = 0.
            !--------------------------------------------------------------------------------------------------------
            if (preconditioner == 0) then
                do i = SS, NN
                    phi(i) = phi(i) + omega*rr(i) + alpha*pp(i)
                    rr (i) = rr (i) - omega*Ar(i)
                    norm_inf = MAX(ABS(rr(i)),norm_inf)
                end do
            else
                do i = SS, NN
                    phi(i) = phi(i) + omega*z2(i) + alpha*z1(i)
                    rr (i) = rr (i) - omega*Ar(i)
                    norm_inf = MAX(ABS(rr(i)),norm_inf)
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            call MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            norm_inf = norm_inf_global
           !========================================================================================================
     
        end do ITERATE
  
  
    end subroutine BiCGstab2_1D
  
  
  
  
  

 
  
  
  
  
    subroutine product_div_grad_transp_1D(m,bbL,bbU,Nmax,phi,Lap) ! TEST!!! woanders hin ...
  
        implicit none
  
        integer, intent(in   ) ::  m
        integer, intent(in   ) ::  bbL, bbU, Nmax
        real   , intent(inout) ::  phi (bbL:(Nmax+bbU))
        real   , intent(out  ) ::  Lap (bbL:(Nmax+bbU))
        real                   ::  work(bbL:(Nmax+bbU))
  
  
        call divergence_transp_1D      (m,bbL,bbU,Nmax,phi,work)
        call bc_extrapolation_transp_1D(m,bbL,bbU,Nmax,work    )
        call gradient_transp_1D        (m,bbL,bbU,Nmax,work,Lap)
  
  
    end subroutine product_div_grad_transp_1D
  
  
  
  
  
  
  
  
  
  
  
    subroutine bc_extrapolation_transp_1D(m,bbL,bbU,Nmax,phi) ! TEST!!! woanders hin ...
  
        implicit none
  
        integer, intent(in   ) ::  m
        integer, intent(in   ) ::  bbL, bbU, Nmax
        real   , intent(inout) ::  phi(bbL:(Nmax+bbU))
  
        integer                ::  i, ii
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (m == 1) then
            if (BC_1L > 0) then
                i = 0
                do ii = 0, d1U
                    phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIup(ii,1)/cIup(-1,1)
                end do
                phi(i) = phi(i) / cIup(-1,1)
            end if
            if (BC_1U > 0) then
                i = N1
                do ii = d1L, -1
                    phi(i+ii) = phi(i+ii) - phi(i)*cIup(ii,i)/cIup(0,i)
                end do
                phi(i) = phi(i) / cIup(0,i)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 2) then
            if (BC_2L > 0) then
                i = 0
                do ii = 0, d2U
                    phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIvp(ii,1)/cIvp(-1,1)
                end do
                phi(i) = phi(i) / cIvp(-1,1)
            end if
            if (BC_2U > 0) then
                i = N2
                do ii = d2L, -1
                    phi(i+ii) = phi(i+ii) - phi(i)*cIvp(ii,i)/cIvp(0,i)
                end do
                phi(i) = phi(i) / cIvp(0,i)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (m == 3) then
            if (BC_3L > 0) then
                i = 0
                do ii = 0, d3U
                    phi(i+ii+1) = phi(i+ii+1) - phi(i)*cIwp(ii,1)/cIwp(-1,1)
                end do
                phi(i) = phi(i) / cIwp(-1,1)
            end if
            if (BC_3U > 0) then
                i = N3
                do ii = d3L, -1
                    phi(i+ii) = phi(i+ii) - phi(i)*cIwp(ii,i)/cIwp(0,i)
                end do
                phi(i) = phi(i) / cIwp(0,i)
            end if
        end if
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine bc_extrapolation_transp_1D
  
  
  
  
  
  
  
  
  
  
  
    subroutine divergence_transp_1D(m,bbL,bbU,Nmax,phi,div) ! TEST!!! woanders hin ...
  
        implicit none
  
        integer, intent(in   ) ::  m
        integer, intent(in   ) ::  bbL, bbU, Nmax
        real   , intent(inout) ::  phi(bbL:(Nmax+bbU))
        real   , intent(inout) ::  div(bbL:(Nmax+bbU))
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
  
        !===========================================================================================================
        if (m == 1) then
            !--------------------------------------------------------------------------------------------------------
            if (comp_div_yes) then
        
                j = S2p
                k = S3p
        
                do i = S1p, N1p
                    com(i,j,k) = phi(i)*dx1DM(i)
                end do
        
                call apply_compact_transp(1,0,S1p,j,k,N1p,j,k,S11B,j,k,N11B,j,k,N1,ndL,ndR,dimS1,cDu1CLT,cDu1CLT_LU,cDu1CRT,WDu1T,SDu1T,com,psi) ! ACHTUNG!! "psi"
        
                do i = S11B, N11B
                    div(i) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
                do i = S11B, N11B
                    div(i) = cDu1T(g1L,i)*phi(i+g1L)
                    do ii = g1L+1, g1U
                        div(i) = div(i) + cDu1T(ii,i)*phi(i+ii)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
        end if
        !===========================================================================================================
        if (m == 2) then
            !--------------------------------------------------------------------------------------------------------
            if (comp_div_yes) then
        
                i = S1p
                k = S3p
        
                do j = S2p, N2p
                    com(i,j,k) = phi(j)*dx2DM(j)
                end do
        
                call apply_compact_transp(2,0,i,S2p,k,i,N2p,k,i,S22B,k,i,N22B,k,N2,ndL,ndR,dimS2,cDv2CLT,cDv2CLT_LU,cDv2CRT,WDv2T,SDv2T,com,psi)
        
                do j = S22B, N22B
                    div(j) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
                do j = S22B, N22B
                    div(j) = cDv2T(g2L,j)*phi(j+g2L)
                    do jj = g2L+1, g2U
                        div(j) = div(j) + cDv2T(jj,j)*phi(j+jj)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
        end if
        !===========================================================================================================
        if (m == 3) then
            !--------------------------------------------------------------------------------------------------------
            if (comp_div_yes) then
        
                i = S1p
                j = S2p
        
                do k = S3p, N3p
                    com(i,j,k) = phi(k)*dx3DM(k)
                end do
        
                call apply_compact_transp(3,0,i,j,S3p,i,j,N3p,i,j,S33B,i,j,N33B,N3,ndL,ndR,dimS3,cDw3CLT,cDw3CLT_LU,cDw3CRT,WDw3T,SDw3T,com,psi)
        
                do k = S33B, N33B
                    div(k) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,0,bbL,bbU,Nmax,phi)
        
                do k = S33B, N33B
                    div(k) = cDw3T(g3L,k)*phi(k+g3L)
                    do kk = g3L+1, g3U
                        div(k) = div(k) + cDw3T(kk,k)*phi(k+kk)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
        end if
    !===========================================================================================================
  
  
    end subroutine divergence_transp_1D
  
  
  
  
  
  
  
  
  
  
  
    subroutine gradient_transp_1D(m,bbL,bbU,Nmax,phi,grad) ! TEST!!! woanders hin ...
  
        implicit none
  
        integer, intent(in   ) ::  m
        integer, intent(in   ) ::  bbL, bbU, Nmax
        real   , intent(inout) ::  phi (bbL:(Nmax+bbU))
        real   , intent(out  ) ::  grad(bbL:(Nmax+bbU))
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
  
        !===========================================================================================================
        if (m == 1) then
     
            !--- Randbedingungen ------------------------------------------------------------------------------------
            if (BC_1L > 0) phi(0 ) = 0.
            if (BC_1U > 0) phi(N1) = 0.
     
            !--------------------------------------------------------------------------------------------------------
            if (comp_grad_yes) then
        
                j = S2p
                k = S3p
        
                if (BC_1L > 0) com(0 ,j,k) = 0.
                if (BC_1U > 0) com(N1,j,k) = 0.
        
                do i = S11, N11
                    com(i,j,k) = phi(i)*dx1GM(i)
                end do
        
                call apply_compact_transp(1,1,S11,j,k,N11,j,k,S1p,j,k,N1p,j,k,N1,ndL,ndR,dimS1,cGp1CLT,cGp1CLT_LU,cGp1CRT,WGp1T,SGp1T,com,psi)
        
                do i = S1p, N1p
                    grad(i) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
                do i = S1p, N1p
                    grad(i) = cGp1T(d1L,i)*phi(i+d1L)
                    do ii = d1L+1, d1U
                        grad(i) = grad(i) + cGp1T(ii,i)*phi(i+ii)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
     
        end if
        !===========================================================================================================
        if (m == 2) then
     
            !--- Randbedingungen ------------------------------------------------------------------------------------
            if (BC_2L > 0) phi(0 ) = 0.
            if (BC_2U > 0) phi(N2) = 0.
     
            !--------------------------------------------------------------------------------------------------------
            if (comp_grad_yes) then
        
                i = S1p
                k = S3p
        
                if (BC_2L > 0) com(i,0 ,k) = 0.
                if (BC_2U > 0) com(i,N2,k) = 0.
        
                do j = S22, N22
                    com(i,j,k) = phi(j)*dx2GM(j)
                end do
        
                call apply_compact_transp(2,2,i,S22,k,i,N22,k,i,S2p,k,i,N2p,k,N2,ndL,ndR,dimS2,cGp2CLT,cGp2CLT_LU,cGp2CRT,WGp2T,SGp2T,com,psi)
        
                do j = S2p, N2p
                    grad(j) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
                do j = S2p, N2p
                    grad(j) = cGp2T(d2L,j)*phi(j+d2L)
                    do jj = d2L+1, d2U
                        grad(j) = grad(j) + cGp2T(jj,j)*phi(j+jj)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
     
        end if
        !===========================================================================================================
        if (m == 3) then
     
            !--- Randbedingungen ------------------------------------------------------------------------------------
            if (BC_3L > 0) phi(0 ) = 0.
            if (BC_3U > 0) phi(N3) = 0.
     
            !--------------------------------------------------------------------------------------------------------
            if (comp_grad_yes) then
        
                i = S1p
                j = S2p
        
                if (BC_3L > 0) com(i,j,0 ) = 0.
                if (BC_3U > 0) com(i,j,N3) = 0.
        
                do k = S33, N33
                    com(i,j,k) = phi(k)*dx3GM(k)
                end do
        
                call apply_compact_transp(3,3,i,j,S33,i,j,N33,i,j,S3p,i,j,N3p,N3,ndL,ndR,dimS3,cGp3CLT,cGp3CLT_LU,cGp3CRT,WGp3T,SGp3T,com,psi)
        
                do k = S3p, N3p
                    grad(k) = psi(i,j,k)
                end do
            !--------------------------------------------------------------------------------------------------------
            else
                call exchange_1D(m,m,bbL,bbU,Nmax,phi) ! Muss nach Randbedingungen kommen!
        
                do k = S3p, N3p
                    grad(k) = cGp3T(d3L,k)*phi(k+d3L)
                    do kk = d3L+1, d3U
                        grad(k) = grad(k) + cGp3T(kk,k)*phi(k+kk)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
     
        end if
    !===========================================================================================================
  
  
    end subroutine gradient_transp_1D
  
  
  
  
  
  
  
  
  
  
  
    subroutine exchange_1D(dir,vel_dir,bbL,bbU,Nmax,phi) ! TEST!!! woanders hin ...
  
        implicit none
  
        integer, intent(in   ) ::  dir
        integer, intent(in   ) ::  vel_dir
        integer, intent(in   ) ::  bbL, bbU, Nmax
  
        real   , intent(inout) ::  phi(bbL:(Nmax+bbU))
  
        integer                ::  status(MPI_STATUS_SIZE) ! TEST!!! Warum steht das eigentlich hier lokal? gleiches gilt auch fuer die anderen Module ...
  
        real                   ::  ghostLR( 1:bbU)
        real                   ::  ghostLS( 1:bbU)
        real                   ::  ghostUR(bbL:-1)
        real                   ::  ghostUS(bbL:-1)
  
        integer                ::  lengthL, lengthU
        integer                ::  rankL, rankU
        integer                ::  BCL, BCU
        integer                ::  SS, NN
  
  
        lengthL = ABS(bbL)
        lengthU = ABS(bbU)
  
        if      (dir == 1) then
            SS    = S1p
            NN    = N1p
            BCL   = BC_1L
            BCU   = BC_1U
            rankL = rank1L
            rankU = rank1U
        else if (dir == 2) then
            SS    = S2p
            NN    = N2p
            BCL   = BC_2L
            BCU   = BC_2U
            rankL = rank2L
            rankU = rank2U
        else if (dir == 3) then
            SS    = S3p
            NN    = N3p
            BCL   = BC_3L
            BCU   = BC_3U
            rankL = rank3L
            rankU = rank3U
        end if
  
  
        !======================================================================================================================
        if (BCL == 0) call MPI_IRECV(ghostUR,lengthL,MPI_REAL8,rankL,1,COMM_CART,req1L,merror) ! TEST!!! Was ist mit req1L, etc.?
        if (BCU == 0) call MPI_IRECV(ghostLR,lengthU,MPI_REAL8,rankU,2,COMM_CART,req1U,merror)
  
        if (BCU == 0) ghostUS(bbL:-1) = phi((NN+1+bbL):NN)
        if (BCL == 0) ghostLS(1:bbU ) = phi(SS:(SS-1+bbU))
  
        if (BCU == 0) call MPI_SEND(ghostUS,lengthL,MPI_REAL8,rankU,1,COMM_CART,merror)
        if (BCL == 0) call MPI_SEND(ghostLS,lengthU,MPI_REAL8,rankL,2,COMM_CART,merror)
  
        if (BCL == 0) call MPI_WAIT(req1L,status,merror)
        if (BCU == 0) call MPI_WAIT(req1U,status,merror)
  
        if (BCL == 0) call pseudocall(ghostUR) ! Soll den Compiler daran hindern, das Umspeichern mit MPI_WAIT zu vertauschen.
        if (BCU == 0) call pseudocall(ghostLR)
  
        if (BCL == 0) phi((SS+bbL):(SS-1)) = ghostUR(bbL:-1)
        if (BCU == 0) phi((NN+1):(NN+bbU)) = ghostLR(1:bbU )
        !-------------------------------------------------------------------------------------------------------------------
        if (BCL == -1) then
            phi((SS+bbL):(SS-1)) = phi((NN+1+bbL):NN)
            phi((NN+1):(NN+bbU)) = phi(SS:(SS-1+bbU))
        end if
        !-------------------------------------------------------------------------------------------------------------------
        if (vel_dir == dir) then
            if (BCL > 0) phi( bbL    : -1       ) = 0.
            if (BCU > 0) phi((Nmax+1):(Nmax+bbU)) = 0.
     
            if (BCL  == -2) phi( bbL    :  0       ) = 0.
            if (BCU  == -2) phi( Nmax   :(Nmax+bbU)) = 0.
        else
            if (BCL > 0 .or. BCL == -2) phi( bbL    : 0        ) = 0.
            if (BCU > 0 .or. BCU == -2) phi((Nmax+1):(Nmax+bbU)) = 0.
        end if
    !======================================================================================================================
  
  
    end subroutine exchange_1D
  
  
  
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
    subroutine solve_Helmholtz( &
        m,                      &
        mulL,                   &
        epsU,                   &
        n_it_max,               &
!        init_yes,               &
        bb,                     &
        phi,                    &
        cquiet_yes1,             &
        cquiet_yes2 ) bind ( c, name='OP_SolveHelmholtz' )
  
        implicit none
  
        integer(c_int), intent(in   ) ::  m
  
        !        real(c_double), intent(in   ) ::  multI
        real(c_double), intent(in   ) ::  mulL

        real(c_double), intent(in   ) ::  epsU
        integer(c_int), intent(in   ) ::  n_it_max
  
!        logical(c_bool), intent(in  ) ::  init_yes
  
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
  
  
    end subroutine solve_Helmholtz
  
  
  
    subroutine solve_conc(quiet_yes1,quiet_yes2)
  
        implicit none
  
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
        integer                ::  m
  
        integer                ::  i, ii, imax, iimax, di
        integer                ::  j, jj, jmax, jjmax, dj
        integer                ::  k, kk, kmax, kkmax, dk
  
        integer                ::  g
  
  
        m = conc_nu
  
        if (rank == 0 .and. log_iteration_yes                     ) write(10,'(a     )') 'solving Helmholtz equations for concentration field ...'
        if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(* ,'(a,i2,a)') 'concentration field', m, ':'
  
  
        !===========================================================================================================
        if (BCc_1L(m) == 3) then
            do k = S3p, N3p
                do j = S2p, N2p
                    velReSc1(j,k,1,1) = usReSc(1,m) + bc11(j,k,1)
                end do
            end do
        end if
        if (BCc_1U(m) == 3) then
            do k = S3p, N3p
                do j = S2p, N2p
                    velReSc1(j,k,2,1) = usReSc(1,m) + bc11(j,k,2)
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
            do k = S3p, N3p
                do i = S1p, N1p
                    velReSc2(i,k,1,1) = usReSc(2,m) + bc22(i,k,1)
                end do
            end do
        end if
        if (BCc_2U(m) == 3) then
            do k = S3p, N3p
                do i = S1p, N1p
                    velReSc2(i,k,2,1) = usReSc(2,m) + bc22(i,k,2)
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
            do j = S2p, N2p
                do i = S1p, N1p
                    velReSc3(i,j,1,1) = usReSc(3,m) + bc33(i,j,1)
                end do
            end do
        end if
        if (BCc_3U(m) == 3) then
            do j = S2p, N2p
                do i = S1p, N1p
                    velReSc3(i,j,2,1) = usReSc(3,m) + bc33(i,j,2)
                end do
            end do
        end if
        !===========================================================================================================
  
  
  
        !===========================================================================================================
        do g = 1, n_grids-1

            if (N1 == 1) then
                imax  = 1
                iimax = 1
            else
                imax  = NN(1,g)
                iimax = NN(1,g+1)
            end if

            if (N2 == 1) then
                jmax  = 1
                jjmax = 1
            else
                jmax  = NN(2,g)
                jjmax = NN(2,g+1)
            end if

            if (N3 == 1) then
                kmax  = 1
                kkmax = 1
            else
                kmax  = NN(3,g)
                kkmax = NN(3,g+1)
            end if
     
            di = (imax-1)/(iimax-1)
            dj = (jmax-1)/(jjmax-1)
            dk = (kmax-1)/(kkmax-1)
     
            !--------------------------------------------------------------------------------------------------------
            if (BCc_1L(m) == 3) then
                do kk = 1, kkmax
                    k = 1+dk*(kk-1)
                    !pgi$ unroll = n:8
                    do jj = 1, jjmax
                        j = 1+dj*(jj-1)
                        velReSc1(jj,kk,1,g+1) = velReSc1(j,k,1,g)
                    end do
                end do
            end if
            if (BCc_1U(m) == 3) then
                do kk = 1, kkmax
                    k = 1+dk*(kk-1)
                    !pgi$ unroll = n:8
                    do jj = 1, jjmax
                        j = 1+dj*(jj-1)
                        velReSc1(jj,kk,2,g+1) = velReSc1(j,k,2,g)
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (BCc_2L(m) == 3) then
                do kk = 1, kkmax
                    k = 1+dk*(kk-1)
                    !pgi$ unroll = n:8
                    do ii = 1, iimax
                        i = 1+di*(ii-1)
                        velReSc2(ii,kk,1,g+1) = velReSc2(i,k,1,g)
                    end do
                end do
            end if
            if (BCc_2U(m) == 3) then
                do kk = 1, kkmax
                    k = 1+dk*(kk-1)
                    !pgi$ unroll = n:8
                    do ii = 1, iimax
                        i = 1+di*(ii-1)
                        velReSc2(ii,kk,2,g+1) = velReSc2(i,k,2,g)
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (BCc_3L(m) == 3) then
                do jj = 1, jjmax
                    j = 1+dj*(jj-1)
                    !pgi$ unroll = n:8
                    do ii = 1, iimax
                        i = 1+di*(ii-1)
                        velReSc3(ii,jj,1,g+1) = velReSc3(i,j,1,g)
                    end do
                end do
            end if
            if (BCc_3U(m) == 3) then
                do jj = 1, jjmax
                    j = 1+dj*(jj-1)
                    !pgi$ unroll = n:8
                    do ii = 1, iimax
                        i = 1+di*(ii-1)
                        velReSc3(ii,jj,2,g+1) = velReSc3(i,j,2,g)
                    end do
                end do
            end if
           !--------------------------------------------------------------------------------------------------------
     
        end do
        !===========================================================================================================

        if (corner_yes) call handle_corner_rhs_conc(res)

        call BiCGstab(epsU,n_it_Helmh_conc,init_conc(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,conc(b1L,b2L,b3L,m),3,quiet_yes1,quiet_yes2,precond_Helmh_conc)
        !CALL Richardson(epsU,n_it_Helmh_conc,init_conc(substep),S1p,S2p,S3p,N1p,N2p,N3p,res,conc(b1L,b2L,b3L,m),3,quiet_yes1,quiet_yes2,precond_Helmh_conc)

        if (rank == 0 .and. log_iteration_yes) write(10,'(a)')
  
  
    end subroutine solve_conc
  
  
  
  
  
  
  
  
  
  
  
    subroutine solve_conc_explicit
  
        implicit none
  
        real                   ::  velReSc
  
        integer                ::  m
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
  
        m = conc_nu
  
        ! ACHTUNG!!! Kopieren ist hier etwas ungl�cklich, erh�lt daf�r aber die �bersichtliche Struktur in anderen Routinen:
        conc(S1p:N1p,S2p:N2p,S3p:N3p,m) = res(S1p:N1p,S2p:N2p,S3p:N3p)
  
  
        !===========================================================================================================
        if (BCc_1L(m) == 3) then
            i = 1
            do k = S3p, N3p
                do j = S2p, N2p
                    velReSc = usReSc(1,m) + bc11(j,k,1)
                    !pgi$ unroll = n:8
                    do ii = 1, b1U
                        conc(i,j,k,m) = conc(i,j,k,m) - cc1(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc1(0,i,m) - velReSc)
                end do
            end do
        end if
        if (BCc_1U(m) == 3) then
            i = N1
            do k = S3p, N3p
                do j = S2p, N2p
                    velReSc = usReSc(1,m) + bc11(j,k,2)
                    !pgi$ unroll = n:8
                    do ii = b1L, -1
                        conc(i,j,k,m) = conc(i,j,k,m) - cc1(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc1(0,i,m) - velReSc)
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
            j = 1
            do k = S3p, N3p
                do i = S1p, N1p
                    velReSc = usReSc(2,m) + bc22(i,k,1)
                    !pgi$ unroll = n:8
                    do jj = 1, b2U
                        conc(i,j,k,m) = conc(i,j,k,m) - cc2(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc2(0,j,m) - velReSc)
                end do
            end do
        end if
        if (BCc_2U(m) == 3) then
            j = N2
            do k = S3p, N3p
                do i = S1p, N1p
                    velReSc = usReSc(2,m) + bc22(i,k,2)
                    !pgi$ unroll = n:8
                    do jj = b2L, -1
                        conc(i,j,k,m) = conc(i,j,k,m) - cc2(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc2(0,j,m) - velReSc)
                end do
            end do
        end if
        !-----------------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
            k = 1
            do j = S2p, N2p
                do i = S1p, N1p
                    velReSc = usReSc(3,m) + bc33(i,j,1)
                    !pgi$ unroll = n:8
                    do kk = 1, b3U
                        conc(i,j,k,m) = conc(i,j,k,m) - cc3(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc3(0,k,m) - velReSc)
                end do
            end do
        end if
        if (BCc_3U(m) == 3) then
            k = N3
            do j = S2p, N2p
                do i = S1p, N1p
                    velReSc = usReSc(3,m) + bc33(i,j,2)
                    !pgi$ unroll = n:8
                    do kk = b3L, -1
                        conc(i,j,k,m) = conc(i,j,k,m) - cc3(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                    conc(i,j,k,m) = conc(i,j,k,m) / (cc3(0,k,m) - velReSc)
                end do
            end do
        end if
        !===========================================================================================================
  
  
        if (corner_yes) call handle_corner_rhs_conc(conc(b1L,b2L,b3L,m))
  
  
    end subroutine solve_conc_explicit
  
  
  
  
  
  
  
  
  
  
  
    subroutine relax_restrict(init_yes,nullspace_yes,g,psi,bb,phi,work,coarse,problem_type)
  
        implicit none
  
        logical, intent(in   ) ::  init_yes
        logical, intent(in   ) ::  nullspace_yes
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  psi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
            if (nullspace_yes) call apply_nullspace2(g,psi,bb,problem_type)
            call relaxation_div_grad(init_yes,n_relax_down,g,bb,phi)
            call product_div_grad_relax(g,phi,work)
            call restrict(.true.,g,coarse,bb,work)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 1) then
            if (g == 1) then
                call relaxation_Helmholtz(init_yes,n_relax_down,bb,phi)
                call product_Helmholtz_relax(phi,work)
                call restrict_Helmholtz(coarse,bb,work)
            else
                call relaxation_Helmholtz_coarse(init_yes,n_relax_down,g,bb,phi)
                call product_Helmholtz_relax_coarse(g,phi,work)
                call restrict(.true.,g,coarse,bb,work)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 3) then
            call relaxation_Helmholtz_conc(init_yes,n_relax_down,g,bb,phi)
            call product_Helmholtz_conc_relax(g,phi,work)
            call restrict(.true.,g,coarse,bb   ,work)
        end if
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine relax_restrict
  
  
  
  
  
  
  
  
  
  
  
    subroutine interpolate_relax(g,bb,phi,work,coarse,problem_type)
  
        implicit none
  
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
            call interpolate(.true.,g,coarse,phi,work)
            call relaxation_div_grad_inv(.false.,n_relax_up,g,bb,phi)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 1) then
            if (g == 1) then
                call interpolate_Helmholtz(coarse,phi,work)
                call relaxation_Helmholtz       (.false.,n_relax_up,  bb,phi)
            else
                call interpolate(.true.,g,coarse,phi,work)
                call relaxation_Helmholtz_coarse(.false.,n_relax_up,g,bb,phi)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 3) then
            call interpolate(.true.,g,coarse,phi  ,work)
            call relaxation_Helmholtz_conc(.false.,n_relax_up,g,bb,phi)
        end if
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine interpolate_relax
  
  
  
  
  
  
  
  
  
  
  
    subroutine relax_bottom(init_yes,nullspace_yes,g,psi,bb,phi,problem_type)
  
        implicit none
  
        logical, intent(in   ) ::  init_yes
        logical, intent(in   ) ::  nullspace_yes
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  psi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(inout) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
            if (nullspace_yes) call apply_nullspace2(g,psi,bb,problem_type)
            call relaxation_div_grad    (init_yes,n_relax_bottom,g,bb,phi)
            call relaxation_div_grad_inv(.false. ,n_relax_bottom,g,bb,phi)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 1) then
            if (g  ==  1) call relaxation_Helmholtz       (init_yes,n_relax_bottom,  bb,phi)
            if (g > 1) call relaxation_Helmholtz_coarse(init_yes,n_relax_bottom,g,bb,phi)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 3) call relaxation_Helmholtz_conc(init_yes,n_relax_bottom,g,bb,phi)
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine relax_bottom
  
  
  
  
  
  
  
  
  
  
  
  
    subroutine plain_restrict(nullspace_yes,g,psi,bb,work,coarse,problem_type)
  
        implicit none
  
        logical, intent(in   ) ::  nullspace_yes
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  psi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
            if (nullspace_yes) call apply_nullspace2(g,psi,bb,problem_type)
            call restrict(.false.,g,coarse,bb,work)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 1) then
            if (g == 1) then
                call restrict_Helmholtz(coarse,bb,work)
            else
                call restrict(.false.,g,coarse,bb,work)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 3) call restrict(.false.,g,coarse,bb,work)
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine plain_restrict
  
  
  
  
  
  
  
  
  
  
  
  
    subroutine interpolate_mg(g,bb,phi,work,coarse,problem_type)
  
        implicit none
  
        integer, intent(in   ) ::  g
  
        real   , intent(inout) ::  bb    (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  phi   (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
  
        integer, intent(in   ) ::  problem_type
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
            call interpolate(.false.,g,coarse,work,phi)
            call multigridV (.false.,g,bb,phi,problem_type)
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 1) then
            if (g == 1) then
                call interpolate_Helmholtz(coarse,phi,work)
               !!CALL relaxation_Helmholtz(.FALSE.,n_relax_up,  bb,phi)
               !CALL multigridV(.FALSE.,g,bb,phi,problem_type) ! TEST!!! Test schreiben, damit das nicht aufgerufen werden kann ...
            else
                call interpolate(.false.,g,coarse,work,phi)
                call multigridV (.false.,g,bb,phi,problem_type)
            end if
        end if
        !-----------------------------------------------------------------------------------------------------------
  
  
        !-----------------------------------------------------------------------------------------------------------
        if (problem_type == 3) then
            call interpolate(.false.,g,coarse,work,phi)
            call multigridV (.false.,g,bb,phi,problem_type)
        end if
    !-----------------------------------------------------------------------------------------------------------
  
  
    end subroutine interpolate_mg
  
  
  
  
  
  
  
  
  
  
  
    subroutine multigridV(init_yes,gstart,bb,phi,problem_type)
  
        implicit none
  
        logical, intent(in   ) ::  init_yes
        integer, intent(in   ) ::  gstart
        integer, intent(in   ) ::  problem_type
  
        real   , intent(inout) ::  bb (b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
        real   , intent(inout) ::  phi(b1L:(NN(1,gstart)+b1U),b2L:(NN(2,gstart)+b2U),b3L:(NN(3,gstart)+b3U))
  
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
  
  
    end subroutine multigridV
  
  
  
  
  
  
  
  
  
  
  
    subroutine multigridF(init_yes,gstart,bb,phi,problem_type)
  
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
  
  
    end subroutine multigridF
  
  
  
  
  
  
  
  
  
  
  
    ! TEST!!! Hier lassen sich evtl. Operationen einsparen! (fine1-fine2) nur einmal berechnen, ist aber fraglich, ob das auch schneller ist ...
    subroutine restrict(add_yes,g,coarse,fine1,fine2) ! TEST!!! aufräumen ...
  
        implicit none
  
        logical, intent(in   ) ::  add_yes
        integer, intent(in   ) ::  g
  
        real   , intent(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
        real   , intent(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
        integer                ::  i, ii, di, imax, iimax
        integer                ::  j, jj, dj, jmax, jjmax
        integer                ::  k, kk, dk, kmax, kkmax
  
        integer                ::  sizsg(1:3), offsg(1:3), dispg
        real   , allocatable   ::  sendbuf(:,:,:)
        real                   ::  recvbuf(1:NN(1,g+1)*NN(2,g+1)*NN(3,g+1))
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
        !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
        !                aber im Prinzip redundant (genauer: coarse(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
        !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
        !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
        !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
        !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
        !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
        !----------------------------------------------------------------------------------------------------------!
  
  
        imax    = NN(1,g)
        jmax    = NN(2,g)
        kmax    = NN(3,g)
  
        iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
        jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
        kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
        di      = (NN(1,g)-1)/(iimax-1)
        dj      = (NN(2,g)-1)/(jjmax-1)
        dk      = (NN(3,g)-1)/(kkmax-1)
  
  
        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  1      ,  1      ,1:NN(3,g)) = (fine1(2        ,1      ,1:NN(3,g)) + fine1(1      ,2        ,1:NN(3,g)))/2.
        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  1      ,  NN(2,g),1:NN(3,g)) = (fine1(2        ,NN(2,g),1:NN(3,g)) + fine1(1      ,NN(2,g)-1,1:NN(3,g)))/2.
        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  NN(1,g),  1      ,1:NN(3,g)) = (fine1(NN(1,g)-1,1      ,1:NN(3,g)) + fine1(NN(1,g),2        ,1:NN(3,g)))/2.
        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine1(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine1(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  1      ,1:NN(2,g),  1      ) = (fine1(2        ,1:NN(2,g),1      ) + fine1(1      ,1:NN(2,g),2        ))/2.
        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  1      ,1:NN(2,g),  NN(3,g)) = (fine1(2        ,1:NN(2,g),NN(3,g)) + fine1(1      ,1:NN(2,g),NN(3,g)-1))/2.
        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  1      ) = (fine1(NN(1,g)-1,1:NN(2,g),1      ) + fine1(NN(1,g),1:NN(2,g),2        ))/2.
        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine1(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine1(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  1      ,  1      ) = (fine1(1:NN(1,g),2        ,1      ) + fine1(1:NN(1,g),1      ,2        ))/2.
        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  1      ,  NN(3,g)) = (fine1(1:NN(1,g),2        ,NN(3,g)) + fine1(1:NN(1,g),1      ,NN(3,g)-1))/2.
        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  1      ) = (fine1(1:NN(1,g),NN(2,g)-1,1      ) + fine1(1:NN(1,g),NN(2,g),2        ))/2.
        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine1(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine1(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) fine2(  1      ,  1      ,1:NN(3,g)) = (fine2(2        ,1      ,1:NN(3,g)) + fine2(1      ,2        ,1:NN(3,g)))/2.
        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) fine2(  1      ,  NN(2,g),1:NN(3,g)) = (fine2(2        ,NN(2,g),1:NN(3,g)) + fine2(1      ,NN(2,g)-1,1:NN(3,g)))/2.
        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) fine2(  NN(1,g),  1      ,1:NN(3,g)) = (fine2(NN(1,g)-1,1      ,1:NN(3,g)) + fine2(NN(1,g),2        ,1:NN(3,g)))/2.
        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) fine2(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine2(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine2(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.
  
        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) fine2(  1      ,1:NN(2,g),  1      ) = (fine2(2        ,1:NN(2,g),1      ) + fine2(1      ,1:NN(2,g),2        ))/2.
        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) fine2(  1      ,1:NN(2,g),  NN(3,g)) = (fine2(2        ,1:NN(2,g),NN(3,g)) + fine2(1      ,1:NN(2,g),NN(3,g)-1))/2.
        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) fine2(  NN(1,g),1:NN(2,g),  1      ) = (fine2(NN(1,g)-1,1:NN(2,g),1      ) + fine2(NN(1,g),1:NN(2,g),2        ))/2.
        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) fine2(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine2(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine2(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) fine2(1:NN(1,g),  1      ,  1      ) = (fine2(1:NN(1,g),2        ,1      ) + fine2(1:NN(1,g),1      ,2        ))/2.
        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) fine2(1:NN(1,g),  1      ,  NN(3,g)) = (fine2(1:NN(1,g),2        ,NN(3,g)) + fine2(1:NN(1,g),1      ,NN(3,g)-1))/2.
        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) fine2(1:NN(1,g),  NN(2,g),  1      ) = (fine2(1:NN(1,g),NN(2,g)-1,1      ) + fine2(1:NN(1,g),NN(2,g),2        ))/2.
        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) fine2(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine2(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine2(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.
  
  
        if (ls1 ==  0 .and. (BC(1,1,g) == 0 .or. BC(1,1,g) == -1)) iimax = iimax-1
        if (ls1 == -1 .and. (BC(2,1,g) == 0 .or. BC(2,1,g) == -1)) iimax = iimax-1
  
        if (ls2 ==  0 .and. (BC(1,2,g) == 0 .or. BC(1,2,g) == -1)) jjmax = jjmax-1
        if (ls2 == -1 .and. (BC(2,2,g) == 0 .or. BC(2,2,g) == -1)) jjmax = jjmax-1
  
        if (ls3 ==  0 .and. (BC(1,3,g) == 0 .or. BC(1,3,g) == -1)) kkmax = kkmax-1
        if (ls3 == -1 .and. (BC(2,3,g) == 0 .or. BC(2,3,g) == -1)) kkmax = kkmax-1
  
  
        if (1 == 2) then ! TEST!!!
            if (add_yes) then
                !pgi$ unroll = n:8
                coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk) - fine2(1:imax:di,1:jmax:dj,1:kmax:dk)
            else
                coarse(1:iimax,1:jjmax,1:kkmax) = fine1(1:imax:di,1:jmax:dj,1:kmax:dk)
            end if
     
        else
     
            if (add_yes) then ! TEST!!! etwas seriöser einbauen ...
     
                call exchange_relax(g,0,0,0,0,.true.,fine1)
                call exchange_relax(g,0,0,0,0,.true.,fine2)
     
                if (dimens == 3) then
                    do kk = 1, kkmax
                        k = dk*(kk-1)+1
                        do jj = 1, jjmax
                            j = dj*(jj-1)+1
                            do ii = 1, iimax
                                i = di*(ii-1)+1
                                coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                                    &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                                    &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                                    &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                                    &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k)) +  &
                                    &              cR3(-1,kk,g+1)*(fine1(i,j,k-1)-fine2(i,j,k-1)) +  &
                                    &              cR3( 1,kk,g+1)*(fine1(i,j,k+1)-fine2(i,j,k+1))) / 3.
                            end do
                        end do
                    end do
                else
                    k  = 1
                    kk = 1
                    do jj = 1, jjmax
                        j = dj*(jj-1)+1
                        do ii = 1, iimax
                            i = di*(ii-1)+1
                            coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*(fine1(i,j,k)-fine2(i,j,k)) +   &
                                &              cR1(-1,ii,g+1)*(fine1(i-1,j,k)-fine2(i-1,j,k)) +  &
                                &              cR1( 1,ii,g+1)*(fine1(i+1,j,k)-fine2(i+1,j,k)) +  &
                                &              cR2(-1,jj,g+1)*(fine1(i,j-1,k)-fine2(i,j-1,k)) +  &
                                &              cR2( 1,jj,g+1)*(fine1(i,j+1,k)-fine2(i,j+1,k))) / 2.
                        end do
                    end do
                end if
     
            else
     
                call exchange_relax(g,0,0,0,0,.true.,fine1)
     
                if (dimens == 3) then
                    do kk = 1, kkmax
                        k = dk*(kk-1)+1
                        do jj = 1, jjmax
                            j = dj*(jj-1)+1
                            do ii = 1, iimax
                                i = di*(ii-1)+1
                                coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1)+cR3(0,kk,g+1))*fine1(i,j,k) +   &
                                    &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                                    &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                                    &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                                    &              cR2( 1,jj,g+1)*fine1(i,j+1,k) +  &
                                    &              cR3(-1,kk,g+1)*fine1(i,j,k-1) +  &
                                    &              cR3( 1,kk,g+1)*fine1(i,j,k+1)) / 3.
                            end do
                        end do
                    end do
                else
                    k  = 1
                    kk = 1
                    do jj = 1, jjmax
                        j = dj*(jj-1)+1
                        do ii = 1, iimax
                            i = di*(ii-1)+1
                            coarse(ii,jj,kk) = ((cR1(0,ii,g+1)+cR2(0,jj,g+1))*fine1(i,j,k) +   &
                                &              cR1(-1,ii,g+1)*fine1(i-1,j,k) +  &
                                &              cR1( 1,ii,g+1)*fine1(i+1,j,k) +  &
                                &              cR2(-1,jj,g+1)*fine1(i,j-1,k) +  &
                                &              cR2( 1,jj,g+1)*fine1(i,j+1,k)) / 2.
                        end do
                    end do
                end if
     
            end if
     
        end if
  
  
  
        if (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) > 1) then
     
            allocate(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
            do kk = 1, kkmax
                do jj = 1, jjmax
                    do ii = 1, iimax
                        sendbuf(ii,jj,kk) = coarse(ii,jj,kk)
                    end do
                end do
            end do
     
            call MPI_GATHERv(sendbuf,iimax*jjmax*kkmax,MPI_REAL8,recvbuf,recvR(1,g+1),dispR(1,g+1),MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
            deallocate(sendbuf)
     
     
            if (participate_yes(g+1)) then
                do k = 1, n_gather(3,g+1)
                    do j = 1, n_gather(2,g+1)
                        do i = 1, n_gather(1,g+1)

                            sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                            offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                            dispg      = dispR(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                 
                            do kk = 1, sizsg(3)
                                do jj = 1, sizsg(2)
                                    do ii = 1, sizsg(1)
                                        coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
                                    end do
                                end do
                            end do
                 
                        end do
                    end do
                end do
            end if
     
        end if
  
  
    end subroutine restrict
  
  
  
  
  
  
  
  
  
  
  
    subroutine restrict_Helmholtz(coarse,fine1,fine2)
  
        implicit none
  
        integer, parameter     ::  g = 1
  
        real   , intent(out  ) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
        real   , intent(inout) ::  fine1 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(inout) ::  fine2 (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
        integer                ::  i, ii, di
        integer                ::  j, jj, dj
        integer                ::  k, kk, dk
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
        !              - Da nur in Richtung der jeweiligen Geschwindigkeitskomponente gemittelt wird, muss nicht   !
        !                die spezialisierte Helmholtz-Variante aufgerufen werden.                                  !
        !              - Austauschrichtung ist invers zu ex1, ex2, ex3. Bei mehreren Blöcken wird auch der jeweils !
        !                redundante "überlappende" Punkt aufgrund der zu grossen Intervallgrenzen (1:iimax) zwar   !
        !                berechnet, aber aufgrund des Einweg-Austauschs falsch berechnet! Dieses Vorgehen wurde    !
        !                bislang aus übersichtsgründen vorgezogen, macht aber eine Initialisierung notwendig.      !
        !                Dabei werden Intervalle der Form 0:imax anstelle von 1:imax bearbeitet, da hier nur die   !
        !                das feinste Geschwindigkeitsgitter behandelt wird!                                        !
        !              - INTENT(inout) ist bei den feinen Gittern notwendig, da Ghost-Werte ausgetauscht werden    !
        !                müssen.                                                                                   !
        !              - Zuviele Daten werden ausgetauscht; eigentlich müsste in der Grenzfläche nur jeder 4.      !
        !                Punkt behandelt werden (4x zuviel!). Leider etwas unschön, könnte aber durch eine         !
        !                spezialisierte Austauschroutine behandelt werden, da das übergeben von Feldern mit        !
        !                Intervallen von b1L:(iimax+b1U) nur sehr schlecht funktionieren würde (d.h. mit Um-       !
        !                kopieren).                                                                                !
        !----------------------------------------------------------------------------------------------------------!
  
        ! TEST!!! Test schreiben, um n_gather(:,2) .GT. 1 hier zu vermeiden! Gleiches gilt natürlich für die Interpolation.
  
  
        di = (NN(1,g)-1)/(NN(1,g+1)-1)
        dj = (NN(2,g)-1)/(NN(2,g+1)-1)
        dk = (NN(3,g)-1)/(NN(3,g+1)-1)
  
        call exchange_relax(g,-ex1,-ex2,-ex3,direction,.false.,fine1)
        call exchange_relax(g,-ex1,-ex2,-ex3,direction,.false.,fine2)
  
  
        if (ls1 ==  0 .and. (BC(1,1,g) == 0 .or. BC(1,1,g) == -1)) fine1(1      ,0:NN(2,g),0:NN(3,g)) = 0.
        if (ls1 ==  0 .and. (BC(1,1,g) == 0 .or. BC(1,1,g) == -1)) fine2(1      ,0:NN(2,g),0:NN(3,g)) = 0.
        if (ls1 == -1 .and. (BC(2,1,g) == 0 .or. BC(2,1,g) == -1)) fine1(NN(1,g),0:NN(2,g),0:NN(3,g)) = 0.
        if (ls1 == -1 .and. (BC(2,1,g) == 0 .or. BC(2,1,g) == -1)) fine2(NN(1,g),0:NN(2,g),0:NN(3,g)) = 0.
  
        if (ls2 ==  0 .and. (BC(1,2,g) == 0 .or. BC(1,2,g) == -1)) fine1(0:NN(1,g),1      ,0:NN(3,g)) = 0.
        if (ls2 ==  0 .and. (BC(1,2,g) == 0 .or. BC(1,2,g) == -1)) fine2(0:NN(1,g),1      ,0:NN(3,g)) = 0.
        if (ls2 == -1 .and. (BC(2,2,g) == 0 .or. BC(2,2,g) == -1)) fine1(0:NN(1,g),NN(2,g),0:NN(3,g)) = 0.
        if (ls2 == -1 .and. (BC(2,2,g) == 0 .or. BC(2,2,g) == -1)) fine2(0:NN(1,g),NN(2,g),0:NN(3,g)) = 0.
  
        if (ls3 ==  0 .and. (BC(1,3,g) == 0 .or. BC(1,3,g) == -1)) fine1(0:NN(1,g),0:NN(2,g),1      ) = 0.
        if (ls3 ==  0 .and. (BC(1,3,g) == 0 .or. BC(1,3,g) == -1)) fine2(0:NN(1,g),0:NN(2,g),1      ) = 0.
        if (ls3 == -1 .and. (BC(2,3,g) == 0 .or. BC(2,3,g) == -1)) fine1(0:NN(1,g),0:NN(2,g),NN(3,g)) = 0.
        if (ls3 == -1 .and. (BC(2,3,g) == 0 .or. BC(2,3,g) == -1)) fine2(0:NN(1,g),0:NN(2,g),NN(3,g)) = 0.
  
  
        !===========================================================================================================
        if (direction == 1) then
     
            do kk = 1, NN(3,g+1)
                k = dk*kk-1
                do jj = 1, NN(2,g+1)
                    j = dj*jj-1
                    do ii = 1, NN(1,g+1)
                        i = di*ii-1
                        coarse(ii,jj,kk) = cRH1(1,i)*(fine1(i-1,j,k)-fine2(i-1,j,k)) + cRH1(2,i)*(fine1(i,j,k)-fine2(i,j,k))
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
        if (direction == 2) then
     
            do kk = 1, NN(3,g+1)
                k = dk*kk-1
                do jj = 1, NN(2,g+1)
                    j = dj*jj-1
                    do ii = 1, NN(1,g+1)
                        i = di*ii-1
                        coarse(ii,jj,kk) = cRH2(1,j)*(fine1(i,j-1,k)-fine2(i,j-1,k)) + cRH2(2,j)*(fine1(i,j,k)-fine2(i,j,k))
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
        if (direction == 3) then
     
            do kk = 1, NN(3,g+1)
                k = dk*kk-1
                do jj = 1, NN(2,g+1)
                    j = dj*jj-1
                    do ii = 1, NN(1,g+1)
                        i = di*ii-1
                        coarse(ii,jj,kk) = cRH3(1,k)*(fine1(i,j,k-1)-fine2(i,j,k-1)) + cRH3(2,k)*(fine1(i,j,k)-fine2(i,j,k))
                    end do
                end do
            end do
     
        end if
    !===========================================================================================================
  
  
    end subroutine restrict_Helmholtz
  
  
  
  
  
  
  
  
  
  
  
    subroutine interpolate(add_yes,g,coarse,fine,work)
  
        implicit none
  
        logical, intent(in   ) ::  add_yes
        integer, intent(in   ) ::  g
        real   , intent(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
        real   , intent(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
        integer                ::  i, di, imax, iimax, iiShift
        integer                ::  j, dj, jmax, jjmax, jjShift
        integer                ::  k, dk, kmax, kkmax, kkShift
  
        !******************************************************************
        integer                ::  ii, jj, kk
        integer                ::  dispg, offsg(1:3)
        real   , allocatable   ::  sendbuf(:,:,:)
        real                   ::  recvbuf(1:(NN(1,g+1)+NB(1,g)-1)*(NN(2,g+1)+NB(2,g)-1)*(NN(3,g+1)+NB(3,g)-1)) ! TEST!!! Ist das richtig so?
        !******************************************************************
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - f�r allgemeine di, dj, dk geeignet                                                        !
        !              - di /= 1 <===> N1 /= 1                                                                     !
        !              - es wird nur in eine Richung ausgetauscht                                                  !
        !              - Null-Setzen am Rand nicht notwendig                                                       !
        !              - Es wird sequentiell �ber alle Raumrichtungen interpoliert, um keinen individuellen        !
        !                Interpolationsstencil f�r jeden Punkt im Raum speichern zu m�ssen.                        !
        !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
        !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
        !                ist dabei das zus�tzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
        !                Routine entliehen).                                               .                       !
        !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
        !                auf die entsprechenden Indizes des gr�beren Gitters umrechnen zu m�ssen.                  !
        !              - Die Block-�berlappenden Stirnfl�chen werden ebenfalls mitverarbeitet, aber eigentlich     !
        !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
        !                notwendig w�re. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
        !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
        !----------------------------------------------------------------------------------------------------------!
  
  
        ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
        !            das Prefetching! Andererseits ist der Geschwindigkeitsgewinn nur sehr gering (schon getestet).
        !
        ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
        ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
        imax    = NN(1,g)
        jmax    = NN(2,g)
        kmax    = NN(3,g)
  
        iimax   = (NN(1,g+1)-1)/n_gather(1,g+1)+1
        jjmax   = (NN(2,g+1)-1)/n_gather(2,g+1)+1
        kkmax   = (NN(3,g+1)-1)/n_gather(3,g+1)+1
  
        di      = (NN(1,g)-1)/(iimax-1)
        dj      = (NN(2,g)-1)/(jjmax-1)
        dk      = (NN(3,g)-1)/(kkmax-1)
  
        iiShift = (iimax-1)*MOD(iB(1,g)-1,n_gather(1,g+1))
        jjShift = (jjmax-1)*MOD(iB(2,g)-1,n_gather(2,g+1))
        kkShift = (kkmax-1)*MOD(iB(3,g)-1,n_gather(3,g+1))
  
  
        if (BC(1,1,g+1) > 0 .and. BC(1,2,g+1) > 0) coarse(1        ,1        ,1:NN(3,g+1)) = (coarse(2          ,1        ,1:NN(3,g+1)) + coarse(1        ,2          ,1:NN(3,g+1)) + coarse(2          ,2          ,1:NN(3,g+1)))/3.
        if (BC(1,1,g+1) > 0 .and. BC(2,2,g+1) > 0) coarse(1        ,NN(2,g+1),1:NN(3,g+1)) = (coarse(2          ,NN(2,g+1),1:NN(3,g+1)) + coarse(1        ,NN(2,g+1)-1,1:NN(3,g+1)) + coarse(2          ,NN(2,g+1)-1,1:NN(3,g+1)))/3.
        if (BC(2,1,g+1) > 0 .and. BC(1,2,g+1) > 0) coarse(NN(1,g+1),1        ,1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,1        ,1:NN(3,g+1)) + coarse(NN(1,g+1),2          ,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,2          ,1:NN(3,g+1)))/3.
        if (BC(2,1,g+1) > 0 .and. BC(2,2,g+1) > 0) coarse(NN(1,g+1),NN(2,g+1),1:NN(3,g+1)) = (coarse(NN(1,g+1)-1,NN(2,g+1),1:NN(3,g+1)) + coarse(NN(1,g+1),NN(2,g+1)-1,1:NN(3,g+1)) + coarse(NN(1,g+1)-1,NN(2,g+1)-1,1:NN(3,g+1)))/3.
  
        if (BC(1,1,g+1) > 0 .and. BC(1,3,g+1) > 0) coarse(1        ,1:NN(2,g+1),1        ) = (coarse(2          ,1:NN(2,g+1),1        ) + coarse(1        ,1:NN(2,g+1),2          ) + coarse(2          ,1:NN(2,g+1),2          ))/3.
        if (BC(1,1,g+1) > 0 .and. BC(2,3,g+1) > 0) coarse(1        ,1:NN(2,g+1),NN(3,g+1)) = (coarse(2          ,1:NN(2,g+1),NN(3,g+1)) + coarse(1        ,1:NN(2,g+1),NN(3,g+1)-1) + coarse(2          ,1:NN(2,g+1),NN(3,g+1)-1))/3.
        if (BC(2,1,g+1) > 0 .and. BC(1,3,g+1) > 0) coarse(NN(1,g+1),1:NN(2,g+1),1        ) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),1        ) + coarse(NN(1,g+1),1:NN(2,g+1),2          ) + coarse(NN(1,g+1)-1,1:NN(2,g+1),2          ))/3.
        if (BC(2,1,g+1) > 0 .and. BC(2,3,g+1) > 0) coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)) = (coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)) + coarse(NN(1,g+1),1:NN(2,g+1),NN(3,g+1)-1) + coarse(NN(1,g+1)-1,1:NN(2,g+1),NN(3,g+1)-1))/3.

        if (BC(1,2,g+1) > 0 .and. BC(1,3,g+1) > 0) coarse(1:NN(1,g+1),1        ,1        ) = (coarse(1:NN(1,g+1),2          ,1        ) + coarse(1:NN(1,g+1),1        ,2          ) + coarse(1:NN(1,g+1),2          ,2          ))/3.
        if (BC(1,2,g+1) > 0 .and. BC(2,3,g+1) > 0) coarse(1:NN(1,g+1),1        ,NN(3,g+1)) = (coarse(1:NN(1,g+1),2          ,NN(3,g+1)) + coarse(1:NN(1,g+1),1        ,NN(3,g+1)-1) + coarse(1:NN(1,g+1),2          ,NN(3,g+1)-1))/3.
        if (BC(2,2,g+1) > 0 .and. BC(1,3,g+1) > 0) coarse(1:NN(1,g+1),NN(2,g+1),1        ) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,1        ) + coarse(1:NN(1,g+1),NN(2,g+1),2          ) + coarse(1:NN(1,g+1),NN(2,g+1)-1,2          ))/3.
        if (BC(2,2,g+1) > 0 .and. BC(2,3,g+1) > 0) coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)) = (coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)) + coarse(1:NN(1,g+1),NN(2,g+1),NN(3,g+1)-1) + coarse(1:NN(1,g+1),NN(2,g+1)-1,NN(3,g+1)-1))/3.
  
  
        call exchange_relax(g+1,ex1,ex2,ex3,0,.false.,coarse) ! Anmerkung: .FALSE. ist ok ...
  
  
        !***********************************************************************************************************
        if (n_gather(1,g+1)*n_gather(2,g+1)*n_gather(3,g+1) > 1) then
     
            if (participate_yes(g+1)) then
                do k = 1, n_gather(3,g+1)
                    do j = 1, n_gather(2,g+1)
                        do i = 1, n_gather(1,g+1)
                 
                            dispg      = dispI(    i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)
                            offsg(1:3) = offsI(1:3,i+(j-1)*n_gather(1,g+1)+(k-1)*n_gather(1,g+1)*n_gather(2,g+1),g+1)

                            do kk = 1, kkmax
                                do jj = 1, jjmax
                                    do ii = 1, iimax
                                        recvbuf(dispg+ii+(jj-1)*iimax+(kk-1)*iimax*jjmax) = coarse(ii+offsg(1),jj+offsg(2),kk+offsg(3))
                                    end do
                                end do
                            end do
                 
                        end do
                    end do
                end do
            end if
     
     
            allocate(sendbuf(1:iimax,1:jjmax,1:kkmax)) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "coarse" verwenden!
     
            call MPI_SCATTER(recvbuf,iimax*jjmax*kkmax,MPI_REAL8,sendbuf,iimax*jjmax*kkmax,MPI_REAL8,rankc2(g+1),comm2(g+1),merror)
     
            do kk = 1, kkmax
                k = dk*(kk-1)+1
                do jj = 1, jjmax
                    j = dj*(jj-1)+1
                    !pgi$ unroll = n:8
                    do ii = 1, iimax
                        i = di*(ii-1)+1
                        work(i,j,k) = sendbuf(ii,jj,kk)
                    end do
                end do
            end do
     
            deallocate(sendbuf)
     
        else
     
            !pgi$ unroll = n:8
            work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse((1+iiShift):(iimax+iiShift),(1+jjShift):(jjmax+jjShift),(1+kkShift):(kkmax+kkShift))
     
        end if
        !***********************************************************************************************************
  
  
  
        !===========================================================================================================
        if (dk /= 1) then ! (dimens == 2) <==> (dk == 1) automatisch erf�llt!
     
            do k = 2, kmax-1, dk
                do j = 1, jmax, dj
                    !pgi$ unroll = n:8
                    do i = 1, imax, di
                        work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
        if (dj /= 1) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (kmax == 2??)
     
            do k = 1, kmax
                do j = 2, jmax-1, dj
                    !pgi$ unroll = n:8
                    do i = 1, imax, di
                        work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
        if (di /= 1) then
     
            do k = 1, kmax
                do j = 1, jmax
                    !pgi$ unroll = n:8
                    do i = 2, imax-1, di
                        work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
  
  
        if (add_yes) then
            !pgi$ unroll = n:8
            fine(1:imax,1:jmax,1:kmax) = fine(1:imax,1:jmax,1:kmax) + work(1:imax,1:jmax,1:kmax)
        end if
  
  
    end subroutine interpolate
  
  
  
  
  
  
  
  
  
  
  
    subroutine interpolate_Helmholtz(coarse,fine,work)
  
        implicit none
  
        integer, parameter     ::  g = 1
  
        real   , intent(inout) ::  coarse(b1L:(NN(1,g+1)+b1U),b2L:(NN(2,g+1)+b2U),b3L:(NN(3,g+1)+b3U))
        real   , intent(inout) ::  fine  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
        real   , intent(out  ) ::  work  (b1L:(NN(1,g  )+b1U),b2L:(NN(2,g  )+b2U),b3L:(NN(3,g  )+b3U))
  
        integer                ::  i, di, imax, iimax
        integer                ::  j, dj, jmax, jjmax
        integer                ::  k, dk, kmax, kkmax
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
        !              - di /= 1 <===> N1 /= 1                                                                     !
        !              - Da imax=N1 usw., k�nnten auf dem feinen Gitter auch problemlos die bekannten engeren      !
        !                Intervallgrenzen S11:N11 usw. benutzt werden. Wurde bislang aus �bersichtsgr�nden nicht   !
        !                vollzogen.                                                                                !
        !              - Die Block-�berlappenden Stirnfl�chen werden ebenfalls mitverarbeitet, aber eigentlich     !
        !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
        !                notwendig w�re. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
        !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
        !              - Es wird sequentiell �ber alle Raumrichtungen interpoliert, um keinen individuellen        !
        !                Interpolationsstencil f�r jeden Punkt im Raum speichern zu m�ssen.                        !
        !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
        !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
        !                ist dabei das zus�tzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
        !                Routine entliehen).                                               .                       !
        !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
        !                auf die entsprechenden Indizes des gr�beren Gitters umrechnen zu m�ssen.                  !
        !----------------------------------------------------------------------------------------------------------!
  
  
        iimax = NN(1,g+1) ! TEST!!! iimax, etc. substituieren ...
        jjmax = NN(2,g+1)
        kkmax = NN(3,g+1)
  
        imax  = NN(1,g)
        jmax  = NN(2,g)
        kmax  = NN(3,g)
  
        di = (imax-1)/(iimax-1)
        dj = (jmax-1)/(jjmax-1)
        dk = (kmax-1)/(kkmax-1)
  
  
        ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
        !            das Prefetching!
  
        ! Anmerkung: .FALSE. ist ok!
        call exchange_relax(g+1,ex1,ex2,ex3,0,.false.,coarse)
  
  
        !pgi$ unroll = n:8
        work(1:imax:di,1:jmax:dj,1:kmax:dk) = coarse(1:iimax,1:jjmax,1:kkmax)
  
  
        !===========================================================================================================
        if (direction == 1) then
     
            if (dk /= 1) then ! (dimens == 2) <==> (dk == 1) automatisch erf�llt!
                do k = 2, kmax-1, dk
                    do j = 1, jmax, dj
                        !pgi$ unroll = n:8
                        do i = 1, imax, di
                            work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
                        end do
                    end do
                end do
            end if
     
            !--------------------------------------------------------------------------------------------------------
            if (dj /= 1) then
                do k = 1, kmax
                    do j = 2, jmax-1, dj
                        !pgi$ unroll = n:8
                        do i = 1, imax, di
                            work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
                        end do
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
     
            if (di /= 1) then
                do k = 1, kmax
                    do j = 1, jmax
                        !pgi$ unroll = n:8
                        do i = 1, imax-2, di
                            fine(i  ,j,k) = fine(i  ,j,k) + cIH1(1,i  )*work(i,j,k) + cIH1(2,i  )*work(i+2,j,k)
                            fine(i+1,j,k) = fine(i+1,j,k) + cIH1(1,i+1)*work(i,j,k) + cIH1(2,i+1)*work(i+2,j,k)
                        end do
                    end do
                end do
        
                !--- Extrapolation am Rand ---
                if (BC_1L > 0) then
                    do k = 1, kmax
                        !pgi$ unroll = n:8
                        do j = 1, jmax
                            fine(0   ,j,k) = fine(0   ,j,k) + cIH1(1,0   )*work(1      ,j,k) + cIH1(2,0   )*work(1+di,j,k)
                        end do
                    end do
                end if
                if (BC_1U > 0) then
                    do k = 1, kmax
                        !pgi$ unroll = n:8
                        do j = 1, jmax
                            fine(imax,j,k) = fine(imax,j,k) + cIH1(1,imax)*work(imax-di,j,k) + cIH1(2,imax)*work(imax,j,k)
                        end do
                    end do
                end if
            else
                if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
                call MPI_FINALIZE(merror)
                stop
            end if

        end if
        !===========================================================================================================
        if (direction == 2) then
     
            if (dk /= 1) then ! (dimens == 2) <==> (dk == 1) automatisch erf�llt!
                do k = 2, kmax-1, dk
                    do j = 1, jmax, dj
                        !pgi$ unroll = n:8
                        do i = 1, imax, di
                            work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
                        end do
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (di /= 1) then
                do k = 1, kmax
                    do j = 1, jmax, dj
                        !pgi$ unroll = n:8
                        do i = 2, imax-1, di
                            work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
                        end do
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (dj /= 1) then
                do k = 1, kmax
                    do i = 1, imax
                        !pgi$ unroll = n:8
                        do j = 1, jmax-2, dj
                            fine(i,j  ,k) = fine(i,j  ,k) + cIH2(1,j  )*work(i,j,k) + cIH2(2,j  )*work(i,j+2,k)
                            fine(i,j+1,k) = fine(i,j+1,k) + cIH2(1,j+1)*work(i,j,k) + cIH2(2,j+1)*work(i,j+2,k)
                        end do
                    end do
                end do
        
                !--- Extrapolation am Rand ---
                if (BC_2L > 0) then
                    do k = 1, kmax
                        !pgi$ unroll = n:8
                        do i = 1, imax
                            fine(i,0   ,k) = fine(i,0   ,k) + cIH2(1,0   )*work(i,1      ,k) + cIH2(2,0   )*work(i,1+dj,k)
                        end do
                    end do
                end if
                if (BC_2U > 0) then
                    do k = 1, kmax
                        !pgi$ unroll = n:8
                        do i = 1, imax
                            fine(i,jmax,k) = fine(i,jmax,k) + cIH2(1,jmax)*work(i,jmax-dj,k) + cIH2(2,jmax)*work(i,jmax,k)
                        end do
                    end do
                end if
            else
                if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
                call MPI_FINALIZE(merror)
                stop
            end if

        end if
        !===========================================================================================================
        if (direction == 3) then ! (dimens == 2) <==> (direction /= 3) automatisch erf�llt!
     
            if (dj /= 1) then
                do k = 1, kmax, dk
                    do j = 2, jmax-1, dj
                        !pgi$ unroll = n:8
                        do i = 1, imax, di
                            work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
                        end do
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (di /= 1) then
                do k = 1, kmax, dk
                    do j = 1, jmax
                        !pgi$ unroll = n:8
                        do i = 2, imax-1, di
                            work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
                        end do
                    end do
                end do
            end if
            !--------------------------------------------------------------------------------------------------------
            if (dk /= 1) then
                do j = 1, jmax
                    do i = 1, imax
                        !pgi$ unroll = n:8
                        do k = 1, kmax-2, dk
                            fine(i,j,k  ) = fine(i,j,k  ) + cIH3(1,k  )*work(i,j,k) + cIH3(2,k  )*work(i,j,k+2)
                            fine(i,j,k+1) = fine(i,j,k+1) + cIH3(1,k+1)*work(i,j,k) + cIH3(2,k+1)*work(i,j,k+2)
                        end do
                    end do
                end do

                !--- Extrapolation am Rand ---
                if (BC_3L > 0) then
                    do j = 1, jmax
                        !pgi$ unroll = n:8
                        do i = 1, imax
                            fine(i,j,0   ) = fine(i,j,0   ) + cIH3(1,0   )*work(i,j,1      ) + cIH3(2,0   )*work(i,j,1+dk)
                        end do
                    end do
                end if
                if (BC_3U > 0) then
                    do j = 1, jmax
                        !pgi$ unroll = n:8
                        do i = 1, imax
                            fine(i,j,kmax) = fine(i,j,kmax) + cIH3(1,kmax)*work(i,j,kmax-dk) + cIH3(2,kmax)*work(i,j,kmax)
                        end do
                    end do
                end if
            else
                if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
                call MPI_FINALIZE(merror)
                stop
            end if
     
        end if
    !===========================================================================================================


    end subroutine interpolate_Helmholtz

  
  
  
  
  
  
  
  
  
    !> \brief BICGstap
    !! \param[in] eps
    !! \param[in] n_it_max
    !! \param[in] init_yes
    !! \param[in] SS1 start index in 1-direction
    !! \param[in] SS2 start index in 2-direction
    !! \param[in] SS3 start index in 3-direction
    !! \param[in] NN1 end index in 1-direction
    !! \param[in] NN2 end index in 2-direction
    !! \param[in] NN3 end index in 3-direction
    !! \param[in] preconditioner
    !! \param[in] problem_type problem types:
    !!                                       - 1: Helmholtz
    !!                                       - 2: div grad
    !!                                       - 3: Helmholtz_conc
    !!                                       - 4: div grad transposted
    !! \param[in] bb
    !! \param[inout] phi
    !! \param[in] quiet_yes1
    !! \param[in] quiet_yes2
    subroutine BiCGstab(    &
        eps,                &
        n_it_max,           &
        init_yes,           &
        SS1,SS2,SS3,        &
        NN1,NN2,NN3,        &
        bb,                 &
        phi,                &
        problem_type,       &
        quiet_yes1,         &
        quiet_yes2,         &
        preconditioner)
  
        implicit none
  
        real   , intent(in)    ::  eps
        integer, intent(in)    ::  n_it_max
        logical, intent(in)    ::  init_yes
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        integer, intent(in)    ::  preconditioner
        integer, intent(in)    ::  problem_type
  
        real   , intent(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
        integer                ::  counter
        integer                ::  i, j, k
  
        real                   ::  norm2   !, norm2_global
        real                   ::  norm_inf, norm_inf_global, norm_inf_prev
        real                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
        real                   ::  scalar_global(1:2)
        real                   ::  alpha, beta, omega
  
        logical                ::  exit_yes
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
        !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
        !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
        !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
        !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
        !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
        !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
        !                genommen wird.                                                                            !
        !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
        !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
        !              - rhr muss nicht initialisiert werden!                                                      !
        !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
        !                hinter "CALL status_iteration" eingeführt.                                                !
        !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
        !----------------------------------------------------------------------------------------------------------!
  
  
        ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
        ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
        !===========================================================================================================
        if (.not. init_yes) then
            if (problem_type == 1) call product_Helmholtz(phi,rr)
            if (problem_type == 2) call product_div_grad (phi,rr)
            if (problem_type == 3) call product_Helmholtz_conc(phi,rr)
            if (problem_type == 4) call product_div_grad_transp(phi,rr)
     
            rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
            rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        end if
        !===========================================================================================================
  
        !===========================================================================================================
        !=== Residuum ==============================================================================================
        !===========================================================================================================
        if (init_yes) then
            call get_norms(SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.true.,.true. ,normInf=norm_inf,normTwo=norm2)
        else
            call get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.true.,.false.,normInf=norm_inf)
        end if
        norm_inf_prev = norm_inf
        !===========================================================================================================
  
        counter = 0
  
  
        ITERATE: do

            !========================================================================================================
            !=== überprüfen des Konvergenzkriteriums ================================================================
            !========================================================================================================
            call status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)

            ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
            !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
            !    zu rechnen.
            if (problem_type == 1 .and.                           counter < 1 .and. norm_inf /= 0.) exit_yes = .false.
            if (problem_type == 2 .and. number_poisson == 1 .and. counter < 1 .and. norm_inf /= 0.) exit_yes = .false.
            if (problem_type == 2 .and. number_poisson == 2 .and. counter < 0 .and. norm_inf /= 0.) exit_yes = .false.

            if (exit_yes .and. counter == 0 .and. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
            if (exit_yes) exit ITERATE
            !========================================================================================================


            !========================================================================================================
            !=== nächster Durchlauf =================================================================================
            !========================================================================================================
            counter = counter + 1
            !========================================================================================================


            !========================================================================================================
            rhr_prev = rhr
            if (init_yes) then
                if (counter == 1) then
                    rhr = norm2
                else
                    call product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
                end if
            else
                call product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
            end if
            !========================================================================================================
            if (ABS(rhr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhr =', rhr
                if (rank == 0) write(10,'(a,E13.5)') 'rhr =', rhr

                rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
                call product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
            end if
            !========================================================================================================
            if (counter >= 2) then
                if (omega == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'omega =', omega
                    if (rank == 0) write(10,'(a,E13.5)') 'omega =', omega
                    if (problem_type == 2 .or. problem_type == 4) then
                        exit ITERATE
                    else
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                if (rhr_prev == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    if (rank == 0) write(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    if (problem_type == 2 .or. problem_type == 4) then
                        exit ITERATE
                    else
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                beta = (alpha/omega)*(rhr/rhr_prev)
                omega = -beta*omega
                pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
            else
                if (init_yes .and. counter == 1) then
                    pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
                else
                    pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
                end if
            end if
            !========================================================================================================
            if (preconditioner == 0) then
                if (problem_type == 1) call product_Helmholtz(pp,Ap)
                if (problem_type == 2) call product_div_grad (pp,Ap)
                if (problem_type == 3) call product_Helmholtz_conc(pp,Ap)
                if (problem_type == 4) call product_div_grad_transp(pp,Ap)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                if (preconditioner == 1) call multigridV(.true.,1,pp,z1,problem_type)
                if (preconditioner == 2) call multigridF(.true.,1,pp,z1,problem_type)
                if (problem_type == 1) call product_Helmholtz(z1,Ap)
                if (problem_type == 2) call product_div_grad (z1,Ap)
                if (problem_type == 3) call product_Helmholtz_conc(z1,Ap)
                if (problem_type == 4) call product_div_grad_transp(z1,Ap)
            else
                if (rank == 0) write(*,'(a)') 'ERROR! Specify valid preconditioner!'
                call MPI_FINALIZE(merror)
                stop
            end if
            !========================================================================================================
            if (init_yes) then
                call product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
            else
                call product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
            end if
            !========================================================================================================
            if (ABS(rhAp) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhAp =', rhAp
                if (rank == 0) write(10,'(a,E13.5)') 'rhAp =', rhAp
                if (problem_type == 2 .or. problem_type == 4) then
                    exit ITERATE
                else
                    if (rhr /= 0.) then
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                alpha = 0.
            else
                alpha = rhr / rhAp
            end if
            !========================================================================================================
            if (init_yes .and. counter == 1) then
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
                        end do
                    end do
                end do
            else
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
                        end do
                    end do
                end do
            end if
            !========================================================================================================
            if (preconditioner == 0) then
                if (problem_type == 1) call product_Helmholtz(rr,Ar)
                if (problem_type == 2) call product_div_grad (rr,Ar)
                if (problem_type == 3) call product_Helmholtz_conc(rr,Ar)
                if (problem_type == 4) call product_div_grad_transp(rr,Ar)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                if (preconditioner == 1) call multigridV(.true.,1,rr,z2,problem_type)
                if (preconditioner == 2) call multigridF(.true.,1,rr,z2,problem_type)
                if (problem_type == 1) call product_Helmholtz(z2,Ar)
                if (problem_type == 2) call product_div_grad (z2,Ar)
                if (problem_type == 3) call product_Helmholtz_conc(z2,Ar)
                if (problem_type == 4) call product_div_grad_transp(z2,Ar)
            end if
            !========================================================================================================
            rAr  = 0.
            ArAr = 0.
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
                        ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            rAr  = scalar_global(1)
            ArAr = scalar_global(2)
     
            ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
            !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
            !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
            !========================================================================================================
            if (ABS(rAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rAr =', rAr
                if (rank == 0) write(10,'(a,E13.5)') 'rAr =', rAr
            end if
            if (ABS(ArAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'ArAr =', ArAr
                if (rank == 0) write(10,'(a,E13.5)') 'ArAr =', ArAr
                if (problem_type == 2 .or. problem_type == 4) then
                    exit ITERATE
                else
                    if (rAr /= 0.) then
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                omega = 0.
            else
                omega = rAr / ArAr
            end if
            !========================================================================================================
            norm_inf = 0.
            !--------------------------------------------------------------------------------------------------------
            if (counter == 1 .and. init_yes) then
                if (preconditioner == 0) then
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                !-----------------------------------------------------------------------------------------------------
                else
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                end if
             !--------------------------------------------------------------------------------------------------------
            else
                if (preconditioner == 0) then
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                !-----------------------------------------------------------------------------------------------------
                else
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                end if
            end if
            !--------------------------------------------------------------------------------------------------------
            call MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            norm_inf = norm_inf_global
            !========================================================================================================
     
     
            !========================================================================================================
            !=== Konvergenzstatistik ================================================================================
            !========================================================================================================
            if (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf_prev)
            if (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
            if (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
            if (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
            norm_inf_prev = norm_inf
           !========================================================================================================
     
        end do ITERATE
  
  
    end subroutine BiCGstab
  
  
  
  
  
  
  
  
  
    ! TEST!!! gleiches Spiel; g und SNB uebergeben ....
    ! ACHTUNG!!! Routine sollte nur bei g .GE. 2 aufgerufen werden, weil sonst viele grosse Felder allociert werden ...
    ! TEST!!! weight in Residuum?
    ! N1, N2, N3, g werden übergeben, um z.B. psi auch auf gröberen Gittern berechnen zu können ...
    subroutine BiCGstab2(eps,n_it_max,init_yes,N1,N2,N3,g,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
        implicit none
  
        real   , intent(in)    ::  eps
        integer, intent(in)    ::  n_it_max
        logical, intent(in)    ::  init_yes
  
        integer, intent(in)    ::  N1
        integer, intent(in)    ::  N2
        integer, intent(in)    ::  N3
  
        integer, intent(in)    ::  g
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        integer, intent(in)    ::  preconditioner
        integer, intent(in)    ::  problem_type
  
        real   , intent(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        real                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) ! TEST!!! ! FELD!!!
        real                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
        integer                ::  counter
        integer                ::  i, j, k
  
        real                   ::  norm2   !, norm2_global
        real                   ::  norm_inf, norm_inf_global, norm_inf_prev
        real                   ::  rhr, rhr_prev, ArAr, rAr, rhAp
        real                   ::  scalar_global(1:2)
        real                   ::  alpha, beta, omega
  
        logical                ::  exit_yes
  
  
        !*****************************************************************************************
        !SNB(:,:,g) = (/SS1,NN1/), etc. ! TEST!!!
        !*****************************************************************************************
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Gitterstreckung ist ein Problem bei der Gewichtung der Residuen, was beim Helmholtz-      !
        !                Problem aufgrund der Diagonaldominanz weit weniger ausgeprägt sein sollte als beim reinen !
        !                Poisson-Problem. Da der Multigrid-Glätter (Jacobi, Gauss-Seidel) davon generell unbeein-  !
        !                flusst sind, kann dieses Problem hier mehr oder weniger vernachlässigt werden. Alternativ !
        !                könnten die Stencils analog zu "weight" auch für das gesamte Feld gespeichert werden, was !
        !                allerdings einen sehr grossen Speicheraufwand und lange Speicherzugriffszeiten erwarten   !
        !                lassen, so dass hier eine evtl. geringere Konvergenzrate mit gutem Gewissen in Kauf       !
        !                genommen wird.                                                                            !
        !              - Es wird hier generell rh = rr angenommen, da somit immer <rh,rr> /= 0 gilt!               !
        !              - dAXPY ist nicht wirklich schneller als die direkt programmierten Operationen.             !
        !              - rhr muss nicht initialisiert werden!                                                      !
        !              - Um pp nur bei Bedarf berechnen zu müssen, wird dazu eine Fallunterscheidung nach erst     !
        !                hinter "CALL status_iteration" eingeführt.                                                !
        !              - z2 könnte zu Ungunsten der Geschwindigkeit eingespaart werden.                            !
        !----------------------------------------------------------------------------------------------------------!
  
  
        ! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
        ! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
  
  
        !===========================================================================================================
        if (.not. init_yes) then
            !IF (problem_type == 1) CALL product_Helmholtz(phi,rr)
            !IF (problem_type == 2) CALL product_div_grad (phi,rr)
            !IF (problem_type == 3) CALL product_Helmholtz_conc(phi,rr)
            !IF (problem_type == 4) CALL product_div_grad_transp(phi,rr)
            if (problem_type == 5) call product_div_grad_relax(g,phi,rr)
     
            rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
            rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
        end if
        !===========================================================================================================
  
        !===========================================================================================================
        !=== Residuum ==============================================================================================
        !===========================================================================================================
        if (init_yes) then
            call get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,bb,problem_type,.true.,.true. ,normInf=norm_inf,normTwo=norm2)
        else
            call get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.true.,.false.,normInf=norm_inf)
        end if
        norm_inf_prev = norm_inf
        !===========================================================================================================
  
        counter = 0
  
  
        ITERATE: do
     
            !========================================================================================================
            !=== überprüfen des Konvergenzkriteriums ================================================================
            !========================================================================================================
            call status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
            ! soll Konvergenz sicherstellen (oder offsetprec runtersetzen ...).
            !    norm_inf == 0. << eps ist ein Sonderfall, bei dem es keinen Sinn macht, eine weitere Iteration
            !    zu rechnen.
            if (problem_type == 1 .and.                           counter < 1 .and. norm_inf /= 0.) exit_yes = .false.
            if (problem_type == 2 .and. number_poisson == 1 .and. counter < 1 .and. norm_inf /= 0.) exit_yes = .false.
            if (problem_type == 2 .and. number_poisson == 2 .and. counter < 0 .and. norm_inf /= 0.) exit_yes = .false.
     
            if (exit_yes .and. counter == 0 .and. init_yes) phi(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
            if (exit_yes) exit ITERATE
            !========================================================================================================
     
     
            !========================================================================================================
            !=== nächster Durchlauf =================================================================================
            !========================================================================================================
            counter = counter + 1
            !========================================================================================================
     
     
            !========================================================================================================
            rhr_prev = rhr
            if (init_yes) then
                if (counter == 1) then
                    rhr = norm2
                else
                    call product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,bb,rhr)
                end if
            else
                call product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
            end if
            !========================================================================================================
            if (ABS(rhr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhr =', rhr
                if (rank == 0) write(10,'(a,E13.5)') 'rhr =', rhr
        
                rh(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) ! Neuer Referenzvektor ...
                call product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,rr,rh,rhr)
            end if
            !========================================================================================================
            if (counter >= 2) then
                if (omega == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'omega =', omega
                    if (rank == 0) write(10,'(a,E13.5)') 'omega =', omega
                    if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
                        exit ITERATE
                    else
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                if (rhr_prev == 0.) then
                    if (rank == 0) write(* ,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    if (rank == 0) write(10,'(a,E13.5)') 'rhr_prev =', rhr_prev
                    if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
                        exit ITERATE
                    else
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                beta = (alpha/omega)*(rhr/rhr_prev)
                omega = -beta*omega
                pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3) + beta*pp(SS1:NN1,SS2:NN2,SS3:NN3) + omega*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
            else
                if (init_yes .and. counter == 1) then
                    pp(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
                else
                    pp(SS1:NN1,SS2:NN2,SS3:NN3) = rr(SS1:NN1,SS2:NN2,SS3:NN3)
                end if
            end if
            !========================================================================================================
            if (preconditioner == 0) then
                !IF (problem_type == 1) CALL product_Helmholtz(pp,Ap)
                !IF (problem_type == 2) CALL product_div_grad (pp,Ap)
                !IF (problem_type == 3) CALL product_Helmholtz_conc(pp,Ap)
                !IF (problem_type == 4) CALL product_div_grad_transp(pp,Ap)
                if (problem_type == 5) call product_div_grad_relax(g,pp,Ap)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                if (preconditioner == 1) call multigridV(.true.,g,pp,z1,problem_type)
                if (preconditioner == 2) call multigridF(.true.,g,pp,z1,problem_type)
                !IF (problem_type == 1) CALL product_Helmholtz(z1,Ap)
                !IF (problem_type == 2) CALL product_div_grad (z1,Ap)
                !IF (problem_type == 3) CALL product_Helmholtz_conc(z1,Ap)
                !IF (problem_type == 4) CALL product_div_grad_transp(z1,Ap)
                if (problem_type == 5) call product_div_grad_relax(g,z1,Ap)
            else
                if (rank == 0) write(*,'(a)') 'ERROR! Specify valid preconditioner!'
                call MPI_FINALIZE(merror)
                stop
            end if
            !========================================================================================================
            if (init_yes) then
                call product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,Ap,bb,rhAp)
            else
                call product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,Ap,rh,rhAp)
            end if
            !========================================================================================================
            if (ABS(rhAp) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rhAp =', rhAp
                if (rank == 0) write(10,'(a,E13.5)') 'rhAp =', rhAp
                if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
                    exit ITERATE
                else
                    if (rhr /= 0.) then
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                alpha = 0.
            else
                alpha = rhr / rhAp
            end if
            !========================================================================================================
            if (init_yes .and. counter == 1) then
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            rr(i,j,k) = bb(i,j,k) - alpha*Ap(i,j,k)
                        end do
                    end do
                end do
            else
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            rr(i,j,k) = rr(i,j,k) - alpha*Ap(i,j,k)
                        end do
                    end do
                end do
            end if
            !========================================================================================================
            if (preconditioner == 0) then
                !IF (problem_type == 1) CALL product_Helmholtz(rr,Ar)
                !IF (problem_type == 2) CALL product_div_grad (rr,Ar)
                !IF (problem_type == 3) CALL product_Helmholtz_conc(rr,Ar)
                !IF (problem_type == 4) CALL product_div_grad_transp(rr,Ar)
                if (problem_type == 5) call product_div_grad_relax(g,rr,Ar)
            else if (preconditioner == 1 .or. preconditioner == 2) then
                if (preconditioner == 1) call multigridV(.true.,g,rr,z2,problem_type)
                if (preconditioner == 2) call multigridF(.true.,g,rr,z2,problem_type)
                !IF (problem_type == 1) CALL product_Helmholtz(z2,Ar)
                !IF (problem_type == 2) CALL product_div_grad (z2,Ar)
                !IF (problem_type == 3) CALL product_Helmholtz_conc(z2,Ar)
                !IF (problem_type == 4) CALL product_div_grad_transp(z2,Ar)
                if (problem_type == 5) call product_div_grad_relax(g,z2,Ar)
            end if
            !========================================================================================================
            rAr  = 0.
            ArAr = 0.
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        rAr  = rAr  + rr(i,j,k)*Ar(i,j,k)
                        ArAr = ArAr + Ar(i,j,k)*Ar(i,j,k)
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE((/rAr,ArAr/),scalar_global,2,MPI_REAL8,MPI_SUM,comm1(g),merror)
            rAr  = scalar_global(1)
            ArAr = scalar_global(2)
     
            ! ACHTUNG!!! Zu teuer im Vergleich zur obigen Variante:
            !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,rr,rAr)
            !CALL product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,Ar,Ar,ArAr)
            !========================================================================================================
            if (ABS(rAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'rAr =', rAr
                if (rank == 0) write(10,'(a,E13.5)') 'rAr =', rAr
            end if
            if (ABS(ArAr) == 0.) then
                if (rank == 0) write(* ,'(a,E13.5)') 'ArAr =', ArAr
                if (rank == 0) write(10,'(a,E13.5)') 'ArAr =', ArAr
                if (problem_type == 2 .or. problem_type == 4 .or. problem_type == 5) then
                    exit ITERATE
                else
                    if (rAr /= 0.) then
                        call MPI_FINALIZE(merror)
                        stop
                    end if
                end if
                omega = 0.
            else
                omega = rAr / ArAr
            end if
            !========================================================================================================
            norm_inf = 0.
            !--------------------------------------------------------------------------------------------------------
            if (counter == 1 .and. init_yes) then
                if (preconditioner == 0) then
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*pp(i,j,k) + omega*rr(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                !-----------------------------------------------------------------------------------------------------
                else
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = alpha*z1(i,j,k) + omega*z2(i,j,k)
                                    rr (i,j,k) =       rr(i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                end if
            !--------------------------------------------------------------------------------------------------------
            else
                if (preconditioner == 0) then
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*rr(i,j,k) + alpha*pp(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                !-----------------------------------------------------------------------------------------------------
                else
                    if (problem_type == 2 .and. weighting_yes) then
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)*weight(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    else
                        do k = SS3, NN3
                            do j = SS2, NN2
                                !pgi$ unroll = n:8
                                do i = SS1, NN1
                                    phi(i,j,k) = phi(i,j,k) + omega*z2(i,j,k) + alpha*z1(i,j,k)
                                    rr (i,j,k) = rr (i,j,k) - omega*Ar(i,j,k)
                                    norm_inf = MAX(ABS(rr(i,j,k)),norm_inf)
                                end do
                            end do
                        end do
                    end if
                end if
            end if
            !--------------------------------------------------------------------------------------------------------
            call MPI_ALLREDUCE(norm_inf,norm_inf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            norm_inf = norm_inf_global
            !========================================================================================================
     
     
            !========================================================================================================
            !=== Konvergenzstatistik ================================================================================
            !========================================================================================================
            if (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf_prev)
            if (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf_prev)
     
            if (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
            if (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
            norm_inf_prev = norm_inf
           !========================================================================================================
     
        end do ITERATE
  
  
    end subroutine BiCGstab2
  
  
  
  
  
  
  
  
  
  
  
    subroutine Richardson(eps,n_it_max,init_yes,SS1,SS2,SS3,NN1,NN2,NN3,bb,phi,problem_type,quiet_yes1,quiet_yes2,preconditioner)
  
        implicit none
  
        real   , intent(in)    ::  eps
        integer, intent(in)    ::  n_it_max
        logical, intent(in)    ::  init_yes
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        integer, intent(in)    ::  preconditioner
        integer, intent(in)    ::  problem_type
  
        real   , intent(in)    ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
        integer                ::  counter
        real                   ::  norm_inf, norm_inf0
        logical                ::  exit_yes
  
  
  
        !------------------------------------------------------------
        if (init_yes) then
            !pgi$ unroll = n:8
            rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3)
            phi = 0.
        else
            if (problem_type == 1) call product_Helmholtz(phi,rr)
            if (problem_type == 2) call product_div_grad (phi,rr)
            if (problem_type == 3) call product_Helmholtz_conc(phi,rr)
            if (problem_type == 4) call product_div_grad_transp(phi,rr)
            !pgi$ unroll = n:8
            rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - rr(SS1:NN1,SS2:NN2,SS3:NN3)
        end if
  
        !--- grösstes Residuum --------------------------------------
        call get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.true.,.false.,normInf=norm_inf)
        norm_inf0 = norm_inf
  
        !------------------------------------------------------------
        counter = 0
        !------------------------------------------------------------
        Ap(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
        pp(SS1:NN1,SS2:NN2,SS3:NN3) = 0.
        !------------------------------------------------------------
  
        ITERATE: do
     
            !--- überprüfen des Konvergenzkriteriums -----------------
            call status_iteration(eps,norm_inf,counter,n_it_max,exit_yes,quiet_yes1,quiet_yes2)
     
            if (exit_yes) exit ITERATE
     
     
            !--- nåchster Restart ------------------------------------
            counter = counter + 1
     
     
            if (preconditioner == 0) then
                pp = 0.25*rr ! TEST!!! nicht sehr schön ...
               !pp = 1.0*rr
            else if (preconditioner == 1 .or. preconditioner == 2) then
                if (preconditioner == 1) call multigridV(.true.,1,rr,pp,problem_type)
                if (preconditioner == 2) call multigridF(.true.,1,rr,pp,problem_type)
            else
                if (rank == 0) write(*,'(a)') 'ERROR! Specify valid preconditioner!'
                call MPI_FINALIZE(merror)
                stop
            end if
     
     
            if (1 == 1) then
        
                if (problem_type == 1) call product_Helmholtz      (pp,Ap)
                if (problem_type == 2) call product_div_grad       (pp,Ap)
                if (problem_type == 3) call product_Helmholtz_conc (pp,Ap)
                if (problem_type == 4) call product_div_grad_transp(pp,Ap)
        
                rr (SS1:NN1,SS2:NN2,SS3:NN3) = rr (SS1:NN1,SS2:NN2,SS3:NN3) - Ap(SS1:NN1,SS2:NN2,SS3:NN3)
                phi(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3) + pp(SS1:NN1,SS2:NN2,SS3:NN3)
        
            else
        
                phi(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3) + pp(SS1:NN1,SS2:NN2,SS3:NN3)
        
                if (problem_type == 1) call product_Helmholtz      (phi,Ap)
                if (problem_type == 2) call product_div_grad       (phi,Ap)
                if (problem_type == 3) call product_Helmholtz_conc (phi,Ap)
                if (problem_type == 4) call product_div_grad_transp(phi,Ap)
        
                rr(SS1:NN1,SS2:NN2,SS3:NN3) = bb(SS1:NN1,SS2:NN2,SS3:NN3) - Ap(SS1:NN1,SS2:NN2,SS3:NN3)
            end if
     
     
            !--- grösstes Residuum --------------------------------------
            call get_norms(SS1,SS2,SS3,NN1,NN2,NN3,rr,problem_type,.true.,.false.,normInf=norm_inf)
     
     
            !--- Konvergenzstatistik ---------------------------------
            if (problem_type == 1) ratioH(substep,direction     ) = ratioH(substep,direction     ) + LOG10(norm_inf/norm_inf0)
            if (problem_type == 2) ratioP(substep,number_poisson) = ratioP(substep,number_poisson) + LOG10(norm_inf/norm_inf0)
     
            if (problem_type == 1) countH(substep,direction     ) = countH(substep,direction     ) + 1
            if (problem_type == 2) countP(substep,number_poisson) = countP(substep,number_poisson) + 1
     
            norm_inf0 = norm_inf
     
        end do ITERATE
  
  
    end subroutine Richardson
  
  
  
  
  
  
  
  
  
  
  
    subroutine product_scalar(SS1,SS2,SS3,NN1,NN2,NN3,phi1,phi2,scalar)
  
        implicit none
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        real   , intent(out)   ::  scalar
        real                   ::  scalar_global
        integer                ::  i, j, k
  
  
        scalar = 0.
  
        do k = SS3, NN3
            do j = SS2, NN2
                !pgi$ unroll = n:8
                do i = SS1, NN1
                    scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
                end do
            end do
        end do
  
        call MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
        scalar = scalar_global
  
  
    end subroutine product_scalar
  
  
  
  
  
  
  
  
  
  
    ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
    subroutine product_scalar2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,phi1,phi2,scalar)
  
        implicit none
  
        integer, intent(in)    ::  g
  
        integer, intent(in)    ::  N1
        integer, intent(in)    ::  N2
        integer, intent(in)    ::  N3
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        real   , intent(out)   ::  scalar
        real                   ::  scalar_global
        integer                ::  i, j, k
  
  
        scalar = 0.
  
        do k = SS3, NN3
            do j = SS2, NN2
                !pgi$ unroll = n:8
                do i = SS1, NN1
                    scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
                end do
            end do
        end do
  
        call MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
        scalar = scalar_global
  
  
    end subroutine product_scalar2
  
  
  
  
  
  
  
  
  
  
  
    subroutine multadd1(SS1,SS2,SS3,NN1,NN2,NN3,mult,vec,phi)
  
        implicit none
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  mult
  
        real   , intent(inout) ::  vec(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer                ::  i, j, k
  
  
        do k = SS3, NN3
            do j = SS2, NN2
                !pgi$ unroll = n:8
                do i = SS1, NN1
                    vec(i,j,k) = vec(i,j,k) + mult*phi(i,j,k)
                end do
            end do
        end do
  
  
    end subroutine multadd1
  
  
  
  
  
  
  
  
  
  
  
    subroutine multadd2(SS1,SS2,SS3,NN1,NN2,NN3,mult1,mult2,vec,phi1,phi2,init_yes)
  
        implicit none
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  mult1
        real   , intent(in)    ::  mult2
  
        real   , intent(inout) ::  vec (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical, intent(in)    ::  init_yes
  
        integer                ::  i, j, k
  
  
        if (init_yes) then
     
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        vec(i,j,k) = mult1*phi1(i,j,k) + mult2*phi2(i,j,k)
                    end do
                end do
            end do
     
        else
     
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        vec(i,j,k) = vec(i,j,k) + mult1*phi1(i,j,k) + mult2*phi2(i,j,k)
                    end do
                end do
            end do
     
        end if
  
  
    end subroutine multadd2
  
  
  
  
  
  
  
  
  
  
  
    !> brief computes two or infinity norm( get is misleading)
    !!
    !! ??? shouldn't it be defined in an own module
    !! \param[in] SS1 start index in 1-direction
    !! \param[in] SS2 start index in 1-direction
    !! \param[in] SS3 start index in 1-direction
    !! \param[in] NN1 end index in 1-direction
    !! \param[in] NN2 end index in 2-direction
    !! \param[in] NN3 end index in 3-direction
    !! \param[in] phi vector, from which the norm is taken
    !! \param[in] problem_type if problem_type==2(Helmholtz) and weighted_yes the weighting is used
    !! \param[in] inf_yes if true infinity norm is computed
    !! \param[in] two_yes if trhue two norm is computed
    !! \param[out] normInf gets the infinity norm of phi
    !! \param[out] normTwo get the two norm of phi
    subroutine get_norms(SS1,SS2,SS3,NN1,NN2,NN3,phi,problem_type,inf_yes,two_yes,normInf,normTwo)
  
        implicit none
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer, intent(in)    ::  problem_type
        logical, intent(in)    ::  inf_yes
        logical, intent(in)    ::  two_yes
  
        real   , optional, intent(out) ::  normInf
        real   , optional, intent(out) ::  normTwo
  
        real                   ::  normInf_global, normTwo_global
        integer                ::  i, j, k
  
  
        if (inf_yes .and. two_yes) then
     
            normInf = 0.
            normTwo = 0.
     
            if (problem_type == 2 .and. weighting_yes) then
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                            normTwo = normTwo + phi(i,j,k)**2
                        end do
                    end do
                end do
        
            else
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)),normInf)
                            normTwo = normTwo + phi(i,j,k)**2
                        end do
                    end do
                end do
        
            end if
     
            ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normInf = normInf_global
            normTwo = normTwo_global
     
        else if (inf_yes) then
     
            normInf = 0.
     
            if (problem_type == 2 .and. weighting_yes) then
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                        end do
                    end do
                end do
        
            else
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)),normInf)
                        end do
                    end do
                end do
        
            end if
     
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            normInf = normInf_global
     
        else if (two_yes) then
     
            normTwo = 0.
     
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        normTwo = normTwo + phi(i,j,k)**2
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normTwo = normTwo_global
     
        end if
  
  
    end subroutine get_norms
  
  
  
  
  
  
  
  
  
  
    ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
    subroutine get_norms2(g,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,phi,problem_type,inf_yes,two_yes,normInf,normTwo)
  
        implicit none
  
        integer, intent(in)    ::  g
  
        integer, intent(in)    ::  N1
        integer, intent(in)    ::  N2
        integer, intent(in)    ::  N3
  
        integer, intent(in)    ::  SS1
        integer, intent(in)    ::  SS2
        integer, intent(in)    ::  SS3
  
        integer, intent(in)    ::  NN1
        integer, intent(in)    ::  NN2
        integer, intent(in)    ::  NN3
  
        real   , intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer, intent(in)    ::  problem_type
        logical, intent(in)    ::  inf_yes
        logical, intent(in)    ::  two_yes
  
        real   , optional, intent(out) ::  normInf
        real   , optional, intent(out) ::  normTwo
  
        real                   ::  normInf_global, normTwo_global
        integer                ::  i, j, k
  
  
        if (inf_yes .and. two_yes) then
     
            normInf = 0.
            normTwo = 0.
     
            if (problem_type == 2 .and. weighting_yes) then
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                            normTwo = normTwo + phi(i,j,k)**2
                        end do
                    end do
                end do
        
            else
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)),normInf)
                            normTwo = normTwo + phi(i,j,k)**2
                        end do
                    end do
                end do
        
            end if
     
            ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror)
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
            normInf = normInf_global
            normTwo = normTwo_global
     
        else if (inf_yes) then
     
            normInf = 0.
     
            if (problem_type == 2 .and. weighting_yes) then
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                        end do
                    end do
                end do
        
            else
        
                do k = SS3, NN3
                    do j = SS2, NN2
                        !pgi$ unroll = n:8
                        do i = SS1, NN1
                            normInf = MAX(ABS(phi(i,j,k)),normInf)
                        end do
                    end do
                end do
        
            end if
     
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm1(g),merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            normInf = normInf_global
     
        else if (two_yes) then
     
            normTwo = 0.
     
            do k = SS3, NN3
                do j = SS2, NN2
                    !pgi$ unroll = n:8
                    do i = SS1, NN1
                        normTwo = normTwo + phi(i,j,k)**2
                    end do
                end do
            end do
     
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm1(g),merror)
            normTwo = normTwo_global
     
        end if
  
  
    end subroutine get_norms2
  
  
  
  
  
  
  
  
  
  
  
    subroutine status_iteration(eps,norm,counter,N_restarts,exit_yes,quiet_yes1,quiet_yes2)
  
        implicit none
  
        real   , intent(in)    ::  eps
        real   , intent(in)    ::  norm
        integer, intent(in)    ::  counter
        integer, intent(in)    ::  N_restarts
        logical, intent(out)   ::  exit_yes
        logical, intent(in)    ::  quiet_yes1
        logical, intent(in)    ::  quiet_yes2
  
  
        exit_yes = .false.
  
        if (norm < eps) then
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
            if (rank == 0 .and. log_iteration_yes                     ) write(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  (Termination criterion satisfied)'
            exit_yes = .true.
        end if
  
        if ((.not. exit_yes) .and. counter == N_restarts) then
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes2) write(* ,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
            if (rank == 0 .and. log_iteration_yes                     ) write(10,'(a,i5,a,E24.17,a)') '  Iteration',counter,'; ||res|| =',norm,'  WARNING! Too many iterations!'
            exit_yes = .true.
        end if
  
        if (.not. exit_yes) then
            if (rank == 0 .and. write_stout_yes .and. .not. quiet_yes1) write(* ,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
            if (rank == 0 .and. log_iteration_yes                     ) write(10,'(a,i5,a,E24.17  )') '  Iteration',counter,'; ||res|| =',norm
        end if
  
  
    end subroutine status_iteration
  
  
  
  
  
  
  
  
  
  
  
    !> Anmerkungen: Siehe Subroutine "handle_corner_Lap"!                                                     !
    subroutine handle_corner_rhs(g,phi) ! TEST!!! Ersetzen mit einheitlicher Routine?? ! TEST!!! validieren ...
  
        implicit none
  
        integer, intent(in   ) ::  g
        real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Siehe Subroutine "handle_corner_Lap"!                                                     !
        !----------------------------------------------------------------------------------------------------------!
  
  
        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) phi(1      ,1      ,1:NN(3,g)) = 0.
        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) phi(1      ,NN(2,g),1:NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) phi(NN(1,g),1      ,1:NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) phi(NN(1,g),NN(2,g),1:NN(3,g)) = 0.
  
        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) phi(1      ,1:NN(2,g),1      ) = 0.
        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) phi(1      ,1:NN(2,g),NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) phi(NN(1,g),1:NN(2,g),1      ) = 0.
        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) phi(NN(1,g),1:NN(2,g),NN(3,g)) = 0.
  
        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) phi(1:NN(1,g),1      ,1      ) = 0.
        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) phi(1:NN(1,g),1      ,NN(3,g)) = 0.
        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) phi(1:NN(1,g),NN(2,g),1      ) = 0.
        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) phi(1:NN(1,g),NN(2,g),NN(3,g)) = 0.
  
  
    end subroutine handle_corner_rhs
  
  
  
  
  
  
  
  
  
  
  
    subroutine handle_corner_rhs_conc(phi)
  
        implicit none
  
        real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Siehe Subroutines "handle_corner_Lap" und "handle_corner_rhs"!                            !
        !----------------------------------------------------------------------------------------------------------!
  
  
        if (BCc_1L(conc_nu) > 0 .and. BCc_2L(conc_nu) > 0) phi(1 ,1 ,1:N3) = 0.
        if (BCc_1L(conc_nu) > 0 .and. BCc_2U(conc_nu) > 0) phi(1 ,N2,1:N3) = 0.
        if (BCc_1U(conc_nu) > 0 .and. BCc_2L(conc_nu) > 0) phi(N1,1 ,1:N3) = 0.
        if (BCc_1U(conc_nu) > 0 .and. BCc_2U(conc_nu) > 0) phi(N1,N2,1:N3) = 0.
  
        if (BCc_1L(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1 ,1:N2,1 ) = 0.
        if (BCc_1L(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1 ,1:N2,N3) = 0.
        if (BCc_1U(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(N1,1:N2,1 ) = 0.
        if (BCc_1U(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(N1,1:N2,N3) = 0.
  
        if (BCc_2L(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1:N1,1 ,1 ) = 0.
        if (BCc_2L(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1:N1,1 ,N3) = 0.
        if (BCc_2U(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1:N1,N2,1 ) = 0.
        if (BCc_2U(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1:N1,N2,N3) = 0.
  
  
    end subroutine handle_corner_rhs_conc
  
    !> brief computes two or infinity norm( get is misleading)
    !!
    !! ??? shouldn't it be defined in an own module
    !! \param[in] phi velocity vector, from which the norm is taken
    !! \param[in] inf_yes if true infinity norm is computed
    !! \param[in] two_yes if trhue two norm is computed
    !! \param[out] normInf gets the infinity norm of phi
    !! \param[out] normTwo get the two norm of phi
    subroutine get_norms_vel(phi,inf_yes,two_yes,normInf,normTwo)

        implicit none

        real   , intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)

        logical, intent(in)    ::  inf_yes
        logical, intent(in)    ::  two_yes

        real   , optional, intent(out) ::  normInf
        real   , optional, intent(out) ::  normTwo

        real                   ::  normInf_global, normTwo_global
        integer                ::  i, j, k


        if (inf_yes .and. two_yes) then

            normInf = 0.
            normTwo = 0.

            do k = S31, N31
                do j = S21, N21
                    !pgi$ unroll = n:8
                    do i = S11, N11
                        normInf = MAX(ABS(phi(i,j,k,1)),normInf)
                        normTwo = normTwo + phi(i,j,k,1)**2
                    end do
                end do
            end do
            do k = S32, N32
                do j = S22, N22
                    !pgi$ unroll = n:8
                    do i = S12, N12
                        normInf = MAX(ABS(phi(i,j,k,2)),normInf)
                        normTwo = normTwo + phi(i,j,k,2)**2
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = S33, N33
                    do j = S23, N23
                        !pgi$ unroll = n:8
                        do i = S13, N13
                            normInf = MAX(ABS(phi(i,j,k,3)),normInf)
                            normTwo = normTwo + phi(i,j,k,3)**2
                        end do
                    end do
                end do
            end if

            ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normInf = normInf_global
            normTwo = normTwo_global

        else if (inf_yes) then

            normInf = 0.

            do k = S31, N31
                do j = S21, N21
                    !pgi$ unroll = n:8
                    do i = S11, N11
                        normInf = MAX(ABS(phi(i,j,k,1)),normInf)
                    end do
                end do
            end do
            do k = S32, N32
                do j = S22, N22
                    !pgi$ unroll = n:8
                    do i = S12, N12
                        normInf = MAX(ABS(phi(i,j,k,2)),normInf)
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = S33, N33
                    do j = S23, N23
                        !pgi$ unroll = n:8
                        do i = S13, N13
                            normInf = MAX(ABS(phi(i,j,k,3)),normInf)
                        end do
                    end do
                end do
            end if


            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            normInf = normInf_global

        else if (two_yes) then

            normTwo = 0.

            do k = S31, N31
                do j = S21, N21
                    !pgi$ unroll = n:8
                    do i = S11, N11
                        normTwo = normTwo + phi(i,j,k,1)**2
                    end do
                end do
            end do
            do k = S32, N32
                do j = S22, N22
                    !pgi$ unroll = n:8
                    do i = S12, N12
                        normTwo = normTwo + phi(i,j,k,2)**2
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = S33, N33
                    do j = S23, N23
                        !pgi$ unroll = n:8
                        do i = S13, N13
                            normTwo = normTwo + phi(i,j,k,3)**2
                        end do
                    end do
                end do
            end if

            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normTwo = normTwo_global

        end if


    end subroutine get_norms_vel
  
  
end module cmod_solvers
