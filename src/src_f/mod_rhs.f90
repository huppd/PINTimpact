!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief provides routines for computing the right hand side for the velocity and the concetration
module mod_rhs
  
  
  use mod_dims
  use mod_vars
  use mod_diff
  use mod_particles
  use mod_les
  
  
  private
  
  public rhs_vel, rhs_conc
  
  INCLUDE 'mpif.h'
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> \brief computes right hand side for the velocity component
  !!
  !! \f[ \mathrm{ rhs = vel - (-(1-thetaL)(aRK(substep) + bRK(substep))\frac{\Delta t}{Re} \Delta vel } \f]
  !! \f[ \mathrm{ rhs = rhs -\Delta t \, bRK(substep) nl_{old} } \f]
  !!
  !! \f[ \mathrm{ rhs = vel - (-(1-thetaL)(aRK(substep) + bRK(substep))\frac{\Delta t}{Re} \Delta vel -\Delta t \, bRK(substep) nl_{old} - \Delta t \, aRK(substep) nl} \f]
  subroutine rhs_vel
  
  implicit none
  
  integer                ::  i, j, k, m
  real                   ::  mult
  real                   ::  flux, flux_global
  
  
  ! TEST!!! generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(vel) ============================================================================
  !===========================================================================================================
  if (timeint_mode == 1 .or. thetaL == 1.) then
     do m = 1, dimens
        if (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
        if (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        if (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
#ifdef NONBOUSSINESQ
        ! TEST!!! s.o.: generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
        ! ==> doppelte Arbeit an der Stelle!
        if (nonBoussinesq_yes) then
           if (m == 1) then
              do k = S31, N31
                 do j = S21, N21
!pgi$ unroll = n:8
                    do i = S11, N11
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    end do
                 end do
              end do
           end if
           if (m == 2) then
              do k = S32, N32
                 do j = S22, N22
!pgi$ unroll = n:8
                    do i = S12, N12
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    end do
                 end do
              end do
           end if
           if (m == 3) then
              do k = S33, N33
                 do j = S23, N23
!pgi$ unroll = n:8
                    do i = S13, N13
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    end do
                 end do
              end do
           end if
        end if
#endif
     end do
  else
     multL = -(1.-thetaL)*(aRK(substep)+bRK(substep))*dtime/Re
     do m = 1, dimens
        call Helmholtz(m,.false.,vel(b1L,b2L,b3L,m),rhs(b1L,b2L,b3L,m))
     end do
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  if (substep /= 1) then
     mult = dtime*bRK(substep)
     
     do m = 1, dimens
        if (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
        if (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
        if (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
     end do
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- advective terms ---------------------------------------------------------------------------------------
  if (Stokes_yes) then
                      nl(S11B:N11B,S21B:N21B,S31B:N31B,1) = 0.
                      nl(S12B:N12B,S22B:N22B,S32B:N32B,2) = 0.
     if (dimens == 3) nl(S13B:N13B,S23B:N23B,S33B:N33B,3) = 0.
  else
     call nonlinear(.false.)
  end if
  
#ifdef NONBOUSSINESQ
  !--- multiply by total fluid density -----------------------------------------------------------------------
  if (.not. Stokes_yes .and. nonBoussinesq_yes) then
     do k = S31B, N31B
        do j = S21B, N21B
!pgi$ unroll = n:8
           do i = S11B, N11B
              nl(i,j,k,1) = nl(i,j,k,1)*dens(i,j,k,1)
           end do
        end do
     end do
     do k = S32B, N32B
        do j = S22B, N22B
!pgi$ unroll = n:8
           do i = S12B, N12B
              nl(i,j,k,2) = nl(i,j,k,2)*dens(i,j,k,2)
           end do
        end do
     end do
     if (dimens == 3) then
     do k = S33B, N33B
        do j = S23B, N23B
!pgi$ unroll = n:8
           do i = S13B, N13B
              nl(i,j,k,3) = nl(i,j,k,3)*dens(i,j,k,3)
           end do
        end do
     end do
     end if
  end if
  
  !--- transport of concentration diffusion ------------------------------------------------------------------
  if (timeint_mode == 1 .and. (.not. Euler_yes) .and. nonBoussinesq_yes) then
     res = 0.
     do conc_nu = 1, n_conc
        m = conc_nu
        multL = Fro*Ric(m)/(Re*Sc(m))
        call Helmholtz_conc_explicit(.false.,res)
     end do
     
     call interpolate2_pre_vel(.true.,1,res,gpre)
     do k = S31B, N31B
        do j = S21B, N21B
!pgi$ unroll = n:8
           do i = S11B, N11B
              nl(i,j,k,1) = nl(i,j,k,1) - vel(i,j,k,1)*gpre(i,j,k)
           end do
        end do
     end do
     
     call interpolate2_pre_vel(.true.,2,res,gpre)
     do k = S32B, N32B
        do j = S22B, N22B
!pgi$ unroll = n:8
           do i = S12B, N12B
              nl(i,j,k,2) = nl(i,j,k,2) - vel(i,j,k,2)*gpre(i,j,k)
           end do
        end do
     end do
     
     if (dimens == 3) then
     call interpolate2_pre_vel(.true.,3,res,gpre)
     do k = S33B, N33B
        do j = S23B, N23B
!pgi$ unroll = n:8
           do i = S13B, N13B
              nl(i,j,k,3) = nl(i,j,k,3) - vel(i,j,k,3)*gpre(i,j,k)
           end do
        end do
     end do
     end if
  end if
#endif
  
  !--- viscose terms -----------------------------------------------------------------------------------------
  if (timeint_mode == 1 .and. (.not. Euler_yes)) then
     multL = 1./Re
     call Helmholtz_explicit(.false.)
  end if
  
  !--- feedback-forces concentrations ------------------------------------------------------------------------
  if (concentration_yes) call interpolate_conc(.false.)
  
  !--- feedback-forces Lagrangian particles ------------------------------------------------------------------
  if (particles_yes) call move_particles(.false.) ! TEST!!! Eine dauerhaftere Loesung suchen, die auch die Kopplung mit conc elegant zulaesst! TEST!!! umbenennen ...
  
  !--- LES closure -------------------------------------------------------------------------------------------
  if (LES_mode == 1) then
     call admrt(1,0,S11B,S21B,S31B,N11B,N21B,N31B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,1),nl(b1L,b2L,b3L,1))
     call admrt(2,0,S12B,S22B,S32B,N12B,N22B,N32B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,2),nl(b1L,b2L,b3L,2))
     call admrt(3,0,S13B,S23B,S33B,N13B,N23B,N33B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,3),nl(b1L,b2L,b3L,3))
  end if
  if (LES_mode == 2) then
     call smag_vel(.false.)
  end if
  
  !--- forcing -----------------------------------------------------------------------------------------------
  call forcing_vel
  
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  
  do m = 1, dimens
     if (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
     if (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
     if (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
  end do
  !===========================================================================================================
  
  
  
  ! TEST!!! BC-Abfrage unklar bzw. inkonsistent ...
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Momentan sind nur Dirichlet-Randbedingungen für Geschwindigkeiten wählbar.                !
  !              - "vel" darf NICHT überschrieben werden, da sonst die Konzentrationen nicht korrekt berech- !
  !                net werden können.                                                                        !
  !----------------------------------------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------------------------------------!
  ! ACHTUNG: Es ist sehr genau auf die "Stetigkeit" der RB über Ecken/Kanten zu achten, bzw. die Reihenfolge !
  !          der Operationen und die Intervallgrenzen müssen UNBEDINGT sorgfaltig gewählt werden (Grund:     !
  !          Divergenz wird auch auf dem Rand gebildet). Siehe dazu auch die Reihenfolge der RB in           !
  !          product_Helmholtz!                                                                              !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Randbedingungen (Zeitintegration) =====================================================================
  !===========================================================================================================
  if (substep /= 1) then
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl_old ---------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*bRK(substep)
     !--------------------------------------------------------------------------------------------------------
     if (outlet(1,1,1) .and. (BC_1L == 1 .or. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
     if (outlet(1,2,1) .and. (BC_1U == 1 .or. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
     
     if (outlet(1,1,2) .and. (BC_1L == 1 .or. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
     if (outlet(1,2,2) .and. (BC_1U == 1 .or. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
     
     if (outlet(1,1,3) .and. (BC_1L == 1 .or. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
     if (outlet(1,2,3) .and. (BC_1U == 1 .or. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     if (outlet(2,1,1) .and. (BC_2L == 1 .or. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
     if (outlet(2,2,1) .and. (BC_2U == 1 .or. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
     
     if (outlet(2,1,2) .and. (BC_2L == 1 .or. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
     if (outlet(2,2,2) .and. (BC_2U == 1 .or. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
     
     if (outlet(2,1,3) .and. (BC_2L == 1 .or. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
     if (outlet(2,2,3) .and. (BC_2U == 1 .or. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     if (outlet(3,1,1) .and. (BC_3L == 1 .or. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
     if (outlet(3,2,1) .and. (BC_3U == 1 .or. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
     
     if (outlet(3,1,2) .and. (BC_3L == 1 .or. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
     if (outlet(3,2,2) .and. (BC_3U == 1 .or. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
     
     if (outlet(3,1,3) .and. (BC_3L == 1 .or. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
     if (outlet(3,2,3) .and. (BC_3U == 1 .or. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
     !--------------------------------------------------------------------------------------------------------
  end if
  
  !--- Ausfluss-RB -------------------------------------------------------------------------------------------
  call outflow_bc
  
  !--- andere RB (instationaer + Zeitintegration) ------------------------------------------------------------
  call boundary_vel_tint
  
  !-----------------------------------------------------------------------------------------------------------
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(1,1,1) .and. (BC_1L == 1 .or. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
  if (outlet(1,2,1) .and. (BC_1U == 1 .or. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
  
  if (outlet(1,1,2) .and. (BC_1L == 1 .or. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
  if (outlet(1,2,2) .and. (BC_1U == 1 .or. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
  
  if (outlet(1,1,3) .and. (BC_1L == 1 .or. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
  if (outlet(1,2,3) .and. (BC_1U == 1 .or. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(2,1,1) .and. (BC_2L == 1 .or. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
  if (outlet(2,2,1) .and. (BC_2U == 1 .or. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
  
  if (outlet(2,1,2) .and. (BC_2L == 1 .or. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
  if (outlet(2,2,2) .and. (BC_2U == 1 .or. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
  
  if (outlet(2,1,3) .and. (BC_2L == 1 .or. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
  if (outlet(2,2,3) .and. (BC_2U == 1 .or. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,1) .and. (BC_3L == 1 .or. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
  if (outlet(3,2,1) .and. (BC_3U == 1 .or. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
  
  if (outlet(3,1,2) .and. (BC_3L == 1 .or. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
  if (outlet(3,2,2) .and. (BC_3U == 1 .or. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
  
  if (outlet(3,1,3) .and. (BC_3L == 1 .or. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
  if (outlet(3,2,3) .and. (BC_3U == 1 .or. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen (instantan) ===========================================================================
  !===========================================================================================================
  call boundary_vel_stat
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== RB einsetzen ==========================================================================================
  !===========================================================================================================
  ! Anmerkung: Reihenfolge bestimmt Konvention (Linien, Kanten)!
  !-----------------------------------------------------------------------------------------------------------
  if (BC_2L > 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
  if (BC_2U > 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
  
  if (BC_3L > 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
  if (BC_3U > 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
  
  if (BC_1L > 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
  if (BC_1U > 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
  !-----------------------------------------------------------------------------------------------------------
  if (BC_1L > 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
  if (BC_1U > 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
  
  if (BC_3L > 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
  if (BC_3U > 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
  
  if (BC_2L > 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
  if (BC_2U > 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
  !-----------------------------------------------------------------------------------------------------------
  if (BC_1L > 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
  if (BC_1U > 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
  
  if (BC_2L > 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
  if (BC_2U > 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
  
  if (BC_3L > 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
  if (BC_3U > 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== Fluss-Korrektur =======================================================================================
  !===========================================================================================================
  if (nullspace_yes) then
     flux = 0.
     
     m = 1
     do k = S31B, N31B
        do j = S21B, N21B
!pgi$ unroll = n:8
           do i = S11B, N11B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           end do
        end do
     end do
     
     m = 2
     do k = S32B, N32B
        do j = S22B, N22B
!pgi$ unroll = n:8
           do i = S12B, N12B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           end do
        end do
     end do
     
     if (dimens == 3) then
     m = 3
     do k = S33B, N33B
        do j = S23B, N23B
!pgi$ unroll = n:8
           do i = S13B, N13B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           end do
        end do
     end do
     end if
     
     
     call MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     flux = flux_global
     
     if (rank == 0) write(10,'(a,E25.17)') 'flux =', flux
     
     
     ! d/dt ist unstetig in der Zeit ==> kein Versuch unternommen, in Zeitintegration einzubauen ...
     if (nullspace_ortho_yes) then
        
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L > 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,1        ,S31B:N31B,1)
        if (BC_2U > 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*psi_vel(S11B:N11B,N2       ,S31B:N31B,1)
        
        if (BC_3L > 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,1        ,1)
        if (BC_3U > 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*psi_vel(S11B:N11B,S21B:N21B,N3       ,1)
        
        if (BC_1L > 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*psi_vel(0        ,S21B:N21B,S31B:N31B,1)
        if (BC_1U > 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*psi_vel(N1       ,S21B:N21B,S31B:N31B,1)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*psi_vel(1        ,S22B:N22B,S32B:N32B,2)
        if (BC_1U > 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*psi_vel(N1       ,S22B:N22B,S32B:N32B,2)
        
        if (BC_3L > 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*psi_vel(S12B:N12B,S22B:N22B,1        ,2)
        if (BC_3U > 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,N3       ,2)
        
        if (BC_2L > 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*psi_vel(S12B:N12B,0        ,S32B:N32B,2)
        if (BC_2U > 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,N2       ,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*psi_vel(1        ,S23B:N23B,S33B:N33B,3)
        if (BC_1U > 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*psi_vel(N1       ,S23B:N23B,S33B:N33B,3)
        
        if (BC_2L > 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*psi_vel(S13B:N13B,1        ,S33B:N33B,3)
        if (BC_2U > 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*psi_vel(S13B:N13B,N2       ,S33B:N33B,3)
        
        if (BC_3L > 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*psi_vel(S13B:N13B,S23B:N23B,0        ,3)
        if (BC_3U > 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*psi_vel(S13B:N13B,S23B:N23B,N3       ,3)
        !-----------------------------------------------------------------------------------------------------
        
                         rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
                         rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        if (dimens == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - flux*psi_vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
     else
        
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L > 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*th12(S11B:N11B,S31B:N31B,1)
        if (BC_2U > 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*th12(S11B:N11B,S31B:N31B,2)
        
        if (BC_3L > 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*th13(S11B:N11B,S21B:N21B,1)
        if (BC_3U > 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*th13(S11B:N11B,S21B:N21B,2)
        
        if (BC_1L > 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*th11(S21B:N21B,S31B:N31B,1)
        if (BC_1U > 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*th11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*th21(S22B:N22B,S32B:N32B,1)
        if (BC_1U > 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*th21(S22B:N22B,S32B:N32B,2)
        
        if (BC_3L > 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*th23(S12B:N12B,S22B:N22B,1)
        if (BC_3U > 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*th23(S12B:N12B,S22B:N22B,2)
        
        if (BC_2L > 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*th22(S12B:N12B,S32B:N32B,1)
        if (BC_2U > 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*th22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*th31(S23B:N23B,S33B:N33B,1)
        if (BC_1U > 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*th31(S23B:N23B,S33B:N33B,2)
        
        if (BC_2L > 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*th32(S13B:N13B,S33B:N33B,1)
        if (BC_2U > 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*th32(S13B:N13B,S33B:N33B,2)
        
        if (BC_3L > 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*th33(S13B:N13B,S23B:N23B,1)
        if (BC_3U > 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*th33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
        
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L > 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
        if (BC_2U > 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
        
        if (BC_3L > 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
        if (BC_3U > 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
        
        if (BC_1L > 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
        if (BC_1U > 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
        if (BC_1U > 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
        
        if (BC_3L > 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
        if (BC_3U > 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
        
        if (BC_2L > 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
        if (BC_2U > 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
        if (BC_1U > 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
        
        if (BC_2L > 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
        if (BC_2U > 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
        
        if (BC_3L > 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
        if (BC_3U > 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
     end if
     
  end if
  !===========================================================================================================
  
  
  end subroutine rhs_vel
  
  
  
  
  
  
  
  
  
  
  
  subroutine rhs_conc
  
  implicit none
  
  integer                ::  m
  real                   ::  mult
  
  
  m = conc_nu
  
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(conc) ===========================================================================
  !===========================================================================================================
  if (timeint_mode == 1 .or. thetaL == 1.) then
     res(S1p:N1p,S2p:N2p,S3p:N3p) = conc(S1p:N1p,S2p:N2p,S3p:N3p,m)
  else
     multL = -(1.-thetaL)*(aRK(substep)+bRK(substep))*dtime / (Re*Sc(m))
     call Helmholtz_conc(.false.,conc(b1L,b2L,b3L,m),res)
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  if (substep /= 1) then
     mult = dtime*bRK(substep)
     res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - mult*nlco(S1p:N1p,S2p:N2p,S3p:N3p,m)
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- Konvektions-Term --------------------------------------------------------------------------------------
  ! TEST!!! Stokes flow?? Evtl. analog fluid auch hier einbauen ...
  call nonlinear_conc(.false.)
  
  !--- diffusive term ----------------------------------------------------------------------------------------
  if (timeint_mode == 1 .and. (.not. Euler_yes)) then
     multL = 1./(Re*Sc(m))
     call Helmholtz_conc_explicit(.false.,nlco(b1L,b2L,b3L,m))
  end if
  
  !--- LES closure -------------------------------------------------------------------------------------------
  if (LES_mode == 1) call admrt(0,m,S1p,S2p,S3p,N1p,N2p,N3p,n_lp_conc(m),n_hp_conc(m),chi_conc(m),conc(b1L,b2L,b3L,m),nlco(b1L,b2L,b3L,m))
  if (LES_mode == 2) call smag_conc(.false.)
  
  !--- forcing -----------------------------------------------------------------------------------------------
  call forcing_conc
  
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - mult*nlco(S1p:N1p,S2p:N2p,S3p:N3p,m)
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== Randbedingungen (Zeitintegration) =====================================================================
  !===========================================================================================================
  ! - sed1L & Co. k�nnten auch in nlco gespeichert werden. So ist es allerdings sicherer!
  ! - evtl. bc-Felder einf�hren und analog zu den Geschwindigkeiten verfahren ...
  
  if (substep == 1) then
     !--- Ausfluss-RB ----------------------------------------------------------------------------------------
     call sediment_bc
     
     !--- andere RB (instationaer + Zeitintegration) ---------------------------------------------------------
     call boundary_conc_tint
     
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl -------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*aRK(substep)
     
     if (BCc_1L(m) == 1 .and. (.not. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = conc(1 ,S2p:N2p,S3p:N3p,m) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     if (BCc_1U(m) == 1 .and. (.not. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = conc(N1,S2p:N2p,S3p:N3p,m) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     if (BCc_2L(m) == 1 .and. (.not. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = conc(S1p:N1p,1 ,S3p:N3p,m) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     if (BCc_2U(m) == 1 .and. (.not. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = conc(S1p:N1p,N2,S3p:N3p,m) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     if (BCc_3L(m) == 1 .and. (.not. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = conc(S1p:N1p,S2p:N2p,1 ,m) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     if (BCc_3U(m) == 1 .and. (.not. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = conc(S1p:N1p,S2p:N2p,N3,m) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
  else
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl_old ---------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*bRK(substep)
     
     if (BCc_1L(m) == 1 .and. (.not. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = conc(1 ,S2p:N2p,S3p:N3p,m) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     if (BCc_1U(m) == 1 .and. (.not. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = conc(N1,S2p:N2p,S3p:N3p,m) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     if (BCc_2L(m) == 1 .and. (.not. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = conc(S1p:N1p,1 ,S3p:N3p,m) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     if (BCc_2U(m) == 1 .and. (.not. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = conc(S1p:N1p,N2,S3p:N3p,m) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     if (BCc_3L(m) == 1 .and. (.not. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = conc(S1p:N1p,S2p:N2p,1 ,m) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     if (BCc_3U(m) == 1 .and. (.not. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = conc(S1p:N1p,S2p:N2p,N3,m) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
     
     !--- Ausfluss-RB ----------------------------------------------------------------------------------------
     call sediment_bc
     
     !--- andere RB (instationaer + Zeitintegration) ---------------------------------------------------------
     call boundary_conc_tint
     
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl -------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*aRK(substep)
     
     if (BCc_1L(m) == 1 .and. (.not. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = res(1 ,S2p:N2p,S3p:N3p) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     if (BCc_1U(m) == 1 .and. (.not. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = res(N1,S2p:N2p,S3p:N3p) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     if (BCc_2L(m) == 1 .and. (.not. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = res(S1p:N1p,1 ,S3p:N3p) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     if (BCc_2U(m) == 1 .and. (.not. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = res(S1p:N1p,N2,S3p:N3p) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     if (BCc_3L(m) == 1 .and. (.not. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = res(S1p:N1p,S2p:N2p,1 ) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     if (BCc_3U(m) == 1 .and. (.not. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = res(S1p:N1p,S2p:N2p,N3) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Randbedingungen (instantan) ===========================================================================
  !===========================================================================================================
  call boundary_conc_stat
  !===========================================================================================================
  
  
  end subroutine rhs_conc
  
  
  
  
end module mod_rhs
