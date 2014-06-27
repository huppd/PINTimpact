!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief provides routines for computing the right hand side for the velocity and the concetration
MODULE mod_rhs
  
  
  USE mod_dims
  USE mod_vars
  USE mod_diff
  USE mod_particles
  USE mod_les
  
  
  PRIVATE
  
  PUBLIC rhs_vel, rhs_conc
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> \brief computes right hand side for the velocity component
  !!
  !! \f[ \mathrm{ rhs = vel - (-(1-thetaL)(aRK(substep) + bRK(substep))\frac{\Delta t}{Re} \Delta vel } \f]
  !! \f[ \mathrm{ rhs = rhs -\Delta t \, bRK(substep) nl_{old} } \f]
  !!
  !! \f[ \mathrm{ rhs = vel - (-(1-thetaL)(aRK(substep) + bRK(substep))\frac{\Delta t}{Re} \Delta vel -\Delta t \, bRK(substep) nl_{old} - \Delta t \, aRK(substep) nl} \f]
  SUBROUTINE rhs_vel
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL                   ::  mult
  REAL                   ::  flux, flux_global
  
  
  ! TEST!!! generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(vel) ============================================================================
  !===========================================================================================================
  IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
#ifdef NONBOUSSINESQ
        ! TEST!!! s.o.: generell muesste eigentlich nur der "innere" Bereich (d.h. ohne Raender) von rhs" belegt werden!
        ! ==> doppelte Arbeit an der Stelle!
        IF (nonBoussinesq_yes) THEN
           IF (m == 1) THEN
              DO k = S31, N31
                 DO j = S21, N21
!pgi$ unroll = n:8
                    DO i = S11, N11
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    END DO
                 END DO
              END DO
           END IF
           IF (m == 2) THEN
              DO k = S32, N32
                 DO j = S22, N22
!pgi$ unroll = n:8
                    DO i = S12, N12
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    END DO
                 END DO
              END DO
           END IF
           IF (m == 3) THEN
              DO k = S33, N33
                 DO j = S23, N23
!pgi$ unroll = n:8
                    DO i = S13, N13
                       rhs(i,j,k,m) = vel(i,j,k,m)*dens(i,j,k,m)
                    END DO
                 END DO
              END DO
           END IF
        END IF
#endif
     END DO
  ELSE
     multL = -(1.-thetaL)*(aRK(substep)+bRK(substep))*dtime/Re
     DO m = 1, dimens
        CALL Helmholtz(m,.FALSE.,vel(b1L,b2L,b3L,m),rhs(b1L,b2L,b3L,m))
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  IF (substep /= 1) THEN
     mult = dtime*bRK(substep)
     
     DO m = 1, dimens
        IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
        IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- advective terms ---------------------------------------------------------------------------------------
  IF (Stokes_yes) THEN
                      nl(S11B:N11B,S21B:N21B,S31B:N31B,1) = 0.
                      nl(S12B:N12B,S22B:N22B,S32B:N32B,2) = 0.
     IF (dimens == 3) nl(S13B:N13B,S23B:N23B,S33B:N33B,3) = 0.
  ELSE
     CALL nonlinear(.FALSE.)
  END IF
  
#ifdef NONBOUSSINESQ
  !--- multiply by total fluid density -----------------------------------------------------------------------
  IF (.NOT. Stokes_yes .AND. nonBoussinesq_yes) THEN
     DO k = S31B, N31B
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              nl(i,j,k,1) = nl(i,j,k,1)*dens(i,j,k,1)
           END DO
        END DO
     END DO
     DO k = S32B, N32B
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              nl(i,j,k,2) = nl(i,j,k,2)*dens(i,j,k,2)
           END DO
        END DO
     END DO
     IF (dimens == 3) THEN
     DO k = S33B, N33B
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              nl(i,j,k,3) = nl(i,j,k,3)*dens(i,j,k,3)
           END DO
        END DO
     END DO
     END IF
  END IF
  
  !--- transport of concentration diffusion ------------------------------------------------------------------
  IF (timeint_mode == 1 .AND. (.NOT. Euler_yes) .AND. nonBoussinesq_yes) THEN
     res = 0.
     DO conc_nu = 1, n_conc
        m = conc_nu
        multL = Fro*Ric(m)/(Re*Sc(m))
        CALL Helmholtz_conc_explicit(.FALSE.,res)
     END DO
     
     CALL interpolate2_pre_vel(.TRUE.,1,res,gpre)
     DO k = S31B, N31B
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              nl(i,j,k,1) = nl(i,j,k,1) - vel(i,j,k,1)*gpre(i,j,k)
           END DO
        END DO
     END DO
     
     CALL interpolate2_pre_vel(.TRUE.,2,res,gpre)
     DO k = S32B, N32B
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              nl(i,j,k,2) = nl(i,j,k,2) - vel(i,j,k,2)*gpre(i,j,k)
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     CALL interpolate2_pre_vel(.TRUE.,3,res,gpre)
     DO k = S33B, N33B
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              nl(i,j,k,3) = nl(i,j,k,3) - vel(i,j,k,3)*gpre(i,j,k)
           END DO
        END DO
     END DO
     END IF
  END IF
#endif
  
  !--- viscose terms -----------------------------------------------------------------------------------------
  IF (timeint_mode == 1 .AND. (.NOT. Euler_yes)) THEN
     multL = 1./Re
     CALL Helmholtz_explicit(.FALSE.)
  END IF
  
  !--- feedback-forces concentrations ------------------------------------------------------------------------
  IF (concentration_yes) CALL interpolate_conc(.FALSE.)
  
  !--- feedback-forces Lagrangian particles ------------------------------------------------------------------
  IF (particles_yes) CALL move_particles(.FALSE.) ! TEST!!! Eine dauerhaftere Loesung suchen, die auch die Kopplung mit conc elegant zulaesst! TEST!!! umbenennen ...
  
  !--- LES closure -------------------------------------------------------------------------------------------
  IF (LES_mode == 1) THEN
     CALL admrt(1,0,S11B,S21B,S31B,N11B,N21B,N31B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,1),nl(b1L,b2L,b3L,1))
     CALL admrt(2,0,S12B,S22B,S32B,N12B,N22B,N32B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,2),nl(b1L,b2L,b3L,2))
     CALL admrt(3,0,S13B,S23B,S33B,N13B,N23B,N33B,n_lp_vel,n_hp_vel,chi_vel,vel(b1L,b2L,b3L,3),nl(b1L,b2L,b3L,3))
  END IF
  IF (LES_mode == 2) THEN
     CALL smag_vel(.FALSE.)
  END IF
  
  !--- forcing -----------------------------------------------------------------------------------------------
  CALL forcing_vel
  
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  
  DO m = 1, dimens
     IF (m == 1) rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - mult*nl(S11B:N11B,S21B:N21B,S31B:N31B,1)
     IF (m == 2) rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - mult*nl(S12B:N12B,S22B:N22B,S32B:N32B,2)
     IF (m == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - mult*nl(S13B:N13B,S23B:N23B,S33B:N33B,3)
  END DO
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
  IF (substep /= 1) THEN
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl_old ---------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*bRK(substep)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
     IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
     
     IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
     IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
     
     IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
     IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
     IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
     
     IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
     IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
     
     IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
     IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
     !--------------------------------------------------------------------------------------------------------
     IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
     IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
     
     IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
     IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
     
     IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
     IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
     !--------------------------------------------------------------------------------------------------------
  END IF
  
  !--- Ausfluss-RB -------------------------------------------------------------------------------------------
  CALL outflow_bc
  
  !--- andere RB (instationaer + Zeitintegration) ------------------------------------------------------------
  CALL boundary_vel_tint
  
  !-----------------------------------------------------------------------------------------------------------
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - mult*nlbc11(S21B:N21B,S31B:N31B,1)
  IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - mult*nlbc11(S21B:N21B,S31B:N31B,2)
  
  IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - mult*nlbc21(S22B:N22B,S32B:N32B,1)
  IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - mult*nlbc21(S22B:N22B,S32B:N32B,2)
  
  IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - mult*nlbc31(S23B:N23B,S33B:N33B,1)
  IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - mult*nlbc31(S23B:N23B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - mult*nlbc12(S11B:N11B,S31B:N31B,1)
  IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - mult*nlbc12(S11B:N11B,S31B:N31B,2)
  
  IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - mult*nlbc22(S12B:N12B,S32B:N32B,1)
  IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - mult*nlbc22(S12B:N12B,S32B:N32B,2)
  
  IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - mult*nlbc32(S13B:N13B,S33B:N33B,1)
  IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - mult*nlbc32(S13B:N13B,S33B:N33B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - mult*nlbc13(S11B:N11B,S21B:N21B,1)
  IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - mult*nlbc13(S11B:N11B,S21B:N21B,2)
  
  IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - mult*nlbc23(S12B:N12B,S22B:N22B,1)
  IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - mult*nlbc23(S12B:N12B,S22B:N22B,2)
  
  IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - mult*nlbc33(S13B:N13B,S23B:N23B,1)
  IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - mult*nlbc33(S13B:N13B,S23B:N23B,2)
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen (instantan) ===========================================================================
  !===========================================================================================================
  CALL boundary_vel_stat
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== RB einsetzen ==========================================================================================
  !===========================================================================================================
  ! Anmerkung: Reihenfolge bestimmt Konvention (Linien, Kanten)!
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L > 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
  IF (BC_2U > 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
  
  IF (BC_3L > 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
  IF (BC_3U > 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
  
  IF (BC_1L > 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
  IF (BC_1U > 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L > 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
  IF (BC_1U > 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
  
  IF (BC_3L > 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
  IF (BC_3U > 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
  
  IF (BC_2L > 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
  IF (BC_2U > 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L > 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
  IF (BC_1U > 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
  
  IF (BC_2L > 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
  IF (BC_2U > 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
  
  IF (BC_3L > 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
  IF (BC_3U > 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== Fluss-Korrektur =======================================================================================
  !===========================================================================================================
  IF (nullspace_yes) THEN
     flux = 0.
     
     m = 1
     DO k = S31B, N31B
        DO j = S21B, N21B
!pgi$ unroll = n:8
           DO i = S11B, N11B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     
     m = 2
     DO k = S32B, N32B
        DO j = S22B, N22B
!pgi$ unroll = n:8
           DO i = S12B, N12B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     m = 3
     DO k = S33B, N33B
        DO j = S23B, N23B
!pgi$ unroll = n:8
           DO i = S13B, N13B
              flux = flux + psi_vel(i,j,k,m)*rhs(i,j,k,m)
           END DO
        END DO
     END DO
     END IF
     
     
     CALL MPI_ALLREDUCE(flux,flux_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     flux = flux_global
     
     IF (rank == 0) WRITE(10,'(a,E25.17)') 'flux =', flux
     
     
     ! d/dt ist unstetig in der Zeit ==> kein Versuch unternommen, in Zeitintegration einzubauen ...
     IF (nullspace_ortho_yes) THEN
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L > 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,1        ,S31B:N31B,1)
        IF (BC_2U > 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*psi_vel(S11B:N11B,N2       ,S31B:N31B,1)
        
        IF (BC_3L > 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,1        ,1)
        IF (BC_3U > 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*psi_vel(S11B:N11B,S21B:N21B,N3       ,1)
        
        IF (BC_1L > 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*psi_vel(0        ,S21B:N21B,S31B:N31B,1)
        IF (BC_1U > 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*psi_vel(N1       ,S21B:N21B,S31B:N31B,1)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*psi_vel(1        ,S22B:N22B,S32B:N32B,2)
        IF (BC_1U > 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*psi_vel(N1       ,S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L > 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*psi_vel(S12B:N12B,S22B:N22B,1        ,2)
        IF (BC_3U > 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,N3       ,2)
        
        IF (BC_2L > 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*psi_vel(S12B:N12B,0        ,S32B:N32B,2)
        IF (BC_2U > 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,N2       ,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*psi_vel(1        ,S23B:N23B,S33B:N33B,3)
        IF (BC_1U > 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*psi_vel(N1       ,S23B:N23B,S33B:N33B,3)
        
        IF (BC_2L > 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*psi_vel(S13B:N13B,1        ,S33B:N33B,3)
        IF (BC_2U > 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*psi_vel(S13B:N13B,N2       ,S33B:N33B,3)
        
        IF (BC_3L > 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*psi_vel(S13B:N13B,S23B:N23B,0        ,3)
        IF (BC_3U > 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*psi_vel(S13B:N13B,S23B:N23B,N3       ,3)
        !-----------------------------------------------------------------------------------------------------
        
                         rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) = rhs(S11B:N11B,S21B:N21B,S31B:N31B,1) - flux*psi_vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
                         rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) = rhs(S12B:N12B,S22B:N22B,S32B:N32B,2) - flux*psi_vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
        IF (dimens == 3) rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) = rhs(S13B:N13B,S23B:N23B,S33B:N33B,3) - flux*psi_vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
        
     ELSE
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L > 0) bc12(S11B:N11B,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1) - flux*th12(S11B:N11B,S31B:N31B,1)
        IF (BC_2U > 0) bc12(S11B:N11B,S31B:N31B,2) = bc12(S11B:N11B,S31B:N31B,2) - flux*th12(S11B:N11B,S31B:N31B,2)
        
        IF (BC_3L > 0) bc13(S11B:N11B,S21B:N21B,1) = bc13(S11B:N11B,S21B:N21B,1) - flux*th13(S11B:N11B,S21B:N21B,1)
        IF (BC_3U > 0) bc13(S11B:N11B,S21B:N21B,2) = bc13(S11B:N11B,S21B:N21B,2) - flux*th13(S11B:N11B,S21B:N21B,2)
        
        IF (BC_1L > 0) bc11(S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1) - flux*th11(S21B:N21B,S31B:N31B,1)
        IF (BC_1U > 0) bc11(S21B:N21B,S31B:N31B,2) = bc11(S21B:N21B,S31B:N31B,2) - flux*th11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) bc21(S22B:N22B,S32B:N32B,1) = bc21(S22B:N22B,S32B:N32B,1) - flux*th21(S22B:N22B,S32B:N32B,1)
        IF (BC_1U > 0) bc21(S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2) - flux*th21(S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L > 0) bc23(S12B:N12B,S22B:N22B,1) = bc23(S12B:N12B,S22B:N22B,1) - flux*th23(S12B:N12B,S22B:N22B,1)
        IF (BC_3U > 0) bc23(S12B:N12B,S22B:N22B,2) = bc23(S12B:N12B,S22B:N22B,2) - flux*th23(S12B:N12B,S22B:N22B,2)
        
        IF (BC_2L > 0) bc22(S12B:N12B,S32B:N32B,1) = bc22(S12B:N12B,S32B:N32B,1) - flux*th22(S12B:N12B,S32B:N32B,1)
        IF (BC_2U > 0) bc22(S12B:N12B,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2) - flux*th22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) bc31(S23B:N23B,S33B:N33B,1) = bc31(S23B:N23B,S33B:N33B,1) - flux*th31(S23B:N23B,S33B:N33B,1)
        IF (BC_1U > 0) bc31(S23B:N23B,S33B:N33B,2) = bc31(S23B:N23B,S33B:N33B,2) - flux*th31(S23B:N23B,S33B:N33B,2)
        
        IF (BC_2L > 0) bc32(S13B:N13B,S33B:N33B,1) = bc32(S13B:N13B,S33B:N33B,1) - flux*th32(S13B:N13B,S33B:N33B,1)
        IF (BC_2U > 0) bc32(S13B:N13B,S33B:N33B,2) = bc32(S13B:N13B,S33B:N33B,2) - flux*th32(S13B:N13B,S33B:N33B,2)
        
        IF (BC_3L > 0) bc33(S13B:N13B,S23B:N23B,1) = bc33(S13B:N13B,S23B:N23B,1) - flux*th33(S13B:N13B,S23B:N23B,1)
        IF (BC_3U > 0) bc33(S13B:N13B,S23B:N23B,2) = bc33(S13B:N13B,S23B:N23B,2) - flux*th33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
        
        !-----------------------------------------------------------------------------------------------------
        IF (BC_2L > 0) rhs(S11B:N11B,1        ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,1)
        IF (BC_2U > 0) rhs(S11B:N11B,N2       ,S31B:N31B,1) = bc12(S11B:N11B,S31B:N31B,2)
        
        IF (BC_3L > 0) rhs(S11B:N11B,S21B:N21B,1        ,1) = bc13(S11B:N11B,S21B:N21B,1)
        IF (BC_3U > 0) rhs(S11B:N11B,S21B:N21B,N3       ,1) = bc13(S11B:N11B,S21B:N21B,2)
        
        IF (BC_1L > 0) rhs(0        ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,1)
        IF (BC_1U > 0) rhs(N1       ,S21B:N21B,S31B:N31B,1) = bc11(S21B:N21B,S31B:N31B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) rhs(1        ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,1)
        IF (BC_1U > 0) rhs(N1       ,S22B:N22B,S32B:N32B,2) = bc21(S22B:N22B,S32B:N32B,2)
        
        IF (BC_3L > 0) rhs(S12B:N12B,S22B:N22B,1        ,2) = bc23(S12B:N12B,S22B:N22B,1)
        IF (BC_3U > 0) rhs(S12B:N12B,S22B:N22B,N3       ,2) = bc23(S12B:N12B,S22B:N22B,2)
        
        IF (BC_2L > 0) rhs(S12B:N12B,0        ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,1)
        IF (BC_2U > 0) rhs(S12B:N12B,N2       ,S32B:N32B,2) = bc22(S12B:N12B,S32B:N32B,2)
        !-----------------------------------------------------------------------------------------------------
        IF (BC_1L > 0) rhs(1        ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,1)
        IF (BC_1U > 0) rhs(N1       ,S23B:N23B,S33B:N33B,3) = bc31(S23B:N23B,S33B:N33B,2)
        
        IF (BC_2L > 0) rhs(S13B:N13B,1        ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,1)
        IF (BC_2U > 0) rhs(S13B:N13B,N2       ,S33B:N33B,3) = bc32(S13B:N13B,S33B:N33B,2)
        
        IF (BC_3L > 0) rhs(S13B:N13B,S23B:N23B,0        ,3) = bc33(S13B:N13B,S23B:N23B,1)
        IF (BC_3U > 0) rhs(S13B:N13B,S23B:N23B,N3       ,3) = bc33(S13B:N13B,S23B:N23B,2)
        !-----------------------------------------------------------------------------------------------------
        
     END IF
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE rhs_vel
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE rhs_conc
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  REAL                   ::  mult
  
  
  m = conc_nu
  
  
  !===========================================================================================================
  !=== rhs = rhs + Helmholtz(conc) ===========================================================================
  !===========================================================================================================
  IF (timeint_mode == 1 .OR. thetaL == 1.) THEN
     res(S1p:N1p,S2p:N2p,S3p:N3p) = conc(S1p:N1p,S2p:N2p,S3p:N3p,m)
  ELSE
     multL = -(1.-thetaL)*(aRK(substep)+bRK(substep))*dtime / (Re*Sc(m))
     CALL Helmholtz_conc(.FALSE.,conc(b1L,b2L,b3L,m),res)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl_old ====================================================================================
  !===========================================================================================================
  IF (substep /= 1) THEN
     mult = dtime*bRK(substep)
     res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - mult*nlco(S1p:N1p,S2p:N2p,S3p:N3p,m)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== rhs = rhs - nl ========================================================================================
  !===========================================================================================================
  !--- Konvektions-Term --------------------------------------------------------------------------------------
  ! TEST!!! Stokes flow?? Evtl. analog fluid auch hier einbauen ...
  CALL nonlinear_conc(.FALSE.)
  
  !--- diffusive term ----------------------------------------------------------------------------------------
  IF (timeint_mode == 1 .AND. (.NOT. Euler_yes)) THEN
     multL = 1./(Re*Sc(m))
     CALL Helmholtz_conc_explicit(.FALSE.,nlco(b1L,b2L,b3L,m))
  END IF
  
  !--- LES closure -------------------------------------------------------------------------------------------
  IF (LES_mode == 1) CALL admrt(0,m,S1p,S2p,S3p,N1p,N2p,N3p,n_lp_conc(m),n_hp_conc(m),chi_conc(m),conc(b1L,b2L,b3L,m),nlco(b1L,b2L,b3L,m))
  IF (LES_mode == 2) CALL smag_conc(.FALSE.)
  
  !--- forcing -----------------------------------------------------------------------------------------------
  CALL forcing_conc
  
  !--- rhs = rhs - nl ----------------------------------------------------------------------------------------
  mult = dtime*aRK(substep)
  res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) - mult*nlco(S1p:N1p,S2p:N2p,S3p:N3p,m)
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== Randbedingungen (Zeitintegration) =====================================================================
  !===========================================================================================================
  ! - sed1L & Co. k�nnten auch in nlco gespeichert werden. So ist es allerdings sicherer!
  ! - evtl. bc-Felder einf�hren und analog zu den Geschwindigkeiten verfahren ...
  
  IF (substep == 1) THEN
     !--- Ausfluss-RB ----------------------------------------------------------------------------------------
     CALL sediment_bc
     
     !--- andere RB (instationaer + Zeitintegration) ---------------------------------------------------------
     CALL boundary_conc_tint
     
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl -------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*aRK(substep)
     
     IF (BCc_1L(m) == 1 .AND. (.NOT. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = conc(1 ,S2p:N2p,S3p:N3p,m) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     IF (BCc_1U(m) == 1 .AND. (.NOT. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = conc(N1,S2p:N2p,S3p:N3p,m) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     IF (BCc_2L(m) == 1 .AND. (.NOT. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = conc(S1p:N1p,1 ,S3p:N3p,m) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     IF (BCc_2U(m) == 1 .AND. (.NOT. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = conc(S1p:N1p,N2,S3p:N3p,m) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     IF (BCc_3L(m) == 1 .AND. (.NOT. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = conc(S1p:N1p,S2p:N2p,1 ,m) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     IF (BCc_3U(m) == 1 .AND. (.NOT. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = conc(S1p:N1p,S2p:N2p,N3,m) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
  ELSE
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl_old ---------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*bRK(substep)
     
     IF (BCc_1L(m) == 1 .AND. (.NOT. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = conc(1 ,S2p:N2p,S3p:N3p,m) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     IF (BCc_1U(m) == 1 .AND. (.NOT. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = conc(N1,S2p:N2p,S3p:N3p,m) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     IF (BCc_2L(m) == 1 .AND. (.NOT. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = conc(S1p:N1p,1 ,S3p:N3p,m) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     IF (BCc_2U(m) == 1 .AND. (.NOT. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = conc(S1p:N1p,N2,S3p:N3p,m) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     IF (BCc_3L(m) == 1 .AND. (.NOT. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = conc(S1p:N1p,S2p:N2p,1 ,m) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     IF (BCc_3U(m) == 1 .AND. (.NOT. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = conc(S1p:N1p,S2p:N2p,N3,m) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
     
     !--- Ausfluss-RB ----------------------------------------------------------------------------------------
     CALL sediment_bc
     
     !--- andere RB (instationaer + Zeitintegration) ---------------------------------------------------------
     CALL boundary_conc_tint
     
     !--------------------------------------------------------------------------------------------------------
     !--- rhs = rhs - nl -------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mult = dtime*aRK(substep)
     
     IF (BCc_1L(m) == 1 .AND. (.NOT. isopycnal(1,1,m))) res(1 ,S2p:N2p,S3p:N3p) = res(1 ,S2p:N2p,S3p:N3p) - mult*sed1L(S2p:N2p,S3p:N3p,m)
     IF (BCc_1U(m) == 1 .AND. (.NOT. isopycnal(1,2,m))) res(N1,S2p:N2p,S3p:N3p) = res(N1,S2p:N2p,S3p:N3p) - mult*sed1U(S2p:N2p,S3p:N3p,m)
     
     IF (BCc_2L(m) == 1 .AND. (.NOT. isopycnal(2,1,m))) res(S1p:N1p,1 ,S3p:N3p) = res(S1p:N1p,1 ,S3p:N3p) - mult*sed2L(S1p:N1p,S3p:N3p,m)
     IF (BCc_2U(m) == 1 .AND. (.NOT. isopycnal(2,2,m))) res(S1p:N1p,N2,S3p:N3p) = res(S1p:N1p,N2,S3p:N3p) - mult*sed2U(S1p:N1p,S3p:N3p,m)
     
     IF (BCc_3L(m) == 1 .AND. (.NOT. isopycnal(3,1,m))) res(S1p:N1p,S2p:N2p,1 ) = res(S1p:N1p,S2p:N2p,1 ) - mult*sed3L(S1p:N1p,S2p:N2p,m)
     IF (BCc_3U(m) == 1 .AND. (.NOT. isopycnal(3,2,m))) res(S1p:N1p,S2p:N2p,N3) = res(S1p:N1p,S2p:N2p,N3) - mult*sed3U(S1p:N1p,S2p:N2p,m)
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Randbedingungen (instantan) ===========================================================================
  !===========================================================================================================
  CALL boundary_conc_stat
  !===========================================================================================================
  
  
  END SUBROUTINE rhs_conc
  
  
  
  
END MODULE mod_rhs
