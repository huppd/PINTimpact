!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module containing the important timeintegration
MODULE mod_timeint
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE mod_inout
  USE mod_coeffs
  USE mod_lib
  USE mod_test
  USE mod_particles
  USE mod_solvers
  USE mod_rhs
  
  PRIVATE
  
  PUBLIC timeintegration
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !> \brief routine for the timeintegration
  !!
  !!  start up:
  !!    - Startvektor
  !!    - Nullräume bestimmen
  !!    - RB initialisieren
  !!    - Divergenzfreiheit bestimmen
  !!    - diverse Files öffnen
  !!    - File fuer zwischenzeitlichen Abbruch der Zeitintegration neu erstellen
  !!    - Zeitmessung starten
  !!    - Ausschreiben
  !!
  !!  Zeitintegration:
  !!    - compute timestep
  !!    - ghost cell update (fuer RHS)
  !!    - interpolate advection velocity + update ghost cells
  !!    - interpolate density on velocity grids
  !!    - rhs (ggf. Neumann-RB überschreiben)
  !!    - Helmholtz-Multiplikator
  !!    - compute right hand side
  !!    - Umskalieren (Effizienz, initial guess)
  !!    - Löser( timeint=1 and thetaL=1 -> ecplicit else outer_iteration/twostep)
  SUBROUTINE timeintegration
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  ! --- needed for stopping criterion -----------------------------------------------------------------------
  REAL                   ::  norm_velp
  REAL                   ::  norm_vel
  REAL                   ::  norm_dummy
  REAL                   ::  periodic
  REAL                   ::  periodic0

  
  !--- Startvektor -------------------------------------------------------------------------------------------
  ! Hier bereits notwendig, um weitere Konzentrationen zuschalten zu koennen (BC werden aber weiter unten initialisiert!).
                         CALL initial_conditions_vel
  IF (concentration_yes) CALL initial_conditions_conc
  IF (particles_yes    ) CALL initial_conditions_part
  
  
  IF (restart == 0) THEN
     time          = time_start
     time_out_vect = time_start
     time_out_scal = time_start
     
     dtime         = 0.
     timestep      = 0
     
     new_dtime      = .TRUE.
     write_out_vect = .TRUE.
     write_out_scal = .TRUE.
     
     write_count = 0
     
     IF (dtime_out_vect == 0.) write_out_vect = .FALSE.
     IF (dtime_out_scal == 0.) write_out_scal = .FALSE.
  ELSE
                               CALL read_restart
     IF (dtime_out_scal /= 0.) CALL read_restart_stats
     
     IF (rank == 0 .AND. write_stout_yes) THEN
        WRITE(*,'(a,1E13.5)') '             time =', time
        WRITE(*,'(a,1i5   )') '         timestep =', timestep
        
        WRITE(*,'(a,1E13.5)') '            dtime =', dtime
        WRITE(*,'(a,1E13.5)') '    time_out_vect =', time_out_vect
        WRITE(*,'(a,1E13.5)') '    time_out_scal =', time_out_scal
        
        WRITE(*,'(a,1L2   )') '        new_dtime =', new_dtime
        WRITE(*,'(a,1L2   )') '   write_out_vect =', write_out_vect
        WRITE(*,'(a,1L2   )') '   write_out_scal =', write_out_scal
        
        WRITE(*,'(a,1i5   )') '      write_count =', write_count
     END IF
  END IF
  
  timestep_old  = timestep
  dtime_old     = dtime
  dtime_average = 0.
  finish_yes    = .FALSE.
  
  !--- Null-Raeume bestimmen ---------------------------------------------------------------------------------
  ! Steht hier, weil Korrekturvektor "th" nach "configuration" erst alloziert und bestimmt werden muss
  ! ("initial_conditions_" werden danach als nächstes gerufen, s.o.)
  ! Alternative: eigene Subroutine für "th" kreieren ...
  IF (nullspace_yes) THEN
     CALL get_stencil_transp
     CALL get_nullspace
     CALL get_stencil ! TEST!!! Unschoen! besser zu impact.f90 packen und Reihenfolge abaendern ...
  END IF
  
  !--- RB initialisieren -------------------------------------------------------------------------------------
  CALL init_BC
  
  !--- Divergenz-Freiheit testen -----------------------------------------------------------------------------
  CALL test_divergence
  
  !--- diverse Files öffnen ----------------------------------------------------------------------------------
  CALL open_stats
  
  !--- File fuer zwischenzeitlichen Abbruch der Zeitintegration neu erstellen --------------------------------
  IF (rank == 0) THEN
     OPEN (10,FILE='send_signal.txt',STATUS='UNKNOWN')
     WRITE(10,'(a)') '0'
     WRITE(10,*) n_timesteps
     WRITE(10,*) time_end
     CALL flush(10)
     CLOSE(10)
  END IF
  
  
  IF (rank == 0 .AND. write_stout_yes) THEN
     WRITE(*,'(a)')
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
     WRITE(*,'(a)') '---------------------------- START TIME-INTEGRATION -----------------------------'
     WRITE(*,'(a)')
     WRITE(*,'(a,E12.5)')                 '                     Re =',Re
     WRITE(*,'(a,i4,a,i4,a,i4)')          '          box resolution:',M1,' x',M2,' x',M3
     WRITE(*,'(a,E12.5,a,E12.5,a,E12.5)') '          box dimension :',L1,' x',L2,' x',L3
     WRITE(*,'(a,E12.5)')                 '                   epsU =',epsU
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)') '================================================================================='
     WRITE(*,'(a)')
  END IF
  
  
  !--- Zeitmessung starten -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     day  = ctime(3)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     elatime = msec+1000*(sec+60*(minu+60*hour))
     OPEN(99,FILE='test_wallclocktime_restart'//restart_char//'.txt',STATUS='UNKNOWN')
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Begin time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &      ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     CALL flush(99)
  END IF
  
  
  !--- Ausschreiben ------------------------------------------------------------------------------------------
  IF (write_out_scal) CALL compute_stats
  IF (write_out_vect) CALL write_fields
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  velp = vel
  IF(rank == 0) OPEN(42,FILE='log_periodicity.txt')
  IF(rank == 0) WRITE(42,*) '  time step    ','time                         ','period                       ','periodicity'
  timeint: DO

     CALL get_dtime

      ! ------- stopping criterion: checking if period state is reached --------------------------------------
      IF( abs(time-tp - 1.0/freq) < dtime ) THEN
        dtime = time-tp - 1.0/freq
        IF(dtime < 0 ) THEN
          IF(rank == 0) WRITE(42,*) 'negative periodicity'
          CALL get_dtime
          EXIT timeint
        END IF
        IF(dtime ==0 ) THEN
          new_dtime = .TRUE.
          CALL get_dtime
          new_dtime = .FALSE.
        END IF
        velp = velp-vel
        CALL get_norms_vel(velp, .FALSE., .TRUE. ,norm_dummy, norm_velp)
!        call get_norms_vel(vel , .false., .true. ,norm_dummy, norm_vel )
        periodic = norm_velp
!        IF( time <= 1.0/freq ) THEN
!          periodic0 = periodic
!        END IF
!        periodic = periodic / periodic0


!        if(rank == 0) WRITE(42,'(a,i8,a,E25.17,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; period =',time*freq,'; periodicity =',periodic
!        if(rank == 0) WRITE(42,'(a,*,a,*,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; period =',time*freq,'; periodicity =',periodic
        IF(rank == 0) WRITE(42,*) timestep,' ',time,'  ',time*freq,'  ',periodic

        IF( periodic < periodic_tol ) THEN
!          dtime_out_vect = dtime/2.
!          time_end = time + 1./freq + 1.E-11
!          periodic_tol = periodic_tol/100.
          EXIT timeint
!          tp = time
!          velp = vel
        ELSE
          tp = time
          velp = vel
        END IF

      END IF

!      IF( time > 0.125/freq ) EXIT timeint
      ! ------- stopping criterion: end ----------------------------------------------------------------------
     
     
     IF (rank == 0 .AND. write_stout_yes) THEN
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)')
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
        WRITE(*,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,' ; time =',time,' ; dtime =',dtime
        IF (timeint_mode == 0) THEN
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a)') '================================================================================='
        END IF
     END IF
     
     IF (rank == 0 .AND. log_iteration_yes) THEN
        OPEN(10,FILE='log_iterations.txt', STATUS='UNKNOWN')
        WRITE(10,'(a)')
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; dtime =',dtime
        WRITE(10,'(a)') '================================================================================='
        WRITE(10,'(a)') '================================================================================='
     END IF
     
     !========================================================================================================
     
     RK_loop: DO substep = 1, RK_steps
        
        IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) THEN
           WRITE(*,'(a)') 
           WRITE(*,'(a)') '================================================================================='
           WRITE(*,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(*,'(a)') '================================================================================='
        END IF
        IF (rank == 0 .AND. log_iteration_yes) THEN
           WRITE(10,'(a)') 
           WRITE(10,'(a)') '================================================================================='
           WRITE(10,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           WRITE(10,'(a)') '================================================================================='
        END IF
        
        
        !--- time --------------------------------------------------------------------------------------------
        IF (substep == 1) subtime = time + dtime* aRK(1)
        IF (substep == 2) subtime = time + dtime*(aRK(1)+aRK(2)+bRK(2))
        IF (substep == 3) subtime = time + dtime
        
        
        !--- ghost cell update (fuer RHS) --------------------------------------------------------------------
        CALL exchange_all_all(.TRUE.,vel)
        
        IF (concentration_yes) THEN
           DO m = 1, n_conc
              CALL exchange(1,0,conc(b1L,b2L,b3L,m))
              CALL exchange(2,0,conc(b1L,b2L,b3L,m))
              CALL exchange(3,0,conc(b1L,b2L,b3L,m))
           END DO
        END IF
        
        
        !--- interpolate advection velocity + update ghost cells ---------------------------------------------
        ! vel(:,:,:,i) --> worki(:,:,:)
        CALL interpolate_vel(.FALSE.) ! TEST!!! Wurde teilweise schon bei Zeitschritt-Bestimmung erledigt!
        
        
#ifdef NONBOUSSINESQ
        !--- interpolate density on velocity grids -----------------------------------------------------------
        IF (nonBoussinesq_yes) THEN
           res(S1p:N1p,S2p:N2p,S3p:N3p) = 1.
           DO m = 1, n_conc
              res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) + Fro*Ric(m)*conc(S1p:N1p,S2p:N2p,S3p:N3p,m)
           END DO
           
                            CALL interpolate2_pre_vel(.TRUE.,1,res,dens(b1L,b2L,b3L,1))
                            CALL interpolate2_pre_vel(.TRUE.,2,res,dens(b1L,b2L,b3L,2))
           IF (dimens == 3) CALL interpolate2_pre_vel(.TRUE.,3,res,dens(b1L,b2L,b3L,3))
        END IF
#endif
        
        !--- rhs (ggf. Neumann-RB überschreiben) -------------------------------------------------------------
        ! Muss vor Konzentrationen kommen, weil
        !  - bcii für die no-flux-RB verwendet werden,
        !  - die Konzentrationen die Eddy-Viscosity des Geschwindigkeitsfeldes benoetigen.
        CALL rhs_vel
        
        
        !--- Neue Partikel hinzufuegen -----------------------------------------------------------------------
        IF (particles_yes .AND. substep == RK_steps) CALL add_particles
        
        
        !--- Konzentrationsfelder lösen ----------------------------------------------------------------------
        ! Anmerkung: Um rhs_conc 1-dimensional wählen zu können, muss dieser Abschnitt zwischen rhs_NS und
        !            Solver stehen! (evtl. ist thetaL == 0. ...)
        IF (concentration_yes) THEN
           DO conc_nu = 1, n_conc
              
              CALL rhs_conc
              
              multL = thetaL*(aRK(substep)+bRK(substep))*dtime / (Re*Sc(conc_nu))
              
              IF (timeint_mode == 1) THEN
                 CALL solve_conc_explicit
              ELSE
                 CALL solve_conc(.FALSE.,.FALSE.)
              END IF
              
           END DO
        END IF
        
        
        !--- Helmholtz-Multiplikator -------------------------------------------------------------------------
        multL = thetaL*(aRK(substep)+bRK(substep))*dtime / Re
        
        
        !--- Umskalieren (Effizienz, initial guess) ----------------------------------------------------------
        IF (.NOT. init_pre(substep)) pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) * (aRK(substep)+bRK(substep)) * dtime
        
        
        !--- Löser -------------------------------------------------------------------------------------------
!        IF (timeint_mode == 1 .OR. thetaL == 1.) THEN !< thetaL == 1 means implicite ?
        IF (timeint_mode == 1 .OR. thetaL == 0.) THEN !< thetaL == 1 means implicite ? danielized doesnot really mather because when thetaL=0 timeintmod=1 is automatically switched on
#ifdef NONBOUSSINESQ
           IF (nonBoussinesq_yes) THEN
              CALL non_Boussinesq
           ELSE
              CALL explicit
           END IF
#else
           CALL explicit
#endif
        ELSE
           IF (twostep_yes) THEN
              CALL twostep
           ELSE
              CALL outer_iteration
           END IF
        END IF
        
        !--- physikalischer Druck ----------------------------------------------------------------------------
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) / (aRK(substep)+bRK(substep)) / dtime
        
        !--- Undefinierte Ecken / Kanten auffüllen -----------------------------------------------------------
        CALL fill_corners(pre)
        IF (corner_yes) THEN
           DO m = 1, n_conc
              CALL fill_corners(conc(b1L,b2L,b3L,m))
           END DO
        END IF
        
     END DO RK_loop
     
     !========================================================================================================
     timestep = timestep + 1
     time     = time + dtime
     
     !--- send_signal.txt lesen ------------------------------------------------------------------------------
     CALL check_signal
     
     
     !--- Druck-Niveau festhalten ----------------------------------------------------------------------------
     CALL level_pressure
     
     
     !--- Ausschreiben ---------------------------------------------------------------------------------------
     IF (write_out_scal) CALL compute_stats
     IF (write_out_vect) CALL write_fields
     
     
     !--------------------------------------------------------------------------------------------------------
     IF (rank == 0 .AND. log_iteration_yes) CLOSE(10)
     
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     IF (rank == 0 .AND. write_stout_yes .AND. timeint_mode == 0) WRITE(*,*)
     
     CALL MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror) ! notwendig fuer "check_alarm"
     IF (finish_yes) EXIT
     
  END DO timeint
  !===========================================================================================================
  
  
  !--- Zeitmessung beenden -----------------------------------------------------------------------------------
  IF (rank == 0) THEN
     CALL DATE_AND_TIME(values=ctime)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     IF (ctime(3) /= day) THEN
        ! Anmerkung: Gilt nur für Jobs <= 24h
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime + 24*60*60*1000
     ELSE
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime
     END IF
     
     WRITE(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Finish time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &       ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     WRITE(99,'(a,E13.5)') 'elapsed time [sec]', REAL(elatime)/1000.
     CLOSE(99)
  END IF
  
  !--- Restart schreiben -------------------------------------------------------------------------------------
  restart = restart + 1
  CALL write_restart
  CALL write_restart_stats
  
  
  !--- Iterationsstatistiken auswerten -----------------------------------------------------------------------
  CALL iteration_stats
  
  !--- diverse Files schliessen ------------------------------------------------------------------------------
  CALL close_stats
  
  CLOSE(42)

  
  END SUBROUTINE timeintegration
  
  
  
END MODULE mod_timeint
