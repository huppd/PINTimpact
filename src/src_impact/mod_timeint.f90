!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module containing the important timeintegration
module mod_timeint
  
  use iso_c_binding
  
  use mod_dims
  use mod_vars
  use mod_exchange
  use mod_diff
  use mod_inout
  use mod_coeffs
  use mod_lib
  use mod_test
  use mod_particles
  use mod_solvers
  use mod_rhs
  
  private
  
  public timeintegration
  
  INCLUDE 'mpif.h'
  
  contains
  
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
  subroutine timeintegration() bind(c,name='ftimeintegration')
  
  implicit none
  
  integer                ::  m
  
  ! --- needed for stopping criterion -----------------------------------------------------------------------
  real                   ::  norm_velp
  real                   ::  norm_vel
  real                   ::  norm_dummy
  real                   ::  periodic
  real                   ::  periodic0

  
  !--- Startvektor -------------------------------------------------------------------------------------------
  ! Hier bereits notwendig, um weitere Konzentrationen zuschalten zu koennen (BC werden aber weiter unten initialisiert!).
                         call initial_conditions_vel
  if (concentration_yes) call initial_conditions_conc
  if (particles_yes    ) call initial_conditions_part
  
  
  if (restart == 0) then
     time          = time_start
     time_out_vect = time_start
     time_out_scal = time_start
     
     dtime         = 0.
     timestep      = 0
     
     new_dtime      = .true.
     write_out_vect = .true.
     write_out_scal = .true.
     
     write_count = 0
     
     if (dtime_out_vect == 0.) write_out_vect = .false.
     if (dtime_out_scal == 0.) write_out_scal = .false.
  else
                               call read_restart
     if (dtime_out_scal /= 0.) call read_restart_stats
     
     if (rank == 0 .and. write_stout_yes) then
        write(*,'(a,1E13.5)') '             time =', time
        write(*,'(a,1i5   )') '         timestep =', timestep
        
        write(*,'(a,1E13.5)') '            dtime =', dtime
        write(*,'(a,1E13.5)') '    time_out_vect =', time_out_vect
        write(*,'(a,1E13.5)') '    time_out_scal =', time_out_scal
        
        write(*,'(a,1L2   )') '        new_dtime =', new_dtime
        write(*,'(a,1L2   )') '   write_out_vect =', write_out_vect
        write(*,'(a,1L2   )') '   write_out_scal =', write_out_scal
        
        write(*,'(a,1i5   )') '      write_count =', write_count
     end if
  end if
  
  timestep_old  = timestep
  dtime_old     = dtime
  dtime_average = 0.
  finish_yes    = .false.
  
  !--- Null-Raeume bestimmen ---------------------------------------------------------------------------------
  ! Steht hier, weil Korrekturvektor "th" nach "configuration" erst alloziert und bestimmt werden muss
  ! ("initial_conditions_" werden danach als nächstes gerufen, s.o.)
  ! Alternative: eigene Subroutine für "th" kreieren ...
  if (nullspace_yes) then
     call get_stencil_transp
     call get_nullspace
     call get_stencil ! TEST!!! Unschoen! besser zu impact.f90 packen und Reihenfolge abaendern ...
  end if
  
  !--- RB initialisieren -------------------------------------------------------------------------------------
  call init_BC
  
  !--- Divergenz-Freiheit testen -----------------------------------------------------------------------------
  call test_divergence
  
  !--- diverse Files öffnen ----------------------------------------------------------------------------------
  call open_stats
  
  !--- File fuer zwischenzeitlichen Abbruch der Zeitintegration neu erstellen --------------------------------
  if (rank == 0) then
     open (10,file='send_signal.txt',status='UNKNOWN')
     write(10,'(a)') '0'
     write(10,*) n_timesteps
     write(10,*) time_end
     call flush(10)
     close(10)
  end if
  
  
  if (rank == 0 .and. write_stout_yes) then
     write(*,'(a)')
     write(*,'(a)') '================================================================================='
     write(*,'(a)') '================================================================================='
     write(*,'(a)') '================================================================================='
     write(*,'(a)')
     write(*,'(a)') '---------------------------- START TIME-INTEGRATION -----------------------------'
     write(*,'(a)')
     write(*,'(a,E12.5)')                 '                     Re =',Re
     write(*,'(a,i4,a,i4,a,i4)')          '          box resolution:',M1,' x',M2,' x',M3
     write(*,'(a,E12.5,a,E12.5,a,E12.5)') '          box dimension :',L1,' x',L2,' x',L3
     write(*,'(a,E12.5)')                 '                   epsU =',epsU
     write(*,'(a)') '================================================================================='
     write(*,'(a)') '================================================================================='
     write(*,'(a)') '================================================================================='
     write(*,'(a)')
  end if
  
  
  !--- Zeitmessung starten -----------------------------------------------------------------------------------
  if (rank == 0) then
     call DATE_AND_TIME(values=ctime)
     day  = ctime(3)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     elatime = msec+1000*(sec+60*(minu+60*hour))
     open(99,file='test_wallclocktime_restart'//restart_char//'.txt',status='UNKNOWN')
     write(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Begin time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &      ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     call flush(99)
  end if
  
  
  !--- Ausschreiben ------------------------------------------------------------------------------------------
  if (write_out_scal) call compute_stats
  if (write_out_vect) call write_fields
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  velp = vel
  if(rank == 0) open(42,file='log_periodicity.txt')
  if(rank == 0) write(42,*) '  time step    ','time                         ','period                       ','periodicity'
  timeint: do

     call get_dtime
      ! ------- stopping criterion: checking if period state is reached --------------------------------------

      if( abs(time-tp - 1.0/freq) < dtime ) then
        dtime = time-tp - 1.0/freq
        if(dtime < 0 ) then
          if(rank == 0) write(42,*) 'negative periodicity'
          call get_dtime
          exit timeint
        end if
        if(dtime ==0 ) then
          new_dtime = .true.
          call get_dtime
          new_dtime = .false.
        end if
        velp = velp-vel
        call get_norms_vel(velp, .false., .true. ,norm_dummy, norm_velp)
!        call get_norms_vel(vel , .false., .true. ,norm_dummy, norm_vel )
        periodic = norm_velp
!        IF( time <= 1.0/freq ) THEN
!          periodic0 = periodic
!        END IF
!        periodic = periodic / periodic0


!        if(rank == 0) WRITE(42,'(a,i8,a,E25.17,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; period =',time*freq,'; periodicity =',periodic
!        if(rank == 0) WRITE(42,'(a,*,a,*,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; period =',time*freq,'; periodicity =',periodic
        if(rank == 0) write(42,*) timestep,' ',time,'  ',time*freq,'  ',periodic

        if( periodic < periodic_tol ) then
!          dtime_out_vect = dtime/2.
!          time_end = time + 1./freq + 1.E-11
!          periodic_tol = periodic_tol/100.
          exit timeint
!          tp = time
!          velp = vel
        else
          tp = time
          velp = vel
        end if

      end if

!      IF( time > 0.125/freq ) EXIT timeint
      ! ------- stopping criterion: end ----------------------------------------------------------------------
     
     
     if (rank == 0 .and. write_stout_yes) then
        if (timeint_mode == 0) then
           write(*,'(a)')
           write(*,'(a)') '================================================================================='
           write(*,'(a)') '================================================================================='
        end if
        write(*,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,' ; time =',time,' ; dtime =',dtime
        if (timeint_mode == 0) then
           write(*,'(a)') '================================================================================='
           write(*,'(a)') '================================================================================='
        end if
     end if
     
     if (rank == 0 .and. log_iteration_yes) then
        open(10,file='log_iterations.txt', status='UNKNOWN')
        write(10,'(a)')
        write(10,'(a)') '================================================================================='
        write(10,'(a)') '================================================================================='
        write(10,'(a,i8,a,E25.17,a,E25.17)') 'time step = ',timestep,'; time =',time,'; dtime =',dtime
        write(10,'(a)') '================================================================================='
        write(10,'(a)') '================================================================================='
     end if
     
     !========================================================================================================
     
     RK_loop: do substep = 1, RK_steps
        
        if (rank == 0 .and. write_stout_yes .and. timeint_mode == 0) then
           write(*,'(a)') 
           write(*,'(a)') '================================================================================='
           write(*,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           write(*,'(a)') '================================================================================='
        end if
        if (rank == 0 .and. log_iteration_yes) then
           write(10,'(a)') 
           write(10,'(a)') '================================================================================='
           write(10,'(a,i2,a)') 'Runge-Kutta sub-step',substep,':'
           write(10,'(a)') '================================================================================='
        end if
        
        
        !--- Zeit --------------------------------------------------------------------------------------------
        if (substep == 1) subtime = time + dtime* aRK(1)
        if (substep == 2) subtime = time + dtime*(aRK(1)+aRK(2)+bRK(2))
        if (substep == 3) subtime = time + dtime
        
        
        !--- ghost cell update (fuer RHS) --------------------------------------------------------------------
        call exchange_all_all(.true.,vel)
        
        if (concentration_yes) then
           do m = 1, n_conc
              call exchange(1,0,conc(b1L,b2L,b3L,m))
              call exchange(2,0,conc(b1L,b2L,b3L,m))
              call exchange(3,0,conc(b1L,b2L,b3L,m))
           end do
        end if
        
        
        !--- interpolate advection velocity + update ghost cells ---------------------------------------------
        ! vel(:,:,:,i) --> worki(:,:,:)
        call interpolate_vel(.false.) ! TEST!!! Wurde teilweise schon bei Zeitschritt-Bestimmung erledigt!
        
        
#ifdef NONBOUSSINESQ
        !--- interpolate density on velocity grids -----------------------------------------------------------
        if (nonBoussinesq_yes) then
           res(S1p:N1p,S2p:N2p,S3p:N3p) = 1.
           do m = 1, n_conc
              res(S1p:N1p,S2p:N2p,S3p:N3p) = res(S1p:N1p,S2p:N2p,S3p:N3p) + Fro*Ric(m)*conc(S1p:N1p,S2p:N2p,S3p:N3p,m)
           end do
           
                            call interpolate2_pre_vel(.true.,1,res,dens(b1L,b2L,b3L,1))
                            call interpolate2_pre_vel(.true.,2,res,dens(b1L,b2L,b3L,2))
           if (dimens == 3) call interpolate2_pre_vel(.true.,3,res,dens(b1L,b2L,b3L,3))
        end if
#endif
        
        !--- rhs (ggf. Neumann-RB überschreiben) -------------------------------------------------------------
        ! Muss vor Konzentrationen kommen, weil
        !  - bcii für die no-flux-RB verwendet werden,
        !  - die Konzentrationen die Eddy-Viscosity des Geschwindigkeitsfeldes benoetigen.
        call rhs_vel
        
        
        !--- Neue Partikel hinzufuegen -----------------------------------------------------------------------
        if (particles_yes .and. substep == RK_steps) call add_particles
        
        
        !--- Konzentrationsfelder lösen ----------------------------------------------------------------------
        ! Anmerkung: Um rhs_conc 1-dimensional wählen zu können, muss dieser Abschnitt zwischen rhs_NS und
        !            Solver stehen! (evtl. ist thetaL == 0. ...)
        if (concentration_yes) then
           do conc_nu = 1, n_conc
              
              call rhs_conc
              
              multL = thetaL*(aRK(substep)+bRK(substep))*dtime / (Re*Sc(conc_nu))
              
              if (timeint_mode == 1) then
                 call solve_conc_explicit
              else
                 call solve_conc(.false.,.false.)
              end if
              
           end do
        end if
        
        
        !--- Helmholtz-Multiplikator -------------------------------------------------------------------------
        multL = thetaL*(aRK(substep)+bRK(substep))*dtime / Re
        
        
        !--- Umskalieren (Effizienz, initial guess) ----------------------------------------------------------
        if (.not. init_pre(substep)) pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) * (aRK(substep)+bRK(substep)) * dtime
        
        
        !--- Löser -------------------------------------------------------------------------------------------
!        IF (timeint_mode == 1 .OR. thetaL == 1.) THEN !< thetaL == 1 means implicite ?
        if (timeint_mode == 1 .or. thetaL == 0.) then !< thetaL == 1 means implicite ? danielized doesnot really mather because when thetaL=0 timeintmod=1 is automatically switched on
#ifdef NONBOUSSINESQ
           if (nonBoussinesq_yes) then
              call non_Boussinesq
           else
              call explicit
           end if
#else
           call explicit
#endif
        else
           if (twostep_yes) then
              call twostep
           else
              call outer_iteration
           end if
        end if
        
        !--- physikalischer Druck ----------------------------------------------------------------------------
        pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) / (aRK(substep)+bRK(substep)) / dtime
        
        !--- Undefinierte Ecken / Kanten auffüllen -----------------------------------------------------------
        call fill_corners(pre)
        if (corner_yes) then
           do m = 1, n_conc
              call fill_corners(conc(b1L,b2L,b3L,m))
           end do
        end if
        
     end do RK_loop
     
     !========================================================================================================
     timestep = timestep + 1
     time     = time + dtime
     
     !--- send_signal.txt lesen ------------------------------------------------------------------------------
     call check_signal
     
     
     !--- Druck-Niveau festhalten ----------------------------------------------------------------------------
     call level_pressure
     
     
     !--- Ausschreiben ---------------------------------------------------------------------------------------
     if (write_out_scal) call compute_stats
     if (write_out_vect) call write_fields
     
     
     !--------------------------------------------------------------------------------------------------------
     if (rank == 0 .and. log_iteration_yes) close(10)
     
     if (rank == 0 .and. write_stout_yes .and. timeint_mode == 0) write(*,*)
     if (rank == 0 .and. write_stout_yes .and. timeint_mode == 0) write(*,*)
     
     call MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror) ! notwendig fuer "check_alarm"
     if (finish_yes) exit
     
  end do timeint
  !===========================================================================================================
  
  
  !--- Zeitmessung beenden -----------------------------------------------------------------------------------
  if (rank == 0) then
     call DATE_AND_TIME(values=ctime)
     hour = ctime(5)
     minu = ctime(6)
     sec  = ctime(7)
     msec = ctime(8)
     
     if (ctime(3) /= day) then
        ! Anmerkung: Gilt nur für Jobs <= 24h
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime + 24*60*60*1000
     else
        elatime = msec+1000*(sec+60*(minu+60*hour)) - elatime
     end if
     
     write(99,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,a,i3)') 'Finish time integration at ', ctime(3),'.',ctime(2),'.',ctime(1),    &
                                                           &       ', ',ctime(5),':',ctime(6),':',ctime(7),'.',ctime(8)
     write(99,'(a,E13.5)') 'elapsed time [sec]', REAL(elatime)/1000.
     close(99)
  end if
  
  !--- Restart schreiben -------------------------------------------------------------------------------------
  restart = restart + 1
  call write_restart
  call write_restart_stats
  
  
  !--- Iterationsstatistiken auswerten -----------------------------------------------------------------------
  call iteration_stats
  
  !--- diverse Files schliessen ------------------------------------------------------------------------------
  call close_stats
  
  
  end subroutine timeintegration
  
  
  
end module mod_timeint
