!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE open_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  
  IMPLICIT NONE
  
  
  !--- open files before time integration starts ---
  ! this is just an example/template
  !
  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
                            OPEN(21,FILE='test_TKE_mod_restart'       //restart_char//'.txt',STATUS='UNKNOWN')
                            OPEN(22,FILE='test_TKE_spectrum_restart'  //restart_char//'.txt',STATUS='UNKNOWN')
                            OPEN(23,FILE='test_turb_stat_restart'     //restart_char//'.txt',STATUS='UNKNOWN')
     IF (concentration_yes) OPEN(24,FILE='test_conc_stat_restart'     //restart_char//'.txt',STATUS='UNKNOWN')
  END IF
  
  
  END SUBROUTINE open_stats
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE close_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_func
  USE usr_vars
  
  IMPLICIT NONE
  
  
  !--- close files after time integration ---
  ! this is just an example/template
  !
  IF (rank == 0 .AND. dtime_out_scal /= 0.) THEN
                            CLOSE(21)
                            CLOSE(22)
                            CLOSE(23)
     IF (concentration_yes) CLOSE(24)
  END IF
  
  
  END SUBROUTINE close_stats
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE compute_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  
  !--- compute your statistics ---
  ! this is just an example/template
  !
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  INTEGER                ::  a, b, m, s
  
  REAL                   ::  sinus, cosinus
  REAL                   ::  mult, dd
  
  REAL                   ::  vel_dist(1:3)
  
  REAL                   ::  u_mod       (1:2,1:3,-amax:amax,0:bmax)
  REAL                   ::  u_mod_global(1:2,1:3,-amax:amax,0:bmax)
  
  REAL                   ::  TKE_mod       (-amax:amax,0:bmax)
  REAL                   ::  TKE_mod_global(-amax:amax,0:bmax)
  
  REAL                   ::  TKE(1:2), TKE_global(1:2)
  
  REAL                   ::  mean       (1:3,S1p:N1p), covar       (1:6,S1p:N1p), diss_visc       (S1p:N1p)
  REAL                   ::  mean_global(1:3,S1p:N1p), covar_global(1:6,S1p:N1p), diss_visc_global(S1p:N1p)
  REAL                   ::  mean_write (1:3,1  :M1 ), covar_write (1:6,1  :M1 ), diss_visc_write (1  :M1 )
  REAL                   ::  deriv      (1:3,1:3)
  
  REAL                   ::  mass       (1:n_conc,1:2)
  REAL                   ::  mass_global(1:n_conc,1:2)
  
  REAL                   ::  mass_flux       (1:n_conc)
  REAL                   ::  mass_flux_vgl   (1:n_conc)
  REAL                   ::  mass_flux_global(1:n_conc)
  
  REAL                   ::  diss, diss_viscInt, diss_viscInt_global
  
  REAL                   ::  diss_conc       (1:n_conc)
  REAL                   ::  diss_conc_global(1:n_conc)
  
  
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  IF (task == 1) CALL interpolate_vel(.TRUE.)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== discrete volume =======================================================================================
  !===========================================================================================================
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           res(i,j,k) = dx1p(i)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
        END DO
     END DO
  END DO
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== integral energy =======================================================================================
  !===========================================================================================================
  TKE = 0.
  !-----------------------------------------------------------------------------------------------------------
  IF (wallnorm_dir == 1 .AND. bulkflow_dir == 2) THEN
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              vel_dist(1) = work1(i,j,k)
              vel_dist(2) = work2(i,j,k) - (L1*x1p(i) - x1p(i)**2)
              vel_dist(3) = work3(i,j,k)
              
              TKE(1) = TKE(1) + res(i,j,k)*(work1(i,j,k)**2 + work2(i,j,k)**2 + work3(i,j,k)**2)
              TKE(2) = TKE(2) + res(i,j,k)*(vel_dist (1)**2 + vel_dist (2)**2 + vel_dist (3)**2)
           END DO
        END DO
     END DO
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              
              TKE(1) = TKE(1) + res(i,j,k)*(work1(i,j,k)**2 + work2(i,j,k)**2 + work3(i,j,k)**2)
              
           END DO
        END DO
     END DO
  !-----------------------------------------------------------------------------------------------------------
  END IF
  
  CALL MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE = TKE_global / 2.
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== modal energies ========================================================================================
  !===========================================================================================================
  TKE_mod  = 0.
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L_global == -1 .AND. BC_3L_global == -1) THEN
     
     DO i = S1p, N1p
        
        u_mod = 0.
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              
              DO b = 0, bmax
                 DO a = -amax, amax
                    cosinus = COS(2.*pi*(REAL(a)*x2p(j)/L2 + REAL(b)*x3p(k)/L3)) * dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
                    sinus   = SIN(2.*pi*(REAL(a)*x2p(j)/L2 + REAL(b)*x3p(k)/L3)) * dx2p(j)*dx3p(k)
                    
                    u_mod(1,1,a,b) = u_mod(1,1,a,b) + work1(i,j,k) * cosinus
                    u_mod(1,2,a,b) = u_mod(1,2,a,b) + work2(i,j,k) * cosinus
                    u_mod(1,3,a,b) = u_mod(1,3,a,b) + work3(i,j,k) * cosinus
                    
                    u_mod(2,1,a,b) = u_mod(2,1,a,b) - work1(i,j,k) * sinus
                    u_mod(2,2,a,b) = u_mod(2,2,a,b) - work2(i,j,k) * sinus
                    u_mod(2,3,a,b) = u_mod(2,3,a,b) - work3(i,j,k) * sinus
                 END DO
              END DO
              
           END DO
        END DO
        
        CALL MPI_ALLREDUCE(u_mod,u_mod_global,6*(2*amax+1)*(bmax+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
        
        DO b = 0, bmax
           DO a = -amax, amax
              TKE_mod_global(a,b) =  u_mod_global(1,1,a,b)**2 + u_mod_global(2,1,a,b)**2      &
                          &        + u_mod_global(1,2,a,b)**2 + u_mod_global(2,2,a,b)**2      &
                          &        + u_mod_global(1,3,a,b)**2 + u_mod_global(2,3,a,b)**2
           END DO
        END DO
        TKE_mod = TKE_mod + TKE_mod_global * dx1p(i)
     END DO
     
     CALL MPI_ALLREDUCE(TKE_mod,TKE_mod_global,(2*amax+1)*(bmax+1),MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
     TKE_mod = TKE_mod_global / (2.*L1*(L2*L3)**2)
     
     IF (rank == 0) THEN
        WRITE(21,'(100E25.17)') time, TKE_mod(-amax:amax,0:bmax)
        CALL flush(21)
        
        DO j = 0, bmax
           DO i = -amax, amax
              WRITE(22,'(2i4,E25.17)') i, j, TKE_mod(i,j)
           END DO
           WRITE(22,*)
        END DO
        
        WRITE(22,*)
        CALL flush(22)
     END IF
     
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== turbulence statistics =================================================================================
  !===========================================================================================================
  IF (wallnorm_dir == 1) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- <u>_xy(z) ------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mean = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              mean(1,i) = mean(1,i) + work1(i,j,k)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              mean(2,i) = mean(2,i) + work2(i,j,k)*dx2p(j)*dx3p(k)
              mean(3,i) = mean(3,i) + work3(i,j,k)*dx2p(j)*dx3p(k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(mean,mean_global,3*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     mean = mean_global / (L2*L3)
     
     CALL MPI_ALLGATHERv(mean,3*(N1p-S1p+1),MPI_REAL8,mean_write,3*bar1_size,3*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- <u'i u'j>_xy(z) ------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     covar = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              covar(1,i) = covar(1,i) + (work1(i,j,k)-mean(1,i))*(work1(i,j,k)-mean(1,i))*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              covar(2,i) = covar(2,i) + (work2(i,j,k)-mean(2,i))*(work2(i,j,k)-mean(2,i))*dx2p(j)*dx3p(k)
              covar(3,i) = covar(3,i) + (work3(i,j,k)-mean(3,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
              
              covar(4,i) = covar(4,i) + (work1(i,j,k)-mean(1,i))*(work2(i,j,k)-mean(2,i))*dx2p(j)*dx3p(k)
              covar(5,i) = covar(5,i) + (work1(i,j,k)-mean(1,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
              covar(6,i) = covar(6,i) + (work2(i,j,k)-mean(2,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
           END DO
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(covar,covar_global,6*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     covar = covar_global / (L2*L3)
     
     CALL MPI_ALLGATHERv(covar,6*(N1p-S1p+1),MPI_REAL8,covar_write,6*bar1_size,6*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- dissipation(z), dissipation_total ------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     diss_visc    = 0.
     diss_viscInt = 0.
     
     IF (M2 == 2) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 
                 deriv = 0.
                 DO ii = b1L, b1U
                    deriv(1,1) = deriv(1,1) + cp1(ii,i)*work1(i+ii,j,k)
                    deriv(2,1) = deriv(2,1) + cp1(ii,i)*work2(i+ii,j,k)
                 END DO
                 DO jj = b2L, b2U
                    deriv(1,2) = deriv(1,2) + cp2(jj,j)*work1(i,j+jj,k)
                    deriv(2,2) = deriv(2,2) + cp2(jj,j)*work2(i,j+jj,k)
                 END DO
                 
                 diss =    2.*deriv(1,1)**2 +    deriv(1,2)**2    &
                     &  +     deriv(2,1)**2 + 2.*deriv(2,2)**2    &
                     &  +  2.*deriv(1,2)*deriv(2,1)
                 
                 diss_viscInt = diss_viscInt + diss*res(i,j,k)
                 diss_visc(i) = diss_visc(i) + diss*dx2p(j)
              END DO
           END DO
        END DO
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 
                 deriv = 0.
                 DO ii = b1L, b1U
                    deriv(1,1) = deriv(1,1) + cp1(ii,i)*work1(i+ii,j,k)
                    deriv(2,1) = deriv(2,1) + cp1(ii,i)*work2(i+ii,j,k)
                    deriv(3,1) = deriv(3,1) + cp1(ii,i)*work3(i+ii,j,k)
                 END DO
                 DO jj = b2L, b2U
                    deriv(1,2) = deriv(1,2) + cp2(jj,j)*work1(i,j+jj,k)
                    deriv(2,2) = deriv(2,2) + cp2(jj,j)*work2(i,j+jj,k)
                    deriv(3,2) = deriv(3,2) + cp2(jj,j)*work3(i,j+jj,k)
                 END DO
                 DO kk = b3L, b3U
                    deriv(1,3) = deriv(1,3) + cp3(kk,k)*work1(i,j,k+kk)
                    deriv(2,3) = deriv(2,3) + cp3(kk,k)*work2(i,j,k+kk)
                    deriv(3,3) = deriv(3,3) + cp3(kk,k)*work3(i,j,k+kk)
                 END DO
                 
                 diss =    2.*deriv(1,1)**2 +    deriv(1,2)**2 +    deriv(1,3)**2                      &
                     &  +     deriv(2,1)**2 + 2.*deriv(2,2)**2 +    deriv(2,3)**2                      &
                     &  +     deriv(3,1)**2 +    deriv(3,2)**2 + 2.*deriv(3,3)**2                      &
                     &  + 2.*(deriv(1,2)*deriv(2,1) + deriv(1,3)*deriv(3,1) + deriv(2,3)*deriv(3,2))
                 
                 diss_viscInt = diss_viscInt + diss*res(i,j,k)
                 diss_visc(i) = diss_visc(i) + diss*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              END DO
           END DO
        END DO
     END IF
     
     CALL MPI_REDUCE(diss_viscInt,diss_viscInt_global,1,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     
     IF (rank == 0) THEN
        diss_viscInt = - diss_viscInt_global / Re
        mult = 0.5*dtime_out_scal
        energy_visc      = energy_visc + mult*(diss_viscInt + diss_viscInt_old)
        diss_viscInt_old = diss_viscInt
     END IF
     
     CALL MPI_ALLREDUCE(diss_visc,diss_visc_global,1*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     diss_visc = - diss_visc_global / (Re*L2*L3)
     CALL MPI_ALLGATHERv(diss_visc,1*(N1p-S1p+1),MPI_REAL8,diss_visc_write,1*bar1_size,1*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- output ---------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (rank == 0) THEN
        WRITE(23,'(a,1i10,1E25.17)') '#', timestep, time
        DO i = 1, M1
           WRITE(23,'(11E25.17)') y1p(i), mean_write(1:3,i), covar_write(1:6,i), diss_visc_write(i)
        END DO
        WRITE(23,*)
        CALL flush(23)
     END IF
     
  END IF
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== concentration statistics ==============================================================================
  !===========================================================================================================
  IF (wallnorm_dir == 1 .AND. concentration_yes) THEN
     
     !--------------------------------------------------------------------------------------------------------
     !--- ghost cell exchange --------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     DO m = 1, n_conc
        CALL exchange(1,0,conc(b1L,b2L,b3L,m))
        CALL exchange(2,0,conc(b1L,b2L,b3L,m))
        CALL exchange(3,0,conc(b1L,b2L,b3L,m))
     END DO
     
     !--------------------------------------------------------------------------------------------------------
     !--- mass / potential energy / Stokes dissipation -------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mass = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              mass(:,1) = mass(:,1) + conc(i,j,k,:)*res(i,j,k)
              mass(:,2) = mass(:,2) + conc(i,j,k,:)*res(i,j,k)*x1p(i)
           END DO
        END DO
     END DO
     
     CALL MPI_REDUCE(mass,mass_global,2*n_conc,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     mass = mass_global
     
     IF (rank == 0) THEN
        DO m = 1, n_conc
           mass(m,:) = Ric(m)*mass(m,:)
           mult = 0.5*dtime_out_scal * us_vec(1,m)
           
           energy_stokes(m) = energy_stokes(m) + mult*(mass(m,1) + mass_old(m))
           mass_flux_vgl(m) = (mass_old(m) - mass(m,1))/dtime_out_scal
           
           mass_old(m) = mass(m,1)
        END DO
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     !--- potential energy (diffusion) -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     diss_conc = 0.
     
     DO m = 1, n_conc
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
                 DO ii = b1L+1, b1U
                    dd = dd + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 END DO
                 DO jj = b2L, b2U
                    dd = dd + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 END DO
                 IF (M2 /= 2) THEN
                    DO kk = b3L, b3U
                       dd = dd + cc33(kk,k,m)*conc(i,j,k+kk,m)
                    END DO
                 END IF
                 diss_conc(m) = diss_conc(m) + dd*res(i,j,k)*x1p(i)
              END DO
           END DO
        END DO
        
        diss_conc(m) = Ric(m)*diss_conc(m)/(Re*Sc(m))
     END DO
     
     CALL MPI_REDUCE(diss_conc,diss_conc_global,n_conc,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     
     IF (rank == 0) THEN
        diss_conc = diss_conc_global
        mult = 0.5*dtime_out_scal
        DO m = 1, n_conc
           Epot_diff(m) = Epot_diff(m) + mult*(diss_conc(m) + diss_conc_old(m))
           diss_conc_old(m) = diss_conc(m)
        END DO
     END IF
     
     !--------------------------------------------------------------------------------------------------------
     !--- deposit and mass flux ------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF ((dtime_out_scal == 0. .AND. dtime_out_vect > 0.) .OR.   &
       & (dtime_out_scal /= 0. .AND. dtime_out_vect > 0. .AND. MOD(dtime_out_vect,dtime_out_scal) > 10.**(-15))) THEN
        IF (rank == 0) WRITE(*,*) 'ERROR! For time integration of particle deposit is MOD(time_out_vect,time_out_scal) == 0. mandatory!'
        CALL MPI_FINALIZE(merror)
        STOP
     END IF
     
     DO m = 1, n_conc
        mult = dtime_out_scal*us_vec(1,m)/2.
        
        dep1L_conc(S2p:N2p,S3p:N3p,m) = dep1L_conc(S2p:N2p,S3p:N3p,m) + mult*(conc(S1p,S2p:N2p,S3p:N3p,m) + conc1L(S2p:N2p,S3p:N3p,m))
        conc1L    (S2p:N2p,S3p:N3p,m) = conc  (S1p,S2p:N2p,S3p:N3p,m)
     END DO
     
     mass_flux_global = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           mass_flux_global(:) = mass_flux_global(:) - Ric(:)*us_vec(1,:)*conc(S1p,j,k,:)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
        END DO
     END DO
     
     CALL MPI_ALLREDUCE(mass_flux_global,mass_flux,n_conc,MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- output ---------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (rank == 0) THEN
        WRITE(24,'(50E25.17)') time, TKE(1), energy_visc, energy_stokes(:), mass(:,2), mass(:,1), mass_flux(:), mass_flux_vgl(:), Epot_diff(:)
        CALL flush(24)
     END IF
     
  END IF
  !===========================================================================================================
  
  
  write_out_scal = .FALSE.
  time_out_scal = time_out_scal + dtime_out_scal
  
  
  END SUBROUTINE compute_stats
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE write_restart_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE mod_lib
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  
  !--- write time-dependent data after time integration for next restart ---
  ! this is just an example/template
  !
  CHARACTER(LEN=3)       ::  next_restart_char
  
  
  IF (write_restart_yes .AND. rank == 0) THEN
     
     !--- new restart number as string for file names ---
     CALL num_to_string(3,restart,next_restart_char)
     
     OPEN (10,FILE='stats_restart'//next_restart_char//'.txt',STATUS='UNKNOWN')
     
                            WRITE(10,'(a,10E25.17)') 'energy_visc:      ', energy_visc 
                            WRITE(10,'(a,10E25.17)') 'diss_viscInt_old: ', diss_viscInt_old
     IF (concentration_yes) WRITE(10,'(a,10E25.17)') 'mass_old:         ', mass_old     (1:n_conc) 
     IF (concentration_yes) WRITE(10,'(a,10E25.17)') 'energy_stokes:    ', energy_stokes(1:n_conc)
     IF (concentration_yes) WRITE(10,'(a,10E25.17)') 'Epot_diff:        ', Epot_diff    (1:n_conc)    
     IF (concentration_yes) WRITE(10,'(a,10E25.17)') 'diss_conc_old:    ', diss_conc_old(1:n_conc)
     
     CLOSE(10)
  END IF
  
  
  END SUBROUTINE write_restart_stats
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE read_restart_stats
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  
  !--- read time-dependent data from previous restart before time integration starts ---
  ! this is just an example/template
  !
  CHARACTER(LEN=80)        ::  dummy
  INTEGER                  ::  ios
  
  
  OPEN(10,FILE='stats_restart'//restart_char//'.txt',ACTION='read',STATUS='old',IOSTAT=ios)
  
  IF (ios /= 0) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Cannot open the stats_restart'//restart_char//'.txt file!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
                         READ(10,*) dummy, energy_visc
                         READ(10,*) dummy, diss_viscInt_old
  IF (concentration_yes) READ(10,*) dummy, mass_old     (1:n_conc_old)
  IF (concentration_yes) READ(10,*) dummy, energy_stokes(1:n_conc_old)
  IF (concentration_yes) READ(10,*) dummy, Epot_diff    (1:n_conc_old)
  IF (concentration_yes) READ(10,*) dummy, diss_conc_old(1:n_conc_old)
  
  CLOSE(10)
  
  
  IF (rank == 0) THEN
                            WRITE(*,'(a,10E25.17)') 'energy_visc:      ', energy_visc
                            WRITE(*,'(a,10E25.17)') 'diss_viscInt_old: ', diss_viscInt_old
     IF (concentration_yes) WRITE(*,'(a,10E25.17)') 'mass_old:         ', mass_old     (1:n_conc_old)
     IF (concentration_yes) WRITE(*,'(a,10E25.17)') 'energy_stokes:    ', energy_stokes(1:n_conc_old)
     IF (concentration_yes) WRITE(*,'(a,10E25.17)') 'Epot_diff:        ', Epot_diff    (1:n_conc_old)
     IF (concentration_yes) WRITE(*,'(a,10E25.17)') 'diss_conc_old:    ', diss_conc_old(1:n_conc_old)
     WRITE(*,*)
     WRITE(*,*)
  END IF
  
  
  END SUBROUTINE read_restart_stats
