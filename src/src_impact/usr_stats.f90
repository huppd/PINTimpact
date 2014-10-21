!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine open_stats
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_func
  use usr_vars
  
  implicit none
  
  
  !--- open files before time integration starts ---
  ! this is just an example/template
  !
  if (rank == 0 .and. dtime_out_scal /= 0.) then
                            open(21,file='test_TKE_mod_restart'       //restart_char//'.txt',status='UNKNOWN')
                            open(22,file='test_TKE_spectrum_restart'  //restart_char//'.txt',status='UNKNOWN')
                            open(23,file='test_turb_stat_restart'     //restart_char//'.txt',status='UNKNOWN')
     if (concentration_yes) open(24,file='test_conc_stat_restart'     //restart_char//'.txt',status='UNKNOWN')
  end if
  
  
  end subroutine open_stats
  
  
  
  
  
  
  
  
  
  
  
  subroutine close_stats
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_func
  use usr_vars
  
  implicit none
  
  
  !--- close files after time integration ---
  ! this is just an example/template
  !
  if (rank == 0 .and. dtime_out_scal /= 0.) then
                            close(21)
                            close(22)
                            close(23)
     if (concentration_yes) close(24)
  end if
  
  
  end subroutine close_stats
  
  
  
  
  
  
  
  
  
  
  
  subroutine compute_stats
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use mod_exchange
  use mod_diff
  use usr_vars
  use usr_func
  
  implicit none
  
  INCLUDE 'mpif.h'
  
  
  !--- compute your statistics ---
  ! this is just an example/template
  !
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  integer                ::  a, b, m, s
  
  real                   ::  sinus, cosinus
  real                   ::  mult, dd
  
  real                   ::  vel_dist(1:3)
  
  real                   ::  u_mod       (1:2,1:3,-amax:amax,0:bmax)
  real                   ::  u_mod_global(1:2,1:3,-amax:amax,0:bmax)
  
  real                   ::  TKE_mod       (-amax:amax,0:bmax)
  real                   ::  TKE_mod_global(-amax:amax,0:bmax)
  
  real                   ::  TKE(1:2), TKE_global(1:2)
  
  real                   ::  mean       (1:3,S1p:N1p), covar       (1:6,S1p:N1p), diss_visc       (S1p:N1p)
  real                   ::  mean_global(1:3,S1p:N1p), covar_global(1:6,S1p:N1p), diss_visc_global(S1p:N1p)
  real                   ::  mean_write (1:3,1  :M1 ), covar_write (1:6,1  :M1 ), diss_visc_write (1  :M1 )
  real                   ::  deriv      (1:3,1:3)
  
  real                   ::  mass       (1:n_conc,1:2)
  real                   ::  mass_global(1:n_conc,1:2)
  
  real                   ::  mass_flux       (1:n_conc)
  real                   ::  mass_flux_vgl   (1:n_conc)
  real                   ::  mass_flux_global(1:n_conc)
  
  real                   ::  diss, diss_viscInt, diss_viscInt_global
  
  real                   ::  diss_conc       (1:n_conc)
  real                   ::  diss_conc_global(1:n_conc)
  
  
  !===========================================================================================================
  !=== interpolation of velocity to pressure grid points =====================================================
  !===========================================================================================================
  ! vel(:,:,:,i) --> worki(:,:,:)
  if (task == 1) call interpolate_vel(.true.)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== discrete volume =======================================================================================
  !===========================================================================================================
  do k = S3p, N3p
     do j = S2p, N2p
        do i = S1p, N1p
           res(i,j,k) = dx1p(i)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
        end do
     end do
  end do
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== integral energy =======================================================================================
  !===========================================================================================================
  TKE = 0.
  !-----------------------------------------------------------------------------------------------------------
  if (wallnorm_dir == 1 .and. bulkflow_dir == 2) then
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              vel_dist(1) = work1(i,j,k)
              vel_dist(2) = work2(i,j,k) - (L1*x1p(i) - x1p(i)**2)
              vel_dist(3) = work3(i,j,k)
              
              TKE(1) = TKE(1) + res(i,j,k)*(work1(i,j,k)**2 + work2(i,j,k)**2 + work3(i,j,k)**2)
              TKE(2) = TKE(2) + res(i,j,k)*(vel_dist (1)**2 + vel_dist (2)**2 + vel_dist (3)**2)
           end do
        end do
     end do
  !-----------------------------------------------------------------------------------------------------------
  else
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              
              TKE(1) = TKE(1) + res(i,j,k)*(work1(i,j,k)**2 + work2(i,j,k)**2 + work3(i,j,k)**2)
              
           end do
        end do
     end do
  !-----------------------------------------------------------------------------------------------------------
  end if
  
  call MPI_ALLREDUCE(TKE,TKE_global,2,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  TKE = TKE_global / 2.
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== modal energies ========================================================================================
  !===========================================================================================================
  TKE_mod  = 0.
  !-----------------------------------------------------------------------------------------------------------
  if (BC_2L_global == -1 .and. BC_3L_global == -1) then
     
     do i = S1p, N1p
        
        u_mod = 0.
        
        do k = S3p, N3p
           do j = S2p, N2p
              
              do b = 0, bmax
                 do a = -amax, amax
                    cosinus = COS(2.*pi*(REAL(a)*x2p(j)/L2 + REAL(b)*x3p(k)/L3)) * dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
                    sinus   = SIN(2.*pi*(REAL(a)*x2p(j)/L2 + REAL(b)*x3p(k)/L3)) * dx2p(j)*dx3p(k)
                    
                    u_mod(1,1,a,b) = u_mod(1,1,a,b) + work1(i,j,k) * cosinus
                    u_mod(1,2,a,b) = u_mod(1,2,a,b) + work2(i,j,k) * cosinus
                    u_mod(1,3,a,b) = u_mod(1,3,a,b) + work3(i,j,k) * cosinus
                    
                    u_mod(2,1,a,b) = u_mod(2,1,a,b) - work1(i,j,k) * sinus
                    u_mod(2,2,a,b) = u_mod(2,2,a,b) - work2(i,j,k) * sinus
                    u_mod(2,3,a,b) = u_mod(2,3,a,b) - work3(i,j,k) * sinus
                 end do
              end do
              
           end do
        end do
        
        call MPI_ALLREDUCE(u_mod,u_mod_global,6*(2*amax+1)*(bmax+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
        
        do b = 0, bmax
           do a = -amax, amax
              TKE_mod_global(a,b) =  u_mod_global(1,1,a,b)**2 + u_mod_global(2,1,a,b)**2      &
                          &        + u_mod_global(1,2,a,b)**2 + u_mod_global(2,2,a,b)**2      &
                          &        + u_mod_global(1,3,a,b)**2 + u_mod_global(2,3,a,b)**2
           end do
        end do
        TKE_mod = TKE_mod + TKE_mod_global * dx1p(i)
     end do
     
     call MPI_ALLREDUCE(TKE_mod,TKE_mod_global,(2*amax+1)*(bmax+1),MPI_REAL8,MPI_SUM,COMM_BAR1,merror)
     TKE_mod = TKE_mod_global / (2.*L1*(L2*L3)**2)
     
     if (rank == 0) then
        write(21,'(100E25.17)') time, TKE_mod(-amax:amax,0:bmax)
        call flush(21)
        
        do j = 0, bmax
           do i = -amax, amax
              write(22,'(2i4,E25.17)') i, j, TKE_mod(i,j)
           end do
           write(22,*)
        end do
        
        write(22,*)
        call flush(22)
     end if
     
  end if
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== turbulence statistics =================================================================================
  !===========================================================================================================
  if (wallnorm_dir == 1) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- <u>_xy(z) ------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mean = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              mean(1,i) = mean(1,i) + work1(i,j,k)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              mean(2,i) = mean(2,i) + work2(i,j,k)*dx2p(j)*dx3p(k)
              mean(3,i) = mean(3,i) + work3(i,j,k)*dx2p(j)*dx3p(k)
           end do
        end do
     end do
     
     call MPI_ALLREDUCE(mean,mean_global,3*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     mean = mean_global / (L2*L3)
     
     call MPI_ALLGATHERv(mean,3*(N1p-S1p+1),MPI_REAL8,mean_write,3*bar1_size,3*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- <u'i u'j>_xy(z) ------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     covar = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              covar(1,i) = covar(1,i) + (work1(i,j,k)-mean(1,i))*(work1(i,j,k)-mean(1,i))*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              covar(2,i) = covar(2,i) + (work2(i,j,k)-mean(2,i))*(work2(i,j,k)-mean(2,i))*dx2p(j)*dx3p(k)
              covar(3,i) = covar(3,i) + (work3(i,j,k)-mean(3,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
              
              covar(4,i) = covar(4,i) + (work1(i,j,k)-mean(1,i))*(work2(i,j,k)-mean(2,i))*dx2p(j)*dx3p(k)
              covar(5,i) = covar(5,i) + (work1(i,j,k)-mean(1,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
              covar(6,i) = covar(6,i) + (work2(i,j,k)-mean(2,i))*(work3(i,j,k)-mean(3,i))*dx2p(j)*dx3p(k)
           end do
        end do
     end do
     
     call MPI_ALLREDUCE(covar,covar_global,6*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     covar = covar_global / (L2*L3)
     
     call MPI_ALLGATHERv(covar,6*(N1p-S1p+1),MPI_REAL8,covar_write,6*bar1_size,6*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     !--------------------------------------------------------------------------------------------------------
     !--- dissipation(z), dissipation_total ------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     diss_visc    = 0.
     diss_viscInt = 0.
     
     if (M2 == 2) then
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 
                 deriv = 0.
                 do ii = b1L, b1U
                    deriv(1,1) = deriv(1,1) + cp1(ii,i)*work1(i+ii,j,k)
                    deriv(2,1) = deriv(2,1) + cp1(ii,i)*work2(i+ii,j,k)
                 end do
                 do jj = b2L, b2U
                    deriv(1,2) = deriv(1,2) + cp2(jj,j)*work1(i,j+jj,k)
                    deriv(2,2) = deriv(2,2) + cp2(jj,j)*work2(i,j+jj,k)
                 end do
                 
                 diss =    2.*deriv(1,1)**2 +    deriv(1,2)**2    &
                     &  +     deriv(2,1)**2 + 2.*deriv(2,2)**2    &
                     &  +  2.*deriv(1,2)*deriv(2,1)
                 
                 diss_viscInt = diss_viscInt + diss*res(i,j,k)
                 diss_visc(i) = diss_visc(i) + diss*dx2p(j)
              end do
           end do
        end do
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 
                 deriv = 0.
                 do ii = b1L, b1U
                    deriv(1,1) = deriv(1,1) + cp1(ii,i)*work1(i+ii,j,k)
                    deriv(2,1) = deriv(2,1) + cp1(ii,i)*work2(i+ii,j,k)
                    deriv(3,1) = deriv(3,1) + cp1(ii,i)*work3(i+ii,j,k)
                 end do
                 do jj = b2L, b2U
                    deriv(1,2) = deriv(1,2) + cp2(jj,j)*work1(i,j+jj,k)
                    deriv(2,2) = deriv(2,2) + cp2(jj,j)*work2(i,j+jj,k)
                    deriv(3,2) = deriv(3,2) + cp2(jj,j)*work3(i,j+jj,k)
                 end do
                 do kk = b3L, b3U
                    deriv(1,3) = deriv(1,3) + cp3(kk,k)*work1(i,j,k+kk)
                    deriv(2,3) = deriv(2,3) + cp3(kk,k)*work2(i,j,k+kk)
                    deriv(3,3) = deriv(3,3) + cp3(kk,k)*work3(i,j,k+kk)
                 end do
                 
                 diss =    2.*deriv(1,1)**2 +    deriv(1,2)**2 +    deriv(1,3)**2                      &
                     &  +     deriv(2,1)**2 + 2.*deriv(2,2)**2 +    deriv(2,3)**2                      &
                     &  +     deriv(3,1)**2 +    deriv(3,2)**2 + 2.*deriv(3,3)**2                      &
                     &  + 2.*(deriv(1,2)*deriv(2,1) + deriv(1,3)*deriv(3,1) + deriv(2,3)*deriv(3,2))
                 
                 diss_viscInt = diss_viscInt + diss*res(i,j,k)
                 diss_visc(i) = diss_visc(i) + diss*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
              end do
           end do
        end do
     end if
     
     call MPI_REDUCE(diss_viscInt,diss_viscInt_global,1,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     
     if (rank == 0) then
        diss_viscInt = - diss_viscInt_global / Re
        mult = 0.5*dtime_out_scal
        energy_visc      = energy_visc + mult*(diss_viscInt + diss_viscInt_old)
        diss_viscInt_old = diss_viscInt
     end if
     
     call MPI_ALLREDUCE(diss_visc,diss_visc_global,1*(N1p-S1p+1),MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     diss_visc = - diss_visc_global / (Re*L2*L3)
     call MPI_ALLGATHERv(diss_visc,1*(N1p-S1p+1),MPI_REAL8,diss_visc_write,1*bar1_size,1*bar1_offset,MPI_REAL8,COMM_BAR1,merror)
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- output ---------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (rank == 0) then
        write(23,'(a,1i10,1E25.17)') '#', timestep, time
        do i = 1, M1
           write(23,'(11E25.17)') y1p(i), mean_write(1:3,i), covar_write(1:6,i), diss_visc_write(i)
        end do
        write(23,*)
        call flush(23)
     end if
     
  end if
  !===========================================================================================================
  
  
  
  
  !===========================================================================================================
  !=== concentration statistics ==============================================================================
  !===========================================================================================================
  if (wallnorm_dir == 1 .and. concentration_yes) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- ghost cell exchange --------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     do m = 1, n_conc
        call exchange(1,0,conc(b1L,b2L,b3L,m))
        call exchange(2,0,conc(b1L,b2L,b3L,m))
        call exchange(3,0,conc(b1L,b2L,b3L,m))
     end do
     
     !--------------------------------------------------------------------------------------------------------
     !--- mass / potential energy / Stokes dissipation -------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     mass = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              mass(:,1) = mass(:,1) + conc(i,j,k,:)*res(i,j,k)
              mass(:,2) = mass(:,2) + conc(i,j,k,:)*res(i,j,k)*x1p(i)
           end do
        end do
     end do
     
     call MPI_REDUCE(mass,mass_global,2*n_conc,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     mass = mass_global
     
     if (rank == 0) then
        do m = 1, n_conc
           mass(m,:) = Ric(m)*mass(m,:)
           mult = 0.5*dtime_out_scal * us_vec(1,m)
           
           energy_stokes(m) = energy_stokes(m) + mult*(mass(m,1) + mass_old(m))
           mass_flux_vgl(m) = (mass_old(m) - mass(m,1))/dtime_out_scal
           
           mass_old(m) = mass(m,1)
        end do
     end if
     
     !--------------------------------------------------------------------------------------------------------
     !--- potential energy (diffusion) -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     diss_conc = 0.
     
     do m = 1, n_conc
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 dd = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
                 do ii = b1L+1, b1U
                    dd = dd + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 end do
                 do jj = b2L, b2U
                    dd = dd + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 end do
                 if (M2 /= 2) then
                    do kk = b3L, b3U
                       dd = dd + cc33(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                 end if
                 diss_conc(m) = diss_conc(m) + dd*res(i,j,k)*x1p(i)
              end do
           end do
        end do
        
        diss_conc(m) = Ric(m)*diss_conc(m)/(Re*Sc(m))
     end do
     
     call MPI_REDUCE(diss_conc,diss_conc_global,n_conc,MPI_REAL8,MPI_SUM,0,COMM_CART,merror)
     
     if (rank == 0) then
        diss_conc = diss_conc_global
        mult = 0.5*dtime_out_scal
        do m = 1, n_conc
           Epot_diff(m) = Epot_diff(m) + mult*(diss_conc(m) + diss_conc_old(m))
           diss_conc_old(m) = diss_conc(m)
        end do
     end if
     
     !--------------------------------------------------------------------------------------------------------
     !--- deposit and mass flux ------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if ((dtime_out_scal == 0. .and. dtime_out_vect > 0.) .or.   &
       & (dtime_out_scal /= 0. .and. dtime_out_vect > 0. .and. MOD(dtime_out_vect,dtime_out_scal) > 10.**(-15))) then
        if (rank == 0) write(*,*) 'ERROR! For time integration of particle deposit is MOD(time_out_vect,time_out_scal) == 0. mandatory!'
        call MPI_FINALIZE(merror)
        stop
     end if
     
     do m = 1, n_conc
        mult = dtime_out_scal*us_vec(1,m)/2.
        
        dep1L_conc(S2p:N2p,S3p:N3p,m) = dep1L_conc(S2p:N2p,S3p:N3p,m) + mult*(conc(S1p,S2p:N2p,S3p:N3p,m) + conc1L(S2p:N2p,S3p:N3p,m))
        conc1L    (S2p:N2p,S3p:N3p,m) = conc  (S1p,S2p:N2p,S3p:N3p,m)
     end do
     
     mass_flux_global = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           mass_flux_global(:) = mass_flux_global(:) - Ric(:)*us_vec(1,:)*conc(S1p,j,k,:)*dx2p(j)*dx3p(k) ! (dx3p = 1. for 2D)
        end do
     end do
     
     call MPI_ALLREDUCE(mass_flux_global,mass_flux,n_conc,MPI_REAL8,MPI_SUM,COMM_SLICE1,merror)
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- output ---------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (rank == 0) then
        write(24,'(50E25.17)') time, TKE(1), energy_visc, energy_stokes(:), mass(:,2), mass(:,1), mass_flux(:), mass_flux_vgl(:), Epot_diff(:)
        call flush(24)
     end if
     
  end if
  !===========================================================================================================
  
  
  write_out_scal = .false.
  time_out_scal = time_out_scal + dtime_out_scal
  
  
  end subroutine compute_stats
  
  
  
  
  
  
  
  
  
  
  
  subroutine write_restart_stats
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use mod_lib
  use usr_vars
  use usr_func
  
  implicit none
  
  
  !--- write time-dependent data after time integration for next restart ---
  ! this is just an example/template
  !
  character(len=3)       ::  next_restart_char
  
  
  if (write_restart_yes .and. rank == 0) then
     
     !--- new restart number as string for file names ---
     call num_to_string(3,restart,next_restart_char)
     
     open (10,file='stats_restart'//next_restart_char//'.txt',status='UNKNOWN')
     
                            write(10,'(a,10E25.17)') 'energy_visc:      ', energy_visc 
                            write(10,'(a,10E25.17)') 'diss_viscInt_old: ', diss_viscInt_old
     if (concentration_yes) write(10,'(a,10E25.17)') 'mass_old:         ', mass_old     (1:n_conc) 
     if (concentration_yes) write(10,'(a,10E25.17)') 'energy_stokes:    ', energy_stokes(1:n_conc)
     if (concentration_yes) write(10,'(a,10E25.17)') 'Epot_diff:        ', Epot_diff    (1:n_conc)    
     if (concentration_yes) write(10,'(a,10E25.17)') 'diss_conc_old:    ', diss_conc_old(1:n_conc)
     
     close(10)
  end if
  
  
  end subroutine write_restart_stats
  
  
  
  
  
  
  
  
  
  
  
  subroutine read_restart_stats
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  INCLUDE 'mpif.h'
  
  
  !--- read time-dependent data from previous restart before time integration starts ---
  ! this is just an example/template
  !
  character(len=80)        ::  dummy
  integer                  ::  ios
  
  
  open(10,file='stats_restart'//restart_char//'.txt',action='read',status='old',iostat=ios)
  
  if (ios /= 0) then
     if (rank == 0) write(*,*) 'ERROR! Cannot open the stats_restart'//restart_char//'.txt file!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
                         read(10,*) dummy, energy_visc
                         read(10,*) dummy, diss_viscInt_old
  if (concentration_yes) read(10,*) dummy, mass_old     (1:n_conc_old)
  if (concentration_yes) read(10,*) dummy, energy_stokes(1:n_conc_old)
  if (concentration_yes) read(10,*) dummy, Epot_diff    (1:n_conc_old)
  if (concentration_yes) read(10,*) dummy, diss_conc_old(1:n_conc_old)
  
  close(10)
  
  
  if (rank == 0) then
                            write(*,'(a,10E25.17)') 'energy_visc:      ', energy_visc
                            write(*,'(a,10E25.17)') 'diss_viscInt_old: ', diss_viscInt_old
     if (concentration_yes) write(*,'(a,10E25.17)') 'mass_old:         ', mass_old     (1:n_conc_old)
     if (concentration_yes) write(*,'(a,10E25.17)') 'energy_stokes:    ', energy_stokes(1:n_conc_old)
     if (concentration_yes) write(*,'(a,10E25.17)') 'Epot_diff:        ', Epot_diff    (1:n_conc_old)
     if (concentration_yes) write(*,'(a,10E25.17)') 'diss_conc_old:    ', diss_conc_old(1:n_conc_old)
     write(*,*)
     write(*,*)
  end if
  
  
  end subroutine read_restart_stats
