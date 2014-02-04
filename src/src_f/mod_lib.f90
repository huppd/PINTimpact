!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module mod_lib
  
    use iso_c_binding
  
  use mod_dims
  use mod_vars
  use mod_exchange
  use mod_diff
  use mod_coeffs
  
  
  private
  
  public get_dtime, get_beta
  public init_BC
  public fill_corners, level_pressure
  public init_alarm, check_alarm, check_signal
  public num_to_string
  public iteration_stats
  
  
  INCLUDE 'mpif.h'
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  !> \brief computes time step according to CFL
  !! \test TEST!!! Noch nicht 100% durchgetested ...
  subroutine get_dtime
  
  implicit none
  
  integer                ::  m
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  real                   ::  alpha_re, alpha_im
  real                   ::  beta(1:2), beta_global(1:2)
  real                   ::  angle, radius
  
  real                   ::  kdx_max_conv(1:3)
  real                   ::  kdx_max_visc(1:3)
  
  real                   ::  array_real(1:2)
  logical                ::  array_log (1:4)
  
  real                   ::  newtime
  real                   ::  pi
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Die Auswertung bezieht sich generell nur auf den Feldbereich. Bei Ausfluss-Rand-          !
  !                bedingungen ist die CFL-Bedingung naturgemäss irrelevant, Einfluss-Randbedingungen        !
  !                sollten generell implizit in der Zeit diskretisiert werden (Dirichlet-, Neumann- oder     !
  !                Robin-Randbedingungen).                                                                   !
  !              - Bei den Geschwindigkeiten werden der Einfachheit halber nur die Druck-Gitterpunkte be-    !
  !                trachtet, nicht die eigentlich Geschwindigkeitsgitter. Die Abschaetzung sollte aber auf   !
  !                der sicheren Seite liegen.                                                                !
  !              - Explizite Differenzen: k^2_mod_max <= ||L||_inf                                           !
  !              - Implizite Differenzen: k^2_mod_max <= max((kdx_max/dx)^2)  (angenommen, nicht exakt!)     !
  !              - Imaginärteile (z.B. bei schiefen Randstencils) werden der Einfachheit halber vernach-     !
  !                lässigt.                                                                                  !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  dtime_old = dtime
  
  if ((Int_dtime /= 0 .and. MOD(timestep,Int_dtime) == 0) .or. new_dtime) then
     
     if (rank == 0 .and. write_stout_yes) write(*,'(a)') 'new dtime ...'
     
     new_dtime = .false.
     
     
     kdx_max_conv = pi
     kdx_max_visc = pi
     
     
     ! vel(:,:,:,i) --> worki(:,:,:)
     call interpolate_vel(.true.)
     
     beta = 0.
     
     
     !========================================================================================================
     do k = S31, N31 ! Nur Feldbereich! (Intervallgrenzen von Tangentialgeschwindigkeiten werden verwendet)
        do j = S21, N21
           do i = S12, N12
              !-----------------------------------------------------------------------------------------------
              if (comp_visc_yes) then
                 if (dimens == 3) then
                    alpha_re = 1./Re * ((kdx_max_visc(1)*dx1pM(i))**2 +  &
                          &             (kdx_max_visc(2)*dx2pM(j))**2 +  &
                          &             (kdx_max_visc(3)*dx3pM(k))**2)
                 else
                    alpha_re = 1./Re * ((kdx_max_visc(1)*dx1pM(i))**2 +  &
                          &             (kdx_max_visc(2)*dx2pM(j))**2)
                 end if
              else
                 alpha_re = 0.
                 do ii = b1L, b1U
                    alpha_re = alpha_re + ABS(cp11(ii,i))
                 end do
                 do jj = b2L, b2U
                    alpha_re = alpha_re + ABS(cp22(jj,j))
                 end do
                 if (dimens == 3) then
                    do kk = b3L, b3U
                       alpha_re = alpha_re + ABS(cp33(kk,k))
                    end do
                 end if
                 alpha_re = alpha_re / Re
              end if
              !-----------------------------------------------------------------------------------------------
              if (comp_conv_yes) then ! TEST!!! Betraege waeren evtl. gar nicht notwendig!? Argument-Berechnung muesste man evtl. anpassen ...
                 if (dimens == 3) then
                    alpha_im = ABS(work1(i,j,k)*kdx_max_conv(1)*dx1pM(i)) +  &
                          &    ABS(work2(i,j,k)*kdx_max_conv(2)*dx2pM(j)) +  &
                          &    ABS(work3(i,j,k)*kdx_max_conv(3)*dx3pM(k))
                 else
                    alpha_im = ABS(work1(i,j,k)*kdx_max_conv(1)*dx1pM(i)) +  &
                          &    ABS(work2(i,j,k)*kdx_max_conv(2)*dx2pM(j))
                 end if
              else
                 dd1 = 0.
                 do ii = b1L, b1U
                    dd1 = dd1 + ABS(cp1(ii,i))
                 end do
                 alpha_im = dd1*ABS(work1(i,j,k))
                 
                 dd1 = 0.
                 do jj = b2L, b2U
                    dd1 = dd1 + ABS(cp2(jj,j))
                 end do
                 alpha_im = alpha_im + dd1*ABS(work2(i,j,k))
                 
                 if (dimens == 3) then
                    dd1 = 0.
                    do kk = b3L, b3U
                       dd1 = dd1 + ABS(cp3(kk,k))
                    end do
                    alpha_im = alpha_im + dd1*ABS(work3(i,j,k))
                 end if
              end if
              !-----------------------------------------------------------------------------------------------
              if (ABS(alpha_im) < 10.**(-16)) then
                 angle = pi/2.
              else
                 angle = ATAN(alpha_re/alpha_im) ! angle = 0,pi/2
              end if
              
              radius = SQRT(alpha_re**2 + alpha_im**2)
              
              beta(1) = MAX(beta(1),radius  /stabilitylimit(NINT(2.*REAL(n_stab-1)*angle/pi)))
              beta(2) = MAX(beta(2),alpha_im/stabilitylimit(0                               ))
           end do
        end do
     end do
     !========================================================================================================
     if (concentration_yes) then
        do m = 1, n_conc
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
                 do i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    if (comp_visc_yes) then
                       if (dimens == 3) then
                          alpha_re = 1./(Re*Sc(m)) * ((kdx_max_visc(1)*dx1pM(i))**2 +  &
                                          &           (kdx_max_visc(2)*dx2pM(j))**2 +  &
                                          &           (kdx_max_visc(3)*dx3pM(k))**2)
                       else
                          alpha_re = 1./(Re*Sc(m)) * ((kdx_max_visc(1)*dx1pM(i))**2 +  &
                                          &           (kdx_max_visc(2)*dx2pM(j))**2)
                       end if
                    else
                       alpha_re = 0.
                       do ii = b1L, b1U
                          alpha_re = alpha_re + ABS(cc11(ii,i,m))
                       end do
                       do jj = b2L, b2U
                          alpha_re = alpha_re + ABS(cc22(jj,j,m))
                       end do
                       if (dimens == 3) then
                          do kk = b3L, b3U
                             alpha_re = alpha_re + ABS(cc33(kk,k,m))
                          end do
                       end if
                       alpha_re = alpha_re / (Re*Sc(m))
                    end if
                    !-----------------------------------------------------------------------------------------
                    if (comp_conv_yes) then
                       if (dimens == 3) then
                          alpha_im = ABS((work1(i,j,k)+us_vec(1,m))*kdx_max_conv(1)*dx1pM(i)) +  &
                                &    ABS((work2(i,j,k)+us_vec(2,m))*kdx_max_conv(2)*dx2pM(j)) +  &
                                &    ABS((work3(i,j,k)+us_vec(3,m))*kdx_max_conv(3)*dx3pM(k))
                       else
                          alpha_im = ABS((work1(i,j,k)+us_vec(1,m))*kdx_max_conv(1)*dx1pM(i)) +  &
                                &    ABS((work2(i,j,k)+us_vec(2,m))*kdx_max_conv(2)*dx2pM(j))
                       end if
                    else
                       dd1 = 0.
                       do ii = b1L, b1U
                          dd1 = dd1 + ABS(cc1(ii,i,m))
                       end do
                       alpha_im = dd1*ABS(work1(i,j,k)+us_vec(1,m))
                       
                       dd1 = 0.
                       do jj = b2L, b2U
                          dd1 = dd1 + ABS(cc2(jj,j,m))
                       end do
                       alpha_im = alpha_im + dd1*ABS(work2(i,j,k)+us_vec(2,m))
                       
                       if (dimens == 3) then
                          dd1 = 0.
                          do kk = b3L, b3U
                             dd1 = dd1 + ABS(cc3(kk,k,m))
                          end do
                          alpha_im = alpha_im + dd1*ABS(work3(i,j,k)+us_vec(3,m))
                       end if
                    end if
                    !-----------------------------------------------------------------------------------------
                    if (ABS(alpha_im) < 10.**(-16)) then
                       angle = pi/2.
                    else
                       angle = ATAN(alpha_re/alpha_im) ! angle = 0,pi/2
                    end if
                    
                    radius = SQRT(alpha_re**2 + alpha_im**2)
                    
                    beta(1) = MAX(beta(1),radius  /stabilitylimit(NINT(2.*REAL(n_stab-1)*angle/pi)))
                    beta(2) = MAX(beta(2),alpha_im/stabilitylimit(0                               ))
                 end do
              end do
           end do
        end do
     end if
     !========================================================================================================
     
     !CALL MPI_ALLREDUCE(beta,beta_global,2,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     call MPI_REDUCE(beta,beta_global,2,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!! ok?
     
     if (rank == 0) then ! TEST!!! ok?
        if (timeint_mode == 0) then
           dtime = beta_global(2) / CFL ! Definition: CFL = 0,1
        else
           dtime = beta_global(1) / CFL ! Definition: CFL = 0,1
        end if
        
        
        if (ABS(dtime) < 1./dtime_max) then
           dtime = dtime_max
        else
           dtime = 1./ dtime
        end if
        
        if (timestep == 0 .and. dtime > dtime0) dtime = dtime0
     end if
     
  end if
  
  !===========================================================================================================
  
  array_real = 0.
  array_log  = .false.
  
  if (rank == 0) then
     ! sicherstellen, dass Zeitschritte nicht zu klein werden:
     newtime = time+dtime + 10.**(-12)
     
     if (newtime >= time_out_vect .and. dtime_out_vect /= 0.) then
        dtime          = time_out_vect - time
        write_out_vect = .true.
        new_dtime      = .true.
     end if
     
     if (newtime >= time_out_scal .and. dtime_out_scal /= 0.) then
        dtime          = time_out_scal - time
        write_out_scal = .true.
        new_dtime      = .true.
     end if
     
     if (newtime >= time_end) then
        dtime          = time_end - time
        finish_yes     = .true.
     end if
     
     if (timestep+1 >= n_timesteps) then
        finish_yes     = .true.
     end if
     
     array_real = (/time,dtime/)
     array_log  = (/finish_yes,new_dtime,write_out_vect,write_out_scal/)
  end if
  
  call MPI_BCAST(array_real,2,MPI_REAL8  ,0,COMM_CART,merror)
  call MPI_BCAST(array_log ,4,MPI_LOGICAL,0,COMM_CART,merror)
  
  time  = array_real(1)
  dtime = array_real(2)
  
  finish_yes     = array_log(1)
  new_dtime      = array_log(2)
  write_out_vect = array_log(3)
  write_out_scal = array_log(4)
  
  
  if (timestep == 0) dtime_old = dtime
  
  dtime_average = dtime_average + dtime
  
  
  end subroutine get_dtime
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! Noch nicht 100% getestet ...
  subroutine get_beta() bind(c,name='fget_beta')
  
  implicit none
  
  integer                ::  m
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1, beta, beta_global
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
  if (comp_visc_yes) then
     
     ! Bei kompakten Differenzen ist die Analyse nicht so einfach ...
     if (rank == 0) write(*,'(a)') 'Note: Cannot determine beta (compact differences)!'
     
  else
     
     if (rank == 0) write(*,'(a)') 'determining beta ...'
     
     beta = 0.
     
     !========================================================================================================
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              dd1 = 0.
              do ii = b1L, b1U
                 dd1 = dd1 + ABS(cu11(ii,i))
              end do
              do jj = b2L, b2U
                 dd1 = dd1 + ABS(cp22(jj,j))
              end do
              if (dimens == 3) then
                 do kk = b3L, b3U
                    dd1 = dd1 + ABS(cp33(kk,k))
                 end do
              end if
              beta = MAX(beta,dd1)
           end do
        end do
     end do
     
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              dd1 = 0.
              do ii = b1L, b1U
                 dd1 = dd1 + ABS(cp11(ii,i))
              end do
              do jj = b2L, b2U
                 dd1 = dd1 + ABS(cv22(jj,j))
              end do
              if (dimens == 3) then
                 do kk = b3L, b3U
                    dd1 = dd1 + ABS(cp33(kk,k))
                 end do
              end if
              beta = MAX(beta,dd1)
           end do
        end do
     end do
     
     if (dimens == 3) then
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = 0.
                 do ii = b1L, b1U
                    dd1 = dd1 + ABS(cp11(ii,i))
                 end do
                 do jj = b2L, b2U
                    dd1 = dd1 + ABS(cp22(jj,j))
                 end do
                 do kk = b3L, b3U
                    dd1 = dd1 + ABS(cw33(kk,k))
                 end do
                 beta = MAX(beta,dd1)
              end do
           end do
        end do
     end if
     
     if (concentration_yes) then
        do m = 1, n_conc
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
                 do i = S1c(m), N1c(m)
                    dd1 = 0.
                    do ii = b1L, b1U
                       dd1 = dd1 + ABS(cc11(ii,i,m))
                    end do
                    do jj = b2L, b2U
                       dd1 = dd1 + ABS(cc22(jj,j,m))
                    end do
                    if (dimens == 3) then
                       do kk = b3L, b3U
                          dd1 = dd1 + ABS(cc33(kk,k,m))
                       end do
                    end if
                    beta = MAX(beta,dd1/Sc(m))
                 end do
              end do
           end do
        end do
     end if
     !========================================================================================================
     
     beta = beta / Re
     
     !CALL MPI_ALLREDUCE(beta,beta_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     call MPI_REDUCE(beta,beta_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     if (rank == 0) then
        open(10,file='test_beta.txt',status='UNKNOWN')
        write(10,'(a,E25.17)') '||L||_infty  = ', beta_global
        close(10)
     end if
     
  end if
  
  
  end subroutine get_beta
  
  
  
  
  
  
  
  
  
  
  
  subroutine init_BC ! TEST!!! Name nicht ganz passend ...
  
  implicit none
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  integer                ::  m
  
  
  !bc11 = 0. ! TEST!!! nach usr_initcond.f90 bzw. usr_boundcond.f90 verschoben ...
  !bc12 = 0.
  !bc13 = 0.
  !
  !bc21 = 0.
  !bc22 = 0.
  !bc23 = 0.
  !
  !bc31 = 0.
  !bc32 = 0.
  !bc33 = 0.
  
  
  !===========================================================================================================
  !=== Ausfluss-RB initialisieren ============================================================================
  !===========================================================================================================
  if (outlet(1,1,2) .and. (BC_1L == 1 .or. BC_1L == 3)) bc21(S22B:N22B,S32B:N32B,1) = vel(1, S22B:N22B,S32B:N32B,2)
  if (outlet(1,2,2) .and. (BC_1U == 1 .or. BC_1U == 3)) bc21(S22B:N22B,S32B:N32B,2) = vel(N1,S22B:N22B,S32B:N32B,2)
  
  if (outlet(1,1,3) .and. (BC_1L == 1 .or. BC_1L == 3)) bc31(S23B:N23B,S33B:N33B,1) = vel(1 ,S23B:N23B,S33B:N33B,3)
  if (outlet(1,2,3) .and. (BC_1U == 1 .or. BC_1U == 3)) bc31(S23B:N23B,S33B:N33B,2) = vel(N1,S23B:N23B,S33B:N33B,3)
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(2,1,1) .and. (BC_2L == 1 .or. BC_2L == 3)) bc12(S11B:N11B,S31B:N31B,1) = vel(S11B:N11B,1 ,S31B:N31B,1)
  if (outlet(2,2,1) .and. (BC_2U == 1 .or. BC_2U == 3)) bc12(S11B:N11B,S31B:N31B,2) = vel(S11B:N11B,N2,S31B:N31B,1)
  
  if (outlet(2,1,3) .and. (BC_2L == 1 .or. BC_2L == 3)) bc32(S13B:N13B,S33B:N33B,1) = vel(S13B:N13B,1 ,S33B:N33B,3)
  if (outlet(2,2,3) .and. (BC_2U == 1 .or. BC_2U == 3)) bc32(S13B:N13B,S33B:N33B,2) = vel(S13B:N13B,N2,S33B:N33B,3)
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,1) .and. (BC_3L == 1 .or. BC_3L == 3)) bc13(S11B:N11B,S21B:N21B,1) = vel(S11B:N11B,S21B:N21B,1 ,1)
  if (outlet(3,2,1) .and. (BC_3U == 1 .or. BC_3U == 3)) bc13(S11B:N11B,S21B:N21B,2) = vel(S11B:N11B,S21B:N21B,N3,1)
  
  if (outlet(3,1,2) .and. (BC_3L == 1 .or. BC_3L == 3)) bc23(S12B:N12B,S22B:N22B,1) = vel(S12B:N12B,S22B:N22B,1 ,2)
  if (outlet(3,2,2) .and. (BC_3U == 1 .or. BC_3U == 3)) bc23(S12B:N12B,S22B:N22B,2) = vel(S12B:N12B,S22B:N22B,N3,2)
  !===========================================================================================================
  if (outlet(1,1,1) .and. (BC_1L == 1 .or. BC_1L == 2)) then
     i = 1
     do k = S31B, N31B
        do j = S21B, N21B
           bc11(j,k,1) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           do ii = d1L+1, d1U
              bc11(j,k,1) = bc11(j,k,1) + cIup(ii,i)*vel(i+ii,j,k,1)
           end do
        end do
     end do
  end if
  if (outlet(1,2,1) .and. (BC_1U == 1 .or. BC_1U == 2)) then
     i = N1
     do k = S31B, N31B
        do j = S21B, N21B
           bc11(j,k,2) = cIup(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           do ii = d1L+1, d1U
              bc11(j,k,2) = bc11(j,k,2) + cIup(ii,i)*vel(i+ii,j,k,1)
           end do
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(2,1,2) .and. (BC_2L == 1 .or. BC_2L == 2)) then
     j = 1
     do k = S32B, N32B
        do i = S12B, N12B
           bc22(i,k,1) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           do jj = d2L+1, d2U
              bc22(i,k,1) = bc22(i,k,1) + cIvp(jj,j)*vel(i,j+jj,k,2)
           end do
        end do
     end do
  end if
  if (outlet(2,2,2) .and. (BC_2U == 1 .or. BC_2U == 2)) then
     j = N2
     do k = S32B, N32B
        do i = S12B, N12B
           bc22(i,k,2) = cIvp(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           do jj = d2L+1, d2U
              bc22(i,k,2) = bc22(i,k,2) + cIvp(jj,j)*vel(i,j+jj,k,2)
           end do
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,3) .and. (BC_3L == 1 .or. BC_3L == 2)) then
     k = 1
     do j = S23B, N23B
        do i = S13B, N13B
           bc33(i,j,1) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           do kk = d3L+1, d3U
              bc33(i,j,1) = bc33(i,j,1) + cIwp(kk,k)*vel(i,j,k+kk,3)
           end do
        end do
     end do
  end if
  if (outlet(3,2,3) .and. (BC_3U == 1 .or. BC_3U == 2)) then
     k = N3
     do j = S23B, N23B
        do i = S13B, N13B
           bc33(i,j,2) = cIwp(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           do kk = d3L+1, d3U
              bc33(i,j,2) = bc33(i,j,2) + cIwp(kk,k)*vel(i,j,k+kk,3)
           end do
        end do
     end do
  end if
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !=== Konzentrationsfeld ====================================================================================
  !===========================================================================================================
  if (concentration_yes) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- nicht-Durchfluss-RB berücksichtigen ----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     do m = 1, n_conc
        
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1L(m) == 3) then
           i = 1
           do k = S3p, N3p
              do j = S2p, N2p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do ii = 1, b1U
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i+ii,j,k,m) * cc1(ii,i,m) / (cc1(0,i,m) - usReSc(1,m) - bc11(j,k,1))
                 end do
              end do
           end do
        end if
        
        if (BCc_1U(m) == 3) then
           i = N1
           do k = S3p, N3p
              do j = S2p, N2p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do ii = b1L, -1
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i+ii,j,k,m) * cc1(ii,i,m) / (cc1(0,i,m) - usReSc(1,m) - bc11(j,k,2))
                 end do
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
           j = 1
           do k = S3p, N3p
              do i = S1p, N1p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do jj = 1, b2U
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i,j+jj,k,m) * cc2(jj,j,m) / (cc2(0,j,m) - usReSc(2,m) - bc22(i,k,1))
                 end do
              end do
           end do
        end if
        
        if (BCc_2U(m) == 3) then
           j = N2
           do k = S3p, N3p
              do i = S1p, N1p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do jj = b2L, -1
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i,j+jj,k,m) * cc2(jj,j,m) / (cc2(0,j,m) - usReSc(2,m) - bc22(i,k,2))
                 end do
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
           k = 1
           do j = S2p, N2p
              do i = S1p, N1p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do kk = 1, b3U
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i,j,k+kk,m) * cc3(kk,k,m) / (cc3(0,k,m) - usReSc(3,m) - bc33(i,j,1))
                 end do
              end do
           end do
        end if
        
        if (BCc_3U(m) == 3) then
           k = N3
           do j = S2p, N2p
              do i = S1p, N1p
                 conc(i,j,k,m) = 0.
!pgi$ unroll = n:8
                 do kk = b3L, -1
                    conc(i,j,k,m) = conc(i,j,k,m) - conc(i,j,k+kk,m) * cc3(kk,k,m) / (cc3(0,k,m) - usReSc(3,m) - bc33(i,j,2))
                 end do
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        
     end do
     !--------------------------------------------------------------------------------------------------------
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- Ausfluss-RB initialisieren (für Statistik) ---------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     conc1L(S2p:N2p,S3p:N3p,1:n_conc) = conc(S1p,S2p:N2p,S3p:N3p,1:n_conc)
     conc1U(S2p:N2p,S3p:N3p,1:n_conc) = conc(N1p,S2p:N2p,S3p:N3p,1:n_conc)
     
     conc2L(S1p:N1p,S3p:N3p,1:n_conc) = conc(S1p:N1p,S2p,S3p:N3p,1:n_conc)
     conc2U(S1p:N1p,S3p:N3p,1:n_conc) = conc(S1p:N1p,N2p,S3p:N3p,1:n_conc)
     
     conc3L(S1p:N1p,S2p:N2p,1:n_conc) = conc(S1p:N1p,S2p:N2p,S3p,1:n_conc)
     conc3U(S1p:N1p,S2p:N2p,1:n_conc) = conc(S1p:N1p,S2p:N2p,N3p,1:n_conc)
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  
  
  
  end subroutine init_BC
  
  
  
  
  
  
  
  
  
  
  
  subroutine level_pressure
  
  implicit none
  
  integer                ::  i, j, k
  real                   ::  pre0, pre0_global
  
  
  if ((Int_lev_pre /= 0 .and. MOD(timestep,Int_lev_pre) == 0) .or. write_out_vect) then
     
     pre0 = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              pre0 = pre0 + pre(i,j,k)
           end do
        end do
     end do
     
     call MPI_ALLREDUCE(pre0,pre0_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
     
     pre0 = pre0_global/REAL(dim1)/REAL(dim2)/REAL(dim3) ! TEST!!! wegen i4 gefaehrlich!
     
     pre(S1p:N1p,S2p:N2p,S3p:N3p) = pre(S1p:N1p,S2p:N2p,S3p:N3p) - pre0
     
  end if
  
  
  end subroutine level_pressure
  
  
  
  
  
  
  
  
  
  
  
  !> fills corners/edges for pressure/concentration, with
  subroutine fill_corners(phi)
  
  implicit none
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei, so dass der Druck hier unbestimmt ist.
  ! Hier Mittelung der Nachbarpunkte
  
  if (BC_1L > 0 .and. BC_2L > 0) phi(1 ,1 ,1:N3) = (phi(2   ,1 ,1:N3) + phi(1 ,2   ,1:N3) + phi(2   ,2   ,1:N3))/3.
  if (BC_1L > 0 .and. BC_2U > 0) phi(1 ,N2,1:N3) = (phi(2   ,N2,1:N3) + phi(1 ,N2-1,1:N3) + phi(2   ,N2-1,1:N3))/3.
  if (BC_1U > 0 .and. BC_2L > 0) phi(N1,1 ,1:N3) = (phi(N1-1,1 ,1:N3) + phi(N1,2   ,1:N3) + phi(N1-1,2   ,1:N3))/3.
  if (BC_1U > 0 .and. BC_2U > 0) phi(N1,N2,1:N3) = (phi(N1-1,N2,1:N3) + phi(N1,N2-1,1:N3) + phi(N1-1,N2-1,1:N3))/3.

  if (BC_1L > 0 .and. BC_3L > 0) phi(1 ,1:N2,1 ) = (phi(2   ,1:N2,1 ) + phi(1 ,1:N2,2   ) + phi(2   ,1:N2,2   ))/3.
  if (BC_1L > 0 .and. BC_3U > 0) phi(1 ,1:N2,N3) = (phi(2   ,1:N2,N3) + phi(1 ,1:N2,N3-1) + phi(2   ,1:N2,N3-1))/3.
  if (BC_1U > 0 .and. BC_3L > 0) phi(N1,1:N2,1 ) = (phi(N1-1,1:N2,1 ) + phi(N1,1:N2,2   ) + phi(N1-1,1:N2,2   ))/3.
  if (BC_1U > 0 .and. BC_3U > 0) phi(N1,1:N2,N3) = (phi(N1-1,1:N2,N3) + phi(N1,1:N2,N3-1) + phi(N1-1,1:N2,N3-1))/3.

  if (BC_2L > 0 .and. BC_3L > 0) phi(1:N1,1 ,1 ) = (phi(1:N1,2   ,1 ) + phi(1:N1,1 ,2   ) + phi(1:N1,2   ,2   ))/3.
  if (BC_2L > 0 .and. BC_3U > 0) phi(1:N1,1 ,N3) = (phi(1:N1,2   ,N3) + phi(1:N1,1 ,N3-1) + phi(1:N1,2   ,N3-1))/3.
  if (BC_2U > 0 .and. BC_3L > 0) phi(1:N1,N2,1 ) = (phi(1:N1,N2-1,1 ) + phi(1:N1,N2,2   ) + phi(1:N1,N2-1,2   ))/3.
  if (BC_2U > 0 .and. BC_3U > 0) phi(1:N1,N2,N3) = (phi(1:N1,N2-1,N3) + phi(1:N1,N2,N3-1) + phi(1:N1,N2-1,N3-1))/3.
  
  
  end subroutine fill_corners
  
  
  
  
  
  
  
  
  
  
  
  subroutine init_alarm()  bind(c, name='finit_alarm')
  
  implicit none
  
  integer                ::  wtime, margin ! TEST!!! margin nach jobscript verlagert!
  integer                ::  ios
  
  
  !*************
  margin = 0 ! [seconds] ! TEST!!! margin nach jobscript verlagert!
  !*************
  
  
  if (rank == 0) then
     open(10,file='queue.txt',action='read',status='old',iostat=ios)
     
     if (ios == 0) then
        read (10,*) wtime
        close(10)
        
        if (wtime > margin) wtime = wtime - margin ! TEST!!! margin nach jobscript verlagert!
        
        write(*,*)
        write(*,'(a,i6,a)') '... found queue.txt, exiting after', wtime, ' seconds ...'
        write(*,*)
        
        call ALARM(wtime,check_alarm)
     else
        write(*,*)
        write(*,*) '... found no queue.txt ...'
        write(*,*)
     end if
  end if
  
  
  end subroutine init_alarm
  
  
  
  
  
  
  
  
  
  
  
  subroutine check_alarm
  
  implicit none
  
  
  if (rank == 0) then
     
     write(*,*) '+++++++++++++++++++++++++++++++++++'
     write(*,*) '+ ... received the alarm sign ... +'
     write(*,*) '+++++++++++++++++++++++++++++++++++'
     
     finish_yes = .true.
     
  end if
  
  
  end subroutine check_alarm
  
  
  
  
  
  
  
  
  
  
  
  subroutine check_signal
  
  implicit none
  
  integer                ::  signal
  integer                ::  timestep_end
  real                   ::  time_end
  
  
  if (rank == 0) then
     
     open (10, file='send_signal.txt', status='UNKNOWN')
     read (10,*) signal
     read (10,*) timestep_end
     read (10,*) time_end
     close(10)
     
     if (signal == 1) finish_yes = .true.
     if (signal == 2 .and. (timestep == timestep_end .or. time >= time_end)) finish_yes = .true.
     
  end if
  
  call MPI_BCAST(finish_yes,1,MPI_LOGICAL,0,COMM_CART,merror)
  
  
  end subroutine check_signal
  
  
  
  
  
  
  
  
  
  
  
  subroutine num_to_string(n_digits,num,num_char)
  
  implicit none
  
  integer     , intent(in ) ::  n_digits
  integer     , intent(in ) ::  num
  character(len=n_digits), intent(out) ::  num_char
  
  integer                   ::  tenthousands, thousands, hundreds, tens, ones
  
  
  !===========================================================================================================
  if (num < 0) then
     if (rank == 0) write(*,*) 'ERROR! Cannot convert negative integers to a string.'
     call MPI_FINALIZE(merror)
     stop
  !===========================================================================================================
  else if (n_digits >= 6) then
     if (rank == 0) write(*,*) 'ERROR! Cannot convert integers > 99999 to a string.'
     call MPI_FINALIZE(merror)
     stop
  !===========================================================================================================
  else if (n_digits == 1) then
     write(num_char,'(i1.1)') num
  !===========================================================================================================
  else if (n_digits == 2) then
     write(num_char,'(i2.2)') num
  !===========================================================================================================
  else if (n_digits == 3) then
     write(num_char,'(i3.3)') num
  !===========================================================================================================
  else if (n_digits == 4) then
     write(num_char,'(i4.4)') num
  !===========================================================================================================
  else if (n_digits == 5) then
     !tenthousands =  num / 10000
     !thousands    = (num - 10000 * tenthousands) / 1000
     !hundreds     = (num - 10000 * tenthousands - 1000 * thousands) / 100
     !tens         = (num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds) / 10
     !ones         =  num - 10000 * tenthousands - 1000 * thousands - 100 * hundreds - tens*10
     !
     !num_char(1:1) = CHAR(ICHAR('0')+tenthousands)
     !num_char(2:2) = CHAR(ICHAR('0')+thousands)
     !num_char(3:3) = CHAR(ICHAR('0')+hundreds)
     !num_char(4:4) = CHAR(ICHAR('0')+tens)
     !num_char(5:5) = CHAR(ICHAR('0')+ones)
     write(num_char,'(i5.5)') num
  end if
  !===========================================================================================================
  
  
  end subroutine num_to_string
  
  
  
  
  
  
  
  
  
  
  
  subroutine iteration_stats
  
  ! revised: 24.10.07
  
  implicit none
  
  real                ::  ratioO_tot
  real                ::  ratioH_tot
  real                ::  ratioP_tot
  
  integer             ::  countO_tot
  integer             ::  countP_tot
  integer             ::  countH_tot
  
  real                ::  ratioP_all
  real                ::  ratioH_all
  
  integer             ::  countP_all
  integer             ::  countH_all
  
  real                ::  countO_av
  real                ::  countP_av(1:2)
  real                ::  countH_av(1:3)
  real                ::  countP_all_av
  real                ::  countH_all_av
  
  integer             ::  tsteps
  
  
  tsteps = timestep - timestep_old
  
  
  if (rank == 0 .and. tsteps >= 1) then
     open(10, file='test_iteration_stats_restart'//restart_char//'.txt', status='UNKNOWN')
     
     dtime_average = dtime_average / REAL(tsteps)
     
     countO_tot = 0
     countP_tot = 0
     countH_tot = 0
     
     ratioO_tot = 0.
     ratioP_tot = 0.
     ratioH_tot = 0.
     
     write(10,'(a       )') ''
     write(10,'(a, G13.5)') '                        Re =', Re
     !WRITE(10,'(a, G13.5)') '                      freq =', freq
     write(10,'(a,3E13.5)') '                Dimensions =', L1 , L2 , L3
     write(10,'(a,3i5   )') '                Resolution =', M1 , M2 , M3
     write(10,'(a,3i5   )') '                     Procs =', NB1, NB2, NB3
     write(10,'(a, E13.5)') '                   <dtime> =', dtime_average
     write(10,'(a, i7   )') '               n_timesteps =', tsteps
     write(10,'(a, E13.5)') '        elapsed time [sec] =', REAL(elatime)/1000.
     write(10,'(a, E13.5)') 'el. time/n_timesteps [sec] =', REAL(elatime)/1000./REAL(tsteps)
     write(10,'(a,2E13.5)') '               div(vel(0)) =', max_div_init(1:2)
     write(10,'(a       )') ''
     write(10,'(a       )') 'number of iterations (number/timestep, <conv. ratio>):'
     write(10,'(a       )') ''
     
     do substep = 1, RK_steps
        
        countP_all = countP(substep,1) + countP(substep,2)
        countH_all = countH(substep,1) + countH(substep,2) + countH(substep,3)
        
        ratioP_all = ratioP(substep,1) + ratioP(substep,2)
        ratioH_all = ratioH(substep,1) + ratioH(substep,2) + ratioH(substep,3)
        
        !-----------------------------------------------------------------------------------------------------
        
        countO_tot = countO_tot + countO(substep)
        countP_tot = countP_tot + countP_all
        countH_tot = countH_tot + countH_all
        
        ratioO_tot = ratioO_tot + ratioO(substep)
        ratioP_tot = ratioP_tot + ratioP_all
        ratioH_tot = ratioH_tot + ratioH_all
        
        !-----------------------------------------------------------------------------------------------------
        
        countO_av      = REAL(countO(substep))     / REAL(tsteps)
        countP_av(1:2) = REAL(countP(substep,1:2)) / REAL(tsteps)
        countH_av(1:3) = REAL(countH(substep,1:3)) / REAL(tsteps)
        countP_all_av  = REAL(countP_all)          / REAL(tsteps)
        countH_all_av  = REAL(countH_all)          / REAL(tsteps)
        
        
        if (countO(substep)   /= 0) ratioO(substep)   = ratioO(substep)   / REAL(countO(substep)  )
        
        if (countP(substep,1) /= 0) ratioP(substep,1) = ratioP(substep,1) / REAL(countP(substep,1))
        if (countP(substep,2) /= 0) ratioP(substep,2) = ratioP(substep,2) / REAL(countP(substep,2))
        
        if (countH(substep,1) /= 0) ratioH(substep,1) = ratioH(substep,1) / REAL(countH(substep,1))
        if (countH(substep,2) /= 0) ratioH(substep,2) = ratioH(substep,2) / REAL(countH(substep,2))
        if (countH(substep,3) /= 0) ratioH(substep,3) = ratioH(substep,3) / REAL(countH(substep,3))
        
        if (countP_all        /= 0) ratioP_all        = ratioP_all        / REAL(countP_all       )
        if (countH_all        /= 0) ratioH_all        = ratioH_all        / REAL(countH_all       )
        
        !-----------------------------------------------------------------------------------------------------
        
        write(10,'(a,i2,a  )') '   substep', substep, ':'
        write(10,'(a       )') '      '
        write(10,'(a,2G12.4)') '      outer Iterations:', countO_av    , ratioO    (substep)
        write(10,'(a       )') '         '
        write(10,'(a,2G12.4)') '         Poisson   (1):', countP_av(1) , ratioP    (substep,1)
        write(10,'(a,2G12.4)') '         Poisson   (2):', countP_av(2) , ratioP    (substep,2)
        write(10,'(a,2G12.4)') '         all          :', countP_all_av, ratioP_all
        write(10,'(a       )') '         '
        write(10,'(a,2G12.4)') '         Helmholtz (1):', countH_av(1) , ratioH    (substep,1)
        write(10,'(a,2G12.4)') '         Helmholtz (2):', countH_av(2) , ratioH    (substep,2)
        write(10,'(a,2G12.4)') '         Helmholtz (3):', countH_av(3) , ratioH    (substep,3)
        write(10,'(a,2G12.4)') '         all          :', countH_all_av, ratioH_all
        write(10,'(a       )') '         '
        
     end do
     
     if (countO_tot /= 0) ratioO_tot = ratioO_tot / REAL(countO_tot)
     if (countP_tot /= 0) ratioP_tot = ratioP_tot / REAL(countP_tot)
     if (countH_tot /= 0) ratioH_tot = ratioH_tot / REAL(countH_tot)
     
     write(10,'(a       )') '   total:'
     write(10,'(a       )') '      '
     write(10,'(a,2G12.4)') '      outer Iterations:', REAL(countO_tot)/REAL(tsteps), ratioO_tot
     write(10,'(a,2G12.4)') '         Poisson      :', REAL(countP_tot)/REAL(tsteps), ratioP_tot
     write(10,'(a,2G12.4)') '         Helmholtz    :', REAL(countH_tot)/REAL(tsteps), ratioH_tot
     write(10,'(a       )') '         '
     
     close(10)
     
  end if
  
  end subroutine iteration_stats
  
  
  
  
end module mod_lib
