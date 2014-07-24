!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module cmod_test

  use iso_c_binding
  use mpi
  
  
  use mod_dims
  use mod_vars
  use mod_exchange
  use mod_diff      ! CALL divergence (analyze_matrix), CALL gradient (analyze_matrix)
  use mod_helmholtz ! CALL Helmholtz (analyze_matrix)
  use mod_laplace   ! CALL product_div_grad (analyze_matrix)
  
  
  private
  
  public test_parameter
  public test_divergence, test_momentum
  public analyze_matrix
  
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine test_parameter() bind(c,name='ftest_parameterlight')
  
  ! revised: 26.05.08
  
  implicit none
  
  logical                ::  error(1:4)
  integer                ::  m
  integer                ::  comm_size
  
  
  !*******************************************************************************
!  if (M3 == 2) then ! TEST!!! wird auch in init_general nochmals ausgefuehrt!
!     dimens = 2
!  else
!     dimens = 3
!  end if
  
  if (M3 == 2) then ! TEST!!! Notwendig? Pruefen, ggf. aendern und diese Fallunterscheidung loeschen ...
     BC_3L_global = -1
     BC_3U_global = -1
     
     impl_dir(3) = 0
     
     outlet   (3,:,:) = .false.
     isopycnal(3,:,:) = .false.
     
     gravity(3) = 0.
  end if
  !*******************************************************************************
  
  error(1:4) = .false.
  
  
  ! TEST!!! - Testroutine schreiben für dtime_out_scal .LT. 0. & dtime_out_vect .LT. 0. ...
  !         - noch einen Test schreiben um ndL <= n1U zu testen ...
  !         - dimXX > 1 ebenfalls testen ...
  !         - ndR =< b1U, ndR =< b2U, ndR =< b3U
  if (rank == 0) then
     
     !========================================================================================================
     !=== Dimensionen ========================================================================================
     !========================================================================================================
     ! TEST!!! unschoen ...
!#ifdef ALLOC
!     N1 = 1+(M1-1)/NB1
!     N2 = 1+(M2-1)/NB2
!     N3 = 1+(M3-1)/NB3
!#endif
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 2 .and. N3 /= 2) then
        write(*,*) 'ERROR! Choose N3 = 2 for 2D!'
        error(1) = .true.
     end if
     !--------------------------------------------------------------------------------------------------------
     if (n_conc < 1) then
        write(*,'(a)') 'ERROR! n_conc < 1'
        error(1) = .true.
     end if
     !--------------------------------------------------------------------------------------------------------
!     if (N1 < 3) then
!        write(*,*) 'ERROR! N1 < 3!'
!        error(1) = .true.
!     end if
!     if (N2 < 3) then
!        write(*,*) 'ERROR! N2 < 3!'
!        error(1) = .true.
!     end if
     if (N3 < 3 .and. dimens == 3) then
        write(*,*) 'ERROR! N3 < 3!'
        error(1) = .true.
     end if
     !--------------------------------------------------------------------------------------------------------
!     if (MOD(N1-1,2) /= 0) then
!        write(*,'(a)') 'ERROR! Dimension N1 cannot be used for multigrid!'
!        error(1) = .true.
!     end if
!     if (MOD(N2-1,2) /= 0) then
!        write(*,'(a)') 'ERROR! Dimension N2 cannot be used for multigrid!'
!        error(1) = .true.
!     end if
     if (MOD(N3-1,2) /= 0 .and. dimens == 3) then
        write(*,'(a)') 'ERROR! Dimension N3 cannot be used for multigrid!'
        error(1) = .true.
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Bockaufteilung =====================================================================================
     !========================================================================================================
     if (dimens == 2 .and. NB3 /= 1) then
        write(*,*) 'ERROR! Choose NB3 = 1 for 2D!'
        error(1) = .true.
     end if
     !--------------------------------------------------------------------------------------------------------
     if (ls1 /= 0 .and. ls1 /= -1) then
        write(*,*) 'ERROR! ls1 /= 0,-1!'
        error(1) = .true.
     end if
     if (ls2 /= 0 .and. ls2 /= -1) then
        write(*,*) 'ERROR! ls2 /= 0,-1!'
        error(1) = .true.
     end if
     if (ls3 /= 0 .and. ls3 /= -1) then
        write(*,*) 'ERROR! ls3 /= 0,-1!'
        error(1) = .true.
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== FD-Stencil-Dimensionen =============================================================================
     !========================================================================================================
     if (b1L > -1) then
        write(*,*) 'ERROR! b1L > -1!'
        error(1) = .true.
     end if
     if (b2L > -1) then
        write(*,*) 'ERROR! b2L > -1!'
        error(1) = .true.
     end if
     if (b3L > -1) then
        write(*,*) 'ERROR! b3L > -1!'
        error(1) = .true.
     end if
     
     if (b1U < 1) then
        write(*,*) 'ERROR! b1U < 1!'
        error(1) = .true.
     end if
     if (b2U < 1) then
        write(*,*) 'ERROR! b2U < 1!'
        error(1) = .true.
     end if
     if (b3U < 1) then
        write(*,*) 'ERROR! b3U < 1!'
        error(1) = .true.
     end if
     !--------------------------------------------------------------------------------------------------------
     if ((ndR > b1U .or. ndR > -b1L) .and. (comp_visc_yes .or. comp_conv_yes .or. comp_inter_yes .or. comp_div_yes .or. comp_grad_yes)) then ! TEST!!!
        write(*,*) 'ERROR! ndR > b1U or ndR > -b1L!'
        error(1) = .true.
     end if
     if ((ndR > b2U .or. ndR > -b2L) .and. (comp_visc_yes .or. comp_conv_yes .or. comp_inter_yes .or. comp_div_yes .or. comp_grad_yes)) then
        write(*,*) 'ERROR! ndR > b2U or ndR > -b2L!'
        error(1) = .true.
     end if
     if ((ndR > b3U .or. ndR > -b3L) .and. (comp_visc_yes .or. comp_conv_yes .or. comp_inter_yes .or. comp_div_yes .or. comp_grad_yes)) then
        write(*,*) 'ERROR! ndR > b3U or ndR > -b3L!'
        error(1) = .true.
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Periodische Raender prüfen =========================================================================
     !========================================================================================================
     ! 1L/1U-Rand:
     if ((BC_1L_global == -1 .and. BC_1U_global /= -1) .or. (BC_1L_global /= -1 .and. BC_1U_global == -1)) error(2) = .true.
     if ((BC_2L_global == -1 .and. BC_2U_global /= -1) .or. (BC_2L_global /= -1 .and. BC_2U_global == -1)) error(2) = .true.
     if ((BC_3L_global == -1 .and. BC_3U_global /= -1) .or. (BC_3L_global /= -1 .and. BC_3U_global == -1)) error(2) = .true.
     
     if (error(2)) write(*,*) 'ERROR! A periodic boundary cannot be combined with a boundary condition!'
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Ausflussrand prüfen ================================================================================
     !========================================================================================================
     ! 1L/1U-Rand:
     if (outlet(1,1,1) .and. BC_1L_global /= 1 .and. BC_1L_global /= 3) error(3) = .true.
     if (outlet(1,2,1) .and. BC_1U_global /= 1 .and. BC_1U_global /= 3) error(3) = .true.
     if (outlet(1,1,2) .and. BC_1L_global /= 1 .and. BC_1L_global /= 2) error(3) = .true.
     if (outlet(1,2,2) .and. BC_1U_global /= 1 .and. BC_1U_global /= 2) error(3) = .true.
     if (outlet(1,1,3) .and. BC_1L_global /= 1 .and. BC_1L_global /= 3) error(3) = .true.
     if (outlet(1,2,3) .and. BC_1U_global /= 1 .and. BC_1U_global /= 3) error(3) = .true.
     
     ! 2L/2U-Rand:
     if (outlet(2,1,1) .and. BC_2L_global /= 1 .and. BC_2L_global /= 2) error(3) = .true.
     if (outlet(2,2,1) .and. BC_2U_global /= 1 .and. BC_2U_global /= 2) error(3) = .true.
     if (outlet(2,1,2) .and. BC_2L_global /= 1 .and. BC_2L_global /= 3) error(3) = .true.
     if (outlet(2,2,2) .and. BC_2U_global /= 1 .and. BC_2U_global /= 3) error(3) = .true.
     if (outlet(2,1,3) .and. BC_2L_global /= 1 .and. BC_2L_global /= 3) error(3) = .true.
     if (outlet(2,2,3) .and. BC_2U_global /= 1 .and. BC_2U_global /= 3) error(3) = .true.
     
     ! 3L/3U-Rand:
     if (outlet(3,1,1) .and. BC_3L_global /= 1 .and. BC_3L_global /= 3) error(3) = .true.
     if (outlet(3,2,1) .and. BC_3U_global /= 1 .and. BC_3U_global /= 3) error(3) = .true.
     if (outlet(3,1,2) .and. BC_3L_global /= 1 .and. BC_3L_global /= 3) error(3) = .true.
     if (outlet(3,2,2) .and. BC_3U_global /= 1 .and. BC_3U_global /= 3) error(3) = .true.
     if (outlet(3,1,3) .and. BC_3L_global /= 1 .and. BC_3L_global /= 2) error(3) = .true.
     if (outlet(3,2,3) .and. BC_3U_global /= 1 .and. BC_3U_global /= 2) error(3) = .true.
     
     if (error(3)) write(*,*) 'ERROR! Choice of outlet configuartion is not suitable!'
     !========================================================================================================
     
     
     !========================================================================================================
     if (concentration_yes .and. n_conc >= 1) then
        do m = 1, n_conc
           if (BC_1L == -2 .and. us_vec(1,m) < 0.) error(4) = .true.
           if (BC_1U == -2 .and. us_vec(1,m) > 0.) error(4) = .true.
           if (BC_2L == -2 .and. us_vec(2,m) < 0.) error(4) = .true.
           if (BC_2U == -2 .and. us_vec(2,m) > 0.) error(4) = .true.
           if (BC_3L == -2 .and. us_vec(3,m) < 0.) error(4) = .true.
           if (BC_3U == -2 .and. us_vec(3,m) > 0.) error(4) = .true.
        end do
        if (error(4)) write(*,*) 'ERROR! Symmetry boundary conditions and settling velocity are contradictory!'
     end if
     !========================================================================================================
     
     
     !========================================================================================================
     !=== Sonstiges prüfen ===================================================================================
     !========================================================================================================
     if (time_start > time_end) then 
        write(*,*) 'ERROR! time_start > time_end!'
        error(1) = .true.
     end if
     
#ifdef NONBOUSSINESQ
     if (nullspace_yes .and. nonBoussinesq_yes) then
        write(*,*) 'ERROR! Flux correction not applicable in non-Boussinesq simulations!'
        error(1) = .true.
     end if
     
     if (timeint_mode /= 1 .and. nonBoussinesq_yes) then
        write(*,*) 'ERROR! Non-Boussinesq simulations cannot be performed with (semi-)implicit time integration!'
        error(1) = .true.
     end if
#endif
     !========================================================================================================
     
  end if
  
  
  call MPI_BCAST(error(1:4),4,MPI_LOGICAL,0,MPI_COMM_WORLD,merror) ! MPI_COMM_CART ist hier noch nicht initialisiert ...
  
  if (error(1) .or. error(2) .or. error(3) .or. error(4)) then
     if (rank == 0) write(*,*) 'Exiting ...'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  
  !===========================================================================================================
  !=== Zeitintegration =======================================================================================
  !===========================================================================================================
  ! Euler_yes ==> thetaL ==> timeint_mode ==> ...
  
  if (thetaL == 0. .and. timeint_mode == 0) then
     write(*,*) 'WARNING! Setting timeint_mode = 1 for thetaL == 0. ...'
     timeint_mode = 1
  end if
  
  if (Euler_yes .and. timeint_mode == 0) then
     write(*,*) 'WARNING! Setting timeint_mode = 1 for Euler_yes == .TRUE. ...'
     timeint_mode = 1
  end if
  !===========================================================================================================
  
  
  end subroutine test_parameter
  
  
  
  
  
  
  
  
  
  
  
!  subroutine test_divergence
!
!  implicit none
!
!  integer                ::  i, j, k, m
!  real                   ::  max_div(1:2), max_div_global(1:2)
!
!
!  ! TEST!!! handle_corner_rhs(res) evtl. zu Divergenz hineinziehen?
!  call divergence2(vel,res)
!
!  !IF (corner_yes) CALL handle_corner_rhs(res) ! TEST!!! mod_solvers wird z.Zt. nicht gelinked ...
!  ! weight ist in den Ecken allerdings auch zu Null gesetzt ...
!
!  max_div = 0.
!
!  do k = S3p, N3p
!     do j = S2p, N2p
!        do i = S1p, N1p
!           if (ABS(res(i,j,k)) >= max_div(1)) then
!              max_div(1) = ABS(res(i,j,k))
!           end if
!           if (ABS(res(i,j,k)*weight(i,j,k)) >= max_div(2)) then
!              max_div(2) = ABS(res(i,j,k)*weight(i,j,k))
!           end if
!        end do
!     end do
!  end do
!
!  call MPI_REDUCE(max_div,max_div_global,2,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
!
!  if (rank == 0) then
!     write(*,'(a)')
!     write(*,'(a,E25.17)') 'MAX(div(u))        =', max_div_global(1)
!     write(*,'(a,E25.17)') 'MAX(div(u)*weight) =', max_div_global(2)
!  end if
!
!  ! Für Auswertung der Iterationsstatistiken:
!  if (timestep == 0) max_div_init = max_div_global
!
!
!  end subroutine test_divergence
  
  
  
  
  
  
  
  
  
  
  
!  subroutine test_momentum
!
!  implicit none
!
!  integer                ::  i, j, k, m
!  real                   ::  max_err, max_err_global
!  real                   ::  dd1, dd2
!
!
!  dd1 = aRK(RK_steps)+bRK(RK_steps)
!
!  !===========================================================================================================
!  m         = 1
!  direction = 1
!
!  call gradient(m,pre,gpre)
!  call product_Helmholtz(vel(b1L,b2L,b3L,m),res)
!
!  max_err = 0.
!
!  do k = S31B, N31B
!     do j = S21B, N21B
!        do i = S11B, N11B
!           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
!           if (dd2 > max_err) max_err = dd2
!        end do
!     end do
!  end do
!
!  call MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
!
!  if (rank == 0) write(*,'(a      )')
!  if (rank == 0) write(*,'(a,E13.5)') 'MAX(res(vel_1)) =', max_err_global
!  !-----------------------------------------------------------------------------------------------------------
!  m         = 2
!  direction = 2
!
!  call gradient(m,pre,gpre)
!  call product_Helmholtz(vel(b1L,b2L,b3L,m),res)
!
!  max_err = 0.
!
!  do k = S32B, N32B
!     do j = S22B, N22B
!        do i = S12B, N12B
!           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
!           if (dd2 > max_err) max_err = dd2
!        end do
!     end do
!  end do
!
!  call MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
!
!  if (rank == 0) write(*,'(a      )')
!  if (rank == 0) write(*,'(a,E13.5)') 'MAX(res(vel_2)) =', max_err_global
!  !-----------------------------------------------------------------------------------------------------------
!  m         = 3
!  direction = 3
!
!  call gradient(m,pre,gpre)
!  call product_Helmholtz(vel(b1L,b2L,b3L,m),res)
!
!  max_err = 0.
!
!  do k = S33B, N33B
!     do j = S23B, N23B
!        do i = S13B, N13B
!           dd2 = ABS(rhs(i,j,k,m)-res(i,j,k)-dd1*gpre(i,j,k))
!           if (dd2 > max_err) max_err = dd2
!        end do
!     end do
!  end do
!
!  call MPI_REDUCE(max_err,max_err_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
!
!  if (rank == 0) write(*,'(a      )')
!  if (rank == 0) write(*,'(a,E13.5)') 'MAX(res(vel_3)) =', max_err_global
!  !===========================================================================================================
!
!
!  end subroutine test_momentum
!
  

  
!  subroutine analyze_matrix(gridtype) bind(c,name='fanalyze_matrix')! TEST!!! ggf. in anderes Modul??
!
!  implicit none
!
!  integer(c_int), intent(in)    ::  gridtype
!
!  integer                ::  i, j, k
!  integer                ::  ii, jj, kk
!  integer                ::  iM, jM
!
!  integer                ::  SS1, SS2, SS3
!  integer                ::  NN1, NN2, NN3
!
!  integer                ::  Nmax
!
!  complex, allocatable   ::  MM(:,:)
!  complex, allocatable   ::  NN(:,:)
!
!  integer, allocatable   ::  PP(:)
!  integer                ::  lerror
!
!  complex, allocatable   ::  EV  (:)
!  complex, allocatable   ::  work(:)
!  real   , allocatable   ::  rwork(:)
!
!  integer                ::  iEV
!  real                   ::  maxRe, minRe, diag
!  real                   ::  beta, ssum, maxdiag
!
!
!  if (gridtype == 0) then
!     SS1 = S1p
!     SS2 = S2p
!     SS3 = S3p
!     NN1 = N1p
!     NN2 = N2p
!     NN3 = N3p
!  else if (gridtype == 1) then
!     SS1 = S11B
!     SS2 = S21B
!     SS3 = S31B
!     NN1 = N11B
!     NN2 = N21B
!     NN3 = N31B
!  else if (gridtype == 2) then
!     SS1 = S12B
!     SS2 = S22B
!     SS3 = S32B
!     NN1 = N12B
!     NN2 = N22B
!     NN3 = N32B
!  else if (gridtype == 3) then
!     SS1 = S13B
!     SS2 = S23B
!     SS3 = S33B
!     NN1 = N13B
!     NN2 = N23B
!     NN3 = N33B
!  end if
!
!  multL = 1./Re
!
!  Nmax = (NN1-SS1+1)*(NN2-SS2+1)*(NN3-SS3+1)
!
!
!  allocate(MM(1:Nmax,1:Nmax))
!  allocate(NN(1:Nmax,1:Nmax))
!
!  allocate(EV(1:Nmax))
!
!  allocate(work (1:2*Nmax))
!  allocate(rwork(1:2*Nmax))
!
!  MM = 0.
!  NN = 0.
!
!  EV = 0.
!
!  work = 0.
!  rwork = 0.
!
!  jM = 1
!
!
!  do k = SS3, NN3
!     do j = SS2, NN2
!        do i = SS1, NN1
!
!           pre  = 0.
!           gpre = 0.
!           pre(i,j,k) = 1.
!
!           if (gridtype == 0) then
!              call product_div_grad (pre,gpre)
!           else if (gridtype == 1) then
!              direction = 1
!              call product_helmholtz(pre,gpre)
!           else if (gridtype == 2) then
!              direction = 2
!              call product_helmholtz(pre,gpre)
!           else if (gridtype == 3) then
!              direction = 3
!              call product_helmholtz(pre,gpre)
!           end if
!
!           iM = 1
!
!           do kk = SS3, NN3
!              do jj = SS2, NN2
!                 do ii = SS1, NN1
!
!                    MM(iM,jM) = CMPLX(gpre(ii,jj,kk),0.)
!                    iM = iM + 1
!                 end do
!              end do
!           end do
!
!           jM = jM + 1
!        end do
!     end do
!  end do
!
!  !*****************************************************
!  if (1 == 1) then
!     do k = S3p, N3p
!        do j = S2p, N2p
!           do i = S1p, N1p
!
!              pre  = 0.
!              nl   = 0.
!              vel  = 0.
!              pre(i,j,k) = 1.
!
!              call gradient(1,pre,nl(b1L,b2L,b3L,1))
!              call gradient(2,pre,nl(b1L,b2L,b3L,2))
!
!              call divergence_transp(1,pre,vel(b1L,b2L,b3L,1))
!              call divergence_transp(2,pre,vel(b1L,b2L,b3L,2))
!
!              write(57,'(5000E25.17)') nl (S11B:N11B,S21B:N21B,S31B:N31B,1), nl (S12B:N12B,S22B:N22B,S32B:N32B,2)
!              write(58,'(5000E25.17)') vel(S11B:N11B,S21B:N21B,S31B:N31B,1), vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
!
!           end do
!        end do
!     end do
!  end if
!  !*****************************************************
!
!
!  ! umskalieren auf diag(MM)=I:
!  !DO iM = 1, Nmax
!  !   diag = REAL(MM(iM,iM))
!  !   IF (diag /= 0.) THEN
!  !      DO jM = 1, Nmax
!  !         MM(iM,jM) = MM(iM,jM) / CMPLX(diag,0.)
!  !      END DO
!  !   END IF
!  !END DO
!
!
!
!  do iM = 1, Nmax
!     write(47,'(500E25.17)') REAL(MM(iM,1:Nmax))
!  end do
!
!  beta = 0.
!
!  do iM = 1, Nmax
!     ssum = 0.
!     do jM = 1, Nmax
!        ssum = ssum + ABS(REAL(MM(iM,jM)))
!     end do
!     ssum = ssum - ABS(REAL(MM(iM,iM)))
!     maxdiag = MAX(ABS(REAL(MM(iM,iM))),maxdiag)
!     beta = MAX(beta,ssum)
!  end do
!
!
!  call ZGEEV('N','V',Nmax,MM,Nmax,EV,NN,Nmax,NN,Nmax,work,2*Nmax,rwork,lerror)
!
!
!  write(*,*) lerror
!  write(*,*)
!
!  maxRe = -10.**16
!  minRe =  10.**16
!
!  do i = 1, Nmax
!     !WRITE(*,*) i, EV(i)
!     write(48,*) i, REAL(EV(i)), AIMAG(EV(i))
!
!     if (REAL(EV(i)) /= 0.) maxRe = MAX(maxRe,REAL(EV(i)))
!     minRe = MIN(minRe,REAL(EV(i)))
!  end do
!
!
!  !iEV = 1
!  !
!  !iM = 1
!  !DO j = SS2, NN2
!  !   DO i = SS1, NN1
!  !      WRITE(49,'(500E25.17)') x1u(i), x2p(j), REAL(NN(iM,1:Nmax))
!  !      iM = iM+1
!  !   END DO
!  !   WRITE(49,*)
!  !END DO
!
!  write(*,*) maxRe, minRe
!  write(*,*) 2.*beta+1.
!  write(*,*) beta, maxdiag
!  write(*,*) beta+maxdiag
!
!
!  end subroutine analyze_matrix
  
  
  
end module cmod_test
