!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module cmod_setup
  
    use iso_c_binding
  
    use mod_dims
    use mod_vars
    use mod_lib ! (num_to_string)

!    use mpi
  
    private
  
    public init_general, init_parallel, init_boundaries, init_limits

    INCLUDE 'mpif.h'
  
contains
  
    !pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
  
  
  
    subroutine init_general() bind(c,name='finit_generallight')
  
        implicit none
  
        integer                ::  b, g, m
        integer                ::  stride(1:3)
        integer                ::  MM(1:3,1:n_grids_max)
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Da bei den Helmholtz-Problemen für die Geschwindigkeiten auf dem feinsten Gitter          !
        !                mindestens soviele Punkte in Wand-normaler Richtung vorhanden sind wie bei dem Poisson-   !
        !                Problem bzw. die selben groben Gittern verwendet werden, kann es hier zu keinem Konflikt  !
        !                kommen.                                                                                   !
        !----------------------------------------------------------------------------------------------------------!
  
  
        !===========================================================================================================
        !=== Zeitintegration =======================================================================================
        !===========================================================================================================
        ! - Euler-explizit ist nicht erlaubt.
        ! - Euler_yes (== Re --> infty) ==> ist immer explizit ==> timeint_mode = 1.
        ! - thetaL ist nur bei timeint_mode == 0 von Bedeutung.
  
        if (Euler_yes) timeint_mode = 1
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Dimension =============================================================================================
        !===========================================================================================================
        if (M3 == 2) then
            dimens = 2
        else
            dimens = 3
        end if
        !===========================================================================================================
  
  
#ifdef ALLOC
        !===========================================================================================================
        !=== Domain- und Blockspezifikationen ======================================================================
        !===========================================================================================================
        !  N1 = 1+(M1-1)/NB1
        !  N2 = 1+(M2-1)/NB2
        !  N3 = 1+(M3-1)/NB3
        !
        dimS1 = 2*ndL*NB1
        dimS2 = 2*ndL*NB2
        dimS3 = 2*ndL*NB3
        !===========================================================================================================
#endif
  
  
        !===========================================================================================================
        !=== Multigrid-Gitterspezifikationen =======================================================================
        !===========================================================================================================
        ! TEST!!! Test schreiben, dass Preconditioning von Helmholtz-Gleichungen (vel+conc) momentan nicht voll unterstützt ist!
  
        MM(1:3,1) = (/M1 ,M2 ,M3 /)
        NN(1:3,1) = (/N1 ,N2 ,N3 /)
        NB(1:3,1) = (/NB1,NB2,NB3/)
  
        n_grids  = 1
        n_gather = 1
        stride   = 1
  
        ! (alte Version geloescht am 24.01.2010)
        do g = 2, n_grids_max
            do m = 1, 3
                !--- global number of grid points ---
                if (MOD(MM(m,g-1)-1,2) == 0 .and. MM(m,g-1) >= 5) then ! 5=2*(Nmin-1)+1, Nmin=3 for Laplacian
                    MM(m,g) = (MM(m,g-1)-1) / 2 + 1
                    n_grids = g
                else
                    MM(m,g) = MM(m,g-1)
                end if
        
                !--- find largest possible number of blocks ---
                NB(m,g) = 1
                do b = 2, NB(m,g-1)
                    if (MOD(MM(m,g)-1,b) == 0 .and. MOD(NB(m,g-1),b) == 0) then
                        ! prolongation to next coarser grid must be possible on each block, independent of the other blocks!
                        ! NN(m,g+1)-1 .GE. 1 (next coarser grid)
                        if (MOD((MM(m,g)-1)/b,2) == 0 .and. (MM(m,g)-1)/b >= 2) NB(m,g) = b
                    end if
                end do
        
                !--- enforce gathering of coarsest grid to one processor ---
                if ((MOD(MM(m,g)-1,2) /= 0 .or. MM(m,g) < 5) .and. NB(m,g) > 1) NB(m,g) = 1
        
                !--- number of gathered blocks ---
                n_gather(m,g) = NB(m,g-1)/NB(m,g)
        
                !--- local grid size ---
                NN(m,g) = (MM(m,g)-1)/NB(m,g)+1
            end do
     
        end do
  
  
        if (n_grids > n_grids_limit) n_grids = n_grids_limit
  
  
        if (rank == 0) then
            open(10,file='test_multigrid_properties.txt', status='UNKNOWN')
            write(10,'(a     )') 'grid dimensions (total):'
            write(10,'(a,20i6)') '         M1: ', (NN(1,1:n_grids)-1)*NB(1,1:n_grids)+1
            write(10,'(a,20i6)') '         M2: ', (NN(2,1:n_grids)-1)*NB(2,1:n_grids)+1
            write(10,'(a,20i6)') '         M3: ', (NN(3,1:n_grids)-1)*NB(3,1:n_grids)+1
            write(10,*)
            write(10,'(a     )') 'grid dimensions (per processor block):'
            write(10,'(a,20i6)') '         N1: ', NN(1,1:n_grids)
            write(10,'(a,20i6)') '         N2: ', NN(2,1:n_grids)
            write(10,'(a,20i6)') '         N3: ', NN(3,1:n_grids)
            write(10,*)
            write(10,'(a     )') 'Processor blocks:'
            write(10,'(a,20i6)') '        NB1: ', NB(1,1:n_grids)
            write(10,'(a,20i6)') '        NB2: ', NB(2,1:n_grids)
            write(10,'(a,20i6)') '        NB3: ', NB(3,1:n_grids)
            write(10,*)
            write(10,'(a     )') 'Processor block gathering:'
            write(10,'(a,20i6)') 'n_gather(1): ', n_gather(1,1:n_grids)
            write(10,'(a,20i6)') 'n_gather(2): ', n_gather(2,1:n_grids)
            write(10,'(a,20i6)') 'n_gather(3): ', n_gather(3,1:n_grids)
            close(10)
        end if
  
  
        !===========================================================================================================
        !=== Felder allocieren =====================================================================================
        !===========================================================================================================
        ! TEST!!! alloc.f90 ganz hineinziehen?
#ifdef ALLOC
  INCLUDE 'alloc.f90'
#ifdef NONBOUSSINESQ
        allocate(dens(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)) ! checkt der Compiler sonst nicht, wurde aus alloc.f90 herausgezogen ...
#endif
#endif
  
  
        if (n_grids >=  2) allocate(vec2A (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
        if (n_grids >=  2) allocate(vec2B (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
        if (n_grids >=  2) allocate(vec2C (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
  
        if (n_grids >=  3) allocate(vec3A (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
        if (n_grids >=  3) allocate(vec3B (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
        if (n_grids >=  3) allocate(vec3C (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
  
        if (n_grids >=  4) allocate(vec4A (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
        if (n_grids >=  4) allocate(vec4B (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
        if (n_grids >=  4) allocate(vec4C (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
  
        if (n_grids >=  5) allocate(vec5A (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
        if (n_grids >=  5) allocate(vec5B (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
        if (n_grids >=  5) allocate(vec5C (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
  
        if (n_grids >=  6) allocate(vec6A (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
        if (n_grids >=  6) allocate(vec6B (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
        if (n_grids >=  6) allocate(vec6C (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
  
        if (n_grids >=  7) allocate(vec7A (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
        if (n_grids >=  7) allocate(vec7B (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
        if (n_grids >=  7) allocate(vec7C (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
  
        if (n_grids >=  8) allocate(vec8A (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
        if (n_grids >=  8) allocate(vec8B (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
        if (n_grids >=  8) allocate(vec8C (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
  
        if (n_grids >=  9) allocate(vec9A (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
        if (n_grids >=  9) allocate(vec9B (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
        if (n_grids >=  9) allocate(vec9C (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
  
        if (n_grids >= 10) allocate(vec10A(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
        if (n_grids >= 10) allocate(vec10B(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
        if (n_grids >= 10) allocate(vec10C(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
  
        if (n_grids >= 11) allocate(vec11A(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
        if (n_grids >= 11) allocate(vec11B(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
        if (n_grids >= 11) allocate(vec11C(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
  
        if (n_grids >= 12) allocate(vec12A(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
        if (n_grids >= 12) allocate(vec12B(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
        if (n_grids >= 12) allocate(vec12C(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
  
        if (n_grids >= 13) allocate(vec13A(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
        if (n_grids >= 13) allocate(vec13B(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
        if (n_grids >= 13) allocate(vec13C(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
  
        if (n_grids >= 14) allocate(vec14A(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
        if (n_grids >= 14) allocate(vec14B(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
        if (n_grids >= 14) allocate(vec14C(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
  
        if (n_grids >= 15) allocate(vec15A(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
        if (n_grids >= 15) allocate(vec15B(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
  
  
        if (n_grids >=  2) allocate(psi_rel2 (b1L:(NN(1,2 )+b1U),b2L:(NN(2,2 )+b2U),b3L:(NN(3,2 )+b3U)))
        if (n_grids >=  3) allocate(psi_rel3 (b1L:(NN(1,3 )+b1U),b2L:(NN(2,3 )+b2U),b3L:(NN(3,3 )+b3U)))
        if (n_grids >=  4) allocate(psi_rel4 (b1L:(NN(1,4 )+b1U),b2L:(NN(2,4 )+b2U),b3L:(NN(3,4 )+b3U)))
        if (n_grids >=  5) allocate(psi_rel5 (b1L:(NN(1,5 )+b1U),b2L:(NN(2,5 )+b2U),b3L:(NN(3,5 )+b3U)))
        if (n_grids >=  6) allocate(psi_rel6 (b1L:(NN(1,6 )+b1U),b2L:(NN(2,6 )+b2U),b3L:(NN(3,6 )+b3U)))
        if (n_grids >=  7) allocate(psi_rel7 (b1L:(NN(1,7 )+b1U),b2L:(NN(2,7 )+b2U),b3L:(NN(3,7 )+b3U)))
        if (n_grids >=  8) allocate(psi_rel8 (b1L:(NN(1,8 )+b1U),b2L:(NN(2,8 )+b2U),b3L:(NN(3,8 )+b3U)))
        if (n_grids >=  9) allocate(psi_rel9 (b1L:(NN(1,9 )+b1U),b2L:(NN(2,9 )+b2U),b3L:(NN(3,9 )+b3U)))
        if (n_grids >= 10) allocate(psi_rel10(b1L:(NN(1,10)+b1U),b2L:(NN(2,10)+b2U),b3L:(NN(3,10)+b3U)))
        if (n_grids >= 11) allocate(psi_rel11(b1L:(NN(1,11)+b1U),b2L:(NN(2,11)+b2U),b3L:(NN(3,11)+b3U)))
        if (n_grids >= 12) allocate(psi_rel12(b1L:(NN(1,12)+b1U),b2L:(NN(2,12)+b2U),b3L:(NN(3,12)+b3U)))
        if (n_grids >= 13) allocate(psi_rel13(b1L:(NN(1,13)+b1U),b2L:(NN(2,13)+b2U),b3L:(NN(3,13)+b3U)))
        if (n_grids >= 14) allocate(psi_rel14(b1L:(NN(1,14)+b1U),b2L:(NN(2,14)+b2U),b3L:(NN(3,14)+b3U)))
        if (n_grids >= 15) allocate(psi_rel15(b1L:(NN(1,15)+b1U),b2L:(NN(2,15)+b2U),b3L:(NN(3,15)+b3U)))
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Konzentration =========================================================================================
        !===========================================================================================================
        if (concentration_yes) then
            do m = 1, n_conc
                usReSc(1,m) = gravity(1)*usc(m)*Re*Sc(m)
                usReSc(2,m) = gravity(2)*usc(m)*Re*Sc(m)
                usReSc(3,m) = gravity(3)*usc(m)*Re*Sc(m)
        
                us_vec(1,m) = gravity(1)*usc(m)
                us_vec(2,m) = gravity(2)*usc(m)
                us_vec(3,m) = gravity(3)*usc(m)
            end do
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Konvergenzstatistik ===================================================================================
        !===========================================================================================================
        countO = 0
        countP = 0
        countH = 0
  
        ratioO = 0.
        ratioH = 0.
        ratioP = 0.
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Abbruchkriterium der äusseren Iteration ===============================================================
        !===========================================================================================================
        allocate(precOffset(1:RK_steps,1:n_it_outer))
        allocate(precRatio (1:RK_steps,1:n_it_outer))
  
        do substep = 1, RK_steps
            precOffset(substep,:) = precOffset0(substep)
            precRatio (substep,:) = precRatio0 (substep)
        end do
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Restart-Nr. als String für File-Namen =================================================================
        !===========================================================================================================
        call num_to_string(3,restart,restart_char)
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Feld-I/O ==============================================================================================
        !===========================================================================================================
        write_large = .false.
        write_med   = .false.
        write_small = .false.
  
        if (stride_large(1) /= 0 .and. stride_large(2) /= 0 .and. stride_large(3) /= 0) write_large = .true.
        if (stride_med  (1) /= 0 .and. stride_med  (2) /= 0 .and. stride_med  (3) /= 0) write_med   = .true.
        if (stride_small(1) /= 0 .and. stride_small(2) /= 0 .and. stride_small(3) /= 0) write_small = .true.
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Smagorinsky LES-Modell ================================================================================
        !===========================================================================================================
        if (LES_mode == 2) allocate(smag_diffus(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:dimens))
        !===========================================================================================================
  
  
        epsU0 = epsU ! TEST!!! Woanders hin?
  
  
    !IF (rank == 0) CALL SYSTEM('echo "IMPACT is run by $USER on $HOST which is a $HOSTTYPE" | mail -s "IMPACT usage notification" henniger@ifd.mavt.ethz.ch')
  
  
    end subroutine init_general
  
  
  
  
  
  
  
  
  
    !> \brief intializes comunicator, COMM_CART, slice
    !! \deprecated
    subroutine init_parallel() bind(c,name='finit_parallel')
  
        implicit none
  
        integer                ::  ijkB(1:3)      ! mpi grid coordinates
        integer                ::  periodic(1:3)  ! array for mpir to signal where periodic grid is, from manual should be bool
  
  
        ! init periodic
        if (BC_1L_global == -1) then
            periodic(1) = 1
        else
            periodic(1) = 0
        end if
  
        if (BC_2L_global == -1) then
            periodic(2) = 1
        else
            periodic(2) = 0
        end if
  
        if (BC_3L_global == -1) then
            periodic(3) = 1
        else
            periodic(3) = 0
        end if
  
  
        ! Macht keinen messbaren Unterschied (falls doch irgendwann, dann sollte der CALL auch auf die Grobgitter-Kommunikatoren auch angewandt werden!):
        !CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,(/NB1,NB2,NB3/),periodic,.FALSE.,COMM_CART,merror)
        !CALL MPI_CART_MAP(COMM_CART,3,(/NB1,NB2,NB3/),periodic,rank,merror)

        ! .true. means ranking may be reorderd
        ! comm_cart comm with cartesian grid informations
        call MPI_CART_CREATE(MPI_COMM_WORLD,3,(/NB1,NB2,NB3/),periodic,.true.,COMM_CART,merror)
        ! gets rank from COMM_CART
        call MPI_COMM_RANK  (COMM_CART,rank,merror)
        ! gets coordinates in xyz direction from rank and comm_cart
        call MPI_CART_COORDS(COMM_CART,rank,3,ijkB,merror)
  
        ! stores coordinates in a fortran fasion?
        iB(1,1) = ijkB(1)+1
        iB(2,1) = ijkB(2)+1
        iB(3,1) = ijkB(3)+1
  
        ! computes index ofset
        iShift = (iB(1,1)-1)*(N1-1)
        jShift = (iB(2,1)-1)*(N2-1)
        kShift = (iB(3,1)-1)*(N3-1)
  
        call MPI_CART_SHIFT(COMM_CART,0,1,rank1L,rank1U,merror)
        call MPI_CART_SHIFT(COMM_CART,1,1,rank2L,rank2U,merror)
        call MPI_CART_SHIFT(COMM_CART,2,1,rank3L,rank3U,merror)
        !                             ^ ^   ^      ^
        !                             | |   |      |
        !                             d d   r      r
        !                             i i   a      a
        !                             r s   n      n
        !                             e p   k      k
        !                             c l   s      d
        !                             t a   o      e
        !                             i c   u      s
        !                             o m   r      t
        !                             n e   c
        !                               n   e
        !                               t
  
  
        call MPI_CART_SUB(COMM_CART,(/0,1,1/),COMM_SLICE1,merror)
        call MPI_CART_SUB(COMM_CART,(/1,0,1/),COMM_SLICE2,merror)
        call MPI_CART_SUB(COMM_CART,(/1,1,0/),COMM_SLICE3,merror)
  
        call MPI_CART_SUB(COMM_CART,(/1,0,0/),COMM_BAR1  ,merror)
        call MPI_CART_SUB(COMM_CART,(/0,1,0/),COMM_BAR2  ,merror)
        call MPI_CART_SUB(COMM_CART,(/0,0,1/),COMM_BAR3  ,merror)
  
  
        call MPI_COMM_RANK(COMM_SLICE1,rank_slice1,merror)
        call MPI_COMM_RANK(COMM_SLICE2,rank_slice2,merror)
        call MPI_COMM_RANK(COMM_SLICE3,rank_slice3,merror)
  
        call MPI_COMM_RANK(COMM_BAR1  ,rank_bar1  ,merror)
        call MPI_COMM_RANK(COMM_BAR2  ,rank_bar2  ,merror)
        call MPI_COMM_RANK(COMM_BAR3  ,rank_bar3  ,merror)
  
  
        ! Evtl. überflüssig ...
        call MPI_ERRHANDLER_SET(COMM_CART  ,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_SLICE1,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_SLICE2,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_SLICE3,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_BAR1  ,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_BAR2  ,MPI_ERRORS_ARE_FATAL,merror)
        call MPI_ERRHANDLER_SET(COMM_BAR3  ,MPI_ERRORS_ARE_FATAL,merror)
  
  
    end subroutine init_parallel
  
  
  
  
  
  
  
  
  
  
  
    !> \deprecated
    subroutine init_boundaries() bind(c,name='finit_boundaries')
  
        implicit none
  
        integer                ::  m
  
  
        ! Default: Nachbarblöcke
        BC_1L = 0
        BC_1U = 0
  
        BC_2L = 0
        BC_2U = 0
  
        BC_3L = 0
        BC_3U = 0
  
        ! Spezialfall periodische RB bei einem Block:
        if (BC_1L_global == -1 .and. NB1 == 1) then
            BC_1L = -1
            BC_1U = -1
        end if
        if (BC_2L_global == -1 .and. NB2 == 1) then
            BC_2L = -1
            BC_2U = -1
        end if
        if (BC_3L_global == -1 .and. NB3 == 1) then
            BC_3L = -1
            BC_3U = -1
        end if
  
        ! Wände oder Symmetrie-RB:
        if (rank1L < 0) BC_1L = BC_1L_global
        if (rank1U < 0) BC_1U = BC_1U_global
  
        if (rank2L < 0) BC_2L = BC_2L_global
        if (rank2U < 0) BC_2U = BC_2U_global
  
        if (rank3L < 0) BC_3L = BC_3L_global
        if (rank3U < 0) BC_3U = BC_3U_global
  
  
        !===========================================================================================================
        !=== Konzentration =========================================================================================
        !===========================================================================================================
        if (concentration_yes) then
            do m = 1, n_conc
        
                !-----------------------------------------------------------------------------------------------------
                if (BC_1L <= 0) then
                    if (BC_1L == -2 .and. us_vec(1,m) > 0.) then
                        BCc_1L(m) = 3
                    else
                        BCc_1L(m) = BC_1L
                    end if
                else
                    if (us_vec(1,m) < 0. .or. outlet(1,1,1) .or. isopycnal(1,1,m)) then
                        BCc_1L(m) = 1
                    else
                        BCc_1L(m) = 3
                    end if
                end if
        
                if (BC_1U <= 0) then
                    if (BC_1U == -2 .and. us_vec(1,m) < 0.) then
                        BCc_1U(m) = 3
                    else
                        BCc_1U(m) = BC_1U
                    end if
                else
                    if (us_vec(1,m) > 0. .or. outlet(1,2,1) .or. isopycnal(1,2,m)) then
                        BCc_1U(m) = 1
                    else
                        BCc_1U(m) = 3
                    end if
                end if
                !-----------------------------------------------------------------------------------------------------
                if (BC_2L <= 0) then
                    if (BC_2L == -2 .and. us_vec(2,m) > 0.) then
                        BCc_2L(m) = 3
                    else
                        BCc_2L(m) = BC_2L
                    end if
                else
                    if (us_vec(2,m) < 0. .or. outlet(2,1,2) .or. isopycnal(2,1,m)) then
                        BCc_2L(m) = 1
                    else
                        BCc_2L(m) = 3
                    end if
                end if
        
                if (BC_2U <= 0) then
                    if (BC_2U == -2 .and. us_vec(2,m) < 0.) then
                        BCc_2U(m) = 3
                    else
                        BCc_2U(m) = BC_2U
                    end if
                else
                    if (us_vec(2,m) > 0. .or. outlet(2,2,2) .or. isopycnal(2,2,m)) then
                        BCc_2U(m) = 1
                    else
                        BCc_2U(m) = 3
                    end if
                end if
                !-----------------------------------------------------------------------------------------------------
                if (BC_3L <= 0) then
                    if (BC_3L == -2 .and. us_vec(3,m) > 0.) then
                        BCc_3L(m) = 3
                    else
                        BCc_3L(m) = BC_3L
                    end if
                else
                    if (us_vec(3,m) < 0. .or. outlet(3,1,3) .or. isopycnal(3,1,m)) then
                        BCc_3L(m) = 1
                    else
                        BCc_3L(m) = 3
                    end if
                end if
        
                if (BC_3U <= 0) then
                    if (BC_3U == -2 .and. us_vec(3,m) < 0.) then
                        BCc_3U(m) = 3
                    else
                        BCc_3U(m) = BC_3U
                    end if
                else
                    if (us_vec(3,m) > 0. .or. outlet(3,2,3) .or. isopycnal(3,2,m)) then
                        BCc_3U(m) = 1
                    else
                        BCc_3U(m) = 3
                    end if
                end if
               !-----------------------------------------------------------------------------------------------------
            end do
        end if
    !===========================================================================================================
  
  
    end subroutine init_boundaries
  
  
  
  
  
  
  
  
  
  
  
    subroutine init_limits() bind(c,name='finit_limits')
  
        implicit none
  
        integer                ::  bar1_out(1:NB1)
        integer                ::  bar2_out(1:NB2)
        integer                ::  bar3_out(1:NB3)
  
        integer                ::  g, m
  
        integer                ::  rank_bar1
        integer                ::  rank_bar2
        integer                ::  rank_bar3
  
        !*************************************************************************
        integer                ::  base_group, new_group, comm_temp
        integer                ::  stride(1:3), counter
        integer                ::  rank_comm1(1:n_grids_max), rank_comm2(1:n_grids_max)
        integer, allocatable   ::  new_ranks(:,:,:)
        integer                ::  i, ii, iii, iimax, iiShift
        integer                ::  j, jj, jjj, jjmax, jjShift
        integer                ::  k, kk, kkk, kkmax, kkShift
        integer                ::  offs_global(1:3,NB1*NB2*NB3)
        integer                ::  sizs_global(1:3,NB1*NB2*NB3)
        integer                ::  count_comm
        !*************************************************************************
  
        integer                ::  periodic(1:3)
        logical                ::  member_yes
  
        if (BC_1L_global == -1) then
            periodic(1) = 1
        else
            periodic(1) = 0
        end if
  
        if (BC_2L_global == -1) then
            periodic(2) = 1
        else
            periodic(2) = 0
        end if
  
        if (BC_3L_global == -1) then
            periodic(3) = 1
        else
            periodic(3) = 0
        end if
  
  
  
        !-----------------------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - S1p = S1R, N1p = N1R (in Relaxationsroutinen)                                                          !
        !              - Default sind periodische / Nachbarblock-RB                                                             !
        !-----------------------------------------------------------------------------------------------------------------------!
  
        !-----------------------------------------------------------------------------------------------------------------------!
        ! ls =  0: erster  Block ab Wand hat 2 Schichten wand-normale Geschwindigkeiten mehr als die folgenden Blöcke           !
        !                         -||-       1 Schicht   wand-parallele Geschw. / Druck mehr als die folgenden Blöcke           !
        !                                                                                                                       !
        ! ls = -1: erster / letzter Block ab Wand haben 1 Schicht wand-normale Geschwindigkeiten mehr als die Blöcke dazwischen !
        !                   letzter Block ab Wand hat   1 Schicht wand-parallele Geschw. / Druck mehr als die Blöcke davor      !
        !-----------------------------------------------------------------------------------------------------------------------!
  
  
        if (ls1 ==  0) ex1 =  1
        if (ls1 == -1) ex1 = -1
  
        if (ls2 ==  0) ex2 =  1
        if (ls2 == -1) ex2 = -1
  
        if (ls3 ==  0) ex3 =  1
        if (ls3 == -1) ex3 = -1
  
  
        S1R  = 2 + ls1
        S2R  = 2 + ls2
        S3R  = 2 + ls3
  
        S11R = 2 + ls1
        S22R = 2 + ls2
        S33R = 2 + ls3
  
        d1R  = ls1
        d2R  = ls2
        d3R  = ls3
  
        d11R = ls1
        d22R = ls2
        d33R = ls3
  
  
        ! Druck / Konzentrationen (MIT Rand):
!        S1p  = 2 + ls1
!        S2p  = 2 + ls2
!        S3p  = 2 + ls3
!
!        N1p  = N1 + ls1
!        N2p  = N2 + ls2
!        N3p  = N3 + ls3
  
  
        ! Geschwindigkeiten (Feld MIT Rand)
        S11B = 2 + ls1
        S21B = 2 + ls2
        S31B = 2 + ls3
  
        S12B = 2 + ls1
        S22B = 2 + ls2
        S32B = 2 + ls3
  
        S13B = 2 + ls1
        S23B = 2 + ls2
        S33B = 2 + ls3
  
        N11B = N1 + ls1
        N21B = N2 + ls2
        N31B = N3 + ls3
  
        N12B = N1 + ls1
        N22B = N2 + ls2
        N32B = N3 + ls3
  
        N13B = N1 + ls1
        N23B = N2 + ls2
        N33B = N3 + ls3
  
        ! Geschwindigkeiten (Feld OHNE Rand)
        S11 = 2 + ls1
        S21 = 2 + ls2
        S31 = 2 + ls3
  
        S12 = 2 + ls1
        S22 = 2 + ls2
        S32 = 2 + ls3
  
        S13 = 2 + ls1
        S23 = 2 + ls2
        S33 = 2 + ls3
  
        N11 = N1 + ls1
        N21 = N2 + ls2
        N31 = N3 + ls3
  
        N12 = N1 + ls1
        N22 = N2 + ls2
        N32 = N3 + ls3
  
        N13 = N1 + ls1
        N23 = N2 + ls2
        N33 = N3 + ls3
  
  
        ! Konzentrationen (Feld OHNE Rand):
        do m = 1, n_conc
            S1c(m)  = 2 + ls1
            S2c(m)  = 2 + ls2
            S3c(m)  = 2 + ls3
     
            N1c(m)  = N1 + ls1
            N2c(m)  = N2 + ls2
            N3c(m)  = N3 + ls3
        end do
  
        !===========================================================================================================
        if (BC_1L > 0) then
!            S1p  = 1
            S1R  = 1
            S11R = 2
     
            S11B = 0
            S12B = 1
            S13B = 1
     
            S11  = 1
            S12  = 2
            S13  = 2
        end if
        if (BC_1L == -2) then
!            S1p  = 1
            S1R  = 1
            S11R = 1
     
            S11B = 1
            S12B = 1
            S13B = 1
     
            S11  = 1
            S12  = 1
            S13  = 1
        end if
        do m = 1, n_conc
            if (BCc_1L(m) > 0) S1c(m) = 2
            if (BCc_1L(m)  == -2) S1c(m) = 1
        end do
        !-----------------------------------------------------------------------------------------------------------
        if (BC_1U > 0) then
!            N1p  = N1
            d1R  =  0
            d11R = -1
     
            N11B = N1
            N12B = N1
            N13B = N1
     
            N11  = N1-1
            N12  = N1-1
            N13  = N1-1
        end if
        if (BC_1U == -2) then
!            N1p  = N1
            d1R  =  0
            d11R =  0
     
            N11B = N1-1
            N12B = N1
            N13B = N1
     
            N11  = N1-1
            N12  = N1
            N13  = N1
        end if
        do m = 1, n_conc
            if (BCc_1U(m) > 0) N1c(m) = N1-1
            if (BCc_1U(m)  == -2) N1c(m) = N1
        end do
        !===========================================================================================================
        if (BC_2L > 0) then
!            S2p  = 1
            S2R  = 1
            S22R = 2
     
            S21B = 1
            S22B = 0
            S23B = 1
     
            S21  = 2
            S22  = 1
            S23  = 2
        end if
        if (BC_2L == -2) then
!            S2p  = 1
            S2R  = 1
            S22R = 1
     
            S21B = 1
            S22B = 1
            S23B = 1
     
            S21  = 1
            S22  = 1
            S23  = 1
        end if
        do m = 1, n_conc
            if (BCc_2L(m) > 0) S2c(m) = 2
            if (BCc_2L(m)  == -2) S2c(m) = 1
        end do
        !-----------------------------------------------------------------------------------------------------------
        if (BC_2U > 0) then
!            N2p  = N2
            d2R  =  0
            d22R = -1
     
            N21B = N2
            N22B = N2
            N23B = N2
     
            N21  = N2-1
            N22  = N2-1
            N23  = N2-1
        end if
        if (BC_2U == -2) then
!            N2p  = N2
            d2R  =  0
            d22R =  0
     
            N21B = N2
            N22B = N2-1
            N23B = N2
     
            N21  = N2
            N22  = N2-1
            N23  = N2
        end if
        do m = 1, n_conc
            if (BCc_2U(m) > 0) N2c(m) = N2-1
            if (BCc_2U(m)  == -2) N2c(m) = N2
        end do
        !===========================================================================================================
        if (BC_3L > 0) then
!            S3p  = 1
            S3R  = 1
            S33R = 2
     
            S31B = 1
            S32B = 1
            S33B = 0
     
            S31  = 2
            S32  = 2
            S33  = 1
        end if
        if (BC_3L == -2) then
!            S3p  = 1
            S3R  = 1
            S33R = 1
     
            S31B = 1
            S32B = 1
            S33B = 1
     
            S31  = 1
            S32  = 1
            S33  = 1
        end if
        do m = 1, n_conc
            if (BCc_3L(m) > 0) S3c(m) = 2
            if (BCc_3L(m)  == -2) S3c(m) = 1
        end do
        !-----------------------------------------------------------------------------------------------------------
        if (BC_3U > 0) then
!            N3p  = N3
            d3R  =  0
            d33R = -1
     
            N31B = N3
            N32B = N3
            N33B = N3
     
            N31  = N3-1
            N32  = N3-1
            N33  = N3-1
        end if
        if (BC_3U == -2) then
!            N3p  = N3
            d3R  =  0
            d33R =  0
     
            N31B = N3
            N32B = N3
            N33B = N3-1
     
            N31  = N3
            N32  = N3
            N33  = N3-1
        end if
        do m = 1, n_conc
            if (BCc_3U(m) > 0) N3c(m) = N3-1
            if (BCc_3U(m)  == -2) N3c(m) = N3
        end do
  
  
  
        !*********************************************************************
        participate_yes = .false.
  
        stride     = 1
        count_comm = 0
  
        do g = 1, n_grids
     
            !--------------------------------------------------------------
            stride(1) = stride(1)*n_gather(1,g)
            stride(2) = stride(2)*n_gather(2,g)
            stride(3) = stride(3)*n_gather(3,g)
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            iB(1,g) = (iB(1,1)-1)*NB(1,g)/NB1 + 1
            iB(2,g) = (iB(2,1)-1)*NB(2,g)/NB2 + 1
            iB(3,g) = (iB(3,1)-1)*NB(3,g)/NB3 + 1
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            ! Analog "init_boundaries"
            BC(1:2,1,g) = (/BC_1L,BC_1U/)
            BC(1:2,2,g) = (/BC_2L,BC_2U/)
            BC(1:2,3,g) = (/BC_3L,BC_3U/)
     
            if (iB(1,g) == 1       .and. .not. BC_1L_global == -1) BC(1,1,g) = BC_1L_global
            if (iB(1,g) == NB(1,g) .and. .not. BC_1U_global == -1) BC(2,1,g) = BC_1U_global
     
            if (iB(2,g) == 1       .and. .not. BC_2L_global == -1) BC(1,2,g) = BC_2L_global
            if (iB(2,g) == NB(2,g) .and. .not. BC_2U_global == -1) BC(2,2,g) = BC_2U_global
     
            if (iB(3,g) == 1       .and. .not. BC_3L_global == -1) BC(1,3,g) = BC_3L_global
            if (iB(3,g) == NB(3,g) .and. .not. BC_3U_global == -1) BC(2,3,g) = BC_3U_global
     
     
            if (NB(1,g) == 1 .and. BC_1L_global == -1) BC(1:2,1,g) = BC_1L_global
            if (NB(2,g) == 1 .and. BC_2L_global == -1) BC(1:2,2,g) = BC_2L_global
            if (NB(3,g) == 1 .and. BC_3L_global == -1) BC(1:2,3,g) = BC_3L_global
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            SNB(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
            SNB(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
            SNB(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)
     
            if (BC(1,1,g) > 0) SNB(1,1,g) = 1
            if (BC(1,1,g)  == -2) SNB(1,1,g) = 1
            if (BC(2,1,g) > 0) SNB(2,1,g) = NN(1,g)
            if (BC(2,1,g)  == -2) SNB(2,1,g) = NN(1,g)
     
            if (BC(1,2,g) > 0) SNB(1,2,g) = 1
            if (BC(1,2,g)  == -2) SNB(1,2,g) = 1
            if (BC(2,2,g) > 0) SNB(2,2,g) = NN(2,g)
            if (BC(2,2,g)  == -2) SNB(2,2,g) = NN(2,g)
     
            if (BC(1,3,g) > 0) SNB(1,3,g) = 1
            if (BC(1,3,g)  == -2) SNB(1,3,g) = 1
            if (BC(2,3,g) > 0) SNB(2,3,g) = NN(3,g)
            if (BC(2,3,g)  == -2) SNB(2,3,g) = NN(3,g)
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            SNF(1:2,1,g) = (/2+ls1,NN(1,g)+ls1/)
            SNF(1:2,2,g) = (/2+ls2,NN(2,g)+ls2/)
            SNF(1:2,3,g) = (/2+ls3,NN(3,g)+ls3/)
     
            if (BC(1,1,g) > 0) SNF(1,1,g) = 2
            if (BC(1,1,g)  == -2) SNF(1,1,g) = 1
            if (BC(2,1,g) > 0) SNF(2,1,g) = NN(1,g)-1
            if (BC(2,1,g)  == -2) SNF(2,1,g) = NN(1,g)
     
            if (BC(1,2,g) > 0) SNF(1,2,g) = 2
            if (BC(1,2,g)  == -2) SNF(1,2,g) = 1
            if (BC(2,2,g) > 0) SNF(2,2,g) = NN(2,g)-1
            if (BC(2,2,g)  == -2) SNF(2,2,g) = NN(2,g)
     
            if (BC(1,3,g) > 0) SNF(1,3,g) = 2
            if (BC(1,3,g)  == -2) SNF(1,3,g) = 1
            if (BC(2,3,g) > 0) SNF(2,3,g) = NN(3,g)-1
            if (BC(2,3,g)  == -2) SNF(2,3,g) = NN(3,g)
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            call MPI_CART_SHIFT(COMM_CART,0,stride(1),ngb(1,1,g),ngb(2,1,g),merror)
            call MPI_CART_SHIFT(COMM_CART,1,stride(2),ngb(1,2,g),ngb(2,2,g),merror)
            call MPI_CART_SHIFT(COMM_CART,2,stride(3),ngb(1,3,g),ngb(2,3,g),merror)
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            !--- coarse grid communicators (comm1) ---
            ! siehe auch Kommentare in alten Versionen bis 19.01.2010!
            if (g >= 2 .and. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) > 1) then
     
                allocate(new_ranks(1:NB(1,g),1:NB(2,g),1:NB(3,g)))
        
                !-----------------------------------------------------------------------------------
                call MPI_COMM_GROUP(COMM_CART,base_group,merror)
                participate_yes(g) = .false.
                do k = 1, NB(3,g)
                    do j = 1, NB(2,g)
                        do i = 1, NB(1,g)
                            call MPI_CART_RANK(COMM_CART,(/MOD((i-1)*stride(1),NB1),  &
                                &          MOD((j-1)*stride(2),NB2),  &
                                &          MOD((k-1)*stride(3),NB3)/),new_ranks(i,j,k),merror)
                            if (rank == new_ranks(i,j,k)) participate_yes(g) = .true.
                        end do
                    end do
                end do
        
                call MPI_GROUP_INCL(base_group,NB(1,g)*NB(2,g)*NB(3,g),new_ranks,new_group,merror)
                call MPI_COMM_CREATE(COMM_CART,new_group,comm_temp,merror)
                call MPI_GROUP_FREE(base_group,merror)
                call MPI_GROUP_FREE(new_group ,merror)
        
                if (participate_yes(g)) then
                    call MPI_CART_CREATE(comm_temp,3,NB(1:3,g),periodic,.true.,comm1(g),merror)
                    call MPI_COMM_FREE(comm_temp,merror)
                end if
                !-----------------------------------------------------------------------------------
        
                deallocate(new_ranks)
            else
                if (g == 1) then
                    comm1          (g) = COMM_CART
                    rank_comm1     (g) = rank
                    participate_yes(g) = .true.
                else
                    comm1          (g) = comm1          (g-1)
                    rank_comm1     (g) = rank_comm1     (g-1)
                    participate_yes(g) = participate_yes(g-1)
                end if
            end if
            !--------------------------------------------------------------
     
     
            !--------------------------------------------------------------
            !--- gathering communicators (comm2) ---
            if (g >= 2 .and. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) > 1) then
        
                allocate(new_ranks(1:n_gather(1,g),1:n_gather(2,g),1:n_gather(3,g)))
        
                !-----------------------------------------------------------------------------------------------
                call MPI_COMM_GROUP(COMM_CART,base_group,merror)
        
                do kk = 1, NB(3,g)
                    do jj = 1, NB(2,g)
                        do ii = 1, NB(1,g)
                 
                            member_yes = .false.
                            do k = 1, n_gather(3,g)
                                do j = 1, n_gather(2,g)
                                    do i = 1, n_gather(1,g)
                                        ! TEST!!! Das sieht etwas komisch aus ...
                                        call MPI_CART_RANK(COMM_CART,(/MOD(((ii-1)*n_gather(1,g)+i-1)*NB1/NB(1,g-1),NB1),                             &
                                            &          MOD(((jj-1)*n_gather(2,g)+j-1)*NB2/NB(2,g-1),NB2),                             &
                                            &          MOD(((kk-1)*n_gather(3,g)+k-1)*NB3/NB(3,g-1),NB3)/),new_ranks(i,j,k),merror)
                                        if (rank == new_ranks(i,j,k)) member_yes = .true.
                                    end do
                                end do
                            end do
                 
                            call MPI_GROUP_INCL(base_group,n_gather(1,g)*n_gather(2,g)*n_gather(3,g),new_ranks,new_group,merror)
                            call MPI_COMM_CREATE(COMM_CART,new_group,comm_temp,merror)
                            call MPI_GROUP_FREE(new_group,merror)
                 
                            if (member_yes) then
                                call MPI_CART_CREATE(comm_temp,3,n_gather(1:3,g),periodic,.true.,comm2(g),merror)
                                call MPI_COMM_FREE(comm_temp,merror)
                    
                                !--- comm1 und comm2 synchronisieren ---
                                rank_comm2(g) = 0
                                if (participate_yes(g)) call MPI_COMM_RANK(comm2(g),rank_comm2(g),merror) ! Nur fuer "rankc2(g)" gilt "participate_yes(g)"
                                call MPI_ALLREDUCE(rank_comm2(g),rankc2(g),1,MPI_INTEGER,MPI_SUM,comm2(g),merror)
                                call MPI_COMM_RANK(comm2(g),rank_comm2(g),merror)
                            end if
                 
                        end do
                    end do
                end do
        
                call MPI_GROUP_FREE(base_group,merror)
                !-----------------------------------------------------------------------------------------------
        
                deallocate(new_ranks)
        
            end if
            !--------------------------------------------------------------
     
     
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (g >= 2 .and. n_gather(1,g)*n_gather(2,g)*n_gather(3,g) > 1 .and. participate_yes(g-1)) then
        
                iimax   = (NN(1,g)-1)/n_gather(1,g)+1
                jjmax   = (NN(2,g)-1)/n_gather(2,g)+1
                kkmax   = (NN(3,g)-1)/n_gather(3,g)+1
        
                iiShift = (iimax-1)*MOD(iB(1,g-1)-1,n_gather(1,g))
                jjShift = (jjmax-1)*MOD(iB(2,g-1)-1,n_gather(2,g))
                kkShift = (kkmax-1)*MOD(iB(3,g-1)-1,n_gather(3,g))
        
                call MPI_COMM_RANK(comm2(g),rank_comm2(g),merror)
        
                !--- Interpolation ---------------------------------------------------------------------
                offs_global = 0
                sizs_global = 0
        
                offs_global(1:3,rank_comm2(g)+1) = (/iiShift,jjShift,kkShift/)
                sizs_global(1:3,rank_comm2(g)+1) = (/iimax  ,jjmax  ,kkmax  /)
        
                call MPI_ALLREDUCE(offs_global,offsI(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
                call MPI_ALLREDUCE(sizs_global,sizsI(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        
                counter = 0
                do k = 1, n_gather(3,g)
                    do j = 1, n_gather(2,g)
                        do i = 1, n_gather(1,g)
                            recvI                (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &          = sizsI(1,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &           *sizsI(2,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &           *sizsI(3,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                 
                            dispI                    (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g) = counter
                            counter = counter + recvI(  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                        end do
                    end do
                end do
        
        
                !--- Restriction -----------------------------------------------------------------------
                if (ls1 ==  0 .and. (BC(1,1,g-1) == 0 .or. BC(1,1,g-1) == -1)) iimax = iimax-1
                if (ls1 == -1 .and. (BC(2,1,g-1) == 0 .or. BC(2,1,g-1) == -1)) iimax = iimax-1
        
                if (ls2 ==  0 .and. (BC(1,2,g-1) == 0 .or. BC(1,2,g-1) == -1)) jjmax = jjmax-1
                if (ls2 == -1 .and. (BC(2,2,g-1) == 0 .or. BC(2,2,g-1) == -1)) jjmax = jjmax-1
        
                if (ls3 ==  0 .and. (BC(1,3,g-1) == 0 .or. BC(1,3,g-1) == -1)) kkmax = kkmax-1
                if (ls3 == -1 .and. (BC(2,3,g-1) == 0 .or. BC(2,3,g-1) == -1)) kkmax = kkmax-1
        
                if (ls1 ==  0 .and. (BC(1,1,g-1) == 0 .or. BC(1,1,g-1) == -1)) iiShift = iiShift+1
                if (ls2 ==  0 .and. (BC(1,2,g-1) == 0 .or. BC(1,2,g-1) == -1)) jjShift = jjShift+1
                if (ls3 ==  0 .and. (BC(1,3,g-1) == 0 .or. BC(1,3,g-1) == -1)) kkShift = kkShift+1
        
                offs_global = 0
                sizs_global = 0
        
                offs_global(1:3,rank_comm2(g)+1) = (/iiShift,jjShift,kkShift/)
                sizs_global(1:3,rank_comm2(g)+1) = (/iimax  ,jjmax  ,kkmax  /)
        
                call MPI_ALLREDUCE(offs_global,offsR(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
                call MPI_ALLREDUCE(sizs_global,sizsR(1,1,g),3*n_gather(1,g)*n_gather(2,g)*n_gather(3,g),MPI_INTEGER,MPI_SUM,comm2(g),merror)
        
                counter = 0
                do k = 1, n_gather(3,g)
                    do j = 1, n_gather(2,g)
                        do i = 1, n_gather(1,g)
                            recvR                (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &          = sizsR(1,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &           *sizsR(2,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)      &
                                &           *sizsR(3,i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                 
                            dispR                    (  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g) = counter
                            counter = counter + recvR(  i+(j-1)*n_gather(1,g)+(k-1)*n_gather(1,g)*n_gather(2,g),g)
                        end do
                    end do
                end do
               !---------------------------------------------------------------------------------------
            end if
           !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        end do
        !*********************************************************************
  
        !===========================================================================================================
  
        dim1 = M1
        dim2 = M2
        dim3 = M3
  
        if (BC_1L_global == -1) dim1 = M1-1
        if (BC_2L_global == -1) dim2 = M2-1
        if (BC_3L_global == -1) dim3 = M3-1
  
        !===========================================================================================================
  
        bar1_size = 0
        bar2_size = 0
        bar3_size = 0

        bar1_offset = 0
        bar2_offset = 0
        bar3_offset = 0
  
        call MPI_COMM_RANK(COMM_BAR1,rank_bar1,merror)
        call MPI_COMM_RANK(COMM_BAR2,rank_bar2,merror)
        call MPI_COMM_RANK(COMM_BAR3,rank_bar3,merror)
  
        bar1_size(rank_bar1+1) = N1p-S1p+1
        bar2_size(rank_bar2+1) = N2p-S2p+1
        bar3_size(rank_bar3+1) = N3p-S3p+1
  
        call MPI_ALLREDUCE(bar1_size,bar1_out,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
        call MPI_ALLREDUCE(bar2_size,bar2_out,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
        call MPI_ALLREDUCE(bar3_size,bar3_out,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
        bar1_size = bar1_out
        bar2_size = bar2_out
        bar3_size = bar3_out
  
  
        bar1_offset(rank_bar1+1) = iShift + 1 + ls1
        bar2_offset(rank_bar2+1) = jShift + 1 + ls2
        bar3_offset(rank_bar3+1) = kShift + 1 + ls3
  
        call MPI_ALLREDUCE(bar1_offset,bar1_out,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
        call MPI_ALLREDUCE(bar2_offset,bar2_out,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
        call MPI_ALLREDUCE(bar3_offset,bar3_out,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
        bar1_offset = bar1_out
        bar2_offset = bar2_out
        bar3_offset = bar3_out
  
  
    end subroutine init_limits
  
  
  
end module cmod_setup
