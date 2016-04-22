!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing functions to initiliaze stencil arrays
module cmod_FDcoeff

  use iso_c_binding
  use mpi

  implicit none

contains

  !> \brief calculates finite differences coefficients.
  !!
  !! \param[in] Nmax dimension of cc and xc,xs corresponds to \c GridSizeLocal
  !! \param[in] bl  coordinate lower bound \c FieldSpace
  !! \param[in] bu  coordinate upper bound \c FieldSpace
  !! \param[in] cl  stencil lower bound \c FieldSpace
  !! \param[in] cu  stencil upper bound \c FieldSpace
  !! \param[in] BCL local lower Boundary conditions \c BoundaryConditionsLocal
  !! \param[in] BCU local upper Boundary conditions \c BoundaryConditionsLocal
  !! \param[in] SShift shift in \c ProcGrid
  !! \param[in] grid_type (5: p, 1: u, 2: v, 3:w ) used weridly
  !! \param[in] dir direction ( 1:x, 2:y, 3:z )
  !! \param[in] abl degree of derivative
  !! \param[in] upwind (0: central, 1: up, -1: low)
  !! \param[in] mapping_yes computation of finite difference coefficients (relevant only for explicit finite differences), note affects the accuracy and temporal stability of the discretization
  !!            -mapping_yes = .TRUE.  (computation on computational, equidistant grid and subsequent mapping on
  !!                        physical grid)
  !!            -mapping_yes = .FALSE. (direct computation on physical grid)
  !! \param[in] dim_ncb dimension of ??
  !! \param[in] n_coeff_bound ??
  !! \param[in] xC coordinate in of varaiable
  !! \param[in] xE coordinate out of equations
  !! \param[out] cc coefficients
  !!
  !! todo move error handling else where
  !!
  !! \note: - Upwinding wird auf dem Rand unterdrücken, da Differenzenstencils dort ohnehin schief sind.
  !!        - Generell koennte/sollte die Initialisierung der Index-Grenzen auch schon vor dieser       
  !!          Routine ausgefuehrt werden, um hier die Uebersichtlichkeit zu verbessern.                 
  subroutine FD_getDiffCoeff( &
      Nmax,                   &
      bL,                     &
      bU,                     &
      cL,                     &
      cU,                     &
      BCL,                    &
      BCU,                    &
      SShift,                 &
      grid_type,              &
      dir,                    &
      abl,                    &
      upwind,                 &
      mapping_yes,            &
      dim_ncb,                &
      n_coeff_bound,          &
      xC,                     &
      xE,                     &
      cc ) bind(c,name='FD_getDiffCoeff')

    implicit none


    integer(c_int),  intent(in)   :: Nmax

    integer(c_int),  intent(in)   :: bL
    integer(c_int),  intent(in)   :: bU

    integer(c_int),  intent(in)   :: cL
    integer(c_int),  intent(in)   :: cU

    integer(c_int),  intent(in)   :: BCL
    integer(c_int),  intent(in)   :: BCU

    integer(c_int),  intent(in)   :: SShift

    integer(c_int),  intent(in)   :: grid_type

    integer(c_int),  intent(in)   :: dir
    integer(c_int),  intent(in)   :: abl

    integer(c_int),  intent(in)   :: upwind
    logical(c_bool), intent(in)   :: mapping_yes

    integer(c_int),  intent(in)   :: dim_ncb
    integer(c_int),  intent(in)   :: n_coeff_bound(1:dim_ncb)

    real(c_double),  intent(in)   :: xC(bL:(Nmax+bU))
    real(c_double),  intent(in)   :: xE(bL:(Nmax+bU))

    real(c_double),  intent(out)  :: cc(cL:cU,0:Nmax)

    integer(c_int)        :: n_coeff
    integer(c_int)        :: dim_n_coeff_bound
    integer(c_int)        :: i, ii, iC, iStart
    integer(c_int)        :: k!, kk

    integer(c_int)        :: left, right

    real(c_double)        :: dxi(1:2)

    real(c_double)        :: dxL, dxU ! Für Integrationskoeffizienten

    real(c_double), allocatable  :: cc_xi(:,:)
    real(c_double), allocatable  :: deltaX(:)

    logical(c_bool)              :: filter_yes

    integer(c_int)               :: merror


    if (mapping_yes .and. abl > 2) then
      !if (rank == 0)
      write(*,*) 'ERROR! Can`t handle derivatives > 2 combined with mapping ...'
      !call MPI_FINALIZE(merror)
      !stop
    end if

    !===========================================================================================
    !=== Startindex für Entwicklungspunkte =====================================================
    !===========================================================================================
    if( grid_type == 5 .or. grid_type == 3 ) then
      iStart = 1
    else
      iStart = 0
    end if
    !===========================================================================================

    filter_yes = .false.
    if( abl==0 .and. (grid_type==5 .or. grid_type==1) ) filter_yes = .true.

    dim_n_coeff_bound = SIZE(n_coeff_bound)

    cc = 0.

    right = 0 ! to please compiler

    do i = iStart, Nmax

      !=========================================================================================
      !=== Stencil-Breite auslesen =============================================================
      !=========================================================================================
      if( BCL>0 .and. i<=( iStart - 1 + dim_n_coeff_bound ) ) then
        n_coeff = n_coeff_bound(i + 1 - iStart)
        !IF (upwind /= 0 .AND. i == iStart) n_coeff = 0
      else if( BCU>0 .and. i >=( Nmax + 1 - dim_n_coeff_bound ) ) then
        n_coeff = n_coeff_bound(Nmax + 1 - i)
        !IF (upwind /= 0 .AND. i == Nmax  ) n_coeff = 0
      else
        n_coeff = n_coeff_bound(dim_n_coeff_bound)
      end if

      !=========================================================================================

      if (n_coeff > 0 .and. n_coeff > abl) then

        !=======================================================================================
        !=== Anzahl der Koeffizienten RECHTS vom Entwicklungspunkt bestimmen ===================
        !=======================================================================================
        if (grid_type == 5 .or. grid_type == 1) then
          !=====================================================================================
          !=== "normales" Gitter, keine Zwischengitterpunkte ===================================
          !=====================================================================================
          if      (BCL > 0 .and. i <= (iStart - 1 + n_coeff/2)) then
            right = n_coeff - i + iStart - 1
          else if (BCU > 0 .and. i >= (Nmax   + 1 - n_coeff/2)) then
            right = Nmax - i
          else
            ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= ungerade Anzahl Koeffizienten)!
            if (MOD(n_coeff,2) == 0) then
              !if (rank == 0) then
              write(*,'(a   )') 'ERROR! Choose odd number of coefficients!'
              write(*,'(a,i4)') '    direction =', dir
              write(*,'(a,i4)') '    grid_type =', grid_type
              write(*,'(a,i4)') '            i =', i+SShift
              !end if
            end if
            right = (n_coeff-1)/2
          end if
        else if( grid_type==2 ) then
          !=====================================================================================
          !=== Druck ==> Impuls ================================================================
          !=====================================================================================
          if( BCL>0 .and. i<=n_coeff/2 ) then
            right = n_coeff - i - iStart
          else if( BCU>0 .and. i>=( Nmax + 0 - n_coeff/2 ) ) then
            right = Nmax - i
          else
            ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
            if (MOD(n_coeff,2) /= 0) then
              !if (rank == 0) then
              write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
              write(*,'(a,i4)') '    direction =', dir
              write(*,'(a,i4)') '    grid_type =', grid_type
              write(*,'(a,i4)') '            i =', i+SShift
              !end if
              call MPI_FINALIZE(merror)
              stop
            end if
            right = n_coeff/2
          end if
        else if( grid_type==3 ) then
          !=====================================================================================
          !=== Geschwindigkeit ==> Konti =======================================================
          !=====================================================================================
          if( BCL>0 .and. i <= n_coeff/2 ) then
            right = n_coeff - i - iStart
          else if( BCU>0 .and. i>=( Nmax + 1 - n_coeff/2) ) then
            right = Nmax - i
          else
            ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
            if (MOD(n_coeff,2) /= 0) then
              !if (rank == 0) then
              write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
              write(*,'(a,i4)') '    direction =', dir
              write(*,'(a,i4)') '    grid_type =', grid_type
              write(*,'(a,i4)') '            i =', i+SShift
              !end if
              call MPI_FINALIZE(merror)
              stop
            end if
            right = n_coeff/2 - 1
          end if
        end if
        !=======================================================================================


        !=======================================================================================
        !=== Anzahl der Koeffizienten LINKS vom Entwicklungspunkt bestimmen ====================
        !=======================================================================================
        ! (0-ter bzw. 1-ter Koeffizient links vom Entwicklungspunkt
        !          == "zentraler" Koeffizient im Speicher)
        left = right - n_coeff + 1
        !=======================================================================================


        !=======================================================================================
        if( upwind == -1 ) then
          if( .not. (BCL > 0 .and. i < (iStart + n_coeff/2)) ) then
            n_coeff = n_coeff - 1
            left    = left    + 1
          end if
        end if
        if( upwind ==  1 ) then
          if( .not. (BCU > 0 .and. i > (Nmax   - n_coeff/2)) ) then
            n_coeff = n_coeff - 1
            right   = right   - 1
          end if
        end if
        !=======================================================================================


        !=======================================================================================
        !=== Stencilanordnung testen ===========================================================
        !=======================================================================================
        if( right > cU .or. left < cL ) then
          !if (rank == 0) then
          WRITE(*,'(a   )') 'WARNING! The FD-Stencil does probably not fit into provided array!'
          write(*,'(a   )') 'ERROR! Stencil doesn`t fit into provided array!'
          write(*,'(a,i4)') '    direction =', dir
          write(*,'(a,i4)') '    grid_type =', grid_type
          write(*,'(a,i4)') '            i =', i+SShift
          !end if
          call MPI_FINALIZE(merror)
          stop
        end if
        !=======================================================================================


        allocate(deltaX(1:n_coeff))


        !=======================================================================================
        !=== räumliche Abstände zum Entwicklungspunkt =========================================+
        !=======================================================================================
        do ii = 1, n_coeff
          ! Koeffizienten-Punkte ("iC"):
          iC = i + left + (ii-1)

          if (mapping_yes) then
            if( grid_type == 2 ) then
              deltaX(ii) = REAL(iC-i)-0.5
            else if (grid_type == 3) then
              deltaX(ii) = REAL(iC-i)+0.5
            else
              deltaX(ii) = REAL(iC-i)
            end if
          else
            deltaX(ii) = xC(iC) - xE(i)
          end if
        end do
        !=======================================================================================



        !=======================================================================================
        !=== Integrations-Intervall ============================================================
        !=======================================================================================
        ! Anmerkungen: - Kein Mapping aus Genauigkeitsgründen vorgesehen.
        !              - Nur für Druckgitter vorgesehen.
        if (abl == -1) then

          if (BCL > 0 .and. i == iStart) then
            dxL = 0.
          else
            dxL = (xC(i-1) - xE(i)) / 2. ! der Einfachheit halber anstelle der Geschwindigkeitspunkte ...
          end if
          if (BCU > 0 .and. i == Nmax  ) then
            dxU = 0.
          else
            dxU = (xC(i+1) - xE(i)) / 2.
          end if

        end if
        !=======================================================================================


        !=======================================================================================
        !=== Bestimmung der Koeffizienten ======================================================
        !=======================================================================================
        ! Anmerkung: Filter und Integratoren sollen nicht gemappt werden (siehe z.B. Stolz).
        if (mapping_yes .and. abl /= 0 .and. abl /= -1) then

          !--- derivatives (mapped) ---

          allocate(cc_xi(left:right,1:abl))

          cc_xi = 0.
          dxi   = 0.

          do k = 1, abl
            if (n_coeff <= 7) then ! TEST!!!
              !kk = k
              !INCLUDE 'FD_coeffs_expl_map.f90'
              call diff_coeffs_exact(k,n_coeff,deltaX(1),cc_xi(left:right,k)) ! TEST!!!
            else
              call FD_coeffs_solver(k,filter_yes,n_coeff,deltaX,cc_xi(left:right,k))
            end if

            do ii = left, right
              dxi(k) = dxi(k) + cc_xi(ii,k)*xC(i+ii)
            end do
          end do

          !---------------------------------------------------!
          ! Mapping:                                          !
          ! d/dx     = (xi') * d/dxi                          !
          ! d^2/dx^2 = (xi") * d/dxi + (xi')**2 * d^2/dxi^2   !
          !                                                   !
          ! Mapping-Faktoren:                                 !
          ! xi' =   1  /  x'                                  !
          ! xi" = - x" / (x')**3                              !
          !---------------------------------------------------!

          if (abl == 1) cc(left:right,i) = cc_xi(left:right,1)/dxi(1)
          if (abl == 2) cc(left:right,i) = cc_xi(left:right,2)/dxi(1)**2 - cc_xi(left:right,1)*dxi(2)/dxi(1)**3

          deallocate(cc_xi)

        else if( mapping_yes .and. abl==0 .and. .not. filter_yes ) then ! TEST!!!

          !--- interpolation (mapped) ---

          !k  = 1
          !kk = 0
          !ALLOCATE(cc_xi(left:right,k:k))
          !INCLUDE 'FD_coeffs_expl_map.f90'
          !cc(left:right,i) = cc_xi(left:right,k)
          !DEALLOCATE(cc_xi)

          if( n_coeff <= 7 ) then ! TEST!!!
            call diff_coeffs_exact(abl,n_coeff,deltaX(1),cc(left:right,i))
          else
            call FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
          end if

        else

          if (abl == -1) then
            !--- integration coefficients ---
            call FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc(left:right,i))
          else
            !--- interpolation and derivatives ---
            call FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
          end if

        end if
        !=======================================================================================


        deallocate(deltaX)


        !=======================================================================================
        !=== Symmetrie =========================================================================
        !=======================================================================================
        if (BCL == -2 .and. abl == -1 .and. i == iStart) then
          cc(0,i) = 0.5*cc(0,i)
          do ii = cL, -1
            cc(ii,i) = 0.
          end do
        end if
        if (BCU == -2 .and. abl == -1 .and. i == Nmax  ) then
          cc(0,i) = 0.5*cc(0,i)
          do ii = 1, cU
            cc(ii,i) = 0.
          end do
        end if
        !---------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------
        ! TEST!!! abl == -1 ???
        if (BCL == -2 .and. (i+cL) < 1) then
          if (grid_type == 5 .or. (grid_type == 2 .and. i >= 1)) then
            do ii = 1, 1-(i+cL)
              cc(1-i+ii,i) = cc(1-i+ii,i) + cc(1-i-ii,i)
              cc(1-i-ii,i) = 0.
            end do
          end if
          if ((grid_type == 1 .and. i >= 1) .or. grid_type == 3) then
            do ii = 1, 1-(i+cL)
              cc(0-i+ii,i) = cc(0-i+ii,i) - cc(1-i-ii,i)
              cc(1-i-ii,i) = 0.
            end do
          end if
        end if
        !---------------------------------------------------------------------------------------
        if (BCU == -2 .and. (i+cU) > (Nmax+0)) then
          if (grid_type == 5 .or. (grid_type == 2 .and. i <= Nmax-1)) then
            do ii = 1, i+cU-(Nmax-0)
              cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) + cc(Nmax-i  +ii,i)
              cc(Nmax-i  +ii,i) = 0.
            end do
          end if
        end if
        if (BCU == -2 .and. (i+cU) > (Nmax-1)) then
          if ((grid_type == 1 .and. i <= Nmax-1) .or. grid_type == 3) then
            do ii = 1, i+cU-(Nmax-1)
              cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) - cc(Nmax-i-1+ii,i)
              cc(Nmax-i-1+ii,i) = 0.
            end do
          end if
        end if
        !=======================================================================================

      end if

    end do


  end subroutine FD_getDiffCoeff





  subroutine diff_coeffs_exact( abl, n_coeff, deltaX, cc ) ! TEST!!!

    implicit none

    integer(c_int), intent(in   ) ::  abl
    integer(c_int), intent(in   ) ::  n_coeff
    real(c_double), intent(in   ) ::  deltaX
    real(c_double), intent(out  ) ::  cc(1:n_coeff)

    ! TEST!!!
    ! - korrekt ausgerichtet?
    ! - Vorzeichen ok?


    if (n_coeff == 2) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,-1./)/2.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 1., 1./)/2.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/-1., 3./)/2.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.

      if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
    end if
    if (n_coeff == 3) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/15.,-10., 3./)/8.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,  6.,-1./)/8.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/-1.,  6., 3./)/8.
      if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,-10.,15./)/8.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-2.,  3.,-1./)/1.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1.,  1., 0./)/1.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/ 0., -1., 1./)/1.
      if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/ 1., -3., 2./)/1.

      if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-3.,  4.,-1./)/2.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/-1.,  0., 1./)/2.
      if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/ 1., -4., 3./)/2.

      if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
      if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
      if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
    end if
    if (n_coeff == 4) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/ 35.,-35.,  21., -5./)/16.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/  5., 15.,  -5.,  1./)/16.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/ -1.,  9.,   9., -1./)/16.
      if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/  1., -5.,  15.,  5./)/16.
      if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/ -5., 21., -35., 35./)/16.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-71.,141., -93., 23./)/24.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-23., 21.,   3., -1./)/24.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/  1.,-27.,  27., -1./)/24.
      if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/  1., -3., -21., 23./)/24.
      if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/-23., 93.,-141., 71./)/24.

      if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-11., 18.,  -9.,  2./)/6.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/ -2., -3.,   6., -1./)/6.
      if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/  1., -6.,   3.,  2./)/6.
      if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/ -2.,  9., -18., 11./)/6.

      if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/  2., -5.,   4., -1./)/1.
      if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/  1., -2.,   1.,  0./)/1.
      if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/  0.,  1.,  -2.,  1./)/1.
      if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,  4.,  -5.,  2./)/1.
    end if
    if (n_coeff == 5) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/315.,-420., 378.,-180., 35./)/128.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 35., 140., -70.,  28., -5./)/128.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/ -5.,  60.,  90., -20.,  3./)/128.
      if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/  3., -20.,  90.,  60., -5./)/128.
      if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/ -5.,  28., -70., 140., 35./)/128.
      if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/ 35.,-180., 378.,-420.,315./)/128.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-93., 229.,-225., 111.,-22./)/24.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-22.,  17.,   9.,  -5.,  1./)/24.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/  1., -27.,  27.,  -1.,  0./)/24.
      if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/  0.,   1., -27.,  27., -1./)/24.
      if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/ -1.,   5.,  -9., -17., 22./)/24.
      if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/ 22.,-111., 225.,-229., 93./)/24.

      if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-25.,  48., -36.,  16., -3./)/12.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/ -3., -10.,  18.,  -6.,  1./)/12.
      if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/  1.,  -8.,   0.,   8., -1./)/12.
      if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/ -1.,   6., -18.,  10.,  3./)/12.
      if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/  3., -16.,  36., -48., 25./)/12.

      if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/ 35.,-104., 114., -56., 11./)/12.
      if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/ 11., -20.,   6.,   4., -1./)/12.
      if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,  16., -30.,  16., -1./)/12.
      if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,   4.,   6., -20., 11./)/12.
      if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/ 11., -56., 114.,-104., 35./)/12.
    end if
    if (n_coeff == 6) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/  693.,-1155.,  1386., -990.,   385., -63./)/256.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/   63.,  315.,  -210.,  126.,   -45.,   7./)/256.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/   -7.,  105.,   210.,  -70.,    21.,  -3./)/256.
      if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/    3.,  -25.,   150.,  150.,   -25.,   3./)/256.
      if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/   -3.,   21.,   -70.,  210.,   105.,  -7./)/256.
      if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/    7.,  -45.,   126., -210.,   315.,  63./)/256.
      if (deltaX == -5.5 .and. abl == 0) cc(1:n_coeff) = (/  -63.,  385.,  -990., 1386., -1155., 693./)/256.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-9129.,26765.,-34890.,25770.,-10205.,1689./)/1920.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1689., 1005.,  1430.,-1110.,   435., -71./)/1920.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/   71.,-2115.,  2070.,   10.,   -45.,   9./)/1920.
      if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/   -9.,  125., -2250., 2250.,  -125.,   9./)/1920.
      if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/   -9.,   45.,   -10.,-2070.,  2115., -71./)/1920.
      if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/   71., -435.,  1110.,-1430., -1005.,1689./)/1920.
      if (deltaX == -5.5 .and. abl == 1) cc(1:n_coeff) = (/-1689.,10205.,-25770.,34890.,-26765.,9129./)/1920.

      if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/ -137.,  300.,  -300.,  200.,   -75.,  12./)/60.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/  -12.,  -65.,   120.,  -60.,    20.,  -3./)/60.
      if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/    3.,  -30.,   -20.,   60.,   -15.,   2./)/60.
      if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/   -2.,   15.,   -60.,   20.,    30.,  -3./)/60.
      if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/    3.,  -20.,    60., -120.,    65.,  12./)/60.
      if (deltaX == -5.0 .and. abl == 1) cc(1:n_coeff) = (/  -12.,   75.,  -200.,  300.,  -300., 137./)/60.

      if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/   45., -154.,   214., -156.,    61., -10./)/12.
      if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/   10.,  -15.,    -4.,   14.,    -6.,   1./)/12.
      if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/   -1.,   16.,   -30.,   16.,    -1.,   0./)/12.
      if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/    0.,   -1.,    16.,  -30.,    16.,  -1./)/12.
      if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/    1.,   -6.,    14.,   -4.,   -15.,  10./)/12.
      if (deltaX == -5.0 .and. abl == 2) cc(1:n_coeff) = (/  -10.,   61.,  -156.,  214.,  -154.,  45./)/12.
    end if
    if (n_coeff == 7) then
      if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/  3003., -6006.,  9009.,-8580.,   5005., -1638.,  231./)/1024.
      if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/   231.,  1386., -1155.,  924.,   -495.,   154.,  -21./)/1024.
      if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/   -21.,   378.,   945., -420.,    189.,   -54.,    7./)/1024.
      if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/     7.,   -70.,   525.,  700.,   -175.,    42.,   -5./)/1024.
      if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/    -5.,    42.,  -175.,  700.,    525.,   -70.,    7./)/1024.
      if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/     7.,   -54.,   189., -420.,    945.,   378.,  -21./)/1024.
      if (deltaX == -5.5 .and. abl == 0) cc(1:n_coeff) = (/   -21.,   154.,  -495.,  924.,  -1155.,  1386.,  231./)/1024.
      if (deltaX == -6.5 .and. abl == 0) cc(1:n_coeff) = (/   231., -1638.,  5005.,-8580.,   9009., -6006., 3003./)/1024.

      if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-10756., 36527.,-59295., 58310.,-34610., 11451.,-1627./)/1920.
      if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/ -1627.,   633.,  2360., -2350.,  1365.,  -443.,   62./)/1920.
      if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/    62., -2061.,  1935.,   190.,  -180.,    63.,   -9./)/1920.
      if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/    -9.,   125., -2250.,  2250.,  -125.,     9.,    0./)/1920.
      if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/     0.,    -9.,   125., -2250.,  2250.,  -125.,    9./)/1920.
      if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/     9.,   -63.,   180.,  -190., -1935.,  2061.,  -62./)/1920.
      if (deltaX == -5.5 .and. abl == 1) cc(1:n_coeff) = (/   -62.,   443., -1365.,  2350., -2360.,  -633., 1627./)/1920.
      if (deltaX == -6.5 .and. abl == 1) cc(1:n_coeff) = (/  1627.,-11451., 34610.,-58310., 59295.,-36527.,10756./)/1920.

      if (deltaX == -0.0 .and. abl == 1) cc(1:n_coeff) = (/  -147.,   360.,  -450.,   400.,  -225.,    72.,  -10./)/60.
      if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/   -10.,   -77.,   150.,  -100.,    50.,   -15.,    2./)/60.
      if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/     2.,   -24.,   -35.,    80.,   -30.,     8.,   -1./)/60.
      if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/    -1.,     9.,   -45.,     0.,    45.,    -9.,    1./)/60.
      if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/     1.,    -8.,    30.,   -80.,    35.,    24.,   -2./)/60.
      if (deltaX == -5.0 .and. abl == 1) cc(1:n_coeff) = (/    -2.,    15.,   -50.,   100.,  -150.,    77.,   10./)/60.
      if (deltaX == -6.0 .and. abl == 1) cc(1:n_coeff) = (/    10.,   -72.,   225.,  -400.,   450.,  -360.,  147./)/60.

      if (deltaX == -0.0 .and. abl == 2) cc(1:n_coeff) = (/   812., -3132.,  5265., -5080.,  2970.,  -972.,  137./)/180.
      if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/   137.,  -147.,  -225.,   470.,  -285.,    93.,  -13./)/180.
      if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/   -13.,   228.,  -420.,   200.,    15.,   -12.,    2./)/180.
      if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/     2.,   -27.,   270.,  -490.,   270.,   -27.,    2./)/180.
      if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/     2.,   -12.,    15.,   200.,  -420.,   228.,  -13./)/180.
      if (deltaX == -5.0 .and. abl == 2) cc(1:n_coeff) = (/   -13.,    93.,  -285.,   470.,  -225.,  -147.,  137./)/180.
      if (deltaX == -6.0 .and. abl == 2) cc(1:n_coeff) = (/   137.,  -972.,  2970., -5080.,  5265., -3132.,  812./)/180.
    end if


  end subroutine diff_coeffs_exact











  !    subroutine test_diff(dir,abl,xC,xE,Nmax,bL,bU,cc,BCL,grid_type,name)
  !
  !        implicit none
  !
  !        integer, intent(in)   ::  dir
  !        integer, intent(in)   ::  abl
  !        integer, intent(in)   ::  Nmax
  !        integer, intent(in)   ::  bL
  !        integer, intent(in)   ::  bU
  !        real   , intent(in)   ::  xC(bL:(Nmax+bU))
  !        real   , intent(in)   ::  xE(bL:(Nmax+bU))
  !        real   , intent(in)   ::  cc(bL:bU,0:Nmax)
  !        integer, intent(in)   ::  BCL
  !        integer, intent(in)   ::  grid_type
  !
  !        character(*), intent(in) ::  name
  !
  !        integer               ::  i, ii, iStartC, iStartE, SShift
  !        integer               ::  nw, nw_max, n_dx
  !        real                  ::  wave, wave_mod(1:2), dx
  !
  !        character(len=1)      ::  part
  !        real                  ::  pi
  !
  !
  !        !----------------------------------------------------------------------------------------------------------!
  !        ! Anmerkungen: - Hier sollen alle Differenzen-Stencils innerhalb eines Blockes getestet werden, analog zu  !
  !        !                ihrer Berechnung. Daher wird bewusst nicht mittels MPI in eine Datei geschrieben, sondern !
  !        !                in verschiedene Blöcke, um Fehler bei der MPI-Programmierung auszuschliessen.             !
  !        !----------------------------------------------------------------------------------------------------------!
  !
  !
  !        !--- pi ---
  !        pi = 2.*ABS(ACOS(0.))
  !
  !
  !        if (dir == 1) then
  !            write(part,'(i1.1)') iB(1,1)
  !            SShift = iShift
  !        else if (dir == 2) then
  !            write(part,'(i1.1)') iB(2,1)
  !            SShift = jShift
  !        else if (dir == 3) then
  !            write(part,'(i1.1)') iB(3,1)
  !            SShift = kShift
  !        else
  !            if (rank == 0) write(*,*) 'ERROR! Wrong input at subroutine `test_diff`!'
  !            call MPI_FINALIZE(merror)
  !            stop
  !        end if
  !
  !
  !        open(10,file='test_'//name//'_transfer_block'//part//'.txt',status='UNKNOWN')
  !        open(11,file='test_'//name//'_coeffs_block'  //part//'.txt',status='UNKNOWN')
  !
  !
  !        !===========================================================================================================
  !        !=== Startindex für Entwicklungspunkte =====================================================================
  !        !===========================================================================================================
  !        if (BCL <= 0 .or. grid_type == 0 .or. grid_type == 3) then
  !            iStartE = 1
  !        else
  !            iStartE = 0
  !        end if
  !
  !
  !        !===========================================================================================================
  !        !=== Startindex für Koeffizienten ==========================================================================
  !        !===========================================================================================================
  !        if (BCL > 0 .and. (grid_type == 1 .or. grid_type == 3)) then
  !            iStartC = 0
  !        else
  !            iStartC = 1
  !        end if
  !
  !
  !        !===========================================================================================================
  !        !=== Test der Koeffizienten ================================================================================
  !        !===========================================================================================================
  !        nw_max = 100
  !
  !        do i = iStartE, Nmax
  !
  !            do nw = 1, nw_max
  !
  !                wave     = pi*REAL(nw)/REAL(nw_max)
  !                wave_mod = 0.
  !
  !                !=== Referenz-Gitterweite bestimmen ==================================================================
  !                dx   = 0.
  !                n_dx = 0
  !
  !                do ii = bL, bU-1
  !                    !IF (i+ii .GE. iStartC .AND. i+ii .LT. Nmax .AND. ABS(xC(i+ii+1)-xC(i+ii)) .GT. dx) THEN
  !                    if (i+ii >= iStartC .and. i+ii < Nmax) then
  !                        dx   = dx + ABS(xC(i+ii+1)-xC(i+ii))
  !                        n_dx = n_dx + 1
  !                    end if
  !                end do
  !
  !                dx = dx / REAL(n_dx)
  !
  !
  !                !=== Transferfunktion plotten ========================================================================
  !                if (abl == 0) then
  !                    do ii = bL, bU
  !                        wave_mod(1) = wave_mod(1) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))
  !                        wave_mod(2) = wave_mod(2) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))
  !                    end do
  !                else if (abl == 1) then
  !                    do ii = bL, bU
  !                        wave_mod(1) = wave_mod(1) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx
  !                        wave_mod(2) = wave_mod(2) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx/wave**2
  !                    end do
  !                else if (abl == 2) then
  !                    do ii = bL, bU
  !                        wave_mod(1) = wave_mod(1) - cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx**2
  !                        wave_mod(2) = wave_mod(2) - cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx**2/wave
  !                    end do
  !                end if
  !
  !                write(10,'(i7,3E26.17)') i+SShift, wave, wave_mod(1:2)
  !
  !            end do
  !
  !            write(10,*)
  !            write(11,'(i7,100E26.17)') i+SShift, cc(:,i)
  !
  !        end do
  !
  !        close(10)
  !        close(11)
  !
  !
  !    end subroutine test_diff











  subroutine FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc)

    implicit none

    !integer(c_int) , intent(in)   ::  rank
    integer(c_int) , intent(in)   ::  abl
    logical(c_bool), intent(in)   ::  filter_yes
    integer(c_int) , intent(in)   ::  n_coeff
    real(c_double) , intent(in)   ::  deltaX(1:n_coeff)
    real(c_double) , intent(out)  ::  cc    (1:n_coeff)

    integer(c_int)               ::  i, j

    real(c_double)                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
    real(c_double)                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)

    real(c_double)                  ::  const


    !===========================================================================================
    !=== Aufstellen des Gleichungssystems ======================================================
    !===========================================================================================
    do i = 1, n_coeff
      do j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
      end do
    end do


    !===========================================================================================
    !=== Zusatzbedingungen =====================================================================
    !===========================================================================================
    if (filter_yes) then
      ! G(pi) = 0.
      do i = 1, n_coeff
        polyn_vals(i,n_coeff) = (-1.)**i
      end do
    end if


    !===========================================================================================
    !=== Lösen des Gleichungssystems ===========================================================
    !===========================================================================================
    call Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)


    !===========================================================================================
    !=== Funktionswert am Entwicklungspunkt (deltaX = 0) =======================================
    !==========================================================================================
    const = 1.
    do i = 1, abl-1
      const = const*REAL(1+i)
    end do


    !===========================================================================================
    !=== Koeffizienten bestimmen (explizite Differenzen) =======================================
    !===========================================================================================
    cc(1:n_coeff) = const*polyn_vals_inv(1+abl,1:n_coeff)


  end subroutine FD_coeffs_solver











  subroutine FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc)

    implicit none

    !integer(c_int), intent(in)   ::  rank

    integer(c_int), intent(in)   ::  n_coeff
    real(c_double), intent(in)   ::  deltaX(1:n_coeff)
    real(c_double), intent(out)  ::  cc    (1:n_coeff)
    real(c_double), intent(in)   ::  dxL, dxU

    integer(c_int)               ::  i, j!, k

    real(c_double)                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
    real(c_double)                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)

    !real(c_double)                  ::  const


    !===========================================================================================
    !=== Aufstellen des Gleichungssystems ======================================================
    !===========================================================================================
    do i = 1, n_coeff
      do j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
      end do
    end do


    !===========================================================================================
    !=== Lösen des Gleichungssystems ===========================================================
    !===========================================================================================
    call Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)


    !===========================================================================================
    !=== Koeffizienten bestimmen (explizite Differenzen) =======================================
    !===========================================================================================
    !  (Matrix-Vektor-Multiplikation)
    cc = 0.
    do j = 1, n_coeff
      do i = 1, n_coeff
        cc(j) = cc(j) + (dxU**i - dxL**i)/REAL(i)*polyn_vals_inv(i,j)
      end do
    end do

  end subroutine FD_coeffs_solver_integral



  subroutine Matrix_invert(N,matrix,matrix_inv)

    implicit none

    !integer(c_int), intent(in)  ::  rank
    integer(c_int), intent(in)  ::  N
    real(c_double), intent(in)  ::  matrix     (1:N,1:N)
    real(c_double), intent(out) ::  matrix_inv (1:N,1:N)

    real(c_double)              ::  matrix_left(1:N,1:N)
    real(c_double)              ::  mult1, mult2
    real(c_double)              ::  eps
    integer(c_int)              ::  i, j, k

    real(c_double)              ::  store
    integer(c_int)              ::  maxValue, maxValuePos



    eps = 10.**(-20) ! double precision
    !eps = 10.**(-??) ! single precision

    matrix_left = matrix


    ! Rechte Seite initialisieren (= Einheitsmatrix):
    matrix_inv = 0.

    do i = 1, N
      matrix_inv(i,i) = 1.
    end do


    !--- Vorwaertsschritt -------------------------------------
    ! linke Seite umformen in obere Dreiecksmatrix
    ! rechte Seite umformen in untere Dreiecksmatrix

    do j = 1, N-1

      ! Pivoting == Umsortieren der aktuellen Untermatrix (j:N,j:N)
      ! (Diagonalelemente der linken Dreiecksmatrix sind betragsmaessig zu maximieren)
      ! Groesster Wert in Spalte j = zukuenftiges Diagonalelement j,j
      maxValue    = int( abs(matrix_left(j,j)), c_int )
      maxValuePos = j
      do i = j+1, N
        if (ABS(matrix_left(i,j)) > maxValue) then
          maxValue    = int( ABS(matrix_left(i,j)), c_int )
          maxValuePos = i
        end if
      end do

      ! Zeilen vertauschen:
      if (maxValuePos /= j) then
        do i = 1, N
          store                      = matrix_left(maxValuePos,i)
          matrix_left(maxValuePos,i) = matrix_left(j,i)
          matrix_left(j,i)           = store

          store                      = matrix_inv(maxValuePos,i)
          matrix_inv(maxValuePos,i)  = matrix_inv(j,i)
          matrix_inv(j,i)            = store
        end do
      end if

      mult1 = matrix_left(j,j)
      if (ABS(mult1) < eps ) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (1)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', j
      end if

      do i = j+1, N
        mult2 = matrix_left(i,j)/mult1

        ! Aufaddieren von Zeile j auf aktuelle Zeile i (linke Seite):
        do k = j, N
          ! To avoid underflow:
          if (.not. (ABS(mult2) < eps .and. ABS(matrix_left(j,k)) < eps)) then
            matrix_left(i,k) = matrix_left(i,k) - matrix_left(j,k)*mult2
          end if
        end do

        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        do k = 1, N
          ! To avoid underflow:
          if (.not. (ABS(mult2) < eps .and. ABS(matrix_inv(j,k)) < eps)) then
            matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
          end if
        end do

        ! Komponente i,j explizit zu Null setzen:
        matrix_left(i,j) = 0.
      end do

    end do


    !--- Rueckwaertsschritt -----------------------------------
    ! linke Seite umformen in Einheitsmatrix
    ! rechte Seite umformen in gesuchte inverse Matrix

    do j = N, 2, -1

      ! Multiplikator:
      mult1 = matrix_left(j,j)
      if (ABS(mult1) < eps ) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (2)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', j
      end if

      do i = 1, j-1
        mult2 = matrix_left(i,j)/mult1

        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        do k = 1, N
          ! To avoid underflow:
          if (.not. (ABS(mult2) < eps .and. ABS(matrix_inv(j,k)) < eps)) then
            matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
          end if
        end do

        ! nicht-Diagonalelement explizit zu 0 setzen:
        matrix_left(i,j) = 0.
      end do

    end do


    ! linke Matrix auf Einheitsmatrix bringen:
    do i = 1, N
      mult1 = matrix_left(i,i)
      if (ABS(mult1) < eps ) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (3)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', i
      end if
      matrix_inv(i,:) = matrix_inv(i,:)/mult1
    end do


  end subroutine Matrix_invert


end module cmod_FDcoeff
