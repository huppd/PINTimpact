!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing functions to initiliaze stencil arrays
module cmod_FDcoeff

    use iso_c_binding
    use mpi

    private

    public  ::  FD_getDiffCoeff
  
contains
  
    !!pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
  
  
  
  
!    subroutine myDGBTRF(N,ndL,AB,LDAB,pivot)
!
!        implicit none
!
!        integer, intent(in   )  ::  ndL, LDAB, N
!        integer, intent(out  )  ::  pivot( * )
!        real   , intent(inout)  ::  AB( LDAB, * )
!        integer, parameter      ::  NBMAX = 64
!        integer, parameter      ::  LDWORK = NBMAX+1
!        integer                 ::  I, I2, I3, II, IP, J, J2, J3, JJ, JM, JP, JU, K2, KM, KV, NW
!        real                    ::  TEMP
!        real                    ::  WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX )
!        integer                 ::  IDAMAX
!
!        external IDAMAX
!        external DCOPY, DGEMM, DGER, DLASWP, DSCAL, DSWAP, DTRSM
!
!
!        KV = ndL + ndL
!
!        do J = ndL + 2, MIN( KV, N )
!            do I = KV - J + 2, ndL
!                AB( I, J ) = 0.
!            end do
!        end do
!
!        JU = 1
!
!        do J = 1, N
!            I2 = MIN(ndL-1,N-J)
!            I3 = MIN(1,N-J-ndL+1)
!
!            if(J+KV <= N) then
!                do I = 1, ndL
!                    AB(I,J+KV) = 0.
!                end do
!            end if
!
!            KM = MIN(ndL,N-J)
!
!            !--- Pivoting ---
!            !JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
!            JP = 1 ! TEST!!!
!
!            pivot(J) = JP
!
!            if (AB(KV+JP,J) /= 0.) then
!                JU = MAX( JU, MIN(J+ndL+JP-1,N))
!                if (JP /= 1) then
!                    if (JP-1 < ndL) then
!                        call DSWAP( 1, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
!                    else
!                        call DSWAP( 0, AB( KV+1, J ), LDAB-1,WORK31( JP-ndL, 1 ), LDWORK )
!                        call DSWAP( 1, AB( KV+1, J ), LDAB-1,AB( KV+JP, J ), LDAB-1 )
!                    end if
!                end if
!                call DSCAL( KM, 1. / AB( KV+1, J ), AB( KV+2, J ),1 )
!                JM = MIN( JU, J )
!                if (JM > J) call DGER( KM, JM-J, -1., AB( KV+2, J ), 1, AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), LDAB-1 )
!            end if
!            NW = MIN( 1, I3 )
!            if(NW > 0) call DCOPY( NW, AB( KV+ndL+1, J ), 1, WORK31( 1, 1 ), 1 )
!
!
!
!            if (J+1 <= N) then
!
!                J2 = MIN( JU-J+1, KV ) - 1
!                J3 = MAX( 0, JU-J-KV+1 )
!
!                call DLASWP( J2, AB( KV, J+1 ), LDAB-1, 1, 1, pivot( J ), 1 )
!
!                pivot( J ) = pivot( J ) + J - 1
!
!                K2 = J + J2
!                do I = 1, J3
!                    JJ = K2 + I
!                    do II = J + I - 1, J
!                        IP = pivot( II )
!                        if (IP /= II) then
!                            TEMP = AB( KV+1+II-JJ, JJ )
!                            AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
!                            AB( KV+1+IP-JJ, JJ ) = TEMP
!                        end if
!                    end do
!                end do
!
!                if (J2 > 0) then
!                    call DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J2, 1., AB( KV+1, J ), LDAB-1, AB( KV, J+1 ), LDAB-1 )
!                    if (I2 > 0) call DGEMM( 'No transpose', 'No transpose', I2, J2, 1, -1., AB( KV+2, J ), LDAB-1, AB( KV, J+1 ), LDAB-1, 1., AB( KV+1, J+1 ), LDAB-1 )
!                    if (I3 > 0) call DGEMM( 'No transpose', 'No transpose', I3, J2, 1, -1., WORK31, LDWORK, AB( KV, J+1 ), LDAB-1, 1., AB( KV+ndL, J+1 ), LDAB-1 )
!                end if
!
!                if (J3 > 0) then
!                    do JJ = 1, J3
!                        do II = JJ, 1
!                            WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
!                        end do
!                    end do
!
!                    call DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J3, 1., AB( KV+1, J ), LDAB-1, WORK13, LDWORK )
!                    if( I2>0 ) call DGEMM( 'No transpose', 'No transpose', I2, J3, 1, -1., AB( KV+2, J ), LDAB-1, WORK13, LDWORK, 1., AB( 2, J+KV ), LDAB-1 )
!                    if( I3>0 ) call DGEMM( 'No transpose', 'No transpose', I3, J3, 1, -1., WORK31, LDWORK, WORK13, LDWORK, 1., AB( 1+ndL, J+KV ), LDAB-1 )
!
!                    do JJ = 1, J3
!                        do II = JJ, 1
!                            AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
!                        end do
!                    end do
!                end if
!
!            else
!                pivot( J ) = pivot( J ) + J - 1
!            end if
!
!            JP = pivot( J ) - J + 1
!            if (JP /= 1) then
!                if (JP-1 < ndL) then
!                    call DSWAP( 0, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
!                else
!                    call DSWAP( 0, AB( KV+1, J ), LDAB-1, WORK31( JP-ndL, 1 ), LDWORK )
!                end if
!            end if
!            NW = MIN( I3, 1 )
!            if( NW>0 ) call DCOPY( NW, WORK31( 1, 1 ), 1, AB( KV+ndL+1, J ), 1 )
!
!        end do
!
!
!    end subroutine myDGBTRF
  
  
  
  
  
  
  
  
!    subroutine test_coeffs()
!
!        implicit none
!
!
!        ! ACHTUNG!!! Upwind-Differenzen noch nicht getestet!!!
!        !===========================================================================================================
!        if (iB(2,1) == 1 .and. iB(3,1) == 1) then
!            call test_diff(1,0,x1p              ,x1p              ,N1,b1L,b1U,cFp1 ,BC_1L,0,'cFp1' )
!            call test_diff(1,0,x1u              ,x1u              ,N1,b1L,b1U,cFu1 ,BC_1L,1,'cFu1' )
!            call test_diff(1,1,x1u              ,x1u              ,N1,b1L,b1U,cu1  ,BC_1L,1,'cu1'  )
!            call test_diff(1,1,x1p              ,x1p              ,N1,b1L,b1U,cp1  ,BC_1L,0,'cp1'  )
!            call test_diff(1,2,x1u              ,x1u              ,N1,b1L,b1U,cu11 ,BC_1L,1,'cu11' )
!            call test_diff(1,2,x1p              ,x1p              ,N1,b1L,b1U,cp11 ,BC_1L,0,'cp11' )
!
!            !CALL test_diff(1,1,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cGp1 ,BC_1L,2,'cGp1' ) ! TEST!!! Funktioniert nicht mehr ...
!            !CALL test_diff(1,1,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cDu1 ,BC_1L,3,'cDu1' )
!            !CALL test_diff(1,0,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cIpu ,BC_1L,2,'cIpu' )
!            !CALL test_diff(1,0,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cIup ,BC_1L,3,'cIup' )
!
!            call test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1D,BC_1L,1,'cNu1D')
!            call test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1U,BC_1L,1,'cNu1U')
!            call test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1D,BC_1L,0,'cNp1D')
!            call test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1U,BC_1L,0,'cNp1U')
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        if (iB(1,1) == 1 .and. iB(3,1) == 1) then
!            call test_diff(2,0,x2p              ,x2p              ,N2,b2L,b2U,cFp2 ,BC_2L,0,'cFp2' )
!            call test_diff(2,0,x2v              ,x2v              ,N2,b2L,b2U,cFv2 ,BC_2L,1,'cFv2' )
!            call test_diff(2,1,x2p              ,x2p              ,N2,b2L,b2U,cp2  ,BC_2L,0,'cp2'  )
!            call test_diff(2,1,x2v              ,x2v              ,N2,b2L,b2U,cv2  ,BC_2L,1,'cv2'  )
!            call test_diff(2,2,x2p              ,x2p              ,N2,b2L,b2U,cp22 ,BC_2L,0,'cp22' )
!            call test_diff(2,2,x2v              ,x2v              ,N2,b2L,b2U,cv22 ,BC_2L,1,'cv22' )
!
!            !CALL test_diff(2,1,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cGp2 ,BC_2L,2,'cGp2' )
!            !CALL test_diff(2,1,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cDv2 ,BC_2L,3,'cDv2' )
!            !CALL test_diff(2,0,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cIpv ,BC_2L,2,'cIpv' )
!            !CALL test_diff(2,0,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cIvp ,BC_2L,3,'cIvp' )
!
!            call test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2D,BC_2L,1,'cNv2D')
!            call test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2U,BC_2L,1,'cNv2U')
!            call test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2D,BC_2L,0,'cNp2D')
!            call test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2U,BC_2L,0,'cNp2U')
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. dimens == 3) then
!            call test_diff(3,0,x3p              ,x3p              ,N3,b3L,b3U,cFp3 ,BC_3L,0,'cFp3' )
!            call test_diff(3,0,x3w              ,x3w              ,N3,b3L,b3U,cFw3 ,BC_3L,1,'cFw3' )
!            call test_diff(3,1,x3p              ,x3p              ,N3,b3L,b3U,cp3  ,BC_3L,0,'cp3'  )
!            call test_diff(3,1,x3w              ,x3w              ,N3,b3L,b3U,cw3  ,BC_3L,1,'cw3'  )
!            call test_diff(3,2,x3p              ,x3p              ,N3,b3L,b3U,cp33 ,BC_3L,0,'cp33' )
!            call test_diff(3,2,x3w              ,x3w              ,N3,b3L,b3U,cw33 ,BC_3L,1,'cw33' )
!
!            !CALL test_diff(3,1,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cGp3 ,BC_3L,2,'cGp3' )
!            !CALL test_diff(3,1,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cDw3 ,BC_3L,3,'cDw3' )
!            !CALL test_diff(3,0,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cIpw ,BC_3L,2,'cIpw' )
!            !CALL test_diff(3,0,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cIwp ,BC_3L,3,'cIwp' )
!
!            call test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3D,BC_3L,1,'cNw3D')
!            call test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3U,BC_3L,1,'cNw3U')
!            call test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3D,BC_3L,0,'cNp3D')
!            call test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3U,BC_3L,0,'cNp3U')
!        end if
!    !===========================================================================================================
!
!
!    end subroutine test_coeffs
!
  
  
  
  
    !> \brief calculates finite differences coefficients.
    !!
    !! \param[in] rank only needed for error handling
    !! \param[in] Nmax dimension of cc and xc,xs
    !! \param[in] bl  stencil lower bound
    !! \param[in] bu  stencil upper bound
    !! \param[in] BCL local lower Boundary conditions
    !! \param[in] BCU local upper Boundary conditions
    !! \param[in] SShift shift in procgrid
    !! \param[in] grid_type (5: p, 1: u, 2: v, 3:w )
    !! \param[in] dir direction
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
    subroutine FD_getDiffCoeff( &
        rank,                   &
        Nmax,                   &
        bL, bU,                 &
        BCL, BCU,               &
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
  
        integer(c_int),  intent(in)   :: rank

        integer(c_int),  intent(in)   :: Nmax
        integer(c_int),  intent(in)   :: bL
        integer(c_int),  intent(in)   :: bU

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
  
        real(c_double),  intent(out)  :: cc(bL:bU,0:Nmax)
  
        integer               :: n_coeff
        integer               :: dim_n_coeff_bound
        integer               :: i, ii, iC, iStart
        integer               :: k, kk

        integer               :: left, right
  
        real                  :: dxi(1:2)
  
        real                  :: dxL, dxU ! Für Integrationskoeffizienten
  
        real   , allocatable  :: cc_xi(:,:)
        real   , allocatable  :: deltaX(:)
  
        logical               :: filter_yes

        integer               :: merror
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Upwinding wird auf dem Rand unterdrücken, da Differenzenstencils dort ohnehin schief sind.!
        !              - Generell koennte/sollte die Initialisierung der Index-Grenzen auch schon vor dieser       !
        !                Routine ausgefuehrt werden, um hier die Uebersichtlichkeit zu verbessern.                 !
        !----------------------------------------------------------------------------------------------------------!
  
  
        if (mapping_yes .and. abl > 2) then
            if (rank == 0) write(*,*) 'ERROR! Can`t handle derivatives > 2 combined with mapping ...'
            call MPI_FINALIZE(merror)
            stop
        end if
  
  
        !===========================================================================================================
        !=== Startindex für Entwicklungspunkte =====================================================================
        !===========================================================================================================
        if( grid_type == 5 .or. grid_type == 3 ) then
            iStart = 1
        else
            iStart = 0
        end if
        !===========================================================================================================
  
        filter_yes = .false.
        if (abl == 0 .and. (grid_type == 5 .or. grid_type == 1)) filter_yes = .true.
  
  
        dim_n_coeff_bound = SIZE(n_coeff_bound)
  
        cc = 0.
  
        do i = iStart, Nmax
     
            !========================================================================================================
            !=== Stencil-Breite auslesen ============================================================================
            !========================================================================================================
            if      (BCL > 0 .and. i <= (iStart - 1 + dim_n_coeff_bound)) then
                n_coeff = n_coeff_bound(i + 1 - iStart)
               !IF (upwind /= 0 .AND. i == iStart) n_coeff = 0
            else if (BCU > 0 .and. i >= (Nmax   + 1 - dim_n_coeff_bound)) then
                n_coeff = n_coeff_bound(Nmax + 1 - i)
               !IF (upwind /= 0 .AND. i == Nmax  ) n_coeff = 0
            else
                n_coeff = n_coeff_bound(dim_n_coeff_bound)
            end if
     
            !========================================================================================================
     
            if (n_coeff > 0 .and. n_coeff > abl) then
        
                !=====================================================================================================
                !=== Anzahl der Koeffizienten RECHTS vom Entwicklungspunkt bestimmen =================================
                !=====================================================================================================
                if (grid_type == 5 .or. grid_type == 1) then
                    !==================================================================================================
                    !=== "normales" Gitter, keine Zwischengitterpunkte ================================================
                    !==================================================================================================
                    if      (BCL > 0 .and. i <= (iStart - 1 + n_coeff/2)) then
                        right = n_coeff - i + iStart - 1
                    else if (BCU > 0 .and. i >= (Nmax   + 1 - n_coeff/2)) then
                        right = Nmax - i
                    else
                        ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= ungerade Anzahl Koeffizienten)!
                        if (MOD(n_coeff,2) == 0) then
                            if (rank == 0) then
                                write(*,'(a   )') 'ERROR! Choose odd number of coefficients!'
                                write(*,'(a,i4)') '    direction =', dir
                                write(*,'(a,i4)') '    grid_type =', grid_type
                                write(*,'(a,i4)') '            i =', i+SShift
                            end if
                        end if
                        right = (n_coeff-1)/2
                    end if
                else if (grid_type == 2) then
                    !==================================================================================================
                    !=== Druck ==> Impuls =============================================================================
                    !==================================================================================================
                    if      (BCL > 0 .and. i <= n_coeff/2) then
                        right = n_coeff - i - iStart
                    else if (BCU > 0 .and. i >= (Nmax + 0 - n_coeff/2)) then
                        right = Nmax - i
                    else
                        ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
                        if (MOD(n_coeff,2) /= 0) then
                            if (rank == 0) then
                                write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                                write(*,'(a,i4)') '    direction =', dir
                                write(*,'(a,i4)') '    grid_type =', grid_type
                                write(*,'(a,i4)') '            i =', i+SShift
                            end if
                            call MPI_FINALIZE(merror)
                            stop
                        end if
                        right = n_coeff/2
                    end if
                else if (grid_type == 3) then
                    !==================================================================================================
                    !=== Geschwindigkeit ==> Konti ====================================================================
                    !==================================================================================================
                    if      (BCL > 0 .and. i <= n_coeff/2) then
                        right = n_coeff - i - iStart
                    else if (BCU > 0 .and. i >= (Nmax + 1 - n_coeff/2)) then
                        right = Nmax - i
                    else
                        ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
                        if (MOD(n_coeff,2) /= 0) then
                            if (rank == 0) then
                                write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                                write(*,'(a,i4)') '    direction =', dir
                                write(*,'(a,i4)') '    grid_type =', grid_type
                                write(*,'(a,i4)') '            i =', i+SShift
                            end if
                            call MPI_FINALIZE(merror)
                            stop
                        end if
                        right = n_coeff/2 - 1
                    end if
                end if
                !=====================================================================================================
        
        
                !=====================================================================================================
                !=== Anzahl der Koeffizienten LINKS vom Entwicklungspunkt bestimmen ==================================
                !=====================================================================================================
                ! (0-ter bzw. 1-ter Koeffizient links vom Entwicklungspunkt
                !          == "zentraler" Koeffizient im Speicher)
                left = right - n_coeff + 1
                !=====================================================================================================
        
        
                !=====================================================================================================
                if (upwind == -1) then
                    if (.not. (BCL > 0 .and. i < (iStart + n_coeff/2))) then
                        n_coeff = n_coeff - 1
                        left    = left    + 1
                    end if
                end if
                if (upwind ==  1) then
                    if (.not. (BCU > 0 .and. i > (Nmax   - n_coeff/2))) then
                        n_coeff = n_coeff - 1
                        right   = right   - 1
                    end if
                end if
                !=====================================================================================================
        
        
                !=====================================================================================================
                !=== Stencilanordnung testen =========================================================================
                !=====================================================================================================
                if (right > bU .or. left < bL) then
                    if (rank == 0) then
                        !WRITE(*,'(a   )') 'WARNING! The FD-Stencil does probably not fit into provided array!'
                        write(*,'(a   )') 'ERROR! Stencil doesn`t fit into provided array!'
                        write(*,'(a,i4)') '    direction =', dir
                        write(*,'(a,i4)') '    grid_type =', grid_type
                        write(*,'(a,i4)') '            i =', i+SShift
                    end if
                    call MPI_FINALIZE(merror)
                    stop
                end if
                !=====================================================================================================
        
        
                allocate(deltaX(1:n_coeff))
        
        
                !=====================================================================================================
                !=== räumliche Abstände zum Entwicklungspunkt ========================================================
                !=====================================================================================================
                do ii = 1, n_coeff
                    ! Koeffizienten-Punkte ("iC"):
                    iC = i + left + (ii-1)
           
                    if (mapping_yes) then
                        if      (grid_type == 2) then
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
                !=====================================================================================================
        
        
        
                !=====================================================================================================
                !=== Integrations-Intervall ==========================================================================
                !=====================================================================================================
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
                !=====================================================================================================
        
        
                !=====================================================================================================
                !=== Bestimmung der Koeffizienten ====================================================================
                !=====================================================================================================
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
                            call FD_coeffs_solver(rank,k,filter_yes,n_coeff,deltaX,cc_xi(left:right,k))
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
           
                else if (mapping_yes .and. abl == 0 .and. .not. filter_yes) then ! TEST!!!
           
                    !--- interpolation (mapped) ---
           
                    !k  = 1
                    !kk = 0
                    !ALLOCATE(cc_xi(left:right,k:k))
                    !INCLUDE 'FD_coeffs_expl_map.f90'
                    !cc(left:right,i) = cc_xi(left:right,k)
                    !DEALLOCATE(cc_xi)
           
                    if (n_coeff <= 7) then ! TEST!!!
                        call diff_coeffs_exact(abl,n_coeff,deltaX(1),cc(left:right,i))
                    else
                        call FD_coeffs_solver(rank,abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
                    end if
           
                else
           
                    if (abl == -1) then
                        !--- integration coefficients ---
                        call FD_coeffs_solver_integral(rank,n_coeff,deltaX,dxL,dxU,cc(left:right,i))
                    else
                        !--- interpolation and derivatives ---
                        call FD_coeffs_solver(rank,abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
                    end if
           
                end if
                !=====================================================================================================
        
        
                deallocate(deltaX)
        
        
                !=====================================================================================================
                !=== Symmetrie =======================================================================================
                !=====================================================================================================
                if (BCL == -2 .and. abl == -1 .and. i == iStart) then
                    cc(0,i) = 0.5*cc(0,i)
                    do ii = bL, -1
                        cc(ii,i) = 0.
                    end do
                end if
                if (BCU == -2 .and. abl == -1 .and. i == Nmax  ) then
                    cc(0,i) = 0.5*cc(0,i)
                    do ii = 1, bU
                        cc(ii,i) = 0.
                    end do
                end if
                !-----------------------------------------------------------------------------------------------------
                !-----------------------------------------------------------------------------------------------------
                ! TEST!!! abl == -1 ???
                if (BCL == -2 .and. (i+bL) < 1) then
                    if (grid_type == 5 .or. (grid_type == 2 .and. i >= 1)) then
                        do ii = 1, 1-(i+bL)
                            cc(1-i+ii,i) = cc(1-i+ii,i) + cc(1-i-ii,i)
                            cc(1-i-ii,i) = 0.
                        end do
                    end if
                    if ((grid_type == 1 .and. i >= 1) .or. grid_type == 3) then
                        do ii = 1, 1-(i+bL)
                            cc(0-i+ii,i) = cc(0-i+ii,i) - cc(1-i-ii,i)
                            cc(1-i-ii,i) = 0.
                        end do
                    end if
                end if
                !-----------------------------------------------------------------------------------------------------
                if (BCU == -2 .and. (i+bU) > (Nmax+0)) then
                    if (grid_type == 5 .or. (grid_type == 2 .and. i <= Nmax-1)) then
                        do ii = 1, i+bU-(Nmax-0)
                            cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) + cc(Nmax-i  +ii,i)
                            cc(Nmax-i  +ii,i) = 0.
                        end do
                    end if
                end if
                if (BCU == -2 .and. (i+bU) > (Nmax-1)) then
                    if ((grid_type == 1 .and. i <= Nmax-1) .or. grid_type == 3) then
                        do ii = 1, i+bU-(Nmax-1)
                            cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) - cc(Nmax-i-1+ii,i)
                            cc(Nmax-i-1+ii,i) = 0.
                        end do
                    end if
                end if
               !=====================================================================================================
        
            end if
     
        end do
  
  
    end subroutine FD_getDiffCoeff
  
  
  
  
  
    subroutine diff_coeffs_exact(abl,n_coeff,deltaX,cc) ! TEST!!!

        implicit none

        integer, intent(in   ) ::  abl
        integer, intent(in   ) ::  n_coeff
        real   , intent(in   ) ::  deltaX
        real   , intent(out  ) ::  cc(1:n_coeff)

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
  
  
  
  
  
  
  
  
  
  
  
    subroutine FD_coeffs_solver(rank,abl,filter_yes,n_coeff,deltaX,cc)

        implicit none

        integer, intent(in)   ::  rank
        integer, intent(in)   ::  abl
        logical, intent(in)   ::  filter_yes
        integer, intent(in)   ::  n_coeff
        real   , intent(in)   ::  deltaX(1:n_coeff)
        real   , intent(out)  ::  cc    (1:n_coeff)

        integer               ::  i, j

        real                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
        real                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)

        real                  ::  const


        !===========================================================================================================
        !=== Aufstellen des Gleichungssystems ======================================================================
        !===========================================================================================================
        do i = 1, n_coeff
            do j = 1, n_coeff
                polyn_vals(i,j) = deltaX(i)**(j-1)
            end do
        end do


        !===========================================================================================================
        !=== Zusatzbedingungen =====================================================================================
        !===========================================================================================================
        if (filter_yes) then
            ! G(pi) = 0.
            do i = 1, n_coeff
                polyn_vals(i,n_coeff) = (-1.)**i
            end do
        end if


        !===========================================================================================================
        !=== Lösen des Gleichungssystems ===========================================================================
        !===========================================================================================================
        call Matrix_invert(rank,n_coeff,polyn_vals,polyn_vals_inv)


        !===========================================================================================================
        !=== Funktionswert am Entwicklungspunkt (deltaX = 0) =======================================================
        !===========================================================================================================
        const = 1.
        do i = 1, abl-1
            const = const*REAL(1+i)
        end do


        !===========================================================================================================
        !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
        !===========================================================================================================
        cc(1:n_coeff) = const*polyn_vals_inv(1+abl,1:n_coeff)


    end subroutine FD_coeffs_solver
  
  
  
  
  
  
  
  
  
  
  
    subroutine FD_coeffs_solver_integral(rank,n_coeff,deltaX,dxL,dxU,cc)

        implicit none

        integer, intent(in)   ::  rank

        integer, intent(in)   ::  n_coeff
        real   , intent(in)   ::  deltaX(1:n_coeff)
        real   , intent(out)  ::  cc    (1:n_coeff)
        real   , intent(in)   ::  dxL, dxU

        integer               ::  i, j, k

        real                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
        real                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)

        real                  ::  const


        !===========================================================================================================
        !=== Aufstellen des Gleichungssystems ======================================================================
        !===========================================================================================================
        do i = 1, n_coeff
            do j = 1, n_coeff
                polyn_vals(i,j) = deltaX(i)**(j-1)
            end do
        end do


        !===========================================================================================================
        !=== Lösen des Gleichungssystems ===========================================================================
        !===========================================================================================================
        call Matrix_invert(rank,n_coeff,polyn_vals,polyn_vals_inv)


        !===========================================================================================================
        !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
        !===========================================================================================================
        !  (Matrix-Vektor-Multiplikation)
        cc = 0.
        do j = 1, n_coeff
            do i = 1, n_coeff
                cc(j) = cc(j) + (dxU**i - dxL**i)/REAL(i)*polyn_vals_inv(i,j)
            end do
        end do


    end subroutine FD_coeffs_solver_integral
  
  
  
  
  
  
  
  
  
  
  
!    subroutine get_stencil()
!
!        implicit none
!
!        integer                ::  i, imax
!        integer                ::  j, jmax
!        integer                ::  k, kmax
!
!        integer                ::  g
!
!        real                   ::  cDu1R(-1:0,1:N1)
!        real                   ::  cDv2R(-1:0,1:N2)
!        real                   ::  cDw3R(-1:0,1:N3)
!
!        real                   ::  cGp1R( 0:1,0:N1)
!        real                   ::  cGp2R( 0:1,0:N2)
!        real                   ::  cGp3R( 0:1,0:N3)
!
!
!        cdg1 = 0.
!        cdg2 = 0.
!        cdg3 = 0.
!
!        do g = 1, n_grids
!
!            !========================================================================================================
!            imax = NN(1,g)
!
!            cDu1R = 0.
!            cGp1R = 0.
!
!            !---------------------------------------------------------------------------!
!            ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
!            !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
!            !---------------------------------------------------------------------------!
!
!            do i = 1, imax
!                cDu1R(-1,i) = -1./(x1uR(i  ,g) - x1uR(i-1,g))
!                cDu1R( 0,i) =  1./(x1uR(i  ,g) - x1uR(i-1,g))
!            end do
!
!            do i = 1, imax-1
!                cGp1R( 0,i) = -1./(x1pR(i+1,g) - x1pR(i  ,g))
!                cGp1R( 1,i) =  1./(x1pR(i+1,g) - x1pR(i  ,g))
!            end do
!
!            if (BC(1,1,g) > 0) then
!                cGp1R( :,0   ) =  0.
!                cDu1R( :,1   ) =  0.
!                !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
!                cDu1R( 0,1   ) =  1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
!            else
!                cGp1R( 0,0   ) = -1./(x1pR(1,g) - x1pR(0,g))
!                cGp1R( 1,0   ) =  1./(x1pR(1,g) - x1pR(0,g))
!            end if
!
!            if (BC(2,1,g) > 0) then
!                cGp1R( :,imax) =  0.
!                cDu1R( :,imax) =  0.
!                !cDu1R(-1,imax) = -1./(x1u (N1      ) - x1u (N1-1  )) ! TEST!!! Das geht in die Hose ...
!                cDu1R(-1,imax) = -1./(y1u (M1      ) - y1u (M1-1  ))
!            else
!                cGp1R( 0,imax) = -1./(x1pR(imax+1,g) - x1pR(imax,g))
!                cGp1R( 1,imax) =  1./(x1pR(imax+1,g) - x1pR(imax,g))
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            do i = 1, imax
!                cdg1(-1,i,g) = cDu1R(-1,i)*cGp1R(0,i-1)
!                cdg1( 0,i,g) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
!                cdg1( 1,i,g) =                            cDu1R(0,i)*cGp1R(1,i)
!            end do
!
!            ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
!            if (BC(1,1,g) > 0) cdg1( :,1   ,g) = 2.*cdg1(:,1   ,g) ! TEST!!! Ist das wirklich so optimal?
!            if (BC(2,1,g) > 0) cdg1( :,imax,g) = 2.*cdg1(:,imax,g)
!
!            if (BC(1,1,g) == -2) then
!                cdg1( 1,1   ,g) = cdg1( 1,1   ,g) + cdg1(-1,1   ,g)
!                cdg1(-1,1   ,g) = 0.
!            end if
!            if (BC(2,1,g) == -2) then
!                cdg1(-1,imax,g) = cdg1(-1,imax,g) + cdg1( 1,imax,g)
!                cdg1( 1,imax,g) = 0.
!            end if
!            !========================================================================================================
!            jmax = NN(2,g)
!
!            cDv2R = 0.
!            cGp2R = 0.
!
!            do j = 1, jmax
!                cDv2R(-1,j) = -1./(x2vR(j  ,g) - x2vR(j-1,g))
!                cDv2R( 0,j) =  1./(x2vR(j  ,g) - x2vR(j-1,g))
!            end do
!
!            do j = 1, jmax-1
!                cGp2R( 0,j) = -1./(x2pR(j+1,g) - x2pR(j  ,g))
!                cGp2R( 1,j) =  1./(x2pR(j+1,g) - x2pR(j  ,g))
!            end do
!
!            if (BC(1,2,g) > 0) then
!                cGp2R( :,0   ) =  0.
!                cDv2R( :,1   ) =  0.
!                !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
!                cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
!            else
!                cGp2R( 0,0   ) = -1./(x2pR(1,g) - x2pR(0  ,g))
!                cGp2R( 1,0   ) =  1./(x2pR(1,g) - x2pR(0  ,g))
!            end if
!
!            if (BC(2,2,g) > 0) then
!                cGp2R( :,jmax) =  0.
!                cDv2R( :,jmax) =  0.
!                !cDv2R(-1,jmax) = -1./(x2v (N2      ) - x2v (N2-1  )) ! TEST!!! Das geht in die Hose ...
!                cDv2R(-1,jmax) = -1./(y2v (M2      ) - y2v (M2-1  ))
!            else
!                cGp2R( 0,jmax) = -1./(x2pR(jmax+1,g) - x2pR(jmax,g))
!                cGp2R( 1,jmax) =  1./(x2pR(jmax+1,g) - x2pR(jmax,g))
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            do j = 1, jmax
!                cdg2(-1,j,g) = cDv2R(-1,j)*cGp2R(0,j-1)
!                cdg2( 0,j,g) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
!                cdg2( 1,j,g) =                            cDv2R(0,j)*cGp2R(1,j)
!            end do
!
!            ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
!            if (BC(1,2,g) > 0) cdg2( :,1   ,g) = 2.*cdg2(:,1   ,g)
!            if (BC(2,2,g) > 0) cdg2( :,jmax,g) = 2.*cdg2(:,jmax,g)
!
!            if (BC(1,2,g) == -2) then
!                cdg2( 1,1   ,g) = cdg2( 1,1   ,g) + cdg2(-1,1   ,g)
!                cdg2(-1,1   ,g) = 0.
!            end if
!            if (BC(2,2,g) == -2) then
!                cdg2(-1,jmax,g) = cdg2(-1,jmax,g) + cdg2( 1,jmax,g)
!                cdg2( 1,jmax,g) = 0.
!            end if
!            !========================================================================================================
!            if (dimens == 3) then
!
!                kmax = NN(3,g)
!
!                cDw3R = 0.
!                cGp3R = 0.
!
!                do k = 1, kmax
!                    cDw3R(-1,k) = -1./(x3wR(k  ,g) - x3wR(k-1,g))
!                    cDw3R( 0,k) =  1./(x3wR(k  ,g) - x3wR(k-1,g))
!                end do
!
!                do k = 1, kmax-1
!                    cGp3R( 0,k) = -1./(x3pR(k+1,g) - x3pR(k  ,g))
!                    cGp3R( 1,k) =  1./(x3pR(k+1,g) - x3pR(k  ,g))
!                end do
!
!                if (BC(1,3,g) > 0) then
!                    cGp3R( :,0   ) =  0.
!                    cDw3R( :,1   ) =  0.
!                    !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
!                    cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
!                else
!                    cGp3R( 0,0   ) = -1./(x3pR(1,g) - x3pR(0,g))
!                    cGp3R( 1,0   ) =  1./(x3pR(1,g) - x3pR(0,g))
!                end if
!
!                if (BC(2,3,g) > 0) then
!                    cGp3R( :,kmax) =  0.
!                    cDw3R( :,kmax) =  0.
!                    !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
!                    cDw3R(-1,kmax) = -1./(y3w (M3      ) - y3w (M3-1  ))
!                else
!                    cGp3R( 0,kmax) = -1./(x3pR(kmax+1,g) - x3pR(kmax,g))
!                    cGp3R( 1,kmax) =  1./(x3pR(kmax+1,g) - x3pR(kmax,g))
!                end if
!                !--------------------------------------------------------------------------------------------------------
!                do k = 1, kmax
!                    cdg3(-1,k,g) = cDw3R(-1,k)*cGp3R(0,k-1)
!                    cdg3( 0,k,g) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
!                    cdg3( 1,k,g) =                            cDw3R(0,k)*cGp3R(1,k)
!                end do
!
!                ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
!                if (BC(1,3,g) > 0) cdg3( :,1   ,g) = 2.*cdg3(:,1   ,g)
!                if (BC(2,3,g) > 0) cdg3( :,kmax,g) = 2.*cdg3(:,kmax,g)
!
!                if (BC(1,3,g) == -2) then
!                    cdg3( 1,1   ,g) = cdg3( 1,1   ,g) + cdg3(-1,1   ,g)
!                    cdg3(-1,1   ,g) = 0.
!                end if
!                if (BC(2,3,g) == -2) then
!                    cdg3(-1,kmax,g) = cdg3(-1,kmax,g) + cdg3( 1,kmax,g)
!                    cdg3( 1,kmax,g) = 0.
!                end if
!
!            end if
!           !========================================================================================================
!
!        end do
!
!
!    end subroutine get_stencil
  
  
  
  
  
  
  
  
  
  
  
!    subroutine get_stencil_transp
!
!        implicit none
!
!        integer                ::  i, j, k, g
!        real   , allocatable   ::  work(:,:)
!
!
!        do g = 1, n_grids
!
!            !========================================================================================================
!            allocate(work(-1:1,0:(NN(1,g)+1)))
!
!            work = 0.
!            do i = 1, NN(1,g)
!                work(-1:1,i) = cdg1(-1:1,i,g)
!            end do
!
!            cdg1(:,:,g) = 0.
!            do i = 1, NN(1,g)
!                cdg1(-1,i,g) = work( 1,i-1)
!                cdg1( 0,i,g) = work( 0,i  )
!                cdg1( 1,i,g) = work(-1,i+1)
!            end do
!
!            deallocate(work)
!
!            if (BC(1,1,g) == 0 .or. BC(1,1,g) == -1) then
!                i = 1
!                cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
!                cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
!            end if
!            if (BC(2,1,g) == 0 .or. BC(2,1,g) == -1) then
!                i = NN(1,g)
!                cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
!                cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
!            end if
!            !========================================================================================================
!            allocate(work(-1:1,0:(NN(2,g)+1)))
!
!            work = 0.
!            do j = 1, NN(2,g)
!                work(-1:1,j) = cdg2(-1:1,j,g)
!            end do
!
!            cdg2(:,:,g) = 0.
!            do j = 1, NN(2,g)
!                cdg2(-1,j,g) = work( 1,j-1)
!                cdg2( 0,j,g) = work( 0,j  )
!                cdg2( 1,j,g) = work(-1,j+1)
!            end do
!
!            deallocate(work)
!
!            if (BC(1,2,g) == 0 .or. BC(1,2,g) == -1) then
!                j = 1
!                cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
!                cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
!            end if
!            if (BC(2,2,g) == 0 .or. BC(2,2,g) == -1) then
!                j = NN(2,g)
!                cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
!                cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
!            end if
!            !========================================================================================================
!            if (dimens == 3) then
!
!                allocate(work(-1:1,0:(NN(3,g)+1)))
!
!                work = 0.
!                do k = 1, NN(3,g)
!                    work(-1:1,k) = cdg3(-1:1,k,g)
!                end do
!
!                cdg3(:,:,g) = 0.
!                do k = 1, NN(3,g)
!                    cdg3(-1,k,g) = work( 1,k-1)
!                    cdg3( 0,k,g) = work( 0,k  )
!                    cdg3( 1,k,g) = work(-1,k+1)
!                end do
!
!                deallocate(work)
!
!                if (BC(1,3,g) == 0 .or. BC(1,3,g) == -1) then
!                    k = 1
!                    cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
!                    cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
!                end if
!                if (BC(2,3,g) == 0 .or. BC(2,3,g) == -1) then
!                    k = NN(3,g)
!                    cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
!                    cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
!                end if
!
!            end if
!           !========================================================================================================
!
!        end do
!
!
!        ! sicherheitshalber (cdg3 = 0. sollte eigentlich schon erf�llt sein):
!        if (dimens == 2) cdg3 = 0.
!
!
!    end subroutine get_stencil_transp
  
  
  
  
  
  
  
  
  
  
  
!    subroutine get_stencil_Helm()
!
!        implicit none
!
!        integer                ::  imax, jmax, kmax
!        integer                ::  g, m
!
!
!        !----------------------------------------------------------------------------------------------------------!
!        ! Anmerkungen: - Bei Symmetrie ergeben sich automatisch auf der Symmetrie-Ebene Dirichlet-RB für die       !
!        !                Normalkomponente der Geschwindigkeit (betrifft cu11R, cv22R, cw33R und tritt aufgrund des !
!        !                Gittertyps nur auf den gröberen Gittern auf). Die zugehörigen Restriktionskoeffizienten   !
!        !                liefern automatisch entsprechend die "Randbedingung" vel_n = 0.                           !
!        !              - Bei Symmetrie-RB werden in "CALL diff_coeffs" die Koeffizienten am Rand gemäss Gittertyp  !
!        !                manipuliert, d.h. bei cp11R, cp22R, cp33R werden Geschwindigkeitskomponenten parallel     !
!        !                zur Symmetrieebene (bzw. Druck) angenommen, so dass cu11R, cv22R, cw33R nach dem Kopieren !
!        !                von cp11R, cp22R, cp33R entsprechend für Geschwindigkeitskomponenten normal zur           !
!        !                Symmetrie-Ebene geändert werden müssen.                                                   !
!        !----------------------------------------------------------------------------------------------------------!
!
!        ! sicherheitshalber (dimens == 2, ...):
!        cp11R = 0.
!        cp22R = 0.
!        cp33R = 0.
!
!        cu11R = 0.
!        cv22R = 0.
!        cw33R = 0.
!
!        cc11R = 0.
!        cc22R = 0.
!        cc33R = 0.
!
!
!        !===========================================================================================================
!        !=== grobe & feinstes Gitter ===============================================================================
!        !===========================================================================================================
!        do g = 1, n_grids
!
!            !========================================================================================================
!            imax = NN(1,g)
!
!            call diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1pR(-1:imax+1,g),x1pR(-1:imax+1,g),BC_1L,BC_1U,0,imax,-1,1,cp11R(-1,0,g))
!
!            if (BC_1L == 1 .or. BC_1L == 3) then
!                cp11R( :,1   ,g) =  0.
!                cp11R( 0,1   ,g) =  1.
!            end if
!            if (BC_1U == 1 .or. BC_1U == 3) then
!                cp11R( :,imax,g) =  0.
!                cp11R( 0,imax,g) =  1.
!            end if
!
!            if (BC_1L == 2 .or. BC_1L == 4) then
!                cp11R(-1,1   ,g) =  0.
!                cp11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
!                cp11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
!            end if
!            if (BC_1U == 2 .or. BC_1U == 4) then
!                cp11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
!                cp11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
!                cp11R( 1,imax,g) =  0.
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            cu11R(:,:,g) = cp11R(:,:,g)
!
!            if (BC_1L ==  2) then
!                cu11R( :,1   ,g) =  0.
!                cu11R( 0,1   ,g) =  1.
!            end if
!            if (BC_1U ==  2) then
!                cu11R( :,imax,g) =  0.
!                cu11R( 0,imax,g) =  1.
!            end if
!
!            if (BC_1L ==  3) then
!                cu11R(-1,1   ,g) =  0.
!                cu11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
!                cu11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
!            end if
!            if (BC_1U ==  3) then
!                cu11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
!                cu11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
!                cu11R( 1,imax,g) =  0.
!            end if
!
!            if (BC_1L == -2) then
!                cu11R(-1,1   ,g) =  0.
!                cu11R( 1,1   ,g) =  0.
!            end if
!            if (BC_1U == -2) then
!                cu11R(-1,imax,g) =  0.
!                cu11R( 1,imax,g) =  0.
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            do m = 1, n_conc
!                cc11R(:,:,g,m) = cp11R(:,:,g)
!
!                if (BCc_1L(m) > 0) then
!                    cc11R(-1,1   ,g,m) =  0.
!                    cc11R( 0,1   ,g,m) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
!                    cc11R( 1,1   ,g,m) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
!                end if
!                if (BCc_1U(m) > 0) then
!                    cc11R(-1,imax,g,m) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
!                    cc11R( 0,imax,g,m) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
!                    cc11R( 1,imax,g,m) =  0.
!                end if
!            end do
!            !--------------------------------------------------------------------------------------------------------
!            !--------------------------------------------------------------------------------------------------------
!            jmax = NN(2,g)
!
!            call diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2pR(-1:jmax+1,g),x2pR(-1:jmax+1,g),BC_2L,BC_2U,0,jmax,-1,1,cp22R(-1,0,g))
!
!            if (BC_2L == 1 .or. BC_2L == 3) then
!                cp22R(:,1    ,g) =  0.
!                cp22R(0,1    ,g) =  1.
!            end if
!            if (BC_2U == 1 .or. BC_2U == 3) then
!                cp22R( :,jmax,g) =  0.
!                cp22R( 0,jmax,g) =  1.
!            end if
!
!            if (BC_2L == 2 .or. BC_2L == 4) then
!                cp22R(-1,1   ,g) =  0.
!                cp22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
!                cp22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
!            end if
!            if (BC_2U == 2 .or. BC_2U == 4) then
!                cp22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                cp22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                cp22R( 1,jmax,g) =  0.
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            cv22R(:,:,g) = cp22R(:,:,g)
!
!            if (BC_2L ==  2) then
!                cv22R( :,1   ,g) =  0.
!                cv22R( 0,1   ,g) =  1.
!            end if
!            if (BC_2U ==  2) then
!                cv22R( :,jmax,g) =  0.
!                cv22R( 0,jmax,g) =  1.
!            end if
!
!            if (BC_2L ==  3) then
!                cv22R(-1,1   ,g) =  0.
!                cv22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
!                cv22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
!            end if
!            if (BC_2U ==  3) then
!                cv22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                cv22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                cv22R( 1,jmax,g) =  0.
!            end if
!
!            if (BC_2L == -2) then
!                cv22R(-1,1   ,g) =  0.
!                cv22R( 1,1   ,g) =  0.
!            end if
!            if (BC_2U == -2) then
!                cv22R(-1,jmax,g) =  0.
!                cv22R( 1,jmax,g) =  0.
!            end if
!            !--------------------------------------------------------------------------------------------------------
!            do m = 1, n_conc
!                cc22R(:,:,g,m) = cp22R(:,:,g)
!
!                if (BCc_2L(m) > 0) then
!                    cc22R(-1,1   ,g,m) =  0.
!                    cc22R( 0,1   ,g,m) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
!                    cc22R( 1,1   ,g,m) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
!                end if
!                if (BCc_2U(m) > 0) then
!                    cc22R(-1,jmax,g,m) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                    cc22R( 0,jmax,g,m) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
!                    cc22R( 1,jmax,g,m) =  0.
!                end if
!            end do
!            !--------------------------------------------------------------------------------------------------------
!            !--------------------------------------------------------------------------------------------------------
!            if (dimens == 3) then
!
!                kmax = NN(3,g)
!
!                call diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3pR(-1:kmax+1,g),x3pR(-1:kmax+1,g),BC_3L,BC_3U,0,kmax,-1,1,cp33R(-1,0,g))
!
!                if (BC_3L == 1 .or. BC_3L == 3) then
!                    cp33R( :,1   ,g) =  0.
!                    cp33R( 0,1   ,g) =  1.
!                end if
!                if (BC_3U == 1 .or. BC_3U == 3) then
!                    cp33R( :,kmax,g) =  0.
!                    cp33R( 0,kmax,g) =  1.
!                end if
!
!                if (BC_3L == 2 .or. BC_3L == 4) then
!                    cp33R(-1,1   ,g) =  0.
!                    cp33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
!                    cp33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
!                end if
!                if (BC_3U == 2 .or. BC_3U == 4) then
!                    cp33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                    cp33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                    cp33R( 1,kmax,g) =  0.
!                end if
!                !--------------------------------------------------------------------------------------------------------
!                cw33R(:,:,g) = cp33R(:,:,g)
!
!                if (BC_3L ==  2) then
!                    cw33R( :,1   ,g) =  0.
!                    cw33R( 0,1   ,g) =  1.
!                end if
!                if (BC_3U ==  2) then
!                    cw33R( :,kmax,g) =  0.
!                    cw33R( 0,kmax,g) =  1.
!                end if
!
!                if (BC_3L ==  3) then
!                    cw33R(-1,1   ,g) =  0.
!                    cw33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
!                    cw33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
!                end if
!                if (BC_3U ==  3) then
!                    cw33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                    cw33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                    cw33R( 1,kmax,g) =  0.
!                end if
!
!                if (BC_3L == -2) then
!                    cw33R(-1,1   ,g) =  0.
!                    cw33R( 1,1   ,g) =  0.
!                end if
!                if (BC_3U == -2) then
!                    cw33R(-1,kmax,g) =  0.
!                    cw33R( 1,kmax,g) =  0.
!                end if
!                !--------------------------------------------------------------------------------------------------------
!                do m = 1, n_conc
!                    cc33R(:,:,g,m) = cp33R(:,:,g)
!
!                    if (BCc_3L(m) > 0) then
!                        cc33R(-1,1   ,g,m) =  0.
!                        cc33R( 0,1   ,g,m) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
!                        cc33R( 1,1   ,g,m) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
!                    end if
!                    if (BCc_3U(m) > 0) then
!                        cc33R(-1,kmax,g,m) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                        cc33R( 0,kmax,g,m) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
!                        cc33R( 1,kmax,g,m) =  0.
!                    end if
!                end do
!
!            end if
!           !========================================================================================================
!
!        end do
!
!        g = 1
!
!        !===========================================================================================================
!        !=== feinstes Gitter (Tangential-Koeffizienten überschreiben) ==============================================
!        !===========================================================================================================
!        call diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1u(-1:N1+1),x1u(-1:N1+1),BC_1L,BC_1U,1,N1,-1,1,cu11R(-1,0,g))
!
!        if (BC_1L == 1 .or. BC_1L == 2) then
!            cu11R(-1,0 ,g) = 0.
!            cu11R( 0,0 ,g) = 1.- (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
!            cu11R( 1,0 ,g) =     (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
!        end if
!        if (BC_1U == 1 .or. BC_1U == 2) then
!            cu11R(-1,N1,g) = 1.- (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
!            cu11R( 0,N1,g) =     (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
!            cu11R( 1,N1,g) = 0.
!        end if
!
!        if (BC_1L == 3 .or. BC_1L == 4) then
!            cu11R(-1,0 ,g) =  0.
!            cu11R( 0,0 ,g) = -1./(x1u(1 )-x1u(0   ))
!            cu11R( 1,0 ,g) =  1./(x1u(1 )-x1u(0   ))
!        end if
!        if (BC_1U == 3 .or. BC_1U == 4) then
!            cu11R(-1,N1,g) = -1./(x1u(N1)-x1u(N1-1))
!            cu11R( 0,N1,g) =  1./(x1u(N1)-x1u(N1-1))
!            cu11R( 1,N1,g) =  0.
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        call diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2v(-1:N2+1),x2v(-1:N2+1),BC_2L,BC_2U,1,N2,-1,1,cv22R(-1,0,g))
!
!        if (BC_2L == 1 .or. BC_2L == 2) then
!            cv22R(-1,0 ,g) = 0.
!            cv22R( 0,0 ,g) = 1.- (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
!            cv22R( 1,0 ,g) =     (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
!        end if
!        if (BC_2U == 1 .or. BC_2U == 2) then
!            cv22R(-1,N2,g) = 1.- (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
!            cv22R( 0,N2,g) =     (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
!            cv22R( 1,N2,g) = 0.
!        end if
!
!        if (BC_2L == 3 .or. BC_2L == 4) then
!            cv22R(-1,0 ,g) =  0.
!            cv22R( 0,0 ,g) = -1./(x2v(1 )-x2v(0   ))
!            cv22R( 1,0 ,g) =  1./(x2v(1 )-x2v(0   ))
!        end if
!        if (BC_2U == 3 .or. BC_2U == 4) then
!            cv22R(-1,N2,g) = -1./(x2v(N2)-x2v(N2-1))
!            cv22R( 0,N2,g) =  1./(x2v(N2)-x2v(N2-1))
!            cv22R( 1,N2,g) =  0.
!        end if
!        !-----------------------------------------------------------------------------------------------------------
!        if (dimens == 3) then
!
!            call diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3w(-1:N3+1),x3w(-1:N3+1),BC_3L,BC_3U,1,N3,-1,1,cw33R(-1,0,g))
!
!            if (BC_3L == 1 .or. BC_3L == 2) then
!                cw33R(-1,0 ,g) = 0.
!                cw33R( 0,0 ,g) = 1.- (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
!                cw33R( 1,0 ,g) =     (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
!            end if
!            if (BC_3U == 1 .or. BC_3U == 2) then
!                cw33R(-1,N3,g) = 1.- (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
!                cw33R( 0,N3,g) =     (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
!                cw33R( 1,N3,g) = 0.
!            end if
!
!            if (BC_3L == 3 .or. BC_3L == 4) then
!                cw33R(-1,0 ,g) =  0.
!                cw33R( 0,0 ,g) = -1./(x3w(1 )-x3w(0   ))
!                cw33R( 1,0 ,g) =  1./(x3w(1 )-x3w(0   ))
!            end if
!            if (BC_3U == 3 .or. BC_3U == 4) then
!                cw33R(-1,N3,g) = -1./(x3w(N3)-x3w(N3-1))
!                cw33R( 0,N3,g) =  1./(x3w(N3)-x3w(N3-1))
!                cw33R( 1,N3,g) =  0.
!            end if
!
!        end if
!    !===========================================================================================================
!
!
!    end subroutine get_stencil_Helm
!


    subroutine Matrix_invert(rank,N,matrix,matrix_inv)

        implicit none

        integer, intent(in)  ::  rank
        integer, intent(in)  ::  N
        real   , intent(in)  ::  matrix     (1:N,1:N)
        real   , intent(out) ::  matrix_inv (1:N,1:N)

        real                 ::  matrix_left(1:N,1:N)
        real                 ::  mult1, mult2
        real                 ::  eps
        integer              ::  i, j, k

        real                 ::  store
        integer              ::  maxValue, maxValuePos



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
            maxValue    = ABS(matrix_left(j,j))
            maxValuePos = j
            do i = j+1, N
                if (ABS(matrix_left(i,j)) > maxValue) then
                    maxValue    = ABS(matrix_left(i,j))
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
            if (ABS(mult1) < eps .and. rank == 0) then
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
            if (ABS(mult1) < eps .and. rank == 0) then
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
            if (ABS(mult1) < eps .and. rank == 0) then
                write(*,'(a)'      ) 'WARNING! Matrix is probably singular (3)!'
                write(*,'(a,E13.5)') '       ... dividing by', mult1
                write(*,'(a,i2)'   ) '                   i =', i
            end if
            matrix_inv(i,:) = matrix_inv(i,:)/mult1
        end do


    end subroutine Matrix_invert

  
end module cmod_FDcoeff
