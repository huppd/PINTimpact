!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module for computing laplace
!! should be moved to cmod_Operator
module cmod_laplace


    use mod_dims
    use mod_vars
    use mod_diff
    use mod_exchange
  
    use iso_c_binding
  

    private
  
    public product_div_grad, product_div_grad_transp
    public product_div_grad_relax, relaxation_div_grad, relaxation_div_grad_inv
    public handle_corner_Lap
  
  
contains
  
    !pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
  
  
  
    !> \brief computes first \c lap = \c div( \c grad( \c phi))
    !! computes realy first \c grad then \c div
    !! dirty hack
    subroutine product_div_grad( ccorner_yes, phi, Lap ) bind (c,name='OP_div_grad')
  
        implicit none
  
        real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(out  ) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        logical(c_bool), intent(in) :: ccorner_yes

        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        integer                ::  m
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Achtung: - Ein direktes Aufstellen des Laplace-Operators funktioniert i.A. nur sehr schlecht als Vor-    !
        !            konditionierer.                                                                               !
        !          - Das Feld "grad" kann nur eliminiert werden, wenn der div-grad-Stencil komplett ausmultipli-   !
        !            ziert wird (aktuell zu aufwendig).                                                            !
        !----------------------------------------------------------------------------------------------------------!
  
        Lap(S1p:N1p,S2p:N2p,S3p:N3p) = 0
  
        do m = 1, dimens
            call gradient        (m,phi,dig)
            call bc_extrapolation(m,dig    )
            call divergence      (m,dig,Lap)
        !     Lap(S1p:N1p,S2p:N2p,S3p:N3p) = com(S1p:N1p,S2p:N2p,S3p:N3p) + Lap(S1p:N1p,S2p:N2p,S3p:N3p)
        end do
  
    !if (ccorner_yes) call handle_corner_Lap(1,Lap)
    !  if (ccorner_yes) call handle_corner_Lap(1,pre)
  
    !  Lap(:,:,:) = 0 + pre(:,:,:)
  
    end subroutine product_div_grad
  
  
  
  
  
  
  
  
  
  
  
    subroutine product_div_grad_transp(phi,Lap)
  
        implicit none
  
        real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real   , intent(out  ) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk
  
        integer                ::  m
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Achtung: - Alle Operatoren werden transponiert und in umgekehrter Reihenfolge angewendet im Vergleich zu !
        !            "product_div_grad".                                                                           !
        !----------------------------------------------------------------------------------------------------------!
  
  
        if (corner_yes) call handle_corner_Lap(1,phi)
  
        do m = 1, dimens
     
            call divergence_transp      (m,phi,dig)
            call bc_extrapolation_transp(m,dig    )
            call gradient_transp        (m,dig,Lap)
     
        end do
  
  
    end subroutine product_div_grad_transp
  
  
  
  
  
  
  
  
  
  
  
    subroutine product_div_grad_relax(g,phi,Lap) ! TEST!!! aufraeumen und Variablen substituieren ...
  
        implicit none
  
        integer, intent(in   ) ::  g
  
        real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(out  ) ::  Lap(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
        integer                ::  i, N1, S1R, N1R, S11R, N11R
        integer                ::  j, N2, S2R, N2R, S22R, N22R
        integer                ::  k, N3, S3R, N3R, S33R, N33R
  
        integer                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
  
        !*****************************************************************************************
        N1 = NN(1,g)
        N2 = NN(2,g)
        N3 = NN(3,g)
  
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
  
        S1R  = SNB(1,1,g)
        S2R  = SNB(1,2,g)
        S3R  = SNB(1,3,g)
  
        N1R  = SNB(2,1,g)
        N2R  = SNB(2,2,g)
        N3R  = SNB(2,3,g)
  
        S11R = SNF(1,1,g)
        S22R = SNF(1,2,g)
        S33R = SNF(1,3,g)
  
        N11R = SNF(2,1,g)
        N22R = SNF(2,2,g)
        N33R = SNF(2,3,g)
        !*****************************************************************************************
  
  
        !-----------------------------------------------------------------------------------------------------!
        ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!        !
        !            - Direktes Verbinden mit Restriktion erscheint nicht sehr attraktiv, da Lap auf dem      !
        !              jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-     !
        !              Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei    !
        !              der Interpolation verwendet (ebenfalls zur Vereinfachung der Programmierung) und ist   !
        !              somit ohnehin vorhanden. Geschwindigkeitsmässig ist vermutlich auch nicht mehr viel zu !
        !              holen.                                                                                 !
        !-----------------------------------------------------------------------------------------------------!
  
  
        call exchange_relax(g,0,0,0,0,.true.,phi)
  
  
        !===========================================================================================================
        if (dimens == 3) then
     
            do k = S33R, N33R
                do j = S22R, N22R
                    !pgi$ unroll = n:8
                    do i = S11R, N11R
                        Lap(i,j,k) =  cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(1,i,g)*phi(i+1,j,k)   &
                            &  +  cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(1,j,g)*phi(i,j+1,k)   &
                            &  +  cdg3(-1,k,g)*phi(i,j,k-1) + cdg3(1,k,g)*phi(i,j,k+1)   &
                            &  + (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))*phi(i,j,k)
                    end do
                end do
            end do
     
        !===========================================================================================================
        else
     
            do k = S33R, N33R
                do j = S22R, N22R
                    !pgi$ unroll = n:8
                    do i = S11R, N11R
                        Lap(i,j,k) =  cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(1,i,g)*phi(i+1,j,k)   &
                            &  +  cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(1,j,g)*phi(i,j+1,k)   &
                            &  + (cdg1(0,i,g) + cdg2(0,j,g))*phi(i,j,k)
                    end do
                end do
            end do
     
        end if
        !===========================================================================================================
  
  
        !===========================================================================================================
        !=== Randbedingungen =======================================================================================
        !===========================================================================================================
  
        if (BC_1L > 0) then
            i = 1
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do j = S2R, N2R
                    Lap(i,j,k) = cdg1(0,i,g)*phi(i,j,k) + cdg1(1,i,g)*phi(i+1,j,k)
                end do
            end do
        end if
  
        if (BC_1U > 0) then
            i = N1
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do j = S2R, N2R
                    Lap(i,j,k) = cdg1(-1,i,g)*phi(i-1,j,k) + cdg1(0,i,g)*phi(i,j,k)
                end do
            end do
        end if
  
        !===========================================================================================================
  
        if (BC_2L > 0) then
            j = 1
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg2(0,j,g)*phi(i,j,k) + cdg2(1,j,g)*phi(i,j+1,k)
                end do
            end do
        end if
  
        if (BC_2U > 0) then
            j = N2
            do k = S3R, N3R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg2(-1,j,g)*phi(i,j-1,k) + cdg2(0,j,g)*phi(i,j,k)
                end do
            end do
        end if
  
        !===========================================================================================================
  
        if (BC_3L > 0) then
            k = 1
            do j = S2R, N2R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg3(0,k,g)*phi(i,j,k) + cdg3(1,k,g)*phi(i,j,k+1)
                end do
            end do
        end if
  
        if (BC_3U > 0) then
            k = N3
            do j = S2R, N2R
                !pgi$ unroll = n:8
                do i = S1R, N1R
                    Lap(i,j,k) = cdg3(-1,k,g)*phi(i,j,k-1) + cdg3(0,k,g)*phi(i,j,k)
                end do
            end do
        end if
  
        !===========================================================================================================
  
  
        if (corner_yes) call handle_corner_Lap(g,Lap)
  
  
    end subroutine product_div_grad_relax
  
  
  
  
  
  
  
  
  
    ! Habe RB-Linienrelaxation herausgeschmissen, weil die Gewinne in der Geschwindigkeit zu gering waren.
    ! Gleichzeitig war der Code viel zu lang und unuebersichtlich.
    subroutine relaxation_div_grad(init_yes,n_relax,g,bb,rel) ! TEST!!! reine 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
        implicit none
  
        !*************************************************************************************************
        integer                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
        !*************************************************************************************************
  
        logical, intent(in   ) ::  init_yes
        integer, intent(in   ) ::  n_relax
  
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
        integer                ::  i, N1, S1R, N1R, S11R, N11R
        integer                ::  j, N2, S2R, N2R, S22R, N22R
        integer                ::  k, N3, S3R, N3R, S33R, N33R
  
        integer                ::  r, ss
        real                   ::  mult
  
        ! TEST!!! herausziehen?
        ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
        ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
        ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
        integer, parameter     ::  RBGS_mode = 2
        logical, parameter     ::  SOR_yes = .false.
        !REAL   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
        real                   ::  omega
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Optimization annotations:                                                                                !
        !    - line relaxation runs on Rosa (Cray XT5, Istanbul hex-core processors) in single core mode with      !
        !      about 5.6-5.8% peak performance in direction 1, 4.2-4.6% in direction 2, 3.9-4.0% in direction 3.   !
        !    - line relaxation in direction 3 is about 50% more expensive than in direction 1; direction 2 is      !
        !      somewhere in between.                                                                               !
        !    - line relaxation in direction 1 is typically 30% of the total execution time(!).                     !
        !    - RB-GS converges slightly worse than standard GS.                                                    !
        !    - RB-GS allows full vectorization (execpt some operations in direction 1), provided that the loops    !
        !      are reordered to let direction 1 be the innermost (speedup was not as good as hoped (<5%) for       !
        !      direction 1 and 2, direction 3 was not tested).                                                     !
        !    - unrolling the non-vectorizable loops 8 times gives roughly the best results.                        !
        !    - IF statements are automatically pulled out of the loops by the compiler.                            !
        !----------------------------------------------------------------------------------------------------------!
  
        ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
        !*****************************************************************************************
        N1 = NN(1,g)
        N2 = NN(2,g)
        N3 = NN(3,g)
  
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
  
        S1R  = SNB(1,1,g)
        S2R  = SNB(1,2,g)
        S3R  = SNB(1,3,g)
  
        N1R  = SNB(2,1,g)
        N2R  = SNB(2,2,g)
        N3R  = SNB(2,3,g)
  
        S11R = SNF(1,1,g)
        S22R = SNF(1,2,g)
        S33R = SNF(1,3,g)
  
        N11R = SNF(2,1,g)
        N22R = SNF(2,2,g)
        N33R = SNF(2,3,g)
        !*****************************************************************************************
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!           !
        !              - Bei der Initialisierung m�ssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
        !                gew�hlt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
        !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder �berschrieben !
        !                werden. Ansonsten w�re eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
        !                kann aber z.B. zu Floating Point Exceptions f�hren (die f�r die eigentliche Rechnung      !
        !                allerdings irrelevant w�ren)!                                                             !
        !              - obere Stirnfl�chen werden bei der Initialisierung ebenfalls ber�cksichtigt, da ggf.       !
        !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
        !                Seite sonst ber�cksichtigt werden m�sste. ==> geringerer Mehraufwand beim Rechnen,        !
        !                weniger Programmierdetails.                                                               !
        !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
        !                �quidistant ist. Der L�sungsaufwand w�rde sich etwa halbieren!!                           !
        !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
        !                hier aber aus Gr�nden der �bersicht bisher nicht ausgef�hrt wurde.                        !
        !----------------------------------------------------------------------------------------------------------!
  
  
        omega = 1.
  
  
        if (Jacobi_yes) then
  
            do r = 1, n_relax
     
                if (r == 1 .and. init_yes) then
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0 .and. dimens == 3) then
                        k = 1
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                            end do
                        end do
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                            end do
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0 .and. dimens == 3) then
                        k = N3
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                            end do
                        end do
                    end if
                   !=====================================================================================================
        
                else
        
                    call exchange_relax(g,0,0,0,0,.true.,rel)
        
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0 .and. dimens == 3) then
                        k = 1
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k)                                                &
                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0 .and. dimens == 3) then
                        k = N3
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !=====================================================================================================
        
                    rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
                end if
     
            end do
  
  
        !***********************************************************************************************************
        !***********************************************************************************************************
        !***********************************************************************************************************
        else
  
  
            if (init_yes) then
                if (BC_1L <= 0) rel(S11R-1 ,S2R:N2R,S3R:N3R) = 0.
                if (BC_1U <= 0) rel(N11R+1 ,S2R:N2R,S3R:N3R) = 0.
                if (BC_2L <= 0) rel(S1R:N1R,S22R-1 ,S3R:N3R) = 0.
                if (BC_2U <= 0) rel(S1R:N1R,N22R+1 ,S3R:N3R) = 0.
                if (BC_3L <= 0) rel(S1R:N1R,S2R:N2R,S33R-1 ) = 0.
                if (BC_3U <= 0) rel(S1R:N1R,S2R:N2R,N33R+1 ) = 0.
     
               !rel = 0. ! TEST!!!  Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
            end if
  
  
            do r = 1, n_relax
     
                if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
     
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. impl_dir(3) == 0) then
        
                    if (RBGS_mode > 0) then
                        !-------------------------------------------------------------------------------------------------!
                        ! Hinweis:                                                                                        !
                        ! Die Red-Black-Sortierung sollte f�r eine gute Performance global sein, d.h. ueber die Bl�cke    !
                        ! hinweg implementiert sein. Die aktuelle Umsetzung macht sich zunutze, dass im Parallelbetrieb   !
                        ! die Anzahl der Gitterpunkte in jeder Raumrichtung automatisch eine gerade Zahl ist (die Tat-    !
                        ! sache, dass in Randbl�cken ggf. eine ungerade Anzahl Gitterpunkte abgearbeitet wird hat damit   !
                        ! nichts zu tun!). Daher kann sich der Shift "ss" in jedem Block einfach auf den lokalen Gitter-  !
                        ! punkt rel(1,1,1) beziehen.                                                                      !
                        !                                                                                                 !
                        ! Die aktuelle Schreibweise in den Schleifen ist �quivalent zu:                                   !
                        ! DO k = S33R, N33R                                                                               !
                        !    DO j = S22R, N22R                                                                            !
                        !       DO i = S11R, N11R                                                                         !
                        !          IF (MOD(i+j+k) == 1) ... ! (red)                                                       !
                        !          ! bzw.                                                                                 !
                        !          IF (MOD(i+j+k) == 0) ... ! (black)                                                     !
                        !       END DO                                                                                    !
                        !    END DO                                                                                       !
                        ! END DO                                                                                          !
                        !-------------------------------------------------------------------------------------------------!
           
                        !==================================================================================================
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2) ! "+k" ist offenbar elementar f�r grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            else
                                                rel(i,j,k) =       bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                                !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Gl�ttung unabh�ngig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                                !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Gl�ttung unabh�ngig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2) ! "+k" ist offenbar elementar f�r grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            else
                                                rel(i,j,k) =       (bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1)                           )    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1)                           )    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                       !==================================================================================================
                       !==================================================================================================
           
                    else
           
                        !==================================================================================================
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_3L > 0) then
                            k = 1
                            if (r == 1 .and. init_yes) then
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                    end do
                                end do
                            else
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        if (SOR_yes) then
                                            rel(i,j,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                                            rel(i,j,k) = omega*rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        else
                                            rel(i,j,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                                            rel(i,j,k) = rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        end if
                                    end do
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        if (SOR_yes) then
                                            rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                        else
                                            rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        end if
                                    end do
                                end do
                            end do
                        end if
                        !==================================================================================================
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_3U > 0) then
                            k = N3
                            if (r == 1 .and. init_yes) then
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                    end do
                                end do
                            else
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                       !==================================================================================================
                       !==================================================================================================
                    end if
        
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(1) == 1) then
        
                    !=====================================================================================================
                    if (BC_2L > 0) then
                        j = 1
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0) then
                        k = 1
                        if (r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do j = S22R, N22R
              
                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                        &      - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    SOR1(i) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB �berschreiben ----------------------------------------------------------------------
                            if (BC_1L > 0) then
                                band1(1,S1R) = cdg1(0,S1R,g)
                                band1(2,S1R) = bb  (S1R,j,k)
                            else
                                band1(2,S1R) = band1(2,S1R) - cdg1(-1,S1R,g)*rel(S1R-1,j,k)
                            end if
              
                            if (BC_1U > 0) then
                                band1(1,N1R) = cdg1(0,N1R,g)
                                band1(2,N1R) = bb  (N1R,j,k)
                            else
                                band1(2,N1R) = band1(2,N1R) - cdg1( 1,N1R,g)*rel(N1R+1,j,k)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do i = S1R+1, N1R
                                mult       = cdg1(-1,i,g) / band1(1,i-1)
                                band1(1,i) = band1(1,i) - mult * cdg1(1,i-1,g)
                                band1(2,i) = band1(2,i) - mult * band1(2,i-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(N1R,j,k) = band1(2,N1R) / band1(1,N1R)
                            !pgi$ unroll = n:8
                            do i = N1R-1, S1R, -1
                                rel(i,j,k) = (band1(2,i) - cdg1(1,i,g)*rel(i+1,j,k)) / band1(1,i)
                            end do

              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then ! TEST!!! Das koennte man auch besser mit der vorherigen Schleife verbinden ...
                                if (r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do i = S1R, N1R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do i = S1R, N1R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR1(i)
                                    end do
                                end if
                            end if
              
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_2U > 0) then
                        j = N2
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0) then
                        k = N3
                        if (r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(2) == 1) then
        
                    if (.not. (impl_dir(1) == 0 .and. r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel) ! TEST!!! Zwischenaktualisierung auch bei Helmholtz-Routinen einbauen!
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0) then
                        k = 1
                        if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do i = S11R, N11R
              
                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                        &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (impl_dir(1) == 0 .and. r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    SOR2(j) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB �berschreiben ----------------------------------------------------------------------
                            if (BC_2L > 0) then
                                band2(1,S2R) = cdg2(0,S2R,g)
                                band2(2,S2R) = bb  (i,S2R,k)
                            else
                                band2(2,S2R) = band2(2,S2R) - cdg2(-1,S2R,g)*rel(i,S2R-1,k)
                            end if
              
                            if (BC_2U > 0) then
                                band2(1,N2R) = cdg2(0,N2R,g)
                                band2(2,N2R) = bb  (i,N2R,k)
                            else
                                band2(2,N2R) = band2(2,N2R) - cdg2( 1,N2R,g)*rel(i,N2R+1,k)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do j = S2R+1, N2R
                                mult       = cdg2(-1,j,g) / band2(1,j-1)
                                band2(1,j) = band2(1,j) - mult * cdg2(1,j-1,g)
                                band2(2,j) = band2(2,j) - mult * band2(2,j-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(i,N2R,k) = band2(2,N2R) / band2(1,N2R)
                            !pgi$ unroll = n:8
                            do j = N2R-1, S2R, -1
                                rel(i,j,k) = (band2(2,j) - cdg2(1,j,g)*rel(i,j+1,k)) / band2(1,j)
                            end do
              
              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then
                                if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do j = S2R, N2R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do j = S2R, N2R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR2(j)
                                    end do
                                end if
                            end if

                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0) then
                        k = N3
                        if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================
        
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(3) == 1) then
        
                    if (.not. (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do j = S22R, N22R
                        do i = S11R, N11R
              
                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)   &
                                        &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    SOR3(k) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB �berschreiben ----------------------------------------------------------------------
                            if (BC_3L > 0) then
                                band3(1,S3R) = cdg3(0,S3R,g)
                                band3(2,S3R) = bb  (i,j,S3R)
                            else
                                band3(2,S3R) = band3(2,S3R) - cdg3(-1,S3R,g)*rel(i,j,S3R-1)
                            end if
              
                            if (BC_3U > 0) then
                                band3(1,N3R) = cdg3(0,N3R,g)
                                band3(2,N3R) = bb  (i,j,N3R)
                            else
                                band3(2,N3R) = band3(2,N3R) - cdg3( 1,N3R,g)*rel(i,j,N3R+1)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do k = S3R+1, N3R
                                mult       = cdg3(-1,k,g) / band3(1,k-1)
                                band3(1,k) = band3(1,k) - mult * cdg3(1,k-1,g)
                                band3(2,k) = band3(2,k) - mult * band3(2,k-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(i,j,N3R) = band3(2,N3R) / band3(1,N3R)
                            !pgi$ unroll = n:8
                            do k = N3R-1, S3R, -1
                                rel(i,j,k) = (band3(2,k) - cdg3(1,k,g)*rel(i,j,k+1)) / band3(1,k)
                            end do
              
              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then
                                if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do k = S3R, N3R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do k = S3R, N3R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR3(k)
                                    end do
                                end if
                            end if

                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================

                end if
               !========================================================================================================
               !========================================================================================================
               !========================================================================================================

            end do
        !***********************************************************************************************************
        !***********************************************************************************************************
        !***********************************************************************************************************

        end if
  
  
        if (corner_yes) call handle_corner_Lap(g,rel)


    end subroutine relaxation_div_grad
  
  
  
  
  
  
  
  
  
  
  
    subroutine relaxation_div_grad_inv(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ... ! TEST!!! aufraeumen und Variablen substituieren ...
  
        implicit none
  
        !*************************************************************************************************
        integer                ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U
        !*************************************************************************************************
  
        logical, intent(in   ) ::  init_yes
        integer, intent(in   ) ::  n_relax
  
        integer, intent(in   ) ::  g
  
        real   , intent(in   ) ::  bb  (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real   , intent(inout) ::  rel (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
        real                   ::  comp(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U)) ! FELD!!!
  
        integer                ::  i, N1, S1R, N1R, S11R, N11R
        integer                ::  j, N2, S2R, N2R, S22R, N22R
        integer                ::  k, N3, S3R, N3R, S33R, N33R
  
        integer                ::  r, ss
        real                   ::  mult
  
        ! TEST!!! herausziehen?
        ! RBGS_mode =< 0 : naturally ordered Gauss-Seidel (slow on CPU, cannot be vectorized)
        ! RBGS_mode == 1 : 3D-Red-Black Gauss-Seidel (faster on CPU due to vectorization, parallelization-independent)
        ! RBGS_mode >= 2 : 2D-Red-Black Gauss-Seidel (normally the fastest variant on CPU due to better cache utilization)
        integer, parameter     ::  RBGS_mode = 2
        logical, parameter     ::  SOR_yes = .false.
        !REAL   , PARAMETER     ::  omega   = 0.8! 1.2 ! 1.27 !1.27  ! omega = 0.9
        real                   ::  omega
  
        ! Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden aus dem Loop herauskopiert, muss man nicht von Hand machen!
  
  
        !*****************************************************************************************
        N1 = NN(1,g)
        N2 = NN(2,g)
        N3 = NN(3,g)
  
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
  
        S1R  = SNB(1,1,g)
        S2R  = SNB(1,2,g)
        S3R  = SNB(1,3,g)
  
        N1R  = SNB(2,1,g)
        N2R  = SNB(2,2,g)
        N3R  = SNB(2,3,g)
  
        S11R = SNF(1,1,g)
        S22R = SNF(1,2,g)
        S33R = SNF(1,3,g)
  
        N11R = SNF(2,1,g)
        N22R = SNF(2,2,g)
        N33R = SNF(2,3,g)
        !*****************************************************************************************
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!           !
        !              - Bei der Initialisierung m�ssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
        !                gew�hlt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
        !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder �berschrieben !
        !                werden. Ansonsten w�re eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
        !                kann aber z.B. zu Floating Point Exceptions f�hren (die f�r die eigentliche Rechnung      !
        !                allerdings irrelevant w�ren)!                                                             !
        !              - obere Stirnfl�chen werden bei der Initialisierung ebenfalls ber�cksichtigt, da ggf.       !
        !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
        !                Seite sonst ber�cksichtigt werden m�sste. ==> geringerer Mehraufwand beim Rechnen,        !
        !                weniger Programmierdetails.                                                               !
        !              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn mindestens eine Richtung     !
        !                �quidistant ist. Der L�sungsaufwand w�rde sich etwa halbieren!!                           !
        !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
        !                hier aber aus Gr�nden der �bersicht bisher nicht ausgef�hrt wurde.                        !
        !----------------------------------------------------------------------------------------------------------!
  
  
        omega = 1.
  
  
        if (Jacobi_yes) then
  
            do r = 1, n_relax
     
                if (r == 1 .and. init_yes) then
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0 .and. dimens == 3) then
                        k = 1
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                            end do
                        end do
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                            end do
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0 .and. dimens == 3) then
                        k = N3
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                            end do
                        end do
                    end if
                   !=====================================================================================================
        
                else
        
                    call exchange_relax(g,0,0,0,0,.true.,rel)
        
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0 .and. dimens == 3) then
                        k = 1
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !=====================================================================================================
                    do k = S33R, N33R
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k)                                                &
                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do j = S22R, N22R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        do k = S33R, N33R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0 .and. dimens == 3) then
                        k = N3
                        do j = S22R, N22R
                            !pgi$ unroll = n:8
                            do i = S11R, N11R
                                comp(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                            end do
                        end do
                    end if
                    !=====================================================================================================
        
                    rel(1:N1,1:N2,1:N3) = comp(1:N1,1:N2,1:N3)
        
                end if
     
            end do
  
  
        !***********************************************************************************************************
        !***********************************************************************************************************
        !***********************************************************************************************************
        else
  
  
            if (init_yes) then
                if (BC_1L <= 0) rel(S11R-1 ,S2R:N2R,S3R:N3R) = 0.
                if (BC_1U <= 0) rel(N11R+1 ,S2R:N2R,S3R:N3R) = 0.
                if (BC_2L <= 0) rel(S1R:N1R,S22R-1 ,S3R:N3R) = 0.
                if (BC_2U <= 0) rel(S1R:N1R,N22R+1 ,S3R:N3R) = 0.
                if (BC_3L <= 0) rel(S1R:N1R,S2R:N2R,S33R-1 ) = 0.
                if (BC_3U <= 0) rel(S1R:N1R,S2R:N2R,N33R+1 ) = 0.
     
               !rel = 0. ! TEST!!!  Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
            end if
  
  
            do r = 1, n_relax
     
                if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
     
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. impl_dir(3) == 0) then
        
                    if (RBGS_mode > 0) then
                        !-------------------------------------------------------------------------------------------------!
                        ! Hinweis:                                                                                        !
                        ! Die Red-Black-Sortierung sollte f�r eine gute Performance global sein, d.h. ueber die Bl�cke    !
                        ! hinweg implementiert sein. Die aktuelle Umsetzung macht sich zunutze, dass im Parallelbetrieb   !
                        ! die Anzahl der Gitterpunkte in jeder Raumrichtung automatisch eine gerade Zahl ist (die Tat-    !
                        ! sache, dass in Randbl�cken ggf. eine ungerade Anzahl Gitterpunkte abgearbeitet wird hat damit   !
                        ! nichts zu tun!). Daher kann sich der Shift "ss" in jedem Block einfach auf den lokalen Gitter-  !
                        ! punkt rel(1,1,1) beziehen.                                                                      !
                        !                                                                                                 !
                        ! Die aktuelle Schreibweise in den Schleifen ist �quivalent zu:                                   !
                        ! DO k = S33R, N33R                                                                               !
                        !    DO j = S22R, N22R                                                                            !
                        !       DO i = S11R, N11R                                                                         !
                        !          IF (MOD(i+j+k) == 1) ... ! (red)                                                       !
                        !          ! bzw.                                                                                 !
                        !          IF (MOD(i+j+k) == 0) ... ! (black)                                                     !
                        !       END DO                                                                                    !
                        !    END DO                                                                                       !
                        ! END DO                                                                                          !
                        !-------------------------------------------------------------------------------------------------!
           
                        !==================================================================================================
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2) ! "+k" ist offenbar elementar f�r grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            else
                                                rel(i,j,k) =       bb(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                                !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Gl�ttung unabh�ngig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                                !CALL exchange_relax(g,0,0,0,0,.TRUE.,rel) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Gl�ttung unabh�ngig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
                                do k = S33R, N33R
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                do k = N33R, S33R, -1
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2) ! "+k" ist offenbar elementar f�r grosse Konvergenzrate. Symmetrie zu relaxation_div_grad_inv ist hier interessanterweise kontraproduktiv!
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            else
                                                rel(i,j,k) =       (bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &                                  - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &                                  - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            else
                                do k = N33R, S33R, -1
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            if (SOR_yes) then
                                                rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                            else
                                                rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                    &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                    &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                    &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                    &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                            end if
                                        end do
                                    end do
                                end do
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(i+k+S22R+1,2)
                                    !pgi$ unroll = n:8
                                    do j = ss+S22R, N22R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    ss = MOD(j+k+S11R+1,2)
                                    !pgi$ unroll = n:8
                                    do i = ss+S11R, N11R, 2
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (RBGS_mode == 1) then
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3U > 0) then
                                k = N3
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        ss = MOD(j+k+S11R+1,2)
                                        !pgi$ unroll = n:8
                                        do i = ss+S11R, N11R, 2
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        else
                            !--------------------------------------------------------------------------------------------------
                            if (BC_3L > 0) then
                                k = 1
                                if (r == 1 .and. init_yes) then
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                        end do
                                    end do
                                else
                                    do j = S22R, N22R
                                        !pgi$ unroll = n:8
                                        do i = S11R, N11R
                                            rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                        end do
                                    end do
                                end if
                            end if
                        !--------------------------------------------------------------------------------------------------
                        end if
                       !==================================================================================================
                       !==================================================================================================

                    else

                        !==================================================================================================
                        !==================================================================================================
                        if (BC_1U > 0) then
                            i = N1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2U > 0) then
                            j = N2
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_3U > 0) then
                            k = N3
                            if (r == 1 .and. init_yes) then
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                    end do
                                end do
                            else
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !==================================================================================================
                        if (r == 1 .and. init_yes) then
                            do k = N33R, S33R, -1
                                do j = N22R, S22R, -1
                                    !pgi$ unroll = n:8
                                    do i = N11R, S11R, -1
                                        if (SOR_yes) then
                                            rel(i,j,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                            rel(i,j,k) = omega*rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        else
                                            rel(i,j,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                            rel(i,j,k) = rel(i,j,k) / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        end if
                                    end do
                                end do
                            end do
                        else
                            do k = N33R, S33R, -1
                                do j = N22R, S22R, -1
                                    !pgi$ unroll = n:8
                                    do i = N11R, S11R, -1
                                        if (SOR_yes) then
                                            rel(i,j,k) = omega*(bb(i,j,k)                                                 &
                                                &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)) + (1.-omega)*rel(i,j,k)
                                        else
                                            rel(i,j,k) =       (bb(i,j,k)                                                 &
                                                &      - cdg1(-1,i,g)*rel(i-1,j,k) - cdg1(1,i,g)*rel(i+1,j,k)     &
                                                &      - cdg2(-1,j,g)*rel(i,j-1,k) - cdg2(1,j,g)*rel(i,j+1,k)     &
                                                &      - cdg3(-1,k,g)*rel(i,j,k-1) - cdg3(1,k,g)*rel(i,j,k+1))    &
                                                &      / (cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g))
                                        end if
                                    end do
                                end do
                            end do
                        end if
                        !==================================================================================================
                        if (BC_1L > 0) then
                            i = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do j = S22R, N22R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_2L > 0) then
                            j = 1
                            if (r == 1 .and. init_yes) then
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                                    end do
                                end do
                            else
                                do k = S33R, N33R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                        !--------------------------------------------------------------------------------------------------
                        if (BC_3L > 0) then
                            k = 1
                            if (r == 1 .and. init_yes) then
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                    end do
                                end do
                            else
                                do j = S22R, N22R
                                    !pgi$ unroll = n:8
                                    do i = S11R, N11R
                                        rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                    end do
                                end do
                            end if
                        end if
                       !==================================================================================================
                       !==================================================================================================
                    end if
        
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(3) == 1) then
        
                    if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2U > 0) then
                        j = N2
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do j = N22R, S22R, -1
                        do i = N11R, S11R, -1
              
                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band3(2,k) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    band3(1,k) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band3(2,k) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)   &
                                        &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do k = S3R, N3R
                                    SOR3(k) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB �berschreiben ----------------------------------------------------------------------
                            if (BC_3L > 0) then
                                band3(1,S3R) = cdg3(0,S3R,g)
                                band3(2,S3R) = bb  (i,j,S3R)
                            else
                                band3(2,S3R) = band3(2,S3R) - cdg3(-1,S3R,g)*rel(i,j,S3R-1)
                            end if
              
                            if (BC_3U > 0) then
                                band3(1,N3R) = cdg3(0,N3R,g)
                                band3(2,N3R) = bb  (i,j,N3R)
                            else
                                band3(2,N3R) = band3(2,N3R) - cdg3( 1,N3R,g)*rel(i,j,N3R+1)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do k = S3R+1, N3R
                                mult       = cdg3(-1,k,g) / band3(1,k-1)
                                band3(1,k) = band3(1,k) - mult * cdg3(1,k-1,g)
                                band3(2,k) = band3(2,k) - mult * band3(2,k-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(i,j,N3R) = band3(2,N3R) / band3(1,N3R)
                            !pgi$ unroll = n:8
                            do k = N3R-1, S3R, -1
                                rel(i,j,k) = (band3(2,k) - cdg3(1,k,g)*rel(i,j,k+1)) / band3(1,k)
                            end do
              
              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then
                                if (r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do k = S3R, N3R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do k = S3R, N3R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR3(k)
                                    end do
                                end if
                            end if

                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_2L > 0) then
                        j = 1
                        if (r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================
        
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(2) == 1) then
        
                    if (.not. (impl_dir(3) == 0 .and. r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
                    !=====================================================================================================
                    if (BC_1U > 0) then
                        i = N1
                        if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0) then
                        k = N3
                        if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do k = N33R, S33R, -1
                        do i = N11R, S11R, -1
              
                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band2(2,j) = bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    band2(1,j) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band2(2,j) = bb(i,j,k) - cdg1(-1,i,g)*rel(i-1,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                        &      - cdg1( 1,i,g)*rel(i+1,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (impl_dir(3) == 0 .and. r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do j = S2R, N2R
                                    SOR2(j) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB �berschreiben ----------------------------------------------------------------------
                            if (BC_2L > 0) then
                                band2(1,S2R) = cdg2(0,S2R,g)
                                band2(2,S2R) = bb  (i,S2R,k)
                            else
                                band2(2,S2R) = band2(2,S2R) - cdg2(-1,S2R,g)*rel(i,S2R-1,k)
                            end if
              
                            if (BC_2U > 0) then
                                band2(1,N2R) = cdg2(0,N2R,g)
                                band2(2,N2R) = bb  (i,N2R,k)
                            else
                                band2(2,N2R) = band2(2,N2R) - cdg2( 1,N2R,g)*rel(i,N2R+1,k)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do j = S2R+1, N2R
                                mult       = cdg2(-1,j,g) / band2(1,j-1)
                                band2(1,j) = band2(1,j) - mult * cdg2(1,j-1,g)
                                band2(2,j) = band2(2,j) - mult * band2(2,j-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(i,N2R,k) = band2(2,N2R) / band2(1,N2R)
                            !pgi$ unroll = n:8
                            do j = N2R-1, S2R, -1
                                rel(i,j,k) = (band2(2,j) - cdg2(1,j,g)*rel(i,j+1,k)) / band2(1,j)
                            end do
              
              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then
                                if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do j = S2R, N2R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do j = S2R, N2R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR2(j)
                                    end do
                                end if
                            end if

                        end do
                    end do
                    !=====================================================================================================
                    if (BC_1L > 0) then
                        i = 1
                        if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do j = S22R, N22R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i,g)*rel(i+1,j,k)) / cdg1(0,i,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0) then
                        k = 1
                        if (impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================
        
                end if
                !========================================================================================================
                !========================================================================================================
                !========================================================================================================
                if (impl_dir(1) == 1) then
        
                    !=====================================================================================================
                    if (BC_2U > 0) then
                        j = N2
                        if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3U > 0) then
                        k = N3
                        if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*bb(i,j,k) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k,g)*rel(i,j,k-1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !=====================================================================================================
                    do k = N33R, S33R, -1
                        do j = N22R, S22R, -1

                            !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
                            if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band1(2,i) = bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            else
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    band1(1,i) = cdg1(0,i,g) + cdg2(0,j,g) + cdg3(0,k,g)
                                    band1(2,i) = bb(i,j,k) - cdg2(-1,j,g)*rel(i,j-1,k) - cdg3(-1,k,g)*rel(i,j,k-1)   &
                                        &      - cdg2( 1,j,g)*rel(i,j+1,k) - cdg3( 1,k,g)*rel(i,j,k+1)
                                end do
                            end if
              
                            !--- SOR (1) -----------------------------------------------------------------------------------
                            if (SOR_yes .and. .not. (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes)) then
                                !pgi$ unroll = n:8
                                do i = S1R, N1R
                                    SOR1(i) = rel(i,j,k)
                                end do
                            end if
              
              
                            !--- mit RB überschreiben ----------------------------------------------------------------------
                            if (BC_1L > 0) then
                                band1(1,S1R) = cdg1(0,S1R,g)
                                band1(2,S1R) = bb  (S1R,j,k)
                            else
                                band1(2,S1R) = band1(2,S1R) - cdg1(-1,S1R,g)*rel(S1R-1,j,k)
                            end if
              
                            if (BC_1U > 0) then
                                band1(1,N1R) = cdg1(0,N1R,g)
                                band1(2,N1R) = bb  (N1R,j,k)
                            else
                                band1(2,N1R) = band1(2,N1R) - cdg1( 1,N1R,g)*rel(N1R+1,j,k)
                            end if
              
              
                            !--- Gauss-Elimination hoch --------------------------------------------------------------------
                            !pgi$ unroll = n:8
                            do i = S1R+1, N1R
                                mult       = cdg1(-1,i,g) / band1(1,i-1)
                                band1(1,i) = band1(1,i) - mult * cdg1(1,i-1,g)
                                band1(2,i) = band1(2,i) - mult * band1(2,i-1)
                            end do
              
              
                            !--- Gauss-Elimination runter ------------------------------------------------------------------
                            rel(N1R,j,k) = band1(2,N1R) / band1(1,N1R)
                            !pgi$ unroll = n:8
                            do i = N1R-1, S1R, -1
                                rel(i,j,k) = (band1(2,i) - cdg1(1,i,g)*rel(i+1,j,k)) / band1(1,i)
                            end do
              
              
                            !--- SOR (2) -----------------------------------------------------------------------------------
                            if (SOR_yes) then
                                if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                                    !pgi$ unroll = n:8
                                    do i = S1R, N1R
                                        rel(i,j,k) = omega*rel(i,j,k)
                                    end do
                                else
                                    !pgi$ unroll = n:8
                                    do i = S1R, N1R
                                        rel(i,j,k) = omega*rel(i,j,k) + (1.-omega)*SOR1(i)
                                    end do
                                end if
                            end if
              
                        end do
                    end do
                    !=====================================================================================================
                    if (BC_2L > 0) then
                        j = 1
                        if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g)
                                end do
                            end do
                        else
                            do k = S33R, N33R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j,g)*rel(i,j+1,k)) / cdg2(0,j,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                    !-----------------------------------------------------------------------------------------------------
                    if (BC_3L > 0) then
                        k = 1
                        if (impl_dir(2) == 0 .and. impl_dir(3) == 0 .and. r == 1 .and. init_yes) then
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g)
                                end do
                            end do
                        else
                            do j = S22R, N22R
                                !pgi$ unroll = n:8
                                do i = S11R, N11R
                                    rel(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k,g)*rel(i,j,k+1)) / cdg3(0,k,g) + (1.-omega)*rel(i,j,k)
                                end do
                            end do
                        end if
                    end if
                   !=====================================================================================================
                end if
               !========================================================================================================
               !========================================================================================================
               !========================================================================================================
     
            end do
        !***********************************************************************************************************
        !***********************************************************************************************************
        !***********************************************************************************************************
  
        end if
  
  
        if (corner_yes) call handle_corner_Lap(g,rel)
  
  
    end subroutine relaxation_div_grad_inv
  
  
  
  
  
  
  
  
  
  
  
    subroutine handle_corner_Lap(g,phi)
  
        implicit none
  
        integer, intent(in   ) ::  g
        real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  
        !----------------------------------------------------------------------------------------------------------!
        ! Anmerkungen: - Diese Routine dient dazu, den unbestimmten Druck in den Ecken und Kanten explizit zu      !
        !                behandeln und gleichzeitig das Konvergenzverhalten der Löser (BiCGstab oder Richardson)   !
        !                möglichst nicht zu beeinträchtigen.                                                       !
        !              - Der Druck wird hier direkt zu Null gesetzt, um innerhalb des iterativen Lösers kein Re-   !
        !                siduum zu erzeugen (anstelle über die RHS).                                               !
        !              - Der Druck wird erst bei Bedarf (z.B. vor einem Ausschrieb) auf einen sinnvollen Wert ge-  !
        !                setzt.                                                                                    !
        !              - Siehe dazu auch die korrespondierende Subroutine "handle_corner_rhs"!                     !
        !----------------------------------------------------------------------------------------------------------!
  
  
        if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) phi(1      ,1      ,1:NN(3,g)) = 0. ! TEST!!! verifizieren ...
        if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) phi(1      ,NN(2,g),1:NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) phi(NN(1,g),1      ,1:NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) phi(NN(1,g),NN(2,g),1:NN(3,g)) = 0.
  
        if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) phi(1      ,1:NN(2,g),1      ) = 0.
        if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) phi(1      ,1:NN(2,g),NN(3,g)) = 0.
        if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) phi(NN(1,g),1:NN(2,g),1      ) = 0.
        if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) phi(NN(1,g),1:NN(2,g),NN(3,g)) = 0.
  
        if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) phi(1:NN(1,g),1      ,1      ) = 0.
        if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) phi(1:NN(1,g),1      ,NN(3,g)) = 0.
        if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) phi(1:NN(1,g),NN(2,g),1      ) = 0.
        if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) phi(1:NN(1,g),NN(2,g),NN(3,g)) = 0.
  
  
    end subroutine handle_corner_Lap
  
  
  
  
end module cmod_laplace
