!> \brief module providing functions to initiliaze and apply RestrictionOp
module cmod_InterpolationOp

    use iso_c_binding

contains

    subroutine MG_getCIS(   &
        N,                  &
        bL, bU,             &
        xs,                 &
        cI ) bind(c,name='MG_getCIS')

        implicit none

        integer(c_int), intent(in)  :: N
        integer(c_int), intent(in)  :: bL, bU

        real(c_double), intent(in)  :: xs(bL:(N+bU))

        real(c_double), intent(out) :: cI(1:2,1:N)

        integer(c_int)              :: i
        real(c_double)              :: Dx10, Dx12

        cI = 0.

        !===========================================================================================================
        !=== Interpolation, linienweise, 1d ========================================================================
        !===========================================================================================================
        ! coarse
        !    xs(i-1)                        xs(i)
        !  ----o------------------------------o------------------------------o-----------------
        !      |------Dx10----|
        !      |------------Dx12--------------|
        !  ----o--------------o---------------o-----------o-----
        !     xs(i-1)        xs(i)           xs(i+1)
        ! fine


        !--------------------------------------------------------------------------------------------------------
        do i = 1, N
            Dx10 = xs(i  )-xs(i-1)
            Dx12 = xs(i+1)-xs(i-1)

            cI(1,i) = 1.- Dx10/Dx12
            cI(2,i) =     Dx10/Dx12
        end do
        !--------------------------------------------------------------------------------------------------------
    end subroutine MG_getCIS



    !> \todo understand from where comes the weith 0.75 / 0.25
    subroutine MG_getCIV(   &
        Nc,                 &
        bLc, bUc,           &
        SSc,                &
        NNc,                &
        BC_L, BC_U,         &
        Nf,                 &
        bLf, bUf,           &
        SSf,                &
        !        NNf,                &
        xc,xf,              &
        cIV ) bind(c,name='MG_getCIV')


        implicit none


        integer(c_int), intent(in)    :: Nc

        integer(c_int), intent(in)    :: bLc,bUc

        integer(c_int), intent(in)    :: SSc
        integer(c_int), intent(in)    :: NNc

        integer(c_int), intent(in)    :: BC_L,BC_U

        integer(c_int), intent(in)    :: Nf

        integer(c_int), intent(in)    :: bLf,bUf

        integer(c_int), intent(in)    :: SSf
        !        integer(c_int), intent(in)    :: NNf


        real(c_double), intent(in)    :: xc( bLc:(Nc+bUc) )

        real(c_double), intent(in)    :: xf( bLf:(Nf+bUf) )

        real(c_double), intent(inout) :: cIV( 1:2,SSf:NNf)

        integer(c_int)                ::  i, ic, dd
        real(c_double)                ::  Dx1a, Dx12


        !===========================================================================================================
        !=== Interpolation, linienweise, 1d ========================================================================
        !===========================================================================================================
        ! fine
        !    xf(i-2)        xf(i-1)        xf(i)       xf(i+1)      xf(i+2)
        !  ---->-------------->--------------->----------->----------->------------
        !                     |--------Dx1a---------|
        !                                     |Dx1a-|
        !                                           |Dx1a-|
        !                                           |-------Dx1a------|
        !             |-----------Dx12--------------|-----------Dx12--------------|
        !  ----------->----------------------------->----------------------------->----------
        !          xc(i-1)                        xc(i)                         xc(i+1)
        ! coarse

        cIV = 0.

        dd = (Nf-1)/(Nc-1)


        do i = SSf, NNf
            ic = (i-SSf)/dd + SSc

            Dx1a = xc(ic)-xf(i -1)
            Dx12 = xc(ic)-xc(ic-1)

            cIV(1,ic) = 1-Dx1a/Dx12

            Dx1a = xc(ic)-xf(i)

            cIV(0,ic) = 1-Dx1a/Dx12

            Dx1a = xf(i +1)-xc(ic)
            Dx12 = xc(ic+1)-xc(ic)

            cIV(1,ic) = 1-Dx1a/Dx12

            Dx1a = xf(i+2)-xc(ic)

            cIV(2,ic) = 1-Dx1a/Dx12

        end do

        ! a little bit shaky, please verify this when IO is ready
        if (BC_L > 0) then
            cIV( :,  0) = 0.
            cIV(-1,SSc) = 0.
            cIV( 0,SSc) = 1.
        end if

        ! maybe false
        if (BC_L == -2) then
            cIV(-1:0,SSc) = 0.
            cIV(  : ,  0) = 0.
        end if

        if (BC_U > 0) then
            cIV(1,NNc) = 1.
            cIV(2,NNc) = 0.

            cIV(:, Nc) = 0.

        end if
        if (BC_U == -2) then
            cIV(1:2,NNc) = 0.
            cIV(  :, Nc) = 0.
        end if

    end subroutine MG_getCIV



    !----------------------------------------------------------------------------------------------------------!
    !> \note:       - für allgemeine dd geeignet                                                        !
    !!              - dd(i) /= 1 <===> N(i) /= 1                                                                     !
    !!              - es wird nur in eine Richung ausgetauscht                                                  !
    !!              - Null-Setzen am Rand nicht notwendig                                                       !
    !!              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
    !!                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
    !!              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
    !!                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)).              !
    !!              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
    !!                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
    !!              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
    !!                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
    !!                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
    !!                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
    !!----------------------------------------------------------------------------------------------------------!
    !! \attention: Verwendung von parametrischen Strides dd in Schleifen verhindert die Vektorisierung bzw.
    !!            das Prefetching! Andererseits ist der Geschwindigkeitsgewinn nur sehr gering (schon getestet).
    !!
    !! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die Speicherzugriffszeit limitiert
    !! - Das wird z.B. deutlich bei Single- vs. Dualcorebetrieb
    subroutine MG_interpolate( &
        dimens,             &
        Nc,                 &
        bLc,bUc,            &
        BCL, BCU,           &
        Nf,                 &
        bLf,bUf,            &
        cI1,cI2,cI3,        &
        phic,               &
        phif ) bind(c,name='MG_interpolate')

        implicit none

        integer(c_int), intent(in)     :: dimens

        integer(c_int), intent(in)     :: Nc(3)

        integer(c_int), intent(in)     :: bLc(3)
        integer(c_int), intent(in)     :: bUc(3)

        integer(c_int), intent(in)     :: BCL(3)
        integer(c_int), intent(in)     :: BCU(3)

        integer(c_int), intent(in)     :: Nf(3)

        integer(c_int), intent(in)     :: bLf(3)
        integer(c_int), intent(in)     :: bUf(3)

        real(c_double),  intent(in)    :: cI1 ( 1:2, 1:Nc(1) )
        real(c_double),  intent(in)    :: cI2 ( 1:2, 1:Nc(2) )
        real(c_double),  intent(in)    :: cI3 ( 1:2, 1:Nc(3) )

        real(c_double),  intent(inout) :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

        real(c_double),  intent(out)   :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))



        integer                ::  i
        integer                ::  j
        integer                ::  k

        integer                :: dd(3)


        dd = 1
        do i=1,dimens
            dd(i) = ( Nf(i)-1 )/( Nc(i)-1 )
        end do


        if (BCL(1) > 0 .and. BCL(2) > 0) phic( 1,     1,     1:Nc(3) ) = ( phic(1+1,     1,     1:Nc(3)) + phic(1    , 1+1,     1:Nc(3)) + phic(1+1,     1+1,     1:Nc(3)) )/3.
        if (BCL(1) > 0 .and. BCU(2) > 0) phic( 1,     Nc(2), 1:Nc(3) ) = ( phic(1+1,     Nc(2), 1:Nc(3)) + phic(1    , Nc(2)-1, 1:Nc(3)) + phic(1+1,     Nc(2)-1, 1:Nc(3)) )/3.
        if (BCU(1) > 0 .and. BCL(2) > 0) phic( Nc(1), 1,     1:Nc(3) ) = ( phic(Nc(1)-1, 1,     1:Nc(3)) + phic(Nc(1), 1+1,     1:Nc(3)) + phic(Nc(1)-1, 1+1,     1:Nc(3)) )/3.
        if (BCU(1) > 0 .and. BCU(2) > 0) phic( Nc(1), Nc(2), 1:Nc(3) ) = ( phic(Nc(1)-1, Nc(2), 1:Nc(3)) + phic(Nc(1), Nc(2)-1, 1:Nc(3)) + phic(Nc(1)-1, Nc(2)-1, 1:Nc(3)) )/3.

        if (BCL(1) > 0 .and. BCL(3) > 0) phic( 1,     1:Nc(2), 1 )     = ( phic(1+1,     1:Nc(2), 1    ) + phic(1,     1:Nc(2), 1+1)     + phic(1+1,     1:Nc(2), 1+1)     )/3.
        if (BCL(1) > 0 .and. BCU(3) > 0) phic( 1,     1:Nc(2), Nc(3) ) = ( phic(1+1,     1:Nc(2), Nc(3)) + phic(1,     1:Nc(2), Nc(3)-1) + phic(1+1,     1:Nc(2), Nc(3)-1) )/3.
        if (BCU(1) > 0 .and. BCL(3) > 0) phic( Nc(1), 1:Nc(2), 1 )     = ( phic(Nc(1)-1, 1:Nc(2), 1    ) + phic(Nc(1), 1:Nc(2), 1+1)     + phic(Nc(1)-1, 1:Nc(2), 1+1)     )/3.
        if (BCU(1) > 0 .and. BCU(3) > 0) phic( Nc(1), 1:Nc(2), Nc(3) ) = ( phic(Nc(1)-1, 1:Nc(2), Nc(3)) + phic(Nc(1), 1:Nc(2), Nc(3)-1) + phic(Nc(1)-1, 1:Nc(2), Nc(3)-1) )/3.

        if (BCL(2) > 0 .and. BCL(3) > 0) phic( 1:Nc(1), 1,     1     ) = ( phic(1:Nc(1), 1+1,     1    ) + phic(1:Nc(1), 1,     1+1    ) + phic(1:Nc(1), 1+1,     1+1)     )/3.
        if (BCL(2) > 0 .and. BCU(3) > 0) phic( 1:Nc(1), 1,     Nc(3) ) = ( phic(1:Nc(1), 1+1,     Nc(3)) + phic(1:Nc(1), 1,     Nc(3)-1) + phic(1:Nc(1), 1+1,     Nc(3)-1) )/3.
        if (BCU(2) > 0 .and. BCL(3) > 0) phic( 1:Nc(1), Nc(2), 1     ) = ( phic(1:Nc(1), Nc(2)-1, 1    ) + phic(1:Nc(1), Nc(2), 1+1    ) + phic(1:Nc(1), Nc(2)-1, 1+1)     )/3.
        if (BCU(2) > 0 .and. BCU(3) > 0) phic( 1:Nc(1), Nc(2), Nc(3) ) = ( phic(1:Nc(1), Nc(2)-1, Nc(3)) + phic(1:Nc(1), Nc(2), Nc(3)-1) + phic(1:Nc(1), Nc(2)-1, Nc(3)-1) )/3.



        !***********************************************************************************************************

        !pgi$ unroll = n:8
        phif( 1:Nf(1):dd(1), 1:Nf(2):dd(2), 1:Nf(3):dd(3) )  =  phic( 1:Nc(1), 1:Nc(2), 1:Nc(3) )

        !***********************************************************************************************************

        !===========================================================================================================
        if( dd(3) /= 1 ) then ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!

            do k = 1+1, Nf(3)-1, dd(3)
                do j = 1, Nf(2), dd(2)
                    !pgi$ unroll = n:8
                    do i = 1, Nf(1), dd(1)
                        phif(i,j,k) = cI3(1,k/dd(3))*phif(i,j,k-1) + cI3(2,k/dd(3))*phif(i,j,k+1)
                    end do
                end do
            end do

        end if
        !===========================================================================================================
        if( dd(2) /= 1 ) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (NNf(3) == 2??)

            do k = 1, Nf(3)
                do j = 1+1, Nf(2)-1, dd(2)
                    !pgi$ unroll = n:8
                    do i = 1, Nf(1), dd(1)
                        phif(i,j,k) = cI2(1,j/dd(2))*phif(i,j-1,k) + cI2(2,j/dd(2))*phif(i,j+1,k)
                    end do
                end do
            end do

        end if
        !===========================================================================================================
        if( dd(1) /= 1 ) then

            do k = 1, Nf(3)
                do j = 1, Nf(2)
                    !pgi$ unroll = n:8
                    do i = 1+1, Nf(1)-1, dd(1)
                        phif(i,j,k) = cI1(1,i/dd(1))*phif(i-1,j,k) + cI1(2,i/dd(1))*phif(i+1,j,k)
                    end do
                end do
            end do

        end if
        !===========================================================================================================


    end subroutine MG_interpolate



    !> \todo solve dirty hack
    !! \todo test 3d
    subroutine MG_interpolateV( &
        dimens,                 &
        dir,                    &
        Nc,                     &
        bLc,bUc,                &
        SSc,NNc,                &
        Nf,                     &
        bLf,bUf,                &
        SSf,NNf,                &
        BCL,BCU,                &
        cIV,                    &
        cI1,                    &
        cI2,                    &
        cI3,                    &
        phic,                   &
        phif ) bind (c,name='MG_interpolateV')

        implicit none

        integer(c_int), intent(in)     :: dimens

        integer(c_int), intent(in)     :: dir

        integer(c_int), intent(in)     :: Nc(1:3)

        integer(c_int), intent(in)     :: bLc(1:3)
        integer(c_int), intent(in)     :: bUc(1:3)

        integer(c_int), intent(in)     :: SSc(1:3)
        integer(c_int), intent(in)     :: NNc(1:3)

        integer(c_int), intent(in)     :: Nf(1:3)

        integer(c_int), intent(in)     :: bLf(1:3)
        integer(c_int), intent(in)     :: bUf(1:3)

        integer(c_int), intent(in)     :: SSf(1:3)
        integer(c_int), intent(in)     :: NNf(1:3)

        integer(c_int), intent(in)     :: BCL(1:3)
        integer(c_int), intent(in)     :: BCU(1:3)

        real(c_double),  intent(in)    :: cIV ( -1:2, 0:Nc(dir) )

        real(c_double),  intent(in)    :: cI1 ( 1:2, 1:Nc(1) )
        real(c_double),  intent(in)    :: cI2 ( 1:2, 1:Nc(2) )
        real(c_double),  intent(in)    :: cI3 ( 1:2, 1:Nc(3) )

        real(c_double),  intent(in)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

        real(c_double),  intent(out)  :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))



        integer                ::  i, ic
        integer                ::  j, jc
        integer                ::  k, kc

        integer                :: l

        integer                :: dd(1:3)

        integer                :: S(1:3)
        integer                :: N(1:3)


        do i = 1,3
            dd(i) = ( Nf(i)-1 )/( Nc(i)-1 )

            if( 0 < BCL(i) ) then
                S(i) = SSf(i)
            else
                S(i) = SSf(i)
            end if

            if( 0 < BCU(i) ) then
                N(i) = NNf(i)
            else
                N(i) = NNf(i)+1
            end if

        end do
        S(dir) = SSf(dir)

        !pgi$ unroll = n:8
        !        phif( SSf(1):NNf(1), SSf(2):NNf(2), SSf(3):NNf(3) ) = 0.

        if( 1==dir ) then
            !            do k = SSf(3)+1, NNf(3)-1, dd(3)
            do k = SSf(3), NNf(3)
                kc = k/dd(3)
                !                if( 2==dimens )
                kc = k
                do j = S(2), N(2), dd(2)
                    jc = ( j )/dd(2)
                    do i = S(1), N(1)
                        ic = ( i )/dd(1)
                        phif(i,j,k) = cIV(1,ic)*phic(ic,jc,kc)+cIV(2,ic+1)*phic(ic+1,jc,kc)
!                        phif(i,j,k) = 0.5*phic(ic,jc,kc)+0.5*phic(ic+1,jc,kc)
                    end do
                end do
            end do
            !            do k = 0, Nc(3)
            !                kc = (k-SSc(3))*dd(3) + SSf(3)
            !!                if( 2==dimens )
            !kc = k
            !                do jc = 1, Nc(2)
            !                    j = ( jc-1 )*dd(2)+1
            !                    do ic = 0, Nc(1)
            !                        i = ( ic-1 )*dd(1)+1
            !                        do l=-1,2
            !                            phif(i+l,j,k) = phif(i+l,j,k) + cIV(l,ic)*phic(ic,jc,kc)
            !                        end do
            !                    end do
            !                end do
            !            end do
            !            do k = SSc(3), NNNc(3)
            !                kc = (k-SSc(3))*dd(3) + SSf(3)
            !                if( 2==dimens ) kc = k
            !                do jc = SSc(2),NNNc(2)
            !                    j = ( jc-SSc(2) )*dd(2) + SSf(2) + 1
            !                    do ic = SSc(1),NNNc(1)
            !                        i = ( ic-SSc(1) )*dd(1) + SSf(1)
            !                        do l=-1,2
            !                            phif(i+l,j,k) = phif(i+l,j,k) + cIV(l,ic)*phic(ic,jc,kc)
            !                        end do
            !                    end do
            !                end do
            !            end do

            !            if( dd(3) /= 1 ) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (NNf(3) == 2??)
            !                do k = SSf(3), NNf(3)
            !                    kc = ( k-SSf(3) )/dd(3) + SSc(3)
            !                    do j = SSf(2), NNf(2), dd(2)
            !                        jc = ( j-SSf(2) )/dd(2) + SSc(2)
            !                        !pgi$ unroll = n:8
            !                        do i = SSf(1), NNf(1)
            !                            phif(i,j,k) = cI3(1,kc)*phif(i,j,k-1) + cI3(2,kc)*phif(i,j,k+1)
            !                        end do
            !                    end do
            !                end do

            !               ! dirty hack
            !               if( BCL(1) > 0 ) then
            !                   phif( SSf(1):NNf(1),SSf(2)::NNf(2),SSf(3) ) = phif( SSf(1):NNf(1),SSf(2):NNf(2),SSf(3)+1 )
            !               end if
            !               if( BCU(1) > 0 ) then
            !                   phif( SSf(1):NNf(1),SSf(2)::NNf(2),NNf(3) ) = phif( SSf(1):NNf(1),SSf(2):NNf(2),NNf(3)-1 )
            !               end if
            !           end if

            if( dd(2) /= 1 ) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (NNf(3) == 2??)

                do k = SSf(3), NNf(3)
                    do j = S(2)+1, N(2)-1, dd(2)
                        jc = j/dd(2)
                        !pgi$ unroll = n:8
                        do i = S(1), N(1)
                            phif(i,j,k) = cI2(1,jc)*phif(i,j-1,k) + cI2(2,jc)*phif(i,j+1,k)
                        end do
                    end do
                end do

                ! dirty hack
                if( BCL(2) > 0 ) then
!                    phif( SSf(1):NNf(1),SSf(2)+1,SSf(3):NNf(3) ) = phif( SSf(1):NNf(1),SSf(2)+2,SSf(3):NNf(3) )
                    phif( SSf(1):NNf(1),SSf(2),SSf(3):NNf(3) ) = phif( SSf(1):NNf(1),SSf(2)+1,SSf(3):NNf(3) )
                end if
                if( BCU(2) > 0 ) then
                    phif( SSf(1):NNf(1),NNf(2),SSf(3):NNf(3) ) = phif( SSf(1):NNf(1),NNf(2)-1,SSf(3):NNf(3) )
                end if

            end if


        else
        end if
    end subroutine MG_interpolateV


end module cmod_InterpolationOp
