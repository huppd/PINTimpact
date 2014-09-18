!> \brief module providing functions to initiliaze and apply RestrictionOp
module cmod_InterpolationOp

    use iso_c_binding

contains

    subroutine MG_getCIS(   &
        N,                  &
        bL, bU,             &
        SS,                 &
        NN,                 &
        xs,                 &
        cI ) bind(c,name='MG_getCIS')

        implicit none

        integer(c_int), intent(in)  :: N
        integer(c_int), intent(in)  :: bL, bU

        integer(c_int), intent(in)  :: SS
        integer(c_int), intent(in)  :: NN

        real(c_double), intent(in)  :: xs(bL:(N+bU))

        real(c_double), intent(out) :: cI(1:2,SS:NN)

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
        do i = SS, NN
            Dx10 = xs(i  )-xs(i-1)
            Dx12 = xs(i+1)-xs(i-1)

            cI(1,i) = 1.- Dx10/Dx12
            cI(2,i) =     Dx10/Dx12
        end do
        !--------------------------------------------------------------------------------------------------------
    end subroutine MG_getCIS



    !> \todo understand from where comes the weith 0.75 / 0.25
    subroutine MG_getCIV(   &
        N,                  &
        bL, bU,             &
        SS,                 &
        NN,                 &
        BC_L, BC_U,         &
        xs,xv,              &
        cIV ) bind(c,name='MG_getCIV')


        implicit none

        integer(c_int),intent(in)    :: N


        integer(c_int),intent(in)    :: bL,bU

        integer(c_int),intent(in)    :: SS
        integer(c_int),intent(in)    :: NN

        integer(c_int),intent(in)    :: BC_L,BC_U

        real(c_double),intent(in)    :: xs(bL:(N+bU))
        real(c_double),intent(in)    :: xv(bL:(N+bU))

        real(c_double),intent(inout) :: cIV(1:2,SS:NN)

        integer(c_int)        ::  i
        real(c_double)        ::  Dx12, Dx1a


        !===========================================================================================================
        !=== Interpolation, linienweise, 1d ========================================================================
        !===========================================================================================================
        ! fine
        !    xs(i-1)        xv(i-1)         xs(i)       xv(i)
        !  ----o-------------->---------------o----------->-----
        !                     |------Dx1a-----|
        !                     |-----------Dx12------------|
        !                                     >
        !                                   xv(i-1)
        ! coarse
        cIV = 0.


        do i = SS, NN
            Dx1a = xs(i)-xv(i-1)
            Dx12 = xv(i)-xv(i-1)

            cIV(1,i) = 1.- Dx1a/Dx12
            cIV(2,i) =     Dx1a/Dx12
        end do

        ! a little bit shaky, please verify this when IO is ready
        if (BC_L > 0) then
            cIV(1,SS) = 1.
            cIV(2,SS) = 0.
            cIV(1,SS+1) = 0.
            cIV(2,SS+1) = 1.
        end if
        if (BC_L == -2) then
            cIV(:,SS) = 0.
        end if

        if (BC_U > 0) then
            cIV(1,NN) = 0.
            cIV(2,NN) = 1.
        end if
        if (BC_U == -2) then
            cIV(1:2,NN) = 0.
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
        SSc,NNc,            &
        BCL, BCU,           &
        Nf,                 &
        bLf,bUf,            &
        SSf,NNf,            &
        cI1,cI2,cI3,        &
        phic,               &
        phif ) bind(c,name='MG_interpolate')

        implicit none

        integer(c_int), intent(in)     :: dimens

        integer(c_int), intent(in)     :: Nc(3)

        integer(c_int), intent(in)     :: bLc(3)
        integer(c_int), intent(in)     :: bUc(3)

        integer(c_int), intent(in)     :: SSc(3)
        integer(c_int), intent(in)     :: NNc(3)

        integer(c_int), intent(in)     :: BCL(3)
        integer(c_int), intent(in)     :: BCU(3)

        integer(c_int), intent(in)     :: Nf(3)

        integer(c_int), intent(in)     :: bLf(3)
        integer(c_int), intent(in)     :: bUf(3)

        integer(c_int), intent(in)     :: SSf(3)
        integer(c_int), intent(in)     :: NNf(3)

        real(c_double),  intent(in)    :: cI1 ( 1:2, SSc(1):Nc(1) )
        real(c_double),  intent(in)    :: cI2 ( 1:2, SSc(2):Nc(2) )
        real(c_double),  intent(in)    :: cI3 ( 1:2, SSc(3):Nc(3) )

        real(c_double),  intent(inout) :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

        real(c_double),  intent(out)   :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))



        integer                ::  i!, ii!, di, imax, iimax
        integer                ::  j!, jj!, dj, jmax, jjmax
        integer                ::  k!, kk!, dk, kmax, kkmax

        integer                :: dd(3)


        dd = 1
        do i=1,dimens
            dd(i) = ( NNf(i)-SSf(i) )/( NNc(i)-SSc(i) )
        end do


        if (BCL(1) > 0 .and. BCL(2) > 0) phic( SSc(1), SSc(2), SSc(3):NNc(3) ) = ( phic(SSc(1)+1, SSc(2), SSc(3):NNc(3)) + phic(SSc(1), SSc(2)+1, SSc(3):NNc(3)) + phic(SSc(1)+1, SSc(2)+1, SSc(3):NNc(3)) )/3.
        if (BCL(1) > 0 .and. BCU(2) > 0) phic( SSc(1), NNc(2), SSc(3):NNc(3) ) = ( phic(SSc(1)+1, NNc(2), SSc(3):NNc(3)) + phic(SSc(1), NNc(2)-1, SSc(3):NNc(3)) + phic(SSc(1)+1, NNc(2)-1, SSc(3):NNc(3)) )/3.
        if (BCU(1) > 0 .and. BCL(2) > 0) phic( NNc(1), SSc(2), SSc(3):NNc(3) ) = ( phic(NNc(1)-1, SSc(2), SSc(3):NNc(3)) + phic(NNc(1), SSc(2)+1, SSc(3):NNc(3)) + phic(NNc(1)-1, SSc(2)+1, SSc(3):NNc(3)) )/3.
        if (BCU(1) > 0 .and. BCU(2) > 0) phic( NNc(1), NNc(2), SSc(3):NNc(3) ) = ( phic(NNc(1)-1, NNc(2), SSc(3):NNc(3)) + phic(NNc(1), NNc(2)-1, SSc(3):NNc(3)) + phic(NNc(1)-1, NNc(2)-1, SSc(3):NNc(3)) )/3.

        if (BCL(1) > 0 .and. BCL(3) > 0) phic( SSc(1), SSc(2):NNc(2), SSc(3) ) = ( phic(SSc(1)+1, SSc(2):NNc(2), SSc(3)) + phic(SSc(1), SSc(2):NNc(2), SSc(3)+1) + phic(SSc(1)+1, SSc(2):NNc(2), SSc(3)+1) )/3.
        if (BCL(1) > 0 .and. BCU(3) > 0) phic( SSc(1), SSc(2):NNc(2), NNc(3) ) = ( phic(SSc(1)+1, SSc(2):NNc(2), NNc(3)) + phic(SSc(1), SSc(2):NNc(2), NNc(3)-1) + phic(SSc(1)+1, SSc(2):NNc(2), NNc(3)-1) )/3.
        if (BCU(1) > 0 .and. BCL(3) > 0) phic( NNc(1), SSc(2):NNc(2), SSc(3) ) = ( phic(NNc(1)-1, SSc(2):NNc(2), SSc(3)) + phic(NNc(1), SSc(2):NNc(2), SSc(3)+1) + phic(NNc(1)-1, SSc(2):NNc(2), SSc(3)+1) )/3.
        if (BCU(1) > 0 .and. BCU(3) > 0) phic( NNc(1), SSc(2):NNc(2), NNc(3) ) = ( phic(NNc(1)-1, SSc(2):NNc(2), NNc(3)) + phic(NNc(1), SSc(2):NNc(2), NNc(3)-1) + phic(NNc(1)-1, SSc(2):NNc(2), NNc(3)-1) )/3.

        if (BCL(2) > 0 .and. BCL(3) > 0) phic( SSc(1):NNc(1), SSc(2), SSc(3) ) = ( phic(SSc(1):NNc(1), SSc(2)+1, SSc(3)) + phic(SSc(1):NNc(1), SSc(2), SSc(3)+1) + phic(SSc(1):NNc(1), SSc(2)+1, SSc(3)+1) )/3.
        if (BCL(2) > 0 .and. BCU(3) > 0) phic( SSc(1):NNc(1), SSc(2), NNc(3) ) = ( phic(SSc(1):NNc(1), SSc(2)+1, NNc(3)) + phic(SSc(1):NNc(1), SSc(2), NNc(3)-1) + phic(SSc(1):NNc(1), SSc(2)+1, NNc(3)-1) )/3.
        if (BCU(2) > 0 .and. BCL(3) > 0) phic( SSc(1):NNc(1), NNc(2), SSc(3) ) = ( phic(SSc(1):NNc(1), NNc(2)-1, SSc(3)) + phic(SSc(1):NNc(1), NNc(2), SSc(3)+1) + phic(SSc(1):NNc(1), NNc(2)-1, SSc(3)+1) )/3.
        if (BCU(2) > 0 .and. BCU(3) > 0) phic( SSc(1):NNc(1), NNc(2), NNc(3) ) = ( phic(SSc(1):NNc(1), NNc(2)-1, NNc(3)) + phic(SSc(1):NNc(1), NNc(2), NNc(3)-1) + phic(SSc(1):NNc(1), NNc(2)-1, NNc(3)-1) )/3.



        !***********************************************************************************************************
        ! huppd: to be save?
!                phif = 0.
        !

        !pgi$ unroll = n:8
!        phif( SSf(1):NNf(1):dd(1), SSf(2):NNf(2):dd(2), SSf(3):NNf(3):dd(3) )   =   &
!        phic( SSc(1):NNc(1),       SSc(2):NNc(2),       SSc(3):NNc(3) )
        phif( SSf(1):Nf(1):dd(1), SSf(2):Nf(2):dd(2), SSf(3):Nf(3):dd(3) )   =   &
        phic( SSc(1):Nc(1),       SSc(2):Nc(2),       SSc(3):Nc(3) )

        !***********************************************************************************************************

!dd(2)=1
!dd = 1

        !===========================================================================================================
        if( dd(3) /= 1 ) then ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!

            do k = SSf(3)+1, Nf(3)-1, dd(3)
                do j = SSf(2), Nf(2), dd(2)
                    !pgi$ unroll = n:8
                    do i = SSf(1), Nf(1), dd(1)
                        phif(i,j,k) = cI3(1,k/dd(3))*phif(i,j,k-1) + cI3(2,k/dd(3))*phif(i,j,k+1)
                    end do
                end do
            end do

        end if
        !===========================================================================================================
        if( dd(2) /= 1 ) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (NNf(3) == 2??)

            do k = SSf(3), Nf(3)
                do j = SSf(2)+1, Nf(2)-1, dd(2)
                    !pgi$ unroll = n:8
                    do i = SSf(1), Nf(1), dd(1)
                        phif(i,j,k) = cI2(1,j/dd(2))*phif(i,j-1,k) + cI2(2,j/dd(2))*phif(i,j+1,k)
                    end do
                end do
            end do

        end if
        !===========================================================================================================
        if( dd(1) /= 1 ) then

            do k = SSf(3), Nf(3)
                do j = SSf(2), Nf(2)
                    !pgi$ unroll = n:8
                    do i = SSf(1)+1, Nf(1)-1, dd(1)
                        phif(i,j,k) = cI1(1,i/dd(1))*phif(i-1,j,k) + cI1(2,i/dd(1))*phif(i+1,j,k)
                    end do
                end do
            end do

        end if
        !===========================================================================================================


    end subroutine MG_interpolate



    !> \todo maybe use SS and NN with boundary instead should be more cleaner
    subroutine MG_interpolateV(    &
        dimens,                 &
        dir,                    &
        Nf,                     &
        bLf,bUf,                &
        SSf,NNf,                &
        Nc,                     &
        bLc,bUc,                &
        SSc,NNc,                &
        cIV,                    &
        phif,                   &
        phic ) bind (c,name='MG_interpolateV')

        implicit none

        integer(c_int), intent(in)     :: dimens

        integer(c_int), intent(in)     :: dir

        integer(c_int), intent(in)     :: Nf(1:3)

        integer(c_int), intent(in)     :: bLf(1:3)
        integer(c_int), intent(in)     :: bUf(1:3)

        integer(c_int), intent(in)     :: SSf(1:3)
        integer(c_int), intent(in)     :: NNf(1:3)

        integer(c_int), intent(in)     :: Nc(1:3)

        integer(c_int), intent(in)     :: bLc(1:3)
        integer(c_int), intent(in)     :: bUc(1:3)

        integer(c_int), intent(in)     :: SSc(1:3)
        integer(c_int), intent(in)     :: NNc(1:3)

        real(c_double),  intent(in)    :: cIV ( 1:2, SSc(dir):NNc(dir) )

        real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

        real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


        integer                ::  i, ii
        integer                ::  j, jj
        integer                ::  k, kk

        integer                :: dd(1:3)


                !----------------------------------------------------------------------------------------------------------!
                ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
                !              - Da nur in Richtung der jeweiligen Geschwindigkeitskomponente gemittelt wird, muss nicht   !
                !                die spezialisierte Helmholtz-Variante aufgerufen werden.                                  !
                !              - Austauschrichtung ist invers zu ex1, ex2, ex3. Bei mehreren Blöcken wird auch der jeweils !
                !                redundante "überlappende" Punkt aufgrund der zu grossen Intervallgrenzen (1:iimax) zwar   !
                !                berechnet, aber aufgrund des Einweg-Austauschs falsch berechnet! Dieses Vorgehen wurde    !
                !                bislang aus übersichtsgründen vorgezogen, macht aber eine Initialisierung notwendig.      !
                !                Dabei werden Intervalle der Form 0:imax anstelle von 1:imax bearbeitet, da hier nur die   !
                !                das feinste Geschwindigkeitsgitter behandelt wird!                                        !
                !              - INTENT(inout) ist bei den feinen Gittern notwendig, da Ghost-Werte ausgetauscht werden    !
                !                müssen.                                                                                   !
                !              - Zuviele Daten werden ausgetauscht; eigentlich müsste in der Grenzfläche nur jeder 4.      !
                !                Punkt behandelt werden (4x zuviel!). Leider etwas unschön, könnte aber durch eine         !
                !                spezialisierte Austauschroutine behandelt werden, da das übergeben von Feldern mit        !
                !                Intervallen von b1L:(iimax+b1U) nur sehr schlecht funktionieren würde (d.h. mit Um-       !
                !                kopieren).                                                                                !
                !----------------------------------------------------------------------------------------------------------!

                ! TEST!!! Test schreiben, um n_gather(:,2) .GT. 1 hier zu vermeiden! Gleiches gilt natürlich für die Interpolation.



        dd=1
        do i=1,dimens
            dd(i) = ( NNf(i)-SSf(i) )/( NNc(i)-SSc(i) )
        end do


        !===========================================================================================================
        if( 1==dir ) then

            if( 2==dimens ) then
                kk = SSc(3)
                k = SSc(3)
                do jj = SSc(2), NNc(2)
                    j = dd(2)*jj-1
                    do ii = SSc(1), NNc(1)
                        i = dd(1)*ii-1
                        phic(ii,jj,kk) = cIV(1,ii)*phif(i-1,j,k) + cIV(2,ii)*phif(i,j,k)
                    end do
                end do
            else
                do kk = SSc(3), NNc(3)
                    k = dd(3)*kk-1
                    do jj = SSc(2), NNc(2)
                        j = dd(2)*jj-1
                        do ii = SSc(1), NNc(1)
                            i = dd(1)*ii-1
                            phic(ii,jj,kk) = cIV(1,ii)*phif(i-1,j,k) + cIV(2,ii)*phif(i,j,k)
                        end do
                    end do
                end do
            end if

        end if
        !===========================================================================================================
        if( dir==2 ) then

            if( 2==dimens ) then
                kk = SSc(3)
                k = SSc(3)
                do jj = SSc(2), NNc(2)
                    j = dd(2)*jj-1
                    do ii = SSc(1), NNc(1)
                        i = dd(1)*ii-1
                        phic(ii,jj,kk) = cIV(1,jj)*phif(i,j-1,k) + cIV(2,jj)*phif(i,j,k)
                    end do
                end do
            else
                do kk = SSc(3), NNc(3)
                    k = dd(3)*kk-1
                    do jj = SSc(2), NNc(2)
                        j = dd(2)*jj-1
                        do ii = SSc(1), NNc(1)
                            i = dd(1)*ii-1
                            phic(ii,jj,kk) = cIV(1,jj)*phif(i,j-1,k) + cIV(2,jj)*phif(i,j,k)
                        end do
                    end do
                end do
            end if

        end if
        !===========================================================================================================
        if( 3==dimens .and. dir==3 ) then

            do kk = SSc(3), NNc(3)
                k = dd(3)*kk-1
                do jj = SSc(2), NNc(2)
                    j = dd(2)*jj-1
                    do ii = SSc(1), NNc(1)
                        i = dd(1)*ii-1
                        phic(ii,jj,kk) = cIV(1,kk)*phif(i,j,k-1) + cIV(2,kk)*phif(i,j,k)
                    end do
                end do
            end do

        end if
    !===========================================================================================================
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig!                                                      !
  !              - dd(i) /= 1 <===> N(i) /= 1                                                                     !
  !              - Da imax=N1 usw., könnten auf dem feinen Gitter auch problemlos die bekannten engeren      !
  !                Intervallgrenzen S11:N11 usw. benutzt werden. Wurde bislang aus übersichtsgründen nicht   !
  !                vollzogen.                                                                                !
  !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
  !                nicht gebraucht (erleichtert die Programmierung), so dass eigentlich eine Initialisierung !
  !                notwendig wäre. Dies wird jedoch zuvor schon in der korrespondierenden Restriktions-      !
  !                Routine erledigt, so dass dies hier nicht mehr notwendig ist.                             !
  !              - Es wird sequentiell über alle Raumrichtungen interpoliert, um keinen individuellen        !
  !                Interpolationsstencil für jeden Punkt im Raum speichern zu müssen.                        !
  !              - Durch das sequentielle Interpolieren kann der Interpolationsstencil klein und damit der   !
  !                Gesamtaufwand minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil (!)). Nachteilig    !
  !                ist dabei das zusätzliche Arbeitsfeld auf dem feineren Gitterniveau (wird der Multigrid-  !
  !                Routine entliehen).                                               .                       !
  !              - Interpolationskoeffizienten werden auf dem jeweils feineren Gitter gespeichert, um nicht  !
  !                auf die entsprechenden Indizes des gröberen Gitters umrechnen zu müssen.                  !
  !----------------------------------------------------------------------------------------------------------!


  iimax = Nc(1) ! TEST!!! iimax, etc. substituieren ...
  jjmax = Nc(2)
  kkmax = Nc(3)

  imax  = Nf(1)
  jmax  = Nf(2)
  kmax  = Nf(3)

  dd(1) = (imax-1)/(iimax-1)
  dd(2) = (jmax-1)/(jjmax-1)
  dd(3) = (kmax-1)/(kkmax-1)


  ! ACHTUNG!!! Verwendung von parametrischen Strides di,dj,dk in Schleifen verhindert die Vektorisierung bzw.
  !            das Prefetching!

  ! Anmerkung: .FALSE. ist ok!
!  call exchange_relax(g+1,ex1,ex2,ex3,0,.false.,coarse)


!pgi$ unroll = n:8
  work(1:Nf(1):dd(1),1:nf(2):dd(2),1:Nf(3):dd(3)) = coarse(1:Nc(1),1:Nc(2),1:Nc(3))


  !===========================================================================================================
  if (direction == 1) then

     if (dd(3) /= 1) then ! (dimens == 2) <==> (dd(3) == 1) automatisch erfüllt!
        do k = 2, Nf(3)-1, dd(3)
           do j = 1, Nf(2), dd(2)
!pgi$ unroll = n:8
              do i = 1, Nf(1), dd(1)
                 work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
              end do
           end do
        end do
     end if

     !--------------------------------------------------------------------------------------------------------
     if (dj /= 1) then
        do k = 1, kmax
           do j = 2, jmax-1, dj
!pgi$ unroll = n:8
              do i = 1, imax, di
                 work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------

     if (di /= 1) then
        do k = 1, kmax
           do j = 1, jmax
!pgi$ unroll = n:8
              do i = 1, imax-2, di
                 fine(i  ,j,k) = fine(i  ,j,k) + cIH1(1,i  )*work(i,j,k) + cIH1(2,i  )*work(i+2,j,k)
                 fine(i+1,j,k) = fine(i+1,j,k) + cIH1(1,i+1)*work(i,j,k) + cIH1(2,i+1)*work(i+2,j,k)
              end do
           end do
        end do

        !--- Extrapolation am Rand ---
        if (BC_1L > 0) then
           do k = 1, kmax
!pgi$ unroll = n:8
              do j = 1, jmax
                 fine(0   ,j,k) = fine(0   ,j,k) + cIH1(1,0   )*work(1      ,j,k) + cIH1(2,0   )*work(1+di,j,k)
              end do
           end do
        end if
        if (BC_1U > 0) then
           do k = 1, kmax
!pgi$ unroll = n:8
              do j = 1, jmax
                 fine(imax,j,k) = fine(imax,j,k) + cIH1(1,imax)*work(imax-di,j,k) + cIH1(2,imax)*work(imax,j,k)
              end do
           end do
        end if
     else
        if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        call MPI_FINALIZE(merror)
        stop
     end if

  end if
  !===========================================================================================================
  if (direction == 2) then

     if (dk /= 1) then ! (dimens == 2) <==> (dk == 1) automatisch erfüllt!
        do k = 2, kmax-1, dk
           do j = 1, jmax, dj
!pgi$ unroll = n:8
              do i = 1, imax, di
                 work(i,j,k) = cI3(1,k,g)*work(i,j,k-1) + cI3(2,k,g)*work(i,j,k+1)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (di /= 1) then
        do k = 1, kmax
           do j = 1, jmax, dj
!pgi$ unroll = n:8
              do i = 2, imax-1, di
                 work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (dj /= 1) then
        do k = 1, kmax
           do i = 1, imax
!pgi$ unroll = n:8
              do j = 1, jmax-2, dj
                 fine(i,j  ,k) = fine(i,j  ,k) + cIH2(1,j  )*work(i,j,k) + cIH2(2,j  )*work(i,j+2,k)
                 fine(i,j+1,k) = fine(i,j+1,k) + cIH2(1,j+1)*work(i,j,k) + cIH2(2,j+1)*work(i,j+2,k)
              end do
           end do
        end do

        !--- Extrapolation am Rand ---
        if (BC_2L > 0) then
           do k = 1, kmax
!pgi$ unroll = n:8
              do i = 1, imax
                 fine(i,0   ,k) = fine(i,0   ,k) + cIH2(1,0   )*work(i,1      ,k) + cIH2(2,0   )*work(i,1+dj,k)
              end do
           end do
        end if
        if (BC_2U > 0) then
           do k = 1, kmax
!pgi$ unroll = n:8
              do i = 1, imax
                 fine(i,jmax,k) = fine(i,jmax,k) + cIH2(1,jmax)*work(i,jmax-dj,k) + cIH2(2,jmax)*work(i,jmax,k)
              end do
           end do
        end if
     else
        if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        call MPI_FINALIZE(merror)
        stop
     end if

  end if
  !===========================================================================================================
  if (direction == 3) then ! (dimens == 2) <==> (direction /= 3) automatisch erfüllt!

     if (dj /= 1) then
        do k = 1, kmax, dk
           do j = 2, jmax-1, dj
!pgi$ unroll = n:8
              do i = 1, imax, di
                 work(i,j,k) = cI2(1,j,g)*work(i,j-1,k) + cI2(2,j,g)*work(i,j+1,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (di /= 1) then
        do k = 1, kmax, dk
           do j = 1, jmax
!pgi$ unroll = n:8
              do i = 2, imax-1, di
                 work(i,j,k) = cI1(1,i,g)*work(i-1,j,k) + cI1(2,i,g)*work(i+1,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (dk /= 1) then
        do j = 1, jmax
           do i = 1, imax
!pgi$ unroll = n:8
              do k = 1, kmax-2, dk
                 fine(i,j,k  ) = fine(i,j,k  ) + cIH3(1,k  )*work(i,j,k) + cIH3(2,k  )*work(i,j,k+2)
                 fine(i,j,k+1) = fine(i,j,k+1) + cIH3(1,k+1)*work(i,j,k) + cIH3(2,k+1)*work(i,j,k+2)
              end do
           end do
        end do

        !--- Extrapolation am Rand ---
        if (BC_3L > 0) then
           do j = 1, jmax
!pgi$ unroll = n:8
              do i = 1, imax
                 fine(i,j,0   ) = fine(i,j,0   ) + cIH3(1,0   )*work(i,j,1      ) + cIH3(2,0   )*work(i,j,1+dk)
              end do
           end do
        end if
        if (BC_3U > 0) then
           do j = 1, jmax
!pgi$ unroll = n:8
              do i = 1, imax
                 fine(i,j,kmax) = fine(i,j,kmax) + cIH3(1,kmax)*work(i,j,kmax-dk) + cIH3(2,kmax)*work(i,j,kmax)
              end do
           end do
        end if
     else
        if (rank == 0) write(*,*) 'ERROR! This choice of sub-grids for the Helmholtz-problem is not intended!'
        call MPI_FINALIZE(merror)
        stop
     end if

  end if
  !===========================================================================================================

    end subroutine MG_interpolateV


end module cmod_InterpolationOp
