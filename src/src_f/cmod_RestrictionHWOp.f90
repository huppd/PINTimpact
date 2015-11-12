!> \brief module providing functions to initiliaze and apply RestrictionHWOp
module cmod_RestrictionHWOp

  use iso_c_binding

  use mpi

  implicit none

contains



  subroutine MG_getCRS( &
      iimax,            &
      BC_L,             &
      BC_U,             &
      dd,               &
      Nf,               &
      bL,               &
      bU,               &
      xf,               &
      cR ) bind( c, name='MG_getCRS' )

    implicit none


    integer(c_int), intent(in)  :: iimax
    integer(c_int), intent(in)  :: BC_L
    integer(c_int), intent(in)  :: BC_U

    integer(c_int), intent(in)  :: dd

    integer(c_int), intent(in)  :: Nf

    integer(c_int), intent(in)  :: bL, bU

    real(c_double), intent(in)  :: xf( bL:(Nf+bU) )

    real(c_double), intent(out) :: cR(-1:1,1:iimax)

    real(c_double)              :: h(1:2)
    integer(c_int)              :: i, ii

    cR = 0.

    !===========================================================================================================
    !=== Restriktion, linienweise, 1d ==========================================================================
    !===========================================================================================================

    do i = 1,iimax 
      if( 1==dd ) then
        cR(-1,i) = 0.
        cR( 0,i) = 1.
        cR( 1,i) = 0.
      else
        ii = dd*(i-1)+1
        h( 1 ) = xf(ii  ) - xf(ii-1)
        h( 2 ) = xf(ii+1) - xf(ii  )
        cR(-1,i) = h(2)/( h(1) + h(2) )/2. ! 1./4. for equi
        cR( 0,i) = 1./2.
        cR( 1,i) = h(1)/( h(1) + h(2) )/2. ! 1./4/ for equi
      end if
    end do

    if( BC_L > 0 ) then
      cR(-1,1) = 0.
      cR( 0,1) = 1.
      cR( 1,1) = 0.

      cR(-1,2) = 0. ! TEST!!! Sollte evtl. noch ergaenzt werden ...
      cR( 0,2) = 1.
      cR( 1,2) = 0.
    end if
    if (BC_L == -2) then
      cR( 1,1) = cR( 1,1) + cR(-1,1)
      cR(-1,1) = 0.
    end if

    if (BC_U > 0) then
      cR(-1,iimax) = 0.
      cR( 0,iimax) = 1.
      cR( 1,iimax) = 0.

      cR(-1,iimax-1) = 0.
      cR( 0,iimax-1) = 1.
      cR( 1,iimax-1) = 0.

    end if

    if (BC_U == -2) then
      cR(-1,iimax) = cR( 1,iimax) + cR(-1,iimax)
      cR( 1,iimax) = 0.
    end if

  end subroutine MG_getCRS



  !> \todo understand from where comes the weith 0.75 / 0.25
  !! \todo change to simple HW or use with a better coarsestrategy
  subroutine MG_getCRV(   &
      N,                  &
      bL, bU,             &
      SS,                 &
      iimax,              &
      BC_L, BC_U,         &
      xs,xv,              &
      cRV ) bind(c,name='MG_getCRV')


    implicit none

    integer(c_int), intent(in)    :: N

    integer(c_int), intent(in)    :: bL,bU

    integer(c_int), intent(in)    :: SS
    integer(c_int), intent(in)    :: iimax

    integer(c_int), intent(in)    :: BC_L,BC_U

    real(c_double), intent(in)    :: xs(bL:(N+bU))
    real(c_double), intent(in)    :: xv(bL:(N+bU))

    real(c_double), intent(inout) :: cRV(1:2,0:iimax)

    integer(c_int)                ::  i
    real(c_double)                ::  Dx12, Dx1a


    !===========================================================================================================
    !=== Restriktion, linienweise, 1d ==========================================================================
    !===========================================================================================================
    ! fine
    !    xs(i-1)        vv(i-1)         xs(i)       xv(i)
    !  ----o-------------->---------------o----------->-----
    !                     |------Dx1x-----|
    !                     |-----------Dx12------------|
    !                                     >
    !                                   xv(i-1)
    ! coarse
    cRV = 0.


    do i = 0, iimax
      Dx1a = xs(i)-xv(i-1)
      Dx12 = xv(i)-xv(i-1)

      cRV(1,i) = 1.- Dx1a/Dx12
      cRV(2,i) =     Dx1a/Dx12
    end do

    ! a little bit shaky, please verify this when IO is ready
    if (BC_L > 0) then
      cRV(1,0) = 0.
      cRV(2,0) = 1.
    end if
    if (BC_L == -2) then
      cRV(:,0) = 0.
    end if

    if (BC_U > 0) then
      cRV(1,iimax) = 1.
      cRV(2,iimax) = 0.
    end if
    if (BC_U == -2) then
      cRV(:,iimax) = 0.
    end if

  end subroutine MG_getCRV



  subroutine MG_RestrictCorners(  &
      Nc,                         &
      bLc,bUc,                    &
      BCL, BCU,                   &
      phif ) bind( c, name='MG_RestrictCorners' )

    implicit none

    integer(c_int), intent(in)     :: Nc(3)

    integer(c_int), intent(in)     :: bLc(3)
    integer(c_int), intent(in)     :: bUc(3)

    integer(c_int), intent(in)     :: BCL(3)
    integer(c_int), intent(in)     :: BCU(3)

    real(c_double), intent(inout)  :: phif (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

    !if (BC(1,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  1      ,  1      ,1:NN(3,g)) = (fine1(2        ,1      ,1:NN(3,g)) + fine1(1      ,2        ,1:NN(3,g)))/2.
    !if (BC(1,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  1      ,  NN(2,g),1:NN(3,g)) = (fine1(2        ,NN(2,g),1:NN(3,g)) + fine1(1      ,NN(2,g)-1,1:NN(3,g)))/2.
    !if (BC(2,1,g) > 0 .and. BC(1,2,g) > 0) fine1(  NN(1,g),  1      ,1:NN(3,g)) = (fine1(NN(1,g)-1,1      ,1:NN(3,g)) + fine1(NN(1,g),2        ,1:NN(3,g)))/2.
    !if (BC(2,1,g) > 0 .and. BC(2,2,g) > 0) fine1(  NN(1,g),  NN(2,g),1:NN(3,g)) = (fine1(NN(1,g)-1,NN(2,g),1:NN(3,g)) + fine1(NN(1,g),NN(2,g)-1,1:NN(3,g)))/2.

    if (BCL(1) > 0 .and. BCL(2) > 0) phif( 1,     1,     1:Nc(3) ) = ( phif(1+1,     1,     1:Nc(3)) + phif(1    , 1+1,     1:Nc(3)) )/2.
    if (BCL(1) > 0 .and. BCU(2) > 0) phif( 1,     Nc(2), 1:Nc(3) ) = ( phif(1+1,     Nc(2), 1:Nc(3)) + phif(1    , Nc(2)-1, 1:Nc(3)) )/2.
    if (BCU(1) > 0 .and. BCL(2) > 0) phif( Nc(1), 1,     1:Nc(3) ) = ( phif(Nc(1)-1, 1,     1:Nc(3)) + phif(Nc(1), 1+1,     1:Nc(3)) )/2.
    if (BCU(1) > 0 .and. BCU(2) > 0) phif( Nc(1), Nc(2), 1:Nc(3) ) = ( phif(Nc(1)-1, Nc(2), 1:Nc(3)) + phif(Nc(1), Nc(2)-1, 1:Nc(3)) )/2.

    !if (BC(1,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  1      ,1:NN(2,g),  1      ) = (fine1(2        ,1:NN(2,g),1      ) + fine1(1      ,1:NN(2,g),2        ))/2.
    !if (BC(1,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  1      ,1:NN(2,g),  NN(3,g)) = (fine1(2        ,1:NN(2,g),NN(3,g)) + fine1(1      ,1:NN(2,g),NN(3,g)-1))/2.
    !if (BC(2,1,g) > 0 .and. BC(1,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  1      ) = (fine1(NN(1,g)-1,1:NN(2,g),1      ) + fine1(NN(1,g),1:NN(2,g),2        ))/2.
    !if (BC(2,1,g) > 0 .and. BC(2,3,g) > 0) fine1(  NN(1,g),1:NN(2,g),  NN(3,g)) = (fine1(NN(1,g)-1,1:NN(2,g),NN(3,g)) + fine1(NN(1,g),1:NN(2,g),NN(3,g)-1))/2.

    if (BCL(1) > 0 .and. BCL(3) > 0) phif( 1,     1:Nc(2), 1 )     = ( phif(1+1,     1:Nc(2), 1    ) + phif(1,     1:Nc(2), 1+1)      )/2.
    if (BCL(1) > 0 .and. BCU(3) > 0) phif( 1,     1:Nc(2), Nc(3) ) = ( phif(1+1,     1:Nc(2), Nc(3)) + phif(1,     1:Nc(2), Nc(3)-1)  )/2.
    if (BCU(1) > 0 .and. BCL(3) > 0) phif( Nc(1), 1:Nc(2), 1 )     = ( phif(Nc(1)-1, 1:Nc(2), 1    ) + phif(Nc(1), 1:Nc(2), 1+1)      )/2.
    if (BCU(1) > 0 .and. BCU(3) > 0) phif( Nc(1), 1:Nc(2), Nc(3) ) = ( phif(Nc(1)-1, 1:Nc(2), Nc(3)) + phif(Nc(1), 1:Nc(2), Nc(3)-1)  )/2.

    !if (BC(1,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  1      ,  1      ) = (fine1(1:NN(1,g),2        ,1      ) + fine1(1:NN(1,g),1      ,2        ))/2.
    !if (BC(1,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  1      ,  NN(3,g)) = (fine1(1:NN(1,g),2        ,NN(3,g)) + fine1(1:NN(1,g),1      ,NN(3,g)-1))/2.
    !if (BC(2,2,g) > 0 .and. BC(1,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  1      ) = (fine1(1:NN(1,g),NN(2,g)-1,1      ) + fine1(1:NN(1,g),NN(2,g),2        ))/2.
    !if (BC(2,2,g) > 0 .and. BC(2,3,g) > 0) fine1(1:NN(1,g),  NN(2,g),  NN(3,g)) = (fine1(1:NN(1,g),NN(2,g)-1,NN(3,g)) + fine1(1:NN(1,g),NN(2,g),NN(3,g)-1))/2.

    if (BCL(2) > 0 .and. BCL(3) > 0) phif( 1:Nc(1), 1,     1     ) = ( phif(1:Nc(1), 1+1,     1    ) + phif(1:Nc(1), 1,     1+1    ) )/2.
    if (BCL(2) > 0 .and. BCU(3) > 0) phif( 1:Nc(1), 1,     Nc(3) ) = ( phif(1:Nc(1), 1+1,     Nc(3)) + phif(1:Nc(1), 1,     Nc(3)-1) )/2.
    if (BCU(2) > 0 .and. BCL(3) > 0) phif( 1:Nc(1), Nc(2), 1     ) = ( phif(1:Nc(1), Nc(2)-1, 1    ) + phif(1:Nc(1), Nc(2), 1+1    ) )/2.
    if (BCU(2) > 0 .and. BCU(3) > 0) phif( 1:Nc(1), Nc(2), Nc(3) ) = ( phif(1:Nc(1), Nc(2)-1, Nc(3)) + phif(1:Nc(1), Nc(2), Nc(3)-1) )/2.

  end subroutine MG_RestrictCorners



  subroutine MG_restrictHW( &
      dimens,             &
      Nf,                 &
      bLf,bUf,            &
      Nc,                 &
      bLc,bUc,            &
      iimax,              &
      dd,                 &
      cR1,cR2,cR3,        &
      phif,               &
      phic ) bind(c,name='MG_restrictHW')

    implicit none

    integer(c_int), intent(in)     :: dimens

    integer(c_int), intent(in)     :: Nf(3)

    integer(c_int), intent(in)     :: bLf(3)
    integer(c_int), intent(in)     :: bUf(3)


    integer(c_int), intent(in)     :: Nc(3)

    integer(c_int), intent(in)     :: bLc(3)
    integer(c_int), intent(in)     :: bUc(3)

    integer(c_int), intent(in)     :: iimax(3)

    integer(c_int), intent(in)     :: dd(3)

    real(c_double),  intent(in)    :: cR1 ( -1:1, 1:iimax(1) )
    real(c_double),  intent(in)    :: cR2 ( -1:1, 1:iimax(2) )
    real(c_double),  intent(in)    :: cR3 ( -1:1, 1:iimax(3) )

    real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

    real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


    integer(c_int)                ::  i, ii
    integer(c_int)                ::  j, jj
    integer(c_int)                ::  k, kk


    !----------------------------------------------------------------------------------------------------------!
    ! Anmerkungen: - für allgemeine di, dj, dk geeignet!                                                       !
    !              - überlappende Schicht in Blöcken wird (der Einfachheit halber) ebenfalls ausgetauscht, ist !
    !                aber im Prinzip redundant (genauer: phiC(S1R:N1R,S2R:N2R,S3R:N3R) = ...).               !
    !              - Motivation für diese kurze Routine ist die Möglichkeit, auch Varianten wie Full-Weighting !
    !                etc. ggf. einzubauen, ansonsten könnte sie auch eingespaart werden.                       !
    !              - Die Block-überlappenden Stirnflächen werden ebenfalls mitverarbeitet, aber eigentlich     !
    !                nicht gebraucht (erleichtert die Programmierung), so dass eine Initialisierung notwendig  !
    !                ist. Dies wiederum bedingt die INTENT(inout)-Deklaration.                                 !
    !----------------------------------------------------------------------------------------------------------!



    if (dimens == 3) then
            do kk = 1, iimax(3)
                k = dd(3)*(kk-1)+1
                do jj = 1, iimax(2) 
                    j = dd(2)*(jj-1)+1
                    do ii = 1, iimax(1)
                        i = dd(1)*(ii-1)+1
                        phic(ii,jj,kk) = ((cR1( 0,ii)+cR2( 0,jj)+cR3( 0,kk))*phif(i  ,j  ,k  ) +  &
                            &              cR1(-1,ii)                       *phif(i-1,j  ,k  ) +  &
                            &              cR1( 1,ii)                       *phif(i+1,j  ,k  ) +  &
                            &                         cR2(-1,jj)            *phif(i  ,j-1,k  ) +  &
                            &                         cR2( 1,jj)            *phif(i  ,j+1,k  ) +  &
                            &                                     cR3(-1,kk)*phif(i  ,j  ,k-1) +  &
                            &                                     cR3( 1,kk)*phif(i  ,j  ,k+1)) / 3.
                    end do
                end do
            end do
        else
            k  = 1
            kk = 1
            do jj = 1, iimax(2)
                j = dd(2)*(jj-1)+1
                do ii = 1, iimax(1)
                    i = dd(1)*(ii-1)+1
                    phic(ii,jj,kk) = ((cR1( 0,ii)+cR2(0,jj))*phif(i  ,j  ,k) +  &
                        &              cR1(-1,ii)           *phif(i-1,j  ,k) +  &
                        &              cR1( 1,ii)           *phif(i+1,j  ,k) +  &
                        &                         cR2(-1,jj)*phif(i  ,j-1,k) +  &
                        &                         cR2( 1,jj)*phif(i  ,j+1,k)) / 2.
                end do
            end do
        end if


    end subroutine MG_restrictHW



    !> \todo maybe use SS and NN with boundary instead should be more cleaner
    !! \todo make restrictHWion proper in 3d not just in 1d!!!
    subroutine MG_restrictHWV(    &
        dimens,                 &
        dir,                    &
        Nf,                     &
        bLf,bUf,                &
        SSf,NNf,                &
        Nc,                     &
        bLc,bUc,                &
        SSc,NNc,                &
        iimax,                  &
        dd,                     &
        cRV,                    &
        phif,                   &
        phic ) bind (c,name='MG_restrictHWV')

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

        integer(c_int), intent(in)     :: iimax(1:3)

        integer(c_int), intent(in)     :: dd(1:3)

        real(c_double),  intent(in)    :: cRV ( 1:2, 0:iimax(dir) )

        real(c_double),  intent(in)    :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

        real(c_double),  intent(out)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


        integer(c_int)               ::  i, ii
        integer(c_int)               ::  j, jj
        integer(c_int)               ::  k, kk



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



        !dd=1
        !do i=1,dimens
          !!if( 0/=NNc(i)-SSc(i) ) dd(i) = ( NNf(i)-SSf(i) )/( NNc(i)-SSc(i) )
          !if( 0/=Nc(i) ) dd(i) = ( Nf(i) )/( Nc(i)-1 )
        !end do
        !!write(*,*) dd


        !===========================================================================================================
        if( 1==dir ) then

          if( 2==dimens ) then

            kk = SSc(3)
            k = SSc(3)
            do jj = SSc(2), iimax(2)
                if( 1==dd(2) ) then
                  j = jj
                else
                  j = dd(2)*jj-1
                endif
              do ii = SSc(1), iimax(1)
                if( 1==dd(1) ) then
                  i = ii
                  phic(ii,jj,kk) = phif(i,j,k)
                else
                  i = dd(1)*ii-1
                  phic(ii,jj,kk) = cRV(1,ii)*phif(i,j,k) + cRV(2,ii)*phif(i+1,j,k)
                end if
              end do
            end do

          else

            do kk = SSc(3), iimax(3)
              if( 1==dd(3) ) then
                k = kk
              else
                k = dd(3)*kk-1
              end if
              do jj = SSc(2), iimax(2)
                if( 1==dd(2) ) then
                  j = jj
                else
                  j = dd(2)*jj-1
                end if
                do ii = SSc(1), iimax(1)
                  if( 1==dd(1) ) then
                    i = ii
                    phic(ii,jj,kk) = phif(i,j,k)
                  else
                    i = dd(1)*ii-1
                    phic(ii,jj,kk) = cRV(1,ii)*phif(i,j,k) + cRV(2,ii)*phif(i+1,j,k)
                  end if
                end do
              end do
            end do

          end if

        end if
        !===========================================================================================================
        if( 2==dir ) then

          if( 2==dimens ) then

            kk = SSc(3)
            k = SSc(3)
            do jj = SSc(2), iimax(2)
              if( 1==dd(2) ) then
                j = jj
              else
                j = dd(2)*jj-1
              end if
              do ii = SSc(1), iimax(1)

                    if( 1==dd(1) ) then
                      i = ii
                    else
                      i = dd(1)*ii-1
                    end if
                    if( 1==dd(2) ) then
                      phic(ii,jj,kk) = phif(i,j,k)
                    else
                      phic(ii,jj,kk) = cRV(1,jj)*phif(i,j,k) + cRV(2,jj)*phif(i,j+1,k)
                    end if

                end do
              end do

            else

              do kk = SSc(3), iimax(3)
                if( 1==dd(3) ) then
                  k = kk
                else
                  k = dd(3)*kk-1
                end if
                do jj = SSc(2), iimax(2)
                  if( 1==dd(2) ) then
                    j = jj
                  else
                    j = dd(2)*jj-1
                  end if
                  do ii = SSc(1), iimax(1)
                    if( 1==dd(1) ) then
                      i = ii
                    else
                      i = dd(1)*ii-1
                    end if
                    if( 1==dd(2) ) then
                      phic(ii,jj,kk) = phif(i,j,k)
                    else
                      phic(ii,jj,kk) = cRV(1,jj)*phif(i,j,k) + cRV(2,jj)*phif(i,j+1,k)
                    end if
                  end do
                end do
              end do

            end if

        end if
        !===========================================================================================================
        if( 3==dimens .and. dir==3 ) then

          do kk = SSc(3), iimax(3)
            if( 1==dd(3) ) then
              k = kk
            else
              k = dd(3)*kk-1
            end if
            do jj = SSc(2), iimax(2)
              if( 1==dd(2) ) then
                j = jj
              else
                j = dd(2)*jj-1
              end if
              do ii = SSc(1), iimax(1)
                if( 1==dd(1) ) then
                  i = ii
                else
                  i = dd(1)*ii-1
                end if
                if( 1==dd(3) ) then
                  phic(ii,jj,kk) = phif(i,j,k)
                else
                  phic(ii,jj,kk) = cRV(1,kk)*phif(i,j,k) + cRV(2,kk)*phif(i,j,k+1)
                end if
              end do
            end do
          end do

        end if
    !===========================================================================================================


    end subroutine MG_restrictHWV

    

    !> \todo ununderstand recevbuf size not coorect
    subroutine MG_RestrictGather( &
        Nc,                       &
        bLc,bUc,                  &
        iimax,                    &
        n_gather,                 &
        participate_yes,          &
        rankc2,                   &
        comm2,                    &
        recvR,                    &
        dispR,                    &
        sizsR,                    &
        offsR,                    &
        phic ) bind(c,name='MG_RestrictGather')

        implicit none


        integer(c_int), intent(in)     :: Nc(3)

        integer(c_int), intent(in)     :: bLc(3)
        integer(c_int), intent(in)     :: bUc(3)

        integer(c_int), intent(in)     :: iimax(3)

        integer(c_int), intent(in)     :: n_gather(3)

        logical(c_bool), intent(in)     :: participate_yes

        integer(c_int), intent(in)     :: rankc2

        integer(c_int), intent(in)     :: comm2

        integer(c_int), intent(in)     :: recvR(     1:n_gather(1)*n_gather(2)*n_gather(3) )
        integer(c_int), intent(in)     :: dispR(     1:n_gather(1)*n_gather(2)*n_gather(3) )
        integer(c_int), intent(in)     :: sizsR(1:3, 1:n_gather(1)*n_gather(2)*n_gather(3) )
        integer(c_int), intent(in)     :: offsR(1:3, 1:n_gather(1)*n_gather(2)*n_gather(3) )


        real(c_double),intent(inout)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


        integer(c_int)                ::  i, ii
        integer(c_int)                ::  j, jj
        integer(c_int)                ::  k, kk

        integer(c_int)                ::  merror
        real(c_double), allocatable   ::  sendbuf(:,:,:)
        integer(c_int)                ::  sizsg(1:3), offsg(1:3), dispg
        real(c_double)                ::  recvbuf( 1:( Nc(1)*Nc(2)*Nc(3)*3 ) )

        !if (n_gather(1)*n_gather(2)*n_gather(3) > 1) then
    
          allocate(sendbuf(1:iimax(1),1:iimax(2),1:iimax(3))) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "phic" verwenden!
    
          do kk = 1, iimax(3)
            do jj = 1, iimax(2)
              do ii = 1, iimax(1)
                sendbuf(ii,jj,kk) = phic(ii,jj,kk)
              end do
            end do
          end do
    
          call MPI_GATHERv(sendbuf,iimax(1)*iimax(2)*iimax(3),MPI_REAL8,recvbuf,recvR,dispR,MPI_REAL8,rankc2,comm2,merror)

          deallocate(sendbuf)
    
    
          if( participate_yes ) then
            do k = 1, n_gather(3)
              do j = 1, n_gather(2)
                do i = 1, n_gather(1)
          
                  sizsg(1:3) = sizsR(1:3,i+(j-1)*n_gather(1)+(k-1)*n_gather(1)*n_gather(2))
                  offsg(1:3) = offsR(1:3,i+(j-1)*n_gather(1)+(k-1)*n_gather(1)*n_gather(2))
                  dispg      = dispR(    i+(j-1)*n_gather(1)+(k-1)*n_gather(1)*n_gather(2))
          
                  do kk = 1, sizsg(3)
                    do jj = 1, sizsg(2)
                      do ii = 1, sizsg(1)
                        phic(ii+offsg(1),jj+offsg(2),kk+offsg(3)) = recvbuf(dispg+ii+(jj-1)*sizsg(1)+(kk-1)*sizsg(1)*sizsg(2))
                      end do
                    end do
                  end do
          
                end do
              end do
            end do
          end if
          
        !end if

  end subroutine MG_RestrictGather


end module cmod_RestrictionHWOp
