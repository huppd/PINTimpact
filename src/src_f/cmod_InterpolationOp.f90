!> \brief module providing functions to initiliaze and apply RestrictionOp
module cmod_InterpolationOp

  use iso_c_binding

  implicit none

contains


  !> \todo fix that, for first entry we get the weights 0.4/0.6
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

      !cI(1,i) = 1.- Dx10/Dx12
      !cI(2,i) =     Dx10/Dx12
      cI(1,i) = 0.5
      cI(2,i) = 0.5
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

    integer(c_int), intent(in)    :: bLc, bUc

    integer(c_int), intent(in)    :: SSc
    integer(c_int), intent(in)    :: NNc

    integer(c_int), intent(in)    :: BC_L,BC_U

    integer(c_int), intent(in)    :: Nf

    integer(c_int), intent(in)    :: bLf, bUf

    integer(c_int), intent(in)    :: SSf
    !        integer(c_int), intent(in)    :: NNf


    real(c_double), intent(in)    :: xc( bLc:(Nc+bUc) )

    real(c_double), intent(in)    :: xf( bLf:(Nf+bUf) )

    real(c_double), intent(inout) :: cIV( 1:2, 0:Nf )

    integer(c_int)                ::  i, ic, dd
    real(c_double)                ::  Dx1a, Dx12


    !===========================================================================================================
    !=== Interpolation, linienweise, 1d ========================================================================
    !===========================================================================================================
    ! fine
    !    xf(i-1)        xf(i)         xf(i+1)       xf(i+1)      xf(i+2)
    !  ---->-------------->--------------->----------->----------->------------
    !             |--Dx1a-|
    !             |-----------Dx12--------------|
    !  ----------->----------------------------->----------------------------->----------
    !          xc(ic)                      xc(ic+1)                        xc(ic+1)
    ! coarse

    cIV = 0.

    dd = (Nf-1)/(Nc-1)

    do i = 0, Nf
      ic = ( i )/dd

      Dx1a = xc(ic)-xf(i   )
      Dx12 = xc(ic)-xc(ic+1)

      cIV(1,i) = 1-Dx1a/Dx12
      cIV(2,i) =   Dx1a/Dx12
    end do

    ! a little bit shaky, please verify this when IO is ready
    if (BC_L > 0) then
      cIV( 1,0) = 1.
      cIV( 2,0) = 0.
    end if

    ! maybe false
    if (BC_L == -2) then
      cIV(1:2,0) = 0.
    end if

    if (BC_U > 0) then
      cIV(1,Nf) = 0.
      cIV(2,Nf) = 1.
    end if
    if (BC_U == -2) then
      cIV(1:2,Nf) = 0.
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
  !! \note right now we have fixed the weights, so that gridscretching is not yet usable
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

    integer(c_int)               ::  i,ic
    integer(c_int)               ::  j
    integer(c_int)               ::  k

    integer(c_int)                :: dd(3)


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
    if( dd(3) /= 1 ) then ! (dimens == 2) <==> (dd(3) == 1) automatisch erfüllt!

      do k = 1+1, Nf(3)-1, dd(3)
        ic = k/dd(3)
        do j = 1, Nf(2), dd(2)
          !pgi$ unroll = n:8
          do i = 1, Nf(1), dd(1)
            phif(i,j,k) = cI3(1,ic)*phif(i,j,k-1) + cI3(2,ic)*phif(i,j,k+1)
          end do
        end do
      end do

    end if
    !===========================================================================================================
    if( dd(2) /= 1 ) then ! TEST!!! in 2D wird hier doppelte Arbeit geleistet! (NNf(3) == 2??)

      do k = 1, Nf(3)
        do j = 1+1, Nf(2)-1, dd(2)
          ic = j/dd(2)
          !pgi$ unroll = n:8
          do i = 1, Nf(1), dd(1)
            phif(i,j,k) = cI2(1,ic)*phif(i,j-1,k) + cI2(2,ic)*phif(i,j+1,k)
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
              ic = i/dd(1)
              phif(i,j,k) = cI1(1,ic)*phif(i-1,j,k) + cI1(2,ic)*phif(i+1,j,k)
            end do
          end do
        end do

      end if
      !===========================================================================================================

    end subroutine MG_interpolate



    !> \brief interpolating 
    !! \todo implement 3d
    !! \test make fallunterscheidung for dd(dir)==0 (simple copy then)
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

      real(c_double),  intent(in)    :: cIV ( 1:2, 0:Nf(dir) )

      real(c_double),  intent(in)    :: cI1 ( 1:2, 1:Nc(1) )
      real(c_double),  intent(in)    :: cI2 ( 1:2, 1:Nc(2) )
      real(c_double),  intent(in)    :: cI3 ( 1:2, 1:Nc(3) )

      real(c_double),  intent(in)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

      real(c_double),  intent(out)  :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

      integer(c_int)                ::  i, ic
      integer(c_int)                ::  j, jc
      integer(c_int)                ::  k, kc

      integer(c_int)                :: l

      integer(c_int)                :: dd(1:3)

      integer(c_int)                :: S(1:3)
      integer(c_int)                :: N(1:3)


      do i = 1,3
        dd(i) = ( Nf(i)-1 )/( Nc(i)-1 )
        !
        !            if( 0 < BCL(i) ) then
        !                S(i) = 0
        !            else
        !                S(i) = SSf(i)
        !            end if
        !
        !            if( 0 < BCU(i) ) then
        !                N(i) = NNf(i)+1
        !            else
        !                N(i) = NNf(i)+1
        !            end if
        !
      end do
      !        S(dir) = SSf(dir)


      if( 1==dir ) then

        if( dd(1) /= 1 ) then

          do k = SSf(3), Nf(3), dd(3)
            kc = ( k+1 )/dd(3)
            do j = SSf(2), Nf(2), dd(2)
              jc = ( j+1 )/dd(2) ! holy shit
              do i = SSf(1), Nf(1) ! zero for dirichlet
                ic = ( i )/dd(1)
                phif(i,j,k) = cIV(1,i)*phic(ic,jc,kc)+cIV(2,i)*phic(ic+1,jc,kc)
              end do
            end do
          end do

        else

          do k = SSf(3), Nf(3), dd(3)
            kc = ( k+1 )/dd(3)
            do j = SSf(2), Nf(2), dd(2)
              jc = ( j+1 )/dd(2) 
              do i = SSf(1), Nf(1)
                phif(i,j,k) = phic(i,jc,kc)
              end do
            end do
          end do

        end if

        if( dd(2) /= 1 ) then

          do k = 1, Nf(3), dd(3)
            do j = SSf(2)+1, Nf(2)-1, dd(2)
              jc = (  j+1  )/dd(2)
              !pgi$ unroll = n:8
              do i = SSf(1), Nf(1)
                phif(i,j,k) = 0.5*phif(i,j-1,k) + 0.5*phif(i,j+1,k)
              end do
            end do
          end do

        end if

        if( dd(3) /= 1 ) then

          do k = SSf(3)+1, Nf(3)-1, dd(3)
            do j = SSf(2), Nf(2)
              !pgi$ unroll = n:8
              do i = SSf(1), Nf(1)
                phif(i,j,k) = 0.5*phif(i,j,k-1) + 0.5*phif(i,j,k+1)
              end do
            end do
          end do

        end if

      end if

      if( 2==dir ) then

        if( dd(2) /= 1 ) then

          do k = SSf(3), Nf(3), dd(3)
            kc = ( k+1 )/dd(3)
            do j = SSf(2), Nf(2) ! zero for dirichlet
              jc = ( j )/dd(2)
              do i = SSf(1), Nf(1), dd(1)
                ic = ( i+1 )/dd(2) ! holy shit
                phif(i,j,k) = cIV(1,j)*phic(ic,jc,kc)+cIV(2,j)*phic(ic,jc+1,kc)
              end do
            end do
          end do

        else

          do k = SSf(3), Nf(3), dd(3)
            kc = ( k+1 )/dd(3)
            do j = SSf(2), Nf(2) 
              do i = SSf(1), Nf(1), dd(1)
                ic = ( i+1 )/dd(2) 
                phif(i,j,k) = phic(ic,j,kc)
              end do
            end do
          end do

        end if

        if( dd(1) /= 1 ) then

          do k = 1, Nf(3), dd(3)
            do j = SSf(2), Nf(2)
              !pgi$ unroll = n:8
              do i = SSf(1)+1, Nf(1)-1, dd(1)
                phif(i,j,k) = 0.5*phif(i-1,j,k) + 0.5*phif(i+1,j,k)
              end do
            end do
          end do

        end if

        if( dd(3) /= 1 ) then

          do k = SSf(3)+1, Nf(3)-1, dd(3)
            do j = SSf(2), Nf(2)
              !pgi$ unroll = n:8
              do i = SSf(1), Nf(1)
                phif(i,j,k) = 0.5*phif(i,j,k-1) + 0.5*phif(i,j,k+1)
              end do
            end do
          end do

        end if

      end if

      if( 3==dir ) then

        if( dd(3) /= 1 ) then

          do k = SSf(3), Nf(3)
            kc = ( k )/dd(3)
            do j = SSf(2), Nf(2), dd(2) ! zero for dirichlet
              jc = ( j+1 )/dd(2)
              do i = SSf(1), Nf(1), dd(1)
                ic = ( i+1 )/dd(2) ! holy shit
                phif(i,j,k) = cIV(1,k)*phic(ic,jc,kc)+cIV(2,k)*phic(ic,jc,kc+1)
              end do
            end do
          end do

        else

          do k = SSf(3), Nf(3)
            do j = SSf(2), Nf(2), dd(2)
              jc = ( j+1 )/dd(2)
              do i = SSf(1), Nf(1), dd(1)
                ic = ( i+1 )/dd(2)
                phif(i,j,k) = phic(ic,jc,k)
              end do
            end do
          end do

        end if

        if( dd(1) /= 1 ) then

          do k = SSf(3), Nf(3)
            do j = 1, Nf(2), dd(2)
              !pgi$ unroll = n:8
              do i = SSf(1)+1, Nf(1)-1, dd(1)
                phif(i,j,k) = 0.5*phif(i-1,j,k) + 0.5*phif(i+1,j,k)
              end do
            end do
          end do

        end if

        if( dd(2) /= 1 ) then

          do k = SSf(3), Nf(3)
            do j = SSf(2)+1, Nf(2)-1, dd(2)
              !pgi$ unroll = n:8
              do i = SSf(1), Nf(1)
                phif(i,j,k) = 0.5*phif(i,j-1,k) + 0.5*phif(i,j+1,k)
              end do
            end do
          end do

        end if

      end if

    end subroutine MG_interpolateV


end module cmod_InterpolationOp
