!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing routine, that initializes the coordinates, which are stored in mod_vars
!! could be easily transfered to c++
!!  by using get_coords
module cmod_geometry

    use iso_c_binding
  
contains

    !> \brief sets xx and dx so that one gets a nice grid streching
    !!
    !! \note if ii0U==ii0L equidastant grid is used
    !! \param[in] Lmax Li
    !! \param[in] iimax Mi
    !! \param[in] ii index
    !! \param[out] xx
    !! \param[out] dx
    subroutine coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
        ! (sample subroutine)
  
        implicit none
  
        real   , intent(in)    ::  Lmax
        real   , intent(in)    ::  iimax
  
        real   , intent(in)    ::  ii0L
        real   , intent(in)    ::  ii0U
  
        real   , intent(in)    ::  ii
        real   , intent(out)   ::  xx
        real   , intent(out)   ::  dx
  
        real                   ::  yy, cmin, cmax
  
  
        if (ii0U == ii0L) then
            ! equidistant grid:
            xx = (ii-1.)*Lmax/(iimax-1.)
            dx =         Lmax/(iimax-1.)
        else
            cmax =  TANH(ii0U)
            cmin = -TANH(ii0L)
     
            ! coordinate transformation (index i=1. is the origin):
            ! y = (i-1.)/(imax-1.)*(cmax-cmin) + cmin
            yy = (ii-1.)/(iimax-1.)*(cmax-cmin) + cmin
     
            ! mapping funktion f(yy)
            ! x = L * (f(y)-f(cmin)) / (f(cmax)-f(cmin))
            xx = Lmax * (atanh(yy)-atanh(cmin)) / (ii0U+ii0L)
     
            ! dx/di = L / (f(cmax)-f(cmin)) * dy/di                 * df(y(i))/dy
            !       = L / (f(cmax)-f(cmin)) * (cmax-cmin)/(imax-1.) * df(y(i))/dy
            dx = Lmax / (atanh(cmax)-atanh(cmin)) * (cmax-cmin)/(iimax-1.) * 1./(1.-yy**2)
        end if
  
  
    end subroutine coord_tanh



    !> \brief here can the user sets the coordinates according to his favorite stretching
    !! \note: - local processor-block coordinates and grid spacings are automatically derived from global grid
    !!        - dy3p = dy3w = 1. for 2D (may simplify programming)
    !!        - ensure that for all i
    !!             y1p(i) < y1p(i+1)
    !!             y1u(i) < y1u(i+1)
    !!             y1p(i) < y1u(i) < y1p(i+1)
    !!                etc.
    !!        - code is tested only for
    !!             y1p(1 ) = 0.
    !!             y1p(M1) = L1
    !!                etc.
    subroutine PI_getGlobalCoordinates(  &
        dimens,             &
        L, M,               &
        y_origin,           &
        y1p,y2p,y3p,        &
        y1u,y2v,y3w,        &
        dy1p,dy2p,dy3p,     &
        dy1u,dy2v,dy3w      &
        ) bind(c,name='PI_getGlobalCoordinates')


        implicit none

        integer(c_int),intent(in)   :: dimens

        real(c_double),intent(in)   :: L(3)
        integer(c_int),intent(in)   :: M(3)

        real(c_double),intent(in)   :: y_origin(3)

        real(c_double),intent(inout)  :: y1p(1:M(1))
        real(c_double),intent(inout)  :: y2p(1:M(2))
        real(c_double),intent(inout)  :: y3p(1:M(3))

        real(c_double),intent(inout)  :: y1u(0:M(1))
        real(c_double),intent(inout)  :: y2v(0:M(2))
        real(c_double),intent(inout)  :: y3w(0:M(3))

        real(c_double),intent(inout)  :: dy1p(1:M(1))
        real(c_double),intent(inout)  :: dy2p(1:M(2))
        real(c_double),intent(inout)  :: dy3p(1:M(3))

        real(c_double),intent(inout)  :: dy1u(0:M(1))
        real(c_double),intent(inout)  :: dy2v(0:M(2))
        real(c_double),intent(inout)  :: dy3w(0:M(3))

        integer     ::  i



        !--- specify global coordinates and grid spacings ---
        ! example:
        !
        do i = 1, M(1)
            call coord_tanh(L(1),REAL(M(1)),0.,0.,REAL(i)    ,y1p(i),dy1p(i))
        end do
        do i = 0, M(1)
            call coord_tanh(L(1),REAL(M(1)),0.,0.,REAL(i)+0.5,y1u(i),dy1u(i))
        end do
        y1p = y1p - y_origin(1)
        y1u = y1u - y_origin(1)

        do i = 1, M(2)
            !     CALL coord_tanh(L(2),REAL(M(2)),0.,1.,REAL(i)    ,y2p(i),dy2p(i))
            call coord_tanh(L(2),REAL(M(2)),0.,0.,REAL(i)    ,y2p(i),dy2p(i))
        end do
        do i = 0, M(2)
            !     CALL coord_tanh(L(2),REAL(M(2)),0.,1.,REAL(i)+0.5,y2v(j),dy2v(j))
            call coord_tanh(L(2),REAL(M(2)),0.,0.,REAL(i)+0.5,y2v(i),dy2v(i))
        end do
        y2p = y2p - y_origin(2)
        y2v = y2v - y_origin(2)


        if( 3==dimens ) then
            do i = 1, M(3)
                call coord_tanh(L(3),REAL(M(3)),0.,0.,REAL(i)    ,y3p(i),dy3p(i))
            end do
            do i = 0, M(3)
                call coord_tanh(L(3),REAL(M(3)),0.,0.,REAL(i)+0.5,y3w(i),dy3w(i))
            end do
            y3p = y3p - y_origin(3)
            y3w = y3w - y_origin(3)
        end if


    end subroutine PI_getGlobalCoordinates


    !pgi$g unroll = n:8
    !!pgi$r unroll = n:8
    !!pgi$l unroll = n:8
    subroutine PI_getLocalCoordinates(  &
        M,                  &
        N,                  &
        bL,bU,              &
        x1p,x2p,x3p,                    &
        x1u,x2v,x3w,                    &
        dx1p,dx2p,dx3p,                 &
        dx1u,dx2v,dx3w                 &
        ) bind(c,name='PI_getLocalCoordinates')

        implicit none

        integer               ::  g
        integer               ::  Nmax, Mmax

        integer               ::  i, ii, di, iiShift
        integer               ::  j, jj, dj, jjShift
        integer               ::  k, kk, dk, kkShift

        real                  ::  y1pR(b1L:(M1+b1U)), y1uR(b1L:(M1+b1U))
        real                  ::  y2pR(b2L:(M2+b2U)), y2vR(b2L:(M2+b2U))
        real                  ::  y3pR(b3L:(M3+b3U)), y3wR(b3L:(M3+b3U))

        character(len=1)      ::  part, grid

        integer(c_int), intent(in)     :: M(3)

        integer(c_int), intent(in)     :: N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)


        real(c_double),intent(inout) ::  x1p(bL(1):(N(1)+bU(1)))
        real(c_double),intent(inout) ::  x2p(bL(2):(N(2)+bU(2)))
        real(c_double),intent(inout) ::  x3p(bL(3):(N(3)+bU(3)))

        real(c_double),intent(inout) ::  x1u(bL(1):(N(1)+bU(1)))
        real(c_double),intent(inout) ::  x2v(bL(2):(N(2)+bU(2)))
        real(c_double),intent(inout) ::  x3w(bL(3):(N(3)+bU(3)))


        x1p   = 0.
        x2p   = 0.
        x3p   = 0.

        x1u   = 0.
        x2v   = 0.
        x3w   = 0.

        dx1p  = 0.
        dx2p  = 0.
        dx3p  = 0.

        dx1u  = 0.
        dx2v  = 0.
        dx3w  = 0.

                !  Smagorinski-model
                !        dx1pS = 0.
                !        dx2pS = 0.
                !        dx3pS = 0.
                !
                !        dx1uS = 0.
                !        dx2vS = 0.
                !        dx3wS = 0.



        g=1

        Nmax = N(1)
        Mmax = M(1)

        di   = (M1-1) / (Mmax-1)

        y1pR = 0.
        y1uR = 0.

            y1pR(1:M1) = y1p(1:M1)
            y1uR(0:M1) = y1u(0:M1)

            do i = 1, Mmax
                y1pR(i) = y1p(1 + di*i - di  )
            end do
            do i = 1, Mmax-1
                y1uR(i) = y1p(1 + di*i - di/2)
            end do
        end if


        !--- periodische RB -------------------------------------------------------------------------------------
        if (BC_1L_global == -1) then
            y1pR(b1L:0) = y1pR((Mmax-1+b1L):(Mmax-1)) - L1
            y1uR(b1L:0) = y1uR((Mmax-1+b1L):(Mmax-1)) - L1

            y1pR(Mmax:(Mmax+b1U)) = y1pR(1:(1+b1U)) + L1
            y1uR(Mmax:(Mmax+b1U)) = y1uR(1:(1+b1U)) + L1
        end if


        !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
        if (BC_1L_global == -2 .or. BC_1L_global > 0) then
            if (BC_1L_global > 0 .and. g == 1) then
                y1pR(b1L: 0) = 2.*y1pR(1) - y1pR((2-b1L):2:-1)
                y1uR(b1L:-1) = 2.*y1pR(1) - y1uR((1-b1L):2:-1)
            else
                y1pR(b1L: 0) = 2.*y1pR(1) - y1pR((2-b1L):2:-1)
                y1uR(b1L: 0) = 2.*y1pR(1) - y1uR((1-b1L):1:-1)
            end if
        end if

        if (BC_1U_global == -2 .or. BC_1U_global > 0) then
            if (BC_1L_global > 0 .and. g == 1) then
                y1pR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1pR((Mmax-1):(Mmax  -b1U):-1)
                y1uR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1uR((Mmax-2):(Mmax-1-b1U):-1)
            else
                y1pR((Mmax+1):(Mmax+b1U)) = 2.*y1pR(Mmax) - y1pR((Mmax-1):(Mmax  -b1U):-1)
                y1uR( Mmax   :(Mmax+b1U)) = 2.*y1pR(Mmax) - y1uR((Mmax-1):(Mmax-1-b1U):-1)
            end if
        end if


        !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
        iiShift = (iB(1,g)-1)*(Nmax-1)

        x1pR(b1L:(Nmax+b1U),g) = y1pR((iiShift+b1L):(iiShift+Nmax+b1U))
        x1uR(b1L:(Nmax+b1U),g) = y1uR((iiShift+b1L):(iiShift+Nmax+b1U))



    x1p(:) = x1pR(:,1)
    x1u(:) = x1uR(:,1)

    ! Ist besser geeignet fuer Auswertung (Integrationsgewichte):
    do i = 1, N1
        dx1p(i) = x1u(i)-x1u(i-1)
    end do
    if (BC_1L > 0 .or. BC_1L == -2) dx1p(1 ) = x1u(1 )-x1p(1   )
    if (BC_1U > 0 .or. BC_1U == -2) dx1p(N1) = x1p(N1)-x1u(N1-1)

    do i = 1, N1-1
        dx1u(i) = x1p(i+1)-x1p(i)
    end do
    dx1u(0 ) = dx1u(1   ) ! Notwendig fuer Partikel-Rückkopplung
    dx1u(N1) = dx1u(N1-1)
    if (BC_1L == 0 .or. BC_1L == -1) dx1u(0 ) = x1p(1   )-x1p(0 )
    if (BC_1U == 0 .or. BC_1U == -1) dx1u(N1) = x1p(N1+1)-x1p(N1)

    !--- Smagorinsky LES model ---
    dx1pS(1:N1) = dy1p((iShift+1):(iShift+N1))
    dx1uS(0:N1) = dy1u((iShift+0):(iShift+N1))
    !===========================================================================================================



    !=== grobe Gitter (global) =================================================================================
    do g = 1, n_grids

        Nmax = NN(2,g)
        Mmax = NB(2,g)*(Nmax-1)+1

        dj   = (M2-1) / (Mmax-1)

        y2pR = 0.
        y2vR = 0.

        if (g == 1) then
            y2pR(1:M2) = y2p(1:M2)
            y2vR(0:M2) = y2v(0:M2)
        else
            do j = 1, Mmax
                y2pR(j) = y2p(1 + dj*j - dj  )
            end do
            do j = 1, Mmax-1
                y2vR(j) = y2p(1 + dj*j - dj/2)
            end do
        end if


        !--- periodische RB -------------------------------------------------------------------------------------
        if (BC_2L_global == -1) then
            y2pR(b2L:0) = y2pR((Mmax-1+b2L):(Mmax-1)) - L2
            y2vR(b2L:0) = y2vR((Mmax-1+b2L):(Mmax-1)) - L2

            y2pR(Mmax:(Mmax+b2U)) = y2pR(1:(1+b2U)) + L2
            y2vR(Mmax:(Mmax+b2U)) = y2vR(1:(1+b2U)) + L2
        end if


        !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
        if (BC_2L_global == -2 .or. BC_2L_global > 0) then
            if (BC_2L_global > 0 .and. g == 1) then
                y2pR(b2L: 0) = 2.*y2pR(1) - y2pR((2-b2L):2:-1)
                y2vR(b2L:-1) = 2.*y2pR(1) - y2vR((1-b2L):2:-1)
            else
                y2pR(b2L: 0) = 2.*y2pR(1) - y2pR((2-b2L):2:-1)
                y2vR(b2L: 0) = 2.*y2pR(1) - y2vR((1-b2L):1:-1)
            end if
        end if

        if (BC_2U_global == -2 .or. BC_2U_global > 0) then
            if (BC_2U_global > 0 .and. g == 1) then
                y2pR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2pR((Mmax-1):(Mmax  -b2U):-1)
                y2vR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2vR((Mmax-2):(Mmax-1-b2U):-1)
            else
                y2pR((Mmax+1):(Mmax+b2U)) = 2.*y2pR(Mmax) - y2pR((Mmax-1):(Mmax  -b2U):-1)
                y2vR( Mmax   :(Mmax+b2U)) = 2.*y2pR(Mmax) - y2vR((Mmax-1):(Mmax-1-b2U):-1)
            end if
        end if


        !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
        jjShift = (iB(2,g)-1)*(Nmax-1)

        x2pR(b2L:(Nmax+b2U),g) = y2pR((jjShift+b2L):(jjShift+Nmax+b2U))
        x2vR(b2L:(Nmax+b2U),g) = y2vR((jjShift+b2L):(jjShift+Nmax+b2U))


        !--- Schreiben ------------------------------------------------------------------------------------------
        if (rank == 0 .and. write_test_yes) then

            write(grid,'(i1.1)') g

            open(10, file='test_coord_p2_grid'//grid//'_restart'//restart_char//'.txt', status='UNKNOWN')
            open(11, file='test_coord_v2_grid'//grid//'_restart'//restart_char//'.txt', status='UNKNOWN')

            do jj = b2L, Nmax+b2U
                j = 1 + dj*(jj-1)
                write(10,'(2E25.17)') REAL(j+jShift)             , x2pR(jj,g)
                write(11,'(2E25.17)') REAL(j+jShift)+0.5*REAL(dj), x2vR(jj,g)
            end do

            close(10)
            close(11)

        end if

    end do

    x2p(:) = x2pR(:,1)
    x2v(:) = x2vR(:,1)

    do j = 1, N2
        dx2p(j) = x2v(j)-x2v(j-1)
    end do
    if (BC_2L > 0 .or. BC_2L == -2) dx2p(1 ) = x2v(1 )-x2p(1   )
    if (BC_2U > 0 .or. BC_2U == -2) dx2p(N2) = x2p(N2)-x2v(N2-1)

    do j = 1, N2-1
        dx2v(j) = x2p(j+1)-x2p(j)
    end do
    dx2v(0 ) = dx2v(1   ) ! Notwendig fuer Partikel-Rückkopplung
    dx2v(N2) = dx2v(N2-1)
    if (BC_2L == 0 .or. BC_2L == -1) dx2v(0 ) = x2p(1   )-x2p(0 )
    if (BC_2U == 0 .or. BC_2U == -1) dx2v(N2) = x2p(N2+1)-x2p(N2)

    !--- Smagorinsky LES model ---
    dx2pS(1:N2) = dy2p((jShift+1):(jShift+N2))
    dx2vS(0:N2) = dy2v((jShift+0):(jShift+N2))
    !===========================================================================================================


    if (dimens == 3) then
        !=== grobe Gitter (global) =================================================================================
        do g = 1, n_grids

            Nmax = NN(3,g)
            Mmax = NB(3,g)*(Nmax-1)+1

            dk   = (M3-1) / (Mmax-1)

            y3pR = 0.
            y3wR = 0.

            if (g == 1) then
                y3pR(1:M3) = y3p(1:M3)
                y3wR(0:M3) = y3w(0:M3)
            else
                do k = 1, Mmax
                    y3pR(k) = y3p(1 + dk*k - dk  )
                end do
                do k = 1, Mmax-1
                    y3wR(k) = y3p(1 + dk*k - dk/2)
                end do
            end if


            !--- periodische RB -------------------------------------------------------------------------------------
            if (BC_3L_global == -1) then
                y3pR(b3L:0) = y3pR((Mmax-1+b3L):(Mmax-1)) - L3
                y3wR(b3L:0) = y3wR((Mmax-1+b3L):(Mmax-1)) - L3

                y3pR(Mmax:(Mmax+b3U)) = y3pR(1:(1+b3U)) + L3
                y3wR(Mmax:(Mmax+b3U)) = y3wR(1:(1+b3U)) + L3
            end if


            !--- Symmetrie-RB oder feste Wände ----------------------------------------------------------------------
            if (BC_3L_global == -2 .or. BC_3L_global > 0) then
                if (BC_3L_global > 0 .and. g == 1) then
                    y3pR(b3L: 0) = 2.*y3pR(1) - y3pR((2-b3L):2:-1)
                    y3wR(b3L:-1) = 2.*y3pR(1) - y3wR((1-b3L):2:-1)
                else
                    y3pR(b3L: 0) = 2.*y3pR(1) - y3pR((2-b3L):2:-1)
                    y3wR(b3L: 0) = 2.*y3pR(1) - y3wR((1-b3L):1:-1)
                end if
            end if

            if (BC_3U_global == -2 .or. BC_3U_global > 0) then
                if (BC_3U_global > 0 .and. g == 1) then
                    y3pR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3pR((Mmax-1):(Mmax  -b3U):-1)
                    y3wR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3wR((Mmax-2):(Mmax-1-b3U):-1)
                else
                    y3pR((Mmax+1):(Mmax+b3U)) = 2.*y3pR(Mmax) - y3pR((Mmax-1):(Mmax  -b3U):-1)
                    y3wR( Mmax   :(Mmax+b3U)) = 2.*y3pR(Mmax) - y3wR((Mmax-1):(Mmax-1-b3U):-1)
                end if
            end if


            !--- Verteilen auf Blöcke -------------------------------------------------------------------------------
            kkShift = (iB(3,g)-1)*(Nmax-1)

            x3pR(b3L:(Nmax+b3U),g) = y3pR((kkShift+b3L):(kkShift+Nmax+b3U))
            x3wR(b3L:(Nmax+b3U),g) = y3wR((kkShift+b3L):(kkShift+Nmax+b3U))


            !--- Schreiben ------------------------------------------------------------------------------------------
            if (rank == 0 .and. write_test_yes) then

                write(grid,'(i1.1)') g

                open(10, file='test_coord_p3_grid'//grid//'_restart'//restart_char//'.txt', status='UNKNOWN')
                open(11, file='test_coord_w3_grid'//grid//'_restart'//restart_char//'.txt', status='UNKNOWN')

                do kk = b3L, Nmax+b3U
                    k = 1 + dk*(kk-1)
                    write(10,'(2E25.17)') REAL(k+kShift)             , x3pR(kk,g)
                    write(11,'(2E25.17)') REAL(k+kShift)+0.5*REAL(dk), x3wR(kk,g)
                end do

                close(10)
                close(11)

            end if

        end do

        x3p(:) = x3pR(:,1)
        x3w(:) = x3wR(:,1)

        do k = 1, N3
            dx3p(k) = x3w(k)-x3w(k-1)
        end do
        if (BC_3L > 0 .or. BC_3L == -2) dx3p(1 ) = x3w(1 )-x3p(1   )
        if (BC_3U > 0 .or. BC_3U == -2) dx3p(N3) = x3p(N3)-x3w(N3-1)

        do k = 1, N3-1
            dx3w(k) = x3p(k+1)-x3p(k)
        end do
        dx3w(0 ) = dx3w(1   ) ! Notwendig fuer Partikel-Rückkopplung
        dx3w(N3) = dx3w(N3-1)
        if (BC_3L == 0 .or. BC_3L == -1) dx3w(0 ) = x3p(1   )-x3p(0 )
        if (BC_3U == 0 .or. BC_3U == -1) dx3w(N3) = x3p(N3+1)-x3p(N3)

        !--- Smagorinsky LES model ---
        dx3pS(1:N3) = dy3p((kShift+1):(kShift+N3))
        dx3wS(0:N3) = dy3w((kShift+0):(kShift+N3))
    !===========================================================================================================
    else
        ! fuer Statistiken: ! TEST!!! ok?
        dy3p  = 1.
        dy3w  = 1.
        dx3p  = 1.
        dx3w  = 1.
        dx3pS = 1.
        dx3wS = 1.
    end if


end subroutine PI_getLocalCoordinates



end module cmod_geometry
