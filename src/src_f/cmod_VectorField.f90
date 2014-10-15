!>  Modul: cmod_VectorField
!!
!! impact functions \c Pimpact::VectorField e.g. scales, norms ...
!! \author huppd
module cmod_VectorField

    use iso_c_binding

!    use mod_vars, only: L1,L2,L3,M1,M2,M3

    implicit none

contains

    !> \brief calculate distance to immersed boundary
    !! \author huppd
    !! \param[in] x
    !! \param[in] y
    !! \param[in] z
    FUNCTION distance2ib( x, y, z, x0, y0, R, dr ) RESULT(dis)

        IMPLICIT NONE

        REAL, INTENT(IN)  ::  x,y,z

        REAL, intent(in)  ::  x0
        REAL, intent(in)  ::  y0

        REAL, intent(in)  ::  R

        REAL, intent(in)  ::  dr

        REAL              ::  z0

        REAL              ::  dis

        INTEGER           ::  m

        ! geometric properties of disc
        !        x0 = L1/2. + L1*amp*SIN(2.*pi*freq*subtime)
        !        y0 = L2/4.
        !        z0 = L3/2.

        !        R = L1/10.
        !  write(*,*) 'dr=',dr

        dis = SQRT( (x0-x)**2 + (y0-y)**2 )

        !        IF( dis <= R) THEN
        !            dis = 1.
        !        ELSE
        IF( Abs(dis-R) < dr) THEN
            dis = ABS(1./( 1. + EXP(1./(-ABS(dis-R)/dr) + 1./(1.-Abs(dis-R)/dr) ) ));
        !    dis = (1-(dis-R)/dr)
        ELSE
            dis = 0.
        END IF
    !    dis = 0.

    END FUNCTION distance2ib






    !> extracts field to field without ghostlayers or boundaries
    subroutine extract_dof(    &
        dimens,                &
        N,                     &
        bL,bU,                 &
        SU,NU,                 &
        SV,NV,                 &
        SW,NW,                 &
        phiiU,phiiV,phiiW,     &
        phioU,phioV,phioW ) bind ( c, name='VF_extract_dof' )

        implicit none

        integer(c_int), intent(in)    ::  dimens

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SU(3)
        integer(c_int), intent(in)    ::  NU(3)

        integer(c_int), intent(in)    ::  SV(3)
        integer(c_int), intent(in)    ::  NV(3)

        integer(c_int), intent(in)    ::  SW(3)
        integer(c_int), intent(in)    ::  NW(3)

        real(c_double),  intent(in)   ::  phiiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phiiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phiiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double),  intent(out)   ::  phioU(SU(1):NU(1),SU(2):NU(2),SU(3):NU(3))
        real(c_double),  intent(out)   ::  phioV(SV(1):NV(1),SV(2):NV(2),SV(3):NV(3))
        real(c_double),  intent(out)   ::  phioW(SW(1):NW(1),SW(2):NW(2),SW(3):NW(3))

        integer                       ::  i, j, k


        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                !pgi$ unroll = n:8
                do i = SU(1), NU(1)
                    phioU(i,j,k) =  phiiU(i,j,k)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                !pgi$ unroll = n:8
                do i = SV(1), NV(1)
                    phioV(i,j,k) =  phiiV(i,j,k)
                end do
            end do
        end do

        if (dimens == 3) then
            do k = SW(3), NW(3)
                do j = SW(2), NW(2)
                    !pgi$ unroll = n:8
                    do i = SW(1), NW(1)
                        phioW(i,j,k) =  phiiW(i,j,k)
                    end do
                end do
            end do
        end if

    end subroutine extract_dof


    subroutine extract_dof_reverse(    &
        dimens,                &
        N,                     &
        bL,bU,                 &
        SU,NU,                 &
        SV,NV,                 &
        SW,NW,                 &
        phiiU,phiiV,phiiW,     &
        phioU,phioV,phioW ) bind ( c, name='VF_extract_dof_reverse' )

        implicit none

        integer(c_int), intent(in)    ::  dimens

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SU(3)
        integer(c_int), intent(in)    ::  NU(3)

        integer(c_int), intent(in)    ::  SV(3)
        integer(c_int), intent(in)    ::  NV(3)

        integer(c_int), intent(in)    ::  SW(3)
        integer(c_int), intent(in)    ::  NW(3)


        real(c_double),  intent(in)  ::  phiiU(SU(1):NU(1),SU(2):NU(2),SU(3):NU(3))
        real(c_double),  intent(in)  ::  phiiV(SV(1):NV(1),SV(2):NV(2),SV(3):NV(3))
        real(c_double),  intent(in)  ::  phiiW(SW(1):NW(1),SW(2):NW(2),SW(3):NW(3))

        real(c_double),  intent(out)   ::  phioU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(out)   ::  phioV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(out)   ::  phioW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                       ::  i, j, k


        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                !pgi$ unroll = n:8
                do i = SU(1), NU(1)
                    phioU(i,j,k) =  phiiU(i,j,k)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                !pgi$ unroll = n:8
                do i = SV(1), NV(1)
                    phioV(i,j,k) =  phiiV(i,j,k)
                end do
            end do
        end do

        if (dimens == 3) then
            do k = SW(3), NW(3)
                do j = SW(2), NW(2)
                    !pgi$ unroll = n:8
                    do i = SW(1), NW(1)
                        phioW(i,j,k) =  phiiW(i,j,k)
                    end do
                end do
            end do
        end if

    end subroutine extract_dof_reverse










    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine VF_init_2DPulsatileXC(   &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        L2,                             &
        x2p,                            &
        re_, om, px,                    &
        phiU, phiV, phiW ) bind ( c, name='VF_init_2DPulsatileXC' )

        implicit none


        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L2

        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om
        real(c_double), intent(in)    ::  px

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: pi
        real :: Lh
        real :: mu
        real :: ny
        real :: c



        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        !  phiU = 0.
        !  phiV = 0.
        !  phiW = 0.
        pi = 4.*atan(1.)
        Lh  = L2/2.
        mu = sqrt( om*re_/2. )*Lh
        c  = px/( om*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    ny = sqrt( om*re_/2. )*( x2p(j)-Lh )
                    phiU(i,j,k) = -c*( -cos(ny)*cosh(ny)*sin(mu)*sinh(mu) +sin(ny)*sinh(ny)*cos(mu)*cosh(mu) )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_2DPulsatileXC


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine VF_init_2DPulsatileYC(   &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        L1,                             &
        x1,                             &
        re_, om_, px,                   &
        phiU, phiV, phiW ) bind ( c, name='VF_init_2DPulsatileYC' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    :: bL(3)
        integer(c_int), intent(in)    :: bU(3)

        integer(c_int), intent(in)    :: SU(3)
        integer(c_int), intent(in)    :: NU(3)

        integer(c_int), intent(in)    :: SV(3)
        integer(c_int), intent(in)    :: NV(3)

        integer(c_int), intent(in)    :: SW(3)
        integer(c_int), intent(in)    :: NW(3)

        real(c_double), intent(in)     :: L1
        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double), intent(in)    :: re_
        real(c_double), intent(in)    :: om_
        real(c_double), intent(in)    :: px

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: pi
        real :: Lh
        real :: mu
        real :: ny
        real :: c

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        !  phiU = 0.
        !  phiV = 0.
        !  phiW = 0.
        pi = 4.*atan(1.)
        Lh  = L1/2.
        mu = sqrt( om_*re_/2. )*Lh
        c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    ny = sqrt( om_*re_/2. )*( x1(i)-Lh )
                    phiV(i,j,k) = -c*( -cos(ny)*cosh(ny)*sin(mu)*sinh(mu) +sin(ny)*sinh(ny)*cos(mu)*cosh(mu) )
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_2DPulsatileYC


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine VF_init_2DPulsatileXS(   &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        L2,                     &
        x2,                     &
        re_, om_, px,           &
        phiU, phiV, phiW ) bind ( c, name='VF_init_2DPulsatileXS' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L2
        real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om_
        real(c_double), intent(in)    ::  px

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: pi
        real :: Lh
        real :: mu
        real :: ny
        real :: c

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        !  phiU = 0.
        !  phiV = 0.
        !  phiW = 0.
        pi = 4.*atan(1.)
        Lh  = L2/2.
        mu = sqrt( om_*re_/2. )*Lh
        c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    ny = sqrt( om_*Re_/2. )*( x2(j)-Lh )
                    phiU(i,j,k) = -c*( cos(ny)*cosh(ny)*cos(mu)*cosh(mu) +sin(ny)*sinh(ny)*sin(mu)*sinh(mu) ) + px/om_
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_2DPulsatileXS


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine VF_init_2DPulsatileYS(   &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        L1,                     &
        x1,                     &
        re_, om_, px,           &
        phiU,phiV,phiW ) bind ( c, name='VF_init_2DPulsatileYS' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L1
        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om_
        real(c_double), intent(in)    ::  px

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: pi
        real :: Lh
        real :: mu
        real :: ny
        real :: c

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        !  phiU = 0.
        !  phiV = 0.
        !  phiW = 0.
        pi = 4.*atan(1.)
        Lh  = L1/2.
        mu = sqrt( om_*re_/2. )*Lh
        c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    ny = sqrt( om_*re_/2. )*( x1(i)-Lh )
                    phiV(i,j,k) = -c*( cos(ny)*cosh(ny)*cos(mu)*cosh(mu) +sin(ny)*sinh(ny)*sin(mu)*sinh(mu) ) + px/om_
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_2DPulsatileYS



    !> \brief \f$ amp*\cos( \frac{2.*pi*x}{L_x} ) \f$
    subroutine VF_init_StreamingC(  &
        N,                          &
        bL,bU,                      &
        SU,NU,                      &
        SV,NV,                      &
        SW,NW,                      &
        L1,                         &
        x1,                         &
        amp,                        &
        phiU,phiV,phiW ) bind ( c, name='VF_init_StreamingC' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L1
        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double), intent(in)    ::  amp

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: pi

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        pi = 4.*atan(1.)

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = amp*cos( 2.*pi*x1(i)/L1 )
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_StreamingC



    !> \brief \f$ amp*\sin( \frac{2.*pi*x}{L_x} ) \f$
    subroutine VF_init_StreamingS(  &
        N,                          &
        bL,bU,                      &
        SU,NU,                      &
        SV,NV,                      &
        SW,NW,                      &
        L1,                         &
        x1,                         &
        amp,                        &
        phiU,phiV,phiW ) bind ( c, name='VF_init_StreamingS' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L1
        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double), intent(in)    ::  amp

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        integer                       ::  i, j, k

        real :: pi


        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        pi = 4.*atan(1.)
        !  Lh  = L1/2.
        !  mu = sqrt( om_*re_/2. )*Lh
        !  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = amp*sin( 2.*pi*x1(i)/L1 )
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_StreamingS



    subroutine VF_init_Vpoint(  &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        L,                      &
        x1u,                    &
        x2p,                    &
        sig,                    &
        phiU,phiV,phiW ) bind ( c, name='VF_init_Vpoint' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L(3)
        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)    :: sig

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !


        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = exp( -((x1u(i)-L(1)/2.)/sig)**2 -((x2p(j)-L(2)/2.)/sig/sig)**2 )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = 0.
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_Vpoint



    subroutine VF_init_Circle(  &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        L,                      &
        x1,                     &
        x2,                     &
        phiU,phiV,phiW ) bind ( c, name='VF_init_Circle' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))



        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = -(x2(j)-L(2)/2)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = x1(i)-L(1)/2
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_Circle



    subroutine VF_init_RankineVortex(   &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        L,                              &
        x1p,x2p,                        &
        x1u,x2v,                        &
        phiU,phiV,phiW ) bind ( c, name='VF_init_RankineVortex' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2v( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: circ
        real :: rad
        real :: r
        real :: pi


        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !
        pi = 4.*atan(1.)
        rad = L(1)/2./2.
        circ = 2.*pi*rad

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    r = sqrt( (x1u(i)-L(1)/2)**2 + (x2p(j)-L(2)/2)**2 )
                    if( r<= rad ) then
                        phiU(i,j,k) = -(x2p(j)-L(2)/2)/rad
                    else
                        phiU(i,j,k) = -(x2p(j)-L(2)/2)*rad/r/r
                    endif
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    r = sqrt( (x1p(i)-L(1)/2)**2 + (x2v(j)-L(2)/2)**2 )
                    if( r<=rad ) then
                        phiV(i,j,k) = (x1p(i)-L(1)/2)/rad
                    else
                        phiV(i,j,k) = (x1p(i)-L(1)/2)*rad/r/r
                    endif
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_RankineVortex



    subroutine VF_init_GaussianForcing1D(   &
        N,                                  &
        bL,bU,                              &
        SU,NU,                              &
        SV,NV,                              &
        SW,NW,                              &
        L1,                                 &
        x1u,                                &
        phiU,phiV,phiW ) bind ( c, name='VF_init_GaussianForcing1D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L1
        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: sig

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        sig = 0.2

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2  )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_GaussianForcing1D




    subroutine VF_init_GaussianForcing2D(   &
        N,                                  &
        bL,bU,                              &
        SU,NU,                              &
        SV,NV,                              &
        SW,NW,                              &
        L,                                  &
        x1p,x2p,                            &
        x1u,x2v,                            &
        phiU,phiV,phiW ) bind ( c, name='VF_init_GaussianForcing2D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2v( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real :: sig

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        sig = 0.2

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = exp( -((x1u(i)   )/sig)**2 - ((x2p(j)   )/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)-L(1))/sig)**2 - ((x2p(j)-L(2))/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)-L(1))/sig)**2 - ((x2p(j)   )/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)   )/sig)**2 - ((x2p(j)-L(2))/sig)**2 ) / sqrt(2.)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = exp( -((x1p(i)   )/sig)**2 -((x2v(j)   )/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)-L(1))/sig)**2 -((x2v(j)-L(2))/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)-L(1))/sig)**2 -((x2v(j)   )/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)   )/sig)**2 -((x2v(j)-L(2))/sig)**2 ) / sqrt(2.)
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_GaussianForcing2D




    subroutine VF_init_BoundaryFilter1D(    &
        N,                                  &
        bL,bU,                              &
        SU,NU,                              &
        SV,NV,                              &
        SW,NW,                              &
        L1,                                 &
        x1u,                                &
        phiU,phiV,phiW ) bind ( c, name='VF_init_BoundaryFilter1D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L1

        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: h
        real :: h2

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        h = 0.1
        h2 = h*h

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = 10*( exp( -(x1u(i)**2)/h2  ) + exp( -((x1u(i)-L1)**2)/h2  ) )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_BoundaryFilter1D



    subroutine VF_init_BoundaryFilter2D(    &
        N,                                  &
        bL,bU,                              &
        SU,NU,                              &
        SV,NV,                              &
        SW,NW,                              &
        L,                                  &
        x1p,x2p,                            &
        x1u,x2v,                            &
        phiU,phiV,phiW ) bind ( c, name='VF_init_BoundaryFilter2D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2v( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real :: h
        real :: h2

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !

        h = 0.1
        h2 = h*h

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = max( 10*( exp( -(x1u(i)**2)/h2  ) + exp( -((x1u(i)-L(1))**2)/h2  ) ) , &
                        10*( exp( -(x2p(j)**2)/h2  ) + exp( -((x2p(j)-L(2))**2)/h2  ) ) )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = max( 10*( exp( -(x1p(i)**2)/h2  ) + exp( -((x1p(i)-L(1))**2)/h2  ) ) , &
                        10*( exp( -(x2v(j)**2)/h2  ) + exp( -((x2v(j)-L(2))**2)/h2  ) ) )
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_BoundaryFilter2D



    !> \f$ u = \min( 4*\exp( -((x-xm)/rad)**2 -((y-ym)/rad)**2 ),1.) \f$
    subroutine VF_init_Disc(    &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        x1p,x2p,x3p,            &
        x1u,x2v,                &
        xm,ym, rad,             &
        phiU,phiV,phiW ) bind ( c, name='VF_init_Disc' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )
        real(c_double), intent(in)     :: x3p( bL(3):(N(3)+bU(3)) )

        real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2v( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double), intent(in)    :: xm
        real(c_double), intent(in)    :: ym
        real(c_double), intent(in)    :: rad


        real    :: dr

        integer ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !


!        dr = SQRT( ( 2*L1/REAL(M1-1) )**2 + ( 2*L2/REAL(M2-1) )**2 )
        dr = SQRT( ( 2*(x1p(2)-x1p(1)) )**2 + ( 2*(x2p(2)-x2p(1)) )**2 )

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = distance2ib( x1u(i),x2p(j),x3p(k),xm,ym,rad,dr )
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = distance2ib( x1p(i),x2v(j),x3p(k),xm,ym,rad,dr )
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_Disc


    subroutine VF_init_RotatingDisc(    &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        x1p,x2p,                        &
        xm,ym, omega,                   &
        phiU,phiV,phiW ) bind ( c, name='VF_init_RotatingDisc' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SU(3)
        integer(c_int), intent(in)     :: NU(3)

        integer(c_int), intent(in)     :: SV(3)
        integer(c_int), intent(in)     :: NV(3)

        integer(c_int), intent(in)     :: SW(3)
        integer(c_int), intent(in)     :: NW(3)

        real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
        real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

        real(c_double), intent(in)    :: xm
        real(c_double), intent(in)    :: ym
        real(c_double), intent(in)    :: omega

        real(c_double),  intent(inout) :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(inout) :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        ! note: - cf. sketch in file "usr_geometry.f90"
        !
        !         grid points in the domain and on the boundary
        !         |         |         |     velocity component
        !         |         |         |     |
        ! vel(S1U:N1U,S2U:N2U,S3U:N3U,1)
        ! vel(S1V:N1V,S2V:N2V,S3V:N3V,2)
        ! vel(S1W:N1W,S2W:N2W,S3W:N3W,3)
        !


        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                do i = SU(1), NU(1)
                    phiU(i,j,k) = -omega*(x2p(j)-ym)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                do i = SV(1), NV(1)
                    phiV(i,j,k) = omega*(x1p(i)-xm)
                end do
            end do
        end do

        do k = SW(3), NW(3)
            do j = SW(2), NW(2)
                do i = SW(1), NW(1)
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine VF_init_RotatingDisc



end module cmod_VectorField
