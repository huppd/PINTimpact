!>  Modul: cmod_VectorField
!!
!! impact functions \c Pimpact::VectorField e.g. scales, norms ...
!! \author huppd
module cmod_VectorField

    use iso_c_binding
    use mpi

    use mod_vars, only: x1p,x1u,x2p,x2v,x3p,L1,L2,L3

    implicit none
!  public get_norms

contains



    !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
    !! \param[in] dimens dimension
    !! \param[in] N ammount of local elements
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] SU start index in 1-direction
    !! \param[in] NU end index in 1-direction
    !! \param[in] SV start index in 1-direction
    !! \param[in] NV end index in 2-direction
    !! \param[in] SW start index in 1-direction
    !! \param[in] NW end index in 3-direction
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    !! \param[out] scalar dot product
    subroutine product_scalar_vel(    &
        dimens,                         &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        !    phi1,                           &
        !    phi2,                           &
        phi1U,phi1V,phi1W,              &
        phi2U,phi2V,phi2W,              &
        scalar ) bind ( c, name='VF_dot' )

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

        real(c_double),  intent(in)   ::  phi1U(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phi1V(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phi1W(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double),  intent(in)   ::  phi2U(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phi2V(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   ::  phi2W(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        !  real(c_double),  intent(in)   ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)),3)
        !  real(c_double),  intent(in)   ::  phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)),3)

        !  real(c_double),               ::  scalars(4)
        real(c_double),  intent(out)  ::  scalar

        integer                       ::  i, j, k

        scalar = 0.

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                !pgi$ unroll = n:8
                do i = SU(1), NU(1)
                    scalar = scalar + phi1U(i,j,k)*phi2U(i,j,k)
                !        scalar = scalar + phi1(i,j,k,1)*phi2(i,j,k,1)
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                !pgi$ unroll = n:8
                do i = SV(1), NV(1)
                    scalar = scalar + phi1V(i,j,k)*phi2V(i,j,k)
                !        scalar = scalar + phi1(i,j,k,2)*phi2(i,j,k,2)
                end do
            end do
        end do

        if (dimens == 3) then
            do k = SW(3), NW(3)
                do j = SW(2), NW(2)
                    !pgi$ unroll = n:8
                    do i = SW(1), NW(1)
                        scalar = scalar + phi1W(i,j,k)*phi2W(i,j,k)
                    !          scalar = scalar + phi1(i,j,k,3)*phi2(i,j,k,3)
                    end do
                end do
            end do
        end if

    end subroutine product_scalar_vel



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
        real(c_double),  intent(out)   ::  phioW(SV(1):NW(1),SW(2):NW(2),SW(3):NW(3))

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


    !> brief computes two or infinity norm
    !!
    !! \param[in] COMM_CART mpi communicator(obsolete) communications should be done level above
    !! \param[in] dimens dimension
    !! \param[in] N ammount of local elements
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] SU start index in 1-direction
    !! \param[in] NU end index in 1-direction
    !! \param[in] SV start index in 1-direction
    !! \param[in] NV end index in 2-direction
    !! \param[in] SW start index in 1-direction
    !! \param[in] NW end index in 3-direction
    !! \param[in] phiU velocity vector, from which the norm is taken
    !! \param[in] phiV velocity vector, from which the norm is taken
    !! \param[in] phiW velocity vector, from which the norm is taken
    !! \param[in] inf_yes if true infinity norm is computed
    !! \param[in] two_yes if trhue two norm is computed
    !! \param[out] normInf gets the infinity norm of phi
    !! \param[out] normTwo get the two norm of phi
    subroutine get_norms_vel( &
        COMM_CART,              &
        dimens,                 &
        N,                      &
        bL,bU,                  &
        SU,NU,                  &
        SV,NV,                  &
        SW,NW,                  &
        phiU,phiV,phiW,         &
        inf_yes,two_yes,        &
        normInf,normTwo ) bind (c,name='VF_compNorm')

        implicit none

        integer(c_int), intent(in)    :: COMM_CART

        integer(c_int), intent(in)    :: dimens

        integer(c_int), intent(in)    :: N(3)

        integer(c_int), intent(in)    :: bL(3)
        integer(c_int), intent(in)    :: bU(3)

        integer(c_int), intent(in)    :: SU(3)
        integer(c_int), intent(in)    :: NU(3)

        integer(c_int), intent(in)    :: SV(3)
        integer(c_int), intent(in)    :: NV(3)

        integer(c_int), intent(in)    :: SW(3)
        integer(c_int), intent(in)    :: NW(3)

        real(c_double),  intent(in)   :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        logical(c_bool), intent(in)   :: inf_yes
        logical(c_bool), intent(in)  ::  two_yes

        real(c_double),  intent(out) ::  normInf
        real(c_double),  intent(out) ::  normTwo

        real(c_double)               ::  normInf_global, normTwo_global
        integer(c_int)               ::  i, j, k

        integer                      ::  merror


        if (inf_yes .and. two_yes) then

            normInf = 0.
            normTwo = 0.

            do k = SU(3), NU(3)
                do j = SU(2), NU(2)
                    !pgi$ unroll = n:8
                    do i = SU(1), NU(1)
                        normInf = MAX(ABS(phiU(i,j,k)),normInf)
                        normTwo = normTwo + phiU(i,j,k)**2
                    end do
                end do
            end do
            do k = SV(3), NV(3)
                do j = SV(2), NV(2)
                    !pgi$ unroll = n:8
                    do i = SV(1), NV(1)
                        normInf = MAX(ABS(phiV(i,j,k)),normInf)
                        normTwo = normTwo + phiV(i,j,k)**2
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = SW(3), NW(3)
                    do j = SW(2), NW(2)
                        !pgi$ unroll = n:8
                        do i = SW(1), NW(1)
                            normInf = MAX(ABS(phiW(i,j,k)),normInf)
                            normTwo = normTwo + phiW(i,j,k)**2
                        end do
                    end do
                end do
            end if

            ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normInf = normInf_global
            normTwo = normTwo_global

        else if (inf_yes) then

            normInf = 0.

            do k = SU(3), NU(3)
                do j = SU(2), NU(2)
                    !pgi$ unroll = n:8
                    do i = SU(1), NU(1)
                        normInf = MAX(ABS(phiU(i,j,k)),normInf)
                    end do
                end do
            end do

            do k = SV(3), NV(3)
                do j = SV(2), NV(2)
                    !pgi$ unroll = n:8
                    do i = SW(1), NW(1)
                        normInf = MAX(ABS(phiV(i,j,k)),normInf)
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = SW(3), NW(3)
                    do j = SW(2), NW(2)
                        !pgi$ unroll = n:8
                        do i = SW(1), NW(1)
                            normInf = MAX(ABS(phiW(i,j,k)),normInf)
                        end do
                    end do
                end do
            end if


            call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
            normInf = normInf_global

        else if (two_yes) then

            normTwo = 0.

            do k = SU(3), NU(3)
                do j = SU(2), NU(2)
                    !pgi$ unroll = n:8
                    do i = SU(1), NU(1)
                        normTwo = normTwo + phiU(i,j,k)**2
                    end do
                end do
            end do
            do k = SV(3), NV(3)
                do j = SV(2), NV(2)
                    !pgi$ unroll = n:8
                    do i = SV(1), NV(1)
                        normTwo = normTwo + phiV(i,j,k)**2
                    end do
                end do
            end do

            if (dimens == 3) then
                do k = SW(3), NW(3)
                    do j = SW(2), NW(2)
                        !pgi$ unroll = n:8
                        do i = SW(1), NW(1)
                            normTwo = normTwo + phiW(i,j,k)**2
                        end do
                    end do
                end do
            end if

            call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
            normTwo = normTwo_global

        end if

    end subroutine get_norms_vel


    !> brief computes weighted two norm
    !!
    !! \param[in] COMM_CART mpi communicator(obsolete) communications should be done level above
    !! \param[in] dimens dimension
    !! \param[in] N ammount of local elements
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] SU start index in 1-direction
    !! \param[in] NU end index in 1-direction
    !! \param[in] SV start index in 1-direction
    !! \param[in] NV end index in 2-direction
    !! \param[in] SW start index in 1-direction
    !! \param[in] NW end index in 3-direction
    !! \param[in] phiU velocity vector, from which the norm is taken
    !! \param[in] phiV velocity vector, from which the norm is taken
    !! \param[in] phiW velocity vector, from which the norm is taken
    !! \param[out] norm gets the weighted two norm of phi
    subroutine get_weighted_norm_vel( &
        COMM_CART,                      &
        dimens,                         &
        N,                              &
        bL,bU,                          &
        SU,NU,                          &
        SV,NV,                          &
        SW,NW,                          &
        phiU,phiV,phiW,                 &
        wU,wV,wW,                       &
        norm ) bind (c,name='VF_weightedNorm')

        implicit none

        integer(c_int), intent(in)    :: COMM_CART

        integer(c_int), intent(in)    :: dimens

        integer(c_int), intent(in)    :: N(3)

        integer(c_int), intent(in)    :: bL(3)
        integer(c_int), intent(in)    :: bU(3)

        integer(c_int), intent(in)    :: SU(3)
        integer(c_int), intent(in)    :: NU(3)

        integer(c_int), intent(in)    :: SV(3)
        integer(c_int), intent(in)    :: NV(3)

        integer(c_int), intent(in)    :: SW(3)
        integer(c_int), intent(in)    :: NW(3)

        real(c_double),  intent(in)   :: phiU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: phiV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: phiW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double),  intent(in)   :: wU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: wV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in)   :: wW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        real(c_double),  intent(out) ::  norm

        real(c_double)               ::  norm_global
        integer                      ::  i, j, k

        integer                      ::  merror


        norm = 0.

        do k = SU(3), NU(3)
            do j = SU(2), NU(2)
                !pgi$ unroll = n:8
                do i = SU(1), NU(1)
                    norm = norm + (wU(i,j,k)*phiU(i,j,k))**2
                end do
            end do
        end do

        do k = SV(3), NV(3)
            do j = SV(2), NV(2)
                !pgi$ unroll = n:8
                do i = SW(1), NW(1)
                    norm = norm + (wV(i,j,k)*phiV(i,j,k))**2
                end do
            end do
        end do

        if (dimens == 3) then
            do k = SW(3), NW(3)
                do j = SW(2), NW(2)
                    !pgi$ unroll = n:8
                    do i = SW(1), NW(1)
                        norm = norm + (wW(i,j,k)*phiW(i,j,k))**2
                    end do
                end do
            end do
        end if

        ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
        call MPI_ALLREDUCE(norm,norm_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
        norm = norm_global

    end subroutine get_weighted_norm_vel




    !> \brief init vector field with 2d Poiseuille flow in x-direction
    subroutine cinit_2DPoiseuilleX(   &
        N,      &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_2DPoiseuilleX' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W

        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))

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
        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = x2p(j)*( L2 - x2p(j) )*4/L2/L2
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPoiseuilleX


    !> \brief init vector field with 2d Poiseuille flow in y-direction
    subroutine cinit_2DPoiseuilleY(   &
        !    N1,N2,N3,                       &
        N,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_2DPoiseuilleY' )

        implicit none

        !  integer(c_int), intent(in)    ::  N1
        !  integer(c_int), intent(in)    ::  N2
        !  integer(c_int), intent(in)    ::  N3
        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W

        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0.
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = x1p(i)*( L1 - x1p(i) )*4/L1/L1
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPoiseuilleY


    !> \brief init vector field with a zero flow
    subroutine cinit_zero(            &
        !    N1,N2,N3,                       &
        N,  &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_Zero' )

        implicit none

        !  integer(c_int), intent(in)    ::  N1
        !  integer(c_int), intent(in)    ::  N2
        !  integer(c_int), intent(in)    ::  N3
        integer(c_int), intent(in)    ::  N(3)


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N(1)+b1U),b2L:(N(2)+b2U),b3L:(N(3)+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0.0
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = 0.0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_zero


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine cinit_2DPulsatileXC(   &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW, re_, om, px ) bind ( c, name='VF_init_2DPulsatileXC' )

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U


        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om
        real(c_double), intent(in)    ::  px

        real :: pi
        real :: Lh
        real :: mu
        real :: nu
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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    nu = sqrt( om*re_/2. )*( x2p(j)-Lh )
                    phiU(i,j,k) = -c*( -cos(nu)*cosh(nu)*sin(mu)*sinh(mu) +sin(nu)*sinh(nu)*cos(mu)*cosh(mu) )
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPulsatileXC


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine cinit_2DPulsatileYC(   &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileYC' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W

        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om_
        real(c_double), intent(in)    ::  px

        real :: pi
        real :: Lh
        real :: mu
        real :: nu
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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    nu = sqrt( om_*re_/2. )*( x1p(i)-Lh )
                    phiV(i,j,k) = -c*( -cos(nu)*cosh(nu)*sin(mu)*sinh(mu) +sin(nu)*sinh(nu)*cos(mu)*cosh(mu) )
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPulsatileYC


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine cinit_2DPulsatileXS(   &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileXS' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om_
        real(c_double), intent(in)    ::  px

        real :: pi
        real :: Lh
        real :: mu
        real :: nu
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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    nu = sqrt( om_*Re_/2. )*( x2p(j)-Lh )
                    phiU(i,j,k) = -c*( cos(nu)*cosh(nu)*cos(mu)*cosh(mu) +sin(nu)*sinh(nu)*sin(mu)*sinh(mu) ) + px/om_
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPulsatileXS


    !> \brief init vector field with 2d pulsatile flow in x-direction
    !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
    !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
    !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
    subroutine cinit_2DPulsatileYS(   &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileYS' )

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(in)    ::  re_
        real(c_double), intent(in)    ::  om_
        real(c_double), intent(in)    ::  px

        real :: pi
        real :: Lh
        real :: mu
        real :: nu
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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    nu = sqrt( om_*re_/2. )*( x1p(i)-Lh )
                    phiV(i,j,k) = -c*( cos(nu)*cosh(nu)*cos(mu)*cosh(mu) +sin(nu)*sinh(nu)*sin(mu)*sinh(mu) ) + px/om_
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_2DPulsatileYS



    subroutine cinit_StreamingC(       &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW,                 &
        amp ) bind ( c, name='VF_init_StreamingC' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double), intent(in)    ::  amp

        real :: pi
        !  real :: Lh
        !  real :: mu
        !  real :: nu
        !  real :: c

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
        !  Lh  = L1/2.
        !  mu = sqrt( om_*re_/2. )*Lh
        !  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = amp*cos( 2.*pi*x1p(i)/L1 )
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_StreamingC



    subroutine cinit_StreamingS(       &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW,                 &
        amp ) bind ( c, name='VF_init_StreamingS' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double), intent(in)    ::  amp

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = 0
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = amp*sin( 2.*pi*x1p(i)/L1 )
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_StreamingS



    subroutine cinit_Vpoint( &
        N1,N2,N3,                         &
        S1U,S2U,S3U, N1U,N2U,N3U,         &
        S1V,S2V,S3V, N1V,N2V,N3V,         &
        S1W,S2W,S3W, N1W,N2W,N3W,         &
        b1L,b2L,b3L, b1U,b2U,b3U,         &
        phiU,phiV,phiW,                   &
        sig ) bind ( c, name='VF_init_Vpoint' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double), intent(in)    :: sig

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


        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = exp( -((x1u(i)-L1/2.)/sig)**2 -((x2p(j)-L2/2.)/sig/sig)**2 )
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = 0.
                !           phiV(i,j,k) = min(2.*exp( -((x1p(i)-L1/2.)/sig)**2 -((x2v(j)-L2/2.)/sig)**2 ),1.)
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_Vpoint



    subroutine cinit_Circle(          &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_Circle' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        !  real :: phi
        !  real :: d
        !  real :: Lh
        !  real :: mu
        !  real :: nu
        !  real :: c

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
        !  pi = 4.*atan(1.)
        !  Lh  = L1/2.
        !  mu = sqrt( om_*re_/2. )*Lh
        !  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )


        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    !           d = sqrt( (x1u(i)-L1/2)**2 + (x2p(j)-L2/2)**2 )
                    !           phi = atan( (x2p(j)-L2/2) / (x1u(i)-L1/2) )
                    !           phiU(i,j,k) = d*sin( phi )
                    phiU(i,j,k) = -(x2p(j)-L2/2)
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    !           d = sqrt( (x1p(i)-L1/2)**2 + (x2v(j)-L2/2)**2 )
                    !           phi = atan( (x2v(j)-L2/2) / (x1p(i)-L1/2) )
                    !           phiV(i,j,k) = d*cos( phi )
                    phiV(i,j,k) = x1p(i)-L1/2
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_Circle



    subroutine cinit_RankineVortex(        &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_RankineVortex' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real :: circ
        real :: rad
        real :: r
        !  real :: phi
        real :: pi

        !  real :: d
        !  real :: Lh
        !  real :: mu
        !  real :: nu
        !  real :: c

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
        rad = L1/2./2.
        circ = 2.*pi*rad
        !  Lh  = L1/2.
        !  mu = sqrt( om_*re_/2. )*Lh
        !  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )


        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    r = sqrt( (x1u(i)-L1/2)**2 + (x2p(j)-L2/2)**2 )
                    if( r<= rad ) then
                        phiU(i,j,k) = -(x2p(j)-L2/2)/rad
                    else
                        phiU(i,j,k) = -(x2p(j)-L2/2)*rad/r/r
                    endif

                !           phi = atan( (x2p(j)-L2/2) / (x1u(i)-L1/2) )
                !           phiU(i,j,k) = d*sin( phi )
                !           r = sqrt( (x1u-

                !           phiU(i,j,k) = -(x2p(j)-L2/2)
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    r = sqrt( (x1p(i)-L1/2)**2 + (x2v(j)-L2/2)**2 )
                    !           d = sqrt( (x1p(i)-L1/2)**2 + (x2v(j)-L2/2)**2 )
                    !           phi = atan( (x2v(j)-L2/2) / (x1p(i)-L1/2) )
                    !           phiV(i,j,k) = d*cos( phi )
                    if( r<=rad ) then
                        phiV(i,j,k) = (x1p(i)-L1/2)/rad
                    else
                        phiV(i,j,k) = (x1p(i)-L1/2)*rad/r/r
                    endif
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_RankineVortex



    subroutine cinit_GaussianForcing1D( &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_GaussianForcing1D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 )
                    phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2  )
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    !           phiV(i,j,k) = exp( -((x1u(i))/sig)**2 -((x2p(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    !           phiV(i,j,k) = exp( -((x1p(i)-L1/2)/sig)**2 -((x2v(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_GaussianForcing1D




    subroutine cinit_GaussianForcing2D( &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_GaussianForcing2D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = exp( -((x1u(i)   )/sig)**2 - ((x2p(j)   )/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)-L1)/sig)**2 - ((x2p(j)-L2)/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)-L1)/sig)**2 - ((x2p(j)   )/sig)**2 ) / sqrt(2.)    &
                        + exp( -((x1u(i)   )/sig)**2 - ((x2p(j)-L2)/sig)**2 ) / sqrt(2.)
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = exp( -((x1p(i)   )/sig)**2 -((x2v(j)   )/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)-L1)/sig)**2 -((x2v(j)-L2)/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)-L1)/sig)**2 -((x2v(j)   )/sig)**2 ) / sqrt(2.) &
                        + exp( -((x1p(i)   )/sig)**2 -((x2v(j)-L2)/sig)**2 ) / sqrt(2.)
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_GaussianForcing2D




    subroutine cinit_BoundaryFilter1D(        &
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_BoundaryFilter1D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 )
                    phiU(i,j,k) = 10*( exp( -(x1u(i)**2)/h2  ) + exp( -((x1u(i)-L1)**2)/h2  ) )
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    !           phiV(i,j,k) = exp( -((x1p(i)-L1/2)/sig)**2 -((x2v(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    phiV(i,j,k) = 0
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_BoundaryFilter1D



    subroutine cinit_BoundaryFilter2D(&
        N1,N2,N3,                       &
        S1U,S2U,S3U, N1U,N2U,N3U,       &
        S1V,S2V,S3V, N1V,N2V,N3V,       &
        S1W,S2W,S3W, N1W,N2W,N3W,       &
        b1L,b2L,b3L, b1U,b2U,b3U,       &
        phiU,phiV,phiW ) bind ( c, name='VF_init_BoundaryFilter2D' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

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

        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    !           phiU(i,j,k) = exp( -((x1u(i)-L1/2)/sig)**2 -((x2p(j)-L2/2)/sig)**2 )
                    phiU(i,j,k) = max( 10*( exp( -(x1u(i)**2)/h2  ) + exp( -((x1u(i)-L1)**2)/h2  ) ) , &
                        10*( exp( -(x2p(j)**2)/h2  ) + exp( -((x2p(j)-L2)**2)/h2  ) ) )
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    !           phiV(i,j,k) = exp( -((x1p(i)-L1/2)/sig)**2 -((x2v(j)-L2/2)/sig)**2 ) / sqrt(2.)
                    !           phiV(i,j,k) = 0
                    phiV(i,j,k) = max( 10*( exp( -(x1p(i)**2)/h2  ) + exp( -((x1p(i)-L1)**2)/h2  ) ) , &
                        10*( exp( -(x2v(j)**2)/h2  ) + exp( -((x2v(j)-L2)**2)/h2  ) ) )
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_BoundaryFilter2D



    subroutine cinit_Disc( &
        N1,N2,N3,                         &
        S1U,S2U,S3U, N1U,N2U,N3U,         &
        S1V,S2V,S3V, N1V,N2V,N3V,         &
        S1W,S2W,S3W, N1W,N2W,N3W,         &
        b1L,b2L,b3L, b1U,b2U,b3U,         &
        phiU,phiV,phiW,                   &
        xm,ym, rad ) bind ( c, name='VF_init_Disc' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double), intent(in)    :: xm
        real(c_double), intent(in)    :: ym
        real(c_double), intent(in)    :: rad

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


        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    phiU(i,j,k) = min(4*exp( -((x1u(i)-xm)/rad)**2 -((x2p(j)-ym)/rad)**2 ),1.)
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    phiV(i,j,k) = min(4*exp( -((x1p(i)-xm)/rad)**2 -((x2v(j)-ym)/rad)**2 ),1.)
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_Disc


    subroutine cinit_RotatingDisc( &
        N1,N2,N3,                         &
        S1U,S2U,S3U, N1U,N2U,N3U,         &
        S1V,S2V,S3V, N1V,N2V,N3V,         &
        S1W,S2W,S3W, N1W,N2W,N3W,         &
        b1L,b2L,b3L, b1U,b2U,b3U,         &
        phiU,phiV,phiW,                   &
        xm,ym, omega ) bind ( c, name='VF_init_RotatingDisc' )
        ! (basic subroutine)

        implicit none

        integer(c_int), intent(in)    ::  N1
        integer(c_int), intent(in)    ::  N2
        integer(c_int), intent(in)    ::  N3


        integer(c_int), intent(in)    ::  S1U
        integer(c_int), intent(in)    ::  S2U
        integer(c_int), intent(in)    ::  S3U

        integer(c_int), intent(in)    ::  N1U
        integer(c_int), intent(in)    ::  N2U
        integer(c_int), intent(in)    ::  N3U

        integer(c_int), intent(in)    ::  S1V
        integer(c_int), intent(in)    ::  S2V
        integer(c_int), intent(in)    ::  S3V

        integer(c_int), intent(in)    ::  N1V
        integer(c_int), intent(in)    ::  N2V
        integer(c_int), intent(in)    ::  N3V

        integer(c_int), intent(in)    ::  S1W
        integer(c_int), intent(in)    ::  S2W
        integer(c_int), intent(in)    ::  S3W

        integer(c_int), intent(in)    ::  N1W
        integer(c_int), intent(in)    ::  N2W
        integer(c_int), intent(in)    ::  N3W


        integer(c_int), intent(in)    ::  b1L
        integer(c_int), intent(in)    ::  b2L
        integer(c_int), intent(in)    ::  b3L

        integer(c_int), intent(in)    ::  b1U
        integer(c_int), intent(in)    ::  b2U
        integer(c_int), intent(in)    ::  b3U

        real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
        real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

        real(c_double), intent(in)    :: xm
        real(c_double), intent(in)    :: ym
        real(c_double), intent(in)    :: omega

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


        do k = S3U, N3U
            do j = S2U, N2U
                do i = S1U, N1U
                    !           phiU(i,j,k) = min(4*exp( -((x1u(i)-xm)/rad)**2 -((x2p(j)-ym)/rad)**2 ),1.)
                    phiU(i,j,k) = -omega*(x2p(j)-ym)
                end do
            end do
        end do

        do k = S3V, N3V
            do j = S2V, N2V
                do i = S1V, N1V
                    !           phiV(i,j,k) = min(4*exp( -((x1p(i)-ym)/rad)**2 -((x2v(j)-ym)/rad)**2 ),1.)
                    phiV(i,j,k) = omega*(x1p(i)-xm)
                end do
            end do
        end do

        do k = S3W, N3W
            do j = S2W, N2W
                do i = S1W, N1W
                    phiW(i,j,k) = 0.0
                end do
            end do
        end do

    end subroutine cinit_RotatingDisc

end module cmod_VectorField
