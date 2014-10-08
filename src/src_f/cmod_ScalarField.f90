!>  Modul: cmod_ScalarVector
!!
!! mature
!! impact functions \c Pimpact::ScalarVector e.g. scales, norms ...
!! \author huppd
module cmod_ScalarField

    use iso_c_binding

    implicit none

contains


    !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi first vector from which is taken the product
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    subroutine SF_add(      &
        N,                  &
        bL,bU,              &
        SS,NN,              &
        phi,phi1,phi2,      &
        scalar1,scalar2 )   &
        bind ( c, name='SF_add' )

        implicit none

        integer(c_int), intent(in)     :: N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)


        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in )   :: phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in )   :: phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        real(c_double),  intent(in)   ::  scalar1
        real(c_double),  intent(in)   ::  scalar2

        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar1*phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))+scalar2*phi2(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

    end subroutine SF_add


    !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi first vector from which is taken the product
    !! \param[in] phi1 first vector from which is taken the product
    subroutine SF_add2( &
        N,                  &
        bL,bU,              &
        SS,NN,              &
        phi,phi1,           &
        scalar1,scalar2 )   &
        bind ( c, name='SF_add2' )

        implicit none

        integer(c_int), intent(in)     :: N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)


        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in )   :: phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        real(c_double),  intent(in)   ::  scalar1
        real(c_double),  intent(in)   ::  scalar2


        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar1*phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))+scalar2*phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

    end subroutine SF_add2


    !> \brief copys absolute value from \c phi1 in \phi
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
    subroutine SF_abs(  &
        N,              &
        bL, bU,         &
        SS, NN,         &
        phi,            &
        phi1 ) bind ( c, name='SF_abs' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
        real(c_double), intent(in)    ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = abs( phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) )

    end subroutine SF_abs


    !> \brief copys reciprocal value from \c phi1 in \phi
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    subroutine SF_reciprocal(   &
        N,                      &
        bL, bU,                 &
        SS, NN,                 &
        phi,                    &
        phi1 ) bind ( c, name='SF_reciprocal' )

        implicit none


        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
        real(c_double), intent(in)    ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        integer                       :: i,j,k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    if( phi1(i,j,k).eq. 0. ) then
                        phi(i,j,k) = 0.
                    else
                        phi(i,j,k) = 1./phi1(i,j,k)
                    endif
                end do
            end do
        end do

    end subroutine SF_reciprocal


    !> \brief scales Field
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    subroutine SF_scale(    &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        phi,                &
        scalar) bind ( c, name='SF_scale' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(in) ::  scalar


        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar*phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

    end subroutine SF_scale


    !> \brief scales \c phi with \phi1.
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
    !! \param[in] inf_yes if true infinity norm is computed
    !! \param[in] two_yes if trhue two norm is computed
    subroutine SF_scale2(   &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        phi, phi1 ) bind ( c, name='SF_scale2' )

        implicit none


        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(inout) ::  phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
        real(c_double), intent(in)    ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))


        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))        &
            = phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))  &
            * phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

    end subroutine SF_scale2


    !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
    subroutine SF_dot(  &
        N,              &
        bL, bU,         &
        SS, NN,         &
        phi1,phi2,      &
        scalar ) bind ( c, name='SF_dot' )

        implicit none


        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in)    ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
        real(c_double), intent(in)    ::  phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(inout) ::  scalar
        integer                       ::  i, j, k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
                end do
            end do
        end do

    end subroutine SF_dot


    !> \brief computes two or infinity norm( get is misleading)
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi vector, from which the norm is taken
    !! \param[out] norm get the two norm of phi
    ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
    subroutine SF_comp1Norm(    &
        N,                      &
        bL, bU,                 &
        SS, NN,                 &
        phi,                    &
        norm ) bind (c, name='SF_comp1Norm')

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in)    ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(inout) ::  norm

        integer                       ::  i, j, k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    norm = norm + abs( phi(i,j,k) )
                end do
            end do
        end do

    end subroutine SF_comp1Norm


    !> \brief computes two or infinity norm( get is misleading)
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi vector, from which the norm is taken
    !! \param[out] norm get the two norm of phi
    ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
    subroutine SF_comp2Norm(    &
        N,                      &
        bL, bU,                 &
        SS, NN,                 &
        phi,                    &
        norm ) bind (c, name='SF_comp2Norm')

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in)    ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(inout) ::  norm

        integer                       ::  i, j, k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    norm = norm + phi(i,j,k)**2
                end do
            end do
        end do

    end subroutine SF_comp2Norm


    !> \brief computes two or infinity norm( get is misleading)
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi vector, from which the norm is taken
    !! \param[out] norm gets the infinity norm of phi
    subroutine SF_compInfNorm(  &
        N,                      &
        bL, bU,                 &
        SS, NN,                 &
        phi,                    &
        norm ) bind (c, name='SF_compInfNorm')

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in)    ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(inout) ::  norm

        integer                       ::  i, j, k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    norm = MAX( ABS(phi(i,j,k)),norm )
                end do
            end do
        end do

    end subroutine SF_compInfNorm


    !> \brief computes two or infinity norm( get is misleading)
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi vector, from which the norm is taken
    !! \param[in] weights weights
    !! \param[out] norm result
    subroutine SF_weightedNorm( &
        N,                      &
        bL, bU,                 &
        SS, NN,                 &
        phi,                    &
        weights,                &
        norm ) bind (c, name='SF_weightedNorm')

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(in)    ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
        real(c_double), intent(in)    ::  weights(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(inout) ::  norm

        integer                     ::  i, j, k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    norm = norm + ( weights(i,j,k)*phi(i,j,k) )**2
                end do
            end do
        end do

    end subroutine SF_weightedNorm


    !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    !! \param[in] phi first vector from which is taken the product
    !! \param[in] phi1 first vector from which is taken the product
    !! \param[in] phi2 second vector from which is taken the product
    subroutine SF_assign(   &
        N,                  &
        bL,bU,              &
        SS,NN,              &
        phi,phi1 )     &
        bind ( c, name='SF_assign' )

        implicit none

        integer(c_int), intent(in)     :: N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
        real(c_double),  intent(in )   :: phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

    end subroutine SF_assign


    !> \brief random values in [-0.5,0.5]
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    subroutine SF_random(   &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        phi) bind ( c, name='SF_random' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)    ::  bL(3)
        integer(c_int), intent(in)    ::  bU(3)

        integer(c_int), intent(in)    ::  SS(3)
        integer(c_int), intent(in)    ::  NN(3)

        real(c_double), intent(out)   ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))


        call Random_number( phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) )

        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) - 0.5

    end subroutine SF_random


    !> \brief inits all with scalar
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    subroutine SF_init( &
        N,              &
        bL, bU,         &
        SS, NN,         &
        phi,            &
        scalar) bind ( c, name='SF_init' )

        implicit none

        integer(c_int), intent(in)      ::  N(3)

        integer(c_int), intent(in)      ::  bL(3)
        integer(c_int), intent(in)      ::  bU(3)

        integer(c_int), intent(in)      ::  SS(3)
        integer(c_int), intent(in)      ::  NN(3)

        real(c_double), intent(out)   ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        real(c_double), intent(in) ::  scalar


        phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar

    end subroutine SF_init



    !> \brief writes Field to std out
    !!
    !! \param[in] N ammount of local elements in 1-direction
    !! \param[in] SS start index
    !! \param[in] NN end index
    !! \param[in] bL start index of storage in 1-direction
    !! \param[in] bU end offset of storage in 1-direction
    subroutine SF_print(    &
        N,                  &
        bL, bU,             &
        SS, NN,             &
        phi) bind ( c, name='SF_print' )

        implicit none

        integer(c_int), intent(in)   ::  N(3)

        integer(c_int), intent(in)   ::  bL(3)
        integer(c_int), intent(in)   ::  bU(3)

        integer(c_int), intent(in)   ::  SS(3)
        integer(c_int), intent(in)   ::  NN(3)

        real(c_double), intent(in)   ::  phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

        integer                      ::  i,j,k


        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    write(*,*) i,j,k,phi(i,j,k)
                end do
            end do
        end do

    end subroutine SF_print


    !> \brief init vector field with 2d Poiseuille flow in x-direction
    !! \f[ u(y) = y*( L_y - y )*4/L_y/L_y \f]
    subroutine SF_init_2DPoiseuilleX(   &
        N,                              &
        bL,bU,                          &
        SS,NN,                          &
        L2,                             &
        x2,                             &
        phi ) bind ( c, name='SF_init_2DPoiseuilleX' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double), intent(in)     :: L2

        real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    phi(i,j,k) = x2(j)*( L2 - x2(j) )*4/L2/L2
                end do
            end do
        end do

    end subroutine SF_init_2DPoiseuilleX


    !> \brief init vector field with 2d Poiseuille flow in y-direction
    !! \f[ v(x) = x*( L_x - x )*4/L_x/L_x \f]
    subroutine SF_init_2DPoiseuilleY(   &
        N,                              &
        bL,bU,                          &
        SS,NN,                          &
        L1,                             &
        x1,                             &
        phi ) bind ( c, name='SF_init_2DPoiseuilleY' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double), intent(in)     :: L1

        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
                    phi(i,j,k) = x1(i)*( L1 - x1(i) )*4/L1/L1
                end do
            end do
        end do

    end subroutine SF_init_2DPoiseuilleY


    !> \brief init vector field with constant gradient in x-direction
    !! \f[ u(y) = y*( L_y - y )*4/L_y/L_y \f]
    subroutine SF_init_2DGradX( &
        N,                      &
        bL,bU,                  &
        SS,NN,                  &
        L,                      &
        x1,                     &
        phi ) bind ( c, name='SF_init_2DGradX' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )

        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
!                    phi(i,j,k) = 8*( x1(i)-L(1)/2 )/L(2)/L(2)/L(1)
                    phi(i,j,k) = ( x1(i)-L(1)/2 )
                end do
            end do
        end do

    end subroutine SF_init_2DGradX


    !> \brief init vector field with constant gradient in x-direction
    !! \f[ u(y) = y*( L_y - y )*4/L_y/L_y \f]
    subroutine SF_init_2DGradY( &
        N,                      &
        bL,bU,                  &
        SS,NN,                  &
        L,                      &
        x2,                     &
        phi ) bind ( c, name='SF_init_2DGradY' )

        implicit none

        integer(c_int), intent(in)    ::  N(3)

        integer(c_int), intent(in)     :: bL(3)
        integer(c_int), intent(in)     :: bU(3)

        integer(c_int), intent(in)     :: SS(3)
        integer(c_int), intent(in)     :: NN(3)

        real(c_double), intent(in)     :: L(3)

        real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )

        real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

        integer                ::  i, j, k

        !--- initial conditions for velocity ---
        do k = SS(3), NN(3)
            do j = SS(2), NN(2)
                !pgi$ unroll = n:8
                do i = SS(1), NN(1)
!                    phi(i,j,k) = 8.*( x2(j)-L(2)/2 )/L(1)/L(1)/L(2)
                    phi(i,j,k) = ( x2(j)-L(2)/2 )
                end do
            end do
        end do

    end subroutine SF_init_2DGradY

end module cmod_ScalarField
