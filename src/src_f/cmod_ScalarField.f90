!>  Modul: cmod_ScalarVector
!!
!! Mature
!! Impact functions \c Pimpact::ScalarVector e.g. scales, norms ...
!! \author huppd
module cmod_ScalarField

  use iso_c_binding
  use mpi

  implicit none

contains


  !> \brief computes scalar product of two scalar fields
  !!
  !! \f[ phi = scalar1 *phi1 + scalar2*phi2 \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[out] phi first vector from which is taken the product
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] scalar1 factor for phi1
  !! \param[in] scalar2 factor for phi2
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


  !> \brief computes scalar product of two scalar fields
  !! 
  !! \f[ phi = scalar1*phi + scalar2*phi1 \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[inout] phi first vector from which is taken the product
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] scalar multiplicator of phi
  !! \param[in] scalar1 multiplicator of phi1
  subroutine SF_add2( &
      N,              &
      bL,bU,          &
      SS,NN,          &
      phi,            &
      phi1,           &
      scalar,         &
      scalar1 )       &
      bind ( c, name='SF_add2' )

    implicit none

    integer(c_int), intent(in)     :: N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)


    real(c_double),  intent(out)   :: phi ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )
    real(c_double),  intent(in )   :: phi1( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )


    real(c_double),  intent(in)   ::  scalar
    real(c_double),  intent(in)   ::  scalar1


    phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar*phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))+scalar1*phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

  end subroutine SF_add2



  !> \brief copys absolute value from \c phi1 in \c phi
  !!
  !! \f[ phi = | phi | \f]
  !! \param[in] N number of local grid points
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] phi field that gets assigned
  !! \param[in] phi1 field thats absolute values are taken
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



  !> \brief copys reciprocal value from \c phi1 in \c phi
  !!
  !! \f[ phi = \frac{1}{phi1} \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] phi field that gets assigned
  !! \param[in] phi1 field from which the reciprocal values are taken
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

    integer(c_int)                :: i,j,k


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
  !! \f[ phi = scalar \cdot phi \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[inout] phi scalar field
  !! \param[in] scalar scaling factor
  subroutine SF_scale(  &
      N,                &
      bL,               &
      bU,               &
      SS,               &
      NN,               &
      phi,              &
      scalar ) bind ( c, name='SF_scale' )

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



  !> \brief scales \c phi with \c phi1.
  !!
  !! \f[ phi = phi * phi1 \f]
  !! \param[in] N number of local grid points
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] phi first vector from which is taken the product
  !! \param[in] phi1 second vector from which is taken the product
  subroutine SF_scale2( &
      N,                &
      bL,               &
      bU,               &
      SS,               &
      NN,               &
      phi,              &
      phi1 ) bind ( c, name='SF_scale2' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)    ::  bL(3)
    integer(c_int), intent(in)    ::  bU(3)

    integer(c_int), intent(in)    ::  SS(3)
    integer(c_int), intent(in)    ::  NN(3)

    real(c_double), intent(inout) ::  phi ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)    ::  phi1( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )


    phi( SS(1):NN(1), SS(2):NN(2), SS(3):NN(3) )       &
      = phi ( SS(1):NN(1), SS(2):NN(2), SS(3):NN(3) )  &
      * phi1( SS(1):NN(1), SS(2):NN(2), SS(3):NN(3) )

  end subroutine SF_scale2



  !> \brief computes dot product of two scalar fields(neccessary condition belong to same sVS)
  !! 
  !! \f[ dot = (phi1 \cdot phi1 ) \f]
  !! \param[in] N number of local elements
  !! \param[in] bL lower bound
  !! \param[in] bU upper bound
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] dot scalar product
  subroutine SF_dot(  &
      N,              &
      bL,             &
      bU,             &
      SS,             &
      NN,             &
      phi1,           &
      phi2,           &
      dot ) bind ( c, name='SF_dot' )

    implicit none


    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)    ::  bL(3)
    integer(c_int), intent(in)    ::  bU(3)

    integer(c_int), intent(in)    ::  SS(3)
    integer(c_int), intent(in)    ::  NN(3)

    real(c_double), intent(in)    ::  phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))
    real(c_double), intent(in)    ::  phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3) ))

    real(c_double), intent(inout) ::  dot
    integer(c_int)                ::  i, j, k


    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          dot = dot + phi1(i,j,k)*phi2(i,j,k)
        end do
      end do
    end do

  end subroutine SF_dot



  !> \brief computes two or infinity norm( get is misleading)
  !! 
  !! \f[ norm = || phi ||_1 \f]
  !! \param[in] N number of local elements
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] phi vector, from which the norm is taken
  !! \param[out] norm get the two norm of phi
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

    integer(c_int)                ::  i, j, k


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
  !! 
  !! \f[ norm = || phi ||_2*^2 \f]
  !! \param[in] N number of local elements
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
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

    integer(c_int)                ::  i, j, k


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
  !!
  !! \f[ norm = || phi ||_\infty \f]
  !! \param[in] N number of local elements
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
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

    integer(c_int)                ::  i, j, k


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
  !!
  !! \f[ norm = ||weights * phi ||_2^2 \f]
  !! \param[in] N number of local elements
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
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

    integer(c_int)                ::  i, j, k


    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          norm = norm + ( weights(i,j,k)*phi(i,j,k) )**2
        end do
      end do
    end do

  end subroutine SF_weightedNorm



  !> \brief assigns field
  !!
  !! \f[ phi = phi1 \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage in
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[out] phi field that gets assigned
  !! \param[in] phi1 field from which phi is assigned
  subroutine SF_assign( &
      N,                &
      bL,               &
      bU,               &
      SS,               &
      NN,               &
      phi,              &
      phi1 )            &
      bind ( c, name='SF_assign' )

    implicit none

    integer(c_int), intent(in)  :: N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(out) :: phi ( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in ) :: phi1( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )


    phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

  end subroutine SF_assign



  !> \brief random values in [-0.5,0.5]
  !!
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage 
  !! \param[in] bU end offset of storage  
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[inout] phi scalar field 
  subroutine SF_random( &
      N,                &
      bL,               &
      bU,               &
      SS,               &
      NN,               &
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
  !!
  !! \f[ phi = scalar \f]
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[out] phi scalar field
  !! \param[in] scalar init value
  subroutine SF_init( &
      N,              &
      bL, bU,         &
      SS, NN,         &
      phi,            &
      scalar) bind ( c, name='SF_init' )

    implicit none

    integer(c_int), intent(in) ::  N(3)

    integer(c_int), intent(in) ::  bL(3)
    integer(c_int), intent(in) ::  bU(3)

    integer(c_int), intent(in) ::  SS(3)
    integer(c_int), intent(in) ::  NN(3)

    real(c_double), intent(out)::  phi( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )

    real(c_double), intent(in) ::  scalar


    phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar

  end subroutine SF_init



  !> \brief writes Field to std out
  !! \deprecated
  !!
  !! \param[in] N number of local elements
  !! \param[in] bL start index of storage
  !! \param[in] bU end offset of storage
  !! \param[in] SS start index
  !! \param[in] NN end index
  !! \param[in] phi scalar field
  subroutine SF_print(  &
      N,                &
      bL,               &
      bU,               &
      SS,               &
      NN,               &
      phi) bind ( c, name='SF_print' )

    implicit none

    integer(c_int), intent(in) ::  N(3)

    integer(c_int), intent(in) ::  bL(3)
    integer(c_int), intent(in) ::  bU(3)

    integer(c_int), intent(in) ::  SS(3)
    integer(c_int), intent(in) ::  NN(3)

    real(c_double), intent(in) ::  phi( bL(1):(N(1)+bU(1)), bL(2):(N(2)+bU(2)), bL(3):(N(3)+bU(3)) )

    integer(c_int)             ::  i,j,k

    write(*,*) "i,           j,           k,           phi(i,j,k)"

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
  !!
  !! \f[ u(x) = x*( L - x )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleX(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleX' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(1):(N(1)+bU(1)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(i)*( L - x(i) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleX



  !> \brief init vector field with 2d Poiseuille flow in y-direction
  !!
  !! \f[ u(y) = y*( L - y )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleY(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleY' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(2):(N(2)+bU(2)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(j)*( L - x(j) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleY



  !> \brief init vector field with 2d Poiseuille flow in y-direction
  !!
  !! \f[ u(z) = y*( L - z )*4/L/L \f]
  subroutine SF_init_2DPoiseuilleZ(   &
      N,                              &
      bL,bU,                          &
      SS,NN,                          &
      L,                              &
      x,                              &
      phi ) bind ( c, name='SF_init_2DPoiseuilleZ' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)

    real(c_double), intent(in)     :: L

    real(c_double), intent(in)     :: x( bL(3):(N(3)+bU(3)) )

    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = x(k)*( L - x(k) )*4/L/L
        end do
      end do
    end do

  end subroutine SF_init_2DPoiseuilleZ



  !> \brief init vector field with constant gradient in x-direction
  !!
  !! \f[ u(x) = \frac{x}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradX( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradX' )

    implicit none

    integer(c_int), intent(in)  ::  N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: L

    real(c_double), intent(in)  :: x( bL(1):(N(1)+bU(1)) )

    real(c_double), intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)  :: alpha

    integer(c_int)              ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) =  alpha*( x(i) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradX



  !> \brief init vector field with constant gradient in y-direction
  !!
  !! \f[ u(y) = \frac{y}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradY( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradY' )

    implicit none

    integer(c_int), intent(in)  ::  N(3)

    integer(c_int), intent(in)  :: bL(3)
    integer(c_int), intent(in)  :: bU(3)

    integer(c_int), intent(in)  :: SS(3)
    integer(c_int), intent(in)  :: NN(3)

    real(c_double), intent(in)  :: L

    real(c_double), intent(in)  :: x( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)  :: alpha

    integer(c_int)                 ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = alpha*( x(j) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradY



  !> \brief init vector field with constant gradient in z-direction
  !!
  !! \f[ u(z) = \frac{z}{L}-\frac{1}{2}  \f]
  subroutine SF_init_2DGradZ( &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      L,                      &
      x,                      &
      phi,                    &
      alpha ) bind ( c, name='SF_init_2DGradZ' )

    implicit none

    integer(c_int), intent(in)   ::  N(3)

    integer(c_int), intent(in)   :: bL(3)
    integer(c_int), intent(in)   :: bU(3)

    integer(c_int), intent(in)   :: SS(3)
    integer(c_int), intent(in)   :: NN(3)

    real(c_double), intent(in)   :: L

    real(c_double), intent(in)   :: x( bL(3):(N(3)+bU(3)) )

    real(c_double),  intent(out) :: phi( bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)   :: alpha

    integer(c_int)               ::  i, j, k

    !--- initial conditions for velocity ---
    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          phi(i,j,k) = alpha*( x(k) - L/2. )
        end do
      end do
    end do

  end subroutine SF_init_2DGradZ



  !> \brief init forcing point
  !!
  !! \f[ phi(x) = amp* \exp( - \frac{||x-xc||_2^2}{sig^2} ) \f]
  subroutine SF_init_Vpoint(  &
      N,                      &
      bL,bU,                  &
      SS,NN,                  &
      x1,                     &
      x2,                     &
      x3,                     &
      xc,                     &
      amp,                    &
      sig,                    &
      phi ) bind ( c, name='SF_init_Vpoint' )

    implicit none

    integer(c_int), intent(in)    ::  N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)


    real(c_double), intent(in)     :: x1( bL(1):(N(1)+bU(1)) )
    real(c_double), intent(in)     :: x2( bL(2):(N(2)+bU(2)) )
    real(c_double), intent(in)     :: x3( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(in)     :: xc(3)

    real(c_double), intent(in)     :: sig(3)
    real(c_double), intent(in)     :: amp

    real(c_double),  intent(inout) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    integer(c_int)               ::  i, j, k


    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        do i = SS(1), NN(1)
          phi(i,j,k) = amp*exp( -( (x1(i)-xc(1) )/sig(1) )**2 -( ( x2(j)-xc(2) )/sig(2) )**2 -( ( x3(k)-xc(3) )/sig(3) )**2 )
        end do
      end do
    end do

  end subroutine SF_init_Vpoint



  !> \brief leves scalar field
  !!
  !! \param[in] COMM_CART spatial communicator
  !! \param[in] M number global dofs
  !! \param[in] N number of stored 
  !! \param[in] bl
  !! \param[in] bu
  !! \param[in] SS start indexes
  !! \param[in] NN end indexes
  !! \param[inout] phi scalar field
  subroutine SF_level(    &
      COMM_CART,          &
      M,                  &
      N,                  &
      bL,bU,              &
      SS,NN,              &
      phi ) bind ( c, name='SF_level' )

    implicit none

    integer(c_int), intent(in)    :: COMM_CART

    integer(c_int), intent(in)    :: M

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: SS(3)
    integer(c_int), intent(in)    :: NN(3)

    real(c_double),  intent(out)  :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                :: i, j, k, merror
    real(c_double)                :: pre0, pre0_global


    pre0 = 0.

    do k = SS(3), NN(3)
      do j = SS(2), NN(2)
        !pgi$ unroll = n:8
        do i = SS(1), NN(1)
          pre0 = pre0 + phi(i,j,k)
        end do
      end do
    end do

    call MPI_ALLREDUCE(pre0,pre0_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)

    pre0 = pre0_global/real(M,c_double) ! TEST!!! wegen i4 gefaehrlich!

    phi( SS(1):NN(1), SS(2):NN(2), SS(3):NN(3) ) = phi( SS(1):NN(1), SS(2):NN(2), SS(3):NN(3) ) - pre0

  end subroutine SF_level



  !> \brief sets corner and edges to zero
  !! 
  !! sets undefined pressure on coreners and edges to zero, to not influence
  !! convergence of iterative solvers
  !! \note Siehe dazu auch die korrespondierende Subroutine "handle_corner_rhs"!
  subroutine SF_handle_corner(  &
      n,                        &
      bL,                       &   
      bU,                       &
      bcl,                      &
      bcu,                      &
      phi) bind( c, name='SF_handle_corner' )

    implicit none

    integer(c_int), intent(in)    :: n(3)

    integer(c_int), intent(in)    :: bl(3)
    integer(c_int), intent(in)    :: bu(3)

    integer(c_int), intent(in)    :: bcl(3)
    integer(c_int), intent(in)    :: bcu(3)

    real(c_double), intent(inout) :: phi (bl(1):(n(1)+bu(1)),bl(2):(n(2)+bu(2)),bl(3):(n(3)+bu(3)))


    !if( bcl(1) > 0 .and. bcl(2) > 0 ) phi(1   ,1   ,1:n(3)) = 0. ! test!!! verifizieren ...
    !if( bcl(1) > 0 .and. bcu(2) > 0 ) phi(1   ,n(2),1:n(3)) = 0.
    !if( bcu(1) > 0 .and. bcl(2) > 0 ) phi(n(1),1   ,1:n(3)) = 0.
    !if( bcu(1) > 0 .and. bcu(2) > 0 ) phi(n(1),n(2),1:n(3)) = 0.

    !if( bcl(1) > 0 .and. bcl(3) > 0 ) phi(1   ,1:n(2),1   ) = 0.
    !if( bcl(1) > 0 .and. bcu(3) > 0 ) phi(1   ,1:n(2),n(3)) = 0.
    !if( bcu(1) > 0 .and. bcl(3) > 0 ) phi(n(1),1:n(2),1   ) = 0.
    !if( bcu(1) > 0 .and. bcu(3) > 0 ) phi(n(1),1:n(2),n(3)) = 0.

    !if( bcl(2) > 0 .and. bcl(3) > 0 ) phi(1:n(1),1   ,1   ) = 0.
    !if( bcl(2) > 0 .and. bcu(3) > 0 ) phi(1:n(1),1   ,n(3)) = 0.
    !if( bcu(2) > 0 .and. bcl(3) > 0 ) phi(1:n(1),n(2),1   ) = 0.
    !if( bcu(2) > 0 .and. bcu(3) > 0 ) phi(1:n(1),n(2),n(3)) = 0.


  end subroutine SF_handle_corner



end module cmod_ScalarField
