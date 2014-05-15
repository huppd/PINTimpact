!>  Modul: cmod_ScalarVector
!!
!! impact functions \c Pimpact::ScalarVector e.g. scales, norms ...
!! \author huppd
module cmod_ScalarVector

  use iso_c_binding
  use mpi


  implicit none
!  public get_norms

  contains

  !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
  !!
  !! \param[in] N ammount of local elements in 1-direction
  !! \param[in] bL start index of storage in 1-direction
  !! \param[in] bU end offset of storage in 1-direction
  !! \param[in] SS start index in 1-direction
  !! \param[in] NN end index in 1-direction
  !! \param[in] phi first vector from which is taken the product
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  subroutine SF_add(    &
    N,                  &
    bL,bU,              &
    SS,NN,              &
    phi,phi1,phi2,      &
    scalar1,scalar2 ) bind ( c, name='SF_add' )

  implicit none

    integer(c_int), intent(in)     :: N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)


    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(in )   :: phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(in )   :: phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

!    real(c_double), intent(out)   ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!    real(c_double), intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
!    real(c_double), intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

    real(c_double),  intent(in)   ::  scalar1
    real(c_double),  intent(in)   ::  scalar2

!    integer                       :: i,j,k


    phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar1*phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))+scalar2*phi2(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

!  do k = SS3, NN3
!    do j = SS2, NN2
!!pgi$ unroll = n:8
!      do i = SS1, NN1
!        phi(i,j,k) = scalar1*phi1(i,j,k)+scalar2*phi2(i,j,k)
!      end do
!    end do
!  end do

  end subroutine SF_add


  !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
  !!
  !! \param[in] N ammount of local elements in 1-direction
  !! \param[in] bL start index of storage in 1-direction
  !! \param[in] bU end offset of storage in 1-direction
  !! \param[in] SS start index in 1-direction
  !! \param[in] NN end index in 1-direction
  !! \param[in] phi first vector from which is taken the product
  !! \param[in] phi1 first vector from which is taken the product
  subroutine SF_add2(   &
    N,                  &
    bL,bU,              &
    SS,NN,              &
    phi,phi1,           &
    scalar1,scalar2 ) bind ( c, name='SF_add2' )

  implicit none

    integer(c_int), intent(in)     :: N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SS(3)
    integer(c_int), intent(in)     :: NN(3)


    real(c_double),  intent(out)   :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(in )   :: phi1(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
!    real(c_double),  intent(in )   :: phi2(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))


    real(c_double),  intent(in)   ::  scalar1
    real(c_double),  intent(in)   ::  scalar2


    phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3)) = scalar1*phi(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))+scalar2*phi1(SS(1):NN(1),SS(2):NN(2),SS(3):NN(3))

!  do k = SS3, NN3
!    do j = SS2, NN2
!!pgi$ unroll = n:8
!      do i = SS1, NN1
!        phi(i,j,k) = scalar1*phi1(i,j,k)+scalar2*phi2(i,j,k)
!      end do
!    end do
!  end do

  end subroutine SF_add2


  !> \brief copys absolute value from \c phi1 in \phi
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
  subroutine SF_abs(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,phi1) bind ( c, name='SF_abs' )

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                       ::  i, j, k

  do k = SS3, NN3
    do j = SS2, NN2
!pgi$ unroll = n:8
      do i = SS1, NN1
        phi(i,j,k) = abs( phi1(i,j,k) )
      end do
    end do
  end do

  end subroutine SF_abs


  !> \brief copys reciprocal value from \c phi1 in \phi
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  subroutine SF_reciprocal(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,phi1) bind ( c, name='SF_reciprocal' )

  implicit none


  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(in   ) ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                       ::  i, j, k

    do k = SS3, NN3
      do j = SS2, NN2
!pgi$ unroll = n:8
        do i = SS1, NN1
          if( phi1(i,j,k) .EQ. 0. ) then
            phi(i,j,k) = 0.
          else
            phi(i,j,k) = 1./phi1(i,j,k)
          end if
        end do
      end do
    end do

  end subroutine SF_reciprocal


  !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
  !! \param[in] COMM_CART communicator belonging to vector
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  !! \test if comm_cart is neccessary, or if comm is enough or even better
  subroutine product_scalar(    &
    N1,N2,N3,                   &
    SS1,SS2,SS3,                &
    NN1,NN2,NN3,                &
    b1L,b2L,b3L,                &
    b1U,b2U,b3U,                &
    phi1,phi2,scalar) bind ( c, name='SF_dot' )

  implicit none

!  integer(c_int), intent(in)    ::  COMM_CART

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double), intent(out)   ::  scalar
  integer                       ::  i, j, k


  scalar = 0.

  do k = SS3, NN3
     do j = SS2, NN2
!pgi$ unroll = n:8
        do i = SS1, NN1
           scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
        end do
     end do
  end do


  end subroutine product_scalar



  !> \brief computes two or infinity norm( get is misleading)
  !! \param[in] comm communicator belonging to vector
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi vector, from which the norm is taken
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  !! \param[out] normInf gets the infinity norm of phi
  !! \param[out] normTwo get the two norm of phi
  ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
  subroutine get_norms(comm,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,inf_yes,two_yes,normInf,normTwo) bind (c, name='SF_compNorm')
  implicit none

  integer(c_int), intent(in)  ::  comm

  integer(c_int), intent(in)  ::  N1
  integer(c_int), intent(in)  ::  N2
  integer(c_int), intent(in)  ::  N3

  integer(c_int), intent(in)  ::  SS1
  integer(c_int), intent(in)  ::  SS2
  integer(c_int), intent(in)  ::  SS3

  integer(c_int), intent(in)  ::  NN1
  integer(c_int), intent(in)  ::  NN2
  integer(c_int), intent(in)  ::  NN3

  integer(c_int), intent(in)  ::  b1L
  integer(c_int), intent(in)  ::  b2L
  integer(c_int), intent(in)  ::  b3L

  integer(c_int), intent(in)  ::  b1U
  integer(c_int), intent(in)  ::  b2U
  integer(c_int), intent(in)  ::  b3U

  real(c_double), intent(in)  ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  logical(c_bool),intent(in)  ::  inf_yes
  logical(c_bool),intent(in)  ::  two_yes

  real(c_double), intent(out) ::  normInf
  real(c_double), intent(out) ::  normTwo

  real(c_double)              ::  normInf_global, normTwo_global
  integer                     ::  i, j, k

  integer                     :: merror

    if (inf_yes .and. two_yes) then

        normInf = 0.
        normTwo = 0.

        do k = SS3, NN3
            do j = SS2, NN2
!pgi$ unroll = n:8
                do i = SS1, NN1
                    normInf = MAX(ABS(phi(i,j,k)),normInf)
                    normTwo = normTwo + phi(i,j,k)**2
                end do
            end do
        end do

        ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
        call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm,merror)
        call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm,merror)

        normInf = normInf_global
        normTwo = normTwo_global

    else if (inf_yes) then

        normInf = 0.

        do k = SS3, NN3
           do j = SS2, NN2
!pgi$ unroll = n:8
              do i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
              end do
           end do
        end do

        call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
        normInf = normInf_global

    else if (two_yes) then

        normTwo = 0.

        do k = SS3, NN3
            do j = SS2, NN2
!pgi$ unroll = n:8
                do i = SS1, NN1
                    normTwo = normTwo + phi(i,j,k)**2
                end do
            end do
        end do

        call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm,merror)
        normTwo = normTwo_global

    end if

  end subroutine get_norms


  !> \brief computes two or infinity norm( get is misleading)
  !! \param[in] comm communicator belonging to vector
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi vector, from which the norm is taken
  !! \param[in] weights weights
  !! \todo weight (easydirty fix: include module) (good persisting fix: move to pimpact)
  ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
  subroutine get_weighted_norms(comm,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,weights,norm) bind (c, name='SF_weightedNorm')
  implicit none

  integer(c_int), intent(in)  ::  comm

  integer(c_int), intent(in)  ::  N1
  integer(c_int), intent(in)  ::  N2
  integer(c_int), intent(in)  ::  N3

  integer(c_int), intent(in)  ::  SS1
  integer(c_int), intent(in)  ::  SS2
  integer(c_int), intent(in)  ::  SS3

  integer(c_int), intent(in)  ::  NN1
  integer(c_int), intent(in)  ::  NN2
  integer(c_int), intent(in)  ::  NN3

  integer(c_int), intent(in)  ::  b1L
  integer(c_int), intent(in)  ::  b2L
  integer(c_int), intent(in)  ::  b3L

  integer(c_int), intent(in)  ::  b1U
  integer(c_int), intent(in)  ::  b2U
  integer(c_int), intent(in)  ::  b3U

  real(c_double), intent(in)  ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(in)  ::  weights(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))


  real(c_double), intent(out) ::  norm

  real(c_double)              ::  norm_global
  integer                     ::  i, j, k

  integer                     :: merror


    norm = 0.

    do k = SS3, NN3
        do j = SS2, NN2
!pgi$ unroll = n:8
            do i = SS1, NN1
                norm = norm + ( weights(i,j,k)*phi(i,j,k) )**2
            end do
        end do
     end do


     call MPI_ALLREDUCE(norm,norm_global,1,MPI_REAL8,MPI_SUM,comm,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden müsste ...
     norm = norm_global

  end subroutine get_weighted_norms


  subroutine scale(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,scalar) bind ( c, name='SF_scale' )

  implicit none

!  integer(c_int), intent(in)    ::  comm

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))


  real(c_double), intent(in) ::  scalar

    phi(SS1:NN1,SS2:NN2,SS3:NN3) = scalar*phi(SS1:NN1,SS2:NN2,SS3:NN3)

  end subroutine scale


  !> \brief scales \c phi with \phi1.
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] SS1 start index in 1-direction
  !! \param[in] SS2 start index in 1-direction
  !! \param[in] SS3 start index in 1-direction
  !! \param[in] NN1 end index in 1-direction
  !! \param[in] NN2 end index in 2-direction
  !! \param[in] NN3 end index in 3-direction
  !! \param[in] b1L start index of storage in 1-direction
  !! \param[in] b2L start index of storage in 2-direction
  !! \param[in] b3L start index of storage in 3-direction
  !! \param[in] b1U end offset of storage in 1-direction
  !! \param[in] b2U end offset of storage in 2-direction
  !! \param[in] b3U end offset of storage in 3-direction
  !! \param[in] phi1 first vector from which is taken the product
  !! \param[in] phi2 second vector from which is taken the product
  !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  subroutine SF_scale2(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,phi1) bind ( c, name='SF_scale2' )

  implicit none


  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                       ::  i, j, k

    do k = SS3, NN3
      do j = SS2, NN2
!pgi$ unroll = n:8
        do i = SS1, NN1
          if( phi1(i,j,k) .EQ. 0. ) then
            phi(i,j,k) = 0.
          else
            phi(i,j,k) = phi(i,j,k)*phi1(i,j,k)
          end if
        end do
      end do
    end do

  end subroutine SF_scale2


  subroutine random(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi) bind ( c, name='SF_random' )

  implicit none

!  integer(c_int), intent(in)    ::  comm

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double),  intent(inout)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))


    call Random_number(phi(SS1:NN1,SS2:NN2,SS3:NN3))
    phi(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3) - 0.5

  end subroutine random


  subroutine init(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,scalar) bind ( c, name='SF_init' )

  implicit none

!  integer(c_int), intent(in)    ::  comm

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double),  intent(inout)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))


  real(c_double), intent(in) ::  scalar

    phi(SS1:NN1,SS2:NN2,SS3:NN3) = scalar

  end subroutine init


  subroutine print(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi) bind ( c, name='SF_print' )

  implicit none

!  integer(c_int), intent(in)    ::  comm

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

  integer(c_int), intent(in)    ::  SS1
  integer(c_int), intent(in)    ::  SS2
  integer(c_int), intent(in)    ::  SS3

  integer(c_int), intent(in)    ::  NN1
  integer(c_int), intent(in)    ::  NN2
  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double),  intent(inout)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  integer                       ::  i,j,k


    do k = SS3, NN3
      do j = SS2, NN2
!pgi$ unroll = n:8
        do i = SS1, NN1
          write(*,*) i,j,k,phi(i,j,k)
        end do
      end do
    end do

  end subroutine print


end module cmod_ScalarVector
