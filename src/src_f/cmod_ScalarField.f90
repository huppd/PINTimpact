!>  Modul: cmod_ScalarVector
!!
!! impact functions \c Pimpact::ScalarVector e.g. scales, norms ...
!! \author huppd
module cmod_ScalarVector

  use iso_c_binding
  use mpi

  use mod_vars, only: weight

  implicit none
!  public get_norms

  contains

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
  subroutine add(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,phi1,phi2,scalar1,scalar2) bind ( c, name='SV_add' )

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

  real(c_double),  intent(out)   ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(in)   ::  scalar1
  real(c_double),  intent(in)   ::  scalar2


    phi(SS1:NN1,SS2:NN2,SS3:NN3) = scalar1*phi1(SS1:NN1,SS2:NN2,SS3:NN3)+scalar2*phi2(SS1:NN1,SS2:NN2,SS3:NN3)

  end subroutine add


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
  subroutine product_scalar(COMM_CART,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi1,phi2,scalar) bind ( c, name='SV_dot' )

  implicit none

  integer(c_int), intent(in)    ::  COMM_CART

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

  real(c_double),  intent(in)    ::  phi1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)    ::  phi2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(out)   ::  scalar
  real(c_double)                 ::  scalar_global
  integer                       ::  i, j, k
  integer                       ::  merror


  scalar = 0.

  do k = SS3, NN3
     do j = SS2, NN2
!pgi$ unroll = n:8
        do i = SS1, NN1
           scalar = scalar + phi1(i,j,k)*phi2(i,j,k)
        end do
     end do
  end do

  call MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  scalar = scalar_global


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
  !! \param[in] weighting_yes if weighting schould be used, using the \c weights from \c mod_vars
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  !! \param[out] normInf gets the infinity norm of phi
  !! \param[out] normTwo get the two norm of phi
  !! \todo weight (easydirty fix: include module) (good persisting fix: move to pimpact)
  ! TEST!!! N1, N2, N3 werden ebenfalls uebergeben ...
  subroutine get_norms(comm,N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,weighting_yes,inf_yes,two_yes,normInf,normTwo) bind (c, name='SV_compNorm')
  implicit none

  integer(c_int), intent(in)    ::  comm

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

  real(c_double),  intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  logical(c_bool), intent(in)   ::  weighting_yes
  logical(c_bool), intent(in)   ::  inf_yes
  logical(c_bool), intent(in)   ::  two_yes

  real(c_double), optional, intent(out) ::  normInf
  real(c_double), optional, intent(out) ::  normTwo

  real(c_double)                        ::  normInf_global, normTwo_global
  integer                              ::  i, j, k

  integer                              :: merror


  if (inf_yes .and. two_yes) then

     normInf = 0.
     normTwo = 0.

     if( weighting_yes ) then

        do k = SS3, NN3
           do j = SS2, NN2
!pgi$ unroll = n:8
              do i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              end do
           end do
        end do

     else

        do k = SS3, NN3
           do j = SS2, NN2
!pgi$ unroll = n:8
              do i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
                 normTwo = normTwo + phi(i,j,k)**2
              end do
           end do
        end do

     end if

     ! Lassen sich wegen MPI_SUM / MPI_MAX nicht zusammenlegen:
     call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,comm,merror)
     call MPI_ALLREDUCE(normTwo,normTwo_global,1,MPI_REAL8,MPI_SUM,comm,merror)
     normInf = normInf_global
     normTwo = normTwo_global

  else if (inf_yes) then

     normInf = 0.

     if( weighting_yes ) then

        do k = SS3, NN3
           do j = SS2, NN2
!pgi$ unroll = n:8
              do i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)*weight(i,j,k)),normInf)
              end do
           end do
        end do

     else

        do k = SS3, NN3
           do j = SS2, NN2
!pgi$ unroll = n:8
              do i = SS1, NN1
                 normInf = MAX(ABS(phi(i,j,k)),normInf)
              end do
           end do
        end do

     end if

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


  subroutine scale(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,scalar) bind ( c, name='SV_scale' )

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

    phi(SS1:NN1,SS2:NN2,SS3:NN3) = scalar*phi(SS1:NN1,SS2:NN2,SS3:NN3)

  end subroutine scale


  subroutine random(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi) bind ( c, name='SV_random' )

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

  end subroutine random


  subroutine init(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi,scalar) bind ( c, name='SV_init' )

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


  subroutine print(N1,N2,N3,SS1,SS2,SS3,NN1,NN2,NN3,b1L,b2L,b3L,b1U,b2U,b3U,phi) bind ( c, name='SV_print' )

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

  real(c_double),  intent(in)    ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
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
