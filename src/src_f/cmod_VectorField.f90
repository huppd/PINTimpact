!>  Modul: cmod_VectorField
!!
!! impact functions \c Pimpact::VectorField e.g. scales, norms ...
!! \author huppd
module cmod_VectorField

  use iso_c_binding
  use mpi

  use mod_vars

  implicit none
!  public get_norms

  contains



  !> \brief computes scalar product of two scalar fields(neccessary condition belong to same sVS)
  !! \param[in] COMM_CART communicator belonging to vector
  !! \param[in] N1 ammount of local elements in 1-direction
  !! \param[in] N2 ammount of local elements in 2-direction
  !! \param[in] N3 ammount of local elements in 3-direction
  !! \param[in] S1U start index in 1-direction
  !! \param[in] S2U start index in 1-direction
  !! \param[in] S3U start index in 1-direction
  !! \param[in] N1U end index in 1-direction
  !! \param[in] N2U end index in 2-direction
  !! \param[in] N3U end index in 3-direction
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
  subroutine product_scalar_vel( COMM_CART, dimens, N1,N2,N3, S1U,S2U,S3U, N1U,N2U,N3U, S1V,S2V,S3V, N1V,N2V,N3V, S1W,S2W,S3W, N1W,N2W,N3W, b1L,b2L,b3L, b1U,b2U,b3U, phi1U,phi1V,phi1W, phi2U,phi2V,phi2W, scalar ) bind ( c, name='VF_dot' )

  implicit none

  integer(c_int), intent(in)    ::  COMM_CART

  integer(c_int), intent(in)    ::  dimens

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

  real(c_double),  intent(in)   ::  phi1U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phi1V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phi1W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(in)   ::  phi2U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phi2V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phi2W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(out)  ::  scalar
  real(c_double)                ::  scalar_global
  integer                       ::  i, j, k
  integer                       ::  merror


  scalar = 0.

  do k = S3U, N3U
    do j = S2U, N2U
!pgi$ unroll = n:8
      do i = S1U, N1U
        scalar = scalar + phi1U(i,j,k)*phi2U(i,j,k)
      end do
   end do
  end do

  do k = S3V, N3V
    do j = S2V, N2V
!pgi$ unroll = n:8
      do i = S1V, N1V
        scalar = scalar + phi1V(i,j,k)*phi2V(i,j,k)
      end do
    end do
  end do

  if (dimens == 3) then
    do k = S3W, N3W
      do j = S2W, N2W
!pgi$ unroll = n:8
        do i = S1W, N1W
          scalar = scalar + phi1W(i,j,k)*phi2W(i,j,k)
        end do
      end do
    end do
  end if

  call MPI_ALLREDUCE(scalar,scalar_global,1,MPI_REAL8,MPI_SUM,COMM_CART,merror)
  scalar = scalar_global

  end subroutine product_scalar_vel



  !> brief computes two or infinity norm
  !! \todo add all parameters to make it more flexible
  !!
  !! \param[in] phi velocity vector, from which the norm is taken
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  !! \param[out] normInf gets the infinity norm of phi
  !! \param[out] normTwo get the two norm of phi
  subroutine get_norms_vel( COMM_CART, dimens, N1,N2,N3, S1U,S2U,S3U, N1U,N2U,N3U, S1V,S2V,S3V, N1V,N2V,N3V, S1W,S2W,S3W, N1W,N2W,N3W, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, inf_yes,two_yes, normInf,normTwo ) bind (c,name='VF_compNorm')

  implicit none

  integer(c_int), intent(in)    ::  COMM_CART

  integer(c_int), intent(in)    ::  dimens

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

  real(c_double),  intent(in)    ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)    ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)    ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  logical(c_bool), intent(in)      ::  inf_yes
  logical(c_bool), intent(in)      ::  two_yes

  real(c_double),  intent(out) ::  normInf
  real(c_double),  intent(out) ::  normTwo

  real(c_double)                   ::  normInf_global, normTwo_global
  integer(c_int)                   ::  i, j, k


  if (inf_yes .and. two_yes) then

     normInf = 0.
     normTwo = 0.

     do k = S3U, N3U
       do j = S2U, N2U
!pgi$ unroll = n:8
         do i = S1U, N1U
           normInf = MAX(ABS(phiU(i,j,k)),normInf)
           normTwo = normTwo + phiU(i,j,k)**2
         end do
       end do
     end do
     do k = S3V, N3V
       do j = S2V, N2V
!pgi$ unroll = n:8
         do i = S1V, N1V
           normInf = MAX(ABS(phiV(i,j,k)),normInf)
           normTwo = normTwo + phiV(i,j,k)**2
         end do
       end do
     end do

     if (dimens == 3) then
       do k = S3W, N3W
         do j = S2W, N2W
  !pgi$ unroll = n:8
           do i = S1W, N1W
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

     do k = S3U, N3U
       do j = S2U, N2U
!pgi$ unroll = n:8
         do i = S1U, N1U
           normInf = MAX(ABS(phiU(i,j,k)),normInf)
         end do
       end do
     end do
     do k = S3V, N3V
       do j = S2V, N2V
!pgi$ unroll = n:8
         do i = S1W, N1W
           normInf = MAX(ABS(phiV(i,j,k)),normInf)
         end do
       end do
     end do

     if (dimens == 3) then
       do k = S3W, N3W
         do j = S2W, N2W
  !pgi$ unroll = n:8
           do i = S1W, N1W
             normInf = MAX(ABS(phiW(i,j,k)),normInf)
           end do
         end do
       end do
     end if


    call MPI_ALLREDUCE(normInf,normInf_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror) ! MPI_REDUCE bringt nichts, weil exit_yes dann mit MPI_BCAST verteilt werden m�sste ...
    normInf = normInf_global

  else if (two_yes) then

     normTwo = 0.

     do k = S3U, N3U
       do j = S2U, N2U
!pgi$ unroll = n:8
         do i = S1U, N1U
           normTwo = normTwo + phiU(i,j,k)**2
         end do
       end do
     end do
     do k = S3V, N3V
       do j = S2V, N2V
!pgi$ unroll = n:8
         do i = S1V, N1V
           normTwo = normTwo + phiV(i,j,k)**2
         end do
       end do
     end do

     if (dimens == 3) then
       do k = S3W, N3W
         do j = S2W, N2W
  !pgi$ unroll = n:8
           do i = S1W, N1W
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
  !! \todo add all parameters to make it more flexible
  !!
  !! \param[in] phi velocity vector, from which the norm is taken
  !! \param[in] inf_yes if true infinity norm is computed
  !! \param[in] two_yes if trhue two norm is computed
  !! \param[out] normInf gets the infinity norm of phi
  !! \param[out] normTwo get the two norm of phi
  subroutine get_weighted_norm_vel( COMM_CART, dimens, N1,N2,N3, S1U,S2U,S3U, N1U,N2U,N3U, S1V,S2V,S3V, N1V,N2V,N3V, S1W,S2W,S3W, N1W,N2W,N3W, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, wU,wV,wW, norm ) bind (c,name='VF_weightedNorm')

  implicit none

  integer(c_int), intent(in)    ::  COMM_CART

  integer(c_int), intent(in)    ::  dimens

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

  real(c_double),  intent(in)   ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(in)   ::  wU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  wV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(in)   ::  wW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))


  real(c_double),  intent(out)   ::  norm

  real(c_double)                 ::  norm_global
  integer                        ::  i, j, k


    norm = 0.

    do k = S3U, N3U
        do j = S2U, N2U
!pgi$ unroll = n:8
            do i = S1U, N1U
                norm = norm + (wU(i,j,k)*phiU(i,j,k))**2
            end do
        end do
    end do

    do k = S3V, N3V
        do j = S2V, N2V
!pgi$ unroll = n:8
            do i = S1V, N1V
                norm = norm + (wV(i,j,k)*phiV(i,j,k))**2
            end do
        end do
    end do

    if (dimens == 3) then
        do k = S3W, N3W
            do j = S2W, N2W
!pgi$ unroll = n:8
                do i = S1W, N1W
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
  subroutine cinit_2DPoiseuilleX( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW ) bind ( c, name='VF_init_2DPoiseuilleX' )

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3
!
  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                ::  i, j, k

  !--- initial conditions for velocity ---
  ! note: - cf. sketch in file "usr_geometry.f90"
  !
  !         grid points in the domain and on the boundary
  !         |         |         |     velocity component
  !         |         |         |     |
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
!  phiU = 0.
!  phiV = 0.
!  phiW = 0.
  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           phiU(i,j,k) = x2p(j)*( L2 - x2p(j) )*4/L2/L2
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = 0
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPoiseuilleX


  !> \brief init vector field with 2d Poiseuille flow in y-direction
  subroutine cinit_2DPoiseuilleY( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW ) bind ( c, name='VF_init_2DPoiseuilleY' )

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3
!
  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                ::  i, j, k

  !--- initial conditions for velocity ---
  ! note: - cf. sketch in file "usr_geometry.f90"
  !
  !         grid points in the domain and on the boundary
  !         |         |         |     velocity component
  !         |         |         |     |
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
          phiU(i,j,k) = 0.
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = x1p(i)*( L1 - x1p(i) )*4/L1/L1
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPoiseuilleY


  !> \brief init vector field with a zero flow
  subroutine cinit_zero( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW ) bind ( c, name='VF_init_Zero' )

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  integer                ::  i, j, k

  !--- initial conditions for velocity ---
  ! note: - cf. sketch in file "usr_geometry.f90"
  !
  !         grid points in the domain and on the boundary
  !         |         |         |     velocity component
  !         |         |         |     |
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
          phiU(i,j,k) = 0.0
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = 0.0
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_zero


  !> \brief init vector field with 2d pulsatile flow in x-direction
  !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
  !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
  !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
  subroutine cinit_2DPulsatileXC( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, re_, om, px ) bind ( c, name='VF_init_2DPulsatileXC' )

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

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
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
!  phiU = 0.
!  phiV = 0.
!  phiW = 0.
  pi = 4.*atan(1.)
  Lh  = L2/2.
  mu = sqrt( om*re_/2. )*Lh
  c  = px/( om*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           nu = sqrt( om*re_/2. )*( x2p(j)-Lh )
           phiU(i,j,k) = -c*( -cos(nu)*cosh(nu)*sin(mu)*sinh(mu) +sin(nu)*sinh(nu)*cos(mu)*cosh(mu) )
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = 0
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPulsatileXC


  !> \brief init vector field with 2d pulsatile flow in x-direction
  !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
  !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
  !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
  subroutine cinit_2DPulsatileYC( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileYC' )
  ! (basic subroutine)

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

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
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
!  phiU = 0.
!  phiV = 0.
!  phiW = 0.
  pi = 4.*atan(1.)
  Lh  = L1/2.
  mu = sqrt( om_*re_/2. )*Lh
  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           phiU(i,j,k) = 0
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           nu = sqrt( om_*re_/2. )*( x1p(i)-Lh )
           phiV(i,j,k) = -c*( -cos(nu)*cosh(nu)*sin(mu)*sinh(mu) +sin(nu)*sinh(nu)*cos(mu)*cosh(mu) )
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPulsatileYC


  !> \brief init vector field with 2d pulsatile flow in x-direction
  !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
  !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
  !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
  subroutine cinit_2DPulsatileXS( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileXS' )
  ! (basic subroutine)

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

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
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
!  phiU = 0.
!  phiV = 0.
!  phiW = 0.
  pi = 4.*atan(1.)
  Lh  = L2/2.
  mu = sqrt( om_*re_/2. )*Lh
  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           nu = sqrt( om_*Re_/2. )*( x2p(j)-Lh )
           phiU(i,j,k) = -c*( cos(nu)*cosh(nu)*cos(mu)*cosh(mu) +sin(nu)*sinh(nu)*sin(mu)*sinh(mu) ) + px/om_
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = 0
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPulsatileXS


  !> \brief init vector field with 2d pulsatile flow in x-direction
  !! \f[ u(y,t) = \hat{u}^+ \exp(i \omega t) + \hat{u}^- \exp(- i \omega t) \f]
  !! \f[ \hat{u}^+(y) = c^+ \left( \exp(+ \lambda_1 y ) + \exp(-\lambda\right) + \frac{p_x i}{\omega \f]
  !! \f[ \hat{u}^-(y) = c^- \left( \exp(+ \lambda_{-1} y ) \right) + \frac{p_x i}{\omega \f]
  subroutine cinit_2DPulsatileYS( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW, re_, om_, px ) bind ( c, name='VF_init_2DPulsatileYS' )
  ! (basic subroutine)

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

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
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
!  phiU = 0.
!  phiV = 0.
!  phiW = 0.
  pi = 4.*atan(1.)
  Lh  = L1/2.
  mu = sqrt( om_*re_/2. )*Lh
  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           phiU(i,j,k) = 0
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           nu = sqrt( om_*re_/2. )*( x1p(i)-Lh )
           phiV(i,j,k) = -c*( cos(nu)*cosh(nu)*cos(mu)*cosh(mu) +sin(nu)*sinh(nu)*sin(mu)*sinh(mu) ) + px/om_
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_2DPulsatileYS



 subroutine cinit_Streaming( N1,N2,N3, b1L,b2L,b3L, b1U,b2U,b3U, phiU,phiV,phiW ) bind ( c, name='VF_init_Streaming' )
  ! (basic subroutine)

  implicit none

  integer(c_int), intent(in)    ::  N1
  integer(c_int), intent(in)    ::  N2
  integer(c_int), intent(in)    ::  N3

!  integer(c_int), intent(in)    ::  SS1
!  integer(c_int), intent(in)    ::  SS2
!  integer(c_int), intent(in)    ::  SS3
!
!  integer(c_int), intent(in)    ::  NN1
!  integer(c_int), intent(in)    ::  NN2
!  integer(c_int), intent(in)    ::  NN3

  integer(c_int), intent(in)    ::  b1L
  integer(c_int), intent(in)    ::  b2L
  integer(c_int), intent(in)    ::  b3L

  integer(c_int), intent(in)    ::  b1U
  integer(c_int), intent(in)    ::  b2U
  integer(c_int), intent(in)    ::  b3U

  real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

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
  ! vel(S11B:N11B,S21B:N21B,S31B:N31B,1)
  ! vel(S12B:N12B,S22B:N22B,S32B:N32B,2)
  ! vel(S13B:N13B,S23B:N23B,S33B:N33B,3)
  !
  pi = 4.*atan(1.)
!  Lh  = L1/2.
!  mu = sqrt( om_*re_/2. )*Lh
!  c  = px/( om_*(cos(mu)**2*cosh(mu)**2 + sin(mu)**2*sinh(mu)**2) )

  do k = S31B, N31B
     do j = S21B, N21B
        do i = S11B, N11B
           phiU(i,j,k) = 0
        end do
     end do
  end do

  do k = S32B, N32B
     do j = S22B, N22B
        do i = S12B, N12B
           phiV(i,j,k) = sin( 2.*pi*x1p(i)/L1 )
        end do
     end do
  end do

  do k = S33B, N33B
     do j = S23B, N23B
        do i = S13B, N13B
           phiW(i,j,k) = 0.0
        end do
     end do
  end do

  end subroutine cinit_Streaming


end module cmod_VectorField
