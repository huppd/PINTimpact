!***********************************************************************************************
!* IMPACT                                                                                      *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)       *
!* Mai 2005 - Dec 2011                                                                         *
!***********************************************************************************************



!> \brief module
module cmod_DivGrad2ndOOp

  use iso_c_binding

  implicit none

contains

  subroutine Op_getCDG(   &
      dimens,             &
      M,                  &
      N,                  &
      BL,                 &
      BU,                 &
      BCL,                &
      BCU,                &
      y1u,                &
      y2v,                &
      y3w,                &
      x1p,                &
      x2p,                &
      x3p,                &
      x1u,                &
      x2v,                &
      x3w,                &
      cdg1,               &
      cdg2,               &
      cdg3 ) bind( c, name='Op_getCDG' )

    implicit none

    integer(c_int), intent(in)    :: dimens

    integer(c_int), intent(in)    :: M(3)

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: BL(3)
    integer(c_int), intent(in)    :: BU(3)

    integer(c_int), intent(in)    :: BCL(3)
    integer(c_int), intent(in)    :: BCU(3)

    real(c_double), intent(in)    :: y1u(0:M(1))
    real(c_double), intent(in)    :: y2v(0:M(2))
    real(c_double), intent(in)    :: y3w(0:M(3))

    real(c_double), intent(in)    :: x1p(bl(1):N(1)+bu(1))
    real(c_double), intent(in)    :: x2p(bl(2):N(2)+bu(2))
    real(c_double), intent(in)    :: x3p(bl(3):N(3)+bu(3))

    real(c_double), intent(in)    :: x1u(bl(1):N(1)+bu(1))
    real(c_double), intent(in)    :: x2v(bl(2):N(2)+bu(2))
    real(c_double), intent(in)    :: x3w(bl(3):N(3)+bu(3))

    real(c_double), intent(inout) :: cdg1(-1:1,1:N(1))
    real(c_double), intent(inout) :: cdg2(-1:1,1:N(2))
    real(c_double), intent(inout) :: cdg3(-1:1,1:N(3))

    real(c_double)                :: cGp1R( 0:1,0:N(1))
    real(c_double)                :: cGp2R( 0:1,0:N(2))
    real(c_double)                :: cGp3R( 0:1,0:N(3))

    real(c_double)                :: cDu1R(-1:0,1:N(1))
    real(c_double)                :: cDv2R(-1:0,1:N(2))
    real(c_double)                :: cDw3R(-1:0,1:N(3))

    integer(c_int)                :: i, j, k
    integer(c_int)                :: jmax, kmax

    cdg1 = 0.
    cdg2 = 0.
    cdg3 = 0.

    cGp1R = 0.
    cDu1R = 0.

    !---------------------------------------------------------------------------!
    ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
    !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
    !---------------------------------------------------------------------------!

    do i = 1, N(1)-1
      cGp1R(0,i) = -1./( x1p(i+1) - x1p(i) )
      cGp1R(1,i) =  1./( x1p(i+1) - x1p(i) )
    end do

    do i = 1, N(1)
      cDu1R(-1,i) = -1./( x1u(i) - x1u(i-1) )
      cDu1R( 0,i) =  1./( x1u(i) - x1u(i-1) )
    end do

    !if (BC(1,1,g) > 0) then
    if( BCL(1) > 0 ) then
      cGp1R(:,0) = 0. !???is set anyway to zero???
      cDu1R(:,1) = 0. 
      !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
      cDu1R(0,1) = 1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
    else
      cGp1R(0,0) = -1./(x1p(1) - x1p(0))
      cGp1R(1,0) =  1./(x1p(1) - x1p(0))
    end if

    !if (BC(2,1,g) > 0) then
    if (BCU(1) > 0) then
      cGp1R(:,N(1)) =  0.
      cDu1R(:,N(1)) =  0.
      !cDu1R(-1,N(1)) = -1./(x1u (N(1)      ) - x1u (N(1)-1  )) ! TEST!!! Das geht in die Hose ...
      cDu1R(-1,N(1)) = -1./( y1u(M(1)) - y1u(M(1)-1) )
    else
      cGp1R( 0,N(1)) = -1./( x1p(N(1)+1) - x1p(N(1)) )
      cGp1R( 1,N(1)) =  1./( x1p(N(1)+1) - x1p(N(1)) )
    end if
    !--------------------------------------------------------------------------------------------
    do i = 1, N(1)
      cdg1(-1,i) = cDu1R(-1,i)*cGp1R(0,i-1)
      cdg1( 0,i) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
      cdg1( 1,i) =                            cDu1R(0,i)*cGp1R(1,i)
    end do

    ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das
    ! Interpolationspolynom am Rand darstellt
    if( BCL(1) > 0 ) cdg1( :,1   ) = 2.*cdg1(:,1   ) ! TEST!!! Ist das wirklich
    ! so optimal?
    if( BCU(1) > 0 ) cdg1( :,N(1)) = 2.*cdg1(:,N(1))

    if (BCL(1) == -2) then
      cdg1( 1,1) = cdg1( 1, 1 ) + cdg1( -1, 1 )
      cdg1(-1,1) = 0.
    end if

    if (BCU(2) == -2) then
      cdg1(-1,N(1)) = cdg1( -1, N(1) ) + cdg1( 1, N(1) )
      cdg1( 1,N(1)) = 0.
    end if

    !============================================================================================   !jmax = NN(2,g)
    jmax = N(2)

    cDv2R = 0.
    cGp2R = 0.

    do j = 1, jmax
      cDv2R(-1,j) = -1./( x2v(j) - x2v(j-1) )
      cDv2R( 0,j) =  1./( x2v(j) - x2v(j-1) )
    end do

    do j = 1, jmax-1
      cGp2R( 0,j) = -1./( x2p(j+1) - x2p(j) )
      cGp2R( 1,j) =  1./( x2p(j+1) - x2p(j) )
    end do

    !if (BC(1,2,g) > 0) then
    if (BCL(2) > 0) then
      cGp2R( :,0   ) =  0. !???is set anyway to zero???
      cDv2R( :,1   ) =  0.
      !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
      cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
    else
      cGp2R( 0,0   ) = -1./( x2p(1) - x2p(0) )
      cGp2R( 1,0   ) =  1./( x2p(1) - x2p(0) )
    end if

    !if (BC(2,2,g) > 0) then
    if (BCU(2) > 0) then
      cGp2R( :,jmax) =  0.
      cDv2R( :,jmax) =  0.
      !cDv2R(-1,jmax) = -1./( x2v (N(2)      ) - x2v (N(2)-1  ) ) ! TEST!!! Das geht in die Hose ...
      cDv2R(-1,jmax) = -1./( y2v (M(2)      ) - y2v (M(2)-1  ) )
    else
      cGp2R( 0,jmax) = -1./(x2p(jmax+1) - x2p(jmax) )
      cGp2R( 1,jmax) =  1./(x2p(jmax+1) - x2p(jmax) )
    end if
    !--------------------------------------------------------------------------------------------
    do j = 1, jmax
      cdg2(-1,j) = cDv2R(-1,j)*cGp2R(0,j-1)
      cdg2( 0,j) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
      cdg2( 1,j) =                            cDv2R(0,j)*cGp2R(1,j)
    end do

    ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
    !if (BC(1,2,g) > 0) cdg2( :,1   ) = 2.*cdg2(:,1   )
    !if (BC(2,2,g) > 0) cdg2( :,jmax) = 2.*cdg2(:,jmax)
    if (BCL(2) > 0) cdg2( :,1   ) = 2.*cdg2(:,1   )
    if (BCU(2) > 0) cdg2( :,jmax) = 2.*cdg2(:,jmax)

    !if (BC(1,2,g) == -2) then
    if (BCL(2) == -2) then
      cdg2( 1,1) = cdg2( 1, 1 ) + cdg2(-1, 1 )
      cdg2(-1,1) = 0.
    end if
    !if (BC(2,2,g) == -2) then
    if (BCU(2) == -2) then
      cdg2(-1,jmax) = cdg2(-1,jmax) + cdg2( 1,jmax)
      cdg2( 1,jmax) = 0.
    end if
    !============================================================================================
    if (dimens == 3) then

      kmax = N(3)

      cDw3R = 0.
      cGp3R = 0.

      do k = 1, kmax
        cDw3R(-1,k) = -1./( x3w(k) - x3w(k-1) )
        cDw3R( 0,k) =  1./( x3w(k) - x3w(k-1) )
      end do

      do k = 1, kmax-1
        cGp3R( 0,k) = -1./( x3p(k+1) - x3p(k) )
        cGp3R( 1,k) =  1./( x3p(k+1) - x3p(k) )
      end do

      !if (BC(1,3,g) > 0) then
      if (BCL(3) > 0) then
        cGp3R( :,0   ) =  0.
        cDw3R( :,1   ) =  0.
        !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
        cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
      else
        cGp3R( 0,0   ) = -1./(x3p(1) - x3p(0))
        cGp3R( 1,0   ) =  1./(x3p(1) - x3p(0))
      end if

      !if (BC(2,3,g) > 0) then
      if (BCU(3) > 0) then
        cGp3R( :,kmax) =  0.
        cDw3R( :,kmax) =  0.
        !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
        cDw3R(-1,kmax) = -1./(y3w (M(3)      ) - y3w (M(3)-1  ))
      else
        cGp3R( 0,kmax) = -1./(x3p(kmax+1) - x3p(kmax))
        cGp3R( 1,kmax) =  1./(x3p(kmax+1) - x3p(kmax))
      end if
      !------------------------------------------------------------------------------------------
      do k = 1, kmax
        cdg3(-1,k) = cDw3R(-1,k)*cGp3R(0,k-1)
        cdg3( 0,k) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
        cdg3( 1,k) =                            cDw3R(0,k)*cGp3R(1,k)
      end do

      ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
      !if (BC(1,3,g) > 0) cdg3( :,1   ) = 2.*cdg3(:,1   )
      !if (BC(2,3,g) > 0) cdg3( :,kmax) = 2.*cdg3(:,kmax)
      if (BCL(3) > 0) cdg3( :,1   ) = 2.*cdg3(:,1   )
      if (BCU(3) > 0) cdg3( :,kmax) = 2.*cdg3(:,kmax)

      if (BCL(3) == -2) then
        cdg3( 1,1) = cdg3( 1,1) + cdg3(-1,1)
        cdg3(-1,1) = 0.
      end if
      if (BCU(3) == -2) then
        cdg3(-1,kmax) = cdg3(-1,kmax) + cdg3( 1,kmax)
        cdg3( 1,kmax) = 0.
      end if

    end if
    !============================================================================================

  end subroutine Op_getCDG




  !>  \brief computes \f$ \mathrm{div = \nabla\cdot\phi } \f$
  !!
  !! \note - not necessary to set bc zeros, as indexes are chosen accordingly 
  !!       - direct merging with restriction is not very attractive, as lap is
  !!         needed anyway on the finer grid.
  !!         jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-
  !!         Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei
  !!         der Interpolation verwendet (ebenfalls zur Vereinfachung der Programmierung) und ist
  !!         somit ohnehin vorhanden. Geschwindigkeitsmässig ist vermutlich auch
  !!         nicht mehr viel zu holen.
  !! \deprecated 
  subroutine OP_DivGradO2Op(  &
      dimens,                 &
      N,                      &
      SR,                     &
      ER,                     &
      bL,                     &
      bU,                     &
      BCL,                    &
      BCU,                    &
      cdg1,                   &
      cdg2,                   &
      cdg3,                   &
      phi,                    &
      Lap ) bind (c,name='OP_DivGradO2Op')

    implicit none

    integer(c_int), intent(in)    :: dimens

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: SR(3)
    integer(c_int), intent(in)    :: ER(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: BCL(3)
    integer(c_int), intent(in)    :: BCU(3)

    real(c_double), intent(in)    :: cdg1(-1:1,1:N(1))
    real(c_double), intent(in)    :: cdg2(-1:1,1:N(2))
    real(c_double), intent(in)    :: cdg3(-1:1,1:N(3))

    real(c_double), intent(in)    :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(inout) :: Lap (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i
    integer(c_int)                ::  j
    integer(c_int)                ::  k

    !===========================================================================================

    if (dimens == 3) then 
      do k = SR(3), ER(3)
        do j = SR(2), ER(2)
          !pgi$ unroll = n:8
          do i = SR(1), ER(1)
            Lap(i,j,k) =  cdg1(-1,i)*phi(i-1,j  ,k  ) + cdg1(1,i)*phi(i+1,j  ,k  )   &
              &        +  cdg2(-1,j)*phi(i  ,j-1,k  ) + cdg2(1,j)*phi(i  ,j+1,k  )   &
              &        +  cdg3(-1,k)*phi(i  ,j  ,k-1) + cdg3(1,k)*phi(i  ,j  ,k+1)   &
              &        + ( cdg1(0,i) + cdg2(0,j) + cdg3(0,k) )*phi(i,j,k)
          end do
        end do
      end do

      !=========================================================================================

    else

      do k = SR(3), ER(3)
        do j = SR(2), ER(2)
          !pgi$ unroll = n:8
          do i = SR(1), ER(1)
            Lap(i,j,k) =  cdg1(-1,i)*phi(i-1,j  ,k) + cdg1(1,i)*phi(i+1,j  ,k)   &
              &        +  cdg2(-1,j)*phi(i  ,j-1,k) + cdg2(1,j)*phi(i  ,j+1,k)   &
              &        + (cdg1(0,i) + cdg2(0,j))*phi(i,j,k)
          end do
        end do
      end do

    end if

    !===========================================================================================


    !===========================================================================================
    !=== boundary conditions ===================================================================
    !===========================================================================================

    if( BCL(1) > 0 ) then
      i = 1
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do j = SR(2), ER(2)
          Lap(i,j,k) = cdg1(0,i)*phi(i,j,k) + cdg1(1,i)*phi(i+1,j,k)
        end do
      end do
    end if

    if( BCU(1) > 0 )then
      i = N(1)
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do j = SR(2), ER(2)
          Lap(i,j,k) = cdg1(-1,i)*phi(i-1,j,k) + cdg1(0,i)*phi(i,j,k)
        end do
      end do
    end if

    !===========================================================================================

    if( BCL(2) > 0 ) then
      j = 1
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          Lap(i,j,k) = cdg2(0,j)*phi(i,j,k) + cdg2(1,j)*phi(i,j+1,k)
        end do
      end do
    end if

    if( BCU(2) > 0 ) then
      j = N(2)
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          Lap(i,j,k) = cdg2(-1,j)*phi(i,j-1,k) + cdg2(0,j)*phi(i,j,k)
        end do
      end do
    end if

    !===========================================================================================

    if( BCL(3) > 0 ) then
      k = 1
      do j = SR(2), ER(2)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          Lap(i,j,k) = cdg3(0,k)*phi(i,j,k) + cdg3(1,k)*phi(i,j,k+1)
        end do
      end do
    end if

    if( BCU(3) > 0 ) then
      k = N(3)
      do j = SR(2), ER(2)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          Lap(i,j,k) = cdg3(-1,k)*phi(i,j,k-1) + cdg3(0,k)*phi(i,j,k)
        end do
      end do
    end if

    !===========================================================================================


  end subroutine OP_DivGradO2Op



  !> \brief applies damped Jaccobi smoother on inner field
  !!
  !! \f$ \mathrm{temp = (1-\omega) phi + \omega D^{-1}(b - N phi) } \f$
  !!
  !! \param[in] dimens spatial dimension
  !! \param[in] N local grid size
  !! \param[in] bL lower storage offset
  !! \param[in] bU upper storage offset
  !! \param[in] BCL local lower boundary conditions
  !! \param[in] BCU local upper boundary conditions
  !! \param[in] cdg1 stencil in x-direction
  !! \param[in] cdg2 stencil in y-direction
  !! \param[in] cdg3 stencil in z-direction
  !! \param[in] omega damping factor
  !! \param[in] b right hand side
  !! \param[in] phi previous 'x'
  !! \param[inout] temp output
  !! 
  !! \note
  !! - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!
  !! - Direktes Verbinden mit Restriktion erscheint nicht sehr attraktiv, da Lap auf dem
  !!   jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-
  !!   Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei
  !!   der Interpolation verwendet (ebenfalls zur Vereinfachung der
  !!   Programmierung) und ist somit ohnehin vorhanden.
  !!   Geschwindigkeitsmässig ist vermutlich auch nicht mehr viel zu holen.
  !!
  !! \deprecated 
  subroutine OP_DivGradO2JSmoother(   &
      dimens,                         &
      N,                              &
      bL,                             &
      bU,                             &
      SR,                             &
      ER,                             &
      cdg1,                           &
      cdg2,                           &
      cdg3,                           &
      omega,                          &
      b,                              &
      phi,                            &
      temp ) bind (c,name='OP_DivGradO2JSmoother')

    implicit none

    integer(c_int), intent(in)    :: dimens

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: SR(3)
    integer(c_int), intent(in)    :: ER(3)

    real(c_double), intent(in)    :: cdg1(-1:1,1:N(1))
    real(c_double), intent(in)    :: cdg2(-1:1,1:N(2))
    real(c_double), intent(in)    :: cdg3(-1:1,1:N(3))

    real(c_double), intent(in)    :: omega

    real(c_double), intent(in)    :: b (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(in)    :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double), intent(inout) :: temp(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i
    integer(c_int)                ::  j
    integer(c_int)                ::  k

    !--------------------------------------------------------------

    if (dimens == 3) then

      do k = SR(3), ER(3)
        do j = SR(2), ER(2)
          !pgi$ unroll = n:8
          do i = SR(1), ER(1)
            temp(i,j,k) = (1-omega)*phi(i,j,k) +                                  &
              &    omega/( cdg1(0,i) + cdg2(0,j) + cdg3(0,k) )*(                  &
              &     b(i,j,k)                                                      &
              &  -  cdg1(-1,i)*phi(i-1,j  ,k  ) - cdg1(1,i)*phi(i+1,j  ,k  )      &
              &  -  cdg2(-1,j)*phi(i  ,j-1,k  ) - cdg2(1,j)*phi(i  ,j+1,k  )      &
              &  -  cdg3(-1,k)*phi(i  ,j  ,k-1) - cdg3(1,k)*phi(i  ,j  ,k+1)      &
              &   )

          end do
        end do
      end do


    else

      do k = SR(3), ER(3)
        do j = SR(2), ER(2)
          !pgi$ unroll = n:8
          do i = SR(1), ER(1)
            temp(i,j,k) = (1-omega)*phi(i,j,k) +                          &
              &   omega/( cdg1(0,i) + cdg2(0,j) )*(                       &
              &   b(i,j,k)                                                &
              &  -  cdg1(-1,i)*phi(i-1,j  ,k) - cdg1(1,i)*phi(i+1,j  ,k)  &
              &  -  cdg2(-1,j)*phi(i  ,j-1,k) - cdg2(1,j)*phi(i  ,j+1,k)  &
              &   )
          end do
        end do
      end do

    end if


  end subroutine OP_DivGradO2JSmoother



  !> \brief applies damped Jaccobi smoother on boundaries 
  !!
  !! \f$ \mathrm{temp = (1-\omega) phi + \omega D^{-1}(b - N phi) } \f$
  !!
  !! \param[in] N local grid size
  !! \param[in] bL lower storage offset
  !! \param[in] bU upper storage offset
  !! \param[in] BCL local lower boundary conditions
  !! \param[in] BCU local upper boundary conditions
  !! \param[in] SR start index for inner field 
  !! \param[in] ER end index for inner field 
  !! \param[in] cdg1 stencil in x-direction
  !! \param[in] cdg2 stencil in y-direction
  !! \param[in] cdg3 stencil in z-direction
  !! \param[in] omBC damping factor
  !! \param[in] b right hand side
  !! \param[in] phi previous 'x'
  !! \param[inout] temp output
  !! 
  !! \note
  !! - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!
  !! - Direktes Verbinden mit Restriktion erscheint nicht sehr attraktiv, da Lap auf dem
  !!   jeweils feineren Gitter ohnehin benötigt wird und coarse nur per MOD( , ) oder IF-
  !!   Abfrage belegt werden könnte. Zum anderen Wird das redundante feinere Feld auch bei
  !!   der Interpolation verwendet (ebenfalls zur Vereinfachung der
  !!   Programmierung) und ist somit ohnehin vorhanden.
  !!   Geschwindigkeitsmässig ist vermutlich auch nicht mehr viel zu holen.
  !! \deprecated
  subroutine OP_DivGradO2JBCSmoother( &
      N,                              &
      bL,                             &
      bU,                             &
      BCL,                            &
      BCU,                            &
      SR,                             &
      ER,                             &
      cdg1,                           &
      cdg2,                           &
      cdg3,                           &
      omBC,                           &
      b,                              &
      phi,                            &
      temp ) bind (c,name='OP_DivGradO2JBCSmoother')

    implicit none

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: BCL(3)
    integer(c_int), intent(in)    :: BCU(3)

    integer(c_int), intent(in)    :: SR(3)
    integer(c_int), intent(in)    :: ER(3)

    real(c_double), intent(in)    :: cdg1(-1:1,1:N(1))
    real(c_double), intent(in)    :: cdg2(-1:1,1:N(2))
    real(c_double), intent(in)    :: cdg3(-1:1,1:N(3))

    real(c_double), intent(in)    :: omBC

    real(c_double), intent(in)    :: b (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(in)    :: phi (bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double), intent(inout) :: temp(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i
    integer(c_int)                ::  j
    integer(c_int)                ::  k

    !--------------------------------------------------------------


    if (BCL(1) > 0) then
      i = 1
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do j = SR(2), ER(2)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg1(0,1)*( b(i,j,k) - cdg1(1,i)*phi(i+1,j,k) )
        end do
      end do
    end if

    if (BCU(1) > 0) then
      i = N(1)
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do j = SR(2), ER(2)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg1(0,i)*( b(i,j,k) - cdg1(-1,i)*phi(i-1,j,k)  )
        end do
      end do
    end if

    !===========================================================================================

    if (BCL(2) > 0) then
      j = 1
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg2(0,j)*( b(i,j,k) - cdg2(1,j)*phi(i,j+1,k) )
        end do
      end do
    end if

    if (BCU(2) > 0) then
      j = N(2)
      do k = SR(3), ER(3)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg2(0,j)*( b(i,j,k) - cdg2(-1,j)*phi(i,j-1,k) )
        end do
      end do
    end if

    !===========================================================================================

    if (BCL(3) > 0) then
      k = 1
      do j = SR(2), ER(2)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg3(0,k)*( b(i,j,k) - cdg3(1,k)*phi(i,j,k+1) )
        end do
      end do
    end if

    if (BCU(3) > 0) then
      k = N(3)
      do j = SR(2), ER(2)
        !pgi$ unroll = n:8
        do i = SR(1), ER(1)
          temp(i,j,k) = (1-omBC)*phi(i,j,k) + omBC/cdg3(0,k)*( b(i,j,k) - cdg3(-1,k)*phi(i,j,k-1) )
        end do
      end do
    end if
    !===========================================================================================

  end subroutine OP_DivGradO2JBCSmoother



  !> \brief applies Gauss--Seidel smoother
  !!
  !! \f$ \mathrm{temp = (1-\omega) phi + \omega L^{-1}(b - U phi) } \f$
  !!
  !! \param[in] N local grid size
  !! \param[in] bL lower storage offset
  !! \param[in] bU upper storage offset
  !! \param[in] BCL local lower boundary conditions
  !! \param[in] BCU local upper boundary conditions
  !! \param[in] cdg1 stencil in x-direction
  !! \param[in] cdg2 stencil in y-direction
  !! \param[in] cdg3 stencil in z-direction
  !! \param[in] omega damping factor
  !! \param[in] RBGS_mode 0:lexi, 1:???, default:2
  !! \param[in] bb right hand side
  !! \param[inout] phi working field
  !! 
  !! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes
  !!                entsprechend gewählt sind
  !!              - Bei der Initialisierung müssen die Intervallgrenzen SiR:NiR
  !!                anstelle von SiiR:N gewählt werden, da beim Aufbau der RHS
  !!                ("vec") innerhalb der Linienrelaxation die Randbereiche SiR
  !!                und NiR aufgebaut aber ggf. mit einer Korrektur wieder übe
  !!                werden. Ansonsten wäre eine Initialisierung allein im
  !!                Feldbereich SiiR:NiiR aus kann aber z.B. zu Floating Point
  !!                Exceptions führen (die für die eigentliche Rec allerdings
  !!                irrelevant wären)
  !!              - obere Stirnflächen werden bei der Initialisierung ebenfalls
  !!                berücksichtigt, da verschiedene Richtungen nacheinander
  !!                bearbeitet werden und dies beim Aufbauen der Seite sonst
  !!                berücksichtigt werden müsste. ==> geringerer Mehraufwand
  !!                beim Rechn weniger Programmierdetails.
  !!              - LU-Zerlegung (d.h. band, mult) kann gespeichert werden, wenn
  !!                mindestens eine Richtung äquidistant ist. Der Lösungsaufwand
  !!                würde sich etwa halbieren
  !!              - "r == 1 .AND. init_yes" sollte idealerweise aus den
  !!                Schleifen herausgezogen wer hier aber aus Gründen der
  !!                Übersicht bisher nicht ausgeführt wurde.
  !!
  !! Optimization annotations:
  !!    - line relaxation runs on Rosa (Cray XT5, Istanbul hex-core processors)
  !!      in single coremode with about 5.6-5.8% peak performance in direction
  !!      1, 4.2-4.6% in direction 2, 3.9-4.0% in direction 3.
  !!    - line relaxation in direction 3 is about 50% more expensive than in
  !!      direction 1; direction 2 is somewhere in between.
  !!    - line relaxation in direction 1 is typically 30% of the total execution time(!).
  !!    - RB-GS converges slightly worse than standard GS.
  !!    - RB-GS allows full vectorization (execpt some operations in direction
  !!      1), provided that the loops are reordered to let direction 1 be the
  !!      innermost (speedup was not as good as hoped (<5%) for direction 1 and
  !!      2, direction 3 was not tested).
  !!    - unrolling the non-vectorizable loops 8 times gives roughly the best
  !!      results.
  !!    - IF statements are automatically pulled out of the loops by the
  !!      compiler.
  !!
  !! \note Allgemein, pgf90: "Invariant if transformation" <==> IF-Statements werden
  !! aus dem Loop herauskopiert, muss man nicht von Hand machen!
  !! \note  Die Red-Black-Sortierung sollte für eine gute Performance global
  !! sein, d.h.  ueber die Bloecke hinweg implementiert sein. Die aktuelle
  !! Umsetzung macht sich zunutze, dass im Parallelbetrieb die Anzahl der
  !! Gitterpunkte in jeder Raumrichtung automatisch eine gerade Zahl ist (die
  !! Tatsache, dass in Randblöcken ggf. eine ungerade Anzahl Gitterpunkte
  !! abgearbeitet wird hat damit  nichts zu tun!). Daher kann sich der Shift
  !! "ss" in jedem Block einfach auf den lokalen Gitterpunkt rel(1,1,1)
  !! beziehen.
  !!
  !! Die aktuelle Schreibweise in den Schleifen ist äquivalent zu:
  !! \code
  !! DO k = S33R, N33R
  !!   DO j = S22R, N22R
  !!     DO i = S11R, N11R
  !!       IF (MOD(i+j+k) == 1) ... ! (red)   
  !!       ! bzw.                               
  !!       IF (MOD(i+j+k) == 0) ... ! (black) 
  !!     END DO                             
  !!   END DO
  !! END DO
  !! \endcode
  !!
  !! RGBS_mode |   
  !! -----------------------------
  !!  0        | lexi
  !!  1        ! ???
  !!  2        | ???
  !!  3        | nothing
  subroutine OP_DivGradO2SORSmoother(  &
      N,                              &
      bL,bU,                          &
      BCL,BCU,                        &
      cdg1,                           &
      cdg2,                           &
      cdg3,                           &
      omega,                          &
      RBGS_mode,                      &
      bb,                             &
      phi ) bind ( c, name='OP_DivGradO2SORSmoother' )

    implicit none

    integer(c_int), intent(in)    :: N(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: BCL(3)
    integer(c_int), intent(in)    :: BCU(3)

    real(c_double), intent(in)    :: cdg1(-1:1,1:N(1))
    real(c_double), intent(in)    :: cdg2(-1:1,1:N(2))
    real(c_double), intent(in)    :: cdg3(-1:1,1:N(3))

    real(c_double), intent(in)    :: omega

    integer(c_int), intent(in)    :: RBGS_mode

    real(c_double), intent(in)    :: bb(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    real(c_double), intent(inout) :: phi(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i, S11R, N11R
    integer(c_int)                ::  j, S22R, N22R
    integer(c_int)                ::  k, S33R, N33R

    integer(c_int)                ::  ss

    logical(c_bool), parameter    ::  SOR_yes   = .false. ! should be equivalent to omega=1

    !--------------------------------------------------------------

    S11R = 1
    N11R = N(1)-1

    if (BCL(1) >   0) S11R = 2
    if (BCL(1) == -2) S11R = 1
    if (BCU(1) >   0) N11R = N(1)-1
    if (BCU(1) == -2) N11R = N(1)

    S22R = 1
    N22R = N(2)-1

    if (BCL(2) >   0) S22R = 2
    if (BCL(2) == -2) S22R = 1
    if (BCU(2) >   0) N22R = N(2)-1
    if (BCU(2) == -2) N22R = N(2)

    S33R = 1
    N33R = N(3)-1

    if (BCL(3) >   0) S33R = 2
    if (BCL(3) == -2) S33R = 1
    if (BCU(3) >   0) N33R = N(3)-1
    if (BCU(3) == -2) N33R = N(3)

    !--------------------------------------------------------------

    if( RBGS_mode>0 ) then
      !==================================================================================
      !==================================================================================
      if( BCL(1)>0 ) then
        i = 1
        do k = S33R, N33R
          ss = MOD(i+k+S22R+1,2)
          !pgi$ unroll = n:8
          do j = ss+S22R, N22R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i)*phi(i+1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(1)>0 ) then
        i = N(1)
        do k = S33R, N33R
          ss = MOD(i+k+S22R+1,2)
          !pgi$ unroll = n:8
          do j = ss+S22R, N22R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i)*phi(i-1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
      if( BCL(2)>0 ) then
        j = 1
        do k = S33R, N33R
          ss = MOD(j+k+S11R+1,2)
          !pgi$ unroll = n:8
          do i = ss+S11R, N11R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j)*phi(i,j+1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(2)>0 ) then
        j = N(2)
        do k = S33R, N33R
          ss = MOD(j+k+S11R+1,2)
          !pgi$ unroll = n:8
          do i = ss+S11R, N11R, 2
            phi(i,j,k) = omega*bb(i,j,k) / cdg2(0,j)
          end do
        end do
        do k = S33R, N33R
          ss = MOD(j+k+S11R+1,2)
          !pgi$ unroll = n:8
          do i = ss+S11R, N11R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j)*phi(i,j-1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
      if( RBGS_mode==1 ) then
        !---------------------------------------------------------------------------------
        if( BCL(3)>0 ) then
          k = 1
          do j = S22R, N22R
            ss = MOD(j+k+S11R+1,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k)*phi(i,j,k+1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
        if( BCU(3)>0 ) then
          k = N(3)
          do j = S22R, N22R
            ss = MOD(j+k+S11R+1,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k)*phi(i,j,k-1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
      else
        !---------------------------------------------------------------------------------
        if( BCL(3)>0 ) then
          k = 1
          do j = S22R, N22R
            !pgi$ unroll = n:8
            do i = S11R, N11R
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k)*phi(i,j,k+1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
      end if
      !===================================================================================
      if( RBGS_mode==1 ) then
        !---------------------------------------------------------------------------------
        do k = S33R, N33R
          do j = S22R, N22R
            ss = MOD(j+k+S11R+1,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              if( SOR_yes ) then
                phi(i,j,k) = omega*(bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k)) + (1.-omega)*phi(i,j,k)
              else
                phi(i,j,k) =       (bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))
              end if
            end do
          end do
        end do
        !CALL exchange_relax(g,0,0,0,0,.TRUE.,phi) ! TEST!!! Austausch bringt praktisch nichts, macht aber die Glättung unabhängig von der Parallelisierung (sonst teilweise Jacobi-Iteration)!
        do k = S33R, N33R
          do j = S22R, N22R
            ss = MOD(j+k+S11R,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              if( SOR_yes ) then
                phi(i,j,k) = omega*(bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k)) + (1.-omega)*phi(i,j,k)
              else
                phi(i,j,k) =       (bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))
              end if
            end do
          end do
        end do
        !---------------------------------------------------------------------------------
      else
        !---------------------------------------------------------------------------------
        do k = S33R, N33R
          do j = S22R, N22R
            ss = MOD(j+k+S11R+1,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              if (SOR_yes) then
                phi(i,j,k) = omega*(bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k)) + (1.-omega)*phi(i,j,k)
              else
                phi(i,j,k) =       (bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))
              end if
            end do
          end do
          do j = S22R, N22R
            ss = MOD(j+k+S11R,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              if (SOR_yes) then
                phi(i,j,k) = omega*(bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k)) + (1.-omega)*phi(i,j,k)
              else
                phi(i,j,k) =       (bb(i,j,k)                                             &
                  &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                  &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                  &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                  &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))
              end if
            end do
          end do
        end do
        !---------------------------------------------------------------------------------
      end if
      !===================================================================================
      if( BCL(1)>0 ) then
        i = 1
        do k = S33R, N33R
          ss = MOD(i+k+S22R,2)
          !pgi$ unroll = n:8
          do j = ss+S22R, N22R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i)*phi(i+1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(1)>0 ) then
        i = N(1)
        do k = S33R, N33R
          ss = MOD(i+k+S22R,2)
          !pgi$ unroll = n:8
          do j = ss+S22R, N22R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i)*phi(i-1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
      if( BCL(2)>0 ) then
        j = 1
        do k = S33R, N33R
          ss = MOD(j+k+S11R,2)
          !pgi$ unroll = n:8
          do i = ss+S11R, N11R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j)*phi(i,j+1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(2)>0 ) then
        j = N(2)
        do k = S33R, N33R
          ss = MOD(j+k+S11R,2)
          !pgi$ unroll = n:8
          do i = ss+S11R, N11R, 2
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j)*phi(i,j-1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
      if( RBGS_mode==1 ) then
        !---------------------------------------------------------------------------------
        if( BCL(3)>0 ) then
          k = 1
          do j = S22R, N22R
            ss = MOD(j+k+S11R,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k)*phi(i,j,k+1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
        if( BCU(3)>0 ) then
          k = N(3)
          do j = S22R, N22R
            ss = MOD(j+k+S11R,2)
            !pgi$ unroll = n:8
            do i = ss+S11R, N11R, 2
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k)*phi(i,j,k-1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
      else
        !---------------------------------------------------------------------------------
        if( BCU(3)>0 ) then
          k = N(3)
          do j = S22R, N22R
            !pgi$ unroll = n:8
            do i = S11R, N11R
              phi(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k)*phi(i,j,k-1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
            end do
          end do
        end if
        !---------------------------------------------------------------------------------
      end if
      !===================================================================================

    else

      !===================================================================================
      if( BCL(1)>0 ) then
        i = 1
        do k = S33R, N33R
          !pgi$ unroll = n:8
          do j = S22R, N22R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1( 1,i)*phi(i+1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCL(2)>0 ) then
        j = 1
        do k = S33R, N33R
          !pgi$ unroll = n:8
          do i = S11R, N11R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2( 1,j)*phi(i,j+1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCL(3)>0 ) then
        k = 1
        do j = S22R, N22R
          !pgi$ unroll = n:8
          do i = S11R, N11R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg3( 1,k)*phi(i,j,k+1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
      do k = S33R, N33R
        do j = S22R, N22R
          !pgi$ unroll = n:8
          do i = S11R, N11R
            if( SOR_yes ) then
              phi(i,j,k) = omega*(bb(i,j,k)                                             &
                &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k)) + (1.-omega)*phi(i,j,k)
            else
              phi(i,j,k) =       (bb(i,j,k)                                             &
                &      - cdg1(-1,i)*phi(i-1,j,k) - cdg1(1,i)*phi(i+1,j,k)               &
                &      - cdg2(-1,j)*phi(i,j-1,k) - cdg2(1,j)*phi(i,j+1,k)               &
                &      - cdg3(-1,k)*phi(i,j,k-1) - cdg3(1,k)*phi(i,j,k+1))              &
                &      / (cdg1(0,i) + cdg2(0,j) + cdg3(0,k))
            end if
          end do
        end do
      end do
      !===================================================================================
      if( BCU(1)>0 ) then
        i = N(1)
        do k = S33R, N33R
          !pgi$ unroll = n:8
          do j = S22R, N22R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg1(-1,i)*phi(i-1,j,k)) / cdg1(0,i) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(2)>0 ) then
        j = N(2)
        do k = S33R, N33R
          !pgi$ unroll = n:8
          do i = S11R, N11R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg2(-1,j)*phi(i,j-1,k)) / cdg2(0,j) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !-----------------------------------------------------------------------------------
      if( BCU(3)>0 ) then
        k = N(3)
        do j = S22R, N22R
          !pgi$ unroll = n:8
          do i = S11R, N11R
            phi(i,j,k) = omega*(bb(i,j,k) - cdg3(-1,k)*phi(i,j,k-1)) / cdg3(0,k) + (1.-omega)*phi(i,j,k)
          end do
        end do
      end if
      !===================================================================================
    end if



  end subroutine OP_DivGradO2SORSmoother


end module cmod_DivGrad2ndOOp
