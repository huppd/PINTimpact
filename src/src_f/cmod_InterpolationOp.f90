!> \brief module providing functions to initiliaze and apply RestrictionOp
module cmod_InterpolationOp

  use iso_c_binding
  use mpi

  implicit none

contains


  subroutine MG_getCIS( &
      Nc,               &
      Nf,               &
      bL,               &
      bU,               &
      xf,               &
      dd,               &
      cI ) bind(c,name='MG_getCIS')

    implicit none

    integer(c_int), intent(in)  :: Nc
    integer(c_int), intent(in)  :: Nf
    integer(c_int), intent(in)  :: bL, bU

    real(c_double), intent(in)  :: xf(bL:(Nf+bU))

    integer(c_int), intent(in)  :: dd 

    real(c_double), intent(out) :: cI(1:2,1:Nc)

    integer(c_int)              :: i,ic
    real(c_double)              :: Dx12

    cI = 0.

    !==========================================================================================
    !=== Interpolation, linienweise, 1d =======================================================
    !==========================================================================================
    ! coarse
    !      ic                             ic+1
    !  ----o------------------------------o------------------------------o-----------------
    !      |------Dx10----|
    !      |------------Dx12--------------|
    !  ----o--------------o---------------o-----------o-----
    !     xf(i-1)        xf(i)           xf(i+1)
    ! fine


    !------------------------------------------------------------------------------------------
    do ic = 1, Nc

      i = dd * ic 
      Dx12 = xf(i+1)-xf(i-1)

      cI(1,ic) = ( xf(i+1)-xf(i  ) )/Dx12
      cI(2,ic) = ( xf(i  )-xf(i-1) )/Dx12

    end do
    !------------------------------------------------------------------------------------------

  end subroutine MG_getCIS



  !> \brief gets coefficents for interpolation of velocity points
  !!
  !! \param[in] Nc grid size on coarse grid
  !! \param[in] bLc lower ghost width on coarse grid
  !! \param[in] bUc upper ghost width on coarse grid
  !! \param[in] SSc start index on coarse grid
  !! \param[in] NNc end index on coarse grid
  !! \param[in] BC_L lower boundary conditions
  !! \param[in] BC_U upper boundary conditions
  !! \param[in] Nf grid size on fine grid
  !! \param[in] bLf lower ghost width on fine grid
  !! \param[in] bUf upper ghost width on fine grid
  !! \param[in] offset offset
  !! \param[in] xc coarse coordinates of velocity
  !! \param[in] xf fine coordinates of velocity
  !! \param[in] dd coarsening factor
  !! \param[out] cIV interpolation coefficients
  subroutine MG_getCIV( &
      Nc,               &
      bLc,              &
      bUc,              &
      SSc,              &
      NNc,              &
      BC_L,             &
      BC_U,             &
      Nf,               &
      bLf,              &
      bUf,              &
      offset,           &
      xc,               &
      xf,               &
      dd,               &
      cIV ) bind(c,name='MG_getCIV')

    implicit none

    integer(c_int), intent(in)    :: Nc

    integer(c_int), intent(in)    :: bLc, bUc

    integer(c_int), intent(in)    :: SSc
    integer(c_int), intent(in)    :: NNc

    integer(c_int), intent(in)    :: BC_L,BC_U

    integer(c_int), intent(in)    :: Nf

    integer(c_int), intent(in)    :: bLf, bUf

    integer(c_int), intent(in)    :: offset

    real(c_double), intent(in)    :: xc( bLc:(Nc+bUc) )

    real(c_double), intent(in)    :: xf( bLf:(Nf+bUf) )

    integer(c_int), intent(in)    :: dd

    real(c_double), intent(inout) :: cIV( 1:2, 0:Nf )

    integer(c_int)                ::  i, ic
    real(c_double)                ::  Dx1a, Dx12


    !==========================================================================================
    !=== Interpolation, linienweise, 1d =======================================================
    !==========================================================================================
    ! fine
    !    xf(i-1)        xf(i)         xf(i+1)       xf(i+1)      xf(i+2)
    !  ---->-------------->--------------->----------->----------->------------
    !             |-----------Dx12--------------|
    !  ----------->----------------------------->----------------------------->----------
    !          xc(ic)                      xc(ic+1)                        xc(ic+1)
    ! coarse

    cIV = 0.

    do i = 0, Nf
      ic = ( i )/dd + offset

      !Dx1a = xc(ic)-xf(i   )
      !Dx12 = xc(ic)-xc(ic+1)

      !cIV(1,i) = 1-Dx1a/Dx12
      !cIV(2,i) =   Dx1a/Dx12

      Dx12 = xc(ic+1) - xc(ic)

      cIV(1,i) = ( xc(ic+1) - xf(i ) )/Dx12
      cIV(2,i) = ( xf(i   ) - xc(ic) )/Dx12

    end do

    ! a little bit shaky, please verify this when IO is ready
    !if (BC_L > 0) then
    !cIV( 1,0) = 0.
    !cIV( 2,0) = 1.
    !end if

    !maybe false
    if( BC_L == -2 )then
      cIV(1:2,0) = 0.
    end if

    !if (BC_U > 0) then
    !cIV(1,Nf) = 1.
    !cIV(2,Nf) = 0.
    !end if
    if( BC_U == -2 )then
      cIV(1:2,Nf) = 0.
    end if

  end subroutine MG_getCIV



  !> \brief corrects values at corners through extrapolation
  !!
  !! \param[in] n
  !! \param[in] bL
  !! \param[in] bU
  !! \param[in] BCL
  !! \param[in] BCU
  !! \param[in] x1
  !! \param[in] x2
  !! \param[in] x3
  !! \param[inout] phi \f$\phi\f$
  !!
  !! \f[ \phi[2] = \frac{ \phi[4]-\phi[3] }{x[4]-x[3]} ( x[2]-x[3] )\f]
  !! \f[ \phi[n-1] = \frac{ \phi[n-2]-\phi[n-3] }{x[n-2]-x[n-3]} ( x[n-1]-x[n-2] )\f]
  subroutine MG_InterpolateCornersPost( &
      n,                                &
      bL,                               &
      bU,                               &
      BCL,                              &
      BCU,                              &
      x1,                               &
      x2,                               &
      x3,                               &
      phi ) bind(c,name='MG_InterpolateCornersPost')

    implicit none

    integer(c_int), intent(in)    :: n(3)

    integer(c_int), intent(in)    :: bL(3)
    integer(c_int), intent(in)    :: bU(3)

    integer(c_int), intent(in)    :: BCL(3)
    integer(c_int), intent(in)    :: BCU(3)

    real(c_double), intent(in)    :: x1( bL(1):(n(1)+bU(1)) )
    real(c_double), intent(in)    :: x2( bL(2):(n(2)+bU(2)) )
    real(c_double), intent(in)    :: x3( bL(3):(n(3)+bU(3)) )

    real(c_double), intent(inout) :: phi( bL(1):(n(1)+bU(1)), bL(2):(n(2)+bU(2)), bL(3):(n(3)+bU(3)) )

    integer(c_int) :: s(3), e(3), i

    do i = 1, 3 

      if( BCL(i)>0 ) then
        s(i) = 2
      else
        s(i) = 1
      end if

      if( BCU(i)>0 ) then
        e(i) = n(i)-1
      else
        e(i) = n(i)
      end if

    end do

    if( BCL(1)>0 .and. BCL(2)>0 ) then
      phi(1, 1, s(3):e(3)) = 0.
      phi(2, 1, s(3):e(3)) = ( phi( 4, 1, s(3):e(3)) - phi( 3, 1, s(3):e(3)) )/( x1(4) - x1(3) )*( x1(2) - x1(3) ) + phi( 3, 1, s(3):e(3))
      phi(1, 2, s(3):e(3)) = ( phi( 1, 4, s(3):e(3)) - phi( 1, 3, s(3):e(3)) )/( x2(4) - x2(3) )*( x2(2) - x2(3) ) + phi( 1, 3, s(3):e(3))
      phi(2, 2, s(3):e(3)) = ( &
        phi( 1, 2, s(3):e(3))*( x1(3)-x1(2) )/( x1(3)-x1(1) ) + &
        phi( 3, 2, s(3):e(3))*( x1(2)-x1(1) )/( x1(3)-x1(1) ) + &
        phi( 2, 1, s(3):e(3))*( x2(3)-x2(2) )/( x2(3)-x2(1) ) + &
        phi( 2, 3, s(3):e(3))*( x2(2)-x2(1) )/( x2(3)-x2(1) )   &
        )/2.
    endif
    if( BCL(1)>0 .and. BCU(2)>0 ) then
      phi(1, n(2)  , s(3):e(3)) =  0.
      phi(2, n(2)  , s(3):e(3)) = ( phi(4, n(2)  , s(3):e(3)) - phi(3, n(2)  , s(3):e(3)) )/( x1(4)      - x1(3)      )*( x1(2)      - x1(3)      ) + phi(3, n(2)  , s(3):e(3))
      phi(1, n(2)-1, s(3):e(3)) = ( phi(1, n(2)-2, s(3):e(3)) - phi(1, n(2)-3, s(3):e(3)) )/( x2(n(2)-2) - x2(n(2)-3) )*( x2(n(2)-1) - x2(n(2)-2) ) + phi(1, n(2)-2, s(3):e(3))
      phi(2, n(2)-1, s(3):e(3)) = ( &
        phi(1, n(2)-1, s(3):e(3))*( x1(3)     -x1(2)      )/( x1(3)-x1(1)         ) + &
        phi(3, n(2)-1, s(3):e(3))*( x1(2)     -x1(1)      )/( x1(3)-x1(1)         ) + &
        phi(2, n(2)  , s(3):e(3))*( x2(n(2)-1)-x2(n(2)-2) )/( x2(n(2))-x2(n(2)-2) ) + &
        phi(2, n(2)-2, s(3):e(3))*( x2(n(2)  )-x2(n(2)-1) )/( x2(n(2))-x2(n(2)-2) )   &
        )/2.
    endif
    if( BCU(1)>0 .and. BCL(2)>0 ) then
      phi(n(1)  , 1, s(3):e(3)) = 0.
      phi(n(1)-1, 1, s(3):e(3)) = ( phi(n(1)-2, 1, s(3):e(3)) - phi(n(1)-3, 1, s(3):e(3)) )/( x1(n(1)-2) - x1(n(1)-3) )*( x1(n(1)-1) - x1(n(1)-2) ) + phi(n(1)-2, 1, s(3):e(3)) 
      phi(n(1)  , 2, s(3):e(3)) = ( phi(n(1)  , 4, s(3):e(3)) - phi(n(1)  , 3, s(3):e(3)) )/( x2(4) - x2(3) )*( x2(2) - x2(3) )                     + phi(n(1)  , 3, s(3):e(3))  
      phi(n(1)-1, 2, s(3):e(3)) = ( &
        phi(n(1)  , 2, s(3):e(3))*( x1(n(1)-1)-x1(n(1)-2) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-2, 2, s(3):e(3))*( x1(n(1)  )-x1(n(1)-1) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-1, 1, s(3):e(3))*( x2(3)-x2(2) )/( x2(3)-x2(1) ) + &
        phi(n(1)-1, 3, s(3):e(3))*( x2(2)-x2(1) )/( x2(3)-x2(1) )   &
        )/2.                     
    endif
    if( BCU(1)>0 .and. BCU(2)>0 ) then
      phi(n(1)  , n(2)  , s(3):e(3)) = 0.
      phi(n(1)-1, n(2)  , s(3):e(3)) = ( phi(n(1)-2, n(2)  , s(3):e(3)) - phi(n(1)-3, n(2)  , s(3):e(3)) )/( x1(n(1)-2) - x1(n(1)-3) )*( x1(n(1)-1) - x1(n(1)-2) ) + phi(n(1)-2, n(2)  , s(3):e(3))
      phi(n(1)  , n(2)-1, s(3):e(3)) = ( phi(n(1)  , n(2)-2, s(3):e(3)) - phi(n(1)  , n(2)-3, s(3):e(3)) )/( x2(n(2)-2) - x2(n(2)-3) )*( x2(n(2)-1) - x2(n(2)-2) ) + phi(n(1)  , n(2)-2, s(3):e(3))
      phi(n(1)-1, n(2)-1, s(3):e(3)) = ( &
        phi(n(1)-1, n(2)  , s(3):e(3))*( x2(n(2)-1)-x2(n(2)-2) )/( x2(n(2))-x2(n(2)-2) ) + &
        phi(n(1)-1, n(2)-2, s(3):e(3))*( x2(n(2)  )-x2(n(2)-1) )/( x2(n(2))-x2(n(2)-2) ) + &
        phi(n(1)  , n(2)-1, s(3):e(3))*( x1(n(1)-1)-x1(n(1)-2) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-2, n(2)-1, s(3):e(3))*( x1(n(1)  )-x1(n(1)-1) )/( x1(n(1))-x1(n(1)-2) )   &
        )/2.                
    endif

    if( BCL(1)>0 .and. BCL(3)>0 ) then
      phi(1, s(2):e(2), 1) = 0. 
      phi(2, s(2):e(2), 1) = ( phi(4, s(2):e(2), 1)-phi(3, s(2):e(2), 1) )/( x1(4)-x1(3) )*( x1(2)-x1(3) ) + phi(3, s(2):e(2), 1)
      phi(1, s(2):e(2), 2) = ( phi(1, s(2):e(2), 4)-phi(1, s(2):e(2), 3) )/( x3(4)-x3(3) )*( x3(2)-x3(3) ) + phi(1, s(2):e(2), 3)
      phi(2, s(2):e(2), 2) = ( &
        phi(2, s(2):e(2), 1)*( x3(3)-x3(2) )/( x3(3)-x3(1) ) + &
        phi(2, s(2):e(2), 3)*( x3(2)-x3(1) )/( x3(3)-x3(1) ) + &
        phi(1, s(2):e(2), 2)*( x1(3)-x1(2) )/( x1(3)-x1(1) ) + &
        phi(3, s(2):e(2), 2)*( x1(2)-x1(1) )/( x1(3)-x1(1) )   &
        )/2. 
    endif                    
    if( BCL(1)>0 .and. BCU(3)>0 ) then
      phi(1, s(2):e(2), n(3)  ) = 0
      phi(2, s(2):e(2), n(3)  ) = ( phi(4, s(2):e(2), n(3)  )-phi(3, s(2):e(2), n(3)  ) )/( x1(4)-x1(3) )*( x1(2)-x1(3) )                     + phi(3, s(2):e(2), n(3)  )
      phi(1, s(2):e(2), n(3)-1) = ( phi(1, s(2):e(2), n(3)-2)-phi(1, s(2):e(2), n(3)-3) )/( x3(n(3)-2)-x3(n(3)-3) )*( x3(n(3)-1)-x3(n(3)-2) ) + phi(1, s(2):e(2), n(3)-2)
      phi(2, s(2):e(2), n(3)-1) = ( &
        phi(2, s(2):e(2), n(3)  )*( x3(n(3)-1)-x3(n(3)-2) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(2, s(2):e(2), n(3)-2)*( x3(n(3)  )-x3(n(3)-1) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(1, s(2):e(2), n(3)-1)*( x1(3)-x1(2) )/( x1(3)-x1(1) ) + &
        phi(3, s(2):e(2), n(3)-1)*( x1(2)-x1(1) )/( x1(3)-x1(1) ) &
        )/2.         
    endif
    if( BCU(1)>0 .and. BCL(3)>0 ) then
      phi(n(1)  , s(2):e(2), 1) = 0.
      phi(n(1)-1, s(2):e(2), 1) = ( phi(n(1)-2, s(2):e(2), 1) - phi(n(1)-3, s(2):e(2), 1) )/( x1(n(1)-2)-x1(n(1)-3) )*( x1(n(1)-1)-x1(n(1)-2) ) + phi(n(1)-2, s(2):e(2), 1)
      phi(n(1)  , s(2):e(2), 2) = ( phi(n(1)  , s(2):e(2), 4) - phi(n(1)  , s(2):e(2), 3) )/( x3(4)     -x3(3)      )*( x3(2)     -x3(3)      ) + phi(n(1)  , s(2):e(2), 3)
      phi(n(1)-1, s(2):e(2), 2) = ( &
        phi(n(1)  , s(2):e(2), 2)*( x1(n(1)-1)-x1(n(1)-2) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-2, s(2):e(2), 2)*( x1(n(1)  )-x1(n(1)-1) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-1, s(2):e(2), 1)*( x3(3)     -x3(2)      )/( x3(3)   -x3(1)      ) + &
        phi(n(1)-1, s(2):e(2), 3)*( x3(2)     -x3(1)      )/( x3(3)   -x3(1)      )   &
        )/2.               
    endif
    if( BCU(1)>0 .and. BCU(3)>0 ) then
      phi(n(1)  , s(2):e(2)  , n(3)  ) = 0.
      phi(n(1)-1, s(2):e(2), n(3)  ) = ( phi(n(1)-2, s(2):e(2), n(3)  )-phi(n(1)-3, s(2):e(2), n(3)  ) )/( x1(n(1)-2)-x1(n(1)-3) )*( x1(n(1)-1)-x1(n(1)-2) ) + phi(n(1)-2, s(2):e(2), n(3)  )
      phi(n(1)  , s(2):e(2), n(3)-1) = ( phi(n(1)  , s(2):e(2), n(3)-2)-phi(n(1)  , s(2):e(2), n(3)-3) )/( x3(n(3)-2)-x3(n(3)-3) )*( x3(n(3)-1)-x3(n(3)-2) ) + phi(n(1)  , s(2):e(2), n(3)-2)
      phi(n(1)-1, s(2):e(2), n(3)-1) = ( &
        phi(n(1)-1, s(2):e(2), n(3)  )*( x3(n(3)-1)-x3(n(3)-2) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(n(1)-1, s(2):e(2), n(3)-2)*( x3(n(3)  )-x3(n(3)-1) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(n(1)  , s(2):e(2), n(3)-1)*( x1(n(1)-1)-x1(n(1)-2) )/( x1(n(1))-x1(n(1)-2) ) + &
        phi(n(1)-2, s(2):e(2), n(3)-1)*( x1(n(1)  )-x1(n(1)-1) )/( x1(n(1))-x1(n(1)-2) )   &
        )/2. 
    endif

    if( BCL(2)>0 .and. BCL(3)>0 ) then
      phi(s(1):e(1), 1, 1 ) = 0.
      phi(s(1):e(1), 2, 1 ) = ( phi(s(1):e(1), 4, 1 )-phi(s(1):e(1), 3, 1 ) )/( x2(4)-x2(3) )*( x2(2)-x2(3) ) + phi(s(1):e(1), 3, 1 )
      phi(s(1):e(1), 1, 2 ) = ( phi(s(1):e(1), 1, 4 )-phi(s(1):e(1), 1, 3 ) )/( x3(4)-x3(3) )*( x3(2)-x3(3) ) + phi(s(1):e(1), 1, 3 )
      phi(s(1):e(1), 2, 2 ) = ( &
        phi(s(1):e(1), 2, 1 )*( x3(3)-x3(2) )/( x3(3)-x3(1) ) + &
        phi(s(1):e(1), 2, 3 )*( x3(2)-x3(1) )/( x3(3)-x3(1) ) + &
        phi(s(1):e(1), 1, 2 )*( x2(3)-x2(2) )/( x2(3)-x2(1) ) + &
        phi(s(1):e(1), 3, 2 )*( x2(2)-x2(1) )/( x2(3)-x2(1) )   &
        )/2.
    endif
    if( BCL(2)>0 .and. BCU(3)>0 ) then
      phi(s(1):e(1), 1, n(3)  ) = 0.
      phi(s(1):e(1), 2, n(3)  ) = ( phi(s(1):e(1), 4, n(3)  )-phi(s(1):e(1), 3, n(3)  ) )/( x2(4)     -x2(3)      )*( x2(2)     -x2(3)      ) + phi(s(1):e(1), 3, n(3)  )
      phi(s(1):e(1), 1, n(3)-1) = ( phi(s(1):e(1), 1, n(3)-2)-phi(s(1):e(1), 1, n(3)-3) )/( x3(n(3)-2)-x3(n(3)-3) )*( x3(n(3)-1)-x3(n(3)-2) ) + phi(s(1):e(1), 1, n(3)-2)
      phi(s(1):e(1), 2, n(3)-1) = ( &
        phi(s(1):e(1), 2, n(3)  )*( x3(n(3)-1)-x3(n(3)-2) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(s(1):e(1), 2, n(3)-2)*( x3(n(3)  )-x3(n(3)-1) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(s(1):e(1), 1, n(3)-1)*( x2(3)     -x2(2)      )/( x2(3)-x2(1)         ) + &
        phi(s(1):e(1), 3, n(3)-1)*( x2(2)     -x2(1)      )/( x2(3)-x2(1)         ) &
        )/2.      
    endif
    if( BCU(2)>0 .and. BCL(3)>0 ) then
      phi(s(1):e(1), n(2)  , 1 ) = 0.
      phi(s(1):e(1), n(2)-1, 1 ) = ( phi(s(1):e(1), n(2)-2, 1)-phi(s(1):e(1), n(2)-3, 1) )/( x2(n(2)-2)-x2(n(2)-3) )*( x2(n(2)-1)-x2(n(2)-2) ) + phi(s(1):e(1), n(2)-2, 1)
      phi(s(1):e(1), n(2)  , 2 ) = ( phi(s(1):e(1), n(2)  , 4)-phi(s(1):e(1), n(2)  , 3) )/( x3(4)     -x3(3)      )*( x3(2)     -x3(3)      ) + phi(s(1):e(1), n(2),   3)
      phi(s(1):e(1), n(2)-1, 2 ) = ( &
        phi(s(1):e(1), n(2)-1, 1)*( x3(3)     -x3(2)      )/( x3(3)-x3(1)         ) + &
        phi(s(1):e(1), n(2)-1, 3)*( x3(2)     -x3(1)      )/( x3(3)-x3(1)         ) + &
        phi(s(1):e(1), n(2)  , 2)*( x2(n(2)-1)-x2(n(2)-2) )/( x2(n(2))-x2(n(2)-2) ) + &
        phi(s(1):e(1), n(2)-2, 2)*( x2(n(2)  )-x2(n(2)-1) )/( x2(n(2))-x2(n(2)-2) ) &
        )/2.     
    endif
    if( BCU(2)>0 .and. BCU(3)>0 ) then
      phi(s(1):e(1), n(2)  , n(3)  ) = 0.
      phi(s(1):e(1), n(2)-1, n(3)  ) = ( phi(s(1):e(1), n(2)-2, n(3)  )-phi(s(1):e(1), n(2)-3, n(3)  ) )/( x2(n(2)-2)-x2(n(2)-3) )*( x2(n(2)-1)-x2(n(2)-2) ) + phi(s(1):e(1), n(2)-2, n(3)  )
      phi(s(1):e(1), n(2)  , n(3)-1) = ( phi(s(1):e(1), n(2)  , n(3)-2)-phi(s(1):e(1), n(2)  , n(3)-3) )/( x3(n(3)-2)-x3(n(3)-3) )*( x3(n(3)-1)-x3(n(3)-2) ) + phi(s(1):e(1), n(2)  , n(3)-2)
      phi(s(1):e(1), n(2)-1, n(3)-1) = (                                                   &
        phi(s(1):e(1), n(2)-1, n(3)  )*( x3(n(3)-1)-x3(n(3)-2) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(s(1):e(1), n(2)-1, n(3)-2)*( x3(n(3)  )-x3(n(3)-1) )/( x3(n(3))-x3(n(3)-2) ) + &
        phi(s(1):e(1), n(2),   n(3)-1)*( x2(n(2)-1)-x2(n(2)-2) )/( x2(n(2))-x2(n(2)-2) ) + &
        phi(s(1):e(1), n(2)-2, n(3)-1)*( x2(n(2)  )-x2(n(2)-1) )/( x2(n(2))-x2(n(2)-2) )   &
        )/2.
    endif

  end subroutine MG_InterpolateCornersPost



  !> \note:       - für allgemeine dd geeignet
  !!              - dd(i) /= 1 <===> N(i) /= 1
  !!              - es wird nur in eine Richung ausgetauscht
  !!              - Null-Setzen am Rand nicht notwendig
  !!              - Es wird sequentiell über alle Raumrichtungen interpoliert,
  !!                um keinen individuellen Interpolationsstencil für jeden
  !!                Punkt im Raum speichern zu müssen.
  !!              - Durch das sequentielle Interpolieren kann der
  !!                Interpolationsstencil klein und damit der Gesamtaufwand
  !!                minimiert werden (Alternative: 8- bzw. 26-Punkt Stencil(!)).
  !!              - Interpolationskoeffizienten werden auf dem jeweils feineren
  !!                Gitter gespeichert, um nicht auf die entsprechenden Indizes
  !!                des gröberen Gitters umrechnen zu müssen.
  !!              - Die Block-überlappenden Stirnflächen werden ebenfalls
  !!                mitverarbeitet, aber eigentlich nicht gebraucht (erleichtert
  !!                die Programmierung), so dass eigentlich eine Initialisierung
  !!                notwendig wäre. Dies wird jedoch zuvor schon in der
  !!                korrespondierenden Restriktions-Routine erledigt, so dass
  !!                dies hier nicht mehr notwendig ist.
  !!
  !! \attention: Verwendung von parametrischen Strides dd in Schleifen
  !! verhindert die Vektorisierung bzw. das Prefetching! Andererseits ist der
  !! Geschwindigkeitsgewinn nur sehr gering (schon getestet).
  !!
  !! - Geschwindigkeit ist trotz Prefetching / Vektorisierung stark durch die
  !!   Speicherzugriffszeit limitiert
  !! - Das wird z.B. deutlich bei Single- vs.  Dualcorebetrieb
  subroutine MG_interpolate(  &
      Nc,                     &
      bLc,                    &
      bUc,                    &
      Nf,                     &
      bLf,                    &
      bUf,                    &
      iimax,                  &
      dd,                     &
      cI1,                    &
      cI2,                    &
      cI3,                    &
      phic,                   &
      phif ) bind(c,name='MG_interpolate')

    implicit none

    integer(c_int), intent(in)     :: Nc(3)

    integer(c_int), intent(in)     :: bLc(3)
    integer(c_int), intent(in)     :: bUc(3)

    integer(c_int), intent(in)     :: Nf(3)

    integer(c_int), intent(in)     :: bLf(3)
    integer(c_int), intent(in)     :: bUf(3)

    integer(c_int), intent(in)     :: iimax(3)
    integer(c_int), intent(in)     :: dd(3)

    real(c_double),  intent(in)    :: cI1 ( 1:2, 1:Nc(1) )
    real(c_double),  intent(in)    :: cI2 ( 1:2, 1:Nc(2) )
    real(c_double),  intent(in)    :: cI3 ( 1:2, 1:Nc(3) )

    real(c_double),  intent(in) :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))

    real(c_double),  intent(out)   :: phif (bLf(1):(Nf(1)+bUf(1)),bLf(2):(Nf(2)+bUf(2)),bLf(3):(Nf(3)+bUf(3)))

    integer(c_int)               ::  i,ic
    integer(c_int)               ::  j
    integer(c_int)               ::  k



    !******************************************************************************************
    !pgi$ unroll = n:8
    phif( 1:Nf(1):dd(1), 1:Nf(2):dd(2), 1:Nf(3):dd(3) ) = phic( 1:iimax(1), 1:iimax(2), 1:iimax(3) )
    !*****************************************************************************************
    !=========================================================================================
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
    !=========================================================================================
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
    !=========================================================================================
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
    !=========================================================================================

  end subroutine MG_interpolate



  !> \brief interpolating velocities
  subroutine MG_interpolateV( &
      dir,                    &
      Nc,                     &
      bLc,bUc,                &
      SSc,NNc,                &
      Nf,                     &
      bLf,bUf,                &
      SSf,NNf,                &
      iimax,                  &
      dd,                     &
      cIV,                    &
      cI1,                    &
      cI2,                    &
      cI3,                    &
      phic,                   &
      phif ) bind (c,name='MG_interpolateV')

    implicit none

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

    !integer(c_int), intent(in)     :: BCL(1:3)
    !integer(c_int), intent(in)     :: BCU(1:3)

    integer(c_int), intent(in)     :: iimax(1:3)
    integer(c_int), intent(in)     :: dd(1:3)

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

    integer(c_int)                :: S(1:3)
    integer(c_int)                :: N(1:3)
    integer(c_int)                :: ds(1:3)


    ds = dd

    do i = 1,3
      if( iimax(i)<Nc(i) ) ds = 1
    end do


    if( 1==dir ) then

      if( dd(1) /= 1 ) then

        do k = 1, Nf(3), dd(3)
          if( 1==dd(3) ) then
            kc = k
          else
            kc = ( k+1 )/dd(3)
          end if
          do j = 1, Nf(2), dd(2)
            if( 1==dd(2) ) then
              jc = j
            else
              jc = ( j+1 )/dd(2) ! holy shit
            end if
            do i = 0, Nf(1) ! zero for dirichlet
              ic = ( i )/dd(1)
              phif(i,j,k) = cIV(1,i)*phic(ic,jc,kc) + cIV(2,i)*phic(ic+1,jc,kc)
            end do
          end do
        end do

      else

        do k = 1, Nf(3), dd(3)
          if( 1==dd(3) ) then
            kc = k
          else
            kc = ( k+1 )/dd(3)
          end if
          do j = 1, Nf(2), dd(2)
            if( 1==dd(2) ) then
              jc = j
            else
              jc = ( j+1 )/dd(2) ! holy shit
            end if
            do i = 0, Nf(1)
              phif(i,j,k) = phic(i,jc,kc)
            end do
          end do
        end do

      end if

      if( dd(2) /= 1 ) then

        do k = 1, Nf(3), dd(3)
          do j = 2, Nf(2)-1, dd(2)
            jc = ( j )/dd(2)
            !pgi$ unroll = n:8
            do i = 0, Nf(1)
              phif(i,j,k) = cI2(1,jc)*phif(i,j-1,k) + cI2(2,jc)*phif(i,j+1,k)
            end do
          end do
        end do

      end if

      if( dd(3) /= 1 ) then

        do k = 2, Nf(3)-1, dd(3)
          kc = k/dd(3)
          do j = 1, Nf(2)
            !pgi$ unroll = n:8
            do i = 0, Nf(1)
              phif(i,j,k) = cI3(1,kc)*phif(i,j,k-1) + cI3(2,kc)*phif(i,j,k+1)
            end do
          end do
        end do

      end if

    end if

    if( 2==dir ) then

      if( dd(2) /= 1 ) then

        do k = 1, Nf(3), dd(3)
          if( 1==dd(3) ) then
            kc = k
          else
            kc = ( k+1 )/dd(3)
          end if
          do j = 0, Nf(2) ! zero for dirichlet
            jc = ( j )/dd(2)
            do i = 1, Nf(1), dd(1)
              if( 1==dd(1) ) then
                ic = i
              else
                ic = ( i+1 )/dd(1) ! holy shit
              end if
              phif(i,j,k) = cIV(1,j)*phic(ic,jc,kc)+cIV(2,j)*phic(ic,jc+1,kc)
            end do
          end do
        end do

      else

        do k = 1, Nf(3), dd(3)
          if( 1==dd(3) ) then
            kc = k
          else
            kc = ( k+1 )/dd(3)
          end if
          do j = 0, Nf(2) 
            do i = 1, Nf(1), dd(1)
              if( 1==dd(1) ) then
                ic = i
              else
                ic = ( i+1 )/dd(1) ! holy shit
              end if
              phif(i,j,k) = phic(ic,j,kc)
            end do
          end do
        end do

      end if

      if( dd(1) /= 1 ) then

        do k = 1, Nf(3), dd(3)
          do j = 0, Nf(2)
            !pgi$ unroll = n:8
            do i = 2, Nf(1)-1, dd(1)
              ic = i/dd(1)
              phif(i,j,k) = cI1(1,ic)*phif(i-1,j,k) + cI1(2,ic)*phif(i+1,j,k)
            end do
          end do
        end do

      end if

      if( dd(3) /= 1 ) then

        do k = 2, Nf(3)-1, dd(3)
          kc = k/dd(3)
          do j = 0, Nf(2)
            !pgi$ unroll = n:8
            do i = 1, Nf(1)
              phif(i,j,k) = cI3(1,kc)*phif(i,j,k-1) + cI3(2,kc)*phif(i,j,k+1)
            end do
          end do
        end do

      end if

    end if

    if( 3==dir ) then

      if( dd(3) /= 1 ) then

        do k = 0, Nf(3)
          kc = ( k )/dd(3)
          do j = 1, Nf(2), dd(2) ! zero for dirichlet
            if( 1==dd(2) ) then
              jc = j
            else
              jc = ( j+1 )/dd(2)
            end if
            do i = 1, Nf(1), dd(1)
              if( 1==dd(1) ) then
                ic = i
              else
                ic = ( i+1 )/dd(1) ! holy shit
              end if
              phif(i,j,k) = cIV(1,k)*phic(ic,jc,kc)+cIV(2,k)*phic(ic,jc,kc+1)
            end do
          end do
        end do

      else

        do k = 0, Nf(3)
          do j = 1, Nf(2), dd(2)
            if( 1==dd(2) ) then
              jc = j
            else
              jc = ( j+1 )/dd(2)
            end if
            do i = 1, Nf(1), dd(1)
              if( 1==dd(1) ) then
                ic = i
              else
                ic = ( i+1 )/dd(1) ! holy shit
              end if
              phif(i,j,k) = phic(ic,jc,k)
            end do
          end do
        end do

      end if

      if( dd(1) /= 1 ) then

        do k = 0, Nf(3)
          do j = 1, Nf(2), dd(2)
            !pgi$ unroll = n:8
            do i = 2, Nf(1)-1, dd(1)
              ic = i/dd(1)
              phif(i,j,k) = cI1(1,ic)*phif(i-1,j,k) + cI1(2,ic)*phif(i+1,j,k)
            end do
          end do
        end do

      end if

      if( dd(2) /= 1 ) then

        do k = 0, Nf(3)
          do j = 2, Nf(2)-1, dd(2)
            jc = j/dd(2)
            !pgi$ unroll = n:8
            do i = 1, Nf(1)
              phif(i,j,k) = cI2(1,jc)*phif(i,j-1,k) + cI2(2,jc)*phif(i,j+1,k)
            end do
          end do
        end do

      end if

    end if

  end subroutine MG_interpolateV



  !> \brief 
  !! \note todo understand recevbuf size not correct
  subroutine MG_InterpolateScatter( &
      Nc,                           &
      bLc,                          &
      bUc,                          &
      np,                           &
      iimax,                        &
      n_gather,                     &
      participate_yes,              &
      rankc2,                       &
      comm2,                        &
      dispI,                        &
      offsI,                        &
      phic ) bind(c,name='MG_InterpolateScatter')

    implicit none


    integer(c_int), intent(in)     :: Nc(3)

    integer(c_int), intent(in)     :: bLc(3)
    integer(c_int), intent(in)     :: bUc(3)

    integer(c_int), intent(in)     :: np(3)

    integer(c_int), intent(in)     :: iimax(3)

    integer(c_int), intent(in)     :: n_gather(3)

    logical(c_bool), intent(in)    :: participate_yes

    integer(c_int), intent(in)     :: rankc2

    integer(c_int), intent(in)     :: comm2

    integer(c_int), intent(in)     :: dispI(     1:n_gather(1)*n_gather(2)*n_gather(3) )
    integer(c_int), intent(in)     :: offsI(1:3, 1:n_gather(1)*n_gather(2)*n_gather(3) )


    real(c_double),intent(inout)   :: phic (bLc(1):(Nc(1)+bUc(1)),bLc(2):(Nc(2)+bUc(2)),bLc(3):(Nc(3)+bUc(3)))


    integer(c_int)                ::  i, ii
    integer(c_int)                ::  j, jj
    integer(c_int)                ::  k, kk

    integer(c_int)                ::  merror
    real(c_double), allocatable   ::  recvbuf(:,:,:)
    integer(c_int)                ::  offsg(1:3), dispg
    !real(c_double)                ::  sendbuf( 1:( ( Nc(1)+np(1)-1 )*( Nc(2)+np(2)-1 )*( Nc(3)+np(3)-1 ) ) ) ! from IMPACT
    !real(c_double)                ::  sendbuf( 0:( ( Nc(1)+np(1)*3 )*( Nc(2)+np(2)*3 )*( Nc(3)+np(3)*3 ) ) )
    real(c_double)                ::  sendbuf( 0:( (iimax(1)+1)*(iimax(2)+1)*(iimax(3)+1)*n_gather(1)*n_gather(2)*n_gather(3) ) )
    !real(c_double)                ::  sendbuf( 0:( ( Nc(1)+12 )*( Nc(2)+12 )*( Nc(3)+12 ) ) )



    if( participate_yes ) then

      do k = 1, n_gather(3)
        do j = 1, n_gather(2)
          do i = 1, n_gather(1)

            dispg      = dispI(    i+(j-1)*n_gather(1)+(k-1)*n_gather(1)*n_gather(2))
            offsg(1:3) = offsI(1:3,i+(j-1)*n_gather(1)+(k-1)*n_gather(1)*n_gather(2))

            do kk = 0, iimax(3)
              do jj = 0, iimax(2)
                do ii = 0, iimax(1)
                  sendbuf( dispg + ii + jj*(iimax(1)+1) + kk*(iimax(1)+1)*(iimax(2)+1) ) = phic(ii+offsg(1),jj+offsg(2),kk+offsg(3))
                end do
              end do
            end do

          end do
        end do
      end do

    end if

    allocate(recvbuf(0:iimax(1),0:iimax(2),0:iimax(3))) ! Anmerkung: Besser nicht fest allocieren um Speicherplatz zu sparen, ODER gleich "phic" verwenden!

    call MPI_SCATTER(sendbuf,(iimax(1)+1)*(iimax(2)+1)*(iimax(3)+1),MPI_REAL8,recvbuf,(iimax(1)+1)*(iimax(2)+1)*(iimax(3)+1),MPI_REAL8,rankc2,comm2,merror)

    phic(0:iimax(1),0:iimax(2),0:iimax(3)) = recvbuf(0:iimax(1),0:iimax(2),0:iimax(3))


    deallocate(recvbuf)
    !------------------------------------------------------------------------------


  end subroutine MG_InterpolateScatter


end module cmod_InterpolationOp
