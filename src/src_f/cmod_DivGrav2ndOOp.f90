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

  subroutine Op_getCDG_dir( &
      M,                    &
      N,                    &
      BL,                   &
      BU,                   &
      BCL,                  &
      BCU,                  &
      yu,                   &
      xp,                   &
      xu,                   &
      cDG ) bind( c, name='Op_getCDG_dir' )

    implicit none

    integer(c_int), intent(in)    :: M

    integer(c_int), intent(in)    :: N

    integer(c_int), intent(in)    :: BL
    integer(c_int), intent(in)    :: BU

    integer(c_int), intent(in)    :: BCL
    integer(c_int), intent(in)    :: BCU

    real(c_double), intent(in)    :: yu(0:M)

    real(c_double), intent(in)    :: xp(bl:N+bu)

    real(c_double), intent(in)    :: xu(bl:N+bu)

    real(c_double), intent(inout) :: cDG(-1:1,1:N)

    real(c_double)                :: cGp( 0:1,0:N)

    real(c_double)                :: cDu(-1:0,1:N)

    integer(c_int)                :: i

    cDG = 0.

    cGp = 0.
    cDu = 0.

    !---------------------------------------------------------------------------!
    ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
    !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
    !---------------------------------------------------------------------------!

    do i = 1, N-1
      cGp(0,i) = -1./( xp(i+1) - xp(i) )
      cGp(1,i) =  1./( xp(i+1) - xp(i) )
    end do

    do i = 1, N
      cDu(-1,i) = -1./( xu(i) - xu(i-1) )
      cDu( 0,i) =  1./( xu(i) - xu(i-1) )
    end do

    if( BCL > 0 ) then
      cGp(:,0) = 0.
      cDu(:,1) = 0.
      !cDu(0,1) = 1./( yu(1) - yu(0) ) ! TEST!!! Der Schoenheit halber, s.u.
      cDu(0,1) =  1./(xu(1) - xu (0)) !! shoudl be the same
      !cDu(0,1) =  0.5/( xu(1) - xp(1) )
    else
      cGp(0,0) = -1./( xp(1) - xp(0) )
      cGp(1,0) =  1./( xp(1) - xp(0) )
    end if

    if( BCU > 0 ) then
      cGp(:,N) =  0.
      cDu(:,N) =  0.
      !cDu(-1,N) = -1./( yu(M) - yu(M-1) )
      cDu(-1,N) = -1./( xu(N) - xu(N-1) ) ! TEST!!! Das geht in die Hose ...
      !cDu(-1,N) = -1./( yu(M+1) - yu(M) )
      !cDu(-1,N) = -0.5/( xu(N) - xp(N) ) ! TEST!!! Das geht in die Hose ...
      ! Aber Warum?
    else
      cGp( 0,N) = -1./( xp(N+1) - xp(N) )
      cGp( 1,N) =  1./( xp(N+1) - xp(N) )
    end if
    !-------------------------------------------------------------------------------------------
    do i = 1, N
      cDG(-1,i) = cDu(-1,i)*cGp(0,i-1)
      cDG( 0,i) = cDu(-1,i)*cGp(1,i-1) + cDu(0,i)*cGp(0,i)
      cDG( 1,i) =                        cDu(0,i)*cGp(1,i)
    end do

    ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das
    ! Interpolationspolynom am Rand darstellt
    if( BCL > 0 ) then
      cDG( :,1) = 2.*cDG(:,1) 
      !cDG(-1,1) = 0. 
      !cDG( 0,1) =  cDG(0,2) 
      !cDG( 1,1) = -cDG(0,2) 
    end if
    if( BCU > 0 ) then 
      cDG( :,N) = 2.*cDG(:,N)
      !cDG(-1,N) =  cDG(0,N-1) 
      !cDG( 0,N) = -cDG(0,N-1) 
      !cDG( 1,N) = 0. 
    end if

    if( BCL == -2) then
      cDG( 1,1) = cDG( 1, 1 ) + cDG( -1, 1 )
      cDG(-1,1) = 0.
    end if

    if( BCU == -2 ) then
      cDG(-1,N) = cDG(-1, N ) + cDG( 1, N )
      cDG( 1,N) = 0.
    end if

  end subroutine Op_getCDG_dir


end module cmod_DivGrad2ndOOp
