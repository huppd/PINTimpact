!>  Modul: cmod_VectorField
!!
!! Impact functions \c Pimpact::VectorField e.g. scales, norms ...
!! \author huppd
module cmod_VectorField

  use iso_c_binding

  implicit none

contains

  !> \brief calculate distance to immersed boundary
  !!
  !! \param[in] x
  !! \param[in] y
  !! \param[in] z
  !! \param[in] x0
  !! \param[in] y0
  !! \param[in] z0
  !! \param[in] R
  !! \param[in] dr
  FUNCTION distance2ib( x, y, z, x0, y0, R, dr ) RESULT(dis)

    IMPLICIT NONE

    REAL(c_double), INTENT(IN)  ::  x,y,z

    REAL(c_double), intent(in)  ::  x0
    REAL(c_double), intent(in)  ::  y0

    REAL(c_double), intent(in)  ::  R

    REAL(c_double), intent(in)  ::  dr

    REAL(c_double)              ::  z0

    REAL(c_double)              ::  dis


    ! geometric properties of disc
    !        x0 = L1/2. + L1*amp*SIN(2.*pi*freq*subtime)
    !        y0 = L2/4.
    !        z0 = L3/2.

    !        R = L1/10.
    !  write(*,*) 'dr=',dr

    dis = SQRT( (x0-x)**2 + (y0-y)**2 )

    if( dis <= R) then
      dis = 1.
    else if( dis <= R+dr) then
      dis = ABS(1./( 1. + EXP(1./(-ABS(dis-R)/dr) + 1./(1.-Abs(dis-R)/dr) ) ));
    else
      !IF( Abs(dis-R) < dr) THEN
      !dis = ABS(1./( 1. + EXP(1./(-ABS(dis-R)/dr) + 1./(1.-Abs(dis-R)/dr) ) ));
      !!    dis = (1-(dis-R)/dr)
      !ELSE
      dis = 0.
    end if
    !    dis = 0.

  END FUNCTION distance2ib



  !> \brief get initial slope/curvature for boundary layer equations as a
  !! function of sweep angle and suction velocity
  !!
  !! \param[in] kappa non-dimensional boundary layer suction velocity
  !! \param[in] sweep_angle_degrees sweep angle (in degrees for ease of interpolation)
  !! \param[out] fxxGuess final guess values (return values)
  !! \param[out] gxGuess
  SUBROUTINE init_value_bl_equation( kappa, sweep_angle_degrees, fxxGuess, gxGuess )

    implicit none

    !----------- mjohn 101111 ------------------------------------------------------------------
    ! - INPUT VARIABLES
    real(c_double), intent(in)    :: kappa
    real(c_double), intent(in)    :: sweep_angle_degrees

    ! - OUTPUT VARIABLES
    real(c_double), intent(out)   ::   fxxGuess, gxGuess

    real(c_double)   :: fxx(1:10,7), gx(1:10,7)    ! initial value matrices for f and g (structure see where declared)
    real(c_double)   ::   fxx1, fxx2, fxx3, fxx4   ! interpolation stencil (2D, since kappa and sweep_angle vary) for f
    real(c_double)   ::   gx1, gx2, gx3, gx4       ! interpolation stencil (2D, since kappa and sweep_angle vary) for g
    real(c_double)   ::   dK, dPhi                 ! interpolation increment remainders (interpolation weights)
    real(c_double)   ::   fxxLo, fxxHi, gxLo, gxHi ! more dummy variables...
    !------------------------------------------------------------------------------------------


    ! v''(0) for kappa = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0) and 
    ! phi =       (     0        10        20        30        40        50        60        70        80      90   )
    fxx(1:10,1) = (/ 1.232578, 1.204596, 1.122776, 0.993369, 0.826410, 0.635208, 0.435783, 0.246543, 0.089191, 0.0 /)
    fxx(1:10,2) = (/ 1.541729, 1.509206, 1.413891, 1.262421, 1.065455, 0.854631, 0.594496, 0.357035, 0.147418, 0.0 /)
    fxx(1:10,3) = (/ 1.889298, 1.851881, 1.742007, 1.566680, 1.337148, 1.068325, 0.778115, 0.486851, 0.217475, 0.0 /)
    fxx(1:10,4) = (/ 2.279667, 2.236852, 2.110903, 1.909177, 1.643485, 1.329463, 0.985815, 0.633578, 0.295836, 0.0 /) !! DUMMY VALUES (mere interpolation of line above and below)
    fxx(1:10,5) = (/ 2.670036, 2.621822, 2.479798, 2.251675, 1.949821, 1.590600, 1.193515, 0.780304, 0.374196, 0.0 /)
    fxx(1:10,6) = (/ 3.098327, 3.044144, 2.884339, 2.626989, 2.285053, 1.884300, 1.419229, 0.938387, 0.457070, 0.0 /) !! DUMMY VALUES (mere interpolation of line above and below)
    fxx(1:10,7) = (/ 3.526617, 3.466466, 3.288879, 3.002302, 2.620284, 2.178000, 1.644944, 1.096470, 0.539944, 0.0 /)
    ! w'(0):
    gx(1:10,1) = (/ 0.570454, 0.566105, 0.552986, 0.530868, 0.499286, 0.457358, 0.403374, 0.333619, 0.237717, 0.0 /)
    gx(1:10,2) = (/ 0.921733, 0.917607, 0.905180, 0.884295, 0.854631, 0.815588, 0.766028, 0.703608, 0.622597, 0.5 /)
    gx(1:10,3) = (/ 1.323674, 1.319995, 1.308946, 1.290479, 1.264492, 1.230778, 1.188932, 1.138139, 1.076642, 1.0 /)
    gx(1:10,4) = (/ 1.767516, 1.764282, 1.754585, 1.738445, 1.715877, 1.686874, 1.651353, 1.609055, 1.559273, 1.5 /) !! DUMMY VALUES (mere interpolation of line above and below)
    gx(1:10,5) = (/ 2.211359, 2.208569, 2.200225, 2.186411, 2.167263, 2.142970, 2.113775, 2.079971, 2.041903, 2.0 /)
    gx(1:10,6) = (/ 2.682235, 2.679764, 2.672382, 2.660189, 2.643351, 2.622097, 2.596728, 2.567613, 2.535191, 2.5 /) !! DUMMY VALUES (mere interpolation of line above and below)
    gx(1:10,7) = (/ 3.153111, 3.150959, 3.144539, 3.133967, 3.119438, 3.101224, 3.079682, 3.055254, 3.028479, 3.0 /)


    ! get interpolation stencil for v'' between starting values in sweep_angle and kappa direction
    fxx1 = fxx(floor(sweep_angle_degrees/10.)+1,  floor(kappa/0.5)+1)
    fxx2 = fxx(floor(sweep_angle_degrees/10.)+1,  ceiling(kappa/0.5)+1)
    fxx3 = fxx(ceiling(sweep_angle_degrees/10.)+1,floor(kappa/0.5)+1)
    fxx4 = fxx(ceiling(sweep_angle_degrees/10.)+1,ceiling(kappa/0.5)+1)
    ! get interpolation stencil for w' between starting values in sweep_angle and kappa direction
    gx1 = gx(floor(sweep_angle_degrees/10.)+1,  floor(kappa/0.5)+1)
    gx2 = gx(floor(sweep_angle_degrees/10.)+1,  ceiling(kappa/0.5)+1)
    gx3 = gx(ceiling(sweep_angle_degrees/10.)+1,floor(kappa/0.5)+1)
    gx4 = gx(ceiling(sweep_angle_degrees/10.)+1,ceiling(kappa/0.5)+1)


    ! get remainders in kappa and sweep_angle dimension
    dK = mod(kappa,0.5)
    dPhi = mod(sweep_angle_degrees,10.0)


    ! interpolate inside 4 point stencil for v'' 
    fxxLo = ((0.5-dK)*fxx1 + dK*fxx2)/0.5
    fxxHi = ((0.5-dK)*fxx3 + dK*fxx4)/0.5
    fxxGuess = ((10-dPhi)*fxxLo + dPhi*fxxHi)/10.0
    ! interpolate inside 4 point stencil for w' 
    gxLo = ((0.5-dK)*gx1 + dK*gx2)/0.5
    gxHi = ((0.5-dK)*gx3 + dK*gx4)/0.5
    gxGuess = ((10-dPhi)*gxLo + dPhi*gxHi)/10.0


    !   ! --- FOR DEBUGGING PURPOSES ONLY
    !    write(*,*) 'Kappa and Phi are', kappa, ', ', sweep_angle_degrees, '.'
    !    write(*,*) 'Interpolation for f between the four values'
    !    write(*,*) fxx1, ', ', fxx2, ', ', fxx3, ', ', fxx4
    !    write(*,*) 'resulted in initial guess', fxxGuess, '.'
    !    write(*,*)
    !    write(*,*) 'Interpolation for g between the four values'
    !    write(*,*) gx1, ', ', gx2, ', ', gx3, ', ', gx4
    !    write(*,*) 'resulted in initial guess', gxGuess, '.'


  END SUBROUTINE init_value_bl_equation



  SUBROUTINE shoot_v(v,grid,n,sweep_angle,angle_attack,kappa,s,r,rank)

    implicit none

    !!!-----------------------------------------------------------------------
    integer(c_int), intent(in)    :: n             !< number of grid points
    real(c_double), intent(out)   :: v(n,3)        !< OUTPUT: solution vector
    real(c_double), intent(in)    :: grid(n)       !< array of grid points
    real(c_double), intent(inout) :: s,r           !< starting and end value for shooting
    real(c_double), intent(in)    :: sweep_angle   !< sweep angle PHI
    real(c_double), intent(in)    :: angle_attack  !< sweep angle ALPHA
    real(c_double), intent(in)    :: kappa         !< boundary suction KAPPA
    integer(c_int), intent(in)    :: rank          !< rank of processor (for write)
    !!!-----------------------------------------------------------------------
    real(c_double) :: rhs(3)        !value of right-hand side
    real(c_double) :: g(3)          !intermediate value in Runge-Kutta integration
    real(c_double) :: t(3)          !intermediate value of v at x
    real(c_double) :: rk1(3), rk2(3)!Runge-Kutta coefficients
    real(c_double) :: dx            !integration step
    real(c_double) :: x             !independent variable
    !integer(c_int) :: i, j, k    !counters
    integer(c_int) :: i, k    !counters
    real(c_double) :: cosPhi        !cos(sweep_angle)
    real(c_double) :: sinAlpha      !sin(angle_attack)
    real(c_double) :: upperBound    !upper interval boundary for integration
    !!!=====================================================================

    !!set Runge-Kutta coefficients
    rk1(1) = 0.
    rk1(2) = -5./9. 
    rk1(3) = -153./128.
    rk2(1) = 1./3.
    rk2(2) = 15./16.
    rk2(3) = 8./15.

    !!set initial condition: v(0) = kappa, v'(0) = 0, initial guess depends on kappa, phi, alpha
    i = 1
    upperBound = 50.0

    cosPhi = cos(sweep_angle)
    sinAlpha = sin(angle_attack)

    t(1) = -kappa
    t(2) = 0.
    t(3) = s
    x = 0.

    do
      ! -- mjohn 111111
      !! catch values in order to prevent divergence collapse (floating point overflow)
      IF (rank .EQ. 0 .AND. abs(t(1)) .ge. 1e6) THEN
        r = t(2) + cosPhi;
        exit
      end if
      ! -- mjohn 111111

      !! save error and integration value
      r = t(2) + cosPhi;
      IF (x .EQ. grid(i)) THEN
        v(i,1:3) = t
        i = i + 1
      end if

      !! exit condition (compute until full grid is integrated AND 30 is reached)
      IF (i .gt. n .and. x .ge. upperBound) EXIT

      !! standard step size
      dx = 1.e-5

      !! reduce step size accordingly if grid point is approached
      IF (i.le.n .and. x+dx .GT. grid(i)) THEN
        dx = grid(i)-x
        x = grid(i)
      else
        x = x+dx
      END IF

      !!three step runge-kutta
      g = 0.

      do k = 1, 3
        if( x .le. upperBound ) then
          rhs(1) = t(2)
          rhs(2) = t(3)
          rhs(3) = t(3)*t(1) + cosPhi**2 - t(2)**2
        else
          rhs(1) = t(2)
          rhs(2) = 0.
          rhs(3) = 0.
        end if

        g = rk1(k)*g + rhs
        t = t + rk2(k)*dx*g

      end do

    end do

  END SUBROUTINE shoot_v



  SUBROUTINE shoot_w(w,v,grid,n,sweep_angle,sweep_angle_degrees,kappa,sv,sw,r,blThick)

    implicit none

    !!!-----------------------------------------------------------------------
    integer(c_int), intent(in) :: n          !number of grid points
    real(c_double), intent(out):: w(n,2)        !OUTPUT: solution vector
    real(c_double), intent(out):: v(n,3)        !OUTPUT: solution vector
    real(c_double), intent(in) :: grid(n)       !array of grid points
    real(c_double), intent(inout):: sv,sw,r       !starting and end value for shooting
    real(c_double), intent(in)   :: sweep_angle   ! sweep angle PHI
    real(c_double), intent(in)   :: sweep_angle_degrees ! sweep angle PHI
    real(c_double), intent(in)   :: kappa         ! boundary suction KAPPA
    real(c_double), intent(out)  :: blThick       ! blThickness (mathematical, non-dimensionalized)
    !!!-----------------------------------------------------------------------
    real(c_double) :: rhs(5)        !value of right-hand side
    real(c_double) :: g(5)          !intermediate value in Runge-Kutta integration
    real(c_double) :: t(5)          !intermediate value of v at x
    real(c_double) :: rk1(3), rk2(3)!Runge-Kutta coefficients
    real(c_double) :: dx            !integration step
    real(c_double) :: x             !independent variable
    !integer(c_int) :: i, j, k          !counters
    integer(c_int) :: i,  k          !counters
    real(c_double) :: cosPhi        !cos(sweep_angle)
    real(c_double) :: upperBound    !upper interval boundary for integration
    !!!=====================================================================

    !!set Runge-Kutta coefficients
    rk1(1) = 0.
    rk1(2) = -5./9. 
    rk1(3) = -153./128.
    rk2(1) = 1./3.
    rk2(2) = 15./16.
    rk2(3) = 8./15.

    !!set initial condition: v(0) = kappa, v'(0) = 0, initial guess depends on kappa, phi, alpha
    cosPhi = cos(sweep_angle)

    i = 1
    t(1) = -kappa
    t(2) = 0.
    t(3) = sv
    t(4) = 0.
    t(5) = sw
    x = 0.

    !! w needs special treatment for phi = 90° and kappa = 0, since solution becomes trivial
    !! activate VERY large integration domain for large phi and small kappa
    if (kappa .LT. 0.01 .AND. sweep_angle_degrees .GT. 85.0) then
      upperBound = min(2/(kappa+epsilon(upperbound)),200.0);
    else
      upperBound = 50.0;
    end if

    ! initialize boundary layer thickness to 0
    blThick = 0;

    !! numerical integration loop (let coordinate x run from 0 to infinity)
    do

      !!save values
      r = t(4) - 1.;
      IF (x .EQ. grid(i)) THEN
        v(i,1:3) = t(1:3)
        w(i,1:2) = t(4:5)
        i = i + 1
      END IF

      !! exit condition
      IF (i .gt. n .and. x .ge. upperBound) EXIT
      dx = 1.e-5

      IF (i.le.n .and. x+dx .GT. grid(i)) THEN
        dx = grid(i)-x
        x = grid(i)
      ELSE
        x = x+dx
      END IF

      !!three step runge-kutta
      g = 0.
      do k = 1, 3
        if (x < upperBound) then
          rhs(1) = t(2)
          rhs(2) = t(3)
          rhs(3) = t(3)*t(1) + cosPhi**2 - t(2)**2
          rhs(4) = t(5)
          rhs(5) = t(5)*t(1)
        else
          rhs(1) = t(2)
          rhs(2) = 0.
          rhs(3) = 0.
          rhs(4) = 0.
          rhs(5) = 0.
        end if
        g = rk1(k)*g + rhs
        t = t + rk2(k)*dx*g
      end do

      !! carry out first order integral: int_0^\inf (g * (1-g)) dx
      blThick = blThick + (t(4)*(1-t(4))) *dx;

      if(blThick .ge. 0.25 * upperBound) then
        write(*,*) 'WARNING. B.L. THICKNESS GREATER THAN 25 % OF INTEGRATION INTERVAL. BASE FLOW PROFILE MIGHT BE WRONG.'
      end if

    end do

  END SUBROUTINE shoot_w



  !> \brief calculate a laminar swept attachment-line flow
  !!
  !! based on calcbasicflow by obristd, 1999/10/11
  !! extended by mjohn, 2011/11/08
  !! Purpose: calculate a generalized 3D boundary layer flow
  !! with SWEEP ANGLE and boundary suction (wall-normal) KAPPA
  !!
  !! 11/11/28: modified: baseflow now contains "real" profile, not only
  !! integration routine, i.e. multiplication with 1/Re, sin(phi), x3
  !! is performed before ub, vb, wb are returned
  !!
  !! 12/10/05: modified: baseflow may now be non-dimensionalized
  !! according to SHBL formalism (1 d.o.f.), wing formalism (2 d.o.f.)
  !! or novel scaling (1 d.o.f.)
  !!
  !! \param[in] rank
  !! \param[in] kappa
  !! \param[in] sweep_angle_degrees
  !! \param[in] sweep_angle
  !! \param[in] angle_attack
  !! \param[in] n
  !! \param[in] grid
  !! \param[in] ub
  !! \param[in] vb
  !! \param[in] wb
  !! \param[in] blthick
  subroutine calcbasicflow( &
      rank,                 &
      kappa,                &
      sweep_angle_degrees,  &
      sweep_angle,          &
      angle_attack,         &
      n,                    &
      grid,                 &
      ub,                   &
      vb,                   &
      wb,                   &
      blthick )

    implicit none

    integer(c_int), intent(in)    :: rank          
    real(c_double), intent(in)    :: kappa         !< boundary suction KAPPA
    real(c_double), intent(in)    :: sweep_angle_degrees !< sweep angle PHI
    real(c_double), intent(in)    :: sweep_angle   !< sweep angle PHI
    real(c_double), intent(in)    :: angle_attack  
    integer(c_int), intent(in)    :: n             !< number of grid points
    real(c_double), intent(in)    :: grid(n)       !< array of grid points
    real(c_double), intent(inout) :: ub(n)         !< chordwise baseflow profile u
    real(c_double), intent(inout) :: vb(n)         !< wall-normal baseflow profile v
    real(c_double), intent(inout) :: wb(n)         !< spanwise baseflow profile w

    real(c_double),allocatable :: v(:,:)       ! dependent variable v and its derivatives
    real(c_double),allocatable :: w(:,:)       ! dependent variable w and its derivatives

    real(c_double), intent(out) :: blThick

    real(c_double) :: initv1, initv2 !starting values for v equation
    real(c_double) :: initw1, initw2 !starting values for w equation
    real(c_double) :: sv, s0, s1, s1old  !starting values
    real(c_double) :: r0, r1         !values of object function
    real(c_double) :: vxxGuess, wxGuess                ! final guess values (return values)  !-- mjohn 101111

    integer(c_int) :: i

    real(c_double) :: pi 


    allocate(v(n,3),w(n,2))   !! allocate local memory
    pi = 4.*atan(1.)    !!set constants


    !!! ---------------------------------------------------------- INITIALIZE SHOOTING ROUTINE (GET STARTING VALUES) ----------------------------------------------------
    !----------- mjohn 101111 
    ! -- get initial slope/curvature for boundary layer equations as a function of sweep angle and suction velocity
    call init_value_bl_equation( kappa, sweep_angle_degrees, vxxguess, wxguess )

    !!!*** solve equation for v
    ! move to VectorField::initField
    if (rank .EQ. 0 .AND. (kappa.GT.3.0 .OR. kappa.LT.0.0 .OR. sweep_angle_degrees.GT.90.0  .OR. sweep_angle_degrees.LT.0.0 )) then
      WRITE(*,*) 'WARNING: kappa or sweep angle outside allowed interval!! Shooting integration might fail.'
    end if

    initv1 = -vxxGuess*(1+1e-3) !value of ddv for first shot
    initv2 = -vxxGuess*(1-1e-3) !value of ddv for second shot
    initw1 =   wxGuess*(1+1e-3) !value of dw for first shot
    initw2 =   wxGuess*(1-1e-3) !value of dw for second shot
    !----------- mjohn 101111 

    !!! ---------------------------------------------------------- SOLVE SHOOTING INTEGRATION PROBLEM FOR VELOCITY V ----------------------------------------------------
    if (rank.eq.0) then
      WRITE(*,*) ' Solving equation v'''''' = cos^2(sweep_angle) - v''^2  + v v'''' for v'
    end if

    !!set initial values
    s0 = initv1
    s1 = initv2

    !! two initial shots
    call shoot_v(v,grid,n,sweep_angle,angle_attack,kappa,s0,r0,rank)
    call shoot_v(v,grid,n,sweep_angle,angle_attack,kappa,s1,r1,rank)

    i = 1
    !! controlled shooting loop
    do

      !!calculate new initial value (only if residuum is large enough!)
      if( abs(r1) .gt. 1E-13 .and. abs(r0-r1) .gt. 1E-13) then
        s1old = s1
        s1 = (r0*s1 - r1*s0)/(r0-r1)
        s0 = s1old
        r0 = r1
      end if

      call shoot_v(v,grid,n,sweep_angle,angle_attack,kappa,s1,r1,rank)

      if(abs(r1) .le. 1e-10) then
        !!goal reached
        exit
      elseif ((s1 .eq. s0) .or. (r1 .eq. r0)) then
        !!no more progress
        exit
      end if
    end do

    !! save starting value
    sv = s1

    !! print output  
    if (rank.eq.0) then
      !      write(*,*) '   Ended shoot_v successfully. ddy(v)|0 = ', v(1,3)
      write(*,*) ' residuum =', r1
      !      if ( ((initv1+initv2)/2.0) .NE. 0.0) then
      !         write(*,*) 'Relative deviation from initial guess:', (v(1,3)-((initv1+initv2)/2))/((initv1+initv2)/2)
      !      else
      !         write(*,*) 'Total deviation from initial guess:   ', v(1,3)
      !      end if
      !      write(*,*) 
    end if

    !!! ---------------------------------------------------------- SOLVE SHOOTING INTEGRATION PROBLEM FOR VELOCITY W ----------------------------------------------------
    if (rank.eq.0) then
      WRITE(*,*) ' Solving equation w'''' = v w'' for w'
    end if

    !! set initial values
    s0 = initw1
    s1 = initw2

    !! two initial shots
    call shoot_w(w,v,grid,n,sweep_angle,sweep_angle_degrees,kappa,sv,s0,r0,blThick)
    call shoot_w(w,v,grid,n,sweep_angle,sweep_angle_degrees,kappa,sv,s1,r1,blThick)

    !! controlled shooting loop
    do

      !!calculate new initial value
      s1old = s1
      s1 = (r0*s1 - r1*s0)/(r0-r1)
      s0 = s1old
      r0 = r1
      !!shoot
      call shoot_w(w,v,grid,n,sweep_angle,sweep_angle_degrees,kappa,sv,s1,r1,blThick)

      !!check result
      IF(abs(r1) .LE. 1E-10) THEN
        !!goal reached
        EXIT
      ELSEIF ((s1 .EQ. s0) .OR. (r1 .EQ. r0)) THEN
        !!no more progress
        EXIT
      END IF
    END DO

    !! print output
    if (rank.eq.0) then
      !      write(*,*) '   Ended shoot_w successfully. dy(w)|0 = ', w(1,2)
      write(*,*) ' residuum =', r1
      write(*,*)
      !      if (((initw1+initw2)/2).NE.0) then
      !         write(*,*) 'Relative deviation from initial guess:', (w(1,2)-((initw1+initw2)/2))/((initw1+initw2)/2)
      !      else
      !         write(*,*) 'Total deviation from initial guess:   ', w(1,2)
      !      end if
      !      write(*,*) 
    end if


    !!! ---------------------------------------------------------- CONVERT INTEGRATION RESULTS TO VERITABLE VELOCITIES --------------------------------------------------
    ub = -v(:,2)
    vb = v(:,1)
    wb = w(:,1)

    !! deallocate local memory
    deallocate(v,w)

  END SUBROUTINE calcbasicflow



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

    integer(c_int)                       ::  i, j, k


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

    integer(c_int)                       ::  i, j, k


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

    real(c_double) :: pi
    real(c_double) :: Lh
    real(c_double) :: mu
    real(c_double) :: ny
    real(c_double) :: c



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

    real(c_double):: pi
    real(c_double):: Lh
    real(c_double):: mu
    real(c_double):: ny
    real(c_double):: c

    integer(c_int)                ::  i, j, k

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

    real(c_double):: pi
    real(c_double):: Lh
    real(c_double):: mu
    real(c_double):: ny
    real(c_double):: c

    integer(c_int)               ::  i, j, k

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

    real(c_double) :: pi
    real(c_double) :: Lh
    real(c_double) :: mu
    real(c_double) :: ny
    real(c_double) :: c

    integer(c_int)                ::  i, j, k

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

    real(c_double) :: pi

    integer(c_int)                ::  i, j, k

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


    integer(c_int)                      ::  i, j, k

    real(c_double) :: pi


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


    integer(c_int)               ::  i, j, k

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

    real(c_double):: circ
    real(c_double):: rad
    real(c_double):: r
    real(c_double):: pi


    integer(c_int)               ::  i, j, k

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

    real(c_double):: sig

    integer(c_int)               ::  i, j, k

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
    real(c_double):: sig

    integer(c_int)               ::  i, j, k

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

    real(c_double):: h
    real(c_double):: h2

    integer(c_int)                ::  i, j, k

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

    real(c_double):: h
    real(c_double):: h2

    integer(c_int)               ::  i, j, k

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
      xm,ym, rad,sca,         &
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
    real(c_double), intent(in)    :: sca


    real(c_double)   :: dr

    integer(c_int) ::  i, j, k

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
          phiU(i,j,k) = sca*distance2ib( x1u(i),x2p(j),x3p(k),xm,ym,rad,dr )
        end do
      end do
    end do

    do k = SV(3), NV(3)
      do j = SV(2), NV(2)
        do i = SV(1), NV(1)
          phiV(i,j,k) = sca*distance2ib( x1p(i),x2v(j),x3p(k),xm,ym,rad,dr )
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
      phiU,                           &
      phiV,                           &
      phiW ) bind ( c, name='VF_init_RotatingDisc' )

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


    integer(c_int)               ::  i, j, k

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



  !> \brief init swept Hiemenz boundary layer
  !! 
  !! \note - baseflow_global(*,1) need to run from 0 to M1, i.e. full axis, on u grid, since interpolated to p grid by init_BC
  !!       - baseflow_global(*,2) and (*,3) are computed on p grid, since they start exactly on p(1) wall
  !!
  !! difficulty: - full u grid serves also negative values below the wall, where shooting integration is invalid
  !! solution:   - compute velocity for y = 0 and y > 0 on u grid, then extrapolate from 0 to first (negative) u value
  !!             - perform extrapolation not on all ranks, but only those which touch the ground
  !!             - since only one value below the wall is relevant for interpolation, extrapolation formula reads:
  !!                      base(0,1)|_u = ( base(0,1)|_p - sum_0^d1U{cIup(j,1)*base(1+j,1)|_u} ) / cIup(-1,1)
  subroutine VF_init_SHBF(  &
      rank,                 &
      iShift,               &
      IB1,                  &
      M,                    &
      N,                    &
      bL,bU,                &
      dL,dU,                &
      SU,NU,                &
      SV,NV,                &
      SW,NW,                &
      y1p,                  &
      y1u,                  &
      x3w,                  &
      cIup,                 &
      Re,                   &
      nonDim,               &
      kappa,                &
      sweep_angle_degrees,  &
      sweep_angle,          &
      angle_attack,         &
      velU,                 &
      velV,                 &
      velW ) bind ( c, name='VF_init_SHBF' )

    implicit none

    integer(c_int), intent(in)     :: rank
    integer(c_int), intent(in)     :: iShift
    integer(c_int), intent(in)     :: IB1

    integer(c_int), intent(in)     :: M(3)
    integer(c_int), intent(in)     :: N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: dL(3)
    integer(c_int), intent(in)     :: dU(3)

    integer(c_int), intent(in)     :: SU(3)
    integer(c_int), intent(in)     :: NU(3)

    integer(c_int), intent(in)     :: SV(3)
    integer(c_int), intent(in)     :: NV(3)

    integer(c_int), intent(in)     :: SW(3)
    integer(c_int), intent(in)     :: NW(3)

    real(c_double), intent(in)     :: y1p( 1:M(1) )

    real(c_double), intent(in)     :: y1u( 0:M(1) )

    !real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
    !real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(in)     :: x3w( bL(3):(N(3)+bU(3)) )

    real(c_double), intent(in)     :: cIup( dU(1):dL(1), 0:N(1) )

    real(c_double), intent(in)     :: Re              
    integer(c_int), intent(in)     :: nonDim
    real(c_double), intent(in)     :: kappa               !< Properties of swept Hiemenz flow
    real(c_double), intent(in)     :: sweep_angle_degrees !< Properties of swept Hiemenz flow
    real(c_double), intent(in)     :: sweep_angle         !< Properties of swept Hiemenz flow
    real(c_double), intent(in)     :: angle_attack        !< Properties of swept Hiemenz flow

    real(c_double),  intent(inout) :: velU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(inout) :: velV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(inout) :: velW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i, ii, j, k

    !! mjohn 051012 - variables required for changing the nondimensionalization
    real(c_double)                ::  blThick
    !! mjohn 051012

    real(c_double)                ::  y1u_temp(0:M(1))


    real(c_double)                ::  baseflow_global(bL(1):(M(1)+bU(1)),1:3)
    real(c_double)                ::  baseflow       (bL(1):(N(1)+bU(1)),1:3)
    !real(c_double)                ::  baseflow1p     (bL(1):(N(1)+bU(1)))


    !---- mjohn 120207  ------------------------------------------------ SWEPT HIEMENZ BASE FLOW -----------------------------------------------------------------------

    ! initialize variables
    velU = 0.
    velV = 0.
    velW = 0.

    baseflow_global = 0.
    baseflow        = 0.

    ! ------------------------- get two tangential base flow components (staggered grid! need information on p grid)
    ! get global basic flow (exactly FULL y axis, over all ranks) on yp grid, store in baseflow_global, for two tangent coordinates
    call calcbasicflow( rank, kappa, sweep_angle_degrees, sweep_angle, angle_attack, M(1),y1p(1),baseflow_global(1,3),baseflow_global(1,1),baseflow_global(1,2), blThick)
    ! --- mjohn 051012
    ! set different velocity fields if integration with Hiemenz is performed
    ! integration routines however are unaltered (see below), because variables are overwritten (see usr_config.f90)
    select case (nonDim)
    case(0)
      if(rank .eq. 0) then
        write(*,'(a)') 'classical non-dimensionalization applied'
      end if
      ! traditinoal SHBL case
      baseflow(bL(1):(N(1)+bU(1)),2) = baseflow_global((bL(1)+iShift):(N(1)+bU(1)+iShift),2)
    case default
      if(rank .eq. 0) then
        write(*,'(a)') 'novel non-dimensionalization applied'
      end if
      ! 2 degrees of freedom and novel formalism (cases 1 and 2)
      baseflow(bL(1):(N(1)+bU(1)),2) = baseflow_global((bL(1)+iShift):(N(1)+bU(1)+iShift),2) * sin(sweep_angle)
    end select
    baseflow(bL(1):(N(1)+bU(1)),3) = baseflow_global((bL(1)+iShift):(N(1)+bU(1)+iShift),3) / Re

    ! ------------------------- get wall-normal base flow component (staggered grid! need information on u grid)
    ! get global basic flow (ENTIRE y axis, excluding points beyond boundaries, over all ranks) on u grid, store in baseflow_global, for wall-normal coordinate
    y1u_temp(1:M(1))  = y1u(1:M(1))
    y1u_temp(0     )  = 0.

    call calcbasicflow(rank, kappa, sweep_angle_degrees, sweep_angle, angle_attack,M(1)+1,y1u_temp(0),baseflow_global(0,3),baseflow_global(0,1),baseflow_global(0,2), blThick)

    if( 1==IB1 ) then
      write(*,*) "blabla"
      do ii = 0, dU(1)
        write(*,*) ii, cIup(ii,1),baseflow_global(0,1)
        baseflow_global(0,1) = baseflow_global(0,1) - cIup(ii,1)*baseflow_global(1+ii,1)
      end do
      baseflow_global(0,1) = baseflow_global(0,1) / cIup(-1,1)
    end if
    baseflow(bL(1):(N(1)+bU(1)),1) = baseflow_global((bL(1)+iShift):(N(1)+bU(1)+iShift),1) / Re

    write(*,*) ii, cIup(-1,1),baseflow_global(0,1)
    ! TEST!!! - debugging purposes only
    write(*,*) "distance, velocity (both in x1-direction)"
    DO i = SU(1), NU(1)
      write(*,*)  baseflow(i,1)
    end do
    write(*,*)
    write(*,*) "distance, velocity (both in x2-direction)"
    DO i = SV(1), NV(1)
      write(*,*)  baseflow(i,2)
    end do
    write(*,*)
    write(*,*) "distance, velocity (both in x3-direction)"
    DO k = SW(1), NW(1)
      write(*,*) , baseflow(i,3)*x3w(k)
    end do

    do k = SU(3), NU(3)
      do j = SU(2), NU(2)
        do i = SU(1), NU(1)
          velU(i,j,k) = velU(i,j,k) + baseflow(i,1)
        END DO
      END DO
    END DO
    do k = SV(3), NV(3)
      do j = SV(2), NV(2)
        do i = SV(1), NV(1)
          velV(i,j,k) = velV(i,j,k) + baseflow(i,2)
        END DO
      END DO
    END DO
    do k = SW(3), NW(3)
      do j = SW(2), NW(2)
        do i = SW(1), NW(1)
          velW(i,j,k) = velW(i,j,k) + baseflow(i,3)*x3w(k)
        END DO
      END DO
    END DO


  end subroutine VF_init_SHBF




  function eddy(x1,x2,dir) result(fn_val)

    implicit none

    real(c_double), intent(in   ) ::  x1
    real(c_double), intent(in   ) ::  x2
    integer(c_int), intent(in   ) ::  dir

    real(c_double)                ::  length, orient, pi
    real(c_double)                ::  fn_val


    pi = 4.*atan(1.)    !!set constants

    if (1 == 2) then

      length = sqrt(x1**2+x2**2)
      orient = atan2(x2,x1)

      fn_val = sin(pi*erf(1.*length))

      !if (dir == 1) fn_val = -fn_val*sin(orient)
      !if (dir == 2) fn_val =  fn_val*cos(orient)
      !if (dir == 3) fn_val =  fn_val*cos(orient)
      if (dir == 1) fn_val =  fn_val*sin(orient)
      if (dir == 2) fn_val = -fn_val*cos(orient)
      if (dir == 3) fn_val = -fn_val*cos(orient)

    else

      if (x1 .ge. -1. .and. x1 .le. 1. .and. x2 .ge. -1. .and. x2 .le. 1.) then
        !! if (dir == 1) fn_val = -sin(pi*x1)*(1.+cos(pi*x2))/2. ! "quatropol" mit sich teilender instabilitaet (eine verbleibt in der symmetriefläche)
        !! if (dir == 2) fn_val =  sin(pi*x2)*(1.+cos(pi*x1))/2.
        !! if (dir == 3) fn_val =  sin(pi*x2)*(1.+cos(pi*x1))/2.
        ! if (dir == 1) fn_val = -sin(pi*x2)*(1.+cos(pi*x1))/2.
        ! if (dir == 2) fn_val =  sin(pi*x1)*(1.+cos(pi*x2))/2.
        ! if (dir == 3) fn_val =  sin(pi*x1)*(1.+cos(pi*x2))/2.
        if (dir == 1) fn_val =  sin(pi*x2)*(1.+cos(pi*x1))/2.
        if (dir == 2) fn_val = -sin(pi*x1)*(1.+cos(pi*x2))/2.
        if (dir == 3) fn_val = -sin(pi*x1)*(1.+cos(pi*x2))/2.
      else
        fn_val = 0.
      end if

      !if (x1 .lt. -1. .or. x1 .gt. 1.) fn_val = 0.
      !if (x2 .lt. -1. .or. x2 .gt. 1.) fn_val = 0.

    end if

    return


  end function eddy

  !>  0 generic GH
  !!  1 counter-rotating primary vortices, even secondary vortices
  !!  2 counter-rotating primary vortices, odd secondary vortices
  !!  3 pair of antisymmetric streaks 
  !!  4 array of periodic vortices 
  !!  5 line of exaggerated spanwise velocity
  !!  6 Goertler-Haemmerlin disturbance in w component (spanwise)
  !!  7 generic noise (OU process)
  !!  8 pair of counter-rotating vortices with d/dz (init cond)
  !!  9 array of counter-rotating vortices with d/dz (init cond)
  subroutine VF_init_Dist(  &
      rank,                 &
      N,                    &
      bL,bU,                &
      SU,NU,                &
      SV,NV,                &
      SW,NW,                &
      BC_3L_global,         &
      x1u,                  &
      x1p,                  &
      x2p,                  &
      x3w,                  &
      x3p,                  &
      dist_type,            &
      vortex_ampli_prim,    &
      vortex_x1pos,         &
      vortex_x3pos,         &
      vortex_radius,        &
      vortex_band,          &
      velU,velV,velW ) bind ( c, name='VF_init_Dist' )

    implicit none

    integer(c_int), intent(in)     :: rank

    integer(c_int), intent(in)     :: N(3)

    integer(c_int), intent(in)     :: bL(3)
    integer(c_int), intent(in)     :: bU(3)

    integer(c_int), intent(in)     :: SU(3)
    integer(c_int), intent(in)     :: NU(3)

    integer(c_int), intent(in)     :: SV(3)
    integer(c_int), intent(in)     :: NV(3)

    integer(c_int), intent(in)     :: SW(3)
    integer(c_int), intent(in)     :: NW(3)

    integer(c_int), intent(in)     :: BC_3L_global


    real(c_double), intent(in)     :: x1u( bL(1):(N(1)+bU(1)) )
    real(c_double), intent(in)     :: x1p( bL(1):(N(1)+bU(1)) )
    real(c_double), intent(in)     :: x2p( bL(2):(N(2)+bU(2)) )

    real(c_double), intent(in)     :: x3p( bL(3):(N(3)+bU(3)) )
    real(c_double), intent(in)     :: x3w( bL(3):(N(3)+bU(3)) )


    integer(c_int), intent(in)     :: dist_type

    real(c_double), intent(in)     :: vortex_ampli_prim
    real(c_double), intent(in)     :: vortex_x1pos, vortex_x3pos
    real(c_double), intent(in)     :: vortex_radius, vortex_band

    real(c_double),  intent(inout) :: velU(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(inout) :: velV(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))
    real(c_double),  intent(inout) :: velW(bL(1):(N(1)+bU(1)),bL(2):(N(2)+bU(2)),bL(3):(N(3)+bU(3)))

    integer(c_int)                ::  i, ii, j, k

    !! mjohn 051012 - variables required for changing the nondimensionalization
    integer(c_int)                ::  No, Ne, aborted, merror
    real(c_double)                ::  pi
    !! mjohn 051012


    velU( SU(1):NU(1),SU(2):NU(2),SU(3):NU(3) ) = 0.
    velV( SV(1):NV(1),SV(2):NV(2),SV(3):NV(3) ) = 0.
    velW( SW(1):NW(1),SW(2):NW(2),SW(3):NW(3) ) = 0.


    aborted = 0

    pi = 4.*atan(1.)    !!set constants


    ! wall-normal v component
    do k = SU(3), NU(3)
      do j = SU(2), NU(2)
        do i = SU(1), NU(1)
          select case (dist_type)
          case (0)                 !----- mjohn 301111 - generic inflow
            aborted = 1
            ! do nothing
          case (1)                 !----- mjohn 301111 - pair of primary vortices, even secondary forcing
            if( BC_3L_global == -2 ) then
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)
            ELSE
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)
              velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)+vortex_x3pos)/vortex_radius,1)
            END IF
          CASE (2)                 !----- mjohn 301111 - odd secondary forcing
            IF (BC_3L_global == -2) THEN
              WRITE(*,*) "ERROR! Wrong disturbance selected. Symmetric setup does not allow for odd secondary forcing. Simulation aborted."
              CALL MPI_FINALIZE(merror)
            ELSE
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)
              velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)+vortex_x3pos)/vortex_radius,1)
            END IF
          CASE (3)                 !----- mjohn 301111 - mixed secondary forcing
            IF (BC_3L_global == -2) THEN
              WRITE(*,*) "ERROR! Wrong disturbance selected. Symmetric setup does not allow for mixed secondary forcing. Simulation aborted."
              CALL MPI_FINALIZE(merror)
            ELSE
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)+vortex_x3pos)/vortex_radius,1)
              velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)             )/vortex_radius,1)
            END IF
          CASE (4)                 !----- mjohn 260912 - array of primary vortices, no secondary forcing (without distinction between symmetric and non-symmetric case / without IF(BC_3L_..) statement)
            ! get number of even and odd vortices to be initialized within domain
            Ne =   floor((vortex_band/vortex_x3pos +1)/4)
            No = ceiling((vortex_band/vortex_x3pos -1)/4)
            Ne = min(Ne,No)
            No = min(Ne,No)
            ! initialize No odd vortices at positions (1*x3pos, 5*x3pos, ...) and Ne even vortices (with opposite sign) at positions (3*x3pos, 7*x3pos, ...)
            do ii = -ceiling((No-1)/2.0), floor((No-1)/2.0)
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k) - vortex_x3pos - 4*ii*vortex_x3pos)/vortex_radius,1)
              velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k) + vortex_x3pos + 4*ii*vortex_x3pos)/vortex_radius,1)
            end do
          CASE (5)
            aborted = 1
          CASE (6)
            aborted = 1
          CASE (7)
            aborted = 1
          CASE (8)                 !----- mjohn 300513 - pair of primary vortices with modulation in 2 (z) direction
            velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
            velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)+vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
          CASE (9)                 !----- mjohn 310513 - array of primary vortices with modulation in 2 (z) direction
            ! get number of even and odd vortices to be initialized within domain
            Ne =   floor((vortex_band/vortex_x3pos +1)/4)
            No = ceiling((vortex_band/vortex_x3pos -1)/4)
            Ne = min(Ne,No)
            No = min(Ne,No)
            ! initialize No odd vortices at positions (1*x3pos, 5*x3pos, ...) and Ne even vortices (with opposite sign) at positions (3*x3pos, 7*x3pos, ...)
            do ii = -ceiling((No-1)/2.0), floor((No-1)/2.0)
              velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k) - vortex_x3pos - 4*ii*vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
              velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k) + vortex_x3pos + 4*ii*vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
            end do

          case default
            aborted = 1
            ! do nothing
          end select
        end do
      end do
    end do


    ! chordwise u component
    do k = SW(3), NW(3)
      do j = SW(2), NW(2)
        do i = SW(1), NW(1)
          select case (dist_type)
          case (0)                 !----- mjohn 301111 - generic inflow conditions of gh-type
            aborted = 1
            ! do nothing
          case (1)                 !----- mjohn 301111 - even secondary forcing
            if (BC_3L_global == -2) then
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)-vortex_x3pos)/vortex_radius,3)
            else
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)-vortex_x3pos)/vortex_radius,3)
              velW(i,j,k) = velW(i,j,k) - vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)+vortex_x3pos)/vortex_radius,3)
            end if
          case (2)                 !----- mjohn 301111 - odd secondary forcing
            IF (BC_3L_global == -2) THEN
              WRITE(*,*) "ERROR! Wrong disturbance selected. Symmetric setup does not allow for odd secondary forcing. Simulation aborted."
              CALL MPI_FINALIZE(merror)
            else
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)-vortex_x3pos)/vortex_radius,3)
              velW(i,j,k) = velW(i,j,k) - vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)+vortex_x3pos)/vortex_radius,3)
            end if
          case (3)                 !----- mjohn 301111 - mixed secondary forcing
            IF (BC_3L_global == -2) THEN
              WRITE(*,*) "ERROR! Wrong disturbance selected. Symmetric setup does not allow for mixed secondary forcing. Simulation aborted."
              CALL MPI_FINALIZE(merror)
            else
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)-vortex_x3pos)/vortex_radius,3)
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)+vortex_x3pos)/vortex_radius,3)
              velW(i,j,k) = velW(i,j,k) - vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k)             )/vortex_radius,3)
            END IF
          case (4)                 !----- mjohn 260912 - array of primary vortices, no secondary forcing (without distinction between symmetric and non-symmetric case / without IF(BC_3L_..) statement)
            ! get number of even and odd vortices to be initialized within domain
            Ne =   floor((vortex_band/vortex_x3pos +1)/4)
            No = ceiling((vortex_band/vortex_x3pos -1)/4)
            Ne = min(Ne,No)
            No = min(Ne,No)
            ! initialize No odd vortices at positions (1*x3pos, 5*x3pos, ...) and Ne even (with opposite sign) vortices at positions (3*x3pos, 7*x3pos, ...)
            do ii = -ceiling((No-1)/2.0), floor((No-1)/2.0)
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k) - vortex_x3pos - 4*ii*vortex_x3pos)/vortex_radius,3)
              velW(i,j,k) = velW(i,j,k) - vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k) + vortex_x3pos + 4*ii*vortex_x3pos)/vortex_radius,3)
            end do
          case (5)
            aborted = 1
          case (6)
            aborted = 1
          case (7)
            aborted = 1
          case (8)                 !----- mjohn 300513 - pair of primary vortices with modulation in 2 (z) direction
            velU(i,j,k) = velU(i,j,k) + vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)-vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
            velU(i,j,k) = velU(i,j,k) - vortex_ampli_prim*eddy((x1u(i)-vortex_x1pos)/vortex_radius,(x3p(k)+vortex_x3pos)/vortex_radius,1)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
          case (9)                 !----- mjohn 310513 - array of primary vortices with modulation in 2 (z) direction
            ! get number of even and odd vortices to be initialized within domain
            Ne =   floor((vortex_band/vortex_x3pos +1)/4)
            No = ceiling((vortex_band/vortex_x3pos -1)/4)
            Ne = min(Ne,No)
            No = min(Ne,No)
            ! initialize No odd vortices at positions (1*x3pos, 5*x3pos, ...) and Ne even (with opposite sign) vortices at positions (3*x3pos, 7*x3pos, ...)
            do ii = -ceiling((No-1)/2.0), floor((No-1)/2.0)
              velW(i,j,k) = velW(i,j,k) + vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k) - vortex_x3pos - 4*ii*vortex_x3pos)/vortex_radius,3)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
              velW(i,j,k) = velW(i,j,k) - vortex_ampli_prim*eddy((x1p(i)-vortex_x1pos)/vortex_radius,(x3w(k) + vortex_x3pos + 4*ii*vortex_x3pos)/vortex_radius,3)*(0.75+0.25*sin(2*pi*x2p(j)/(4*vortex_radius)))
            end do
          case default
            aborted = 1
            ! do nothing
          end select
        end do
      end do
    end do

    !----- mjohn 260912
    IF(rank == 0 .AND. aborted .eq. 1) THEN
      WRITE(*,*) "WARNING! No valid initial disturbance defined. Initial conditions equal to baseflow only."
    END IF

  end subroutine VF_init_Dist

end module cmod_VectorField
