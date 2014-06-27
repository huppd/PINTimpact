!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief user defiend help function, useful for user initial/boundary conditions
MODULE usr_func
  
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  
  
  PRIVATE
  
  PUBLIC interface, erf, atanh
!  PUBLIC parable1
!  PUBLIC distance2ib
!  public vel_ib
  PUBLIC char_ib
  PUBLIC coord_tanh
  PUBLIC init_particles
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  !--- user-specific subroutines and functions ---
  !
  
  
  !> @brief user function, that return some shear profile
  !!
  !! \f[ \mathtt{fn\_val} \gets \left\{
  !!    \begin{array}{cl} 1 & \mathrm{if } \mathtt{ xx} > 3. \\
  !!                      0 & \mathrm{if } \mathtt{ xx} < -3 \\
  !!       \frac12(1+\mathrm{erf}(\sqrt\pi \mathtt{ x}) & \mathtt{else} \end{array} \right. \f]
  !! @param xx inputvale
  !! @return fn_value return value
  FUNCTION interface(xx) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN   ) ::  xx
  REAL                   ::  fn_val
  REAL                   ::  cutoff
  
  
  ! underflow ...
  cutoff = 3.
  
  
  IF      (xx .GT.  cutoff) THEN
     fn_val = 1.
  ELSE IF (xx .LT. -cutoff) THEN
     fn_val = 0.
  ELSE
     fn_val = 0.5*(1.+erf(SQRT(pi)*xx))
  END IF
  
  RETURN
  
  
  END FUNCTION interface
  
!  !> \brief parable function
!  !! returns parable
!  !! \author huppd
!  function parable1(xx) result(fn_val)
!
!  implicit none
!
!  real, intent(in)  :: xx
!  real              :: fn_val
!
!  fn_val = xx*xx
!
!  end function parable1





  
  
  

  !> \brief sets xx and dx so that one gets a nice grid streching
  !! \param[in] Lmax Li
  !! \param[in] iimax Mi
  !! \param[in] ii index
  !! \param[out] xx
  !! \param[out] dx
  SUBROUTINE coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
  ! (sample subroutine)
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)    ::  Lmax
  REAL   , INTENT(IN)    ::  iimax
  
  REAL   , INTENT(IN)    ::  ii0L
  REAL   , INTENT(IN)    ::  ii0U
  
  REAL   , INTENT(IN)    ::  ii
  REAL   , INTENT(OUT)   ::  xx
  REAL   , INTENT(OUT)   ::  dx
  
  REAL                   ::  yy, cmin, cmax
  
  
  IF (ii0U == ii0L) THEN
     ! equidistant grid:
     xx = (ii-1.)*Lmax/(iimax-1.)
     dx =         Lmax/(iimax-1.)
  ELSE
     cmax =  TANH(ii0U)
     cmin = -TANH(ii0L)
     
     ! coordinate transformation (index i=1. is the origin):
     ! y = (i-1.)/(imax-1.)*(cmax-cmin) + cmin
     yy = (ii-1.)/(iimax-1.)*(cmax-cmin) + cmin
     
     ! mapping funktion f(yy)
     ! x = L * (f(y)-f(cmin)) / (f(cmax)-f(cmin))
     xx = Lmax * (atanh(yy)-atanh(cmin)) / (ii0U+ii0L)
     
     ! dx/di = L / (f(cmax)-f(cmin)) * dy/di                 * df(y(i))/dy
     !       = L / (f(cmax)-f(cmin)) * (cmax-cmin)/(imax-1.) * df(y(i))/dy
     dx = Lmax / (atanh(cmax)-atanh(cmin)) * (cmax-cmin)/(iimax-1.) * 1./(1.-yy**2)
  END IF
  
  
  END SUBROUTINE coord_tanh
  
  
  
  
  
  
  
  
  
  
  

  FUNCTION atanh(x) RESULT(fn_val)
  ! (sample function)
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)  :: x
  REAL                 :: fn_val
  
  
  IF (x == 0.) THEN
     fn_val = 0.
     RETURN
  END IF
  
  fn_val = 0.5* LOG((1.+x)/(1.-x))
  RETURN
  
  
  END FUNCTION atanh
  
  
  
  
  
  
  
  
  
  
  
  !> \brief evaluation of the real error function
  FUNCTION erf(x) RESULT(fn_val)
  ! (sample function)
  
  !-----------------------------------------------------------------------
  !             EVALUATION OF THE REAL ERROR FUNCTION
  ! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
  ! Mathematics Library (1993 version).
  ! Adapted by Alan.Miller @ vic.cmis.csiro.au
  !-----------------------------------------------------------------------
  
  IMPLICIT NONE
  
  REAL   , INTENT(IN)  ::  x
  REAL                 ::  fn_val
  
  REAL   , PARAMETER   ::  c = .564189583547756, one = 1.0, half = 0.5, zero = 0.0
  REAL   , PARAMETER   ::  &
           a(5) = (/  7.71058495001320D-05, -1.33733772997339D-03, 3.23076579225834D-02, 4.79137145607681D-02, 1.28379167095513D-01 /),  &
           b(3) = (/  3.01048631703895D-03,  5.38971687740286D-02, 3.75795757275549D-01 /),  &
           p(8) = (/ -1.36864857382717D-07,  5.64195517478974D-01, 7.21175825088309D+00, 4.31622272220567D+01, 1.52989285046940D+02, 3.39320816734344D+02, 4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00,  1.27827273196294D+01, 7.70001529352295D+01, 2.77585444743988D+02, 6.38980264465631D+02, 9.31354094850610D+02, 7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00,  2.62370141675169D+01, 2.13688200555087D+01, 4.65807828718470D+00, 2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01,  1.87114811799590D+02, 9.90191814623914D+01, 1.80124575948747D+01 /)
  REAL                 ::  ax, bot, t, top, x2
  
  
  ax = ABS(x)
  
  IF (ax <= half) THEN
     t = x*x
     top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
     bot = ((b(1)*t + b(2))*t + b(3))*t + one
     fn_val = x*(top/bot)
     RETURN
  END IF
  
  IF (ax <= 4.0) THEN
     top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
           + p(6))*ax + p(7))*ax + p(8)
     bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
           + q(6))*ax + q(7))*ax + q(8)
     fn_val = half + (half - EXP(-x*x)*top/bot)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  IF (ax < 5.8) THEN
     x2 = x*x
     t = one / x2
     top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
     bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
     fn_val = (c - top/(x2*bot)) / ax
     fn_val = half + (half - EXP(-x2)*fn_val)
     IF (x < zero) fn_val = -fn_val
     RETURN
  END IF
  
  fn_val = SIGN(one, x)
  RETURN
  
  
  END FUNCTION erf
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE init_particles(n_part_all,particles)
  ! (sample function)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  n_part_all
  REAL   , INTENT(OUT  ) ::  particles(1:n_args,1:n_part_max)
  REAL                   ::  xx1, xx2, xx3
  REAL                   ::  LL1, LL2, LL3
  REAL                   ::  mult
  REAL                   ::  maxfrac
  INTEGER                ::  n_part_global
  INTEGER                ::  too_many, too_many_global
  INTEGER                ::  i, j, m
  
  
  !--- ratio of maximum initial number and absolute maximum number of particles per subdomain ----------------
  maxfrac = 0.75
  
  
  !===========================================================================================================
  IF (L1c .GT. L1) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! L1c > L1!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (L2c .GT. L2) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! L2c > L2!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  IF (L3c .GT. L3) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! L3c > L3!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  END IF
  !===========================================================================================================
  
  
  ! particle arguments:
  !   particles(1:3,1:n_part) == particle position
  !   particles(4:6,1:n_part) == particle velocity
  !   particles(7:9,1:n_part) == particle acceleration
  !   particles(10 ,1:n_part) == particle feedback force
  !   particles(11 ,1:n_part) == particle settling velocity
  !   particles(12 ,1:n_part) == particle Stokes number
  !   particles(13 ,1:n_part) == particle number
  !   particles(14 ,1:n_part) == particle group
  !   particles(15:,1:n_part) == other particle properties
  
  n_part    = 0
  particles = 0.
  
  LL1 = 0.
  LL2 = 0.
  LL3 = 0.
  
  IF (dimens == 3) THEN
     IF (x1p(S1p) .LT. L1c-y1_origin .AND.  &
      &  x2p(S2p) .LT. L2c-y2_origin .AND.  &
      &  x3p(S3p) .LT. L3c-y3_origin) THEN
        
        IF ((BC_1U == 0 .OR. BC_1L == -1) .AND. x1p(N1p+1) .LT. L1c-y1_origin) THEN
           LL1 = x1p(N1p+1)   -x1p(S1p)
        ELSE
           LL1 = L1c-y1_origin-x1p(S1p)
        END IF
        
        IF ((BC_2U == 0 .OR. BC_2L == -1) .AND. x2p(N2p+1) .LT. L2c-y2_origin) THEN
           LL2 = x2p(N2p+1)   -x2p(S2p)
        ELSE
           LL2 = L2c-y2_origin-x2p(S2p)
        END IF
        
        IF ((BC_3U == 0 .OR. BC_3L == -1) .AND. x3p(N3p+1) .LT. L3c-y3_origin) THEN
           LL3 = x3p(N3p+1)   -x3p(S3p)
        ELSE
           LL3 = L3c-y3_origin-x3p(S3p)
        END IF
     END IF
     n_part = NINT(LL1*LL2*LL3/(L1c*L2c*L3c)*REAL(n_part_all))
  ELSE
     IF (x1p(S1p) .LT. L1c-y1_origin .AND.  &
      &  x2p(S2p) .LT. L2c-y2_origin) THEN
        
        IF ((BC_1U == 0 .OR. BC_1L == -1) .AND. x1p(N1p+1) .LT. L1c-y1_origin) THEN
           LL1 = x1p(N1p+1)   -x1p(S1p)
        ELSE
           LL1 = L1c-y1_origin-x1p(S1p)
        END IF
        
        IF ((BC_2U == 0 .OR. BC_2L == -1) .AND. x2p(N2p+1) .LT. L2c-y2_origin) THEN
           LL2 = x2p(N2p+1)   -x2p(S2p)
        ELSE
           LL2 = L2c-y2_origin-x2p(S2p)
        END IF
     END IF
     n_part = NINT(LL1*LL2    /(L1c*L2c    )*REAL(n_part_all))
  END IF
  
  too_many = 0
  IF (REAL(n_part) .GT. maxfrac*REAL(n_part_max)) too_many = 1
  
  CALL MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  IF (too_many_global .GT. 0) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Too many particles on at least one processor!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  
  IF (n_part .GT. 0) THEN
     
     CALL RANDOM_NUMBER(particles(1:n_args,1:n_part))
     
     DO m = 1, n_part
        
        !--- particle position -------------------------------------------------------------------------------
                         particles(1,m) = x1p(S1p) + LL1*particles(1,m)
                         particles(2,m) = x2p(S2p) + LL2*particles(2,m)
        IF (dimens == 3) particles(3,m) = x3p(S3p) + LL3*particles(3,m)
        
        !--- particle speed ----------------------------------------------------------------------------------
                         particles(4,m) = 0.
                         particles(5,m) = 0.
        IF (dimens == 3) particles(6,m) = 0.
        
        DO i = 10, n_args
           !--- Richardson number ----------------------------------------------------------------------------
           IF (i == 10) particles(i,m) = Rip(1)
           
           !--- settling velocity ----------------------------------------------------------------------------
           IF (i == 11) particles(i,m) = usp(1)
           
           !--- Stokes number --------------------------------------------------------------------------------
           IF (i == 12) particles(i,m) = Stp(1)
           
           !--- particle name --------------------------------------------------------------------------------
           IF (i == 13) particles(i,m) = 0. ! REAL(m+n_part_max-n_part_all) ! TEST!!! needs to be adapted ...
           
           !--- Gruppen werden nach Position von Argument 1 im Intervall [0:1] eingeteilt --------------------
           IF (i == 14) THEN
              DO j = 1, n_groups
                 IF (particles(1,m) .GE. REAL(j-1)/REAL(n_groups)*L1 .AND. particles(1,m) .LE. REAL(j)/REAL(n_groups)*L1) THEN
                    particles(i,m) = REAL(j)
                 END IF
              END DO
           END IF
        END DO
        
     END DO
     
  END IF
  
  
  !===========================================================================================================
  !=== adjust particle weight ================================================================================
  !===========================================================================================================
  CALL MPI_ALLREDUCE(n_part,n_part_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  IF (n_part /= 0) WRITE(*,'(a,i10,a,3i4)') 'n_part        = ', n_part, ', iB(1:3)', iB(1:3,1)
  IF (rank   == 0) WRITE(*,'(a,i10      )') 'n_part_global = ', n_part_global
  
  IF (n_part .GT. 0) THEN
    IF (dimens == 3) THEN
       mult = L1c*L2c*L3c / REAL(n_part_global)
    ELSE
       mult = L1c*L2c     / REAL(n_part_global)
    END IF
    
    i = 10
    DO m = 1, n_part
       particles(i,m) = particles(i,m) * mult
    END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE init_particles
  
END MODULE usr_func
