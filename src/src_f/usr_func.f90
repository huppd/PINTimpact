!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief user defiend help function, useful for user initial/boundary conditions
module usr_func
  
  
  use mod_dims
  use mod_vars
  use usr_vars
  
  
  private
  
  public interface, erf, atanh
!  PUBLIC parable1
!  PUBLIC distance2ib
!  public vel_ib
  public char_ib
  public coord_tanh
  public init_particles
  
  
  INCLUDE 'mpif.h'
  
  contains
  
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
  function interface(xx) result(fn_val)
  ! (sample function)
  
  implicit none
  
  real   , intent(in   ) ::  xx
  real                   ::  fn_val
  real                   ::  cutoff
  
  
  ! underflow ...
  cutoff = 3.
  
  
  if      (xx .gt.  cutoff) then
     fn_val = 1.
  else if (xx .lt. -cutoff) then
     fn_val = 0.
  else
     fn_val = 0.5*(1.+erf(SQRT(pi)*xx))
  end if
  
  return
  
  
  end function interface
  
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
  subroutine coord_tanh(Lmax,iimax,ii0L,ii0U,ii,xx,dx)
  ! (sample subroutine)
  
  implicit none
  
  real   , intent(in)    ::  Lmax
  real   , intent(in)    ::  iimax
  
  real   , intent(in)    ::  ii0L
  real   , intent(in)    ::  ii0U
  
  real   , intent(in)    ::  ii
  real   , intent(out)   ::  xx
  real   , intent(out)   ::  dx
  
  real                   ::  yy, cmin, cmax
  
  
  if (ii0U == ii0L) then
     ! equidistant grid:
     xx = (ii-1.)*Lmax/(iimax-1.)
     dx =         Lmax/(iimax-1.)
  else
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
  end if
  
  
  end subroutine coord_tanh
  
  
  
  
  
  
  
  
  
  
  

  function atanh(x) result(fn_val)
  ! (sample function)
  
  implicit none
  
  real   , intent(in)  :: x
  real                 :: fn_val
  
  
  if (x == 0.) then
     fn_val = 0.
     return
  end if
  
  fn_val = 0.5* LOG((1.+x)/(1.-x))
  return
  
  
  end function atanh
  
  
  
  
  
  
  
  
  
  
  
  !> \brief evaluation of the real error function
  function erf(x) result(fn_val)
  ! (sample function)
  
  !-----------------------------------------------------------------------
  !             EVALUATION OF THE REAL ERROR FUNCTION
  ! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
  ! Mathematics Library (1993 version).
  ! Adapted by Alan.Miller @ vic.cmis.csiro.au
  !-----------------------------------------------------------------------
  
  implicit none
  
  real   , intent(in)  ::  x
  real                 ::  fn_val
  
  real   , parameter   ::  c = .564189583547756, one = 1.0, half = 0.5, zero = 0.0
  real   , parameter   ::  &
           a(5) = (/  7.71058495001320D-05, -1.33733772997339D-03, 3.23076579225834D-02, 4.79137145607681D-02, 1.28379167095513D-01 /),  &
           b(3) = (/  3.01048631703895D-03,  5.38971687740286D-02, 3.75795757275549D-01 /),  &
           p(8) = (/ -1.36864857382717D-07,  5.64195517478974D-01, 7.21175825088309D+00, 4.31622272220567D+01, 1.52989285046940D+02, 3.39320816734344D+02, 4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00,  1.27827273196294D+01, 7.70001529352295D+01, 2.77585444743988D+02, 6.38980264465631D+02, 9.31354094850610D+02, 7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00,  2.62370141675169D+01, 2.13688200555087D+01, 4.65807828718470D+00, 2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01,  1.87114811799590D+02, 9.90191814623914D+01, 1.80124575948747D+01 /)
  real                 ::  ax, bot, t, top, x2
  
  
  ax = ABS(x)
  
  if (ax <= half) then
     t = x*x
     top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
     bot = ((b(1)*t + b(2))*t + b(3))*t + one
     fn_val = x*(top/bot)
     return
  end if
  
  if (ax <= 4.0) then
     top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
           + p(6))*ax + p(7))*ax + p(8)
     bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
           + q(6))*ax + q(7))*ax + q(8)
     fn_val = half + (half - EXP(-x*x)*top/bot)
     if (x < zero) fn_val = -fn_val
     return
  end if
  
  if (ax < 5.8) then
     x2 = x*x
     t = one / x2
     top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
     bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
     fn_val = (c - top/(x2*bot)) / ax
     fn_val = half + (half - EXP(-x2)*fn_val)
     if (x < zero) fn_val = -fn_val
     return
  end if
  
  fn_val = SIGN(one, x)
  return
  
  
  end function erf
  
  
  
  
  
  
  
  
  
  
  subroutine init_particles(n_part_all,particles)
  ! (sample function)
  
  implicit none
  
  integer, intent(in   ) ::  n_part_all
  real   , intent(out  ) ::  particles(1:n_args,1:n_part_max)
  real                   ::  xx1, xx2, xx3
  real                   ::  LL1, LL2, LL3
  real                   ::  mult
  real                   ::  maxfrac
  integer                ::  n_part_global
  integer                ::  too_many, too_many_global
  integer                ::  i, j, m
  
  
  !--- ratio of maximum initial number and absolute maximum number of particles per subdomain ----------------
  maxfrac = 0.75
  
  
  !===========================================================================================================
  if (L1c .gt. L1) then
     if (rank == 0) write(*,*) 'ERROR! L1c > L1!'
     call MPI_FINALIZE(merror)
     stop
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (L2c .gt. L2) then
     if (rank == 0) write(*,*) 'ERROR! L2c > L2!'
     call MPI_FINALIZE(merror)
     stop
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  if (L3c .gt. L3) then
     if (rank == 0) write(*,*) 'ERROR! L3c > L3!'
     call MPI_FINALIZE(merror)
     stop
  end if
  end if
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
  
  if (dimens == 3) then
     if (x1p(S1p) .lt. L1c-y1_origin .and.  &
      &  x2p(S2p) .lt. L2c-y2_origin .and.  &
      &  x3p(S3p) .lt. L3c-y3_origin) then
        
        if ((BC_1U == 0 .or. BC_1L == -1) .and. x1p(N1p+1) .lt. L1c-y1_origin) then
           LL1 = x1p(N1p+1)   -x1p(S1p)
        else
           LL1 = L1c-y1_origin-x1p(S1p)
        end if
        
        if ((BC_2U == 0 .or. BC_2L == -1) .and. x2p(N2p+1) .lt. L2c-y2_origin) then
           LL2 = x2p(N2p+1)   -x2p(S2p)
        else
           LL2 = L2c-y2_origin-x2p(S2p)
        end if
        
        if ((BC_3U == 0 .or. BC_3L == -1) .and. x3p(N3p+1) .lt. L3c-y3_origin) then
           LL3 = x3p(N3p+1)   -x3p(S3p)
        else
           LL3 = L3c-y3_origin-x3p(S3p)
        end if
     end if
     n_part = NINT(LL1*LL2*LL3/(L1c*L2c*L3c)*REAL(n_part_all))
  else
     if (x1p(S1p) .lt. L1c-y1_origin .and.  &
      &  x2p(S2p) .lt. L2c-y2_origin) then
        
        if ((BC_1U == 0 .or. BC_1L == -1) .and. x1p(N1p+1) .lt. L1c-y1_origin) then
           LL1 = x1p(N1p+1)   -x1p(S1p)
        else
           LL1 = L1c-y1_origin-x1p(S1p)
        end if
        
        if ((BC_2U == 0 .or. BC_2L == -1) .and. x2p(N2p+1) .lt. L2c-y2_origin) then
           LL2 = x2p(N2p+1)   -x2p(S2p)
        else
           LL2 = L2c-y2_origin-x2p(S2p)
        end if
     end if
     n_part = NINT(LL1*LL2    /(L1c*L2c    )*REAL(n_part_all))
  end if
  
  too_many = 0
  if (REAL(n_part) .gt. maxfrac*REAL(n_part_max)) too_many = 1
  
  call MPI_ALLREDUCE(too_many,too_many_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (too_many_global .gt. 0) then
     if (rank == 0) write(*,*) 'ERROR! Too many particles on at least one processor!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  
  
  if (n_part .gt. 0) then
     
     call RANDOM_NUMBER(particles(1:n_args,1:n_part))
     
     do m = 1, n_part
        
        !--- particle position -------------------------------------------------------------------------------
                         particles(1,m) = x1p(S1p) + LL1*particles(1,m)
                         particles(2,m) = x2p(S2p) + LL2*particles(2,m)
        if (dimens == 3) particles(3,m) = x3p(S3p) + LL3*particles(3,m)
        
        !--- particle speed ----------------------------------------------------------------------------------
                         particles(4,m) = 0.
                         particles(5,m) = 0.
        if (dimens == 3) particles(6,m) = 0.
        
        do i = 10, n_args
           !--- Richardson number ----------------------------------------------------------------------------
           if (i == 10) particles(i,m) = Rip(1)
           
           !--- settling velocity ----------------------------------------------------------------------------
           if (i == 11) particles(i,m) = usp(1)
           
           !--- Stokes number --------------------------------------------------------------------------------
           if (i == 12) particles(i,m) = Stp(1)
           
           !--- particle name --------------------------------------------------------------------------------
           if (i == 13) particles(i,m) = 0. ! REAL(m+n_part_max-n_part_all) ! TEST!!! needs to be adapted ...
           
           !--- Gruppen werden nach Position von Argument 1 im Intervall [0:1] eingeteilt --------------------
           if (i == 14) then
              do j = 1, n_groups
                 if (particles(1,m) .ge. REAL(j-1)/REAL(n_groups)*L1 .and. particles(1,m) .le. REAL(j)/REAL(n_groups)*L1) then
                    particles(i,m) = REAL(j)
                 end if
              end do
           end if
        end do
        
     end do
     
  end if
  
  
  !===========================================================================================================
  !=== adjust particle weight ================================================================================
  !===========================================================================================================
  call MPI_ALLREDUCE(n_part,n_part_global,1,MPI_INTEGER,MPI_SUM,COMM_CART,merror)
  
  if (n_part /= 0) write(*,'(a,i10,a,3i4)') 'n_part        = ', n_part, ', iB(1:3)', iB(1:3,1)
  if (rank   == 0) write(*,'(a,i10      )') 'n_part_global = ', n_part_global
  
  if (n_part .gt. 0) then
    if (dimens == 3) then
       mult = L1c*L2c*L3c / REAL(n_part_global)
    else
       mult = L1c*L2c     / REAL(n_part_global)
    end if
    
    i = 10
    do m = 1, n_part
       particles(i,m) = particles(i,m) * mult
    end do
  end if
  !===========================================================================================================
  
  
  end subroutine init_particles
  
end module usr_func
