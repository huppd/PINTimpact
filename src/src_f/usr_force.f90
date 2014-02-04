!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  !> \brief allows user to define a custom force to modify the right hand side
  subroutine forcing_vel
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  use usr_imbound
  
  implicit none
  
  integer                ::  i, j, k, m
  real                   ::  mult
  
  mult = dtime*aRK(substep)
  !--- additional volume forces in the momentum equation ---
  ! note: - du/dt = RHS-nl
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !       grid points in the domain
  !       |       |       |    velocity component
  !       |       |       |    |
  ! nl(S11:N11,S21:N21,S31:N31,1)
  ! nl(S12:N12,S22:N22,S32:N32,2)
  ! nl(S13:N13,S23:N23,S33:N33,3)
  !
  do k = S31, N31
     do j = S21, N21
        do i = S11, N11
            do m = 1, dimens
!           nl(i,j,k,1) = nl(i,j,k,1) - 1.*(0. - vel(i,j,k,1))*interface((x2p(j)-0.8*L2)/(0.5*L1)) ! some kind of fringe?
!           nl(i,j,k,1) = (1-distance2ib( x1u(i),x2p(j),x3p(k) ))*nl(i,j,k,1) + (rhs(i,j,k,1)/mult + vel(i,j,k,1))*distance2ib( x1u(i),x2p(j),x3p(k))
!           nl(i,j,k,2) = (1-distance2ib( x1p(i),x2v(j),x3p(k) ))*nl(i,j,k,2) + (rhs(i,j,k,2)/mult + vel(i,j,k,2))*distance2ib( x1p(i),x2v(j),x3p(k))
!           nl(i,j,k,3) = (1-distance2ib( x1p(i),x2p(j),x3w(k) ))*nl(i,j,k,3) + (rhs(i,j,k,3)/mult + vel(i,j,k,3))*distance2ib( x1p(i),x2p(j),x3w(k))
!           nl(i,j,k,1) = nl(i,j,k,1) - 1.e0*distance2ib( x1u(i),x2p(j),x3p(k) )*vel(i,j,k,1)
!           nl(i,j,k,2) = nl(i,j,k,2) - 1.e0*distance2ib( x1p(i),x2v(j),x3p(k) )*vel(i,j,k,2)
!           nl(i,j,k,3) = nl(i,j,k,3) - 1.e0*distance2ib( x1p(i),x2p(j),x3w(k) )*vel(i,j,k,3)
!               WRITE(*,'(a)')  distance2ib( x1u(i),x2p(jWRITE(*,'(a)') ,x3p(k)
               if( m==1 ) rhs(i,j,k,1) = (1-distance2ib( x1u(i),x2p(j),x3p(k) ))*rhs(i,j,k,1) + ( mult*nl(i,j,k,1) - vel(i,j,k,1) + vel_ib(x1u(i),x2p(j),x3p(k),1) )*distance2ib( x1u(i),x2p(j),x3p(k))
               if( m==2 ) rhs(i,j,k,2) = (1-distance2ib( x1p(i),x2v(j),x3p(k) ))*rhs(i,j,k,2) + ( mult*nl(i,j,k,2) - vel(i,j,k,2) + vel_ib(x1p(i),x2v(j),x3p(k),2) )*distance2ib( x1p(i),x2v(j),x3p(k))
               if( m==3 ) rhs(i,j,k,3) = (1-distance2ib( x1p(i),x2p(j),x3w(k) ))*rhs(i,j,k,3) + ( mult*nl(i,j,k,3) - vel(i,j,k,3) + vel_ib(x1p(i),x2p(j),x3w(k),3) )*distance2ib( x1p(i),x2p(j),x3w(k))
!               write(*,*) "vel_ib: ", vel_ib(x1p(i),x2p(j),x3w(k),m)
           end do
        end do
     end do
  end do
  
  
  end subroutine forcing_vel
  
  
  
  
  !> \brief allows user to define a custom force to modify the right hand side
  subroutine forcing_conc
  ! (basic subroutine)
  
  use mod_dims
  use mod_vars
  use usr_vars
  use usr_func
  
  implicit none
  
  integer                ::  i, j, k
  
  
  !--- additional sources in concentration transport equation ---
  ! note: - source terms for concentration index "conc_nu"
  !       - dc/dt = RHS-nl
  !       - cf. sketch in file "usr_geometry.f90"
  !
  !         grid points in the domain and on the boundary
  !         |       |       |     concentration index
  !         |       |       |     |
  ! nlco(S1p:N1p,S2p:N2p,S3p:N3p,1:n_conc)
  !
  do k = S3p, N3p
     do j = S2p, N2p
        do i = S1p, N1p
           nlco(i,j,k,conc_nu) = nlco(i,j,k,conc_nu) - 1.*(1. - conc(i,j,k,conc_nu))*interface((x2p(j)-0.8*L2)/(0.5*L1))
        end do
     end do
  end do
  
  
  end subroutine forcing_conc
  
