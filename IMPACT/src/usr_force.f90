!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  !> \brief allows user to define a custom force to modify the right hand side
  !!
  !!--- additional volume forces in the momentum equation ---
  !! note: - du/dt = RHS-nl
  !!       - cf. sketch in file "usr_geometry.f90"
  !!
  !!       grid points in the domain
  !!       |       |       |    velocity component
  !!       |       |       |    |
  !! nl(S11:N11,S21:N21,S31:N31,1)
  !! nl(S12:N12,S22:N22,S32:N32,2)
  !! nl(S13:N13,S23:N23,S33:N33,3)
  !!
  SUBROUTINE forcing_vel
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  USE usr_imbound
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, m
  REAL                   ::  mult
  
  mult = dtime*aRK(substep)
  DO k = S31, N31
     DO j = S21, N21
        DO i = S11, N11
            DO m = 1, dimens
!           nl(i,j,k,1) = nl(i,j,k,1) - 1.*(0. - vel(i,j,k,1))*interface((x2p(j)-0.8*L2)/(0.5*L1)) ! some kind of fringe?
!           nl(i,j,k,1) = (1-distance2ib( x1u(i),x2p(j),x3p(k) ))*nl(i,j,k,1) + (rhs(i,j,k,1)/mult + vel(i,j,k,1))*distance2ib( x1u(i),x2p(j),x3p(k))
!           nl(i,j,k,2) = (1-distance2ib( x1p(i),x2v(j),x3p(k) ))*nl(i,j,k,2) + (rhs(i,j,k,2)/mult + vel(i,j,k,2))*distance2ib( x1p(i),x2v(j),x3p(k))
!           nl(i,j,k,3) = (1-distance2ib( x1p(i),x2p(j),x3w(k) ))*nl(i,j,k,3) + (rhs(i,j,k,3)/mult + vel(i,j,k,3))*distance2ib( x1p(i),x2p(j),x3w(k))
!           nl(i,j,k,1) = nl(i,j,k,1) - 1.e0*distance2ib( x1u(i),x2p(j),x3p(k) )*vel(i,j,k,1)
!           nl(i,j,k,2) = nl(i,j,k,2) - 1.e0*distance2ib( x1p(i),x2v(j),x3p(k) )*vel(i,j,k,2)
!           nl(i,j,k,3) = nl(i,j,k,3) - 1.e0*distance2ib( x1p(i),x2p(j),x3w(k) )*vel(i,j,k,3)
!               WRITE(*,'(a)')  distance2ib( x1u(i),x2p(jWRITE(*,'(a)') ,x3p(k)
               !IF( m==1 ) rhs(i,j,k,1) = (1-distance2ib( x1u(i),x2p(j),x3p(k) ))*rhs(i,j,k,1) + ( mult*nl(i,j,k,1) - vel(i,j,k,1) + vel_ib(x1u(i),x2p(j),x3p(k),1) )*distance2ib( x1u(i),x2p(j),x3p(k))
               !IF( m==2 ) rhs(i,j,k,2) = (1-distance2ib( x1p(i),x2v(j),x3p(k) ))*rhs(i,j,k,2) + ( mult*nl(i,j,k,2) - vel(i,j,k,2) + vel_ib(x1p(i),x2v(j),x3p(k),2) )*distance2ib( x1p(i),x2v(j),x3p(k))
               !IF( m==3 ) rhs(i,j,k,3) = (1-distance2ib( x1p(i),x2p(j),x3w(k) ))*rhs(i,j,k,3) + ( mult*nl(i,j,k,3) - vel(i,j,k,3) + vel_ib(x1p(i),x2p(j),x3w(k),3) )*distance2ib( x1p(i),x2p(j),x3w(k))
               IF( m==1 ) nl(i,j,k,1) =  nl(i,j,k,1) + 2.*exp( -( (x1u(i)-1.)/0.2 )**2 -( ( x2p(j)-1 )/0.2 )**2 -( ( x3p(k)-2. )/0.2 )**2 )
 
               IF( m==2 ) nl(i,j,k,2) =  nl(i,j,k,2) - 1.*exp( -( (x1p(i)-1.)/0.2 )**2 -( ( x2v(j)-1 )/0.2 )**2 -( ( x3p(k)-2. )/0.2 )**2 )*cos( 2*pi*freq*time )

               !IF( m==3 ) nl(i,j,k,3) =  nl(i,j,k,3) -  
              
           END DO
        END DO
     END DO
  END DO
  
  
  END SUBROUTINE forcing_vel
  
  
  
  
  !> \brief allows user to define a custom force to modify the right hand side
  SUBROUTINE forcing_conc
  ! (basic subroutine)
  
  USE mod_dims
  USE mod_vars
  USE usr_vars
  USE usr_func
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  
  
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
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           nlco(i,j,k,conc_nu) = nlco(i,j,k,conc_nu) - 1.*(1. - conc(i,j,k,conc_nu))*interface((x2p(j)-0.8*L2)/(0.5*L1))
        END DO
     END DO
  END DO
  
  
  END SUBROUTINE forcing_conc
  
