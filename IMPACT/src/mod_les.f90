!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

MODULE mod_les
  
  
  USE mod_dims
  USE mod_vars
  USE mod_diff
  
  
  PRIVATE
  
  
  PUBLIC admrt, smag_vel, smag_conc
  PUBLIC filter3D
  
  
  INCLUDE 'mpif.h'
  
  INTEGER                ::  status(MPI_STATUS_SIZE)
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE admrt(dir,m,SS1,SS2,SS3,NN1,NN2,NN3,gn_max,hn_max,chi,phi,nl)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  dir
  INTEGER, INTENT(IN   ) ::  m
  
  INTEGER, INTENT(IN   ) ::  SS1
  INTEGER, INTENT(IN   ) ::  SS2
  INTEGER, INTENT(IN   ) ::  SS3
  
  INTEGER, INTENT(IN   ) ::  NN1
  INTEGER, INTENT(IN   ) ::  NN2
  INTEGER, INTENT(IN   ) ::  NN3
  
  INTEGER, INTENT(IN   ) ::  gn_max
  INTEGER, INTENT(IN   ) ::  hn_max
  
  REAL   , INTENT(IN   ) ::  chi
  REAL   , INTENT(IN   ) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  nl  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  gn, hn
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - "filter3D" kann nicht verwendet werden, weil "Ap" dann Ein- und Ausgabe-Argument waere    !
  !                und damit innerhalb von "filter3D" zu Problem fuehren wuerde.                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  Ap(SS1:NN1,SS2:NN2,SS3:NN3) = phi(SS1:NN1,SS2:NN2,SS3:NN3)
  
  DO hn = 1, hn_max
     
     pp(SS1:NN1,SS2:NN2,SS3:NN3) = Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     
     !--- Tiefpass -------------------------------------------------------------------------------------------
     DO gn = 1, gn_max
        IF (dimens == 3) THEN
           CALL filter(1,dir,m,Ap,Ar)
           CALL filter(2,dir,m,Ar,rr)
           CALL filter(3,dir,m,rr,Ap)
        ELSE
           CALL filter(1,dir,m,Ap,rr)
           CALL filter(2,dir,m,rr,Ap)
        END IF
     END DO
     
     !--- Hochpass -------------------------------------------------------------------------------------------
     Ap(SS1:NN1,SS2:NN2,SS3:NN3) = pp(SS1:NN1,SS2:NN2,SS3:NN3) - Ap(SS1:NN1,SS2:NN2,SS3:NN3)
     
  END DO
  
  nl(SS1:NN1,SS2:NN2,SS3:NN3) = nl(SS1:NN1,SS2:NN2,SS3:NN3) + chi*Ap(SS1:NN1,SS2:NN2,SS3:NN3)
  
  
  END SUBROUTINE admrt
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE smag_vel(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  INTEGER                ::  i, j, k
  
  REAL                   ::  dd1
  
  REAL                   ::  SSii(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)) ! TEST!!! Felder muessen wieder weg ...
  REAL                   ::  SS12(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  SS13(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  SS23(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  smag_deriv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:6)
  
  REAL                   ::  maxvalue       (0:3), minvalue       (0:3) ! TEST!!! nur zum testen ...
  REAL                   ::  maxvalue_global(0:3), minvalue_global(0:3)
  REAL                   ::  test(0:N1,0:N2,0:N3,0:3)
  REAL                   ::  pi
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
  !****************************
  IF (1 == 2) THEN
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  vel = 0.
  
  DO k = S31B, N31B
     DO j = S21B, N21B
        DO i = S11B, N11B
           !vel(i,j,k,1) = 0.*x1u(i)/SQRT(2.)                  + 0.*x2p(j)                           + 0.*x3p(k)
           vel(i,j,k,1) = 0.*SIN(pi*x1u(i)/L1)*L1/pi/SQRT(2.) + 0.*COS(pi*x2p(j)/L2)*L2/pi          + 0.*COS(pi*x3p(k)/L3)*L3/pi
        END DO
     END DO
  END DO
  DO k = S32B, N32B
     DO j = S22B, N22B
        DO i = S12B, N12B
           !vel(i,j,k,2) = 0.*x1p(i)                           + 0.*x2v(j)                           + 0.*x3p(k)
           vel(i,j,k,2) = 0.*COS(pi*x1p(i)/L1)*L1/pi          + 1.*SIN(pi*x2v(j)/L2)*L2/pi/SQRT(2.) + 0.*COS(pi*x3p(k)/L3)*L3/pi
        END DO
     END DO
  END DO
  DO k = S33B, N33B
     DO j = S23B, N23B
        DO i = S13B, N13B
           !vel(i,j,k,3) = 0.*x1p(i)                           + 0.*x2p(j)                           + 0.*x3w(k)/SQRT(2.)
           vel(i,j,k,3) = 0.*COS(pi*x1p(i)/L1)*L1/pi          + 0.*COS(pi*x2p(j)/L2)*L2/pi          + 0.*SIN(pi*x3w(k)/L3)*L3/pi/SQRT(2.)
        END DO
     END DO
  END DO
  
  
  !exch_yes = .TRUE. ==> diese Subroutine aus "timeintegration" aufrufen!!
  END IF
  !****************************
  
  
  ! TEST!!! Vorsicht ist geboten mit dieser Implementation. Habe bis zum Schluss Bugs gefunden ... (09.03.09)
  ! TEST!!! Validieren:
  !            - v.a. die 2D-Variante!
  !            - Compact / Parallel / unterschiedliche Gitterweiten
  !            - mit usr_stats vergleichen!?
  !            - Gradient und Extrapolation auf dem Rand
  ! TEST!!! Abhaengigkeiten der Felder ueberpruefen!! smag_deriv ==> drei Felder koennen eingespart werden!
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Consistency check of the interpolation / differentiation operators at symmetry planes:                   !
  ! --------------------------------------------------------------------------------------                   !
  !                                                                                                          !
  ! if x=x_symm is a symmetry plane, then                                                                    !
  !    u, du/dy, du/dz, dv/dx, dw/dx           are point-symmetric                                           !
  !    v, w, du/dx, dv/dy, dv/dz, dw/dy, dw/dz are axis-symmetric                                            !
  !    u(x_symm) = du/dy(x_symm) = du/dz(x_symm) = dv/dx(x_symm) = dw/dx(x_symm) = 0                         !
  ! if y=y_symm is a symmetry plane, then                                                                    !
  !    v, dv/dx, dv/dz, du/dy, dw/dy           are point-symmetric                                           !
  !    u, w, dv/dy, du/dx, du/dz, dw/dx, dw/dz are axis-symmetric                                            !
  !    v(y_symm) = dv/dx(y_symm) = dv/dz(y_symm) = du/dy(y_symm) = dw/dy(y_symm) = 0                         !
  ! if z=z_symm is a symmetry plane, then                                                                    !
  !    w, dw/dx, dw/dy, du/dz, dv/dz           are point-symmetric                                           !
  !    u, v, dw/dz, du/dx, du/dy, dv/dx, dv/dy are axis-symmetric                                            !
  !    w(z_symm) = dw/dx(z_symm) = dw/dy(z_symm) = du/dz(z_symm) = dv/dz(z_symm) = 0                         !
  !                                                                                                          !
  !                                                                                                          !
  ! hence,                                                                                                   !
  !    if x=x_symm                                                                                           !
  !       S12, S13  are point-symmetric                                                                      !
  !       Sii, S23  are axis-symmetric                                                                       !
  !    if y=y_symm                                                                                           !
  !       S12, S23  are point-symmetric                                                                      !
  !       Sii, S13  are axis-symmetric                                                                       !
  !    if z=z_symm                                                                                           !
  !       S13, S23  are point-symmetric                                                                      !
  !       Sii, S12  are axis-symmetric                                                                       !
  !                                                                                                          !
  ! Sij**2 and therefore |S| are generally axis-symmetric, such that                                         !
  !    if x=x_symm                                                                                           !
  !       |S|*S12, |S|*S13  are point-symmetric                                                              !
  !       |S|*Sii, |S|*S23  are axis-symmetric                                                               !
  !    if y=y_symm                                                                                           !
  !       |S|*S12, |S|*S23  are point-symmetric                                                              !
  !       |S|*Sii, |S|*S13  are axis-symmetric                                                               !
  !    if z=z_symm                                                                                           !
  !       |S|*S13, |S|*S23  are point-symmetric                                                              !
  !       |S|*Sii, |S|*S12  are axis-symmetric                                                               !
  !                                                                                                          !
  ! finally,                                                                                                 !
  !    if x=x_symm                                                                                           !
  !       d(|S|*Sii)/dx, d(|S|*S12)/dy, d(|S|*S13)/dz  are point-symmetric, just as u                        !
  !       all other                                    are axis-symmetric , just as v, w                     !
  !    if y=y_symm                                                                                           !
  !       d(|S|*Sii)/dy, d(|S|*S12)/dx, d(|S|*S23)/dz  are point-symmetric, just as v                        !
  !       all other                                    are axis-symmetric , just as u, w                     !
  !    if z=z_symm                                                                                           !
  !       d(|S|*Sii)/dz, d(|S|*S13)/dx, d(|S|*S23)/dy  are point-symmetric, just as w                        !
  !       all other                                    are axis-symmetric , just as u, v                     !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== high-pass filter / first derivative ===================================================================
  !===========================================================================================================
  IF (n_hp_vel > 0) THEN
                      CALL filter3D(1,0,vel(b1L,b2L,b3L,1),gpre) ! ok!
                      CALL filter3D(2,0,vel(b1L,b2L,b3L,2),com ) ! ok!
     IF (dimens == 3) CALL filter3D(3,0,vel(b1L,b2L,b3L,3),dig ) ! ok!
     
                      gpre(S11B:N11B,S21B:N21B,S31B:N31B) = vel(S11B:N11B,S21B:N21B,S31B:N31B,1) - gpre(S11B:N11B,S21B:N21B,S31B:N31B)
                      com (S12B:N12B,S22B:N22B,S32B:N32B) = vel(S12B:N12B,S22B:N22B,S32B:N32B,2) - com (S12B:N12B,S22B:N22B,S32B:N32B)
     IF (dimens == 3) dig (S13B:N13B,S23B:N23B,S33B:N33B) = vel(S13B:N13B,S23B:N23B,S33B:N33B,3) - dig (S13B:N13B,S23B:N23B,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
                      CALL first_vel_pre(.TRUE.,1,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,gpre,smag_deriv(b1L,b2L,b3L,1)) ! ok!
                      CALL first_vel_pre(.TRUE.,2,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,com ,smag_deriv(b1L,b2L,b3L,2)) ! ok!
     IF (dimens == 3) CALL first_vel_pre(.TRUE.,3,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,dig ,smag_deriv(b1L,b2L,b3L,3)) ! ok!
     
                      CALL first_pre_vel(.TRUE.,2,S11B,S22B,S3p ,N11B,N22B,N3p ,gpre,pp) ! ok!
                      CALL first_pre_vel(.TRUE.,1,S11B,S22B,S3p ,N11B,N22B,N3p ,com ,rr) ! ok!
                      smag_deriv(S11B:N11B,S22B:N22B,S3p:N3p,4) = 0.5*(pp(S11B:N11B ,S22B:N22B,S3p:N3p) + rr(S11B:N11B,S22B:N22B,S3p:N3p))
     
     IF (dimens == 3) CALL first_pre_vel(.TRUE.,3,S11B,S2p ,S33B,N11B,N2p ,N33B,gpre,pp) ! ok!
     IF (dimens == 3) CALL first_pre_vel(.TRUE.,1,S11B,S2p ,S33B,N11B,N2p ,N33B,dig ,rr) ! ok!
     IF (dimens == 3) smag_deriv(S11B:N11B,S2p:N2p,S33B:N33B,5) = 0.5*(pp(S11B:N11B,S2p:N2p,S33B:N33B) + rr(S11B:N11B,S2p:N2p,S33B:N33B))
     
     IF (dimens == 3) CALL first_pre_vel(.TRUE.,3,S1p ,S22B,S33B,N1p ,N22B,N33B,com,pp) ! ok!
     IF (dimens == 3) CALL first_pre_vel(.TRUE.,2,S1p ,S22B,S33B,N1p ,N22B,N33B,dig,rr) ! ok!
     IF (dimens == 3) smag_deriv(S1p:N1p,S22B:N22B,S33B:N33B,6) = 0.5*(pp(S1p:N1p,S22B:N22B,S33B:N33B) + rr(S1p:N1p,S22B:N22B,S33B:N33B))
  !===========================================================================================================
  ELSE
                      CALL first_vel_pre(exch_yes,1,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,vel(b1L,b2L,b3L,1),smag_deriv(b1L,b2L,b3L,1)) ! ok!
                      CALL first_vel_pre(exch_yes,2,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,vel(b1L,b2L,b3L,2),smag_deriv(b1L,b2L,b3L,2)) ! ok!
     IF (dimens == 3) CALL first_vel_pre(exch_yes,3,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,vel(b1L,b2L,b3L,3),smag_deriv(b1L,b2L,b3L,3)) ! ok!
     
                      CALL first_pre_vel(exch_yes,2,S11B,S22B,S3p ,N11B,N22B,N3p ,vel(b1L,b2L,b3L,1),pp) ! ok!
                      CALL first_pre_vel(exch_yes,1,S11B,S22B,S3p ,N11B,N22B,N3p ,vel(b1L,b2L,b3L,2),rr) ! ok!
                      smag_deriv(S11B:N11B,S22B:N22B,S3p:N3p,4) = 0.5*(pp(S11B:N11B,S22B:N22B,S3p:N3p) + rr(S11B:N11B,S22B:N22B,S3p:N3p))
     
     IF (dimens == 3) CALL first_pre_vel(exch_yes,3,S11B,S2p ,S33B,N11B,N2p ,N33B,vel(b1L,b2L,b3L,1),pp) ! ok!
     IF (dimens == 3) CALL first_pre_vel(exch_yes,1,S11B,S2p ,S33B,N11B,N2p ,N33B,vel(b1L,b2L,b3L,3),rr) ! ok!
     IF (dimens == 3) smag_deriv(S11B:N11B,S2p:N2p,S33B:N33B,5) = 0.5*(pp(S11B:N11B,S2p:N2p,S33B:N33B) + rr(S11B:N11B,S2p:N2p,S33B:N33B))
     
     IF (dimens == 3) CALL first_pre_vel(exch_yes,3,S1p ,S22B,S33B,N1p ,N22B,N33B,vel(b1L,b2L,b3L,2),pp) ! ok!
     IF (dimens == 3) CALL first_pre_vel(exch_yes,2,S1p ,S22B,S33B,N1p ,N22B,N33B,vel(b1L,b2L,b3L,3),rr) ! ok!
     IF (dimens == 3) smag_deriv(S1p:N1p,S22B:N22B,S33B:N33B,6) = 0.5*(pp(S1p:N1p,S22B:N22B,S33B:N33B) + rr(S1p:N1p,S22B:N22B,S33B:N33B))
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- 0-0 ---------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           SSii(i,j,k) = 2.*(smag_deriv(i,j,k,1)**2 + smag_deriv(i,j,k,2)**2 + smag_deriv(i,j,k,3)**2)
        END DO
     END DO
  END DO
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  !--- 1-2 ---------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  DO k = S3p, N3p
     DO j = S22B, N22B
        DO i = S11B, N11B
           SS12(i,j,k) = 2.*smag_deriv(i,j,k,4)
        END DO
     END DO
  END DO
  !-----------------------------------------------------------------------------------------------------------
  
  IF (dimens == 3) THEN
  !-----------------------------------------------------------------------------------------------------------
  !--- 1-3 ---------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  DO k = S33B, N33B
     DO j = S2p, N2p
        DO i = S11B, N11B
           SS13(i,j,k) = 2.*smag_deriv(i,j,k,5)
        END DO
     END DO
  END DO
  !-----------------------------------------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------------------------------------
  !--- 2-3 ---------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  DO k = S33B, N33B
     DO j = S22B, N22B
        DO i = S1p, N1p
           SS23(i,j,k) = 2.*smag_deriv(i,j,k,6)
        END DO
     END DO
  END DO
  !-----------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  
  !===========================================================================================================
  !===========================================================================================================
  !===========================================================================================================
  IF (dimens == 3) THEN
     !========================================================================================================
     !--------------------------------------------------------------------------------------------------------
     !--- 0-0 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_vel_pre(.TRUE.,2,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SS12,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,Ap) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SS13,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,rr) ! ok!
     
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SS12,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,rh) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SS23,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,Ar) ! ok!
     
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SS13,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,z1) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SS23,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,z2) ! ok!
     
     IF (concentration_yes) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = (dx1pS(i)*dx2pS(j)*dx3pS(k))**(2./3.) * SQRT(ABS(SSii(i,j,k) + 0.5*(Ap(i,j,k)**2 + rr(i,j,k)**2 + rh(i,j,k)**2 + Ar(i,j,k)**2 + z1(i,j,k)**2 + z2(i,j,k)**2)))
                 smag_deriv(i,j,k,1:dimens) = dd1*smag_deriv(i,j,k,1:dimens)
                 gpre(i,j,k) = dd1
              END DO
           END DO
        END DO
        ! Anmerkung: Es wird vorzugsweise nu_t selbst interpoliert, anstatt vel nochmals colocated abzuleiten, da nu_t evtl. glatter im Raum ist ...
        CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,gpre,smag_diffus(b1L,b2L,b3L,1)) ! ok!
        CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,gpre,smag_diffus(b1L,b2L,b3L,2)) ! ok!
        CALL interpolate_pre_vel(.TRUE.,3,S1p ,S2p ,S33B,N1p ,N2p ,N33B,gpre,smag_diffus(b1L,b2L,b3L,3)) ! ok!
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = (dx1pS(i)*dx2pS(j)*dx3pS(k))**(2./3.) * SQRT(ABS(SSii(i,j,k) + 0.5*(Ap(i,j,k)**2 + rr(i,j,k)**2 + rh(i,j,k)**2 + Ar(i,j,k)**2 + z1(i,j,k)**2 + z2(i,j,k)**2)))
                 smag_deriv(i,j,k,1:dimens) = dd1*smag_deriv(i,j,k,1:dimens)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- 1-2 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,2,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,Ap) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SS13,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,2,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,rr) ! ok!
     
     CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,rh) ! ok!
     CALL interpolate_vel_pre(.TRUE.,3,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SS23,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,Ar) ! ok!
     
     DO k = S3p, N3p
        DO j = S22B, N22B
           DO i = S11B, N11B
              dd1 = (dx1uS(i)*dx2vS(j)*dx3pS(k))**(2./3.) * SQRT(ABS(SS12(i,j,k)**2 + 0.5*(Ap(i,j,k) + rh(i,j,k)) + rr(i,j,k)**2 + Ar(i,j,k)**2))
              smag_deriv(i,j,k,4) = dd1*smag_deriv(i,j,k,4)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- 1-3 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,3,S11B,S2p ,S33B,N11B,N2p ,N33B,pp  ,Ap) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SS12,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,3,S11B,S2p ,S33B,N11B,N2p ,N33B,pp  ,rr) ! ok!
     
     CALL interpolate_pre_vel(.TRUE.,3,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S33B,N11B,N2p ,N33B,pp  ,rh) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SS23,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S33B,N11B,N2p ,N33B,pp  ,Ar) ! ok!
     
     DO k = S33B, N33B
        DO j = S2p, N2p
           DO i = S11B, N11B
              dd1 = (dx1uS(i)*dx2pS(j)*dx3wS(k))**(2./3.) * SQRT(ABS(SS13(i,j,k)**2 + 0.5*(Ap(i,j,k) + rh(i,j,k)) + rr(i,j,k)**2 + Ar(i,j,k)**2))
              smag_deriv(i,j,k,5) = dd1*smag_deriv(i,j,k,5)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- 2-3 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,3,S1p ,S22B,S33B,N1p ,N22B,N33B,pp  ,Ap) ! ok!
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SS12,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,3,S1p ,S22B,S33B,N1p ,N22B,N33B,pp  ,rr) ! ok!
     
     CALL interpolate_pre_vel(.TRUE.,3,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S33B,N1p ,N22B,N33B,pp  ,rh) ! ok!
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S2p ,S33B,N1p ,N2p ,N33B,SS13,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S33B,N1p ,N22B,N33B,pp  ,Ar) ! ok!
     
     DO k = S33B, N33B
        DO j = S22B, N22B
           DO i = S1p, N1p
              dd1 = (dx1pS(i)*dx2vS(j)*dx3wS(k))**(2./3.) * SQRT(ABS(SS23(i,j,k)**2 + 0.5*(Ap(i,j,k) + rh(i,j,k)) + rr(i,j,k)**2 + Ar(i,j,k)**2))
              smag_deriv(i,j,k,6) = dd1*smag_deriv(i,j,k,6)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     !========================================================================================================
     
     
     ! Achtung: Die Normalen Intervallgrenzen werden in den Routinen anders gewaehlt!
     !========================================================================================================
     CALL first_pre_vel(.TRUE.,1,S11,S21,S31,N11,N21,N31,smag_deriv(b1L,b2L,b3L,1),pp) ! ok!
     CALL first_vel_pre(.TRUE.,2,S11,S21,S31,N11,N21,N31,smag_deriv(b1L,b2L,b3L,4),Ap) ! ok!
     CALL first_vel_pre(.TRUE.,3,S11,S21,S31,N11,N21,N31,smag_deriv(b1L,b2L,b3L,5),rr) ! ok!
     
     nl(S11:N11,S21:N21,S31:N31,1) = nl(S11:N11,S21:N21,S31:N31,1) - 2.*chi_vel**2*(pp(S11:N11,S21:N21,S31:N31) + Ap(S11:N11,S21:N21,S31:N31) + rr(S11:N11,S21:N21,S31:N31))
     !--------------------------------------------------------------------------------------------------------
     CALL first_vel_pre(.TRUE.,1,S12,S22,S32,N12,N22,N32,smag_deriv(b1L,b2L,b3L,4),pp) ! ok!
     CALL first_pre_vel(.TRUE.,2,S12,S22,S32,N12,N22,N32,smag_deriv(b1L,b2L,b3L,2),Ap) ! ok!
     CALL first_vel_pre(.TRUE.,3,S12,S22,S32,N12,N22,N32,smag_deriv(b1L,b2L,b3L,6),rr) ! ok!
     
     nl(S12:N12,S22:N22,S32:N32,2) = nl(S12:N12,S22:N22,S32:N32,2) - 2.*chi_vel**2*(pp(S12:N12,S22:N22,S32:N32) + Ap(S12:N12,S22:N22,S32:N32) + rr(S12:N12,S22:N22,S32:N32))
     !--------------------------------------------------------------------------------------------------------
     CALL first_vel_pre(.TRUE.,1,S13,S23,S33,N13,N23,N33,smag_deriv(b1L,b2L,b3L,5),pp) ! ok!
     CALL first_vel_pre(.TRUE.,2,S13,S23,S33,N13,N23,N33,smag_deriv(b1L,b2L,b3L,6),Ap) ! ok!
     CALL first_pre_vel(.TRUE.,3,S13,S23,S33,N13,N23,N33,smag_deriv(b1L,b2L,b3L,3),rr) ! ok!
     
     nl(S13:N13,S23:N23,S33:N33,3) = nl(S13:N13,S23:N23,S33:N33,3) - 2.*chi_vel**2*(pp(S13:N13,S23:N23,S33:N33) + Ap(S13:N13,S23:N23,S33:N33) + rr(S13:N13,S23:N23,S33:N33))
     !========================================================================================================
     
  !===========================================================================================================
  !===========================================================================================================
  !===========================================================================================================
  ELSE
     !========================================================================================================
     !--------------------------------------------------------------------------------------------------------
     !--- 0-0 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_vel_pre(.TRUE.,2,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SS12,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,Ap) ! ok!
     
     CALL interpolate_vel_pre(.TRUE.,1,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SS12,pp) ! ok!
     CALL interpolate_vel_pre(.TRUE.,2,S1p ,S2p ,S3p ,N1p ,N2p ,N3p ,pp  ,rh) ! ok!
     
     IF (concentration_yes) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = dx1pS(i)*dx2pS(j) * SQRT(ABS(SSii(i,j,k) + 0.5*(Ap(i,j,k)**2 + rh(i,j,k)**2)))
                 smag_deriv(i,j,k,1:dimens) = dd1*smag_deriv(i,j,k,1:dimens)
                 gpre(i,j,k) = dd1
              END DO
           END DO
        END DO
        ! Anmerkung: Es wird vorzugsweise nu_t selbst interpoliert, anstatt vel nochmals colocated abzuleiten, da nu_t evtl. glatter im Raum ist ...
        CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,gpre,smag_diffus(b1L,b2L,b3L,1)) ! ok!
        CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,gpre,smag_diffus(b1L,b2L,b3L,2)) ! ok!
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = dx1pS(i)*dx2pS(j) * SQRT(ABS(SSii(i,j,k) + 0.5*(Ap(i,j,k)**2 + rh(i,j,k)**2)))
                 smag_deriv(i,j,k,1:dimens) = dd1*smag_deriv(i,j,k,1:dimens)
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--------------------------------------------------------------------------------------------------------
     !--- 1-2 ------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,2,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,Ap) ! ok!
     
     CALL interpolate_pre_vel(.TRUE.,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,SSii,pp) ! ok!
     CALL interpolate_pre_vel(.TRUE.,1,S11B,S22B,S3p ,N11B,N22B,N3p ,pp  ,rh) ! ok!
     
     DO k = S3p, N3p
        DO j = S22B, N22B
           DO i = S11B, N11B
              dd1 = dx1uS(i)*dx2vS(j) * SQRT(ABS(SS12(i,j,k)**2 + 0.5*(Ap(i,j,k) + rh(i,j,k))))
              smag_deriv(i,j,k,4) = dd1*smag_deriv(i,j,k,4)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     !========================================================================================================
     
     
     ! Achtung: Die Normalen Intervallgrenzen werden in den Routinen anders gewaehlt!
     !========================================================================================================
     CALL first_pre_vel(.TRUE.,1,S11,S21,S31,N11,N21,N31,smag_deriv(b1L,b2L,b3L,1),pp) ! ok!
     CALL first_vel_pre(.TRUE.,2,S11,S21,S31,N11,N21,N31,smag_deriv(b1L,b2L,b3L,4),Ap) ! ok!
     
     nl(S11:N11,S21:N21,S31:N31,1) = nl(S11:N11,S21:N21,S31:N31,1) - 2.*chi_vel**2*(pp(S11:N11,S21:N21,S31:N31) + Ap(S11:N11,S21:N21,S31:N31))
     !--------------------------------------------------------------------------------------------------------
     CALL first_vel_pre(.TRUE.,1,S12,S22,S32,N12,N22,N32,smag_deriv(b1L,b2L,b3L,4),pp) ! ok!
     CALL first_pre_vel(.TRUE.,2,S12,S22,S32,N12,N22,N32,smag_deriv(b1L,b2L,b3L,2),Ap) ! ok!
     
     nl(S12:N12,S22:N22,S32:N32,2) = nl(S12:N12,S22:N22,S32:N32,2) - 2.*chi_vel**2*(pp(S12:N12,S22:N22,S32:N32) + Ap(S12:N12,S22:N22,S32:N32))
     !========================================================================================================
     
  END IF
  !===========================================================================================================
  !===========================================================================================================
  !===========================================================================================================
  
  
  
 ! DO i = S11B, N11B
 !    IF (rank == 0) WRITE(*,*) test(i,N2/2,N3/2,1), COS(x1u(i)*pi/L1)
 ! END DO
 ! DO j = S22B, N22B
 !    IF (rank == 0) WRITE(*,*) test(N1/2,j,N3/2,3), SIN(x2v(j)*pi/L2)
 ! END DO
  !DO k = S33B, N33B
  !   IF (rank == 0) WRITE(*,*) test(N1/2,N2/2,k,2), SIN(x3w(k)*pi/L3)
  !END DO
  
  
  
  IF (1 == 2) THEN
  !****************************
  maxvalue = 0.
  minvalue = 10000.
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S1p, N1p
           maxvalue(0) = MAX(maxvalue(0),test(i,j,k,0))
           minvalue(0) = MIN(minvalue(0),test(i,j,k,0))
        END DO
     END DO
  END DO
  !DO k = S3p, N3p
  !   DO j = S22B, N22B
  !      DO i = S11B, N11B
  DO k = S31, N31
     DO j = S21, N21
        DO i = S11, N11
           maxvalue(1) = MAX(maxvalue(1),test(i,j,k,1))
           minvalue(1) = MIN(minvalue(1),test(i,j,k,1))
        END DO
     END DO
  END DO
  !DO k = S33B, N33B
  !   DO j = S2p, N2p
  !      DO i = S11B, N11B
  DO k = S32, N32
     DO j = S22, N22
        DO i = S12, N12
           maxvalue(2) = MAX(maxvalue(2),test(i,j,k,2))
           minvalue(2) = MIN(minvalue(2),test(i,j,k,2))
        END DO
     END DO
  END DO
  !DO k = S33B, N33B
  !   DO j = S22B, N22B
  !      DO i = S1p, N1p
  DO k = S33, N33
     DO j = S23, N23
        DO i = S13, N13
           maxvalue(3) = MAX(maxvalue(3),test(i,j,k,3))
           minvalue(3) = MIN(minvalue(3),test(i,j,k,3))
        END DO
     END DO
  END DO
  
  CALL MPI_REDUCE(maxvalue,maxvalue_global,4,MPI_REAL8,MPI_MAX,0,COMM_CART,merror)
  CALL MPI_REDUCE(minvalue,minvalue_global,4,MPI_REAL8,MPI_MIN,0,COMM_CART,merror)
  
  maxvalue = maxvalue_global
  minvalue = minvalue_global
  
  IF (rank == 0) THEN
     WRITE(*,'(5E12.5)') maxvalue(0), maxvalue(1), maxvalue(2), maxvalue(3)
     WRITE(*,'(5E12.5)') minvalue(0), minvalue(1), minvalue(2), minvalue(3)
     !WRITE(*,'(5E12.5)') MAX(maxvalue(0),maxvalue(1),maxvalue(2),maxvalue(3))
     !WRITE(*,'(5E12.5)') MIN(minvalue(0),minvalue(1),minvalue(2),minvalue(3))
     WRITE(*,'(5E12.5)') MAX(maxvalue(1),maxvalue(2),maxvalue(3))
     WRITE(*,'(5E12.5)') MIN(minvalue(1),minvalue(2),minvalue(3))
  END IF
  CALL MPI_FINALIZE(merror)
  STOP
  !****************************
  END IF
  
  
  END SUBROUTINE smag_vel
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE smag_conc(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER                ::  m, i, j, k
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !! Bei Symmetrie sind "CALL first_pre_vel" und "CALL first_vel_pre" nur korrekt, solange us_vect = 0.! TEST!!! warum denn das???
  !IF ((us_vec(1,m) /= 0. .AND. (BC_1L == -2 .OR. BC_1U == -2)) .OR.    &
  !  & (us_vec(2,m) /= 0. .AND. (BC_2L == -2 .OR. BC_2U == -2)) .OR.    &
  !  & (us_vec(3,m) /= 0. .AND. (BC_3L == -2 .OR. BC_3U == -2))) THEN
  !   IF (rank == 0) WRITE(*,*) 'ERROR! Smagorinsky + symmetry + us_vect /= 0 not yet correctly implemented!'
  !   CALL MPI_FINALIZE(merror)
  !   STOP
  !END IF
  
  
  m = conc_nu
  
  
  !===========================================================================================================
  !=== high-pass filter ======================================================================================
  !===========================================================================================================
  IF (n_hp_conc(m) > 0) THEN
     CALL filter3D(0,m,conc(b1L,b2L,b3L,m),dig) ! ok!
     dig(S1p:N1p,S2p:N2p,S3p:N3p) = conc(S1p:N1p,S2p:N2p,S3p:N3p,m) - dig(S1p:N1p,S2p:N2p,S3p:N3p)
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== 0-1 ===================================================================================================
  !===========================================================================================================
  IF (n_hp_conc(m) > 0) THEN
     CALL first_pre_vel(.TRUE.  ,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,dig                ,pp)
  ELSE
     CALL first_pre_vel(exch_yes,1,S11B,S2p ,S3p ,N11B,N2p ,N3p ,conc(b1L,b2L,b3L,m),pp)
  END IF
  
  DO k = S3p, N3p
     DO j = S2p, N2p
        DO i = S11B, N11B
           pp(i,j,k) = pp(i,j,k) * smag_diffus(i,j,k,1)
        END DO
     END DO
  END DO
  
  CALL first_vel_pre(.TRUE.,1,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),pp,rr)
  
  nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) = nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) - chi_conc(m)**2*rr(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m))
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== 0-2 ===================================================================================================
  !===========================================================================================================
  IF (n_hp_conc(m) > 0) THEN
     CALL first_pre_vel(.TRUE.  ,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,dig                ,pp)
  ELSE
     CALL first_pre_vel(exch_yes,2,S1p ,S22B,S3p ,N1p ,N22B,N3p ,conc(b1L,b2L,b3L,m),pp)
  END IF
  
  DO k = S3p, N3p
     DO j = S22B, N22B
        DO i = S1p, N1p
           pp(i,j,k) = pp(i,j,k) * smag_diffus(i,j,k,2)
        END DO
     END DO
  END DO
  
  CALL first_vel_pre(.TRUE.,2,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),pp,rr)
  
  nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) = nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) - chi_conc(m)**2*rr(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m))
  !===========================================================================================================
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== 0-3 ===================================================================================================
  !===========================================================================================================
  IF (n_hp_conc(m) > 0) THEN
     CALL first_pre_vel(.TRUE.  ,3,S1p ,S2p ,S33B,N1p ,N2p ,N33B,dig                ,pp)
  ELSE
     CALL first_pre_vel(exch_yes,3,S1p ,S2p ,S33B,N1p ,N2p ,N33B,conc(b1L,b2L,b3L,m),pp)
  END IF
  
  DO k = S33B, N33B
     DO j = S2p, N2p
        DO i = S1p, N1p
           pp(i,j,k) = pp(i,j,k) * smag_diffus(i,j,k,3)
        END DO
     END DO
  END DO
  
  CALL first_vel_pre(.TRUE.,3,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),pp,rr)
  
  nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) = nlco(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m),m) - chi_conc(m)**2*rr(S1c(m):N1c(m),S2c(m):N2c(m),S3c(m):N3c(m))
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE smag_conc
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE filter3D(dir,m,phi,fil)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)    ::  dir
  INTEGER, INTENT(IN)    ::  m
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  fil(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  
  IF (dimens == 3) THEN
     CALL filter(1,dir,m,phi,fil)
     CALL filter(2,dir,m,fil,rr )
     CALL filter(3,dir,m,rr ,fil)
  ELSE
     CALL filter(1,dir,m,phi,rr )
     CALL filter(2,dir,m,rr ,fil)
  END IF
  
  
  END SUBROUTINE filter3D
  
  
  
END MODULE mod_les
