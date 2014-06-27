!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing many routines to compute derivatives ...
MODULE mod_diff
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  
  
  PRIVATE
  
  
  PUBLIC apply_compact, apply_compact_transp
  PUBLIC divergence, divergence2, divergence_transp
  PUBLIC gradient, gradient_transp
  PUBLIC Helmholtz, Helmholtz_explicit, Helmholtz_conc, Helmholtz_conc_explicit
  PUBLIC nonlinear, nonlinear_conc
  PUBLIC interpolate_vel, interpolate_conc
  PUBLIC outflow_bc, sediment_bc
  PUBLIC filter
  PUBLIC bc_extrapolation, bc_extrapolation_transp
  PUBLIC interpolate_pre_vel , interpolate_vel_pre
  PUBLIC interpolate2_pre_vel, interpolate2_vel_pre
  PUBLIC first_pre_vel, first_vel_pre
  PUBLIC first_adv_pre, first_adv_vel
  PUBLIC Helmholtz_pre_explicit ! TEST!!!
  
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  ! TEST!!! Generell bei Randbehandlung bei ii=0,1 beginnen, bzw. bis ii=N1 rechnen!  
  
  
  ! TEST!!!
  !       - Multiplikation mit 1/dx in Subroutine bzw. in die Koeffizienten hineinziehen (noch besser)
  !       - in dieser Routine gibt es noch Potential für Vektorisierung!!
  !       - Intervalle doch besser erweitern, d.h. S11 --> S11B. Es wird zwar mehr gerechnet, dafuer sind die Operatoren dann allgemeiner (z.B. fuer LES) ...
  !       - ndR nach Unten/Oben unterscheiden ...
  !       - vel_dir ist momentan ueberfluessig ...
!pgi$r unroll = n:8
  SUBROUTINE apply_compact(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,Nmax,ndL,ndR,dimS,ccL,ccL_LU,ccR,WW,Schur,phi,der)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  dir, vel_dir
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  INTEGER, INTENT(IN   ) ::  Nmax, ndL, ndR, dimS
  REAL   , INTENT(IN   ) ::  ccL (-ndL:ndL,0:Nmax)
  REAL   , INTENT(IN   ) ::  ccR (-ndR:ndR,0:Nmax)
  REAL   , INTENT(IN   ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  REAL   , INTENT(IN   ) ::  WW(1:2*ndL,0:Nmax)
  REAL   , INTENT(IN   ) ::  Schur(1:dimS,1:dimS)
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii, i0
  INTEGER                ::  j, jj, j0
  INTEGER                ::  k, kk, k0
  
  INTEGER                ::  dimG
  INTEGER                ::  KD, KO, LM
  
  
  ! Wurde herausgezogen:
  !!CALL exchange (dir,vel_dir,phi)
  !CALL exchange2(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  CALL pseudocall_int(vel_dir) ! TEST!!! this action has no meaning - it shall only to suppress the warnings
                               ! that the argument vel_dir is not used when building the executable ...
  
  KD = ndL+ndL+1
  KO = ndL+ndL
  
  !===========================================================================================================
  IF (dir == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     ! TEST!!! Block-Raender werden nun IMMER per Schur-Komplement behandelt (Pivoting-Problem) ...
     !         Man koennte/sollte aber immerhin das Umspeichern von buffer1A auf buffer1B vermeiden ...
     IF (1 == 2 .AND. NB1 == 1 .AND. BC_1L_global /= -1) THEN
        
        DO k = SS3, NN3
           DO j = SS2, NN2
              
              !--- expliziter Teil ---
              DO i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
                 DO ii = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
              
              !--- L^-1*der ---
              DO i = SS1, NN1-1
                 LM = MIN(ndL,NN1-i)
                 DO ii = 1, LM
                    der(i+ii,j,k) = der(i+ii,j,k) - ccL_LU(KD+ii,i)*der(i,j,k)
                 END DO
              END DO
              
              !--- U^-1*L^-1*der ---
              DO i = NN1, SS1, -1
                 LM = MAX(SS1,i-KO)
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,i)
                 DO ii = i-1, LM, -1
                    der(ii,j,k) = der(ii,j,k) - ccL_LU(KD+ii-i,i)*der(i,j,k)
                 END DO
              END DO
              
           END DO
        END DO
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NN1-SS1+1-2*ndL
        i0   = 2*ndL*(iB(1,1)-1)
        
        DO k = SS3, NN3
           DO j = SS2, NN2
              
              !--- expliziter Teil ---
              DO i = SS1, NN1 ! Anmerkung: Man koennte hier nur den für buffer1A notwendigen Teil von "der" berechnen und den anderen erst nach MPI_ALLGATHERv ...
                 der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
                 DO ii = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
              
              !--- Umspeichern ---
              DO i = 1, ndL
                 buffer1A(j,k,i    ) = der(SS1+i-1  ,j,k)
                 buffer1A(j,k,i+ndL) = der(NN1+i-ndL,j,k)
              END DO
              
              !--- RHS fuer Schur-Komplement ---
              DO i = 1, dimG
                 DO ii = 1, 2*ndL
                    buffer1A(j,k,ii) = buffer1A(j,k,ii) - WW(ii,SS1+i+ndL-1)*der(SS1+i+ndL-1,j,k)
                 END DO
              END DO
              
           END DO
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        ! TEST!!! Evtl. ist das nicht die effizenteste Methode!!
        CALL MPI_ALLGATHERv(buffer1A,2*ndl*(N2+1)*(N3+1),MPI_REAL8,buffer1B,recv1,disp1,MPI_REAL8,COMM_BAR1,merror)
        
        DO k = SS3, NN3
           DO j = SS2, NN2
              
              !--- Loesung des Schur-Komplement-LGS ---
              DO i = 1, ndL
                 der(SS1+i-1  ,j,k) = Schur(i+i0    ,1)*buffer1B(j,k,1)
                 der(NN1+i-ndL,j,k) = Schur(i+i0+ndL,1)*buffer1B(j,k,1)
                 DO ii = 2, dimS
                    der(SS1+i-1  ,j,k) = der(SS1+i-1  ,j,k) + Schur(i+i0    ,ii)*buffer1B(j,k,ii)
                    der(NN1+i-ndL,j,k) = der(NN1+i-ndL,j,k) + Schur(i+i0+ndL,ii)*buffer1B(j,k,ii)
                 END DO
              END DO
              
              !--- RHS des uebrigen Feldes anpassen ---
              DO i = 1, ndL
                 DO ii = i, ndL
                    der(SS1+i-1+ndL,j,k) = der(SS1+i-1+ndL,j,k) - ccL(-ii,SS1+i-1+ndL)*der(SS1+i-1+ndL-ii,j,k)
                    der(NN1-i+1-ndL,j,k) = der(NN1-i+1-ndL,j,k) - ccL( ii,NN1-i+1-ndL)*der(NN1-i+1-ndL+ii,j,k)
                 END DO
              END DO
              
              !--- L^-1*der ---
              DO i = SS1+ndL, NN1-ndL-1
                 LM = MIN(ndL,NN1-ndL-i)
                 DO ii = 1, LM
                    der(i+ii,j,k) = der(i+ii,j,k) - ccL_LU(KD+ii,i)*der(i,j,k)
                 END DO
              END DO
              
              !--- U^-1*L^-1*der ---
              DO i = NN1-ndL, SS1+ndL, -1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,i)
                 LM = MAX(SS1+ndL,i-KO)
                 DO ii = i-1, LM, -1
                    der(ii,j,k) = der(ii,j,k) - ccL_LU(KD+ii-i,i)*der(i,j,k)
                 END DO
              END DO
              
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (1 == 2 .AND. NB2 == 1 .AND. BC_2L_global /= -1) THEN
        
        DO k = SS3, NN3
           
           !--- expliziter Teil ---
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
                 DO jj = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
           
           !--- L^-1*der ---
           DO j = SS2, NN2-1
              LM = MIN(ndL,NN2-j)
              DO i = SS1, NN1
                 DO jj = 1, LM
                    der(i,j+jj,k) = der(i,j+jj,k) - ccL_LU(KD+jj,j)*der(i,j,k)
                 END DO
              END DO
           END DO
           
           !--- U^-1*L^-1*der ---
           DO j = NN2, SS2, -1
              LM = MAX(SS2,j-KO)
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,j)
                 DO jj = j-1, LM, -1
                    der(i,jj,k) = der(i,jj,k) - ccL_LU(KD+jj-j,j)*der(i,j,k)
                 END DO
              END DO
           END DO
           
        END DO
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NN2-SS2+1-2*ndL
        j0   = 2*ndL*(iB(2,1)-1)
        
        DO k = SS3, NN3
           
           !--- expliziter Teil ---
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
                 DO jj = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
           
           !--- Umspeichern ---
           DO j = 1, ndL
              DO i = SS1, NN1
                 buffer2A(i,k,j    ) = der(i,SS2+j-1  ,k)
                 buffer2A(i,k,j+ndL) = der(i,NN2+j-ndL,k)
              END DO
           END DO
           
           !--- RHS fuer Schur-Komplement ---
           DO j = 1, dimG
              DO i = SS1, NN1
                 DO jj = 1, 2*ndL
                    buffer2A(i,k,jj) = buffer2A(i,k,jj) - WW(jj,SS2+j+ndL-1)*der(i,SS2+j+ndL-1,k)
                 END DO
              END DO
           END DO
           
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        CALL MPI_ALLGATHERv(buffer2A,2*ndl*(N1+1)*(N3+1),MPI_REAL8,buffer2B,recv2,disp2,MPI_REAL8,COMM_BAR2,merror)
        
        DO k = SS3, NN3
           
           !--- Loesung des Schur-Komplement-LGS ---
           DO j = 1, ndL
              DO i = SS1, NN1
                 der(i,SS2+j-1  ,k) = Schur(j+j0    ,1)*buffer2B(i,k,1)
                 der(i,NN2+j-ndL,k) = Schur(j+j0+ndL,1)*buffer2B(i,k,1)
                 DO jj = 2, dimS
                    der(i,SS2+j-1  ,k) = der(i,SS2+j-1  ,k) + Schur(j+j0    ,jj)*buffer2B(i,k,jj)
                    der(i,NN2+j-ndL,k) = der(i,NN2+j-ndL,k) + Schur(j+j0+ndL,jj)*buffer2B(i,k,jj)
                 END DO
              END DO
           END DO
           
           !--- RHS des uebrigen Feldes anpassen ---
           DO j = 1, ndL
              DO i = SS1, NN1
                 DO jj = j, ndL
                    der(i,SS2+j-1+ndL,k) = der(i,SS2+j-1+ndL,k) - ccL(-jj,SS2+j-1+ndL)*der(i,SS2+j-1+ndL-jj,k)
                    der(i,NN2-j+1-ndL,k) = der(i,NN2-j+1-ndL,k) - ccL( jj,NN2-j+1-ndL)*der(i,NN2-j+1-ndL+jj,k)
                 END DO
              END DO
           END DO
           
           !--- L^-1*der ---
           DO j = SS2+ndL, NN2-ndL-1
              LM = MIN(ndL,NN2-ndL-j)
              DO i = SS1, NN1
                 DO jj = 1, LM
                    der(i,j+jj,k) = der(i,j+jj,k) - ccL_LU(KD+jj,j)*der(i,j,k)
                 END DO
              END DO
           END DO
           
           !--- U^-1*L^-1*der ---
           DO j = NN2-ndL, SS2+ndL, -1
              LM = MAX(SS2+ndL,j-KO)
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,j)
                 DO jj = j-1, LM, -1
                    der(i,jj,k) = der(i,jj,k) - ccL_LU(KD+jj-j,j)*der(i,j,k)
                 END DO
              END DO
           END DO
           
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (1 == 2 .AND. NB3 == 1 .AND. BC_3L_global /= -1) THEN
        
        !--- expliziter Teil ---
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
                 DO kk = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- L^-1*der ---
        DO k = SS3, NN3-1
           LM = MIN(ndL,NN3-k)
           DO j = SS2, NN2
              DO i = SS1, NN1
                 DO kk = 1, LM
                    der(i,j,k+kk) = der(i,j,k+kk) - ccL_LU(KD+kk,k)*der(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
        !--- U^-1*L^-1*der ---
        DO k = NN3, SS3, -1
           LM = MAX(SS3,k-KO)
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,k)
                 DO kk = k-1, LM, -1
                    der(i,j,kk) = der(i,j,kk) - ccL_LU(KD+kk-k,k)*der(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NN3-SS3+1-2*ndL
        k0   = 2*ndL*(iB(3,1)-1)
        
        !--- expliziter Teil ---
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
                 DO kk = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- Umspeichern ---
        DO k = 1, ndL
           DO j = SS2, NN2
              DO i = SS1, NN1
                 buffer3A(i,j,k    ) = der(i,j,SS3+k-1  )
                 buffer3A(i,j,k+ndL) = der(i,j,NN3+k-ndL)
              END DO
           END DO
        END DO
        
        !--- RHS fuer Schur-Komplement ---
        DO k = 1, dimG
           DO j = SS2, NN2
              DO i = SS1, NN1
                 DO kk = 1, 2*ndL
                    buffer3A(i,j,kk) = buffer3A(i,j,kk) - WW(kk,SS3+k+ndL-1)*der(i,j,SS3+k+ndL-1)
                 END DO
              END DO
           END DO
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        CALL MPI_ALLGATHERv(buffer3A,2*ndl*(N1+1)*(N2+1),MPI_REAL8,buffer3B,recv3,disp3,MPI_REAL8,COMM_BAR3,merror)
        
        !--- Loesung des Schur-Komplement-LGS ---
        DO k = 1, ndL
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,SS3+k-1  ) = Schur(k+k0    ,1)*buffer3B(i,j,1)
                 der(i,j,NN3+k-ndL) = Schur(k+k0+ndL,1)*buffer3B(i,j,1)
                 DO kk = 2, dimS
                    der(i,j,SS3+k-1  ) = der(i,j,SS3+k-1  ) + Schur(k+k0    ,kk)*buffer3B(i,j,kk)
                    der(i,j,NN3+k-ndL) = der(i,j,NN3+k-ndL) + Schur(k+k0+ndL,kk)*buffer3B(i,j,kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- RHS des uebrigen Feldes anpassen ---
        DO k = 1, ndL
           DO j = SS2, NN2
              DO i = SS1, NN1
                 DO kk = k, ndL
                    der(i,j,SS3+k-1+ndL) = der(i,j,SS3+k-1+ndL) - ccL(-kk,SS3+k-1+ndL)*der(i,j,SS3+k-1+ndL-kk)
                    der(i,j,NN3-k+1-ndL) = der(i,j,NN3-k+1-ndL) - ccL( kk,NN3-k+1-ndL)*der(i,j,NN3-k+1-ndL+kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- L^-1*der ---
        DO k = SS3+ndL, NN3-ndL-1
           LM = MIN(ndL,NN3-ndL-k)
           DO j = SS2, NN2
              DO i = SS1, NN1
                 DO kk = 1, LM
                    der(i,j,k+kk) = der(i,j,k+kk) - ccL_LU(KD+kk,k)*der(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
        !--- U^-1*L^-1*der ---
        DO k = NN3-ndL, SS3+ndL, -1
           LM = MAX(SS3+ndL,k-KO)
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,k)
                 DO kk = k-1, LM, -1
                    der(i,j,kk) = der(i,j,kk) - ccL_LU(KD+kk-k,k)*der(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE apply_compact
  
  
  
  
  
  
  
  
  
  
  ! ACTHUNG: phi wird überschrieben (da ohnehin mit dx vormultipliziert werden muss)
  SUBROUTINE apply_compact_transp(dir,vel_dir,SI1,SI2,SI3,NI1,NI2,NI3,SE1,SE2,SE3,NE1,NE2,NE3,Nmax,ndL,ndR,dimS,ccL,ccL_LU,ccR,WW,Schur,phi,der)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  dir, vel_dir
  INTEGER, INTENT(IN   ) ::  SI1, SI2, SI3, NI1, NI2, NI3
  INTEGER, INTENT(IN   ) ::  SE1, SE2, SE3, NE1, NE2, NE3
  INTEGER, INTENT(IN   ) ::  Nmax, ndL, ndR, dimS
  REAL   , INTENT(IN   ) ::  ccL (-ndL:ndL,0:Nmax)
  REAL   , INTENT(IN   ) ::  ccR (-ndR:ndR,0:Nmax)
  REAL   , INTENT(IN   ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  REAL   , INTENT(IN   ) ::  WW(1:2*ndL,0:Nmax)
  REAL   , INTENT(IN   ) ::  Schur(1:dimS,1:dimS)
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii, i0
  INTEGER                ::  j, jj, j0
  INTEGER                ::  k, kk, k0
  
  INTEGER                ::  dimG
  INTEGER                ::  KD, KO, LM
  
  
  KD = ndL+ndL+1
  KO = ndL+ndL
  
  !===========================================================================================================
  IF (dir == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (NB1 == 1 .AND. BC_1L_global /= -1) THEN
        
        DO k = SI3, NI3
           DO j = SI2, NI2
              
              !--- L^-1*der ---
              DO i = SI1, NI1-1
                 LM = MIN(ndL,NI1-i)
                 DO ii = 1, LM
                    phi(i+ii,j,k) = phi(i+ii,j,k) - ccL_LU(KD+ii,i)*phi(i,j,k)
                 END DO
              END DO
              
              !--- U^-1*L^-1*der ---
              DO i = NI1, SI1, -1
                 LM = MAX(SI1,i-KO)
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,i)
                 DO ii = i-1, LM, -1
                    phi(ii,j,k) = phi(ii,j,k) - ccL_LU(KD+ii-i,i)*phi(i,j,k)
                 END DO
              END DO
              
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NI1-SI1+1-2*ndL
        i0   = 2*ndL*(iB(1,1)-1)
        
        DO k = SI3, NI3
           DO j = SI2, NI2
              
              !--- Umspeichern ---
              DO i = 1, ndL
                 buffer1A(j,k,i    ) = phi(SI1+i-1  ,j,k)
                 buffer1A(j,k,i+ndL) = phi(NI1+i-ndL,j,k)
              END DO
              
              !--- RHS fuer Schur-Komplement ---
              DO i = 1, dimG
                 DO ii = 1, 2*ndL
                    buffer1A(j,k,ii) = buffer1A(j,k,ii) - WW(ii,SI1+i+ndL-1)*phi(SI1+i+ndL-1,j,k)
                 END DO
              END DO
              
           END DO
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        CALL MPI_ALLGATHERv(buffer1A,2*ndl*(N2+1)*(N3+1),MPI_REAL8,buffer1B,recv1,disp1,MPI_REAL8,COMM_BAR1,merror)
        
        DO k = SI3, NI3
           DO j = SI2, NI2
              
              !--- Loesung des Schur-Komplement-LGS ---
              DO i = 1, ndL
                 phi(SI1+i-1  ,j,k) = Schur(i+i0    ,1)*buffer1B(j,k,1)
                 phi(NI1+i-ndL,j,k) = Schur(i+i0+ndL,1)*buffer1B(j,k,1)
                 DO ii = 2, dimS
                    phi(SI1+i-1  ,j,k) = phi(SI1+i-1  ,j,k) + Schur(i+i0    ,ii)*buffer1B(j,k,ii)
                    phi(NI1+i-ndL,j,k) = phi(NI1+i-ndL,j,k) + Schur(i+i0+ndL,ii)*buffer1B(j,k,ii)
                 END DO
              END DO
              
              !--- RHS des uebrigen Feldes anpassen ---
              DO i = 1, ndL
                 DO ii = i, ndL
                    phi(SI1+i-1+ndL,j,k) = phi(SI1+i-1+ndL,j,k) - ccL(-ii,SI1+i-1+ndL)*phi(SI1+i-1+ndL-ii,j,k)
                    phi(NI1-i+1-ndL,j,k) = phi(NI1-i+1-ndL,j,k) - ccL( ii,NI1-i+1-ndL)*phi(NI1-i+1-ndL+ii,j,k)
                 END DO
              END DO
              
              !--- L^-1*der ---
              DO i = SI1+ndL, NI1-ndL-1
                 LM = MIN(ndL,NI1-ndL-i)
                 DO ii = 1, LM
                    phi(i+ii,j,k) = phi(i+ii,j,k) - ccL_LU(KD+ii,i)*phi(i,j,k)
                 END DO
              END DO
              
              !--- U^-1*L^-1*der ---
              DO i = NI1-ndL, SI1+ndL, -1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,i)
                 LM = MAX(SI1+ndL,i-KO)
                 DO ii = i-1, LM, -1
                    phi(ii,j,k) = phi(ii,j,k) - ccL_LU(KD+ii-i,i)*phi(i,j,k)
                 END DO
              END DO
              
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     DO k = SE3, NE3
        DO j = SE2, NE2
           DO i = SE1, NE1
              der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
              DO ii = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
              END DO
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (NB2 == 1 .AND. BC_2L_global /= -1) THEN
        
        DO k = SI3, NI3
           
           !--- L^-1*der ---
           DO j = SI2, NI2-1
              LM = MIN(ndL,NI2-j)
              DO i = SI1, NI1
                 DO jj = 1, LM
                    phi(i,j+jj,k) = phi(i,j+jj,k) - ccL_LU(KD+jj,j)*phi(i,j,k)
                 END DO
              END DO
           END DO
           
           !--- U^-1*L^-1*der ---
           DO j = NI2, SI2, -1
              LM = MAX(SI2,j-KO)
              DO i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,j)
                 DO jj = j-1, LM, -1
                    phi(i,jj,k) = phi(i,jj,k) - ccL_LU(KD+jj-j,j)*phi(i,j,k)
                 END DO
              END DO
           END DO
           
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NI2-SI2+1-2*ndL
        j0   = 2*ndL*(iB(2,1)-1)
        
        DO k = SI3, NI3
           
           !--- Umspeichern ---
           DO j = 1, ndL
              DO i = SI1, NI1
                 buffer2A(i,k,j    ) = phi(i,SI2+j-1  ,k)
                 buffer2A(i,k,j+ndL) = phi(i,NI2+j-ndL,k)
              END DO
           END DO
           
           !--- RHS fuer Schur-Komplement ---
           DO j = 1, dimG
              DO i = SI1, NI1
                 DO jj = 1, 2*ndL
                    buffer2A(i,k,jj) = buffer2A(i,k,jj) - WW(jj,SI2+j+ndL-1)*phi(i,SI2+j+ndL-1,k)
                 END DO
              END DO
           END DO
           
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        CALL MPI_ALLGATHERv(buffer2A,2*ndl*(N1+1)*(N3+1),MPI_REAL8,buffer2B,recv2,disp2,MPI_REAL8,COMM_BAR2,merror)
        
        DO k = SI3, NI3
           
           !--- Loesung des Schur-Komplement-LGS ---
           DO j = 1, ndL
              DO i = SI1, NI1
                 phi(i,SI2+j-1  ,k) = Schur(j+j0    ,1)*buffer2B(i,k,1)
                 phi(i,NI2+j-ndL,k) = Schur(j+j0+ndL,1)*buffer2B(i,k,1)
                 DO jj = 2, dimS
                    phi(i,SI2+j-1  ,k) = phi(i,SI2+j-1  ,k) + Schur(j+j0    ,jj)*buffer2B(i,k,jj)
                    phi(i,NI2+j-ndL,k) = phi(i,NI2+j-ndL,k) + Schur(j+j0+ndL,jj)*buffer2B(i,k,jj)
                 END DO
              END DO
           END DO
           
           !--- RHS des uebrigen Feldes anpassen ---
           DO j = 1, ndL
              DO i = SI1, NI1
                 DO jj = j, ndL
                    phi(i,SI2+j-1+ndL,k) = phi(i,SI2+j-1+ndL,k) - ccL(-jj,SI2+j-1+ndL)*phi(i,SI2+j-1+ndL-jj,k)
                    phi(i,NI2-j+1-ndL,k) = phi(i,NI2-j+1-ndL,k) - ccL( jj,NI2-j+1-ndL)*phi(i,NI2-j+1-ndL+jj,k)
                 END DO
              END DO
           END DO
           
           !--- L^-1*der ---
           DO j = SI2+ndL, NI2-ndL-1
              LM = MIN(ndL,NI2-ndL-j)
              DO i = SI1, NI1
                 DO jj = 1, LM
                    phi(i,j+jj,k) = phi(i,j+jj,k) - ccL_LU(KD+jj,j)*phi(i,j,k)
                 END DO
              END DO
           END DO
           
           !--- U^-1*L^-1*der ---
           DO j = NI2-ndL, SI2+ndL, -1
              LM = MAX(SI2+ndL,j-KO)
              DO i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,j)
                 DO jj = j-1, LM, -1
                    phi(i,jj,k) = phi(i,jj,k) - ccL_LU(KD+jj-j,j)*phi(i,j,k)
                 END DO
              END DO
           END DO
           
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     DO k = SE3, NE3
        DO j = SE2, NE2
           DO i = SE1, NE1
              der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
              DO jj = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
              END DO
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (NB3 == 1 .AND. BC_3L_global /= -1) THEN
        
        !--- L^-1*der ---
        DO k = SI3, NI3-1
           LM = MIN(ndL,NI3-k)
           DO j = SI2, NI2
              DO i = SI1, NI1
                 DO kk = 1, LM
                    phi(i,j,k+kk) = phi(i,j,k+kk) - ccL_LU(KD+kk,k)*phi(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
        !--- U^-1*L^-1*der ---
        DO k = NI3, SI3, -1
           LM = MAX(SI3,k-KO)
           DO j = SI2, NI2
              DO i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,k)
                 DO kk = k-1, LM, -1
                    phi(i,j,kk) = phi(i,j,kk) - ccL_LU(KD+kk-k,k)*phi(i,j,k)
                 END DO
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        dimG = NI3-SI3+1-2*ndL
        k0   = 2*ndL*(iB(3,1)-1)
        
        !--- Umspeichern ---
        DO k = 1, ndL
           DO j = SI2, NI2
              DO i = SI1, NI1
                 buffer3A(i,j,k    ) = phi(i,j,SI3+k-1  )
                 buffer3A(i,j,k+ndL) = phi(i,j,NI3+k-ndL)
              END DO
           END DO
        END DO
        
        !--- RHS fuer Schur-Komplement ---
        DO k = 1, dimG
           DO j = SI2, NI2
              DO i = SI1, NI1
                 DO kk = 1, 2*ndL
                    buffer3A(i,j,kk) = buffer3A(i,j,kk) - WW(kk,SI3+k+ndL-1)*phi(i,j,SI3+k+ndL-1)
                 END DO
              END DO
           END DO
        END DO
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        CALL MPI_ALLGATHERv(buffer3A,2*ndl*(N1+1)*(N2+1),MPI_REAL8,buffer3B,recv3,disp3,MPI_REAL8,COMM_BAR3,merror)
        
        !--- Loesung des Schur-Komplement-LGS ---
        DO k = 1, ndL
           DO j = SI2, NI2
              DO i = SI1, NI1
                 phi(i,j,SI3+k-1  ) = Schur(k+k0    ,1)*buffer3B(i,j,1)
                 phi(i,j,NI3+k-ndL) = Schur(k+k0+ndL,1)*buffer3B(i,j,1)
                 DO kk = 2, dimS
                    phi(i,j,SI3+k-1  ) = phi(i,j,SI3+k-1  ) + Schur(k+k0    ,kk)*buffer3B(i,j,kk)
                    phi(i,j,NI3+k-ndL) = phi(i,j,NI3+k-ndL) + Schur(k+k0+ndL,kk)*buffer3B(i,j,kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- RHS des uebrigen Feldes anpassen ---
        DO k = 1, ndL
           DO j = SI2, NI2
              DO i = SI1, NI1
                 DO kk = k, ndL
                    phi(i,j,SI3+k-1+ndL) = phi(i,j,SI3+k-1+ndL) - ccL(-kk,SI3+k-1+ndL)*phi(i,j,SI3+k-1+ndL-kk)
                    phi(i,j,NI3-k+1-ndL) = phi(i,j,NI3-k+1-ndL) - ccL( kk,NI3-k+1-ndL)*phi(i,j,NI3-k+1-ndL+kk)
                 END DO
              END DO
           END DO
        END DO
        
        !--- L^-1*der ---
        DO k = SI3+ndL, NI3-ndL-1
           LM = MIN(ndL,NI3-ndL-k)
           DO j = SI2, NI2
              DO i = SI1, NI1
                 DO kk = 1, LM
                    phi(i,j,k+kk) = phi(i,j,k+kk) - ccL_LU(KD+kk,k)*phi(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
        !--- U^-1*L^-1*der ---
        DO k = NI3-ndL, SI3+ndL, -1
           LM = MAX(SI3+ndL,k-KO)
           DO j = SI2, NI2
              DO i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,k)
                 DO kk = k-1, LM, -1
                    phi(i,j,kk) = phi(i,j,kk) - ccL_LU(KD+kk-k,k)*phi(i,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     CALL exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     DO k = SE3, NE3
        DO j = SE2, NE2
           DO i = SE1, NE1
              der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
              DO kk = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
              END DO
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE apply_compact_transp
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE divergence(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)    ::  m
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(m,m,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phi,div)
        
        DO k = S3p, N3p
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 div(i,j,k) = div(i,j,k)*dx1DM(i)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phi,com)
        
        DO k = S3p, N3p
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 div(i,j,k) = div(i,j,k) + com(i,j,k)*dx2DM(j)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phi,com)
        
        DO k = S3p, N3p
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 div(i,j,k) = div(i,j,k) + com(i,j,k)*dx3DM(k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE divergence2(phi,div)
  
  IMPLICIT NONE
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  REAL   , INTENT(  OUT) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  CALL exchange(1,1,phi(b1L,b2L,b3L,1))
  CALL exchange(2,2,phi(b1L,b2L,b3L,2))
  CALL exchange(3,3,phi(b1L,b2L,b3L,3))
  
  
  !===========================================================================================================
  IF (comp_div_yes) THEN
     !--------------------------------------------------------------------------------------------------------
     CALL apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phi(b1L,b2L,b3L,1),div)
     
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              div(i,j,k) = div(i,j,k)*dx1DM(i)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     CALL apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phi(b1L,b2L,b3L,2),com)
     
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              div(i,j,k) = div(i,j,k) + com(i,j,k)*dx2DM(j)
           END DO
        END DO
     END DO
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     
     CALL apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phi(b1L,b2L,b3L,3),com)
     
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              div(i,j,k) = div(i,j,k) + com(i,j,k)*dx3DM(k)
           END DO
        END DO
     END DO
     
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  !===========================================================================================================
  ELSE
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk,3)
                 END DO
              END DO
           END DO
        END DO
        
     !--------------------------------------------------------------------------------------------------------
     ELSE
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k,2)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence2
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE divergence_transp(m,phi,div)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        
        DO k = S3p, N3p ! Intervall ist ok!
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx1DM(i)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(1,0,S1p,S2p,S3p,N1p,N2p,N3p,S11B,S21B,S31B,N11B,N21B,N31B,N1,ndL,ndR,dimS1,cDu1CLT,cDu1CLT_LU,cDu1CRT,WDu1T,SDu1T,com,div)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,0,phi)
        
        DO k = S31B, N31B ! "B"-Grenzen etwas konsistenter, aber nicht wirklich notwendig, da ohnehin mit bc_extrapolation nachmultipliziert wird ...
           DO j = S21B, N21B
              DO i = S11B, N11B
                 div(i,j,k) = cDu1T(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    div(i,j,k) = div(i,j,k) + cDu1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx2DM(j)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(2,0,S1p,S2p,S3p,N1p,N2p,N3p,S12B,S22B,S32B,N12B,N22B,N32B,N2,ndL,ndR,dimS2,cDv2CLT,cDv2CLT_LU,cDv2CRT,WDv2T,SDv2T,com,div)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,0,phi)
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 div(i,j,k) = cDv2T(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    div(i,j,k) = div(i,j,k) + cDv2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
!pgi$ unroll = n:8
              DO i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx3DM(k)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(3,0,S1p,S2p,S3p,N1p,N2p,N3p,S13B,S23B,S33B,N13B,N23B,N33B,N3,ndL,ndR,dimS3,cDw3CLT,cDw3CLT_LU,cDw3CRT,WDw3T,SDw3T,com,div)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,0,phi)
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 div(i,j,k) = cDw3T(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    div(i,j,k) = div(i,j,k) + cDw3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE divergence_transp
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE gradient(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  
  REAL   , INTENT(INOUT) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(  OUT) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  CALL exchange(m,0,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,phi,grad)
        
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 grad(i,j,k) = grad(i,j,k)*dx1GM(i)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 grad(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 DO ii = g1L+1, g1U
                    grad(i,j,k) = grad(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) grad(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U > 0) grad(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L > 0) grad(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U > 0) grad(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L > 0) grad(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U > 0) grad(S11B:N11B,S21B:N21B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,phi,grad)
        
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 grad(i,j,k) = grad(i,j,k)*dx2GM(j)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 grad(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 DO jj = g2L+1, g2U
                    grad(i,j,k) = grad(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) grad(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U > 0) grad(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L > 0) grad(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U > 0) grad(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L > 0) grad(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U > 0) grad(S12B:N12B,S22B:N22B,N3) = 0.
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,phi,grad)
        
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 grad(i,j,k) = grad(i,j,k)*dx3GM(k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 grad(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 DO kk = g3L+1, g3U
                    grad(i,j,k) = grad(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) grad(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U > 0) grad(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L > 0) grad(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U > 0) grad(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L > 0) grad(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U > 0) grad(S13B:N13B,S23B:N23B,N3) = 0.
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE gradient_transp(m,phi,grad)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  
  REAL   , INTENT(INOUT) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(  OUT) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Umgekehrte Reihenfolge im Vergleich zu Subroutine "gradient".                             !
  !              - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) phi(0 ,S21B:N21B,S31B:N31B) = 0.
     IF (BC_1U > 0) phi(N1,S21B:N21B,S31B:N31B) = 0.
     IF (BC_2L > 0) phi(S11B:N11B,1 ,S31B:N31B) = 0.
     IF (BC_2U > 0) phi(S11B:N11B,N2,S31B:N31B) = 0.
     IF (BC_3L > 0) phi(S11B:N11B,S21B:N21B,1 ) = 0.
     IF (BC_3U > 0) phi(S11B:N11B,S21B:N21B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        
        IF (BC_1L > 0) com(0 ,S21B:N21B,S31B:N31B) = 0.
        IF (BC_1U > 0) com(N1,S21B:N21B,S31B:N31B) = 0.
        IF (BC_2L > 0) com(S11B:N11B,1 ,S31B:N31B) = 0.
        IF (BC_2U > 0) com(S11B:N11B,N2,S31B:N31B) = 0.
        IF (BC_3L > 0) com(S11B:N11B,S21B:N21B,1 ) = 0.
        IF (BC_3U > 0) com(S11B:N11B,S21B:N21B,N3) = 0.
        
        DO k = S31, N31 ! Intervall ist ok!
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 com(i,j,k) = phi(i,j,k)*dx1GM(i)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(1,1,S11,S21,S31,N11,N21,N31,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cGp1CLT,cGp1CLT_LU,cGp1CRT,WGp1T,SGp1T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 grad(i,j,k) = cGp1T(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 DO ii = d1L+1, d1U
                    grad(i,j,k) = grad(i,j,k) + cGp1T(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) phi(1 ,S22B:N22B,S32B:N32B) = 0.
     IF (BC_1U > 0) phi(N1,S22B:N22B,S32B:N32B) = 0.
     IF (BC_2L > 0) phi(S12B:N12B,0 ,S32B:N32B) = 0.
     IF (BC_2U > 0) phi(S12B:N12B,N2,S32B:N32B) = 0.
     IF (BC_3L > 0) phi(S12B:N12B,S22B:N22B,1 ) = 0.
     IF (BC_3U > 0) phi(S12B:N12B,S22B:N22B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        
        IF (BC_1L > 0) com(1 ,S22B:N22B,S32B:N32B) = 0.
        IF (BC_1U > 0) com(N1,S22B:N22B,S32B:N32B) = 0.
        IF (BC_2L > 0) com(S12B:N12B,0 ,S32B:N32B) = 0.
        IF (BC_2U > 0) com(S12B:N12B,N2,S32B:N32B) = 0.
        IF (BC_3L > 0) com(S12B:N12B,S22B:N22B,1 ) = 0.
        IF (BC_3U > 0) com(S12B:N12B,S22B:N22B,N3) = 0.
        
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 com(i,j,k) = phi(i,j,k)*dx2GM(j)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(2,2,S12,S22,S32,N12,N22,N32,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cGp2CLT,cGp2CLT_LU,cGp2CRT,WGp2T,SGp2T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO jj = d2L, d2U
                    grad(i,j,k) = grad(i,j,k) + cGp2T(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (BC_1L > 0) phi(1 ,S23B:N23B,S33B:N33B) = 0.
     IF (BC_1U > 0) phi(N1,S23B:N23B,S33B:N33B) = 0.
     IF (BC_2L > 0) phi(S13B:N13B,1 ,S33B:N33B) = 0.
     IF (BC_2U > 0) phi(S13B:N13B,N2,S33B:N33B) = 0.
     IF (BC_3L > 0) phi(S13B:N13B,S23B:N23B,0 ) = 0.
     IF (BC_3U > 0) phi(S13B:N13B,S23B:N23B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        
        IF (BC_1L > 0) com(1 ,S23B:N23B,S33B:N33B) = 0.
        IF (BC_1U > 0) com(N1,S23B:N23B,S33B:N33B) = 0.
        IF (BC_2L > 0) com(S13B:N13B,1 ,S33B:N33B) = 0.
        IF (BC_2U > 0) com(S13B:N13B,N2,S33B:N33B) = 0.
        IF (BC_3L > 0) com(S13B:N13B,S23B:N23B,0 ) = 0.
        IF (BC_3U > 0) com(S13B:N13B,S23B:N23B,N3) = 0.
        
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 com(i,j,k) = phi(i,j,k)*dx3GM(k)
              END DO
           END DO
        END DO
        
        CALL apply_compact_transp(3,3,S13,S23,S33,N13,N23,N33,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cGp3CLT,cGp3CLT_LU,cGp3CRT,WGp3T,SGp3T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     ELSE
        CALL exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
!pgi$ unroll = n:8
                 DO kk = d3L, d3U
                    grad(i,j,k) = grad(i,j,k) + cGp3T(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
     
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE gradient_transp
  
  
  
  
  
  
  
  
  
  
  !>  \brief computes \f$ \mathrm{lap_m = phi_m - mulL \Delta phi_m} \f$
  !!
  !! used for mod_rhs and for product_Helmholtz
  !! \param[in] m dimension from one to three
  !! \param[in] exch_yes indicates, if phi has to be exchanged(ghost layers)
  !! \param[inout] phi = vel in case of mod_rhs
  !! \param[out] Lap
  !! \todo change for discrete forcing
  SUBROUTINE Helmholtz(m,exch_yes,phi,Lap)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(  OUT) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  IF (exch_yes) THEN
     CALL exchange(1,m,phi)
     CALL exchange(2,m,phi)
     CALL exchange(3,m,phi)
  END IF
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_visc_yes) THEN
        IF (dimens == 3) THEN
           CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,com)
           CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,phi,dig)
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    Lap(i,j,k) = dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i)
                 END DO
              END DO
           END DO
           CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
           CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                 END DO
              END DO
           END DO
           CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
           CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                    Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
                 END DO
              END DO
           END DO
        ELSE
           CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,com)
           CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,phi,dig)
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    Lap(i,j,k) = dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i)
                 END DO
              END DO
           END DO
           CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
           CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
           DO k = S31, N31
              DO j = S21, N21
!pgi$ unroll = n:8
                 DO i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                    Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (dimens == 3) THEN
           DO k = S31, N31
              DO j = S21, N21
                 DO i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    END DO
!pgi$ unroll = n:8
                    DO kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        ELSE
           DO k = S31, N31
              DO j = S21, N21
                 DO i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_visc_yes) THEN
        IF (dimens == 3) THEN
           CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
           CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
                 END DO
              END DO
           END DO
           CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,com)
           CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,phi,dig)
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j)
                 END DO
              END DO
           END DO
           CALL apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
           CALL apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                    Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
                 END DO
              END DO
           END DO
        ELSE
           CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
           CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
                 END DO
              END DO
           END DO
           CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,com)
           CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,phi,dig)
           DO k = S32, N32
              DO j = S22, N22
!pgi$ unroll = n:8
                 DO i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j)
                    Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (dimens == 3) THEN
           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    END DO
!pgi$ unroll = n:8
                    DO kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        ELSE
           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    END DO
!pgi$ unroll = n:8
                    DO jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    END DO
                    Lap(i,j,k) = phi(i,j,k) - multL*dd1
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3 .AND. dimens == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_visc_yes) THEN
        CALL apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
        CALL apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              END DO
           END DO
        END DO
        CALL apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
        CALL apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
              END DO
           END DO
        END DO
        CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,phi,com)
        CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,phi,dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3wM(k)**2 + com(i,j,k)*ddx3wM(k)
                 Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
                 END DO
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Helmholtz_explicit(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange_all_all(.TRUE.,vel)
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (comp_visc_yes) THEN
     CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,vel(b1L,b2L,b3L,1),com)
     CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,vel(b1L,b2L,b3L,1),dig)
     DO k = S31, N31
        DO j = S21, N21
!pgi$ unroll = n:8
           DO i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i))
           END DO
        END DO
     END DO
     CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,1),com)
     CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,vel(b1L,b2L,b3L,1),dig)
     DO k = S31, N31
        DO j = S21, N21
!pgi$ unroll = n:8
           DO i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     
     CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,1),com)
     CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,vel(b1L,b2L,b3L,1),dig)
     DO k = S31, N31
        DO j = S21, N21
!pgi$ unroll = n:8
           DO i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           END DO
        END DO
     END DO
     
     END IF
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     IF (dimens == 3) THEN
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,1)
                 END DO
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 END DO
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (comp_visc_yes) THEN
     CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,2),com)
     CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,vel(b1L,b2L,b3L,2),dig)
     DO k = S32, N32
        DO j = S22, N22
!pgi$ unroll = n:8
           DO i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           END DO
        END DO
     END DO
     CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,vel(b1L,b2L,b3L,2),com)
     CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,vel(b1L,b2L,b3L,2),dig)
     DO k = S32, N32
        DO j = S22, N22
!pgi$ unroll = n:8
           DO i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j))
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     
     CALL apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,2),com)
     CALL apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,vel(b1L,b2L,b3L,2),dig)
     DO k = S32, N32
        DO j = S22, N22
!pgi$ unroll = n:8
           DO i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           END DO
        END DO
     END DO
     
     END IF
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     IF (dimens == 3) THEN
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,2)
                 END DO
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 END DO
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     IF (comp_visc_yes) THEN
        CALL apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,3),com)
        CALL apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,vel(b1L,b2L,b3L,3),dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
              END DO
           END DO
        END DO
        CALL apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,3),com)
        CALL apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,vel(b1L,b2L,b3L,3),dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
              END DO
           END DO
        END DO
        CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,vel(b1L,b2L,b3L,3),com)
        CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,vel(b1L,b2L,b3L,3),dig)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx3wM(k)**2 + com(i,j,k)*ddx3wM(k))
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,3)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,3)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*vel(i,j,k+kk,3)
                 END DO
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_explicit
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! relativ ungetestet ...
  SUBROUTINE Helmholtz_pre_explicit(exch_yes,phi,Lap)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) THEN
     CALL exchange(1,0,phi)
     CALL exchange(2,0,phi)
     CALL exchange(3,0,phi)
  END IF
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (comp_visc_yes) THEN
     CALL apply_compact(1,0,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cp1cL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
     CALL apply_compact(1,0,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           END DO
        END DO
     END DO
     CALL apply_compact(2,0,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
     CALL apply_compact(2,0,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     
     CALL apply_compact(3,0,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
     CALL apply_compact(3,0,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
     DO k = S3p, N3p
        DO j = S2p, N2p
!pgi$ unroll = n:8
           DO i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           END DO
        END DO
     END DO
     
     END IF
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     IF (dimens == 3) THEN
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                 END DO
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 END DO
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_pre_explicit
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Helmholtz_conc(exch_yes,phi,Lap) ! ok (18.02.09)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(  OUT) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  IF (exch_yes) THEN
     CALL exchange(1,0,phi)
     CALL exchange(2,0,phi)
     CALL exchange(3,0,phi)
  END IF
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (comp_visc_yes) THEN
     IF (dimens == 3) THEN
        CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),phi,com)
        CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),phi,dig)
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              DO i = S1c(m), N1c(m)
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              END DO
           END DO
        END DO
        CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),phi,com)
        CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),phi,dig)
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              DO i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
              END DO
           END DO
        END DO
        CALL apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),phi,com)
        CALL apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m),Wc33(1,0,m),Sc33(1,1,m),phi,dig)
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              DO i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                 Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
              END DO
           END DO
        END DO
     ELSE
        CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),phi,com)
        CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),phi,dig)
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              DO i = S1c(m), N1c(m)
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              END DO
           END DO
        END DO
        CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),phi,com)
        CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),phi,dig)
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              DO i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                 Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
              END DO
           END DO
        END DO
     END IF
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     IF (dimens == 3) THEN
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*phi(i,j+jj,k)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cc33(kk,k,m)*phi(i,j,k+kk)
                 END DO
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*phi(i+ii,j,k)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*phi(i,j+jj,k)
                 END DO
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_conc
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Helmholtz_conc_explicit(exch_yes,nlco)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  REAL   , INTENT(INOUT) ::  nlco(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  INTEGER                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  IF (exch_yes) THEN
     CALL exchange(1,0,conc(b1L,b2L,b3L,m))
     CALL exchange(2,0,conc(b1L,b2L,b3L,m))
     CALL exchange(3,0,conc(b1L,b2L,b3L,m))
  END IF
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  IF (comp_visc_yes) THEN
     CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),com)
     CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),conc(b1L,b2L,b3L,m),dig)
     DO k = S3c(m), N3c(m)
        DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           DO i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           END DO
        END DO
     END DO
     CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),com)
     CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),conc(b1L,b2L,b3L,m),dig)
     DO k = S3c(m), N3c(m)
        DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           DO i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           END DO
        END DO
     END DO
     
     IF (dimens == 3) THEN
     
     CALL apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),conc(b1L,b2L,b3L,m),com)
     CALL apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m),Wc33(1,0,m),Sc33(1,1,m),conc(b1L,b2L,b3L,m),dig)
     DO k = S3c(m), N3c(m)
        DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           DO i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           END DO
        END DO
     END DO
     
     END IF
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     IF (dimens == 3) THEN
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 END DO
!pgi$ unroll = n:8
                 DO kk = b3L, b3U
                    dd1 = dd1 + cc33(kk,k,m)*conc(i,j,k+kk,m)
                 END DO
                 !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*dd1
                 nlco(i,j,k) = nlco(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     ELSE
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 END DO
!pgi$ unroll = n:8
                 DO jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 END DO
                 !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*dd1
                 nlco(i,j,k) = nlco(i,j,k) - multL*dd1
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  END SUBROUTINE Helmholtz_conc_explicit
  
  
  
  
  
  
  
  
  
  
  !> \brief computes nonlinear terms
  !!
  !! first interpolates worki to pp, then \f$\partial_i\f$ vel then
  !! \f[ \mathrm{nl = pp*\partial_i vel}\f]
  !! \test Teile davon (nur zentrale Operationen!) koennten jeweils durch interpolate2_pre_vel/interpolate2_vel_pre, first_adv_pre/first_adv_vel
  !!         ersetzt werden (beachte aber Addition von nl!)
  !! \test umbenennen in advect... (?)
  SUBROUTINE nonlinear(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  REAL                   ::  dd1
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]                                         !
  !              - Feld "res" wird mehrfach ausgetauscht, um mit einem statt drei skalaren Feldern arbeiten  !
  !                zu können. Im Prinzip könnte auch rhs(:,:,:,1:3) für die Zwischenspeicherung verwendet    !
  !                werden, jedoch sind dazu einige Umbaumassnahmen notwendig bei geringem Effizienzgewinn.   !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! worki muss bereits ausgetauscht sein!
  IF (exch_yes) CALL exchange_all_all(.TRUE.,vel)
  
  
  !===========================================================================================================
  !=== u*du/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work1,pp) ! work1 -> pp, why?
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNu1U(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1U(ii,i)*vel(i+ii,j,k,1)
                 END DO
              ELSE
                 dd1 = cNu1D(n1L,i)*vel(i+n1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNu1D(ii,i)*vel(i+ii,j,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,vel(b1L,b2L,b3L,1),rr)
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 nl(i,j,k,1) = pp(i,j,k)*rr(i,j,k)*dx1uM(i)
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cu1(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cu1(ii,i)*vel(i+ii,j,k,1)
                 END DO
                 
                 nl(i,j,k,1) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*du/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,1)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,1)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,1),rr)
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 nl(i,j,k,1) = nl(i,j,k,1) + pp(i,j,k)*rr(i,j,k)*dx2pM(j)
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,1)
                 END DO
                 
                 nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*du/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,1,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,1)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,1)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,1)
                 END DO
              END IF
              
              nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,1),rr)
        DO k = S31, N31
           DO j = S21, N21
!pgi$ unroll = n:8
              DO i = S11, N11
                 nl(i,j,k,1) = nl(i,j,k,1) + pp(i,j,k)*rr(i,j,k)*dx3pM(k)
              END DO
           END DO
        END DO
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 dd1 = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,1)
                 END DO
                 
                 nl(i,j,k,1) = nl(i,j,k,1) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  
  !===========================================================================================================
  !=== u*dv/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,2)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,2),rr)
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 nl(i,j,k,2) = pp(i,j,k)*rr(i,j,k)*dx1pM(i)
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,2)
                 END DO
                 
                 nl(i,j,k,2) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dv/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNv2U(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2U(jj,j)*vel(i,j+jj,k,2)
                 END DO
              ELSE
                 dd1 = cNv2D(n2L,j)*vel(i,j+n2L,k,2)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNv2D(jj,j)*vel(i,j+jj,k,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,vel(b1L,b2L,b3L,2),rr)
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 nl(i,j,k,2) = nl(i,j,k,2) + pp(i,j,k)*rr(i,j,k)*dx2vM(j)
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cv2(b2L,j)*vel(i,j+b2L,k,2)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cv2(jj,j)*vel(i,j+jj,k,2)
                 END DO
                 
                 nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== w*dv/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,2,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp3U(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*vel(i,j,k+kk,2)
                 END DO
              ELSE
                 dd1 = cNp3D(n3L,k)*vel(i,j,k+n3L,2)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*vel(i,j,k+kk,2)
                 END DO
              END IF
              
              nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,2),rr)
        DO k = S32, N32
           DO j = S22, N22
!pgi$ unroll = n:8
              DO i = S12, N12
                 nl(i,j,k,2) = nl(i,j,k,2) + pp(i,j,k)*rr(i,j,k)*dx3pM(k)
              END DO
           END DO
        END DO
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 dd1 = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*vel(i,j,k+kk,2)
                 END DO
                 
                 nl(i,j,k,2) = nl(i,j,k,2) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  END IF
  
  
  
  
  
  IF (dimens == 3) THEN
  !===========================================================================================================
  !=== u*dw/dx ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp1U(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*vel(i+ii,j,k,3)
                 END DO
              ELSE
                 dd1 = cNp1D(n1L,i)*vel(i+n1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*vel(i+ii,j,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,3),rr)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = pp(i,j,k)*rr(i,j,k)*dx1pM(i)
              END DO
           END DO
        END DO
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*vel(i+ii,j,k,3)
                 END DO
                 
                 nl(i,j,k,3) = dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dw/dy ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNp2U(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*vel(i,j+jj,k,3)
                 END DO
              ELSE
                 dd1 = cNp2D(n2L,j)*vel(i,j+n2L,k,3)
!pgi$ unroll = n:8
                 DO jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*vel(i,j+jj,k,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,3),rr)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) + pp(i,j,k)*rr(i,j,k)*dx2pM(j)
              END DO
           END DO
        END DO
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*vel(i,j+jj,k,3)
                 END DO
                 
                 nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== w*dw/dz ===============================================================================================
  !===========================================================================================================
  CALL interpolate2_pre_vel(.FALSE.,3,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  IF (upwind_yes) THEN
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (pp(i,j,k) >= 0.) THEN
                 dd1 = cNw3U(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3U(kk,k)*vel(i,j,k+kk,3)
                 END DO
              ELSE
                 dd1 = cNw3D(n3L,k)*vel(i,j,k+n3L,3)
!pgi$ unroll = n:8
                 DO kk = n3L+1, n3U
                    dd1 = dd1 + cNw3D(kk,k)*vel(i,j,k+kk,3)
                 END DO
              END IF
              
              nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     IF (comp_conv_yes) THEN
        CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,vel(b1L,b2L,b3L,3),rr)
        DO k = S33, N33
           DO j = S23, N23
!pgi$ unroll = n:8
              DO i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) + pp(i,j,k)*rr(i,j,k)*dx3wM(k)
              END DO
           END DO
        END DO
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 dd1 = cw3(b3L,k)*vel(i,j,k+b3L,3)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    dd1 = dd1 + cw3(kk,k)*vel(i,j,k+kk,3)
                 END DO
                 
                 nl(i,j,k,3) = nl(i,j,k,3) + dd1*pp(i,j,k)
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE nonlinear
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE nonlinear_conc(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  INTEGER                ::  m
  REAL                   ::  dd1, dd2
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - 
  !----------------------------------------------------------------------------------------------------------!
  
  m = conc_nu
  
  IF (exch_yes) THEN
     CALL exchange(1,0,conc(b1L,b2L,b3L,m))
     CALL exchange(2,0,conc(b1L,b2L,b3L,m))
     CALL exchange(3,0,conc(b1L,b2L,b3L,m))
  END IF
  
  !===========================================================================================================
  IF (dimens == 3) THEN
     
     IF (upwind_conc_yes) THEN
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(1,m) + work1(i,j,k)
                 
                 IF (dd2 >= 0.) THEN
                    dd1 = cNc1U(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNc1U(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                 ELSE
                    dd1 = cNc1D(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNc1D(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                 END IF
                 
                 nlco(i,j,k,m) = dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(2,m) + work2(i,j,k)
                 
                 IF (dd2 >= 0.) THEN
                    dd1 = cNc2U(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNc2U(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                 ELSE
                    dd1 = cNc2D(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNc2D(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                 END IF
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(3,m) + work3(i,j,k)
                 
                 IF (dd2 >= 0.) THEN
                    dd1 = cNc3U(n3L,k,m)*conc(i,j,k+n3L,m)
!pgi$ unroll = n:8
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNc3U(kk,k,m)*conc(i,j,k+kk,m)
                    END DO
                 ELSE
                    dd1 = cNc3D(n3L,k,m)*conc(i,j,k+n3L,m)
!pgi$ unroll = n:8
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNc3D(kk,k,m)*conc(i,j,k+kk,m)
                    END DO
                 END IF
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
              END DO
           END DO
        END DO
     ELSE
        IF (comp_conv_yes) THEN
           CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),pp)
           CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),rr)
           CALL apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),conc(b1L,b2L,b3L,m),Ap)
           DO k = S3c(m), N3c(m)
              DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
                 DO i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    nlco(i,j,k,m) =                 pp(i,j,k)*dx1pM(i)*(us_vec(1,m) + work1(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + rr(i,j,k)*dx2pM(j)*(us_vec(2,m) + work2(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + Ap(i,j,k)*dx3pM(k)*(us_vec(3,m) + work3(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 END DO
              END DO
           END DO
        ELSE
           DO k = S3c(m), N3c(m)
              DO j = S2c(m), N2c(m)
                 DO i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cc1(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                    
                    nlco(i,j,k,m) = dd1*(us_vec(1,m) + work1(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cc2(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cc3(kk,k,m)*conc(i,j,k+kk,m)
                    END DO
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(3,m) + work3(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 END DO
              END DO
           END DO
        END IF
     END IF
  !===========================================================================================================
  ELSE
     IF (upwind_conc_yes) THEN
        DO k = S3c(m), N3c(m)
           DO j = S2c(m), N2c(m)
              DO i = S1c(m), N1c(m)
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(1,m) + work1(i,j,k)
                 
                 IF (dd2 >= 0.) THEN
                    dd1 = cNc1U(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNc1U(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                 ELSE
                    dd1 = cNc1D(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNc1D(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                 END IF
                 
                 nlco(i,j,k,m) = dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(2,m) + work2(i,j,k)
                 
                 IF (dd2 >= 0.) THEN
                    dd1 = cNc2U(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNc2U(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                 ELSE
                    dd1 = cNc2D(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNc2D(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                 END IF
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
              END DO
           END DO
        END DO
     ELSE
        IF (comp_conv_yes) THEN
           CALL apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),pp)
           CALL apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),rr)
           DO k = S3c(m), N3c(m)
              DO j = S2c(m), N2c(m)
!pgi$ unroll = n:8
                 DO i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    nlco(i,j,k,m) =                 pp(i,j,k)*dx1pM(i)*(us_vec(1,m) + work1(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + rr(i,j,k)*dx2pM(j)*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 END DO
              END DO
           END DO
        ELSE
           DO k = S3c(m), N3c(m)
              DO j = S2c(m), N2c(m)
                 DO i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cc1(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                    
                    nlco(i,j,k,m) = dd1*(us_vec(1,m) + work1(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cc2(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 END DO
              END DO
           END DO
        END IF
     END IF
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE nonlinear_conc
  
  
  
  
  
  
  
  
  
  !> \brief interpolates velocity on pressure grid
  !!
  !! vel(:,:,:,i) -> worki(:,:,:)
  !! \test anderes Modul?
  !! \test basiert neu auf interpolate2_vel_pre ... ok??
  SUBROUTINE interpolate_vel(exch_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
                   CALL interpolate2_vel_pre(exch_yes,1,vel(b1L,b2L,b3L,1),work1)
                   CALL interpolate2_vel_pre(exch_yes,2,vel(b1L,b2L,b3L,2),work2)
  IF (dimens == 3) CALL interpolate2_vel_pre(exch_yes,3,vel(b1L,b2L,b3L,3),work3)
  
  
  CALL exchange(1,0,work1)
  CALL exchange(2,0,work1)
  CALL exchange(3,0,work1)
  
  CALL exchange(1,0,work2)
  CALL exchange(2,0,work2)
  CALL exchange(3,0,work2)
  
  IF (dimens == 3) THEN
     CALL exchange(1,0,work3)
     CALL exchange(2,0,work3)
     CALL exchange(3,0,work3)
  END IF
  
  
  END SUBROUTINE interpolate_vel
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN   ) ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(INOUT) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !              - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen (Extrapolation!).                            !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(1,0,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_1L > 0 .AND. SS1 == S11B) THEN
           i = S11B
           DO k = SS3, NN3
              DO j = SS2, NN2
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_1U > 0 .AND. NN1 == N11B) THEN
           i = N11B
           DO k = SS3, NN3
              DO j = SS2, NN2
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11B, N11B ! TEST!!! hier koennte man auch SS1, NN1 nehmen! gilt auch fuer andere Routinen!!
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(2,0,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_2L > 0 .AND. SS2 == S22B) THEN
           j = S22B
           DO k = SS3, NN3
              DO i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_2U > 0 .AND. NN2 == N22B) THEN
           j = N22B
           DO k = SS3, NN3
              DO i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     ELSE
        DO k = SS3, NN3
           DO j = S22B, N22B ! TEST!!! hier koennte man auch SS2, NN2 nehmen!
              DO i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(3,0,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_3L > 0 .AND. SS3 == S33B) THEN
           k = S33B
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_3U > 0 .AND. NN3 == N33B) THEN
           k = N33B
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     ELSE
        DO k = S33B, N33B ! TEST!!! hier koennte man auch SS3, NN3 nehmen!
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_pre_vel
  
  
  
  
  
  
  
  
  
  
  !> \brief interpolates pressure to vel nodes
  !!
  !! \f[ inter = inter + interpolated(phi)
  !! Wie interpolate_pre_vel, allerdings mit fixen Index-Limiten (ohne Rand)
  !! \param[in] exch_yes indicates if fields have exchanged first
  !! \param[in] m dimension
  !! \param[inout] phi input field
  !! \param[inout] inter output field
  SUBROUTINE interpolate2_pre_vel(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN   ) ::  m
  
  REAL   , INTENT(INOUT) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  IF (exch_yes) CALL exchange(m,0,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,phi,inter)
     ELSE
        DO k = S31, N31
           DO j = S21, N21
              DO i = S11, N11
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,phi,inter)
     ELSE
        DO k = S32, N32
           DO j = S22, N22
              DO i = S12, N12
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,phi,inter)
     ELSE
        DO k = S33, N33
           DO j = S23, N23
              DO i = S13, N13
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE interpolate_vel_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN)    ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(INOUT) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(1,1,S1p,SS2,SS3,N1p,NN2,NN3,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,phi,inter)
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(2,2,SS1,S2p,SS3,NN1,N2p,NN3,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,phi,inter)
     ELSE
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(3,3,SS1,SS2,S3p,NN1,NN2,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,phi,inter)
     ELSE
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate_vel_pre
  
  
  
  
  
  
  
  
  
  !> \brief interpolates velocity to pressure grid
  !!
  !! Wie interpolate_vel_pre, allerdings mit fixen Index-Limiten (ohne Rand)
  SUBROUTINE interpolate2_vel_pre(exch_yes,m,phi,inter)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN)    ::  m
  
  REAL   , INTENT(INOUT) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  IF (exch_yes) CALL exchange(m,m,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,phi,inter)
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,phi,inter)
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (comp_inter_yes) THEN
        CALL apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,phi,inter)
     ELSE
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interpolate2_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE interpolate_conc(exch_yes) ! ok (18.02.09)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1, dd2
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -  
  !----------------------------------------------------------------------------------------------------------!
  
  
  DO m = 1, n_conc
     
     !========================================================================================================
     IF (gravity(1) /= 0.) THEN
        
        dd1 = Ric(m)*gravity(1)
        
        IF (exch_yes) CALL exchange(1,0,conc(b1L,b2L,b3L,m))
        
        IF (comp_inter_yes) THEN
           CALL apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIcuCL(-ndL,0,m),cIcuCL_LU(1,0,m),cIcuCR(-ndR,0,m),WIcu(1,0,m),SIcu(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S11:N11,S21:N21,S31:N31,1) = nl(S11:N11,S21:N21,S31:N31,1) - dd1*pp(S11:N11,S21:N21,S31:N31)
        ELSE
           DO k = S31, N31
              DO j = S21, N21
                 DO i = S11, N11
                    dd2 = cIcu(g1L,i,m)*conc(i+g1L,j,k,m)
!pgi$ unroll = n:8
                    DO ii = g1L+1, g1U
                       dd2 = dd2 + cIcu(ii,i,m)*conc(i+ii,j,k,m)
                    END DO
                    nl(i,j,k,1) = nl(i,j,k,1) - dd1*dd2
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (gravity(2) /= 0.) THEN
        
        dd1 = Ric(m)*gravity(2)
        
        IF (exch_yes) CALL exchange(2,0,conc(b1L,b2L,b3L,m))
        
        IF (comp_inter_yes) THEN
           CALL apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIcvCL(-ndL,0,m),cIcvCL_LU(1,0,m),cIcvCR(-ndR,0,m),WIcv(1,0,m),SIcv(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S12:N12,S22:N22,S32:N32,2) = nl(S12:N12,S22:N22,S32:N32,2) - dd1*pp(S12:N12,S22:N22,S32:N32)
        ELSE
           DO k = S32, N32
              DO j = S22, N22
                 DO i = S12, N12
                    dd2 = cIcv(g2L,j,m)*conc(i,j+g2L,k,m)
!pgi$ unroll = n:8
                    DO jj = g2L+1, g2U
                       dd2 = dd2 + cIcv(jj,j,m)*conc(i,j+jj,k,m)
                    END DO
                    nl(i,j,k,2) = nl(i,j,k,2) - dd1*dd2
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (gravity(3) /= 0. .AND. dimens == 3) THEN
        
        dd1 = Ric(m)*gravity(3)
        
        IF (exch_yes) CALL exchange(3,0,conc(b1L,b2L,b3L,m))
        
        IF (comp_inter_yes) THEN
           CALL apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIcwCL(-ndL,0,m),cIcwCL_LU(1,0,m),cIcwCR(-ndR,0,m),WIcw(1,0,m),SIcw(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S13:N13,S23:N23,S33:N33,3) = nl(S13:N13,S23:N23,S33:N33,3) - dd1*pp(S13:N13,S23:N23,S33:N33)
        ELSE
           DO k = S33, N33
              DO j = S23, N23
                 DO i = S13, N13
                    dd2 = cIcw(g3L,k,m)*conc(i,j,k+g3L,m)
!pgi$ unroll = n:8
                    DO kk = g3L+1, g3U
                       dd2 = dd2 + cIcw(kk,k,m)*conc(i,j,k+kk,m)
                    END DO
                    nl(i,j,k,3) = nl(i,j,k,3) - dd1*dd2
                 END DO
              END DO
           END DO
        END IF
     END IF
     !========================================================================================================
     
  END DO
  
  
  END SUBROUTINE interpolate_conc
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN   ) ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen.                                             !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(1,0,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,phi,der)
        
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11, N11
                 der(i,j,k) = der(i,j,k)*dx1GM(i)
              END DO
           END DO
        END DO
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_1L > 0 .AND. SS1 == S11B) THEN
           i = S11B
           DO k = SS3, NN3
              DO j = SS2, NN2
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_1U > 0 .AND. NN1 == N11B) THEN
           i = N11B
           DO k = SS3, NN3
              DO j = SS2, NN2
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11B, N11B
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 DO ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(2,0,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,phi,der)
        
        DO k = SS3, NN3
           DO j = S22, N22
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx2GM(j)
              END DO
           END DO
        END DO
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_2L > 0 .AND. SS2 == S22B) THEN
           j = S22B
           DO k = SS3, NN3
              DO i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_2U > 0 .AND. NN2 == N22B) THEN
           j = N22B
           DO k = SS3, NN3
              DO i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = SS3, NN3
           DO j = S22B, N22B
              DO i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 DO jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_grad_yes) THEN
        CALL apply_compact(3,0,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,phi,der)
        
        DO k = S33, N33
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx3GM(k)
              END DO
           END DO
        END DO
        !+++++++++++++++++++++++++++++++++++++++
        IF (BC_3L > 0 .AND. SS3 == S33B) THEN
           k = S33B
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END IF
        IF (BC_3U > 0 .AND. NN3 == N33B) THEN
           k = N33B
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END IF
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S33B, N33B
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 DO kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_vel_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  
  INTEGER, INTENT(IN   ) ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(1,1,S1p,SS2,SS3,N1p,NN2,NN3,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phi,der)
        
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 der(i,j,k) = der(i,j,k)*dx1DM(i)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 der(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
                 DO ii = d1L+1, d1U
                    der(i,j,k) = der(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(2,2,SS1,S2p,SS3,NN1,N2p,NN3,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phi,der)
        
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx2DM(j)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 der(i,j,k) = cDv2(d2L,j)*phi(i,j+d2L,k)
                 DO jj = d2L+1, d2U
                    der(i,j,k) = der(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (comp_div_yes) THEN
        CALL apply_compact(3,3,SS1,SS2,S3p,NN1,NN2,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phi,der)
        
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx3DM(k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 der(i,j,k) = cDw3(d3L,k)*phi(i,j,k+d3L)
                 DO kk = d3L+1, d3U
                    der(i,j,k) = der(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_adv_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN   ) ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  LOGICAL, INTENT(IN   ) ::  upwind_yes
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(IN   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - die Punkte auf dem Rand der Wand-normalen Komponente werden auch bei kompakter            !
  !                Differenzierung im Feld immer explizit gerechnet, um nur eine Variante kompakter          !
  !                Differenzen abspeichern zu muessen.                                                       !
  !              - fuer die Punkte der Wand-normalen Komponente sind S12, N12 nur Platzhalter fuer S1p+1,    !
  !                N1p-1 usw.                                                                                !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S1p, N1p
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNp1U(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNp1U(ii,i)*phi(i+ii,j,k)
                    END DO
                 ELSE
                    dd1 = cNp1D(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNp1D(ii,i)*phi(i+ii,j,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           CALL apply_compact(1,0,S12,SS2,SS3,N12,NN2,NN3,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,der)
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S12, N12
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx1pM(i)
                 END DO
              END DO
           END DO
           !+++++++++++++++++++++++++++++++++++++++
           IF (BC_1L > 0 .AND. SS1 == S1p) THEN
              i = S1p
              DO k = SS3, NN3
                 DO j = SS2, NN2
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           IF (BC_1U > 0 .AND. NN1 == N1p) THEN
              i = N1p
              DO k = SS3, NN3
                 DO j = SS2, NN2
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           !+++++++++++++++++++++++++++++++++++++++
        ELSE
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S1p, N1p
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = S2p, N2p
              DO i = SS1, NN1
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNp2U(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNp2U(jj,j)*phi(i,j+jj,k)
                    END DO
                 ELSE
                    dd1 = cNp2D(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNp2D(jj,j)*phi(i,j+jj,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           CALL apply_compact(2,0,SS1,S23,SS3,NN1,N23,NN3,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,der)
           DO k = SS3, NN3
              DO j = S23, N23
                 DO i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx2pM(j)
                 END DO
              END DO
           END DO
           !+++++++++++++++++++++++++++++++++++++++
           IF (BC_2L > 0 .AND. SS2 == S2p) THEN
              j = S2p
              DO k = SS3, NN3
                 DO i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           IF (BC_2U > 0 .AND. NN2 == N2p) THEN
              j = N2p
              DO k = SS3, NN3
                 DO i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           !+++++++++++++++++++++++++++++++++++++++
        ELSE
           DO k = SS3, NN3
              DO j = S2p, N2p
                 DO i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = S3p, N3p
           DO j = SS2, NN2
              DO i = SS1, NN1
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNp3U(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNp3U(kk,k)*phi(i,j,k+kk)
                    END DO
                 ELSE
                    dd1 = cNp3D(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNp3D(kk,k)*phi(i,j,k+kk)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           CALL apply_compact(3,0,SS1,SS2,S32,NN1,NN2,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,der)
           DO k = S32, N32
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx3pM(k)
                 END DO
              END DO
           END DO
           !+++++++++++++++++++++++++++++++++++++++
           IF (BC_3L > 0 .AND. SS3 == S3p) THEN
              k = S3p
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           IF (BC_3U > 0 .AND. NN3 == N3p) THEN
              k = N3p
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END IF
           !+++++++++++++++++++++++++++++++++++++++
        ELSE
           DO k = S3p, N3p
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_adv_pre
  
  
  
  
  
  
  
  
  
  ! TEST!!! Routine kann Punkte der Rand-normalen Komponente auf dem Rand NICHT behandeln!! (will sie momentan aber auch gar nicht koennen)
  ! TEST!!! noch relativ ungetestet!
  SUBROUTINE first_adv_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  IMPLICIT NONE
  
  LOGICAL, INTENT(IN   ) ::  exch_yes
  INTEGER, INTENT(IN   ) ::  m
  INTEGER, INTENT(IN   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  LOGICAL, INTENT(IN   ) ::  upwind_yes
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(INOUT) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(IN   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (exch_yes) CALL exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  IF (m == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = SS2, NN2
              DO i = S11, N11
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNu1U(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNu1U(ii,i)*phi(i+ii,j,k)
                    END DO
                 ELSE
                    dd1 = cNu1D(n1L,i)*phi(i+n1L,j,k)
                    DO ii = n1L+1, n1U
                       dd1 = dd1 + cNu1D(ii,i)*phi(i+ii,j,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           CALL apply_compact(1,1,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,der)
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S11, N11
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx1uM(i)
                 END DO
              END DO
           END DO
        ELSE
           DO k = SS3, NN3
              DO j = SS2, NN2
                 DO i = S11, N11
                    dd1 = cu1(b1L,i)*phi(i+b1L,j,k)
                    DO ii = b1L+1, b1U
                       dd1 = dd1 + cu1(ii,i)*phi(i+ii,j,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  END IF
  !===========================================================================================================
  IF (m == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = SS3, NN3
           DO j = S22, N22
              DO i = SS1, NN1
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNv2U(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNv2U(jj,j)*phi(i,j+jj,k)
                    END DO
                 ELSE
                    dd1 = cNv2D(n2L,j)*phi(i,j+n2L,k)
                    DO jj = n2L+1, n2U
                       dd1 = dd1 + cNv2D(jj,j)*phi(i,j+jj,k)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           CALL apply_compact(2,2,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,der)
           DO k = SS3, NN3
              DO j = S22, N22
                 DO i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx2vM(j)
                 END DO
              END DO
           END DO
        ELSE
           DO k = SS3, NN3
              DO j = S22, N22
                 DO i = SS1, NN1
                    dd1 = cv2(b2L,j)*phi(i,j+b2L,k)
                    DO jj = b2L+1, b2U
                       dd1 = dd1 + cv2(jj,j)*phi(i,j+jj,k)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (m == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (upwind_yes) THEN
        DO k = S33, N33
           DO j = SS2, NN2
              DO i = SS1, NN1
                 IF (adv(i,j,k) >= 0.) THEN
                    dd1 = cNw3U(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNw3U(kk,k)*phi(i,j,k+kk)
                    END DO
                 ELSE
                    dd1 = cNw3D(n3L,k)*phi(i,j,k+n3L)
                    DO kk = n3L+1, n3U
                       dd1 = dd1 + cNw3D(kk,k)*phi(i,j,k+kk)
                    END DO
                 END IF
                 der(i,j,k) = dd1*adv(i,j,k)
              END DO
           END DO
        END DO
     !--------------------------------------------------------------------------------------------------------
     ELSE
        IF (comp_conv_yes) THEN
           CALL apply_compact(3,3,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,phi,der)
           DO k = S33, N33
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx3wM(k)
                 END DO
              END DO
           END DO
        ELSE
           DO k = S33, N33
              DO j = SS2, NN2
                 DO i = SS1, NN1
                    dd1 = cw3(b3L,k)*phi(i,j,k+b3L)
                    DO kk = b3L+1, b3U
                       dd1 = dd1 + cw3(kk,k)*phi(i,j,k+kk)
                    END DO
                    der(i,j,k) = dd1*adv(i,j,k)
                 END DO
              END DO
           END DO
        END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE first_adv_vel
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE outflow_bc
  
  ! (revised on 06.08.2009)
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  REAL                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Es wird vorausgesetzt, dass vel(:,:,:,:) zuvor schon ausgetauscht, bzw. an den Raendern   !
  !                zu Null gesetzt wurde (was beides streng genommen aber nichtmal notwentwendig ist).       !
  !              - Kompakte Differenzen nicht notwendig, ist ohnehin nur ein Modell.                         !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== vel(:,:,:,1) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,1) .AND. (BC_2L == 1 .OR. BC_2L == 3)) THEN ! TEST!!! Warum eigentlich BC_2L == 3? Wuerde generell nicht outlet(2,1,1) ausreichen?
     j = 1
     DO k = S31B, N31B
        DO i = S11B, N11B
           nlbc12(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc12(i,k,1) = nlbc12(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,1)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,1)
           END DO
           nlbc12(i,k,1) = dd1*nlbc12(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,1) .AND. (BC_2U == 1 .OR. BC_2U == 3)) THEN
     j = N2
     DO k = S31B, N31B
        DO i = S11B, N11B
           nlbc12(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc12(i,k,2) = nlbc12(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,2)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,2)
           END DO
           nlbc12(i,k,2) = dd1*nlbc12(i,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,1) .AND. (BC_3L == 1 .OR. BC_3L == 3)) THEN
     k = 1
     DO j = S21B, N21B
        DO i = S11B, N11B
           nlbc13(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc13(i,j,1) = nlbc13(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,1)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,1)
           END DO
           nlbc13(i,j,1) = dd1*nlbc13(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,1) .AND. (BC_3U == 1 .OR. BC_3U == 3)) THEN
     k = N3
     DO j = S21B, N21B
        DO i = S11B, N11B
           nlbc13(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc13(i,j,2) = nlbc13(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,1)
           END DO
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,2)
!pgi$ unroll = n:8
           DO ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,2)
           END DO
           nlbc13(i,j,2) = dd1*nlbc13(i,j,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,1) .AND. (BC_1L == 1 .OR. BC_1L == 2)) THEN
     i = 0
     DO k = S31B, N31B
        DO j = S21B, N21B
           nlbc11(j,k,1) = cDu1(d1L,1)*vel(1+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nlbc11(j,k,1) = nlbc11(j,k,1) + cDu1(ii,1)*vel(1+ii,j,k,1)
           END DO
           nlbc11(j,k,1) = drift1(j,k,1)*nlbc11(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,1) .AND. (BC_1U == 1 .OR. BC_1U == 2)) THEN
     i = N1
     DO k = S31B, N31B
        DO j = S21B, N21B
           nlbc11(j,k,2) = cDu1(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           DO ii = d1L+1, d1U
              nlbc11(j,k,2) = nlbc11(j,k,2) + cDu1(ii,i)*vel(i+ii,j,k,1)
           END DO
           nlbc11(j,k,2) = drift1(j,k,2)*nlbc11(j,k,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,2) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,2) .AND. (BC_1L == 1 .OR. BC_1L == 3)) THEN
     i = 1
     DO k = S32B, N32B
        DO j = S22B, N22B
           nlbc21(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc21(j,k,1) = nlbc21(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,1)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,1)
           END DO
           nlbc21(j,k,1) = dd1*nlbc21(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,2) .AND. (BC_1U == 1 .OR. BC_1U == 3)) THEN
     i = N1
     DO k = S32B, N32B
        DO j = S22B, N22B
           nlbc21(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc21(j,k,2) = nlbc21(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,2)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,2)
           END DO
           nlbc21(j,k,2) = dd1*nlbc21(j,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,2) .AND. (BC_3L == 1 .OR. BC_3L == 3)) THEN
     k = 1
     DO j = S22B, N22B
        DO i = S12B, N12B
           nlbc23(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc23(i,j,1) = nlbc23(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,1)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,1)
           END DO
           nlbc23(i,j,1) = dd1*nlbc23(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,2) .AND. (BC_3U == 1 .OR. BC_3U == 3)) THEN
     k = N3
     DO j = S22B, N22B
        DO i = S12B, N12B
           nlbc23(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              nlbc23(i,j,2) = nlbc23(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,2)
           END DO
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,2)
!pgi$ unroll = n:8
           DO jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,2)
           END DO
           nlbc23(i,j,2) = dd1*nlbc23(i,j,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,2) .AND. (BC_2L == 1 .OR. BC_2L == 2)) THEN
     j = 0
     DO k = S32B, N32B
        DO i = S12B, N12B
           nlbc22(i,k,1) = cDv2(d2L,1)*vel(i,1+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nlbc22(i,k,1) = nlbc22(i,k,1) + cDv2(jj,1)*vel(i,1+jj,k,2)
           END DO
           nlbc22(i,k,1) = drift2(i,k,1)*nlbc22(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,2) .AND. (BC_2U == 1 .OR. BC_2U == 2)) THEN
     j = N2
     DO k = S32B, N32B
        DO i = S12B, N12B
           nlbc22(i,k,2) = cDv2(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           DO jj = d2L+1, d2U
              nlbc22(i,k,2) = nlbc22(i,k,2) + cDv2(jj,j)*vel(i,j+jj,k,2)
           END DO
           nlbc22(i,k,2) = drift2(i,k,2)*nlbc22(i,k,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,3) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(1,1,3) .AND. (BC_1L == 1 .OR. BC_1L == 3)) THEN
     i = 1
     DO k = S33B, N33B
        DO j = S23B, N23B
           nlbc31(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc31(j,k,1) = nlbc31(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,1)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,1)
           END DO
           nlbc31(j,k,1) = dd1*nlbc31(j,k,1)
        END DO
     END DO
  END IF
  IF (outlet(1,2,3) .AND. (BC_1U == 1 .OR. BC_1U == 3)) THEN
     i = N1
     DO k = S33B, N33B
        DO j = S23B, N23B
           nlbc31(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              nlbc31(j,k,2) = nlbc31(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,2)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,2)
           END DO
           nlbc31(j,k,2) = dd1*nlbc31(j,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(2,1,3) .AND. (BC_2L == 1 .OR. BC_2L == 3)) THEN
     j = 1
     DO k = S33B, N33B
        DO i = S13B, N13B
           nlbc32(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc32(i,k,1) = nlbc32(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,1)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,1)
           END DO
           nlbc32(i,k,1) = dd1*nlbc32(i,k,1)
        END DO
     END DO
  END IF
  IF (outlet(2,2,3) .AND. (BC_2U == 1 .OR. BC_2U == 3)) THEN
     j = N2
     DO k = S33B, N33B
        DO i = S13B, N13B
           nlbc32(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              nlbc32(i,k,2) = nlbc32(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,3)
           END DO
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,2)
!pgi$ unroll = n:8
           DO kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,2)
           END DO
           nlbc32(i,k,2) = dd1*nlbc32(i,k,2)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  IF (outlet(3,1,3) .AND. (BC_3L == 1 .OR. BC_3L == 2)) THEN
     k = 0
     DO j = S23B, N23B
        DO i = S13B, N13B
           nlbc33(i,j,1) = cDw3(d3L,1)*vel(i,j,1+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nlbc33(i,j,1) = nlbc33(i,j,1) + cDw3(kk,1)*vel(i,j,1+kk,3)
           END DO
           nlbc33(i,j,1) = drift3(i,j,1)*nlbc33(i,j,1)
        END DO
     END DO
  END IF
  IF (outlet(3,2,3) .AND. (BC_3U == 1 .OR. BC_3U == 2)) THEN
     k = N3
     DO j = S23B, N23B
        DO i = S13B, N13B
           nlbc33(i,j,2) = cDw3(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           DO kk = d3L+1, d3U
              nlbc33(i,j,2) = nlbc33(i,j,2) + cDw3(kk,k)*vel(i,j,k+kk,3)
           END DO
           nlbc33(i,j,2) = drift3(i,j,2)*nlbc33(i,j,2)
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE outflow_bc
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE sediment_bc
  
  ! (revised on 06.08.2009)
  
  IMPLICIT NONE
  
  INTEGER                ::  m
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  !===========================================================================================================
  IF (BCc_1L(m) == 1 .AND. (.NOT. isopycnal(1,1,m))) THEN
     i = 1
     DO k = S3p, N3p
        DO j = S2p, N2p
           sed1L(j,k,m) = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              sed1L(j,k,m) = sed1L(j,k,m) + cc1(ii,i,m)*conc(i+ii,j,k,m)
           END DO
           sed1L(j,k,m) = sed1L(j,k,m)*(us_vec(1,m) + work1(i,j,k))
        END DO
     END DO
  END IF
  
  IF (BCc_1U(m) == 1 .AND. (.NOT. isopycnal(1,2,m))) THEN
     i = N1
     DO k = S3p, N3p
        DO j = S2p, N2p
           sed1U(j,k,m) = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
           DO ii = b1L+1, b1U
              sed1U(j,k,m) = sed1U(j,k,m) + cc1(ii,i,m)*conc(i+ii,j,k,m)
           END DO
           sed1U(j,k,m) = sed1U(j,k,m)*(us_vec(1,m) + work1(i,j,k))
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BCc_2L(m) == 1 .AND. (.NOT. isopycnal(2,1,m))) THEN
     j = 1
     DO k = S3p, N3p
        DO i = S1p, N1p
           sed2L(i,k,m) = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              sed2L(i,k,m) = sed2L(i,k,m) + cc2(jj,j,m)*conc(i,j+jj,k,m)
           END DO
           sed2L(i,k,m) = sed2L(i,k,m)*(us_vec(2,m) + work2(i,j,k))
        END DO
     END DO
  END IF
  
  IF (BCc_2U(m) == 1 .AND. (.NOT. isopycnal(2,2,m))) THEN
     j = N2
     DO k = S3p, N3p
        DO i = S1p, N1p
           sed2U(i,k,m) = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
           DO jj = b2L+1, b2U
              sed2U(i,k,m) = sed2U(i,k,m) + cc2(jj,j,m)*conc(i,j+jj,k,m)
           END DO
           sed2U(i,k,m) = sed2U(i,k,m)*(us_vec(2,m) + work2(i,j,k))
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BCc_3L(m) == 1 .AND. (.NOT. isopycnal(3,1,m))) THEN
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           sed3L(i,j,m) = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              sed3L(i,j,m) = sed3L(i,j,m) + cc3(kk,k,m)*conc(i,j,k+kk,m)
           END DO
           sed3L(i,j,m) = sed3L(i,j,m)*(us_vec(3,m) + work3(i,j,k))
        END DO
     END DO
  END IF
  
  IF (BCc_3U(m) == 1 .AND. (.NOT. isopycnal(3,2,m))) THEN
     k = N3
     DO j = S2p, N2p
        DO i = S1p, N1p
           sed3U(i,j,m) = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
           DO kk = b3L+1, b3U
              sed3U(i,j,m) = sed3U(i,j,m) + cc3(kk,k,m)*conc(i,j,k+kk,m)
           END DO
           sed3U(i,j,m) = sed3U(i,j,m)*(us_vec(3,m) + work3(i,j,k))
        END DO
     END DO
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE sediment_bc
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE filter(fil_dir,dir,m,phi,fil)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  fil_dir
  INTEGER, INTENT(IN   ) ::  dir
  INTEGER, INTENT(IN   ) ::  m
  
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , INTENT(  OUT) ::  fil(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Prinzipiell muss das gesamte Feld gefiltert werden; lediglich Dirichlet-Randbedingungen   !
  !                können/sollten unverändert bleiben. Wegen des versetzten Gitters ist das jedoch bei Wand- !
  !                normalen Geschwindigkeitskomponenten nicht möglich und könnte nur durch eine vor- und     !
  !                nachgeschaltete Interpolation auf das Druckgitter erreicht werden.                        !
  !              - Um jede Art der Randbedingungen zuzulassen, wird zunächst auf dem gesamten Feld gefiltert !
  !                und erst anschliessend mit Dirichlet-Randbedingungen überschrieben. Prinzipiell könnten   !
  !                aber auch diese Randbedingungen gefiltert werden.                                         !
  !----------------------------------------------------------------------------------------------------------!
  
  
  CALL exchange(fil_dir,dir,phi)
  
  
  !===========================================================================================================
  !=== Druck (NICHT Konzentrationen) =========================================================================
  !===========================================================================================================
  ! Anmerkung: Nur der Vollstaendigkeit halber aufgefuehrt, Druck sollte eigentlich nicht gefiltert werden.
  IF (dir == 0 .AND. m == 0) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 1) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 2) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 3) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     !   - entfallen, da der Druck ohne Dirichlet-RB auskommen sollte
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Konzentrationen =======================================================================================
  !===========================================================================================================
  IF (dir == 0 .AND. m > 0) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 1) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFc1(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFc1(ii,i,m)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 2) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFc2(b2L,j,m)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFc2(jj,j,m)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 3) THEN
        
        DO k = S3p, N3p
           DO j = S2p, N2p
              DO i = S1p, N1p
                 fil(i,j,k) = cFc3(b3L,k,m)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFc3(kk,k,m)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     ! TEST!!! (siehe Geschwindigkeiten ...)
     IF (filter_BC_yes) THEN ! TEST!!! Name irreführend!
     IF (BCc_1L(m) == 1) fil(1 ,S2p:N2p,S3p:N3p) = phi(1 ,S2p:N2p,S3p:N3p)
     IF (BCc_1U(m) == 1) fil(N1,S2p:N2p,S3p:N3p) = phi(N1,S2p:N2p,S3p:N3p)
     
     IF (BCc_2L(m) == 1) fil(S1p:N1p,1 ,S3p:N3p) = phi(S1p:N1p,1 ,S3p:N3p)
     IF (BCc_2U(m) == 1) fil(S1p:N1p,N2,S3p:N3p) = phi(S1p:N1p,N2,S3p:N3p)
     
     IF (BCc_3L(m) == 1) fil(S1p:N1p,S2p:N2p,1 ) = phi(S1p:N1p,S2p:N2p,1 )
     IF (BCc_3U(m) == 1) fil(S1p:N1p,S2p:N2p,N3) = phi(S1p:N1p,S2p:N2p,N3)
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Geschwindigkeiten =====================================================================================
  !===========================================================================================================
  IF (dir == 1) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 1) THEN
        
        DO k = S31B, N31B ! TEST!!! Grenzen sind nun zu eng, koennen gegen S31, N31, etc. ausgetauscht werden ...
           DO j = S21B, N21B
              DO i = S11B, N11B
                 fil(i,j,k) = cFu1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFu1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 2) THEN
        
        DO k = S31B, N31B
           DO j = S21B, N21B
              DO i = S11B, N11B
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 3) THEN
        
        DO k = S31B, N31B
           DO j = S21B, N21B
              DO i = S11B, N11B
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     ! TEST!!! BC_XX == 2,3 fehlen noch ...
     ! TEST!!! phi durch bc12, bc13 ersetzen?
     ! TEST!!! (Hat sich nicht bewährt. Insbesondere die Extrpolation macht Probleme ...)
     IF (filter_BC_yes) THEN
     IF (BC_2L == 1) fil(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     IF (BC_2U == 1) fil(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     
     IF (BC_3L == 1) fil(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     IF (BC_3U == 1) fil(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     
     IF (BC_1L > 0) THEN ! TEST!!! Funktionert NICHT bei Flussmündung und funktioniert bei Lock-Exchange ...
        fil(0 ,S21B:N21B,S31B:N31B) = bc11(S21B:N21B,S31B:N31B,1)
        CALL bc_extrapolation(dir,fil)
     END IF
     IF (BC_1U > 0) THEN
        fil(N1,S21B:N21B,S31B:N31B) = bc11(S21B:N21B,S31B:N31B,2)
        CALL bc_extrapolation(dir,fil)
     END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 2) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 1) THEN
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 2) THEN
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 fil(i,j,k) = cFv2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFv2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 3) THEN
        
        DO k = S32B, N32B
           DO j = S22B, N22B
              DO i = S12B, N12B
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (filter_BC_yes) THEN
     IF (BC_1L == 1) fil(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     IF (BC_1U == 1) fil(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     
     IF (BC_3L == 1) fil(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     IF (BC_3U == 1) fil(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     
     IF (BC_2L > 0) THEN
        fil(S12B:N12B,0 ,S32B:N32B) = bc22(S12B:N12B,S32B:N32B,1)
        CALL bc_extrapolation(dir,fil)
     END IF
     IF (BC_2U > 0) THEN
        fil(S12B:N12B,N2,S32B:N32B) = bc22(S12B:N12B,S32B:N32B,2)
        CALL bc_extrapolation(dir,fil)
     END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  IF (dir == 3) THEN
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 1) THEN
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 DO ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 2) THEN
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 DO jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (fil_dir == 3) THEN
        
        DO k = S33B, N33B
           DO j = S23B, N23B
              DO i = S13B, N13B
                 fil(i,j,k) = cFw3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 DO kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFw3(kk,k)*phi(i,j,k+kk)
                 END DO
              END DO
           END DO
        END DO
        
     END IF
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     IF (filter_BC_yes) THEN
     IF (BC_1L == 1) fil(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     IF (BC_1U == 1) fil(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     
     IF (BC_2L == 1) fil(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     IF (BC_2U == 1) fil(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     
     IF (BC_3L > 0) THEN
        fil(S13B:N13B,S23B:N23B,0 ) = bc33(S13B:N13B,S23B:N23B,1)
        CALL bc_extrapolation(dir,fil)
     END IF
     IF (BC_3U > 0) THEN
        fil(S13B:N13B,S23B:N23B,N3) = bc33(S13B:N13B,S23B:N23B,2)
        CALL bc_extrapolation(dir,fil)
     END IF
     END IF
     !--------------------------------------------------------------------------------------------------------
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE filter
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE bc_extrapolation(m,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 1) THEN
     IF (BC_1L > 0) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = 0, d1U
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,1)*phi(1+ii,j,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           END DO
        END DO
     END IF
     IF (BC_1U > 0) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = d1L, -1
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,i)*phi(i+ii,j,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (BC_2L > 0) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = 0, d2U
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,1)*phi(i,1+jj,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_2U > 0) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = d2L, -1
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,j)*phi(i,j+jj,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (BC_3L > 0) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = 0, d3U
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,1)*phi(i,j,1+kk)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_3U > 0) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = d3L, -1
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,k)*phi(i,j,k+kk)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE bc_extrapolation
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE bc_extrapolation_transp(m,phi)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  m
  REAL   , INTENT(INOUT) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 1) THEN
     IF (BC_1L > 0) THEN
        i = 0
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = 0, d1U
                 phi(i+ii+1,j,k) = phi(i+ii+1,j,k) - phi(i,j,k)*cIup(ii,1)/cIup(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           END DO
        END DO
     END IF
     IF (BC_1U > 0) THEN
        i = N1
        DO k = S31B, N31B
           DO j = S21B, N21B
!pgi$ unroll = n:8
              DO ii = d1L, -1
                 phi(i+ii,j,k) = phi(i+ii,j,k) - phi(i,j,k)*cIup(ii,i)/cIup(0,i)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 2) THEN
     IF (BC_2L > 0) THEN
        j = 0
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = 0, d2U
                 phi(i,j+jj+1,k) = phi(i,j+jj+1,k) - phi(i,j,k)*cIvp(jj,1)/cIvp(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_2U > 0) THEN
        j = N2
        DO k = S32B, N32B
           DO i = S12B, N12B
!pgi$ unroll = n:8
              DO jj = d2L, -1
                 phi(i,j+jj,k) = phi(i,j+jj,k) - phi(i,j,k)*cIvp(jj,j)/cIvp(0,j)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (m == 3) THEN
     IF (BC_3L > 0) THEN
        k = 0
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = 0, d3U
                 phi(i,j,k+kk+1) = phi(i,j,k+kk+1) - phi(i,j,k)*cIwp(kk,1)/cIwp(-1,1)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           END DO
        END DO
     END IF
     IF (BC_3U > 0) THEN
        k = N3
        DO j = S23B, N23B
           DO i = S13B, N13B
!pgi$ unroll = n:8
              DO kk = d3L, -1
                 phi(i,j,k+kk) = phi(i,j,k+kk) - phi(i,j,k)*cIwp(kk,k)/cIwp(0,k)
              END DO
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           END DO
        END DO
     END IF
  END IF
  !-----------------------------------------------------------------------------------------------------------
  
  
  END SUBROUTINE bc_extrapolation_transp
  
  
  
END MODULE mod_diff
