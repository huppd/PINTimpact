!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing many routines to compute derivatives ...
module cmod_operator
  
  
  use mod_dims
  use mod_vars
  use mod_exchange

  use iso_c_binding
  
  
  private
  
  
  public apply_compact, apply_compact_transp
  public divergence, divergence2, divergence_transp
  public gradient, gradient_transp
  public Helmholtz, Helmholtz_explicit, Helmholtz_conc, Helmholtz_conc_explicit
  public nonlinear, nonlinear_conc
  public interpolate_vel, interpolate_conc
  public outflow_bc, sediment_bc
  public filter
  public bc_extrapolation, bc_extrapolation_transp
  public interpolate_pre_vel , interpolate_vel_pre
  public interpolate2_pre_vel, interpolate2_vel_pre
  public first_pre_vel, first_vel_pre
  public first_adv_pre, first_adv_vel
  public Helmholtz_pre_explicit ! TEST!!!
  
  
  INCLUDE 'mpif.h'
  
  contains
  
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
  subroutine apply_compact(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,Nmax,ndL,ndR,dimS,ccL,ccL_LU,ccR,WW,Schur,phi,der)
  
  implicit none
  
  integer, intent(in   ) ::  dir, vel_dir
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  integer, intent(in   ) ::  Nmax, ndL, ndR, dimS
  real   , intent(in   ) ::  ccL (-ndL:ndL,0:Nmax)
  real   , intent(in   ) ::  ccR (-ndR:ndR,0:Nmax)
  real   , intent(in   ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  real   , intent(in   ) ::  WW(1:2*ndL,0:Nmax)
  real   , intent(in   ) ::  Schur(1:dimS,1:dimS)
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii, i0
  integer                ::  j, jj, j0
  integer                ::  k, kk, k0
  
  integer                ::  dimG
  integer                ::  KD, KO, LM
  
  
  ! Wurde herausgezogen:
  !!CALL exchange (dir,vel_dir,phi)
  !CALL exchange2(dir,vel_dir,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  call pseudocall_int(vel_dir) ! TEST!!! this action has no meaning - it shall only to suppress the warnings
                               ! that the argument vel_dir is not used when building the executable ...
  
  KD = ndL+ndL+1
  KO = ndL+ndL
  
  !===========================================================================================================
  if (dir == 1) then
     !--------------------------------------------------------------------------------------------------------
     ! TEST!!! Block-Raender werden nun IMMER per Schur-Komplement behandelt (Pivoting-Problem) ...
     !         Man koennte/sollte aber immerhin das Umspeichern von buffer1A auf buffer1B vermeiden ...
     if (1 == 2 .and. NB1 == 1 .and. BC_1L_global /= -1) then
        
        do k = SS3, NN3
           do j = SS2, NN2
              
              !--- expliziter Teil ---
              do i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
                 do ii = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
                 end do
              end do
              
              !--- L^-1*der ---
              do i = SS1, NN1-1
                 LM = MIN(ndL,NN1-i)
                 do ii = 1, LM
                    der(i+ii,j,k) = der(i+ii,j,k) - ccL_LU(KD+ii,i)*der(i,j,k)
                 end do
              end do
              
              !--- U^-1*L^-1*der ---
              do i = NN1, SS1, -1
                 LM = MAX(SS1,i-KO)
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,i)
                 do ii = i-1, LM, -1
                    der(ii,j,k) = der(ii,j,k) - ccL_LU(KD+ii-i,i)*der(i,j,k)
                 end do
              end do
              
           end do
        end do
        
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NN1-SS1+1-2*ndL
        i0   = 2*ndL*(iB(1,1)-1)
        
        do k = SS3, NN3
           do j = SS2, NN2
              
              !--- expliziter Teil ---
              do i = SS1, NN1 ! Anmerkung: Man koennte hier nur den für buffer1A notwendigen Teil von "der" berechnen und den anderen erst nach MPI_ALLGATHERv ...
                 der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
                 do ii = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
                 end do
              end do
              
              !--- Umspeichern ---
              do i = 1, ndL
                 buffer1A(j,k,i    ) = der(SS1+i-1  ,j,k)
                 buffer1A(j,k,i+ndL) = der(NN1+i-ndL,j,k)
              end do
              
              !--- RHS fuer Schur-Komplement ---
              do i = 1, dimG
                 do ii = 1, 2*ndL
                    buffer1A(j,k,ii) = buffer1A(j,k,ii) - WW(ii,SS1+i+ndL-1)*der(SS1+i+ndL-1,j,k)
                 end do
              end do
              
           end do
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        ! TEST!!! Evtl. ist das nicht die effizenteste Methode!!
        call MPI_ALLGATHERv(buffer1A,2*ndl*(N2+1)*(N3+1),MPI_REAL8,buffer1B,recv1,disp1,MPI_REAL8,COMM_BAR1,merror)
        
        do k = SS3, NN3
           do j = SS2, NN2
              
              !--- Loesung des Schur-Komplement-LGS ---
              do i = 1, ndL
                 der(SS1+i-1  ,j,k) = Schur(i+i0    ,1)*buffer1B(j,k,1)
                 der(NN1+i-ndL,j,k) = Schur(i+i0+ndL,1)*buffer1B(j,k,1)
                 do ii = 2, dimS
                    der(SS1+i-1  ,j,k) = der(SS1+i-1  ,j,k) + Schur(i+i0    ,ii)*buffer1B(j,k,ii)
                    der(NN1+i-ndL,j,k) = der(NN1+i-ndL,j,k) + Schur(i+i0+ndL,ii)*buffer1B(j,k,ii)
                 end do
              end do
              
              !--- RHS des uebrigen Feldes anpassen ---
              do i = 1, ndL
                 do ii = i, ndL
                    der(SS1+i-1+ndL,j,k) = der(SS1+i-1+ndL,j,k) - ccL(-ii,SS1+i-1+ndL)*der(SS1+i-1+ndL-ii,j,k)
                    der(NN1-i+1-ndL,j,k) = der(NN1-i+1-ndL,j,k) - ccL( ii,NN1-i+1-ndL)*der(NN1-i+1-ndL+ii,j,k)
                 end do
              end do
              
              !--- L^-1*der ---
              do i = SS1+ndL, NN1-ndL-1
                 LM = MIN(ndL,NN1-ndL-i)
                 do ii = 1, LM
                    der(i+ii,j,k) = der(i+ii,j,k) - ccL_LU(KD+ii,i)*der(i,j,k)
                 end do
              end do
              
              !--- U^-1*L^-1*der ---
              do i = NN1-ndL, SS1+ndL, -1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,i)
                 LM = MAX(SS1+ndL,i-KO)
                 do ii = i-1, LM, -1
                    der(ii,j,k) = der(ii,j,k) - ccL_LU(KD+ii-i,i)*der(i,j,k)
                 end do
              end do
              
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (1 == 2 .and. NB2 == 1 .and. BC_2L_global /= -1) then
        
        do k = SS3, NN3
           
           !--- expliziter Teil ---
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
                 do jj = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
           
           !--- L^-1*der ---
           do j = SS2, NN2-1
              LM = MIN(ndL,NN2-j)
              do i = SS1, NN1
                 do jj = 1, LM
                    der(i,j+jj,k) = der(i,j+jj,k) - ccL_LU(KD+jj,j)*der(i,j,k)
                 end do
              end do
           end do
           
           !--- U^-1*L^-1*der ---
           do j = NN2, SS2, -1
              LM = MAX(SS2,j-KO)
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,j)
                 do jj = j-1, LM, -1
                    der(i,jj,k) = der(i,jj,k) - ccL_LU(KD+jj-j,j)*der(i,j,k)
                 end do
              end do
           end do
           
        end do
        
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NN2-SS2+1-2*ndL
        j0   = 2*ndL*(iB(2,1)-1)
        
        do k = SS3, NN3
           
           !--- expliziter Teil ---
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
                 do jj = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
           
           !--- Umspeichern ---
           do j = 1, ndL
              do i = SS1, NN1
                 buffer2A(i,k,j    ) = der(i,SS2+j-1  ,k)
                 buffer2A(i,k,j+ndL) = der(i,NN2+j-ndL,k)
              end do
           end do
           
           !--- RHS fuer Schur-Komplement ---
           do j = 1, dimG
              do i = SS1, NN1
                 do jj = 1, 2*ndL
                    buffer2A(i,k,jj) = buffer2A(i,k,jj) - WW(jj,SS2+j+ndL-1)*der(i,SS2+j+ndL-1,k)
                 end do
              end do
           end do
           
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        call MPI_ALLGATHERv(buffer2A,2*ndl*(N1+1)*(N3+1),MPI_REAL8,buffer2B,recv2,disp2,MPI_REAL8,COMM_BAR2,merror)
        
        do k = SS3, NN3
           
           !--- Loesung des Schur-Komplement-LGS ---
           do j = 1, ndL
              do i = SS1, NN1
                 der(i,SS2+j-1  ,k) = Schur(j+j0    ,1)*buffer2B(i,k,1)
                 der(i,NN2+j-ndL,k) = Schur(j+j0+ndL,1)*buffer2B(i,k,1)
                 do jj = 2, dimS
                    der(i,SS2+j-1  ,k) = der(i,SS2+j-1  ,k) + Schur(j+j0    ,jj)*buffer2B(i,k,jj)
                    der(i,NN2+j-ndL,k) = der(i,NN2+j-ndL,k) + Schur(j+j0+ndL,jj)*buffer2B(i,k,jj)
                 end do
              end do
           end do
           
           !--- RHS des uebrigen Feldes anpassen ---
           do j = 1, ndL
              do i = SS1, NN1
                 do jj = j, ndL
                    der(i,SS2+j-1+ndL,k) = der(i,SS2+j-1+ndL,k) - ccL(-jj,SS2+j-1+ndL)*der(i,SS2+j-1+ndL-jj,k)
                    der(i,NN2-j+1-ndL,k) = der(i,NN2-j+1-ndL,k) - ccL( jj,NN2-j+1-ndL)*der(i,NN2-j+1-ndL+jj,k)
                 end do
              end do
           end do
           
           !--- L^-1*der ---
           do j = SS2+ndL, NN2-ndL-1
              LM = MIN(ndL,NN2-ndL-j)
              do i = SS1, NN1
                 do jj = 1, LM
                    der(i,j+jj,k) = der(i,j+jj,k) - ccL_LU(KD+jj,j)*der(i,j,k)
                 end do
              end do
           end do
           
           !--- U^-1*L^-1*der ---
           do j = NN2-ndL, SS2+ndL, -1
              LM = MAX(SS2+ndL,j-KO)
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,j)
                 do jj = j-1, LM, -1
                    der(i,jj,k) = der(i,jj,k) - ccL_LU(KD+jj-j,j)*der(i,j,k)
                 end do
              end do
           end do
           
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (1 == 2 .and. NB3 == 1 .and. BC_3L_global /= -1) then
        
        !--- expliziter Teil ---
        do k = SS3, NN3
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
                 do kk = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
        !--- L^-1*der ---
        do k = SS3, NN3-1
           LM = MIN(ndL,NN3-k)
           do j = SS2, NN2
              do i = SS1, NN1
                 do kk = 1, LM
                    der(i,j,k+kk) = der(i,j,k+kk) - ccL_LU(KD+kk,k)*der(i,j,k)
                 end do
              end do
           end do
        end do
        
        !--- U^-1*L^-1*der ---
        do k = NN3, SS3, -1
           LM = MAX(SS3,k-KO)
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,k)
                 do kk = k-1, LM, -1
                    der(i,j,kk) = der(i,j,kk) - ccL_LU(KD+kk-k,k)*der(i,j,k)
                 end do
              end do
           end do
        end do
        
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NN3-SS3+1-2*ndL
        k0   = 2*ndL*(iB(3,1)-1)
        
        !--- expliziter Teil ---
        do k = SS3, NN3
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
                 do kk = -ndR+1, ndR
                    der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
        !--- Umspeichern ---
        do k = 1, ndL
           do j = SS2, NN2
              do i = SS1, NN1
                 buffer3A(i,j,k    ) = der(i,j,SS3+k-1  )
                 buffer3A(i,j,k+ndL) = der(i,j,NN3+k-ndL)
              end do
           end do
        end do
        
        !--- RHS fuer Schur-Komplement ---
        do k = 1, dimG
           do j = SS2, NN2
              do i = SS1, NN1
                 do kk = 1, 2*ndL
                    buffer3A(i,j,kk) = buffer3A(i,j,kk) - WW(kk,SS3+k+ndL-1)*der(i,j,SS3+k+ndL-1)
                 end do
              end do
           end do
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        call MPI_ALLGATHERv(buffer3A,2*ndl*(N1+1)*(N2+1),MPI_REAL8,buffer3B,recv3,disp3,MPI_REAL8,COMM_BAR3,merror)
        
        !--- Loesung des Schur-Komplement-LGS ---
        do k = 1, ndL
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,SS3+k-1  ) = Schur(k+k0    ,1)*buffer3B(i,j,1)
                 der(i,j,NN3+k-ndL) = Schur(k+k0+ndL,1)*buffer3B(i,j,1)
                 do kk = 2, dimS
                    der(i,j,SS3+k-1  ) = der(i,j,SS3+k-1  ) + Schur(k+k0    ,kk)*buffer3B(i,j,kk)
                    der(i,j,NN3+k-ndL) = der(i,j,NN3+k-ndL) + Schur(k+k0+ndL,kk)*buffer3B(i,j,kk)
                 end do
              end do
           end do
        end do
        
        !--- RHS des uebrigen Feldes anpassen ---
        do k = 1, ndL
           do j = SS2, NN2
              do i = SS1, NN1
                 do kk = k, ndL
                    der(i,j,SS3+k-1+ndL) = der(i,j,SS3+k-1+ndL) - ccL(-kk,SS3+k-1+ndL)*der(i,j,SS3+k-1+ndL-kk)
                    der(i,j,NN3-k+1-ndL) = der(i,j,NN3-k+1-ndL) - ccL( kk,NN3-k+1-ndL)*der(i,j,NN3-k+1-ndL+kk)
                 end do
              end do
           end do
        end do
        
        !--- L^-1*der ---
        do k = SS3+ndL, NN3-ndL-1
           LM = MIN(ndL,NN3-ndL-k)
           do j = SS2, NN2
              do i = SS1, NN1
                 do kk = 1, LM
                    der(i,j,k+kk) = der(i,j,k+kk) - ccL_LU(KD+kk,k)*der(i,j,k)
                 end do
              end do
           end do
        end do
        
        !--- U^-1*L^-1*der ---
        do k = NN3-ndL, SS3+ndL, -1
           LM = MAX(SS3+ndL,k-KO)
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k) / ccL_LU(KD,k)
                 do kk = k-1, LM, -1
                    der(i,j,kk) = der(i,j,kk) - ccL_LU(KD+kk-k,k)*der(i,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine apply_compact
  
  
  
  
  
  
  
  
  
  
  ! ACTHUNG: phi wird überschrieben (da ohnehin mit dx vormultipliziert werden muss)
  subroutine apply_compact_transp(dir,vel_dir,SI1,SI2,SI3,NI1,NI2,NI3,SE1,SE2,SE3,NE1,NE2,NE3,Nmax,ndL,ndR,dimS,ccL,ccL_LU,ccR,WW,Schur,phi,der)
  
  implicit none
  
  integer, intent(in   ) ::  dir, vel_dir
  integer, intent(in   ) ::  SI1, SI2, SI3, NI1, NI2, NI3
  integer, intent(in   ) ::  SE1, SE2, SE3, NE1, NE2, NE3
  integer, intent(in   ) ::  Nmax, ndL, ndR, dimS
  real   , intent(in   ) ::  ccL (-ndL:ndL,0:Nmax)
  real   , intent(in   ) ::  ccR (-ndR:ndR,0:Nmax)
  real   , intent(in   ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  real   , intent(in   ) ::  WW(1:2*ndL,0:Nmax)
  real   , intent(in   ) ::  Schur(1:dimS,1:dimS)
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii, i0
  integer                ::  j, jj, j0
  integer                ::  k, kk, k0
  
  integer                ::  dimG
  integer                ::  KD, KO, LM
  
  
  KD = ndL+ndL+1
  KO = ndL+ndL
  
  !===========================================================================================================
  if (dir == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (NB1 == 1 .and. BC_1L_global /= -1) then
        
        do k = SI3, NI3
           do j = SI2, NI2
              
              !--- L^-1*der ---
              do i = SI1, NI1-1
                 LM = MIN(ndL,NI1-i)
                 do ii = 1, LM
                    phi(i+ii,j,k) = phi(i+ii,j,k) - ccL_LU(KD+ii,i)*phi(i,j,k)
                 end do
              end do
              
              !--- U^-1*L^-1*der ---
              do i = NI1, SI1, -1
                 LM = MAX(SI1,i-KO)
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,i)
                 do ii = i-1, LM, -1
                    phi(ii,j,k) = phi(ii,j,k) - ccL_LU(KD+ii-i,i)*phi(i,j,k)
                 end do
              end do
              
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NI1-SI1+1-2*ndL
        i0   = 2*ndL*(iB(1,1)-1)
        
        do k = SI3, NI3
           do j = SI2, NI2
              
              !--- Umspeichern ---
              do i = 1, ndL
                 buffer1A(j,k,i    ) = phi(SI1+i-1  ,j,k)
                 buffer1A(j,k,i+ndL) = phi(NI1+i-ndL,j,k)
              end do
              
              !--- RHS fuer Schur-Komplement ---
              do i = 1, dimG
                 do ii = 1, 2*ndL
                    buffer1A(j,k,ii) = buffer1A(j,k,ii) - WW(ii,SI1+i+ndL-1)*phi(SI1+i+ndL-1,j,k)
                 end do
              end do
              
           end do
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        call MPI_ALLGATHERv(buffer1A,2*ndl*(N2+1)*(N3+1),MPI_REAL8,buffer1B,recv1,disp1,MPI_REAL8,COMM_BAR1,merror)
        
        do k = SI3, NI3
           do j = SI2, NI2
              
              !--- Loesung des Schur-Komplement-LGS ---
              do i = 1, ndL
                 phi(SI1+i-1  ,j,k) = Schur(i+i0    ,1)*buffer1B(j,k,1)
                 phi(NI1+i-ndL,j,k) = Schur(i+i0+ndL,1)*buffer1B(j,k,1)
                 do ii = 2, dimS
                    phi(SI1+i-1  ,j,k) = phi(SI1+i-1  ,j,k) + Schur(i+i0    ,ii)*buffer1B(j,k,ii)
                    phi(NI1+i-ndL,j,k) = phi(NI1+i-ndL,j,k) + Schur(i+i0+ndL,ii)*buffer1B(j,k,ii)
                 end do
              end do
              
              !--- RHS des uebrigen Feldes anpassen ---
              do i = 1, ndL
                 do ii = i, ndL
                    phi(SI1+i-1+ndL,j,k) = phi(SI1+i-1+ndL,j,k) - ccL(-ii,SI1+i-1+ndL)*phi(SI1+i-1+ndL-ii,j,k)
                    phi(NI1-i+1-ndL,j,k) = phi(NI1-i+1-ndL,j,k) - ccL( ii,NI1-i+1-ndL)*phi(NI1-i+1-ndL+ii,j,k)
                 end do
              end do
              
              !--- L^-1*der ---
              do i = SI1+ndL, NI1-ndL-1
                 LM = MIN(ndL,NI1-ndL-i)
                 do ii = 1, LM
                    phi(i+ii,j,k) = phi(i+ii,j,k) - ccL_LU(KD+ii,i)*phi(i,j,k)
                 end do
              end do
              
              !--- U^-1*L^-1*der ---
              do i = NI1-ndL, SI1+ndL, -1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,i)
                 LM = MAX(SI1+ndL,i-KO)
                 do ii = i-1, LM, -1
                    phi(ii,j,k) = phi(ii,j,k) - ccL_LU(KD+ii-i,i)*phi(i,j,k)
                 end do
              end do
              
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     call exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     do k = SE3, NE3
        do j = SE2, NE2
           do i = SE1, NE1
              der(i,j,k) = ccR(-ndR,i)*phi(i-ndR,j,k)
              do ii = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (NB2 == 1 .and. BC_2L_global /= -1) then
        
        do k = SI3, NI3
           
           !--- L^-1*der ---
           do j = SI2, NI2-1
              LM = MIN(ndL,NI2-j)
              do i = SI1, NI1
                 do jj = 1, LM
                    phi(i,j+jj,k) = phi(i,j+jj,k) - ccL_LU(KD+jj,j)*phi(i,j,k)
                 end do
              end do
           end do
           
           !--- U^-1*L^-1*der ---
           do j = NI2, SI2, -1
              LM = MAX(SI2,j-KO)
              do i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,j)
                 do jj = j-1, LM, -1
                    phi(i,jj,k) = phi(i,jj,k) - ccL_LU(KD+jj-j,j)*phi(i,j,k)
                 end do
              end do
           end do
           
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NI2-SI2+1-2*ndL
        j0   = 2*ndL*(iB(2,1)-1)
        
        do k = SI3, NI3
           
           !--- Umspeichern ---
           do j = 1, ndL
              do i = SI1, NI1
                 buffer2A(i,k,j    ) = phi(i,SI2+j-1  ,k)
                 buffer2A(i,k,j+ndL) = phi(i,NI2+j-ndL,k)
              end do
           end do
           
           !--- RHS fuer Schur-Komplement ---
           do j = 1, dimG
              do i = SI1, NI1
                 do jj = 1, 2*ndL
                    buffer2A(i,k,jj) = buffer2A(i,k,jj) - WW(jj,SI2+j+ndL-1)*phi(i,SI2+j+ndL-1,k)
                 end do
              end do
           end do
           
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        call MPI_ALLGATHERv(buffer2A,2*ndl*(N1+1)*(N3+1),MPI_REAL8,buffer2B,recv2,disp2,MPI_REAL8,COMM_BAR2,merror)
        
        do k = SI3, NI3
           
           !--- Loesung des Schur-Komplement-LGS ---
           do j = 1, ndL
              do i = SI1, NI1
                 phi(i,SI2+j-1  ,k) = Schur(j+j0    ,1)*buffer2B(i,k,1)
                 phi(i,NI2+j-ndL,k) = Schur(j+j0+ndL,1)*buffer2B(i,k,1)
                 do jj = 2, dimS
                    phi(i,SI2+j-1  ,k) = phi(i,SI2+j-1  ,k) + Schur(j+j0    ,jj)*buffer2B(i,k,jj)
                    phi(i,NI2+j-ndL,k) = phi(i,NI2+j-ndL,k) + Schur(j+j0+ndL,jj)*buffer2B(i,k,jj)
                 end do
              end do
           end do
           
           !--- RHS des uebrigen Feldes anpassen ---
           do j = 1, ndL
              do i = SI1, NI1
                 do jj = j, ndL
                    phi(i,SI2+j-1+ndL,k) = phi(i,SI2+j-1+ndL,k) - ccL(-jj,SI2+j-1+ndL)*phi(i,SI2+j-1+ndL-jj,k)
                    phi(i,NI2-j+1-ndL,k) = phi(i,NI2-j+1-ndL,k) - ccL( jj,NI2-j+1-ndL)*phi(i,NI2-j+1-ndL+jj,k)
                 end do
              end do
           end do
           
           !--- L^-1*der ---
           do j = SI2+ndL, NI2-ndL-1
              LM = MIN(ndL,NI2-ndL-j)
              do i = SI1, NI1
                 do jj = 1, LM
                    phi(i,j+jj,k) = phi(i,j+jj,k) - ccL_LU(KD+jj,j)*phi(i,j,k)
                 end do
              end do
           end do
           
           !--- U^-1*L^-1*der ---
           do j = NI2-ndL, SI2+ndL, -1
              LM = MAX(SI2+ndL,j-KO)
              do i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,j)
                 do jj = j-1, LM, -1
                    phi(i,jj,k) = phi(i,jj,k) - ccL_LU(KD+jj-j,j)*phi(i,j,k)
                 end do
              end do
           end do
           
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     call exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     do k = SE3, NE3
        do j = SE2, NE2
           do i = SE1, NE1
              der(i,j,k) = ccR(-ndR,j)*phi(i,j-ndR,k)
              do jj = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (NB3 == 1 .and. BC_3L_global /= -1) then
        
        !--- L^-1*der ---
        do k = SI3, NI3-1
           LM = MIN(ndL,NI3-k)
           do j = SI2, NI2
              do i = SI1, NI1
                 do kk = 1, LM
                    phi(i,j,k+kk) = phi(i,j,k+kk) - ccL_LU(KD+kk,k)*phi(i,j,k)
                 end do
              end do
           end do
        end do
        
        !--- U^-1*L^-1*der ---
        do k = NI3, SI3, -1
           LM = MAX(SI3,k-KO)
           do j = SI2, NI2
              do i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,k)
                 do kk = k-1, LM, -1
                    phi(i,j,kk) = phi(i,j,kk) - ccL_LU(KD+kk-k,k)*phi(i,j,k)
                 end do
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        
        dimG = NI3-SI3+1-2*ndL
        k0   = 2*ndL*(iB(3,1)-1)
        
        !--- Umspeichern ---
        do k = 1, ndL
           do j = SI2, NI2
              do i = SI1, NI1
                 buffer3A(i,j,k    ) = phi(i,j,SI3+k-1  )
                 buffer3A(i,j,k+ndL) = phi(i,j,NI3+k-ndL)
              end do
           end do
        end do
        
        !--- RHS fuer Schur-Komplement ---
        do k = 1, dimG
           do j = SI2, NI2
              do i = SI1, NI1
                 do kk = 1, 2*ndL
                    buffer3A(i,j,kk) = buffer3A(i,j,kk) - WW(kk,SI3+k+ndL-1)*phi(i,j,SI3+k+ndL-1)
                 end do
              end do
           end do
        end do
        
        !--- Verteilen der RHS fuer Schur-Problem ---
        call MPI_ALLGATHERv(buffer3A,2*ndl*(N1+1)*(N2+1),MPI_REAL8,buffer3B,recv3,disp3,MPI_REAL8,COMM_BAR3,merror)
        
        !--- Loesung des Schur-Komplement-LGS ---
        do k = 1, ndL
           do j = SI2, NI2
              do i = SI1, NI1
                 phi(i,j,SI3+k-1  ) = Schur(k+k0    ,1)*buffer3B(i,j,1)
                 phi(i,j,NI3+k-ndL) = Schur(k+k0+ndL,1)*buffer3B(i,j,1)
                 do kk = 2, dimS
                    phi(i,j,SI3+k-1  ) = phi(i,j,SI3+k-1  ) + Schur(k+k0    ,kk)*buffer3B(i,j,kk)
                    phi(i,j,NI3+k-ndL) = phi(i,j,NI3+k-ndL) + Schur(k+k0+ndL,kk)*buffer3B(i,j,kk)
                 end do
              end do
           end do
        end do
        
        !--- RHS des uebrigen Feldes anpassen ---
        do k = 1, ndL
           do j = SI2, NI2
              do i = SI1, NI1
                 do kk = k, ndL
                    phi(i,j,SI3+k-1+ndL) = phi(i,j,SI3+k-1+ndL) - ccL(-kk,SI3+k-1+ndL)*phi(i,j,SI3+k-1+ndL-kk)
                    phi(i,j,NI3-k+1-ndL) = phi(i,j,NI3-k+1-ndL) - ccL( kk,NI3-k+1-ndL)*phi(i,j,NI3-k+1-ndL+kk)
                 end do
              end do
           end do
        end do
        
        !--- L^-1*der ---
        do k = SI3+ndL, NI3-ndL-1
           LM = MIN(ndL,NI3-ndL-k)
           do j = SI2, NI2
              do i = SI1, NI1
                 do kk = 1, LM
                    phi(i,j,k+kk) = phi(i,j,k+kk) - ccL_LU(KD+kk,k)*phi(i,j,k)
                 end do
              end do
           end do
        end do
        
        !--- U^-1*L^-1*der ---
        do k = NI3-ndL, SI3+ndL, -1
           LM = MAX(SI3+ndL,k-KO)
           do j = SI2, NI2
              do i = SI1, NI1
                 phi(i,j,k) = phi(i,j,k) / ccL_LU(KD,k)
                 do kk = k-1, LM, -1
                    phi(i,j,kk) = phi(i,j,kk) - ccL_LU(KD+kk-k,k)*phi(i,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     call exchange2(dir,vel_dir,SE1,SE2,SE3,NE1,NE2,NE3,phi)
     
     !--- expliziter Teil ---
     do k = SE3, NE3
        do j = SE2, NE2
           do i = SE1, NE1
              der(i,j,k) = ccR(-ndR,k)*phi(i,j,k-ndR)
              do kk = -ndR+1, ndR
                 der(i,j,k) = der(i,j,k) + ccR(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine apply_compact_transp
  
  
  
  
  
  
  
  
  
  
  
  subroutine divergence(m,phi,div)
  
  implicit none
  
  integer, intent(in)    ::  m
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  call exchange(m,m,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phi,div)
        
        do k = S3p, N3p
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 div(i,j,k) = div(i,j,k)*dx1DM(i)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 do ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phi,com)
        
        do k = S3p, N3p
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 div(i,j,k) = div(i,j,k) + com(i,j,k)*dx2DM(j)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
!pgi$ unroll = n:8
                 do jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phi,com)
        
        do k = S3p, N3p
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 div(i,j,k) = div(i,j,k) + com(i,j,k)*dx3DM(k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
!pgi$ unroll = n:8
                 do kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine divergence
  
  
  
  
  
  !> \brief computes divergence
  !! same as Impact \c mod_diff::divergence
  !! side effct phiU, phiV and phiW
  !! \param[inout] phiU
  !! \param[inout] phiV
  !! \param[inout] phiW
  !! \param[out] div
  subroutine divergence2(phiU,phiV,phiW,div) bind (c,name='OP_div')
  
  implicit none
  
  real(c_double), intent(inout) ::  phiU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  phiW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(  out) ::  div (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  call exchange(1,1,phiU(b1L,b2L,b3L))
  call exchange(2,2,phiV(b1L,b2L,b3L))
  call exchange(3,3,phiW(b1L,b2L,b3L))
  
  
  !===========================================================================================================
  if (comp_div_yes) then
     !--------------------------------------------------------------------------------------------------------
     call apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phiU(b1L,b2L,b3L),div)
     
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              div(i,j,k) = div(i,j,k)*dx1DM(i)
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
     call apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phiV(b1L,b2L,b3L),com)
     
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              div(i,j,k) = div(i,j,k) + com(i,j,k)*dx2DM(j)
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     
     call apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phiW(b1L,b2L,b3L),com)
     
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              div(i,j,k) = div(i,j,k) + com(i,j,k)*dx3DM(k)
           end do
        end do
     end do
     
     end if
     !--------------------------------------------------------------------------------------------------------
     
  !===========================================================================================================
  else
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phiU(i+d1L,j,k)
!pgi$ unroll = n:8
                 do ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phiU(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phiV(i,j+jj,k)
                 end do
!pgi$ unroll = n:8
                 do kk = d3L, d3U
                    div(i,j,k) = div(i,j,k) + cDw3(kk,k)*phiW(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     !--------------------------------------------------------------------------------------------------------
     else
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 div(i,j,k) = cDu1(d1L,i)*phiU(i+d1L,j,k)
!pgi$ unroll = n:8
                 do ii = d1L+1, d1U
                    div(i,j,k) = div(i,j,k) + cDu1(ii,i)*phiU(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = d2L, d2U
                    div(i,j,k) = div(i,j,k) + cDv2(jj,j)*phiV(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine divergence2
  
  
  
  
  
  
  
  
  
  
  
  subroutine divergence_transp(m,phi,div)
  
  implicit none
  
  integer, intent(in   ) ::  m
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  div(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        
        do k = S3p, N3p ! Intervall ist ok!
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx1DM(i)
              end do
           end do
        end do
        
        call apply_compact_transp(1,0,S1p,S2p,S3p,N1p,N2p,N3p,S11B,S21B,S31B,N11B,N21B,N31B,N1,ndL,ndR,dimS1,cDu1CLT,cDu1CLT_LU,cDu1CRT,WDu1T,SDu1T,com,div)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,0,phi)
        
        do k = S31B, N31B ! "B"-Grenzen etwas konsistenter, aber nicht wirklich notwendig, da ohnehin mit bc_extrapolation nachmultipliziert wird ...
           do j = S21B, N21B
              do i = S11B, N11B
                 div(i,j,k) = cDu1T(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 do ii = g1L+1, g1U
                    div(i,j,k) = div(i,j,k) + cDu1T(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        
        do k = S3p, N3p
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx2DM(j)
              end do
           end do
        end do
        
        call apply_compact_transp(2,0,S1p,S2p,S3p,N1p,N2p,N3p,S12B,S22B,S32B,N12B,N22B,N32B,N2,ndL,ndR,dimS2,cDv2CLT,cDv2CLT_LU,cDv2CRT,WDv2T,SDv2T,com,div)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,0,phi)
        
        do k = S32B, N32B
           do j = S22B, N22B
              do i = S12B, N12B
                 div(i,j,k) = cDv2T(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 do jj = g2L+1, g2U
                    div(i,j,k) = div(i,j,k) + cDv2T(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        
        do k = S3p, N3p
           do j = S2p, N2p
!pgi$ unroll = n:8
              do i = S1p, N1p
                 com(i,j,k) = phi(i,j,k)*dx3DM(k)
              end do
           end do
        end do
        
        call apply_compact_transp(3,0,S1p,S2p,S3p,N1p,N2p,N3p,S13B,S23B,S33B,N13B,N23B,N33B,N3,ndL,ndR,dimS3,cDw3CLT,cDw3CLT_LU,cDw3CRT,WDw3T,SDw3T,com,div)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,0,phi)
        
        do k = S33B, N33B
           do j = S23B, N23B
              do i = S13B, N13B
                 div(i,j,k) = cDw3T(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 do kk = g3L+1, g3U
                    div(i,j,k) = div(i,j,k) + cDw3T(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine divergence_transp
  
  
  
  
  
  !> \brief applies gradient
  !! to phi, which is a scalar field, and stores the grad in a vector field, the
  !! boundary condtions for Dirichlet and Neumann are set to zero
  !! side effect \c phi is exchanged
  !! \param[in] m direction in which gradient si calculated
  !! \param[inout] phi ScalarField from which the gradient is taken
  !! \param[out] grad gradient in m direction
  subroutine gradient(m,phi,grad) bind(c,name='OP_grad')
  
  implicit none
  
  integer(c_int), intent(in   ) ::  m
  
  real(c_double), intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  call exchange(m,0,phi) ! irrelefant
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,phi,grad)
        
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 grad(i,j,k) = grad(i,j,k)*dx1GM(i)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 grad(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
!pgi$ unroll = n:8
                 do ii = g1L+1, g1U
                    grad(i,j,k) = grad(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) grad(0 ,S21B:N21B,S31B:N31B) = 0.
     if (BC_1U > 0) grad(N1,S21B:N21B,S31B:N31B) = 0.
     if (BC_2L > 0) grad(S11B:N11B,1 ,S31B:N31B) = 0.
     if (BC_2U > 0) grad(S11B:N11B,N2,S31B:N31B) = 0.
     if (BC_3L > 0) grad(S11B:N11B,S21B:N21B,1 ) = 0.
     if (BC_3U > 0) grad(S11B:N11B,S21B:N21B,N3) = 0.
     
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,phi,grad)
        
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 grad(i,j,k) = grad(i,j,k)*dx2GM(j)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 grad(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
!pgi$ unroll = n:8
                 do jj = g2L+1, g2U
                    grad(i,j,k) = grad(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) grad(1 ,S22B:N22B,S32B:N32B) = 0.
     if (BC_1U > 0) grad(N1,S22B:N22B,S32B:N32B) = 0.
     if (BC_2L > 0) grad(S12B:N12B,0 ,S32B:N32B) = 0.
     if (BC_2U > 0) grad(S12B:N12B,N2,S32B:N32B) = 0.
     if (BC_3L > 0) grad(S12B:N12B,S22B:N22B,1 ) = 0.
     if (BC_3U > 0) grad(S12B:N12B,S22B:N22B,N3) = 0.
     
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,phi,grad)
        
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 grad(i,j,k) = grad(i,j,k)*dx3GM(k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 grad(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
!pgi$ unroll = n:8
                 do kk = g3L+1, g3U
                    grad(i,j,k) = grad(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) grad(1 ,S23B:N23B,S33B:N33B) = 0.
     if (BC_1U > 0) grad(N1,S23B:N23B,S33B:N33B) = 0.
     if (BC_2L > 0) grad(S13B:N13B,1 ,S33B:N33B) = 0.
     if (BC_2U > 0) grad(S13B:N13B,N2,S33B:N33B) = 0.
     if (BC_3L > 0) grad(S13B:N13B,S23B:N23B,0 ) = 0.
     if (BC_3U > 0) grad(S13B:N13B,S23B:N23B,N3) = 0.
     
  end if
  !===========================================================================================================
  
  
  end subroutine gradient
  
  
  
  
  
  
  
  subroutine gradient_transp(m,phi,grad)
  
  implicit none
  
  integer, intent(in   ) ::  m
  
  real   , intent(inout) ::  phi (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(  out) ::  grad(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Umgekehrte Reihenfolge im Vergleich zu Subroutine "gradient".                             !
  !              - Randbedingungen könnten nur zum Teil in die Stencils eingebaut werden, so dass sich das   !
  !                vermutlich nicht wirklich lohnt.                                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  if (m == 1) then
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) phi(0 ,S21B:N21B,S31B:N31B) = 0.
     if (BC_1U > 0) phi(N1,S21B:N21B,S31B:N31B) = 0.
     if (BC_2L > 0) phi(S11B:N11B,1 ,S31B:N31B) = 0.
     if (BC_2U > 0) phi(S11B:N11B,N2,S31B:N31B) = 0.
     if (BC_3L > 0) phi(S11B:N11B,S21B:N21B,1 ) = 0.
     if (BC_3U > 0) phi(S11B:N11B,S21B:N21B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        
        if (BC_1L > 0) com(0 ,S21B:N21B,S31B:N31B) = 0.
        if (BC_1U > 0) com(N1,S21B:N21B,S31B:N31B) = 0.
        if (BC_2L > 0) com(S11B:N11B,1 ,S31B:N31B) = 0.
        if (BC_2U > 0) com(S11B:N11B,N2,S31B:N31B) = 0.
        if (BC_3L > 0) com(S11B:N11B,S21B:N21B,1 ) = 0.
        if (BC_3U > 0) com(S11B:N11B,S21B:N21B,N3) = 0.
        
        do k = S31, N31 ! Intervall ist ok!
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 com(i,j,k) = phi(i,j,k)*dx1GM(i)
              end do
           end do
        end do
        
        call apply_compact_transp(1,1,S11,S21,S31,N11,N21,N31,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cGp1CLT,cGp1CLT_LU,cGp1CRT,WGp1T,SGp1T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 grad(i,j,k) = cGp1T(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
                 do ii = d1L+1, d1U
                    grad(i,j,k) = grad(i,j,k) + cGp1T(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (m == 2) then
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) phi(1 ,S22B:N22B,S32B:N32B) = 0.
     if (BC_1U > 0) phi(N1,S22B:N22B,S32B:N32B) = 0.
     if (BC_2L > 0) phi(S12B:N12B,0 ,S32B:N32B) = 0.
     if (BC_2U > 0) phi(S12B:N12B,N2,S32B:N32B) = 0.
     if (BC_3L > 0) phi(S12B:N12B,S22B:N22B,1 ) = 0.
     if (BC_3U > 0) phi(S12B:N12B,S22B:N22B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        
        if (BC_1L > 0) com(1 ,S22B:N22B,S32B:N32B) = 0.
        if (BC_1U > 0) com(N1,S22B:N22B,S32B:N32B) = 0.
        if (BC_2L > 0) com(S12B:N12B,0 ,S32B:N32B) = 0.
        if (BC_2U > 0) com(S12B:N12B,N2,S32B:N32B) = 0.
        if (BC_3L > 0) com(S12B:N12B,S22B:N22B,1 ) = 0.
        if (BC_3U > 0) com(S12B:N12B,S22B:N22B,N3) = 0.
        
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 com(i,j,k) = phi(i,j,k)*dx2GM(j)
              end do
           end do
        end do
        
        call apply_compact_transp(2,2,S12,S22,S32,N12,N22,N32,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cGp2CLT,cGp2CLT_LU,cGp2CRT,WGp2T,SGp2T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
!pgi$ unroll = n:8
                 do jj = d2L, d2U
                    grad(i,j,k) = grad(i,j,k) + cGp2T(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (m == 3) then
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (BC_1L > 0) phi(1 ,S23B:N23B,S33B:N33B) = 0.
     if (BC_1U > 0) phi(N1,S23B:N23B,S33B:N33B) = 0.
     if (BC_2L > 0) phi(S13B:N13B,1 ,S33B:N33B) = 0.
     if (BC_2U > 0) phi(S13B:N13B,N2,S33B:N33B) = 0.
     if (BC_3L > 0) phi(S13B:N13B,S23B:N23B,0 ) = 0.
     if (BC_3U > 0) phi(S13B:N13B,S23B:N23B,N3) = 0.
     
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        
        if (BC_1L > 0) com(1 ,S23B:N23B,S33B:N33B) = 0.
        if (BC_1U > 0) com(N1,S23B:N23B,S33B:N33B) = 0.
        if (BC_2L > 0) com(S13B:N13B,1 ,S33B:N33B) = 0.
        if (BC_2U > 0) com(S13B:N13B,N2,S33B:N33B) = 0.
        if (BC_3L > 0) com(S13B:N13B,S23B:N23B,0 ) = 0.
        if (BC_3U > 0) com(S13B:N13B,S23B:N23B,N3) = 0.
        
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 com(i,j,k) = phi(i,j,k)*dx3GM(k)
              end do
           end do
        end do
        
        call apply_compact_transp(3,3,S13,S23,S33,N13,N23,N33,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cGp3CLT,cGp3CLT_LU,cGp3CRT,WGp3T,SGp3T,com,grad)
     !--------------------------------------------------------------------------------------------------------
     else
        call exchange(m,m,phi) ! Muss nach Randbedingungen kommen!
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
!pgi$ unroll = n:8
                 do kk = d3L, d3U
                    grad(i,j,k) = grad(i,j,k) + cGp3T(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  
  
  end subroutine gradient_transp
  
  
  
  
  
  
  
  
  
  
  !>  \brief computes \f$ \mathrm{lap_m = mulI phi_m - mulL \Delta phi_m} \f$
  !!
  !! used for mod_rhs and for product_Helmholtz(so boundary conditions are included?)
  !! \param[in] m dimension from one to three
  !! \param[in] mulI factor which coresponds to the factor of the identity part
  !! \param[in] mulL factor which coresponds to the factor of the laplace part
  !! \param[inout] phi = vel in case of mod_rhs
  !! \param[out] Lap
  !! \todo change for discrete forcing
  subroutine Helmholtz( m, mulI, mulL, phi, Lap ) bind (c,name='OP_helmholtz')
  
  implicit none
  
  integer(c_int)  , intent(in   ) ::  m
  real(c_double)  , intent(in   ) ::  mulI
  real(c_double)  , intent(in   ) ::  mulL

  real(c_double)  , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double)  , intent(  out) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_visc_yes) then
        if (dimens == 3) then
           call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,com)
           call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,phi,dig)
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    Lap(i,j,k) = dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i)
                 end do
              end do
           end do
           call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
           call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                 end do
              end do
           end do
           call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
           call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*Lap(i,j,k)
                 end do
              end do
           end do
        else
           call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,com)
           call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,phi,dig)
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    Lap(i,j,k) = dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i)
                 end do
              end do
           end do
           call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
           call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*Lap(i,j,k)
                 end do
              end do
           end do
        end if
     !--------------------------------------------------------------------------------------------------------
     else
        if (dimens == 3) then
           do k = S31, N31
              do j = S21, N21
                 do i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    end do
!pgi$ unroll = n:8
                    do jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    end do
!pgi$ unroll = n:8
                    do kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    end do
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                 end do
              end do
           end do
        else
           do k = S31, N31
              do j = S21, N21
                 do i = S11, N11
                    dd1 = cu11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cu11(ii,i)*phi(i+ii,j,k)
                    end do
!pgi$ unroll = n:8
                    do jj = b2L, b2U
                       dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                    end do
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_visc_yes) then
        if (dimens == 3) then
           call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
           call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
                 end do
              end do
           end do
           call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,com)
           call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,phi,dig)
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j)
                 end do
              end do
           end do
           call apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
           call apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*Lap(i,j,k)
                 end do
              end do
           end do
        else
           call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
           call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
                 end do
              end do
           end do
           call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,com)
           call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,phi,dig)
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j)
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*Lap(i,j,k)
                 end do
              end do
           end do
        end if
     !--------------------------------------------------------------------------------------------------------
     else
        if (dimens == 3) then
           do k = S32, N32
              do j = S22, N22
                 do i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    end do
!pgi$ unroll = n:8
                    do jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    end do
!pgi$ unroll = n:8
                    do kk = b3L, b3U
                       dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                    end do
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                 end do
              end do
           end do
        else
           do k = S32, N32
              do j = S22, N22
                 do i = S12, N12
                    dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                    end do
!pgi$ unroll = n:8
                    do jj = b2L, b2U
                       dd1 = dd1 + cv22(jj,j)*phi(i,j+jj,k)
                    end do
                    Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3 .and. dimens == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_visc_yes) then
        call apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
        call apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              end do
           end do
        end do
        call apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
        call apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
              end do
           end do
        end do
        call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,phi,com)
        call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,phi,dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3wM(k)**2 + com(i,j,k)*ddx3wM(k)
                 Lap(i,j,k) = mulI*phi(i,j,k) - mulL*Lap(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*phi(i,j,k+kk)
                 end do
                 Lap(i,j,k) = mulI*phi(i,j,k) - mulL*dd1
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  subroutine Helmholtz_explicit(exch_yes)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange_all_all(.true.,vel)
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (comp_visc_yes) then
     call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,vel(b1L,b2L,b3L,1),com)
     call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,vel(b1L,b2L,b3L,1),dig)
     do k = S31, N31
        do j = S21, N21
!pgi$ unroll = n:8
           do i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx1uM(i)**2 + com(i,j,k)*ddx1uM(i))
           end do
        end do
     end do
     call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,1),com)
     call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,vel(b1L,b2L,b3L,1),dig)
     do k = S31, N31
        do j = S21, N21
!pgi$ unroll = n:8
           do i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           end do
        end do
     end do
     
     if (dimens == 3) then
     
     call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,1),com)
     call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,vel(b1L,b2L,b3L,1),dig)
     do k = S31, N31
        do j = S21, N21
!pgi$ unroll = n:8
           do i = S11, N11
              nl(i,j,k,1) = nl(i,j,k,1) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           end do
        end do
     end do
     
     end if
  !-----------------------------------------------------------------------------------------------------------
  else
     if (dimens == 3) then
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,1)
                 end do
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              end do
           end do
        end do
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 dd1 = cu11(b1L,i)*vel(i+b1L,j,k,1)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cu11(ii,i)*vel(i+ii,j,k,1)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,1)
                 end do
                 nl(i,j,k,1) = nl(i,j,k,1) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (comp_visc_yes) then
     call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,2),com)
     call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,vel(b1L,b2L,b3L,2),dig)
     do k = S32, N32
        do j = S22, N22
!pgi$ unroll = n:8
           do i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           end do
        end do
     end do
     call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,vel(b1L,b2L,b3L,2),com)
     call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,vel(b1L,b2L,b3L,2),dig)
     do k = S32, N32
        do j = S22, N22
!pgi$ unroll = n:8
           do i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx2vM(j)**2 + com(i,j,k)*ddx2vM(j))
           end do
        end do
     end do
     
     if (dimens == 3) then
     
     call apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,vel(b1L,b2L,b3L,2),com)
     call apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,vel(b1L,b2L,b3L,2),dig)
     do k = S32, N32
        do j = S22, N22
!pgi$ unroll = n:8
           do i = S12, N12
              nl(i,j,k,2) = nl(i,j,k,2) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           end do
        end do
     end do
     
     end if
  !-----------------------------------------------------------------------------------------------------------
  else
     if (dimens == 3) then
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*vel(i,j,k+kk,2)
                 end do
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              end do
           end do
        end do
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,2)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cv22(jj,j)*vel(i,j+jj,k,2)
                 end do
                 nl(i,j,k,2) = nl(i,j,k,2) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
     if (comp_visc_yes) then
        call apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,vel(b1L,b2L,b3L,3),com)
        call apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,vel(b1L,b2L,b3L,3),dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
              end do
           end do
        end do
        call apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,vel(b1L,b2L,b3L,3),com)
        call apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,vel(b1L,b2L,b3L,3),dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
              end do
           end do
        end do
        call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,vel(b1L,b2L,b3L,3),com)
        call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,vel(b1L,b2L,b3L,3),dig)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*(dig(i,j,k)*dx3wM(k)**2 + com(i,j,k)*ddx3wM(k))
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = cp11(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*vel(i+ii,j,k,3)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*vel(i,j+jj,k,3)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cw33(kk,k)*vel(i,j,k+kk,3)
                 end do
                 nl(i,j,k,3) = nl(i,j,k,3) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  end subroutine Helmholtz_explicit
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! relativ ungetestet ...
  subroutine Helmholtz_pre_explicit(exch_yes,phi,Lap)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  dd1
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) then
     call exchange(1,0,phi)
     call exchange(2,0,phi)
     call exchange(3,0,phi)
  end if
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (comp_visc_yes) then
     call apply_compact(1,0,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cp1cL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,com)
     call apply_compact(1,0,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,phi,dig)
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           end do
        end do
     end do
     call apply_compact(2,0,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,com)
     call apply_compact(2,0,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,phi,dig)
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           end do
        end do
     end do
     
     if (dimens == 3) then
     
     call apply_compact(3,0,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,com)
     call apply_compact(3,0,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,phi,dig)
     do k = S3p, N3p
        do j = S2p, N2p
!pgi$ unroll = n:8
           do i = S1p, N1p
              Lap(i,j,k) = Lap(i,j,k) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           end do
        end do
     end do
     
     end if
  !-----------------------------------------------------------------------------------------------------------
  else
     if (dimens == 3) then
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cp33(kk,k)*phi(i,j,k+kk)
                 end do
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              end do
           end do
        end do
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 dd1 = cp11(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp11(ii,i)*phi(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cp22(jj,j)*phi(i,j+jj,k)
                 end do
                 Lap(i,j,k) = Lap(i,j,k) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  end subroutine Helmholtz_pre_explicit
  
  
  
  
  
  
  
  
  
  
  
  subroutine Helmholtz_conc(exch_yes,phi,Lap) ! ok (18.02.09)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(  out) ::  Lap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  dd1
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  integer                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  if (exch_yes) then
     call exchange(1,0,phi)
     call exchange(2,0,phi)
     call exchange(3,0,phi)
  end if
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (comp_visc_yes) then
     if (dimens == 3) then
        call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),phi,com)
        call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),phi,dig)
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              do i = S1c(m), N1c(m)
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              end do
           end do
        end do
        call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),phi,com)
        call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),phi,dig)
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              do i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
              end do
           end do
        end do
        call apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),phi,com)
        call apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m),Wc33(1,0,m),Sc33(1,1,m),phi,dig)
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              do i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k)
                 Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
              end do
           end do
        end do
     else
        call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),phi,com)
        call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),phi,dig)
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              do i = S1c(m), N1c(m)
                 Lap(i,j,k) = dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i)
              end do
           end do
        end do
        call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),phi,com)
        call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),phi,dig)
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
              do i = S1c(m), N1c(m)
                 Lap(i,j,k) = Lap(i,j,k) + dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j)
                 Lap(i,j,k) = phi(i,j,k) - multL*Lap(i,j,k)
              end do
           end do
        end do
     end if
  !-----------------------------------------------------------------------------------------------------------
  else
     if (dimens == 3) then
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*phi(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*phi(i,j+jj,k)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cc33(kk,k,m)*phi(i,j,k+kk)
                 end do
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              end do
           end do
        end do
     else
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*phi(i+ii,j,k)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*phi(i,j+jj,k)
                 end do
                 Lap(i,j,k) = phi(i,j,k) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  end subroutine Helmholtz_conc
  
  
  
  
  
  
  
  
  
  
  
  subroutine Helmholtz_conc_explicit(exch_yes,nlco)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  real   , intent(inout) ::  nlco(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real                   ::  dd1
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  integer                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Für den Aufbau der RHS bei expliziter Zeitintegration.                                    !
  !              - Randbedingungen müssen daher nicht berücksichtigt werden.                                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  if (exch_yes) then
     call exchange(1,0,conc(b1L,b2L,b3L,m))
     call exchange(2,0,conc(b1L,b2L,b3L,m))
     call exchange(3,0,conc(b1L,b2L,b3L,m))
  end if
  
  
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  if (comp_visc_yes) then
     call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),com)
     call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m),Wc11(1,0,m),Sc11(1,1,m),conc(b1L,b2L,b3L,m),dig)
     do k = S3c(m), N3c(m)
        do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           do i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx1pM(i)**2 + com(i,j,k)*ddx1pM(i))
           end do
        end do
     end do
     call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),com)
     call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m),Wc22(1,0,m),Sc22(1,1,m),conc(b1L,b2L,b3L,m),dig)
     do k = S3c(m), N3c(m)
        do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           do i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx2pM(j)**2 + com(i,j,k)*ddx2pM(j))
           end do
        end do
     end do
     
     if (dimens == 3) then
     
     call apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),conc(b1L,b2L,b3L,m),com)
     call apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m),Wc33(1,0,m),Sc33(1,1,m),conc(b1L,b2L,b3L,m),dig)
     do k = S3c(m), N3c(m)
        do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
           do i = S1c(m), N1c(m)
              !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
              nlco(i,j,k) = nlco(i,j,k) - multL*(dig(i,j,k)*dx3pM(k)**2 + com(i,j,k)*ddx3pM(k))
           end do
        end do
     end do
     
     end if
  !-----------------------------------------------------------------------------------------------------------
  else
     if (dimens == 3) then
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 end do
!pgi$ unroll = n:8
                 do kk = b3L, b3U
                    dd1 = dd1 + cc33(kk,k,m)*conc(i,j,k+kk,m)
                 end do
                 !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*dd1
                 nlco(i,j,k) = nlco(i,j,k) - multL*dd1
              end do
           end do
        end do
     else
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 dd1 = cc11(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cc11(ii,i,m)*conc(i+ii,j,k,m)
                 end do
!pgi$ unroll = n:8
                 do jj = b2L, b2U
                    dd1 = dd1 + cc22(jj,j,m)*conc(i,j+jj,k,m)
                 end do
                 !nlco(i,j,k,m) = nlco(i,j,k,m) - multL*dd1
                 nlco(i,j,k) = nlco(i,j,k) - multL*dd1
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  
  
  end subroutine Helmholtz_conc_explicit
  
  
  
  
  
  
  
  
  
  !> \brief computes nonlinear terms.
  !! \todo:
  !!    - (parameter in dimension)
  !!    - parameter in work123
  !!    - (parameter2 in work123)
  !!    - exchange vel/work123
  !!    - interpolate_vel
  !!
  !! first interpolates worki to pp, then \f$\partial_i\f$ \c vel then
  !! \f[ \mathrm{nl = pp*\partial_i vel}\f]
  !! \test Teile davon (nur zentrale Operationen!) koennten jeweils durch interpolate2_pre_vel/interpolate2_vel_pre, first_adv_pre/first_adv_vel
  !!         ersetzt werden (beachte aber Addition von nl!)
  !! \test umbenennen in advect... (?)
  !! \relates Pimpact::Nonlinear
  !! \todo extract interpolate vel
  subroutine nonlinear(     &
        phi1U,phi1V,phi1W,  &
        phi2U,phi2V,phi2W,  &
        nlU,nlV,nlW ) bind (c,name='OP_nonlinear')
  
  implicit none
  
  
  real(c_double),  intent(inout) ::  phi1U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(inout) ::  phi1V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(inout) ::  phi1W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(inout) ::  phi2U(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(inout) ::  phi2V(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(inout) ::  phi2W(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real(c_double),  intent(  out) ::  nlU(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(  out) ::  nlV(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double),  intent(  out) ::  nlW(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))

  real                   ::  dd1
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - [Sij,Nij] ist immer eine Untermenge von [Sip,Nip]                                         !
  !              - Feld "res" wird mehrfach ausgetauscht, um mit einem statt drei skalaren Feldern arbeiten  !
  !                zu können. Im Prinzip könnte auch rhs(:,:,:,1:3) für die Zwischenspeicherung verwendet    !
  !                werden, jedoch sind dazu einige Umbaumassnahmen notwendig bei geringem Effizienzgewinn.   !
  !----------------------------------------------------------------------------------------------------------!

                   call interpolate2_vel_pre(.false.,1,phi1U(b1L,b2L,b3L),work1)
                   call interpolate2_vel_pre(.false.,2,phi1V(b1L,b2L,b3L),work2)
  if (dimens == 3) call interpolate2_vel_pre(.false.,3,phi1W(b1L,b2L,b3L),work3)


  call exchange(1,0,work1)
  call exchange(2,0,work1)
  call exchange(3,0,work1)

  call exchange(1,0,work2)
  call exchange(2,0,work2)
  call exchange(3,0,work2)

  if (dimens == 3) then
     call exchange(1,0,work3)
     call exchange(2,0,work3)
     call exchange(3,0,work3)
  end if
  
  
  !===========================================================================================================
  !=== u*du/dx ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,1,work1,pp) ! work1 -> pp, why?
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNu1U(n1L,i)*phi2U(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNu1U(ii,i)*phi2U(i+ii,j,k)
                 end do
              else
                 dd1 = cNu1D(n1L,i)*phi2U(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNu1D(ii,i)*phi2U(i+ii,j,k)
                 end do
              end if
              
              nlU(i,j,k) = dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi2U(b1L,b2L,b3L),rr)
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 nlU(i,j,k) = pp(i,j,k)*rr(i,j,k)*dx1uM(i)
              end do
           end do
        end do
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 dd1 = cu1(b1L,i)*phi2U(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cu1(ii,i)*phi2U(i+ii,j,k)
                 end do
                 
                 nlU(i,j,k) = dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*du/dy ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,1,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp2U(n2L,j)*phi2U(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*phi2U(i,j+jj,k)
                 end do
              else
                 dd1 = cNp2D(n2L,j)*phi2U(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*phi2U(i,j+jj,k)
                 end do
              end if
              
              nlU(i,j,k) = nlU(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi2U(b1L,b2L,b3L),rr)
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 nlU(i,j,k) = nlU(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx2pM(j)
              end do
           end do
        end do
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 dd1 = cp2(b2L,j)*phi2U(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*phi2U(i,j+jj,k)
                 end do
                 
                 nlU(i,j,k) = nlU(i,j,k) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  if (dimens == 3) then
  !===========================================================================================================
  !=== w*du/dz ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,1,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp3U(n3L,k)*phi2U(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*phi2U(i,j,k+kk)
                 end do
              else
                 dd1 = cNp3D(n3L,k)*phi2U(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*phi2U(i,j,k+kk)
                 end do
              end if
              
              nlU(i,j,k) = nlU(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi2U(b1L,b2L,b3L),rr)
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 nlU(i,j,k) = nlU(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx3pM(k)
              end do
           end do
        end do
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 dd1 = cp3(b3L,k)*phi2U(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*phi2U(i,j,k+kk)
                 end do
                 
                 nlU(i,j,k) = nl(i,j,k,1) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  end if
  
  
  
  
  
  
  !===========================================================================================================
  !=== u*dv/dx ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,2,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp1U(n1L,i)*phi2V(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*phi2V(i+ii,j,k)
                 end do
              else
                 dd1 = cNp1D(n1L,i)*phi2V(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*phi2V(i+ii,j,k)
                 end do
              end if
              
              nlV(i,j,k) = dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi2V(b1L,b2L,b3L),rr)
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 nlV(i,j,k) = pp(i,j,k)*rr(i,j,k)*dx1pM(i)
              end do
           end do
        end do
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 dd1 = cp1(b1L,i)*phi2V(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*phi2V(i+ii,j,k)
                 end do
                 
                 nlV(i,j,k) = dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dv/dy ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,2,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNv2U(n2L,j)*phi2V(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNv2U(jj,j)*phi2V(i,j+jj,k)
                 end do
              else
                 dd1 = cNv2D(n2L,j)*phi2V(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNv2D(jj,j)*phi2V(i,j+jj,k)
                 end do
              end if
              
              nlV(i,j,k) = nlV(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi2V(b1L,b2L,b3L),rr)
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 nlV(i,j,k) = nlV(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx2vM(j)
              end do
           end do
        end do
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 dd1 = cv2(b2L,j)*phi2V(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    dd1 = dd1 + cv2(jj,j)*phi2V(i,j+jj,k)
                 end do
                 
                 nlV(i,j,k) = nlV(i,j,k) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  if (dimens == 3) then
  !===========================================================================================================
  !=== w*dv/dz ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,2,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp3U(n3L,k)*phi2V(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNp3U(kk,k)*phi2V(i,j,k+kk)
                 end do
              else
                 dd1 = cNp3D(n3L,k)*phi2V(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNp3D(kk,k)*phi2V(i,j,k+kk)
                 end do
              end if
              
              nlV(i,j,k) = nlV(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(3,0,S12,S22,S32,N12,N22,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi2V(b1L,b2L,b3L),rr)
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 nlV(i,j,k) = nlV(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx3pM(k)
              end do
           end do
        end do
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 dd1 = cp3(b3L,k)*phi2V(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    dd1 = dd1 + cp3(kk,k)*phi2V(i,j,k+kk)
                 end do
                 
                 nlV(i,j,k) = nl(i,j,k,2) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  end if
  
  
  
  
  
  if (dimens == 3) then
  !===========================================================================================================
  !=== u*dw/dx ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,3,work1,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S33, N33
        do j = S23, N23
           do i = S13, N13
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp1U(n1L,i)*phi2W(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNp1U(ii,i)*phi2W(i+ii,j,k)
                 end do
              else
                 dd1 = cNp1D(n1L,i)*phi2W(i+n1L,j,k)
!pgi$ unroll = n:8
                 do ii = n1L+1, n1U
                    dd1 = dd1 + cNp1D(ii,i)*phi2W(i+ii,j,k)
                 end do
              end if
              
              nlW(i,j,k) = dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(1,0,S13,S23,S33,N13,N23,N33,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi2W(b1L,b2L,b3L),rr)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nlW(i,j,k) = pp(i,j,k)*rr(i,j,k)*dx1pM(i)
              end do
           end do
        end do
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = cp1(b1L,i)*phi2W(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    dd1 = dd1 + cp1(ii,i)*phi2W(i+ii,j,k)
                 end do
                 
                 nlW(i,j,k) = dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== v*dw/dy ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,3,work2,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S33, N33
        do j = S23, N23
           do i = S13, N13
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNp2U(n2L,j)*phi2W(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNp2U(jj,j)*phi2W(i,j+jj,k)
                 end do
              else
                 dd1 = cNp2D(n2L,j)*phi2W(i,j+n2L,k)
!pgi$ unroll = n:8
                 do jj = n2L+1, n2U
                    dd1 = dd1 + cNp2D(jj,j)*phi2W(i,j+jj,k)
                 end do
              end if
              
              nlW(i,j,k) = nlW(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(2,0,S13,S23,S33,N13,N23,N33,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi2W(b1L,b2L,b3L),rr)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nlW(i,j,k) = nlW(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx2pM(j)
              end do
           end do
        end do
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = cp2(b2L,j)*phi2W(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    dd1 = dd1 + cp2(jj,j)*phi2W(i,j+jj,k)
                 end do
                 
                 nlW(i,j,k) = nlW(i,j,k) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== w*dw/dz ===============================================================================================
  !===========================================================================================================
  call interpolate2_pre_vel(.false.,3,work3,pp)
  !-----------------------------------------------------------------------------------------------------------
  if (upwind_yes) then
     do k = S33, N33
        do j = S23, N23
           do i = S13, N13
              if (pp(i,j,k) >= 0.) then
                 dd1 = cNw3U(n3L,k)*phi2W(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNw3U(kk,k)*phi2W(i,j,k+kk)
                 end do
              else
                 dd1 = cNw3D(n3L,k)*phi2W(i,j,k+n3L)
!pgi$ unroll = n:8
                 do kk = n3L+1, n3U
                    dd1 = dd1 + cNw3D(kk,k)*phi2W(i,j,k+kk)
                 end do
              end if
              
              nlW(i,j,k) = nlW(i,j,k) + dd1*pp(i,j,k)
           end do
        end do
     end do
  else
     if (comp_conv_yes) then
        call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,phi2W(b1L,b2L,b3L),rr)
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 nlW(i,j,k) = nlW(i,j,k) + pp(i,j,k)*rr(i,j,k)*dx3wM(k)
              end do
           end do
        end do
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 dd1 = cw3(b3L,k)*phi2W(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    dd1 = dd1 + cw3(kk,k)*phi2W(i,j,k+kk)
                 end do
                 
                 nlW(i,j,k) = nlW(i,j,k) + dd1*pp(i,j,k)
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  end if
  
!  nlU(S11:N11,S21:N21,S31:N31) = nl(S11:N11,S21:N21,S31:N31,1)
!  nlV(S12:N12,S22:N22,S32:N32) = nl(S12:N12,S22:N22,S32:N32,2)
!  nlW(S13:N13,S23:N23,S33:N33) = nl(S13:N13,S23:N23,S33:N33,3)
  
  end subroutine nonlinear
  
  
  
  
  
  
  
  
  
  
  
  subroutine nonlinear_conc(exch_yes)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  integer                ::  m
  real                   ::  dd1, dd2
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - 
  !----------------------------------------------------------------------------------------------------------!
  
  m = conc_nu
  
  if (exch_yes) then
     call exchange(1,0,conc(b1L,b2L,b3L,m))
     call exchange(2,0,conc(b1L,b2L,b3L,m))
     call exchange(3,0,conc(b1L,b2L,b3L,m))
  end if
  
  !===========================================================================================================
  if (dimens == 3) then
     
     if (upwind_conc_yes) then
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(1,m) + work1(i,j,k)
                 
                 if (dd2 >= 0.) then
                    dd1 = cNc1U(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNc1U(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                 else
                    dd1 = cNc1D(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNc1D(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                 end if
                 
                 nlco(i,j,k,m) = dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(2,m) + work2(i,j,k)
                 
                 if (dd2 >= 0.) then
                    dd1 = cNc2U(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNc2U(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                 else
                    dd1 = cNc2D(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNc2D(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                 end if
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(3,m) + work3(i,j,k)
                 
                 if (dd2 >= 0.) then
                    dd1 = cNc3U(n3L,k,m)*conc(i,j,k+n3L,m)
!pgi$ unroll = n:8
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNc3U(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                 else
                    dd1 = cNc3D(n3L,k,m)*conc(i,j,k+n3L,m)
!pgi$ unroll = n:8
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNc3D(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                 end if
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
              end do
           end do
        end do
     else
        if (comp_conv_yes) then
           call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),pp)
           call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),rr)
           call apply_compact(3,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N3,ndL,ndR,dimS3,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m),Wc3 (1,0,m),Sc3 (1,1,m),conc(b1L,b2L,b3L,m),Ap)
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
                 do i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    nlco(i,j,k,m) =                 pp(i,j,k)*dx1pM(i)*(us_vec(1,m) + work1(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + rr(i,j,k)*dx2pM(j)*(us_vec(2,m) + work2(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + Ap(i,j,k)*dx3pM(k)*(us_vec(3,m) + work3(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 end do
              end do
           end do
        else
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
                 do i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cc1(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                    
                    nlco(i,j,k,m) = dd1*(us_vec(1,m) + work1(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cc2(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
                    do kk = b3L+1, b3U
                       dd1 = dd1 + cc3(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(3,m) + work3(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 end do
              end do
           end do
        end if
     end if
  !===========================================================================================================
  else
     if (upwind_conc_yes) then
        do k = S3c(m), N3c(m)
           do j = S2c(m), N2c(m)
              do i = S1c(m), N1c(m)
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(1,m) + work1(i,j,k)
                 
                 if (dd2 >= 0.) then
                    dd1 = cNc1U(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNc1U(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                 else
                    dd1 = cNc1D(n1L,i,m)*conc(i+n1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNc1D(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                 end if
                 
                 nlco(i,j,k,m) = dd1*dd2
                 !--------------------------------------------------------------------------------------------
                 dd2 = us_vec(2,m) + work2(i,j,k)
                 
                 if (dd2 >= 0.) then
                    dd1 = cNc2U(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNc2U(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                 else
                    dd1 = cNc2D(n2L,j,m)*conc(i,j+n2L,k,m)
!pgi$ unroll = n:8
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNc2D(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                 end if
                 
                 nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*dd2
                 !--------------------------------------------------------------------------------------------
              end do
           end do
        end do
     else
        if (comp_conv_yes) then
           call apply_compact(1,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N1,ndL,ndR,dimS1,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m),Wc1 (1,0,m),Sc1 (1,1,m),conc(b1L,b2L,b3L,m),pp)
           call apply_compact(2,0,S1c(m),S2c(m),S3c(m),N1c(m),N2c(m),N3c(m),N2,ndL,ndR,dimS2,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m),Wc2 (1,0,m),Sc2 (1,1,m),conc(b1L,b2L,b3L,m),rr)
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
!pgi$ unroll = n:8
                 do i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    nlco(i,j,k,m) =                 pp(i,j,k)*dx1pM(i)*(us_vec(1,m) + work1(i,j,k))
                    nlco(i,j,k,m) = nlco(i,j,k,m) + rr(i,j,k)*dx2pM(j)*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 end do
              end do
           end do
        else
           do k = S3c(m), N3c(m)
              do j = S2c(m), N2c(m)
                 do i = S1c(m), N1c(m)
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cc1(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                    
                    nlco(i,j,k,m) = dd1*(us_vec(1,m) + work1(i,j,k))
                    !-----------------------------------------------------------------------------------------
                    dd1 = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cc2(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                    
                    nlco(i,j,k,m) = nlco(i,j,k,m) + dd1*(us_vec(2,m) + work2(i,j,k))
                    !-----------------------------------------------------------------------------------------
                 end do
              end do
           end do
        end if
     end if
  end if
  !===========================================================================================================
  
  
  end subroutine nonlinear_conc
  
  
  
  
  
  
  
  
  
  !> \brief interpolates velocity on pressure grid
  !!
  !! vel(:,:,:,i) -> worki(:,:,:)
  !! \test anderes Modul?
  !! \test basiert neu auf interpolate2_vel_pre ... ok??
  !! \todo pimpactit
  subroutine interpolate_vel(exch_yes)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
                   call interpolate2_vel_pre(exch_yes,1,vel(b1L,b2L,b3L,1),work1)
                   call interpolate2_vel_pre(exch_yes,2,vel(b1L,b2L,b3L,2),work2)
  if (dimens == 3) call interpolate2_vel_pre(exch_yes,3,vel(b1L,b2L,b3L,3),work3)
  
  
  call exchange(1,0,work1)
  call exchange(2,0,work1)
  call exchange(3,0,work1)
  
  call exchange(1,0,work2)
  call exchange(2,0,work2)
  call exchange(3,0,work2)
  
  if (dimens == 3) then
     call exchange(1,0,work3)
     call exchange(2,0,work3)
     call exchange(3,0,work3)
  end if
  
  
  end subroutine interpolate_vel
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  subroutine interpolate_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,inter)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in   ) ::  m
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  real   , intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !              - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen (Extrapolation!).                            !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     if (comp_inter_yes) then
        call apply_compact(1,0,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_1L > 0 .and. SS1 == S11B) then
           i = S11B
           do k = SS3, NN3
              do j = SS2, NN2
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end if
        if (BC_1U > 0 .and. NN1 == N11B) then
           i = N11B
           do k = SS3, NN3
              do j = SS2, NN2
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     else
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S11B, N11B ! TEST!!! hier koennte man auch SS1, NN1 nehmen! gilt auch fuer andere Routinen!!
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (comp_inter_yes) then
        call apply_compact(2,0,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_2L > 0 .and. SS2 == S22B) then
           j = S22B
           do k = SS3, NN3
              do i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end if
        if (BC_2U > 0 .and. NN2 == N22B) then
           j = N22B
           do k = SS3, NN3
              do i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     else
        do k = SS3, NN3
           do j = S22B, N22B ! TEST!!! hier koennte man auch SS2, NN2 nehmen!
              do i = SS1, NN1
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (comp_inter_yes) then
        call apply_compact(3,0,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,phi,inter)
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_3L > 0 .and. SS3 == S33B) then
           k = S33B
           do j = SS2, NN2
              do i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end if
        if (BC_3U > 0 .and. NN3 == N33B) then
           k = N33B
           do j = SS2, NN2
              do i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     else
        do k = S33B, N33B ! TEST!!! hier koennte man auch SS3, NN3 nehmen!
           do j = SS2, NN2
              do i = SS1, NN1
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  end subroutine interpolate_pre_vel
  
  
  
  
  
  
  
  
  
  
  !> \brief interpolates pressure to vel nodes
  !!
  !! \f[ inter = inter + interpolated(phi) \f]
  !! Wie interpolate_pre_vel, allerdings mit fixen Index-Limiten (ohne Rand)
  !! \param[in] exch_yes indicates if fields have exchanged first
  !! \param[in] m dimension
  !! \param[inout] phi input field
  !! \param[inout] inter output field
  subroutine interpolate2_pre_vel(exch_yes,m,phi,inter)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in   ) ::  m
  
  real   , intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  if (exch_yes) call exchange(m,0,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     if (comp_inter_yes) then
        call apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,phi,inter)
     else
        do k = S31, N31
           do j = S21, N21
              do i = S11, N11
                 inter(i,j,k) = cIpu(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    inter(i,j,k) = inter(i,j,k) + cIpu(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (comp_inter_yes) then
        call apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,phi,inter)
     else
        do k = S32, N32
           do j = S22, N22
              do i = S12, N12
                 inter(i,j,k) = cIpv(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    inter(i,j,k) = inter(i,j,k) + cIpv(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (comp_inter_yes) then
        call apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,phi,inter)
     else
        do k = S33, N33
           do j = S23, N23
              do i = S13, N13
                 inter(i,j,k) = cIpw(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    inter(i,j,k) = inter(i,j,k) + cIpw(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  end subroutine interpolate2_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  subroutine interpolate_vel_pre(   &
    m,                              &
    SS1,SS2,SS3,                    &
    NN1,NN2,NN3,                    &
    phi,inter) bind (c,name='F_interpolateVel2Pre')
  
  implicit none
  
!  logical       , intent(in   ) ::  exch_yes
  integer(c_int), intent(in)    ::  m
  integer(c_int), intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  real(c_double), intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real(c_double), intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
!  if (exch_yes) call exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     if (comp_inter_yes) then
        call apply_compact(1,1,S1p,SS2,SS3,N1p,NN2,NN3,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,phi,inter)
     else
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 do ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (comp_inter_yes) then
        call apply_compact(2,2,SS1,S2p,SS3,NN1,N2p,NN3,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,phi,inter)
     else
        do k = SS3, NN3
           do j = S2p, N2p
              do i = SS1, NN1
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 do jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (comp_inter_yes) then
        call apply_compact(3,3,SS1,SS2,S3p,NN1,NN2,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,phi,inter)
     else
        do k = S3p, N3p
           do j = SS2, NN2
              do i = SS1, NN1
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 do kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  end subroutine interpolate_vel_pre
  
  
  
  
  
  
  
  
  
  !> \brief interpolates velocity to pressure grid
  !!
  !! Wie interpolate_vel_pre, allerdings mit fixen Index-Limiten (ohne Rand)
  subroutine interpolate2_vel_pre(exch_yes,m,phi,inter)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in)    ::  m
  
  real   , intent(inout) ::  phi  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  inter(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  if (exch_yes) call exchange(m,m,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     if (comp_inter_yes) then
        call apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,phi,inter)
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 inter(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
                 do ii = d1L+1, d1U
                    inter(i,j,k) = inter(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (comp_inter_yes) then
        call apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,phi,inter)
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 inter(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
                 do jj = d2L+1, d2U
                    inter(i,j,k) = inter(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (comp_inter_yes) then
        call apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,phi,inter)
     else
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 inter(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
                 do kk = d3L+1, d3U
                    inter(i,j,k) = inter(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
  end if
  !===========================================================================================================
  
  
  end subroutine interpolate2_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  subroutine interpolate_conc(exch_yes) ! ok (18.02.09)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer                ::  m
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1, dd2
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -  
  !----------------------------------------------------------------------------------------------------------!
  
  
  do m = 1, n_conc
     
     !========================================================================================================
     if (gravity(1) /= 0.) then
        
        dd1 = Ric(m)*gravity(1)
        
        if (exch_yes) call exchange(1,0,conc(b1L,b2L,b3L,m))
        
        if (comp_inter_yes) then
           call apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIcuCL(-ndL,0,m),cIcuCL_LU(1,0,m),cIcuCR(-ndR,0,m),WIcu(1,0,m),SIcu(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S11:N11,S21:N21,S31:N31,1) = nl(S11:N11,S21:N21,S31:N31,1) - dd1*pp(S11:N11,S21:N21,S31:N31)
        else
           do k = S31, N31
              do j = S21, N21
                 do i = S11, N11
                    dd2 = cIcu(g1L,i,m)*conc(i+g1L,j,k,m)
!pgi$ unroll = n:8
                    do ii = g1L+1, g1U
                       dd2 = dd2 + cIcu(ii,i,m)*conc(i+ii,j,k,m)
                    end do
                    nl(i,j,k,1) = nl(i,j,k,1) - dd1*dd2
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
     if (gravity(2) /= 0.) then
        
        dd1 = Ric(m)*gravity(2)
        
        if (exch_yes) call exchange(2,0,conc(b1L,b2L,b3L,m))
        
        if (comp_inter_yes) then
           call apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIcvCL(-ndL,0,m),cIcvCL_LU(1,0,m),cIcvCR(-ndR,0,m),WIcv(1,0,m),SIcv(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S12:N12,S22:N22,S32:N32,2) = nl(S12:N12,S22:N22,S32:N32,2) - dd1*pp(S12:N12,S22:N22,S32:N32)
        else
           do k = S32, N32
              do j = S22, N22
                 do i = S12, N12
                    dd2 = cIcv(g2L,j,m)*conc(i,j+g2L,k,m)
!pgi$ unroll = n:8
                    do jj = g2L+1, g2U
                       dd2 = dd2 + cIcv(jj,j,m)*conc(i,j+jj,k,m)
                    end do
                    nl(i,j,k,2) = nl(i,j,k,2) - dd1*dd2
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
     if (gravity(3) /= 0. .and. dimens == 3) then
        
        dd1 = Ric(m)*gravity(3)
        
        if (exch_yes) call exchange(3,0,conc(b1L,b2L,b3L,m))
        
        if (comp_inter_yes) then
           call apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIcwCL(-ndL,0,m),cIcwCL_LU(1,0,m),cIcwCR(-ndR,0,m),WIcw(1,0,m),SIcw(1,1,m),conc(b1L,b2L,b3L,m),pp)
           nl(S13:N13,S23:N23,S33:N33,3) = nl(S13:N13,S23:N23,S33:N33,3) - dd1*pp(S13:N13,S23:N23,S33:N33)
        else
           do k = S33, N33
              do j = S23, N23
                 do i = S13, N13
                    dd2 = cIcw(g3L,k,m)*conc(i,j,k+g3L,m)
!pgi$ unroll = n:8
                    do kk = g3L+1, g3U
                       dd2 = dd2 + cIcw(kk,k,m)*conc(i,j,k+kk,m)
                    end do
                    nl(i,j,k,3) = nl(i,j,k,3) - dd1*dd2
                 end do
              end do
           end do
        end if
     end if
     !========================================================================================================
     
  end do
  
  
  end subroutine interpolate_conc
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  subroutine first_pre_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in   ) ::  m
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - die Punkte auf (bzw. hinter) dem Rand der Wand-normalen Komponente werden auch bei        !
  !                kompakter Differenzierung im Feld immer explizit gerechnet, um nur eine Variante          !
  !                kompakter Differenzen abspeichern zu muessen.                                             !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(1,0,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,phi,der)
        
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S11, N11
                 der(i,j,k) = der(i,j,k)*dx1GM(i)
              end do
           end do
        end do
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_1L > 0 .and. SS1 == S11B) then
           i = S11B
           do k = SS3, NN3
              do j = SS2, NN2
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end if
        if (BC_1U > 0 .and. NN1 == N11B) then
           i = N11B
           do k = SS3, NN3
              do j = SS2, NN2
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     else
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S11B, N11B
                 der(i,j,k) = cGp1(g1L,i)*phi(i+g1L,j,k)
                 do ii = g1L+1, g1U
                    der(i,j,k) = der(i,j,k) + cGp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(2,0,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,phi,der)
        
        do k = SS3, NN3
           do j = S22, N22
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx2GM(j)
              end do
           end do
        end do
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_2L > 0 .and. SS2 == S22B) then
           j = S22B
           do k = SS3, NN3
              do i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end if
        if (BC_2U > 0 .and. NN2 == N22B) then
           j = N22B
           do k = SS3, NN3
              do i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     else
        do k = SS3, NN3
           do j = S22B, N22B
              do i = SS1, NN1
                 der(i,j,k) = cGp2(g2L,j)*phi(i,j+g2L,k)
                 do jj = g2L+1, g2U
                    der(i,j,k) = der(i,j,k) + cGp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_grad_yes) then
        call apply_compact(3,0,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,phi,der)
        
        do k = S33, N33
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx3GM(k)
              end do
           end do
        end do
        !+++++++++++++++++++++++++++++++++++++++
        if (BC_3L > 0 .and. SS3 == S33B) then
           k = S33B
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end if
        if (BC_3U > 0 .and. NN3 == N33B) then
           k = N33B
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end if
        !+++++++++++++++++++++++++++++++++++++++
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S33B, N33B
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = cGp3(g3L,k)*phi(i,j,k+g3L)
                 do kk = g3L+1, g3U
                    der(i,j,k) = der(i,j,k) + cGp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine first_pre_vel
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  subroutine first_vel_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  
  integer, intent(in   ) ::  m
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(1,1,S1p,SS2,SS3,N1p,NN2,NN3,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,phi,der)
        
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S1p, N1p
                 der(i,j,k) = der(i,j,k)*dx1DM(i)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S1p, N1p
                 der(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
                 do ii = d1L+1, d1U
                    der(i,j,k) = der(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(2,2,SS1,S2p,SS3,NN1,N2p,NN3,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,phi,der)
        
        do k = SS3, NN3
           do j = S2p, N2p
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx2DM(j)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = SS3, NN3
           do j = S2p, N2p
              do i = SS1, NN1
                 der(i,j,k) = cDv2(d2L,j)*phi(i,j+d2L,k)
                 do jj = d2L+1, d2U
                    der(i,j,k) = der(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (comp_div_yes) then
        call apply_compact(3,3,SS1,SS2,S3p,NN1,NN2,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,phi,der)
        
        do k = S3p, N3p
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = der(i,j,k)*dx3DM(k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        do k = S3p, N3p
           do j = SS2, NN2
              do i = SS1, NN1
                 der(i,j,k) = cDw3(d3L,k)*phi(i,j,k+d3L)
                 do kk = d3L+1, d3U
                    der(i,j,k) = der(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine first_vel_pre
  
  
  
  
  
  
  
  
  
  
  
  ! TEST!!! noch relativ ungetestet!
  subroutine first_adv_pre(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in   ) ::  m
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  logical, intent(in   ) ::  upwind_yes
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(in   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - die Punkte auf dem Rand der Wand-normalen Komponente werden auch bei kompakter            !
  !                Differenzierung im Feld immer explizit gerechnet, um nur eine Variante kompakter          !
  !                Differenzen abspeichern zu muessen.                                                       !
  !              - fuer die Punkte der Wand-normalen Komponente sind S12, N12 nur Platzhalter fuer S1p+1,    !
  !                N1p-1 usw.                                                                                !
  !              - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange2(m,0,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S1p, N1p
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNp1U(n1L,i)*phi(i+n1L,j,k)
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNp1U(ii,i)*phi(i+ii,j,k)
                    end do
                 else
                    dd1 = cNp1D(n1L,i)*phi(i+n1L,j,k)
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNp1D(ii,i)*phi(i+ii,j,k)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           call apply_compact(1,0,S12,SS2,SS3,N12,NN2,NN3,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,phi,der)
           do k = SS3, NN3
              do j = SS2, NN2
                 do i = S12, N12
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx1pM(i)
                 end do
              end do
           end do
           !+++++++++++++++++++++++++++++++++++++++
           if (BC_1L > 0 .and. SS1 == S1p) then
              i = S1p
              do k = SS3, NN3
                 do j = SS2, NN2
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           if (BC_1U > 0 .and. NN1 == N1p) then
              i = N1p
              do k = SS3, NN3
                 do j = SS2, NN2
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           !+++++++++++++++++++++++++++++++++++++++
        else
           do k = SS3, NN3
              do j = SS2, NN2
                 do i = S1p, N1p
                    dd1 = cp1(b1L,i)*phi(i+b1L,j,k)
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cp1(ii,i)*phi(i+ii,j,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = SS3, NN3
           do j = S2p, N2p
              do i = SS1, NN1
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNp2U(n2L,j)*phi(i,j+n2L,k)
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNp2U(jj,j)*phi(i,j+jj,k)
                    end do
                 else
                    dd1 = cNp2D(n2L,j)*phi(i,j+n2L,k)
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNp2D(jj,j)*phi(i,j+jj,k)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           call apply_compact(2,0,SS1,S23,SS3,NN1,N23,NN3,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,phi,der)
           do k = SS3, NN3
              do j = S23, N23
                 do i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx2pM(j)
                 end do
              end do
           end do
           !+++++++++++++++++++++++++++++++++++++++
           if (BC_2L > 0 .and. SS2 == S2p) then
              j = S2p
              do k = SS3, NN3
                 do i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           if (BC_2U > 0 .and. NN2 == N2p) then
              j = N2p
              do k = SS3, NN3
                 do i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           !+++++++++++++++++++++++++++++++++++++++
        else
           do k = SS3, NN3
              do j = S2p, N2p
                 do i = SS1, NN1
                    dd1 = cp2(b2L,j)*phi(i,j+b2L,k)
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cp2(jj,j)*phi(i,j+jj,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = S3p, N3p
           do j = SS2, NN2
              do i = SS1, NN1
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNp3U(n3L,k)*phi(i,j,k+n3L)
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNp3U(kk,k)*phi(i,j,k+kk)
                    end do
                 else
                    dd1 = cNp3D(n3L,k)*phi(i,j,k+n3L)
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNp3D(kk,k)*phi(i,j,k+kk)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           ! siehe "Anmerkungen" fuer Kommentar zu Index-Limiten
           call apply_compact(3,0,SS1,SS2,S32,NN1,NN2,N32,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,phi,der)
           do k = S32, N32
              do j = SS2, NN2
                 do i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx3pM(k)
                 end do
              end do
           end do
           !+++++++++++++++++++++++++++++++++++++++
           if (BC_3L > 0 .and. SS3 == S3p) then
              k = S3p
              do j = SS2, NN2
                 do i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    do kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           if (BC_3U > 0 .and. NN3 == N3p) then
              k = N3p
              do j = SS2, NN2
                 do i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    do kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end if
           !+++++++++++++++++++++++++++++++++++++++
        else
           do k = S3p, N3p
              do j = SS2, NN2
                 do i = SS1, NN1
                    dd1 = cp3(b3L,k)*phi(i,j,k+b3L)
                    do kk = b3L+1, b3U
                       dd1 = dd1 + cp3(kk,k)*phi(i,j,k+kk)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine first_adv_pre
  
  
  
  
  
  
  
  
  
  ! TEST!!! Routine kann Punkte der Rand-normalen Komponente auf dem Rand NICHT behandeln!! (will sie momentan aber auch gar nicht koennen)
  ! TEST!!! noch relativ ungetestet!
  subroutine first_adv_vel(exch_yes,m,SS1,SS2,SS3,NN1,NN2,NN3,phi,der,adv,upwind_yes)
  
  implicit none
  
  logical, intent(in   ) ::  exch_yes
  integer, intent(in   ) ::  m
  integer, intent(in   ) ::  SS1, SS2, SS3, NN1, NN2, NN3
  logical, intent(in   ) ::  upwind_yes
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  der(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(in   ) ::  adv(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - bei kompakter Differenzierung werden immer ganze Linien behandelt.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (exch_yes) call exchange2(m,m,SS1,SS2,SS3,NN1,NN2,NN3,phi)
  
  
  !===========================================================================================================
  if (m == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = SS3, NN3
           do j = SS2, NN2
              do i = S11, N11
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNu1U(n1L,i)*phi(i+n1L,j,k)
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNu1U(ii,i)*phi(i+ii,j,k)
                    end do
                 else
                    dd1 = cNu1D(n1L,i)*phi(i+n1L,j,k)
                    do ii = n1L+1, n1U
                       dd1 = dd1 + cNu1D(ii,i)*phi(i+ii,j,k)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           call apply_compact(1,1,S11,SS2,SS3,N11,NN2,NN3,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,phi,der)
           do k = SS3, NN3
              do j = SS2, NN2
                 do i = S11, N11
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx1uM(i)
                 end do
              end do
           end do
        else
           do k = SS3, NN3
              do j = SS2, NN2
                 do i = S11, N11
                    dd1 = cu1(b1L,i)*phi(i+b1L,j,k)
                    do ii = b1L+1, b1U
                       dd1 = dd1 + cu1(ii,i)*phi(i+ii,j,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  !===========================================================================================================
  end if
  !===========================================================================================================
  if (m == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = SS3, NN3
           do j = S22, N22
              do i = SS1, NN1
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNv2U(n2L,j)*phi(i,j+n2L,k)
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNv2U(jj,j)*phi(i,j+jj,k)
                    end do
                 else
                    dd1 = cNv2D(n2L,j)*phi(i,j+n2L,k)
                    do jj = n2L+1, n2U
                       dd1 = dd1 + cNv2D(jj,j)*phi(i,j+jj,k)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           call apply_compact(2,2,SS1,S22,SS3,NN1,N22,NN3,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,phi,der)
           do k = SS3, NN3
              do j = S22, N22
                 do i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx2vM(j)
                 end do
              end do
           end do
        else
           do k = SS3, NN3
              do j = S22, N22
                 do i = SS1, NN1
                    dd1 = cv2(b2L,j)*phi(i,j+b2L,k)
                    do jj = b2L+1, b2U
                       dd1 = dd1 + cv2(jj,j)*phi(i,j+jj,k)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (m == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (upwind_yes) then
        do k = S33, N33
           do j = SS2, NN2
              do i = SS1, NN1
                 if (adv(i,j,k) >= 0.) then
                    dd1 = cNw3U(n3L,k)*phi(i,j,k+n3L)
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNw3U(kk,k)*phi(i,j,k+kk)
                    end do
                 else
                    dd1 = cNw3D(n3L,k)*phi(i,j,k+n3L)
                    do kk = n3L+1, n3U
                       dd1 = dd1 + cNw3D(kk,k)*phi(i,j,k+kk)
                    end do
                 end if
                 der(i,j,k) = dd1*adv(i,j,k)
              end do
           end do
        end do
     !--------------------------------------------------------------------------------------------------------
     else
        if (comp_conv_yes) then
           call apply_compact(3,3,SS1,SS2,S33,NN1,NN2,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,phi,der)
           do k = S33, N33
              do j = SS2, NN2
                 do i = SS1, NN1
                    der(i,j,k) = adv(i,j,k)*der(i,j,k)*dx3wM(k)
                 end do
              end do
           end do
        else
           do k = S33, N33
              do j = SS2, NN2
                 do i = SS1, NN1
                    dd1 = cw3(b3L,k)*phi(i,j,k+b3L)
                    do kk = b3L+1, b3U
                       dd1 = dd1 + cw3(kk,k)*phi(i,j,k+kk)
                    end do
                    der(i,j,k) = dd1*adv(i,j,k)
                 end do
              end do
           end do
        end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine first_adv_vel
  
  
  
  
  
  
  
  
  
  
  
  subroutine outflow_bc
  
  ! (revised on 06.08.2009)
  
  implicit none
  
  integer                ::  m
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  real                   ::  dd1
  
  
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
  if (outlet(2,1,1) .and. (BC_2L == 1 .or. BC_2L == 3)) then ! TEST!!! Warum eigentlich BC_2L == 3? Wuerde generell nicht outlet(2,1,1) ausreichen?
     j = 1
     do k = S31B, N31B
        do i = S11B, N11B
           nlbc12(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              nlbc12(i,k,1) = nlbc12(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,1)
           end do
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,1)
!pgi$ unroll = n:8
           do ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,1)
           end do
           nlbc12(i,k,1) = dd1*nlbc12(i,k,1)
        end do
     end do
  end if
  if (outlet(2,2,1) .and. (BC_2U == 1 .or. BC_2U == 3)) then
     j = N2
     do k = S31B, N31B
        do i = S11B, N11B
           nlbc12(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,1)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              nlbc12(i,k,2) = nlbc12(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,1)
           end do
           
           dd1 = cIpu(g1L,i)*drift2(i+g1L,k,2)
!pgi$ unroll = n:8
           do ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift2(i+ii,k,2)
           end do
           nlbc12(i,k,2) = dd1*nlbc12(i,k,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,1) .and. (BC_3L == 1 .or. BC_3L == 3)) then
     k = 1
     do j = S21B, N21B
        do i = S11B, N11B
           nlbc13(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              nlbc13(i,j,1) = nlbc13(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,1)
           end do
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,1)
!pgi$ unroll = n:8
           do ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,1)
           end do
           nlbc13(i,j,1) = dd1*nlbc13(i,j,1)
        end do
     end do
  end if
  if (outlet(3,2,1) .and. (BC_3U == 1 .or. BC_3U == 3)) then
     k = N3
     do j = S21B, N21B
        do i = S11B, N11B
           nlbc13(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,1)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              nlbc13(i,j,2) = nlbc13(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,1)
           end do
           
           dd1 = cIpu(g1L,i)*drift3(i+g1L,j,2)
!pgi$ unroll = n:8
           do ii = g1L+1, g1U
              dd1 = dd1 + cIpu(ii,i)*drift3(i+ii,j,2)
           end do
           nlbc13(i,j,2) = dd1*nlbc13(i,j,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(1,1,1) .and. (BC_1L == 1 .or. BC_1L == 2)) then
     i = 0
     do k = S31B, N31B
        do j = S21B, N21B
           nlbc11(j,k,1) = cDu1(d1L,1)*vel(1+d1L,j,k,1)
!pgi$ unroll = n:8
           do ii = d1L+1, d1U
              nlbc11(j,k,1) = nlbc11(j,k,1) + cDu1(ii,1)*vel(1+ii,j,k,1)
           end do
           nlbc11(j,k,1) = drift1(j,k,1)*nlbc11(j,k,1)
        end do
     end do
  end if
  if (outlet(1,2,1) .and. (BC_1U == 1 .or. BC_1U == 2)) then
     i = N1
     do k = S31B, N31B
        do j = S21B, N21B
           nlbc11(j,k,2) = cDu1(d1L,i)*vel(i+d1L,j,k,1)
!pgi$ unroll = n:8
           do ii = d1L+1, d1U
              nlbc11(j,k,2) = nlbc11(j,k,2) + cDu1(ii,i)*vel(i+ii,j,k,1)
           end do
           nlbc11(j,k,2) = drift1(j,k,2)*nlbc11(j,k,2)
        end do
     end do
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,2) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(1,1,2) .and. (BC_1L == 1 .or. BC_1L == 3)) then
     i = 1
     do k = S32B, N32B
        do j = S22B, N22B
           nlbc21(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              nlbc21(j,k,1) = nlbc21(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,2)
           end do
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,1)
!pgi$ unroll = n:8
           do jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,1)
           end do
           nlbc21(j,k,1) = dd1*nlbc21(j,k,1)
        end do
     end do
  end if
  if (outlet(1,2,2) .and. (BC_1U == 1 .or. BC_1U == 3)) then
     i = N1
     do k = S32B, N32B
        do j = S22B, N22B
           nlbc21(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,2)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              nlbc21(j,k,2) = nlbc21(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,2)
           end do
           
           dd1 = cIpv(g2L,j)*drift1(j+g2L,k,2)
!pgi$ unroll = n:8
           do jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift1(j+jj,k,2)
           end do
           nlbc21(j,k,2) = dd1*nlbc21(j,k,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,2) .and. (BC_3L == 1 .or. BC_3L == 3)) then
     k = 1
     do j = S22B, N22B
        do i = S12B, N12B
           nlbc23(i,j,1) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              nlbc23(i,j,1) = nlbc23(i,j,1) + cp3(kk,k)*vel(i,j,k+kk,2)
           end do
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,1)
!pgi$ unroll = n:8
           do jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,1)
           end do
           nlbc23(i,j,1) = dd1*nlbc23(i,j,1)
        end do
     end do
  end if
  if (outlet(3,2,2) .and. (BC_3U == 1 .or. BC_3U == 3)) then
     k = N3
     do j = S22B, N22B
        do i = S12B, N12B
           nlbc23(i,j,2) = cp3(b3L,k)*vel(i,j,k+b3L,2)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              nlbc23(i,j,2) = nlbc23(i,j,2) + cp3(kk,k)*vel(i,j,k+kk,2)
           end do
           
           dd1 = cIpv(g2L,j)*drift3(i,j+g2L,2)
!pgi$ unroll = n:8
           do jj = g2L+1, g2U
              dd1 = dd1 + cIpv(jj,j)*drift3(i,j+jj,2)
           end do
           nlbc23(i,j,2) = dd1*nlbc23(i,j,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(2,1,2) .and. (BC_2L == 1 .or. BC_2L == 2)) then
     j = 0
     do k = S32B, N32B
        do i = S12B, N12B
           nlbc22(i,k,1) = cDv2(d2L,1)*vel(i,1+d2L,k,2)
!pgi$ unroll = n:8
           do jj = d2L+1, d2U
              nlbc22(i,k,1) = nlbc22(i,k,1) + cDv2(jj,1)*vel(i,1+jj,k,2)
           end do
           nlbc22(i,k,1) = drift2(i,k,1)*nlbc22(i,k,1)
        end do
     end do
  end if
  if (outlet(2,2,2) .and. (BC_2U == 1 .or. BC_2U == 2)) then
     j = N2
     do k = S32B, N32B
        do i = S12B, N12B
           nlbc22(i,k,2) = cDv2(d2L,j)*vel(i,j+d2L,k,2)
!pgi$ unroll = n:8
           do jj = d2L+1, d2U
              nlbc22(i,k,2) = nlbc22(i,k,2) + cDv2(jj,j)*vel(i,j+jj,k,2)
           end do
           nlbc22(i,k,2) = drift2(i,k,2)*nlbc22(i,k,2)
        end do
     end do
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== vel(:,:,:,3) ==========================================================================================
  !===========================================================================================================
  !-----------------------------------------------------------------------------------------------------------
  !--- Tangentialkomponenten ---------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(1,1,3) .and. (BC_1L == 1 .or. BC_1L == 3)) then
     i = 1
     do k = S33B, N33B
        do j = S23B, N23B
           nlbc31(j,k,1) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              nlbc31(j,k,1) = nlbc31(j,k,1) + cp1(ii,i)*vel(i+ii,j,k,3)
           end do
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,1)
!pgi$ unroll = n:8
           do kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,1)
           end do
           nlbc31(j,k,1) = dd1*nlbc31(j,k,1)
        end do
     end do
  end if
  if (outlet(1,2,3) .and. (BC_1U == 1 .or. BC_1U == 3)) then
     i = N1
     do k = S33B, N33B
        do j = S23B, N23B
           nlbc31(j,k,2) = cp1(b1L,i)*vel(i+b1L,j,k,3)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              nlbc31(j,k,2) = nlbc31(j,k,2) + cp1(ii,i)*vel(i+ii,j,k,3)
           end do
           
           dd1 = cIpw(g3L,k)*drift1(j,k+g3L,2)
!pgi$ unroll = n:8
           do kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift1(j,k+kk,2)
           end do
           nlbc31(j,k,2) = dd1*nlbc31(j,k,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(2,1,3) .and. (BC_2L == 1 .or. BC_2L == 3)) then
     j = 1
     do k = S33B, N33B
        do i = S13B, N13B
           nlbc32(i,k,1) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              nlbc32(i,k,1) = nlbc32(i,k,1) + cp2(jj,j)*vel(i,j+jj,k,3)
           end do
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,1)
!pgi$ unroll = n:8
           do kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,1)
           end do
           nlbc32(i,k,1) = dd1*nlbc32(i,k,1)
        end do
     end do
  end if
  if (outlet(2,2,3) .and. (BC_2U == 1 .or. BC_2U == 3)) then
     j = N2
     do k = S33B, N33B
        do i = S13B, N13B
           nlbc32(i,k,2) = cp2(b2L,j)*vel(i,j+b2L,k,3)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              nlbc32(i,k,2) = nlbc32(i,k,2) + cp2(jj,j)*vel(i,j+jj,k,3)
           end do
           
           dd1 = cIpw(g3L,k)*drift2(i,k+g3L,2)
!pgi$ unroll = n:8
           do kk = g3L+1, g3U
              dd1 = dd1 + cIpw(kk,k)*drift2(i,k+kk,2)
           end do
           nlbc32(i,k,2) = dd1*nlbc32(i,k,2)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  !--- Normalkomponente --------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (outlet(3,1,3) .and. (BC_3L == 1 .or. BC_3L == 2)) then
     k = 0
     do j = S23B, N23B
        do i = S13B, N13B
           nlbc33(i,j,1) = cDw3(d3L,1)*vel(i,j,1+d3L,3)
!pgi$ unroll = n:8
           do kk = d3L+1, d3U
              nlbc33(i,j,1) = nlbc33(i,j,1) + cDw3(kk,1)*vel(i,j,1+kk,3)
           end do
           nlbc33(i,j,1) = drift3(i,j,1)*nlbc33(i,j,1)
        end do
     end do
  end if
  if (outlet(3,2,3) .and. (BC_3U == 1 .or. BC_3U == 2)) then
     k = N3
     do j = S23B, N23B
        do i = S13B, N13B
           nlbc33(i,j,2) = cDw3(d3L,k)*vel(i,j,k+d3L,3)
!pgi$ unroll = n:8
           do kk = d3L+1, d3U
              nlbc33(i,j,2) = nlbc33(i,j,2) + cDw3(kk,k)*vel(i,j,k+kk,3)
           end do
           nlbc33(i,j,2) = drift3(i,j,2)*nlbc33(i,j,2)
        end do
     end do
  end if
  !===========================================================================================================
  
  
  end subroutine outflow_bc
  
  
  
  
  
  
  
  
  
  
  
  subroutine sediment_bc
  
  ! (revised on 06.08.2009)
  
  implicit none
  
  integer                ::  m
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: -                                                                                           !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  !===========================================================================================================
  if (BCc_1L(m) == 1 .and. (.not. isopycnal(1,1,m))) then
     i = 1
     do k = S3p, N3p
        do j = S2p, N2p
           sed1L(j,k,m) = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              sed1L(j,k,m) = sed1L(j,k,m) + cc1(ii,i,m)*conc(i+ii,j,k,m)
           end do
           sed1L(j,k,m) = sed1L(j,k,m)*(us_vec(1,m) + work1(i,j,k))
        end do
     end do
  end if
  
  if (BCc_1U(m) == 1 .and. (.not. isopycnal(1,2,m))) then
     i = N1
     do k = S3p, N3p
        do j = S2p, N2p
           sed1U(j,k,m) = cc1(b1L,i,m)*conc(i+b1L,j,k,m)
!pgi$ unroll = n:8
           do ii = b1L+1, b1U
              sed1U(j,k,m) = sed1U(j,k,m) + cc1(ii,i,m)*conc(i+ii,j,k,m)
           end do
           sed1U(j,k,m) = sed1U(j,k,m)*(us_vec(1,m) + work1(i,j,k))
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_2L(m) == 1 .and. (.not. isopycnal(2,1,m))) then
     j = 1
     do k = S3p, N3p
        do i = S1p, N1p
           sed2L(i,k,m) = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              sed2L(i,k,m) = sed2L(i,k,m) + cc2(jj,j,m)*conc(i,j+jj,k,m)
           end do
           sed2L(i,k,m) = sed2L(i,k,m)*(us_vec(2,m) + work2(i,j,k))
        end do
     end do
  end if
  
  if (BCc_2U(m) == 1 .and. (.not. isopycnal(2,2,m))) then
     j = N2
     do k = S3p, N3p
        do i = S1p, N1p
           sed2U(i,k,m) = cc2(b2L,j,m)*conc(i,j+b2L,k,m)
!pgi$ unroll = n:8
           do jj = b2L+1, b2U
              sed2U(i,k,m) = sed2U(i,k,m) + cc2(jj,j,m)*conc(i,j+jj,k,m)
           end do
           sed2U(i,k,m) = sed2U(i,k,m)*(us_vec(2,m) + work2(i,j,k))
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_3L(m) == 1 .and. (.not. isopycnal(3,1,m))) then
     k = 1
     do j = S2p, N2p
        do i = S1p, N1p
           sed3L(i,j,m) = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              sed3L(i,j,m) = sed3L(i,j,m) + cc3(kk,k,m)*conc(i,j,k+kk,m)
           end do
           sed3L(i,j,m) = sed3L(i,j,m)*(us_vec(3,m) + work3(i,j,k))
        end do
     end do
  end if
  
  if (BCc_3U(m) == 1 .and. (.not. isopycnal(3,2,m))) then
     k = N3
     do j = S2p, N2p
        do i = S1p, N1p
           sed3U(i,j,m) = cc3(b3L,k,m)*conc(i,j,k+b3L,m)
!pgi$ unroll = n:8
           do kk = b3L+1, b3U
              sed3U(i,j,m) = sed3U(i,j,m) + cc3(kk,k,m)*conc(i,j,k+kk,m)
           end do
           sed3U(i,j,m) = sed3U(i,j,m)*(us_vec(3,m) + work3(i,j,k))
        end do
     end do
  end if
  !===========================================================================================================
  
  
  end subroutine sediment_bc
  
  
  
  
  
  
  
  
  
  
  
  subroutine filter(fil_dir,dir,m,phi,fil)
  
  implicit none
  
  integer, intent(in   ) ::  fil_dir
  integer, intent(in   ) ::  dir
  integer, intent(in   ) ::  m
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(  out) ::  fil(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Prinzipiell muss das gesamte Feld gefiltert werden; lediglich Dirichlet-Randbedingungen   !
  !                können/sollten unverändert bleiben. Wegen des versetzten Gitters ist das jedoch bei Wand- !
  !                normalen Geschwindigkeitskomponenten nicht möglich und könnte nur durch eine vor- und     !
  !                nachgeschaltete Interpolation auf das Druckgitter erreicht werden.                        !
  !              - Um jede Art der Randbedingungen zuzulassen, wird zunächst auf dem gesamten Feld gefiltert !
  !                und erst anschliessend mit Dirichlet-Randbedingungen überschrieben. Prinzipiell könnten   !
  !                aber auch diese Randbedingungen gefiltert werden.                                         !
  !----------------------------------------------------------------------------------------------------------!
  
  
  call exchange(fil_dir,dir,phi)
  
  
  !===========================================================================================================
  !=== Druck (NICHT Konzentrationen) =========================================================================
  !===========================================================================================================
  ! Anmerkung: Nur der Vollstaendigkeit halber aufgefuehrt, Druck sollte eigentlich nicht gefiltert werden.
  if (dir == 0 .and. m == 0) then
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 1) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 2) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 3) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     !   - entfallen, da der Druck ohne Dirichlet-RB auskommen sollte
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Konzentrationen =======================================================================================
  !===========================================================================================================
  if (dir == 0 .and. m > 0) then
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 1) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFc1(b1L,i,m)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFc1(ii,i,m)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 2) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFc2(b2L,j,m)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFc2(jj,j,m)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 3) then
        
        do k = S3p, N3p
           do j = S2p, N2p
              do i = S1p, N1p
                 fil(i,j,k) = cFc3(b3L,k,m)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFc3(kk,k,m)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     ! TEST!!! (siehe Geschwindigkeiten ...)
     if (filter_BC_yes) then ! TEST!!! Name irreführend!
     if (BCc_1L(m) == 1) fil(1 ,S2p:N2p,S3p:N3p) = phi(1 ,S2p:N2p,S3p:N3p)
     if (BCc_1U(m) == 1) fil(N1,S2p:N2p,S3p:N3p) = phi(N1,S2p:N2p,S3p:N3p)
     
     if (BCc_2L(m) == 1) fil(S1p:N1p,1 ,S3p:N3p) = phi(S1p:N1p,1 ,S3p:N3p)
     if (BCc_2U(m) == 1) fil(S1p:N1p,N2,S3p:N3p) = phi(S1p:N1p,N2,S3p:N3p)
     
     if (BCc_3L(m) == 1) fil(S1p:N1p,S2p:N2p,1 ) = phi(S1p:N1p,S2p:N2p,1 )
     if (BCc_3U(m) == 1) fil(S1p:N1p,S2p:N2p,N3) = phi(S1p:N1p,S2p:N2p,N3)
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Geschwindigkeiten =====================================================================================
  !===========================================================================================================
  if (dir == 1) then
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 1) then
        
        do k = S31B, N31B ! TEST!!! Grenzen sind nun zu eng, koennen gegen S31, N31, etc. ausgetauscht werden ...
           do j = S21B, N21B
              do i = S11B, N11B
                 fil(i,j,k) = cFu1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFu1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 2) then
        
        do k = S31B, N31B
           do j = S21B, N21B
              do i = S11B, N11B
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 3) then
        
        do k = S31B, N31B
           do j = S21B, N21B
              do i = S11B, N11B
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     ! TEST!!! BC_XX == 2,3 fehlen noch ...
     ! TEST!!! phi durch bc12, bc13 ersetzen?
     ! TEST!!! (Hat sich nicht bewährt. Insbesondere die Extrpolation macht Probleme ...)
     if (filter_BC_yes) then
     if (BC_2L == 1) fil(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     if (BC_2U == 1) fil(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     
     if (BC_3L == 1) fil(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     if (BC_3U == 1) fil(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     
     if (BC_1L > 0) then ! TEST!!! Funktionert NICHT bei Flussmündung und funktioniert bei Lock-Exchange ...
        fil(0 ,S21B:N21B,S31B:N31B) = bc11(S21B:N21B,S31B:N31B,1)
        call bc_extrapolation(dir,fil)
     end if
     if (BC_1U > 0) then
        fil(N1,S21B:N21B,S31B:N31B) = bc11(S21B:N21B,S31B:N31B,2)
        call bc_extrapolation(dir,fil)
     end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 2) then
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 1) then
        
        do k = S32B, N32B
           do j = S22B, N22B
              do i = S12B, N12B
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 2) then
        
        do k = S32B, N32B
           do j = S22B, N22B
              do i = S12B, N12B
                 fil(i,j,k) = cFv2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFv2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 3) then
        
        do k = S32B, N32B
           do j = S22B, N22B
              do i = S12B, N12B
                 fil(i,j,k) = cFp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFp3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (filter_BC_yes) then
     if (BC_1L == 1) fil(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     if (BC_1U == 1) fil(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     
     if (BC_3L == 1) fil(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     if (BC_3U == 1) fil(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     
     if (BC_2L > 0) then
        fil(S12B:N12B,0 ,S32B:N32B) = bc22(S12B:N12B,S32B:N32B,1)
        call bc_extrapolation(dir,fil)
     end if
     if (BC_2U > 0) then
        fil(S12B:N12B,N2,S32B:N32B) = bc22(S12B:N12B,S32B:N32B,2)
        call bc_extrapolation(dir,fil)
     end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  if (dir == 3) then
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 1) then
        
        do k = S33B, N33B
           do j = S23B, N23B
              do i = S13B, N13B
                 fil(i,j,k) = cFp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
                 do ii = b1L+1, b1U
                    fil(i,j,k) = fil(i,j,k) + cFp1(ii,i)*phi(i+ii,j,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 2) then
        
        do k = S33B, N33B
           do j = S23B, N23B
              do i = S13B, N13B
                 fil(i,j,k) = cFp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
                 do jj = b2L+1, b2U
                    fil(i,j,k) = fil(i,j,k) + cFp2(jj,j)*phi(i,j+jj,k)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     if (fil_dir == 3) then
        
        do k = S33B, N33B
           do j = S23B, N23B
              do i = S13B, N13B
                 fil(i,j,k) = cFw3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
                 do kk = b3L+1, b3U
                    fil(i,j,k) = fil(i,j,k) + cFw3(kk,k)*phi(i,j,k+kk)
                 end do
              end do
           end do
        end do
        
     end if
     !--------------------------------------------------------------------------------------------------------
     
     
     !--- Randbedingungen ------------------------------------------------------------------------------------
     if (filter_BC_yes) then
     if (BC_1L == 1) fil(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     if (BC_1U == 1) fil(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     
     if (BC_2L == 1) fil(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     if (BC_2U == 1) fil(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     
     if (BC_3L > 0) then
        fil(S13B:N13B,S23B:N23B,0 ) = bc33(S13B:N13B,S23B:N23B,1)
        call bc_extrapolation(dir,fil)
     end if
     if (BC_3U > 0) then
        fil(S13B:N13B,S23B:N23B,N3) = bc33(S13B:N13B,S23B:N23B,2)
        call bc_extrapolation(dir,fil)
     end if
     end if
     !--------------------------------------------------------------------------------------------------------
  end if
  !===========================================================================================================
  
  
  end subroutine filter
  
  
  
  
  
  
  !> \brief extrapolates Direchlet-BC to pressure
  !! \param[in] m direction
  !! \param[inout] phi pressure
  !! \attention Wird in "product_div_grad" und "explicit" verwendet almoste ver after grad
  subroutine bc_extrapolation(m,phi) bind (c,name='OP_bc_extrapolation')
  
  implicit none
  
  integer(c_int), intent(in   ) ::  m
  real(c_double), intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  if (m == 1) then
     if (BC_1L > 0) then
        i = 0
        do k = S31B, N31B
           do j = S21B, N21B
!pgi$ unroll = n:8
              do ii = 0, d1U
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,1)*phi(1+ii,j,k)
              end do
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           end do
        end do
     end if
     if (BC_1U > 0) then
        i = N1
        do k = S31B, N31B
           do j = S21B, N21B
!pgi$ unroll = n:8
              do ii = d1L, -1
                 phi(i,j,k) = phi(i,j,k) - cIup(ii,i)*phi(i+ii,j,k)
              end do
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (BC_2L > 0) then
        j = 0
        do k = S32B, N32B
           do i = S12B, N12B
!pgi$ unroll = n:8
              do jj = 0, d2U
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,1)*phi(i,1+jj,k)
              end do
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           end do
        end do
     end if
     if (BC_2U > 0) then
        j = N2
        do k = S32B, N32B
           do i = S12B, N12B
!pgi$ unroll = n:8
              do jj = d2L, -1
                 phi(i,j,k) = phi(i,j,k) - cIvp(jj,j)*phi(i,j+jj,k)
              end do
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (BC_3L > 0) then
        k = 0
        do j = S23B, N23B
           do i = S13B, N13B
!pgi$ unroll = n:8
              do kk = 0, d3U
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,1)*phi(i,j,1+kk)
              end do
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           end do
        end do
     end if
     if (BC_3U > 0) then
        k = N3
        do j = S23B, N23B
           do i = S13B, N13B
!pgi$ unroll = n:8
              do kk = d3L, -1
                 phi(i,j,k) = phi(i,j,k) - cIwp(kk,k)*phi(i,j,k+kk)
              end do
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  
  
  end subroutine bc_extrapolation
  
  
  
  
  
  
  
  
  
  
  
  subroutine bc_extrapolation_transp(m,phi)
  
  implicit none
  
  integer, intent(in   ) ::  m
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Achtung: - Wird in "product_div_grad" und "explicit" verwendet.                                          !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! bislang nur Dirichlet-RB eingebaut!
  !-----------------------------------------------------------------------------------------------------------
  if (m == 1) then
     if (BC_1L > 0) then
        i = 0
        do k = S31B, N31B
           do j = S21B, N21B
!pgi$ unroll = n:8
              do ii = 0, d1U
                 phi(i+ii+1,j,k) = phi(i+ii+1,j,k) - phi(i,j,k)*cIup(ii,1)/cIup(-1,1)
              end do
              phi(i,j,k) = phi(i,j,k) / cIup(-1,1)
           end do
        end do
     end if
     if (BC_1U > 0) then
        i = N1
        do k = S31B, N31B
           do j = S21B, N21B
!pgi$ unroll = n:8
              do ii = d1L, -1
                 phi(i+ii,j,k) = phi(i+ii,j,k) - phi(i,j,k)*cIup(ii,i)/cIup(0,i)
              end do
              phi(i,j,k) = phi(i,j,k) / cIup(0,i)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 2) then
     if (BC_2L > 0) then
        j = 0
        do k = S32B, N32B
           do i = S12B, N12B
!pgi$ unroll = n:8
              do jj = 0, d2U
                 phi(i,j+jj+1,k) = phi(i,j+jj+1,k) - phi(i,j,k)*cIvp(jj,1)/cIvp(-1,1)
              end do
              phi(i,j,k) = phi(i,j,k) / cIvp(-1,1)
           end do
        end do
     end if
     if (BC_2U > 0) then
        j = N2
        do k = S32B, N32B
           do i = S12B, N12B
!pgi$ unroll = n:8
              do jj = d2L, -1
                 phi(i,j+jj,k) = phi(i,j+jj,k) - phi(i,j,k)*cIvp(jj,j)/cIvp(0,j)
              end do
              phi(i,j,k) = phi(i,j,k) / cIvp(0,j)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (m == 3) then
     if (BC_3L > 0) then
        k = 0
        do j = S23B, N23B
           do i = S13B, N13B
!pgi$ unroll = n:8
              do kk = 0, d3U
                 phi(i,j,k+kk+1) = phi(i,j,k+kk+1) - phi(i,j,k)*cIwp(kk,1)/cIwp(-1,1)
              end do
              phi(i,j,k) = phi(i,j,k) / cIwp(-1,1)
           end do
        end do
     end if
     if (BC_3U > 0) then
        k = N3
        do j = S23B, N23B
           do i = S13B, N13B
!pgi$ unroll = n:8
              do kk = d3L, -1
                 phi(i,j,k+kk) = phi(i,j,k+kk) - phi(i,j,k)*cIwp(kk,k)/cIwp(0,k)
              end do
              phi(i,j,k) = phi(i,j,k) / cIwp(0,k)
           end do
        end do
     end if
  end if
  !-----------------------------------------------------------------------------------------------------------
  
  
  end subroutine bc_extrapolation_transp
  
  
  
end module cmod_operator
