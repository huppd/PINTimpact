!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

module mod_helmholtz
  
  
  use mod_dims
  use mod_vars
  use mod_diff
  use mod_exchange
  
  
  private
  
  public product_Helmholtz, product_Helmholtz_precond, product_Helmholtz_conc
  public product_Helmholtz_relax, product_Helmholtz_relax_coarse, product_Helmholtz_conc_relax
  public relaxation_Helmholtz, relaxation_Helmholtz_coarse, relaxation_Helmholtz_conc
  public handle_corner_conc
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  ! TEST!!! Generell wäre es sinnvoll, für alle Geschwindigkeiten und den Druck eigene Stencils zu speichern,
  !         um die Programmierung übersichtlicher zu halten und die IF-Abfragen zu vermeiden.
  
  
  !> \brief computes \f$ \mathtt{ Hel = (1 - c_m \theta_L \Delta t \mathcal L) phi} \f$
  subroutine product_Helmholtz(phi,Hel)
  
  implicit none
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Routine "Helmholtz" (ohne RB) ist ausgelagert, da sie auch aus RHS heraus aufgerufen wird,!
  !                wo keine RB notwendig sind. Sie werden erst anschliessend aufgeprägt.                     !
  !              - Austausch und Null-Setzen am Rand wird in Helmholtz-Routine erledigt.                     !
  !              - Interpolation in Wand-normaler Richting muss bei i=0 etc. um 1 versetzt werden!           !
  !              - Randbedingungen müssen auch in Ecken und Kanten gerechnet werden, daher die Intervalle    !
  !                S11B:N11B & Co..                                                                          !
  !              - Reihenfolge der Randbedingungen muss jeweils angepasst werden, da sich Ecken und Kanten   !
  !                überschneiden! Für die lid-driven cavity sind die tangentialen Richtungen zuletzt aufzu-  !
  !                zuprägen.                                                                                 !
  !              - Randbedingungen können NICHT in die Stencils mit eingebaut werden, da hier der Laplace-   !
  !                Operator gebildet wird (anders als z.B. bei der Divergenz).                               !
  !----------------------------------------------------------------------------------------------------------!
  
  
#ifdef NONBOUSSINESQ
  if (nonBoussinesq_yes) then
     if (direction == 1) then
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 Hel(i,j,k) = phi(i,j,k)*dens(i,j,k,1)
              end do
           end do
        end do
     end if
     if (direction == 2) then
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 Hel(i,j,k) = phi(i,j,k)*dens(i,j,k,2)
              end do
           end do
        end do
     end if
     if (direction == 3) then
        do k = S33, N33
           do j = S23, N23
!pgi$ unroll = n:8
              do i = S13, N13
                 Hel(i,j,k) = phi(i,j,k)*dens(i,j,k,3)
              end do
           end do
        end do
     end if
#else
  if (thetaL == 0.) then
     if (direction == 1) Hel(S11:N11,S21:N21,S31:N31) = phi(S11:N11,S21:N21,S31:N31)
     if (direction == 2) Hel(S12:N12,S22:N22,S32:N32) = phi(S12:N12,S22:N22,S32:N32)
     if (direction == 3) Hel(S13:N13,S23:N23,S33:N33) = phi(S13:N13,S23:N23,S33:N33)
#endif
  else
     call Helmholtz(direction,.true.,phi,Hel)
  end if
  
  !===========================================================================================================
  if (direction == 1) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S31B, N31B
           do i = S11B, N11B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              do jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S31B, N31B
           do i = S11B, N11B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              do jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S21B, N21B
           do i = S11B, N11B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              do kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S21B, N21B
           do i = S11B, N11B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              do kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 2) then
        i = 0
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = cIup(d1L,1)*phi(1+d1L,j,k)
!pgi$ unroll = n:8
              do ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cIup(ii,1)*phi(1+ii,j,k)
              end do
           end do
        end do
     end if
     if (BC_1U == 1 .or. BC_1U == 2) then
        i = N1
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = cIup(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
              do ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cIup(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 3 .or. BC_1L == 4) then
        i = 0
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = cDu1(d1L,i)*phi(1+d1L,j,k)
!pgi$ unroll = n:8
              do ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cDu1(ii,i)*phi(1+ii,j,k)
              end do
           end do
        end do
     end if
     if (BC_1U == 3 .or. BC_1U == 4) then
        i = N1
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = cDu1(d1L,i)*phi(i+d1L,j,k)
!pgi$ unroll = n:8
              do ii = d1L+1, d1U
                 Hel(i,j,k) = Hel(i,j,k) + cDu1(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (direction == 2) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S32B, N32B
           do j = S22B, N22B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              do ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S32B, N32B
           do j = S22B, N22B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              do ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S22B, N22B
           do i = S12B, N12B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              do kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S22B, N22B
           do i = S12B, N12B
              Hel(i,j,k) = cp3(b3L,k)*phi(i,j,k+b3L)
!pgi$ unroll = n:8
              do kk = b3L+1, b3U
                 Hel(i,j,k) = Hel(i,j,k) + cp3(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 2) then
        j = 0
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = cIvp(d2L,1)*phi(i,1+d2L,k)
!pgi$ unroll = n:8
              do jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cIvp(jj,1)*phi(i,1+jj,k)
              end do
           end do
        end do
     end if
     if (BC_2U == 1 .or. BC_2U == 2) then
        j = N2
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = cIvp(d2L,j)*phi(i,j+d2L,k)
!pgi$ unroll = n:8
              do jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cIvp(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 3 .or. BC_2L == 4) then
        j = 0
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = cDv2(d2L,1)*phi(i,1+d2L,k)
!pgi$ unroll = n:8
              do jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cDv2(jj,1)*phi(i,1+jj,k)
              end do
           end do
        end do
     end if
     if (BC_2U == 3 .or. BC_2U == 4) then
        j = N2
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = cDv2(d2L,j)*phi(i,j+d2L,k)
!pgi$ unroll = n:8
              do jj = d2L+1, d2U
                 Hel(i,j,k) = Hel(i,j,k) + cDv2(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (direction == 3) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S33B, N33B
           do j = S23B, N23B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              do ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S33B, N33B
           do j = S23B, N23B
              Hel(i,j,k) = cp1(b1L,i)*phi(i+b1L,j,k)
!pgi$ unroll = n:8
              do ii = b1L+1, b1U
                 Hel(i,j,k) = Hel(i,j,k) + cp1(ii,i)*phi(i+ii,j,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S33B, N33B
           do i = S13B, N13B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              do jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S33B, N33B
           do i = S13B, N13B
              Hel(i,j,k) = cp2(b2L,j)*phi(i,j+b2L,k)
!pgi$ unroll = n:8
              do jj = b2L+1, b2U
                 Hel(i,j,k) = Hel(i,j,k) + cp2(jj,j)*phi(i,j+jj,k)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 2) then
        k = 0
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = cIwp(d3L,1)*phi(i,j,1+d3L)
!pgi$ unroll = n:8
              do kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cIwp(kk,1)*phi(i,j,1+kk)
              end do
           end do
        end do
     end if
     if (BC_3U == 1 .or. BC_3U == 2) then
        k = N3
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = cIwp(d3L,k)*phi(i,j,k+d3L)
!pgi$ unroll = n:8
              do kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cIwp(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 3 .or. BC_3L == 4) then
        k = 0
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = cDw3(d3L,1)*phi(i,j,1+d3L)
!pgi$ unroll = n:8
              do kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cDw3(kk,1)*phi(i,j,1+kk)
              end do
           end do
        end do
     end if
     if (BC_3U == 3 .or. BC_3U == 4) then
        k = N3
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = cDw3(d3L,k)*phi(i,j,k+d3L)
!pgi$ unroll = n:8
              do kk = d3L+1, d3U
                 Hel(i,j,k) = Hel(i,j,k) + cDw3(kk,k)*phi(i,j,k+kk)
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  
  
  end subroutine product_Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  subroutine product_Helmholtz_precond(phi,Hel) ! TEST!!! Kann weg nach Ende der Diss.!
  
  implicit none
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Siehe Anmerkungen zu Routine "product_Helmholtz"!                                         !
  !              - F�r Vorkonditionierer DH��G ~= (DR�1G)(DHpG)��(DR�1G), wobei R die Randbedingungen und    !
  !                sonst nur "1" auf der Hauptdiagonalen enth�lt, d.h. es soll Hp ~= R�� H gelten. Diese An- !
  !                nahme ist gut, falls der Laplace-Anteil in H gering ist.                                  !
  !              - Da phi hier immer dem Gradienten entspricht, kann direkt phi = 0. auf dem Rand voraus-    !
  !                gesetzt werden (siehe unten).                                                             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (thetaL /= 0.) then
     call Helmholtz(direction,.true.,phi,Hel)
  else
     ! Dieser Fall sollte besser mit dem einfachen DR�1G-Vorkonditionierer behandelt werden:
     if (direction == 1) Hel(S11:N11,S21:N21,S31:N31) = phi(S11:N11,S21:N21,S31:N31)
     if (direction == 2) Hel(S12:N12,S22:N22,S32:N32) = phi(S12:N12,S22:N22,S32:N32)
     if (direction == 3) Hel(S13:N13,S23:N23,S33:N33) = phi(S13:N13,S23:N23,S33:N33)
  end if
  
  !===========================================================================================================
  if (direction == 1) then
     
     if (BC_1L > 0) then
        i = 0
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do ii = 0, d1U
                 Hel(i,j,k) = Hel(i,j,k) - cIup(ii,1)*Hel(1+ii,j,k)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIup(-1,1)
           end do
        end do
     end if
     if (BC_1U > 0) then
        i = N1
        do k = S31B, N31B
           do j = S21B, N21B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do ii = d1L, -1
                 Hel(i,j,k) = Hel(i,j,k) - cIup(ii,i)*Hel(i+ii,j,k)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIup(0,i)
           end do
        end do
     end if
     
     if (BC_2L > 0) Hel(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     if (BC_2U > 0) Hel(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     
     if (BC_3L > 0) Hel(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     if (BC_3U > 0) Hel(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     
  end if
  !===========================================================================================================
  if (direction == 2) then
     
     if (BC_2L > 0) then
        j = 0
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do jj = 0, d2U
                 Hel(i,j,k) = Hel(i,j,k) - cIvp(jj,1)*Hel(i,1+jj,k)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIvp(-1,1)
           end do
        end do
     end if
     if (BC_2U > 0) then
        j = N2
        do k = S32B, N32B
           do i = S12B, N12B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do jj = d2L, -1
                 Hel(i,j,k) = Hel(i,j,k) - cIvp(jj,j)*Hel(i,j+jj,k)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIvp(0,j)
           end do
        end do
     end if
     
     if (BC_1L > 0) Hel(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     if (BC_1U > 0) Hel(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     
     if (BC_3L > 0) Hel(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     if (BC_3U > 0) Hel(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     
  end if
  !===========================================================================================================
  if (direction == 3) then
     
     if (BC_3L > 0) then
        k = 0
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do kk = 0, d3U
                 Hel(i,j,k) = Hel(i,j,k) - cIwp(kk,1)*Hel(i,j,1+kk)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIwp(-1,1)
           end do
        end do
     end if
     if (BC_3U > 0) then
        k = N3
        do j = S23B, N23B
           do i = S13B, N13B
              Hel(i,j,k) = 0. ! == phi(i,j,k)
!pgi$ unroll = n:8
              do kk = d3L, -1
                 Hel(i,j,k) = Hel(i,j,k) - cIwp(kk,k)*Hel(i,j,k+kk)
              end do
              Hel(i,j,k) = Hel(i,j,k) / cIwp(0,k)
           end do
        end do
     end if
     
     if (BC_1L > 0) Hel(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     if (BC_1U > 0) Hel(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     
     if (BC_2L > 0) Hel(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     if (BC_2U > 0) Hel(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     
  end if
  !===========================================================================================================
  
  
  end subroutine product_Helmholtz_precond
  
  
  
  
  
  
  
  
  
  
  
  subroutine product_Helmholtz_conc(phi,Hel) ! TEST!!! Randbedingungen sollten noch in Stencil integriert werden!!! (ggf. auch bei Geschwindigkeiten/Druck?)
  
  implicit none
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Feld ==================================================================================================
  !===========================================================================================================
  call Helmholtz_conc(.true.,phi,Hel)
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  if (BCc_1L(conc_nu) == 1) Hel(1 ,S2p:N2p,S3p:N3p) = phi(1 ,S2p:N2p,S3p:N3p)
  if (BCc_1U(conc_nu) == 1) Hel(N1,S2p:N2p,S3p:N3p) = phi(N1,S2p:N2p,S3p:N3p)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_1L(conc_nu) == 3) then
     i = 1
     do k = S3p, N3p
        do j = S2p, N2p
           Hel(i,j,k) = - velReSc1(j,k,1,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do ii = b1L, b1U
              Hel(i,j,k) = Hel(i,j,k) + cc1(ii,i,conc_nu)*phi(i+ii,j,k)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_1L
        end do
     end do
  end if
  if (BCc_1U(conc_nu) == 3) then
     i = N1
     do k = S3p, N3p
        do j = S2p, N2p
           Hel(i,j,k) = - velReSc1(j,k,2,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do ii = b1L, b1U
              Hel(i,j,k) = Hel(i,j,k) + cc1(ii,i,conc_nu)*phi(i+ii,j,k)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_1U
        end do
     end do
  end if
  !===========================================================================================================
  if (BCc_2L(conc_nu) == 1) Hel(S1p:N1p,1 ,S3p:N3p) = phi(S1p:N1p,1 ,S3p:N3p)
  if (BCc_2U(conc_nu) == 1) Hel(S1p:N1p,N2,S3p:N3p) = phi(S1p:N1p,N2,S3p:N3p)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_2L(conc_nu) == 3) then
     j = 1
     do k = S3p, N3p
        do i = S1p, N1p
           Hel(i,j,k) = - velReSc2(i,k,1,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do jj = b2L, b2U
              Hel(i,j,k) = Hel(i,j,k) + cc2(jj,j,conc_nu)*phi(i,j+jj,k)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_2L
        end do
     end do
  end if
  
  if (BCc_2U(conc_nu) == 3) then
     j = N2
     do k = S3p, N3p
        do i = S1p, N1p
           Hel(i,j,k) = - velReSc2(i,k,2,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do jj = b2L, b2U
              Hel(i,j,k) = Hel(i,j,k) + cc2(jj,j,conc_nu)*phi(i,j+jj,k)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_2U
        end do
     end do
  end if
  !===========================================================================================================
  if (BCc_3L(conc_nu) == 1) Hel(S1p:N1p,S2p:N2p,1 ) = phi(S1p:N1p,S2p:N2p,1 )
  if (BCc_3U(conc_nu) == 1) Hel(S1p:N1p,S2p:N2p,N3) = phi(S1p:N1p,S2p:N2p,N3)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_3L(conc_nu) == 3) then
     k = 1
     do j = S2p, N2p
        do i = S1p, N1p
           Hel(i,j,k) = - velReSc3(i,j,1,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do kk = b3L, b3U
              Hel(i,j,k) = Hel(i,j,k) + cc3(kk,k,conc_nu)*phi(i,j,k+kk)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_3L
        end do
     end do
  end if
  
  if (BCc_3U(conc_nu) == 3) then
     k = N3
     do j = S2p, N2p
        do i = S1p, N1p
           Hel(i,j,k) = - velReSc3(i,j,2,1)*phi(i,j,k)
!pgi$ unroll = n:8
           do kk = b3L, b3U
              Hel(i,j,k) = Hel(i,j,k) + cc3(kk,k,conc_nu)*phi(i,j,k+kk)
           end do
           Hel(i,j,k) = Hel(i,j,k) * dx_3U
        end do
     end do
  end if
  !===========================================================================================================
  
  
  if (corner_yes) call handle_corner_conc(N1,N2,N3,Hel)
  
  
  end subroutine product_Helmholtz_conc






  
  
  
  
  
  subroutine product_Helmholtz_relax(phi,Hel) ! TEST!!! 2D-Variante fehlt noch ...
  
  implicit none
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(out  ) ::  Hel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, j, k, g
  
  
  !-----------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!        !
  !            - Reihenfolge der Randbedingungen analog zu product_Helmholtz!                           !
  !-----------------------------------------------------------------------------------------------------!
  
  ! Alternative: hier direkt auf Druckgitter differenzieren und gleiche Restriktion wie f�r Druck verwenden!
  
  g = 1
  
  call exchange_relax(g,0,0,0,direction,.true.,phi)
  
  !===========================================================================================================
  if (direction == 1) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                         &  +  cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)   &
                         &  + (cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g))*phi(i,j,k))
              end do
           end do
        end do
     else
        do k = S31, N31
           do j = S21, N21
!pgi$ unroll = n:8
              do i = S11, N11
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                         &  + (cu11R(0,i,g) + cp22R(0,j,g))*phi(i,j,k))
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S11B:N11B,1 ,S31B:N31B) = phi(S11B:N11B,1 ,S31B:N31B)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S11B:N11B,N2,S31B:N31B) = phi(S11B:N11B,N2,S31B:N31B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S31B, N31B
!pgi$ unroll = n:8
           do i = S11B, N11B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S31B, N31B
!pgi$ unroll = n:8
           do i = S11B, N11B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S11B:N11B,S21B:N21B,1 ) = phi(S11B:N11B,S21B:N21B,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S11B:N11B,S21B:N21B,N3) = phi(S11B:N11B,S21B:N21B,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S21B, N21B
!pgi$ unroll = n:8
           do i = S11B, N11B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S21B, N21B
!pgi$ unroll = n:8
           do i = S11B, N11B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L > 0) then
        i = 0
        do k = S31B, N31B
!pgi$ unroll = n:8
           do j = S21B, N21B
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     if (BC_1U > 0) then
        i = N1
        do k = S31B, N31B
!pgi$ unroll = n:8
           do j = S21B, N21B
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (direction == 2) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)   &
                         &  +  cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)   &
                         &  + (cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g))*phi(i,j,k))
              end do
           end do
        end do
     else
        do k = S32, N32
           do j = S22, N22
!pgi$ unroll = n:8
              do i = S12, N12
                 Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                         &  +  cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)   &
                         &  + (cp11R(0,i,g) + cv22R(0,j,g))*phi(i,j,k))
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S22B:N22B,S32B:N32B) = phi(1 ,S22B:N22B,S32B:N32B)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S22B:N22B,S32B:N32B) = phi(N1,S22B:N22B,S32B:N32B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S32B, N32B
!pgi$ unroll = n:8
           do j = S22B, N22B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S32B, N32B
!pgi$ unroll = n:8
           do j = S22B, N22B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S12B:N12B,S22B:N22B,1 ) = phi(S12B:N12B,S22B:N22B,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S12B:N12B,S22B:N22B,N3) = phi(S12B:N12B,S22B:N22B,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S22B, N22B
!pgi$ unroll = n:8
           do i = S12B, N12B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S22B, N22B
!pgi$ unroll = n:8
           do i = S12B, N12B
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L > 0) then
        j = 0
        do k = S32B, N32B
!pgi$ unroll = n:8
           do i = S12B, N12B
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     if (BC_2U > 0) then
        j = N2
        do k = S32B, N32B
!pgi$ unroll = n:8
           do i = S12B, N12B
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (direction == 3) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     ! dimens /= 3 durch direction /= 3 trivial erf�llt.
     do k = S33, N33
        do j = S23, N23
!pgi$ unroll = n:8
           do i = S13, N13
              Hel(i,j,k) =  phi(i,j,k) - multL*(                                     &
                      &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)   &
                      &  +  cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)   &
                      &  +  cw33R(-1,k,g)*phi(i,j,k-1) + cw33R(1,k,g)*phi(i,j,k+1)   &
                      &  + (cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g))*phi(i,j,k))
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S23B:N23B,S33B:N33B) = phi(1 ,S23B:N23B,S33B:N33B)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S23B:N23B,S33B:N33B) = phi(N1,S23B:N23B,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S33B, N33B
!pgi$ unroll = n:8
           do j = S23B, N23B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S33B, N33B
!pgi$ unroll = n:8
           do j = S23B, N23B
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S13B:N13B,1 ,S33B:N33B) = phi(S13B:N13B,1 ,S33B:N33B)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S13B:N13B,N2,S33B:N33B) = phi(S13B:N13B,N2,S33B:N33B)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S33B, N33B
!pgi$ unroll = n:8
           do i = S13B, N13B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S33B, N33B
!pgi$ unroll = n:8
           do i = S13B, N13B
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L > 0) then
        k = 0
        do j = S23B, N23B
!pgi$ unroll = n:8
           do i = S13B, N13B
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     if (BC_3U > 0) then
        k = N3
        do j = S23B, N23B
!pgi$ unroll = n:8
           do i = S13B, N13B
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  
  
  end subroutine product_Helmholtz_relax
  
  
  
  
  
  
  
  
  
  
  
  subroutine product_Helmholtz_relax_coarse(g,phi,Hel) ! TEST!!! 2D-Variante fehlt noch ...
  
  implicit none
  
  integer, intent(in   ) ::  g
  
  real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  real   , intent(out  ) ::  Hel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  integer                ::  i, N1, N1R, N11R
  integer                ::  j, N2, N2R, N22R
  integer                ::  k, N3, N3R, N33R
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!             !
  !----------------------------------------------------------------------------------------------------------!
  
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  N1R  = N1 + d1R ! TEST!!! Substituieren ...
  N2R  = N2 + d2R
  N3R  = N3 + d3R
  
  N11R = N1 + d11R
  N22R = N2 + d22R
  N33R = N3 + d33R
  
  
  call exchange_relax(g,0,0,0,0,.true.,phi)
  
  
  !===========================================================================================================
  if (direction == 1) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        do k = S33R, N33R
           do j = S22R, N22R
!pgi$ unroll = n:8
              do i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                         &   + cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)    &
                         &   + phi(i,j,k)*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
              end do
           end do
        end do
     else
        do k = S33R, N33R
           do j = S22R, N22R
!pgi$ unroll = n:8
              do i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cu11R(-1,i,g)*phi(i-1,j,k) + cu11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                         &   + phi(i,j,k)*(cu11R(0,i,g) + cp22R(0,j,g)))
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S1R:N1R,1 ,S3R:N3R) = phi(S1R:N1R,1 ,S3R:N3R)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S1R:N1R,N2,S3R:N3R) = phi(S1R:N1R,N2,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S1R:N1R,S2R:N2R,1 ) = phi(S1R:N1R,S2R:N2R,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S1R:N1R,S2R:N2R,N3) = phi(S1R:N1R,S2R:N2R,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L > 0) then
        i = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_1U > 0) then
        i = N1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cu11R(0,i,g)*phi(i,j,k) + cu11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  if (direction == 2) then
     
     !--------------------------------------------------------------------------------------------------------
     !--- Feld -----------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
        do k = S33R, N33R
           do j = S22R, N22R
!pgi$ unroll = n:8
              do i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)    &
                         &   + cp33R(-1,k,g)*phi(i,j,k-1) + cp33R(1,k,g)*phi(i,j,k+1)    &
                         &   + phi(i,j,k)*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
              end do
           end do
        end do
     else
        do k = S33R, N33R
           do j = S22R, N22R
!pgi$ unroll = n:8
              do i = S11R, N11R
                 Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                         &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                         &   + cv22R(-1,j,g)*phi(i,j-1,k) + cv22R(1,j,g)*phi(i,j+1,k)    &
                         &   + phi(i,j,k)*(cp11R(0,i,g) + cv22R(0,j,g)))
              end do
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S2R:N2R,S3R:N3R) = phi(1 ,S2R:N2R,S3R:N3R)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S2R:N2R,S3R:N3R) = phi(N1,S2R:N2R,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 1 .or. BC_3L == 3) Hel(S1R:N1R,S2R:N2R,1 ) = phi(S1R:N1R,S2R:N2R,1 )
     if (BC_3U == 1 .or. BC_3U == 3) Hel(S1R:N1R,S2R:N2R,N3) = phi(S1R:N1R,S2R:N2R,N3)
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L == 2 .or. BC_3L == 4) then
        k = 1
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        k = N3
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp33R(0,k,g)*phi(i,j,k) + cp33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L > 0) then
        j = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2U > 0) then
        j = N2
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cv22R(0,j,g)*phi(i,j,k) + cv22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  ! dimens /= 3 durch direction /= 3 trivial erf�llt.
  if (direction == 3) then
     do k = S33R, N33R
        do j = S22R, N22R
!pgi$ unroll = n:8
           do i = S11R, N11R
              Hel(i,j,k) = phi(i,j,k) - multL*(                                       &
                      &     cp11R(-1,i,g)*phi(i-1,j,k) + cp11R(1,i,g)*phi(i+1,j,k)    &
                      &   + cp22R(-1,j,g)*phi(i,j-1,k) + cp22R(1,j,g)*phi(i,j+1,k)    &
                      &   + cw33R(-1,k,g)*phi(i,j,k-1) + cw33R(1,k,g)*phi(i,j,k+1)    &
                      &   + phi(i,j,k)*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
           end do
        end do
     end do
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-tangentiale Geschwindigkeiten -----------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 1 .or. BC_1L == 3) Hel(1 ,S2R:N2R,S3R:N3R) = phi(1 ,S2R:N2R,S3R:N3R)
     if (BC_1U == 1 .or. BC_1U == 3) Hel(N1,S2R:N2R,S3R:N3R) = phi(N1,S2R:N2R,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     if (BC_1L == 2 .or. BC_1L == 4) then
        i = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R( 1,i,g)*phi(i+1,j,k)
           end do
        end do
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        i = N1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do j = S2R, N2R
              Hel(i,j,k) = cp11R(0,i,g)*phi(i,j,k) + cp11R(-1,i,g)*phi(i-1,j,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 1 .or. BC_2L == 3) Hel(S1R:N1R,1 ,S3R:N3R) = phi(S1R:N1R,1 ,S3R:N3R)
     if (BC_2U == 1 .or. BC_2U == 3) Hel(S1R:N1R,N2,S3R:N3R) = phi(S1R:N1R,N2,S3R:N3R)
     !--------------------------------------------------------------------------------------------------------
     if (BC_2L == 2 .or. BC_2L == 4) then
        j = 1
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R( 1,j,g)*phi(i,j+1,k)
           end do
        end do
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        j = N2
        do k = S3R, N3R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cp22R(0,j,g)*phi(i,j,k) + cp22R(-1,j,g)*phi(i,j-1,k)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     !--- Wand-normale Geschwindigkeit -----------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (BC_3L > 0) then
        k = 1
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R( 1,k,g)*phi(i,j,k+1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     if (BC_3U > 0) then
        k = N3
        do j = S2R, N2R
!pgi$ unroll = n:8
           do i = S1R, N1R
              Hel(i,j,k) = cw33R(0,k,g)*phi(i,j,k) + cw33R(-1,k,g)*phi(i,j,k-1)
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------------
     
  end if
  !===========================================================================================================
  
  
  end subroutine product_Helmholtz_relax_coarse
  
  
  
  
  
  
  
  
  
  
  
  subroutine product_Helmholtz_conc_relax(g,phi,Hel)
  
  implicit none
  
  integer, intent(in   ) ::  g
  
  real   , intent(inout) ::  phi(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  real   , intent(out  ) ::  Hel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  integer                ::  i, N1, S1Rc, S11Rc, N1Rc, N11Rc
  integer                ::  j, N2, S2Rc, S22Rc, N2Rc, N22Rc
  integer                ::  k, N3, S3Rc, S33Rc, N3Rc, N33Rc
  
  integer                ::  m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkung: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  m = conc_nu
  
  !===========================================================================================================
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  S1Rc  = S1R ! TEST!!! evtl. herausziehen nach init_limits!
  S2Rc  = S2R
  S3Rc  = S3R
  
  S11Rc = S11R
  S22Rc = S22R
  S33Rc = S33R
  
  N1Rc  = N1 + d1R
  N2Rc  = N2 + d2R
  N3Rc  = N3 + d3R
  
  N11Rc = N1 + d11R
  N22Rc = N2 + d22R
  N33Rc = N3 + d33R
  
  !-----------------------------------------------------------------------------------------------------------
  !--- Sonderf�lle (Symmetrie vs. Wand) ----------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_1L(m) > 0) then
     S1Rc  = 1
     S11Rc = 2
  end if
  if (BCc_1U(m) > 0) then
     N1Rc  = N1
     N11Rc = N1-1
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_2L(m) > 0) then
     S2Rc  = 1
     S22Rc = 2
  end if
  if (BCc_2U(m) > 0) then
     N2Rc  = N2
     N22Rc = N2-1
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_3L(m) > 0) then
     S3Rc  = 1
     S33Rc = 2
  end if
  if (BCc_3U(m) > 0) then
     N3Rc  = N3
     N33Rc = N3-1
  end if
  !===========================================================================================================
  
  
  call exchange_relax(g,0,0,0,0,.true.,phi)
  
  !===========================================================================================================
  !=== Feld ==================================================================================================
  !===========================================================================================================
  if (dimens == 3) then
     do k = S33Rc, N33Rc
        do j = S22Rc, N22Rc
!pgi$ unroll = n:8
           do i = S11Rc, N11Rc
              Hel(i,j,k) = phi(i,j,k) - multL*(                                             &
                      &     cc11R(-1,i,g,m)*phi(i-1,j,k) + cc11R(1,i,g,m)*phi(i+1,j,k)      &
                      &   + cc22R(-1,j,g,m)*phi(i,j-1,k) + cc22R(1,j,g,m)*phi(i,j+1,k)      &
                      &   + cc33R(-1,k,g,m)*phi(i,j,k-1) + cc33R(1,k,g,m)*phi(i,j,k+1)      &
                      &   + phi(i,j,k)*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m)))
           end do
        end do
     end do
  else
     do k = S33Rc, N33Rc
        do j = S22Rc, N22Rc
!pgi$ unroll = n:8
           do i = S11Rc, N11Rc
              Hel(i,j,k) = phi(i,j,k) - multL*(                                             &
                      &     cc11R(-1,i,g,m)*phi(i-1,j,k) + cc11R(1,i,g,m)*phi(i+1,j,k)      &
                      &   + cc22R(-1,j,g,m)*phi(i,j-1,k) + cc22R(1,j,g,m)*phi(i,j+1,k)      &
                      &   + phi(i,j,k)*(cc11R(0,i,g,m) + cc22R(0,j,g,m)))
           end do
        end do
     end do
  end if
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  if (BCc_1L(m) == 1) Hel(1 ,S2Rc:N2Rc,S3Rc:N3Rc) = phi(1 ,S2Rc:N2Rc,S3Rc:N3Rc)
  if (BCc_1U(m) == 1) Hel(N1,S2Rc:N2Rc,S3Rc:N3Rc) = phi(N1,S2Rc:N2Rc,S3Rc:N3Rc)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_1L(m) == 3) then
     i = 1
     do k = S3Rc, N3Rc
!pgi$ unroll = n:8
        do j = S2Rc, N2Rc
           Hel(i,j,k) = dx_1L * ((cc11R(0,i,g,m) - velReSc1(j,k,1,g))*phi(i,j,k) + cc11R(1,i,g,m)*phi(i+1,j,k))
        end do
     end do
  end if
  
  if (BCc_1U(m) == 3) then
     i = N1
     do k = S3Rc, N3Rc
!pgi$ unroll = n:8
        do j = S2Rc, N2Rc
           Hel(i,j,k) = dx_1U * ((cc11R(0,i,g,m) - velReSc1(j,k,2,g))*phi(i,j,k) + cc11R(-1,i,g,m)*phi(i-1,j,k))
        end do
     end do
  end if
  !===========================================================================================================
  if (BCc_2L(m) == 1) Hel(S1Rc:N1Rc,1 ,S3Rc:N3Rc) = phi(S1Rc:N1Rc,1 ,S3Rc:N3Rc)
  if (BCc_2U(m) == 1) Hel(S1Rc:N1Rc,N2,S3Rc:N3Rc) = phi(S1Rc:N1Rc,N2,S3Rc:N3Rc)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_2L(m) == 3) then
     j = 1
     do k = S3Rc, N3Rc
!pgi$ unroll = n:8
        do i = S1Rc, N1Rc
           Hel(i,j,k) = dx_2L * ((cc22R(0,j,g,m) - velReSc2(i,k,1,g))*phi(i,j,k) + cc22R(1,j,g,m)*phi(i,j+1,k))
        end do
     end do
  end if
  
  if (BCc_2U(m) == 3) then
     j = N2
     do k = S3Rc, N3Rc
!pgi$ unroll = n:8
        do i = S1Rc, N1Rc
           Hel(i,j,k) = dx_2U * ((cc22R(0,j,g,m) - velReSc2(i,k,2,g))*phi(i,j,k) + cc22R(-1,j,g,m)*phi(i,j-1,k))
        end do
     end do
  end if
  !===========================================================================================================
  if (BCc_3L(m) == 1) Hel(S1Rc:N1Rc,S2Rc:N2Rc,1 ) = phi(S1Rc:N1Rc,S2Rc:N2Rc,1 )
  if (BCc_3U(m) == 1) Hel(S1Rc:N1Rc,S2Rc:N2Rc,N3) = phi(S1Rc:N1Rc,S2Rc:N2Rc,N3)
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_3L(m) == 3) then
     k = 1
     do j = S2Rc, N2Rc
!pgi$ unroll = n:8
        do i = S1Rc, N1Rc
           Hel(i,j,k) = dx_3L * ((cc33R(0,k,g,m) - velReSc3(i,j,1,g))*phi(i,j,k) + cc33R(1,k,g,m)*phi(i,j,k+1))
        end do
     end do
  end if
  
  if (BCc_3U(m) == 3) then
     k = N3
     do j = S2Rc, N2Rc
!pgi$ unroll = n:8
        do i = S1Rc, N1Rc
           Hel(i,j,k) = dx_3U * ((cc33R(0,k,g,m) - velReSc3(i,j,2,g))*phi(i,j,k) + cc33R(-1,k,g,m)*phi(i,j,k-1))
        end do
     end do
  end if
  !===========================================================================================================
  
  
  if (corner_yes) call handle_corner_conc(N1,N2,N3,Hel)
  
  
  end subroutine product_Helmholtz_conc_relax
  
  
  
  
  
  
  
  
  
  
  
  subroutine relaxation_Helmholtz(init_yes,n_relax,bb,rel) ! TEST!!! 2D-Variante fehlt noch ...
  
  implicit none
  
  logical, intent(in   ) ::  init_yes
  integer, intent(in   ) ::  n_relax
  
  real   , intent(in   ) ::  bb (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , intent(inout) ::  rel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  integer                ::  i, j, k, g, r
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gewählt sind!           !
  !              - Randbedingungen müssen auch in Ecken und Kanten gerechnet werden, daher die Intervalle    !
  !                S11B:N11B & Co.                                                                           !
  !              - Initialisierung bezieht sich nur auf Feldbereich, daher die Intervalle S11:N11 usw..      !
  !              - Reihenfolge der Randbedingungen analog zu product_Helmholtz!                              !
  !              - Dirichlet-Randbedingungen sind zur Optimierung der Konvergenzrate z.T. VOR der Relaxation !
  !                eingefügt!                                                                                !
  !----------------------------------------------------------------------------------------------------------!
  
  g = 1
  
  ! Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  !IF (init_yes) rel = 0.
  
  
  do r = 1, n_relax
     
     if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,direction,.true.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 1) then
        
        !=====================================================================================================
        if (r == 1 .and. init_yes) then
           if (BC_1L <= 0) rel(S11-1  ,S21:N21,S31:N31) = 0.
           if (BC_2L <= 0) rel(S11:N11,S21-1  ,S31:N31) = 0.
           if (BC_3L <= 0) rel(S11:N11,S21:N21,S31-1  ) = 0.
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 1 .or. BC_2L == 3) rel(S11B:N11B,1 ,S31B:N31B) = bb(S11B:N11B,1 ,S31B:N31B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 2 .or. BC_2L == 4) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S31B, N31B
!pgi$ unroll = n:8
                 do i = S11B, N11B
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 end do
              end do
           else
              do k = S31B, N31B
!pgi$ unroll = n:8
                 do i = S11B, N11B
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 1 .or. BC_3L == 3) rel(S11B:N11B,S21B:N21B,1 ) = bb(S11B:N11B,S21B:N21B,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 2 .or. BC_3L == 4) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S21B, N21B
!pgi$ unroll = n:8
                 do i = S11B, N11B
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 end do
              end do
           else
              do j = S21B, N21B
!pgi$ unroll = n:8
                 do i = S11B, N11B
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) then
           i = 0
           if (r == 1 .and. init_yes) then
              do k = S31B, N31B
!pgi$ unroll = n:8
                 do j = S21B, N21B
                    rel(i,j,k) = bb(i,j,k) / cu11R(0,i,g)
                 end do
              end do
           else
              do k = S31B, N31B
!pgi$ unroll = n:8
                 do j = S21B, N21B
                    rel(i,j,k) = (bb(i,j,k) - cu11R( 1,i,g)*rel(i+1,j,k)) / cu11R(0,i,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cu11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
           
        else
           
           do k = S31, N31
              do j = S21, N21
!pgi$ unroll = n:8
                 do i = S11, N11
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cu11R(-1,i,g)*rel(i-1,j,k) + cu11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
           
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 1 .or. BC_2U == 3) rel(S11B:N11B,N2,S31B:N31B) = bb(S11B:N11B,N2,S31B:N31B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 2 .or. BC_2U == 4) then
           j = N2
           do k = S31B, N31B
!pgi$ unroll = n:8
              do i = S11B, N11B
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 1 .or. BC_3U == 3) rel(S11B:N11B,S21B:N21B,N3) = bb(S11B:N11B,S21B:N21B,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 2 .or. BC_3U == 4) then
           k = N3
           do j = S21B, N21B
!pgi$ unroll = n:8
              do i = S11B, N11B
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U > 0) then
           i = N1
           do k = S31B, N31B
!pgi$ unroll = n:8
              do j = S21B, N21B
                 rel(i,j,k) = (bb(i,j,k) - cu11R(-1,i,g)*rel(i-1,j,k)) / cu11R(0,i,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 2) then
        
        !=====================================================================================================
        if (r == 1 .and. init_yes) then
           if (BC_1L <= 0) rel(S12-1  ,S22:N22,S32:N32) = 0.
           if (BC_2L <= 0) rel(S12:N12,S22-1  ,S32:N32) = 0.
           if (BC_3L <= 0) rel(S12:N12,S22:N22,S32-1  ) = 0.
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 1 .or. BC_1L == 3) rel(1 ,S22B:N22B,S32B:N32B) = bb(1 ,S22B:N22B,S32B:N32B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 2 .or. BC_1L == 4) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S32B, N32B
!pgi$ unroll = n:8
                 do j = S22B, N22B
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 end do
              end do
           else
              do k = S32B, N32B
!pgi$ unroll = n:8
                 do j = S22B, N22B
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 1 .or. BC_3L == 3) rel(S12B:N12B,S22B:N22B,1 ) = bb(S12B:N12B,S22B:N22B,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 2 .or. BC_3L == 4) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S22B, N22B
!pgi$ unroll = n:8
                 do i = S12B, N12B
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 end do
              end do
           else
              do j = S22B, N22B
!pgi$ unroll = n:8
                 do i = S12B, N12B
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L > 0) then
           j = 0
           if (r == 1 .and. init_yes) then
              do k = S32B, N32B
!pgi$ unroll = n:8
                 do i = S12B, N12B
                    rel(i,j,k) = bb(i,j,k) / cv22R(0,j,g)
                 end do
              end do
           else
              do k = S32B, N32B
!pgi$ unroll = n:8
                 do i = S12B, N12B
                    rel(i,j,k) = (bb(i,j,k) - cv22R( 1,j,g)*rel(i,j+1,k)) / cv22R(0,j,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cv22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
           
        else
           
           do k = S32, N32
              do j = S22, N22
!pgi$ unroll = n:8
                 do i = S12, N12
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cv22R(-1,j,g)*rel(i,j-1,k) + cv22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
           
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 1 .or. BC_1U == 3) rel(N1,S22B:N22B,S32B:N32B) = bb(N1,S22B:N22B,S32B:N32B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 2 .or. BC_1U == 4) then
           i = N1
           do k = S32B, N32B
!pgi$ unroll = n:8
              do j = S22B, N22B
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 1 .or. BC_3U == 3) rel(S12B:N12B,S22B:N22B,N3) = bb(S12B:N12B,S22B:N22B,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 2 .or. BC_3U == 4) then
           k = N3
           do j = S22B, N22B
!pgi$ unroll = n:8
              do i = S12B, N12B
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U > 0) then
           j = N2
           do k = S32, N32
!pgi$ unroll = n:8
              do i = S12, N12
                 rel(i,j,k) = (bb(i,j,k) - cv22R(-1,j,g)*rel(i,j-1,k)) / cv22R(0,j,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 3) then
        
        !=====================================================================================================
        if (r == 1 .and. init_yes) then
           if (BC_1L <= 0) rel(S13-1  ,S23:N23,S33:N33) = 0.
           if (BC_2L <= 0) rel(S13:N13,S23-1  ,S33:N33) = 0.
           if (BC_3L <= 0) rel(S13:N13,S23:N23,S33-1  ) = 0.
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 1 .or. BC_1L == 3) rel(1 ,S23B:N23B,S33B:N33B) = bb(1 ,S23B:N23B,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 2 .or. BC_1L == 4) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S33B, N33B
!pgi$ unroll = n:8
                 do j = S23B, N23B
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 end do
              end do
           else
              do k = S33B, N33B
!pgi$ unroll = n:8
                 do j = S23B, N23B
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 1 .or. BC_2L == 3) rel(S13B:N13B,1 ,S33B:N33B) = bb(S13B:N13B,1 ,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 2 .or. BC_2L == 4) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S33B, N33B
!pgi$ unroll = n:8
                 do i = S13B, N13B
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 end do
              end do
           else
              do k = S33B, N33B
!pgi$ unroll = n:8
                 do i = S13B, N13B
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L > 0) then
           k = 0
           if (r == 1 .and. init_yes) then
              do j = S23B, N23B
!pgi$ unroll = n:8
                 do i = S13B, N13B
                    rel(i,j,k) = bb(i,j,k) / cw33R(0,k,g)
                 end do
              end do
           else
              do j = S23B, N23B
!pgi$ unroll = n:8
                 do i = S13B, N13B
                    rel(i,j,k) = (bb(i,j,k) - cw33R( 1,k,g)*rel(i,j,k+1)) / cw33R(0,k,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           
           do k = S33, N33
              do j = S23, N23
!pgi$ unroll = n:8
                 do i = S13, N13
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cw33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 end do
              end do
           end do
           
        else
           
           do k = S33, N33
              do j = S23, N23
!pgi$ unroll = n:8
                 do i = S13, N13
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cw33R(-1,k,g)*rel(i,j,k-1) + cw33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 end do
              end do
           end do
           
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 1 .or. BC_1U == 3) rel(N1,S23B:N23B,S33B:N33B) = bb(N1,S23B:N23B,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 2 .or. BC_1U == 4) then
           i = N1
           do k = S33B, N33B
!pgi$ unroll = n:8
              do j = S23B, N23B
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 1 .or. BC_2U == 3) rel(S13B:N13B,N2,S33B:N33B) = bb(S13B:N13B,N2,S33B:N33B)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 2 .or. BC_2U == 4) then
           j = N2
           do k = S33B, N33B
!pgi$ unroll = n:8
              do i = S13B, N13B
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U > 0) then
           k = N3
           do j = S23B, N23B
!pgi$ unroll = n:8
              do i = S13B, N13B
                 rel(i,j,k) = (bb(i,j,k) - cw33R(-1,k,g)*rel(i,j,k-1)) / cw33R(0,k,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  end do
  
  
  end subroutine relaxation_Helmholtz
  
  
  
  
  
  
  
  
  
  
  
  subroutine relaxation_Helmholtz_coarse(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ...
  
  implicit none
  
  logical, intent(in   ) ::  init_yes
  integer, intent(in   ) ::  n_relax
  
  integer, intent(in   ) ::  g
  
  real   , intent(in   ) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  real   , intent(inout) ::  rel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  integer                ::  i, N1, N1R, N11R
  integer                ::  j, N2, N2R, N22R
  integer                ::  k, N3, N3R, N33R
  
  integer                ::  r
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Dirichlet-Randbedingungen sind der Konvergenzrate wegen überall VOR der Relaxation        !
  !                eingefügt!                                                                                !
  !----------------------------------------------------------------------------------------------------------!
  
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  N1R  = N1 + d1R
  N2R  = N2 + d2R
  N3R  = N3 + d3R
  
  N11R = N1 + d11R
  N22R = N2 + d22R
  N33R = N3 + d33R
  
  
  ! Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  !IF (init_yes) rel = 0.
  
  
  if (init_yes) then
     if (BC_1L <= 0) rel(S11R-1   ,S22R:N22R,S33R:N33R) = 0.
     if (BC_2L <= 0) rel(S11R:N11R,S22R-1   ,S33R:N33R) = 0.
     if (BC_3L <= 0) rel(S11R:N11R,S22R:N22R,S33R-1   ) = 0.
  end if
  
  
  do r = 1, n_relax
     
     if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 1) then
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 1 .or. BC_2L == 3) rel(S1R:N1R,1 ,S3R:N3R) = bb(S1R:N1R,1 ,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 2 .or. BC_2L == 4) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 1 .or. BC_3L == 3) rel(S1R:N1R,S2R:N2R,1 ) = bb(S1R:N1R,S2R:N2R,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 2 .or. BC_3L == 4) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 end do
              end do
           else
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L > 0) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cu11R(0,i,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cu11R( 1,i,g)*rel(i+1,j,k)) / cu11R(0,i,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           do k = S33R, N33R
              do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cu11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
        else
           do k = S33R, N33R
             do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cu11R(-1,i,g)*rel(i-1,j,k) + cu11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cu11R(0,i,g) + cp22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 1 .or. BC_2U == 3) rel(S1R:N1R,N2,S3R:N3R) = bb(S1R:N1R,N2,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 2 .or. BC_2U == 4) then
           j = N2
           do k = S3R, N3R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 1 .or. BC_3U == 3) rel(S1R:N1R,S2R:N2R,N3) = bb(S1R:N1R,S2R:N2R,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 2 .or. BC_3U == 4) then
           k = N3
           do j = S2R, N2R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U > 0) then
           i = N1
           do k = S3R, N3R
!pgi$ unroll = n:8
              do j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cu11R(-1,i,g)*rel(i-1,j,k)) / cu11R(0,i,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 2) then
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 1 .or. BC_1L == 3) rel(1 ,S2R:N2R,S3R:N3R) = bb(1 ,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 2 .or. BC_1L == 4) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 1 .or. BC_3L == 3) rel(S1R:N1R,S2R:N2R,1 ) = bb(S1R:N1R,S2R:N2R,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L == 2 .or. BC_3L == 4) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp33R(0,k,g)
                 end do
              end do
           else
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp33R( 1,k,g)*rel(i,j,k+1)) / cp33R(0,k,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L > 0) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cv22R(0,j,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cv22R( 1,j,g)*rel(i,j+1,k)) / cv22R(0,j,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           do k = S33R, N33R
              do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cv22R(-1,j,g)*rel(i,j-1,k) + cp33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
        else
           do k = S33R, N33R
             do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cv22R(-1,j,g)*rel(i,j-1,k) + cv22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cp33R(-1,k,g)*rel(i,j,k-1) + cp33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cv22R(0,j,g) + cp33R(0,k,g)))
                 end do
              end do
           end do
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 1 .or. BC_1U == 3) rel(N1,S2R:N2R,S3R:N3R) = bb(N1,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 2 .or. BC_1U == 4) then
           i = N1
           do k = S3R, N3R
!pgi$ unroll = n:8
              do j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 1 .or. BC_3U == 3) rel(S1R:N1R,S2R:N2R,N3) = bb(S1R:N1R,S2R:N2R,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U == 2 .or. BC_3U == 4) then
           k = N3
           do j = S2R, N2R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp33R(-1,k,g)*rel(i,j,k-1)) / cp33R(0,k,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U > 0) then
           j = N2
           do k = S3R, N3R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cv22R(-1,j,g)*rel(i,j-1,k)) / cv22R(0,j,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (direction == 3) then
        
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 1 .or. BC_1L == 3) rel(1 ,S2R:N2R,S3R:N3R) = bb(1 ,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1L == 2 .or. BC_1L == 4) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = bb(i,j,k) / cp11R(0,i,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do j = S2R, N2R
                    rel(i,j,k) = (bb(i,j,k) - cp11R( 1,i,g)*rel(i+1,j,k)) / cp11R(0,i,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 1 .or. BC_2L == 3) rel(S1R:N1R,1 ,S3R:N3R) = bb(S1R:N1R,1 ,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2L == 2 .or. BC_2L == 4) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cp22R(0,j,g)
                 end do
              end do
           else
              do k = S3R, N3R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cp22R( 1,j,g)*rel(i,j+1,k)) / cp22R(0,j,g)
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_3L > 0) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = bb(i,j,k) / cw33R(0,k,g)
                 end do
              end do
           else
              do j = S2R, N2R
!pgi$ unroll = n:8
                 do i = S1R, N1R
                    rel(i,j,k) = (bb(i,j,k) - cw33R( 1,k,g)*rel(i,j,k+1)) / cw33R(0,k,g)
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Feld --------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (r == 1 .and. init_yes) then
           do k = S33R, N33R
              do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                             &
                       &    cp11R(-1,i,g)*rel(i-1,j,k) + cp22R(-1,j,g)*rel(i,j-1,k) + cw33R(-1,k,g)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 end do
              end do
           end do
        else
           do k = S33R, N33R
             do j = S22R, N22R
!pgi$ unroll = n:8
                 do i = S11R, N11R
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cp11R(-1,i,g)*rel(i-1,j,k) + cp11R(1,i,g)*rel(i+1,j,k)     &
                                &      + cp22R(-1,j,g)*rel(i,j-1,k) + cp22R(1,j,g)*rel(i,j+1,k)     &
                                &      + cw33R(-1,k,g)*rel(i,j,k-1) + cw33R(1,k,g)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cp11R(0,i,g) + cp22R(0,j,g) + cw33R(0,k,g)))
                 end do
              end do
           end do
        end if
        !=====================================================================================================
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-tangentiale Geschwindigkeiten --------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 1 .or. BC_1U == 3) rel(N1,S2R:N2R,S3R:N3R) = bb(N1,S2R:N2R,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_1U == 2 .or. BC_1U == 4) then
           i = N1
           do k = S3R, N3R
!pgi$ unroll = n:8
              do j = S2R, N2R
                 rel(i,j,k) = (bb(i,j,k) - cp11R(-1,i,g)*rel(i-1,j,k)) / cp11R(0,i,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 1 .or. BC_2U == 3) rel(S1R:N1R,N2,S3R:N3R) = bb(S1R:N1R,N2,S3R:N3R)
        !-----------------------------------------------------------------------------------------------------
        if (BC_2U == 2 .or. BC_2U == 4) then
           j = N2
           do k = S3R, N3R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cp22R(-1,j,g)*rel(i,j-1,k)) / cp22R(0,j,g)
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !--- Wand-normale Geschwindigkeit --------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        if (BC_3U > 0) then
           k = N3
           do j = S2R, N2R
!pgi$ unroll = n:8
              do i = S1R, N1R
                 rel(i,j,k) = (bb(i,j,k) - cw33R(-1,k,g)*rel(i,j,k-1)) / cw33R(0,k,g)
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  end do
  
  
  end subroutine relaxation_Helmholtz_coarse
  
  
  
  
  
  
  
  
  
  
  
  subroutine relaxation_Helmholtz_conc(init_yes,n_relax,g,bb,rel) ! TEST!!! 2D-Variante fehlt noch ...
  
  implicit none
  
  logical, intent(in   ) ::  init_yes
  integer, intent(in   ) ::  n_relax
  
  integer, intent(in   ) ::  g
  
  real   , intent(in   ) ::  bb (b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  real   , intent(inout) ::  rel(b1L:(NN(1,g)+b1U),b2L:(NN(2,g)+b2U),b3L:(NN(3,g)+b3U))
  
  integer                ::  i, N1, S1Rc, S11Rc, N1Rc, N11Rc
  integer                ::  j, N2, S2Rc, S22Rc, N2Rc, N22Rc
  integer                ::  k, N3, S3Rc, S33Rc, N3Rc, N33Rc
  
  integer                ::  r, m
  
  real                   ::  mul1, mul2, mul3
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Null-Setzen am Rand nicht notwendig, da Startindizes entsprechend gew�hlt sind!           !
  !              - Bei der Initialisierung m�ssen die Intervallgrenzen SiR:NiR anstelle von SiiR:NiiR        !
  !                gew�hlt werden, da beim  Aufbau der RHS ("vec") innerhalb der Linienrelaxation IMMER auch !
  !                die Randbereiche SiR und NiR aufgebaut aber ggf. mit einer Korrektur wieder �berschrieben !
  !                werden. Ansonsten w�re eine Initialisierung allein im Feldbereich SiiR:NiiR ausreichend,  !
  !                kann aber z.B. zu Floating Point Exceptions f�hren (die f�r die eigentliche Rechnung      !
  !                allerdings irrelevant w�ren)!                                                             !
  !              - obere Stirnfl�chen werden bei der Initialisierung ebenfalls ber�cksichtigt, da ggf.       !
  !                verschiedene Richtungen nacheinander bearbeitet werden und dies beim Aufbauen der rechten !
  !                Seite sonst ber�cksichtigt werden m�sste. ==> geringerer Mehraufwand beim Rechnen,        !
  !                weniger Programmierdetails.                                                               !
  !              - LU-Zerlegung (d.h. dia, mult) kann gespeichert werden, wenn mindestens eine Richtung      !
  !                �quidistant ist. Der L�sungsaufwand w�rde sich etwa halbieren!!                           !
  !              - "r == 1 .AND. init_yes" sollte idealerweise aus den Schleifen herausgezogen werden, was   !
  !                hier aber aus Gr�nden der �bersicht bisher nicht ausgef�hrt wurde.                        !
  !----------------------------------------------------------------------------------------------------------!
  
  
  ! TEST!!! dx_ und multL k�nnte/sollte fr�her oder sp�ter direkt in die Stencils eingebaut werden! 
  !         vor jedem RK-Schritt neu berechnen --> insgesamt viel billiger als das hier ...
  
  m = conc_nu
  
  !===========================================================================================================
  N1 = NN(1,g)
  N2 = NN(2,g)
  N3 = NN(3,g)
  
  S1Rc  = S1R ! TEST!!! evtl. herausziehen nach init_limits!
  S2Rc  = S2R ! TEST!!! Substituieren ...
  S3Rc  = S3R
  
  S11Rc = S11R
  S22Rc = S22R
  S33Rc = S33R
  
  N1Rc  = N1 + d1R
  N2Rc  = N2 + d2R
  N3Rc  = N3 + d3R
  
  N11Rc = N1 + d11R
  N22Rc = N2 + d22R
  N33Rc = N3 + d33R
  
  !-----------------------------------------------------------------------------------------------------------
  !--- Sonderf�lle (Symmetrie vs. Wand) ----------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_1L(m) > 0) then
     S1Rc  = 1
     S11Rc = 2
  end if
  if (BCc_1U(m) > 0) then
     N1Rc  = N1
     N11Rc = N1-1
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_2L(m) > 0) then
     S2Rc  = 1
     S22Rc = 2
  end if
  if (BCc_2U(m) > 0) then
     N2Rc  = N2
     N22Rc = N2-1
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BCc_3L(m) > 0) then
     S3Rc  = 1
     S33Rc = 2
  end if
  if (BCc_3U(m) > 0) then
     N3Rc  = N3
     N33Rc = N3-1
  end if
  !===========================================================================================================
  
  
  ! Nur zum testen, falls Initialisierung unten fehlerhaft erscheint:
  !IF (init_yes) rel = 0.
  
  if (init_yes) then
     if (BCc_1L(m) <= 0) rel(S11Rc-1 ,S2Rc:N2Rc,S3Rc:N3Rc) = 0.
     if (BCc_1U(m) <= 0) rel(N11Rc+1 ,S2Rc:N2Rc,S3Rc:N3Rc) = 0.
     if (BCc_2L(m) <= 0) rel(S1Rc:N1Rc,S22Rc-1 ,S3Rc:N3Rc) = 0.
     if (BCc_2U(m) <= 0) rel(S1Rc:N1Rc,N22Rc+1 ,S3Rc:N3Rc) = 0.
     if (BCc_3L(m) <= 0) rel(S1Rc:N1Rc,S2Rc:N2Rc,S33Rc-1 ) = 0.
     if (BCc_3U(m) <= 0) rel(S1Rc:N1Rc,S2Rc:N2Rc,N33Rc+1 ) = 0.
  end if
  
  
  do r = 1, n_relax
     
     if (.not. (r == 1 .and. init_yes)) call exchange_relax(g,0,0,0,0,.true.,rel)
     
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. impl_dir(3) == 0) then
        
        !=====================================================================================================
        if (BCc_1L(m) == 1) rel(1 ,S22Rc:N22Rc,S33Rc:N33Rc) = bb(1 ,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1L(m) == 3) then
           i = 1
           if (r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = idx_1L*bb(i,j,k) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = (idx_1L*bb(i,j,k) - cc11R( 1,i,g,m)*rel(i+1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 1) rel(S11Rc:N11Rc,1 ,S33Rc:N33Rc) = bb(S11Rc:N11Rc,1 ,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_2L*bb(i,j,k) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_2L*bb(i,j,k) - cc22R( 1,j,g,m)*rel(i,j+1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,1 ) = bb(S11Rc:N11Rc,S22Rc:N22Rc,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_3L*bb(i,j,k) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           else
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_3L*bb(i,j,k) - cc33R( 1,k,g,m)*rel(i,j,k+1)) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        if (r == 1 .and. init_yes) then
           
           do k = S33Rc, N33Rc
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                                  &
                       &        cc11R(-1,i,g,m)*rel(i-1,j,k) + cc22R(-1,j,g,m)*rel(i,j-1,k) + cc33R(-1,k,g,m)*rel(i,j,k-1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m)))
                 end do
              end do
           end do
           
        else
           
           do k = S33Rc, N33Rc
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = bb(i,j,k) + multL*(                                                &
                                &        cc11R(-1,i,g,m)*rel(i-1,j,k) + cc11R(1,i,g,m)*rel(i+1,j,k)     &
                                &      + cc22R(-1,j,g,m)*rel(i,j-1,k) + cc22R(1,j,g,m)*rel(i,j+1,k)     &
                                &      + cc33R(-1,k,g,m)*rel(i,j,k-1) + cc33R(1,k,g,m)*rel(i,j,k+1))
                    
                    rel(i,j,k) = rel(i,j,k) / (1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m)))
                 end do
              end do
           end do
           
        end if
        !=====================================================================================================
        if (BCc_1U(m) == 1) rel(N1,S22Rc:N22Rc,S33Rc:N33Rc) = bb(N1,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1U(m) == 3) then
           i = N1
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do j = S22Rc, N22Rc
                 rel(i,j,k) = (idx_1U*bb(i,j,k) - cc11R(-1,i,g,m)*rel(i-1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,2,g))
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2U(m) == 1) rel(S11Rc:N11Rc,N2,S33Rc:N33Rc) = bb(S11Rc:N11Rc,N2,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2U(m) == 3) then
           j = N2
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_2U*bb(i,j,k) - cc22R(-1,j,g,m)*rel(i,j-1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,2,g))
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,N3) = bb(S11Rc:N11Rc,S22Rc:N22Rc,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 3) then
           k = N3
           do j = S22Rc, N22Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_3U*bb(i,j,k) - cc33R(-1,k,g,m)*rel(i,j,k-1)) / (cc33R(0,k,g,m) - velReSc3(i,j,2,g))
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (impl_dir(1) == 1) then
        
        !=====================================================================================================
        if (BCc_2L(m) == 1) rel(S11Rc:N11Rc,1 ,S33Rc:N33Rc) = bb(S11Rc:N11Rc,1 ,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
           j = 1
           if (r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_2L*bb(i,j,k) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_2L*bb(i,j,k) - cc22R( 1,j,g,m)*rel(i,j+1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,1 ) = bb(S11Rc:N11Rc,S22Rc:N22Rc,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
           k = 1
           if (r == 1 .and. init_yes) then
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_3L*bb(i,j,k) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           else
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_3L*bb(i,j,k) - cc33R( 1,k,g,m)*rel(i,j,k+1)) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        do k = S33Rc, N33Rc
           do j = S22Rc, N22Rc
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              if (r == 1 .and. init_yes) then
!pgi$ unroll = n:8
                 do i = S1Rc, N1Rc
                    dia1(i) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec1(i) = bb(i,j,k) + multL*(                                                        &
                                  &     + cc22R(-1,j,g,m)*rel(i,j-1,k) + cc33R(-1,k,g,m)*rel(i,j,k-1))
                 end do
              else
!pgi$ unroll = n:8
                 do i = S1Rc, N1Rc
                    dia1(i) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec1(i) = bb(i,j,k) + multL*(                                                        &
                                  &     + cc22R(-1,j,g,m)*rel(i,j-1,k) + cc33R(-1,k,g,m)*rel(i,j,k-1)    &
                                  &     + cc22R( 1,j,g,m)*rel(i,j+1,k) + cc33R( 1,k,g,m)*rel(i,j,k+1))
                 end do
              end if
              
              
              !--- mit RB �berschreiben ----------------------------------------------------------------------
              if (BCc_1L(m) == 3) then
                 dia1(S1Rc  ) = dx_1L*(cc11R(0,S1Rc,g,m) - velReSc1(j,k,1,g))
                 vec1(S1Rc  ) = bb(S1Rc,j,k)
                 mul1         = multL*cc11R(-1,S1Rc+1,g,m) / dia1(S1Rc)
                 dia1(S1Rc+1) = dia1(S1Rc+1) + mul1*dx_1L*cc11R(1,S1Rc,g,m)
              else if (BCc_1L(m) == 1) then
                 dia1(S1Rc  ) = 1.
                 vec1(S1Rc  ) = bb(S1Rc,j,k)
                 mul1         = multL*cc11R(-1,S1Rc+1,g,m)!/ dia1(S1Rc)
              else
                 vec1(S1Rc  ) = vec1(S1Rc) + multL*cc11R(-1,S1Rc,g,m)*rel(S1Rc-1,j,k)
                 mul1         = multL*cc11R(-1,S1Rc+1,g,m) / dia1(S1Rc)
                 dia1(S1Rc+1) = dia1(S1Rc+1) - mul1*multL*cc11R(1,S1Rc,g,m)
              end if
              vec1(S1Rc+1) = vec1(S1Rc+1) + mul1*vec1(S1Rc)
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              do i = S1Rc+2, N1Rc-1
                 mul1    = multL*cc11R(-1,i,g,m) / dia1(i-1)
                 dia1(i) = dia1(i) - mul1*multL*cc11R(1,i-1,g,m)
                 vec1(i) = vec1(i) + mul1*vec1(i-1)
              end do
              
              
              !--- RB einf�gen -------------------------------------------------------------------------------
              if (BCc_1U(m) == 3) then
                 dia1(N1Rc) = dx_1U*(cc11R(0,N1Rc,g,m) - velReSc1(j,k,2,g))
                 vec1(N1Rc) = bb(N1Rc,j,k)
                 mul1       = dx_1U*cc11R(-1,N1Rc,g,m) / dia1(N1Rc-1)
              else if (BCc_1U(m) == 1) then
                 dia1(N1Rc) = 1.
                 vec1(N1Rc) = bb(N1Rc,j,k)
                 mul1       = 0.
              else
                 vec1(N1Rc) = vec1(N1Rc) + multL*cc11R(1,N1Rc,g,m)*rel(N1Rc+1,j,k)
                 mul1       = -multL*cc11R(-1,N1Rc,g,m) / dia1(N1Rc-1)
              end if
              dia1(N1Rc) = dia1(N1Rc) + mul1 * multL*cc11R(1,N1Rc-1,g,m)
              vec1(N1Rc) = vec1(N1Rc) - mul1 * vec1(N1Rc-1)
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(N1Rc,j,k) = vec1(N1Rc) / dia1(N1Rc)
!pgi$ unroll = n:8
              do i = N1Rc-1, S1Rc+1, -1
                 rel(i,j,k) = (vec1(i) + multL*cc11R(1,i,g,m)*rel(i+1,j,k)) / dia1(i)
              end do
              
              
              if (BCc_1L(m) == 3) then
                 rel(S1Rc,j,k) = (vec1(S1Rc) - dx_1L*cc11R(1,S1Rc,g,m)*rel(S1Rc+1,j,k)) / dia1(S1Rc)
              else if (BCc_1L(m) == 1) then
                 rel(S1Rc,j,k) = vec1(S1Rc)
              else
                 rel(S1Rc,j,k) = (vec1(S1Rc) + multL*cc11R(1,S1Rc,g,m)*rel(S1Rc+1,j,k)) / dia1(S1Rc)
              end if
              
           end do
        end do
        !=====================================================================================================
        if (BCc_2U(m) == 1) rel(S11Rc:N11Rc,N2,S33Rc:N33Rc) = bb(S11Rc:N11Rc,N2,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2U(m) == 3) then
           j = N2
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_2U*bb(i,j,k) - cc22R(-1,j,g,m)*rel(i,j-1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,2,g))
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,N3) = bb(S11Rc:N11Rc,S22Rc:N22Rc,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 3) then
           k = N3
           do j = S22Rc, N22Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_3U*bb(i,j,k) - cc33R(-1,k,g,m)*rel(i,j,k-1)) / (cc33R(0,k,g,m) - velReSc3(i,j,2,g))
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (impl_dir(2) == 1) then
        
        !=====================================================================================================
        if (BCc_1L(m) == 1) rel(1 ,S22Rc:N22Rc,S33Rc:N33Rc) = bb(1 ,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1L(m) == 3) then
           i = 1
           if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = idx_1L*bb(i,j,k) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = (idx_1L*bb(i,j,k) - cc11R( 1,i,g,m)*rel(i+1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,1 ) = bb(S11Rc:N11Rc,S22Rc:N22Rc,1 )
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3L(m) == 3) then
           k = 1
           if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_3L*bb(i,j,k) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           else
              do j = S22Rc, N22Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_3L*bb(i,j,k) - cc33R( 1,k,g,m)*rel(i,j,k+1)) / (cc33R(0,k,g,m) - velReSc3(i,j,1,g))
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        do k = S33Rc, N33Rc
           do i = S11Rc, N11Rc
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              if (impl_dir(1) == 0 .and. r == 1 .and. init_yes) then
!pgi$ unroll = n:8
                 do j = S2Rc, N2Rc
                    dia2(j) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec2(j) = bb(i,j,k) + multL*(                                                       &
                                  &     + cc11R(-1,i,g,m)*rel(i-1,j,k) + cc33R(-1,k,g,m)*rel(i,j,k-1))
                 end do
              else
!pgi$ unroll = n:8
                 do j = S2Rc, N2Rc
                    dia2(j) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec2(j) = bb(i,j,k) + multL*(                                                        &
                                  &     + cc11R(-1,i,g,m)*rel(i-1,j,k) + cc33R(-1,k,g,m)*rel(i,j,k-1)    &
                                  &     + cc11R( 1,i,g,m)*rel(i+1,j,k) + cc33R( 1,k,g,m)*rel(i,j,k+1))
                 end do
              end if
              
              
              !--- mit RB �berschreiben ----------------------------------------------------------------------
              if (BCc_2L(m) == 3) then
                 dia2(S2Rc  ) = dx_2L*(cc22R(0,S2Rc,g,m) - velReSc2(i,k,1,g))
                 vec2(S2Rc  ) = bb(i,S2Rc,k)
                 mul2         = multL*cc22R(-1,S2Rc+1,g,m) / dia2(S2Rc)
                 dia2(S2Rc+1) = dia2(S2Rc+1) + mul2*dx_2L*cc22R(1,S2Rc,g,m)
              else if (BCc_2L(m) == 1) then
                 dia2(S2Rc  ) = 1.
                 vec2(S2Rc  ) = bb(i,S2Rc,k)
                 mul2         = multL*cc22R(-1,S2Rc+1,g,m)!/ dia2(S2Rc)
              else
                 vec2(S2Rc  ) = vec2(S2Rc) + multL*cc22R(-1,S2Rc,g,m)*rel(i,S2Rc-1,k)
                 mul2         = multL*cc22R(-1,S2Rc+1,g,m) / dia2(S2Rc)
                 dia2(S2Rc+1) = dia2(S2Rc+1) - mul2*multL*cc22R(1,S2Rc,g,m)
              end if
              vec2(S2Rc+1) = vec2(S2Rc+1) + mul2*vec2(S2Rc)
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              do j = S2Rc+2, N2Rc-1
                 mul2    = multL*cc22R(-1,j,g,m) / dia2(j-1)
                 dia2(j) = dia2(j) - mul2*multL*cc22R(1,j-1,g,m)
                 vec2(j) = vec2(j) + mul2*vec2(j-1)
              end do
              
              
              !--- RB einf�gen -------------------------------------------------------------------------------
              if (BCc_2U(m) == 3) then
                 dia2(N2Rc) = dx_2U*(cc22R(0,N2Rc,g,m) - velReSc2(i,k,2,g))
                 vec2(N2Rc) = bb(i,N2Rc,k)
                 mul2       = dx_2U*cc22R(-1,N2Rc,g,m) / dia2(N2Rc-1)
              else if (BCc_2U(m) == 1) then
                 dia2(N2Rc) = 1.
                 vec2(N2Rc) = bb(i,N2Rc,k)
                 mul2       = 0.
              else
                 vec2(N2Rc) = vec2(N2Rc) + multL*cc22R(1,N2Rc,g,m)*rel(i,N2Rc+1,k)
                 mul2       = -multL*cc22R(-1,N2Rc,g,m) / dia2(N2Rc-1)
              end if
              dia2(N2Rc) = dia2(N2Rc) + mul2 * multL*cc22R(1,N2Rc-1,g,m)
              vec2(N2Rc) = vec2(N2Rc) - mul2 * vec2(N2Rc-1)
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,N2Rc,k) = vec2(N2Rc) / dia2(N2Rc)
!pgi$ unroll = n:8
              do j = N2Rc-1, S2Rc+1, -1
                 rel(i,j,k) = (vec2(j) + multL*cc22R(1,j,g,m)*rel(i,j+1,k)) / dia2(j)
              end do
              
              
              if (BCc_2L(m) == 3) then
                 rel(i,S2Rc,k) = (vec2(S2Rc) - dx_2L*cc22R(1,S2Rc,g,m)*rel(i,S2Rc+1,k)) / dia2(S2Rc)
              else if (BCc_2L(m) == 1) then
                 rel(i,S2Rc,k) = vec2(S2Rc)
              else
                 rel(i,S2Rc,k) = (vec2(S2Rc) + multL*cc22R(1,S2Rc,g,m)*rel(i,S2Rc+1,k)) / dia2(S2Rc)
              end if
              
           end do
        end do
        !=====================================================================================================
        if (BCc_1U(m) == 1) rel(N1,S22Rc:N22Rc,S33Rc:N33Rc) = bb(N1,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1U(m) == 3) then
           i = N1
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do j = S22Rc, N22Rc
                 rel(i,j,k) = (idx_1U*bb(i,j,k) - cc11R(-1,i,g,m)*rel(i-1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,2,g))
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 1) rel(S11Rc:N11Rc,S22Rc:N22Rc,N3) = bb(S11Rc:N11Rc,S22Rc:N22Rc,N3)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_3U(m) == 3) then
           k = N3
           do j = S22Rc, N22Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_3U*bb(i,j,k) - cc33R(-1,k,g,m)*rel(i,j,k-1)) / (cc33R(0,k,g,m) - velReSc3(i,j,2,g))
              end do
           end do
        end if
        !=====================================================================================================i
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     if (impl_dir(3) == 1) then
        
        !=====================================================================================================
        if (BCc_1L(m) == 1) rel(1 ,S22Rc:N22Rc,S33Rc:N33Rc) = bb(1 ,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1L(m) == 3) then
           i = 1
           if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = idx_1L*bb(i,j,k) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do j = S22Rc, N22Rc
                    rel(i,j,k) = (idx_1L*bb(i,j,k) - cc11R( 1,i,g,m)*rel(i+1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,1,g))
                 end do
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 1) rel(S11Rc:N11Rc,1 ,S33Rc:N33Rc) = bb(S11Rc:N11Rc,1 ,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2L(m) == 3) then
           j = 1
           if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = idx_2L*bb(i,j,k) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           else
              do k = S33Rc, N33Rc
!pgi$ unroll = n:8
                 do i = S11Rc, N11Rc
                    rel(i,j,k) = (idx_2L*bb(i,j,k) - cc22R( 1,j,g,m)*rel(i,j+1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,1,g))
                 end do
              end do
           end if
        end if
        !=====================================================================================================
        do j = S22Rc, N22Rc
           do i = S11Rc, N11Rc
              
              !--- Diagonalelement / rechte Seite aufbauen ---------------------------------------------------
              if (impl_dir(1) == 0 .and. impl_dir(2) == 0 .and. r == 1 .and. init_yes) then
!pgi$ unroll = n:8
                 do k = S3Rc, N3Rc
                    dia3(k) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec3(k) = bb(i,j,k) + multL*(                                                       &
                                  &     + cc11R(-1,i,g,m)*rel(i-1,j,k) + cc22R(-1,j,g,m)*rel(i,j-1,k))
                 end do
              else
!pgi$ unroll = n:8
                 do k = S3Rc, N3Rc
                    dia3(k) = 1. - multL*(cc11R(0,i,g,m) + cc22R(0,j,g,m) + cc33R(0,k,g,m))
                    vec3(k) = bb(i,j,k) + multL*(                                                        &
                                  &     + cc11R(-1,i,g,m)*rel(i-1,j,k) + cc22R(-1,j,g,m)*rel(i,j-1,k)    &
                                  &     + cc11R( 1,i,g,m)*rel(i+1,j,k) + cc22R( 1,j,g,m)*rel(i,j+1,k))
                 end do
              end if
              
              
              !--- mit RB �berschreiben ----------------------------------------------------------------------
              if (BCc_3L(m) == 3) then
                 dia3(S3Rc  ) = dx_3L*(cc33R(0,S3Rc,g,m) - velReSc3(i,j,1,g))
                 vec3(S3Rc  ) = bb(i,j,S3Rc)
                 mul3         = multL*cc33R(-1,S3Rc+1,g,m) / dia3(S3Rc)
                 dia3(S3Rc+1) = dia3(S3Rc+1) + mul3*dx_3L*cc33R(1,S3Rc,g,m)
              else if (BCc_3L(m) == 1) then
                 dia3(S3Rc  ) = 1.
                 vec3(S3Rc  ) = bb(i,j,S3Rc)
                 mul3         = multL*cc33R(-1,S3Rc+1,g,m)!/ dia3(S3Rc)
              else
                 vec3(S3Rc  ) = vec3(S3Rc) + multL*cc33R(-1,S3Rc,g,m)*rel(i,j,S3Rc-1)
                 mul3         = multL*cc33R(-1,S3Rc+1,g,m) / dia3(S3Rc)
                 dia3(S3Rc+1) = dia3(S3Rc+1) - mul3*multL*cc33R(1,S3Rc,g,m)
              end if
              vec3(S3Rc+1) = vec3(S3Rc+1) + mul3*vec3(S3Rc)
              
              
              !--- Gauss-Elimination hoch --------------------------------------------------------------------
!pgi$ unroll = n:8
              do k = S3Rc+2, N3Rc-1
                 mul3    = multL*cc33R(-1,k,g,m) / dia3(k-1)
                 dia3(k) = dia3(k) - mul3*multL*cc33R(1,k-1,g,m)
                 vec3(k) = vec3(k) + mul3*vec3(k-1)
              end do
              
              
              !--- RB einf�gen -------------------------------------------------------------------------------
              if (BCc_3U(m) == 3) then
                 dia3(N3Rc) = dx_3U*(cc33R(0,N3Rc,g,m) - velReSc3(i,j,2,g))
                 vec3(N3Rc) = bb(i,j,N3Rc)
                 mul3       = dx_3U*cc33R(-1,N3Rc,g,m) / dia3(N3Rc-1)
              else if (BCc_3U(m) == 1) then
                 dia3(N3Rc) = 1.
                 vec3(N3Rc) = bb(i,j,N3Rc)
                 mul3       = 0.
              else
                 vec3(N3Rc) = vec3(N3Rc) + multL*cc33R(1,N3Rc,g,m)*rel(i,j,N3Rc+1)
                 mul3       = -multL*cc33R(-1,N3Rc,g,m) / dia3(N3Rc-1)
              end if
              dia3(N3Rc) = dia3(N3Rc) + mul3 * multL*cc33R(1,N3Rc-1,g,m)
              vec3(N3Rc) = vec3(N3Rc) - mul3 * vec3(N3Rc-1)
              
              
              !--- Gauss-Elimination runter ------------------------------------------------------------------
              rel(i,j,N3Rc) = vec3(N3Rc) / dia3(N3Rc)
!pgi$ unroll = n:8
              do k = N3Rc-1, S3Rc+1, -1
                 rel(i,j,k) = (vec3(k) + multL*cc33R(1,k,g,m)*rel(i,j,k+1)) / dia3(k)
              end do
              
              
              if (BCc_3L(m) == 3) then
                 rel(i,j,S3Rc) = (vec3(S3Rc) - dx_3L*cc33R(1,S3Rc,g,m)*rel(i,j,S3Rc+1)) / dia3(S3Rc)
              else if (BCc_3L(m) == 1) then
                 rel(i,j,S3Rc) = vec3(S3Rc)
              else
                 rel(i,j,S3Rc) = (vec3(S3Rc) + multL*cc33R(1,S3Rc,g,m)*rel(i,j,S3Rc+1)) / dia3(S3Rc)
              end if
              
           end do
        end do
        !=====================================================================================================
        if (BCc_1U(m) == 1) rel(N1,S22Rc:N22Rc,S33Rc:N33Rc) = bb(N1,S22Rc:N22Rc,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_1U(m) == 3) then
           i = N1
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do j = S22Rc, N22Rc
                 rel(i,j,k) = (idx_1U*bb(i,j,k) - cc11R(-1,i,g,m)*rel(i-1,j,k)) / (cc11R(0,i,g,m) - velReSc1(j,k,2,g))
              end do
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2U(m) == 1) rel(S11Rc:N11Rc,N2,S33Rc:N33Rc) = bb(S11Rc:N11Rc,N2,S33Rc:N33Rc)
        !-----------------------------------------------------------------------------------------------------
        if (BCc_2U(m) == 3) then
           j = N2
           do k = S33Rc, N33Rc
!pgi$ unroll = n:8
              do i = S11Rc, N11Rc
                 rel(i,j,k) = (idx_2U*bb(i,j,k) - cc22R(-1,j,g,m)*rel(i,j-1,k)) / (cc22R(0,j,g,m) - velReSc2(i,k,2,g))
              end do
           end do
        end if
        !=====================================================================================================
        
     end if
     !========================================================================================================
     !========================================================================================================
     !========================================================================================================
     
  end do
  
  if (corner_yes) call handle_corner_conc(N1,N2,N3,rel)
  
  
  end subroutine relaxation_Helmholtz_conc
  
  
  
  
  
  
  
  
  
  
  
  subroutine handle_corner_conc(N1,N2,N3,phi)
  
  implicit none
  
  integer, intent(in   ) ::  N1
  integer, intent(in   ) ::  N2
  integer, intent(in   ) ::  N3
  
  real   , intent(inout) ::  phi(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Siehe Subroutines "handle_corner_Lap" und "handle_corner_rhs"!                            !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (BCc_1L(conc_nu) > 0 .and. BCc_2L(conc_nu) > 0) phi(1 ,1 ,1:N3) = 0.
  if (BCc_1L(conc_nu) > 0 .and. BCc_2U(conc_nu) > 0) phi(1 ,N2,1:N3) = 0.
  if (BCc_1U(conc_nu) > 0 .and. BCc_2L(conc_nu) > 0) phi(N1,1 ,1:N3) = 0.
  if (BCc_1U(conc_nu) > 0 .and. BCc_2U(conc_nu) > 0) phi(N1,N2,1:N3) = 0.
  
  if (BCc_1L(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1 ,1:N2,1 ) = 0.
  if (BCc_1L(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1 ,1:N2,N3) = 0.
  if (BCc_1U(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(N1,1:N2,1 ) = 0.
  if (BCc_1U(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(N1,1:N2,N3) = 0.
  
  if (BCc_2L(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1:N1,1 ,1 ) = 0.
  if (BCc_2L(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1:N1,1 ,N3) = 0.
  if (BCc_2U(conc_nu) > 0 .and. BCc_3L(conc_nu) > 0) phi(1:N1,N2,1 ) = 0.
  if (BCc_2U(conc_nu) > 0 .and. BCc_3U(conc_nu) > 0) phi(1:N1,N2,N3) = 0.
  
  
  end subroutine handle_corner_conc
  
  
  
end module mod_helmholtz
