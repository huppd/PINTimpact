!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing functions to initiliaze stencil arrays
module mod_coeffs
  
  
  use mod_dims
  use mod_vars
  use mod_exchange
  use mod_diff ! (apply_compact)
  
  
  private
  
  public FD_coeffs_compact, init_compact, diff_coeffs_compact, compact_parallel, test_coeffs_compact
  public FD_coeffs, test_diff, diff_coeffs, diff_coeffs_exact, test_coeffs
  public FD_coeffs_solver, FD_coeffs_solver_integral
  public interp_coeffs, interp_coeffs_Helm, restr_coeffs, restr_coeffs_Helm
  public get_stencil, get_stencil_Helm, get_stencil_transp
  public get_weights
  public Matrix_invert
  
  public myDGBTRF
  
  INCLUDE 'mpif.h'
  
  contains
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  subroutine FD_coeffs_compact()
  
  implicit none
  
  integer                ::  i, j, k
  !----------------------------------------------------
  integer             ::  disp1_global(1:NB1), recv1_global(1:NB1), rank_bar1 ! TEST!!!
  integer             ::  disp2_global(1:NB2), recv2_global(1:NB2), rank_bar2
  integer             ::  disp3_global(1:NB3), recv3_global(1:NB3), rank_bar3
  !----------------------------------------------------
  
  
  ! Anmerkung: So ist das sicherlich noch nicht sehr effizient, sollte aber grundsaetzlich die richtige Strategie sein ...
  
  
  ! TEST!!!
  !=====================================================================================================
  disp1 = 0
  disp2 = 0
  disp3 = 0
  
  call MPI_COMM_RANK(COMM_BAR1,rank_bar1,merror)
  call MPI_COMM_RANK(COMM_BAR2,rank_bar2,merror)
  call MPI_COMM_RANK(COMM_BAR3,rank_bar3,merror)
  
  disp1(rank_bar1+1) = 2*ndL*(iB(1,1)-1)*(N2+1)*(N3+1)
  disp2(rank_bar2+1) = 2*ndL*(iB(2,1)-1)*(N1+1)*(N3+1)
  disp3(rank_bar3+1) = 2*ndL*(iB(3,1)-1)*(N1+1)*(N2+1)
  
  call MPI_ALLREDUCE(disp1,disp1_global,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
  call MPI_ALLREDUCE(disp2,disp2_global,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
  call MPI_ALLREDUCE(disp3,disp3_global,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
  disp1 = disp1_global
  disp2 = disp2_global
  disp3 = disp3_global
  
  recv1 = 2*ndl*(N2+1)*(N3+1)
  recv2 = 2*ndl*(N1+1)*(N3+1)
  recv3 = 2*ndl*(N1+1)*(N2+1)
  !=====================================================================================================
  
  
  
  
  call init_compact(-1,-1)
  
  !=====================================================================================================
  j = 1
  k = 1
  do i = S1p, N1p
     pp(i,j,k) = x1p(i) - (REAL(i+iShift)-1.0)*L1/REAL(M1-1)
  end do
  
  call exchange(1,0,pp)
  call apply_compact(1,0,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,pp,Ap)
  call apply_compact(1,0,S12,1,1,N12,1,1,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ar)
  call apply_compact(1,0,S12,1,1,N12,1,1,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,pp,rh)
  
  do i = S11, N11
     dx1GM (i) = 1./(Ap(i,j,k) + L1/REAL(M1-1))
  end do
  do i = S12, N12
     dx1pM (i) = 1./(Ar(i,j,k) + L1/REAL(M1-1))
     ddx1pM(i) = -rh(i,j,k)*dx1pM(i)**3
  end do
  !-----------------------------------------------------------------------------------------------------
  do i = S11B, N11B
     pp(i,j,k) = x1u(i) - (REAL(i+iShift)-0.5)*L1/REAL(M1-1)
  end do
  
  call exchange(1,1,pp)
  call apply_compact(1,1,S1p,1,1,N1p,1,1,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,pp,Ap)
  call apply_compact(1,1,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ar)
  call apply_compact(1,1,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,pp,rh)
  
  do i = S1p, N1p
     dx1DM (i) = 1./(Ap(i,j,k) + L1/REAL(M1-1))
  end do
  do i = S11, N11
     dx1uM (i) = 1./(Ar(i,j,k) + L1/REAL(M1-1))
     ddx1uM(i) = -rh(i,j,k)*dx1uM(i)**3
  end do
  !=====================================================================================================
  i = 1
  k = 1
  do j = S2p, N2p
     pp(i,j,k) = x2p(j) - (REAL(j+jShift)-1.0)*L2/REAL(M2-1)
  end do
  
  call exchange(2,0,pp)
  call apply_compact(2,0,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,pp,Ap)
  call apply_compact(2,0,1,S21,1,1,N21,1,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ar)
  call apply_compact(2,0,1,S21,1,1,N21,1,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,pp,rh)
  
  do j = S22, N22
     dx2GM (j) = 1./(Ap(i,j,k) + L2/REAL(M2-1))
  end do
  do j = S21, N21
     dx2pM (j) = 1./(Ar(i,j,k) + L2/REAL(M2-1))
     ddx2pM(j) = -rh(i,j,k)*dx2pM(j)**3
  end do
  !-----------------------------------------------------------------------------------------------------
  do j = S22B, N22B
     pp(i,j,k) = x2v(j) - (REAL(j+jShift)-0.5)*L2/REAL(M2-1)
  end do
  
  call exchange(2,2,pp)
  call apply_compact(2,2,1,S2p,1,1,N2p,1,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,pp,Ap)
  call apply_compact(2,2,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ar)
  call apply_compact(2,2,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,pp,rh)
  
  do j = S2p, N2p
     dx2DM (j) = 1./(Ap(i,j,k) + L2/REAL(M2-1))
  end do
  do j = S22, N22
     dx2vM (j) = 1./(Ar(i,j,k) + L2/REAL(M2-1))
     ddx2vM(j) = -rh(i,j,k)*dx2vM(j)**3
  end do
  !=====================================================================================================
  if (dimens == 3) then
  i = 1
  j = 1
  do k = S3p, N3p
     pp(i,j,k) = x3p(k) - (REAL(k+kShift)-1.0)*L3/REAL(M3-1)
  end do
  
  call exchange(3,0,pp)
  call apply_compact(3,0,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,pp,Ap)
  call apply_compact(3,0,1,1,S31,1,1,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ar)
  call apply_compact(3,0,1,1,S31,1,1,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,pp,rh)
  
  do k = S33, N33
     dx3GM (k) = 1./(Ap(i,j,k) + L3/REAL(M3-1))
  end do
  do k = S31, N31
     dx3pM (k) = 1./(Ar(i,j,k) + L3/REAL(M3-1))
     ddx3pM(k) = -rh(i,j,k)*dx3pM(k)**3
  end do
  !-----------------------------------------------------------------------------------------------------
  do k = S33B, N33B
     pp(i,j,k) = x3w(k) - (REAL(k+kShift)-0.5)*L3/REAL(M3-1)
  end do
  
  call exchange(3,3,pp)
  call apply_compact(3,3,1,1,S3p,1,1,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,pp,Ap)
  call apply_compact(3,3,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ar)
  call apply_compact(3,3,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,pp,rh)
  
  do k = S3p, N3p
     dx3DM (k) = 1./(Ap(i,j,k) + L3/REAL(M3-1))
  end do
  do k = S33, N33
     dx3wM (k) = 1./(Ar(i,j,k) + L3/REAL(M3-1))
     ddx3wM(k) = -rh(i,j,k)*dx3wM(k)**3
  end do
  end if
  !=====================================================================================================
  
  call init_compact( 1,-1)
  
  
  end subroutine FD_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  subroutine init_compact(sym_pre,sym_vel)
  
  implicit none
  
  integer, intent(in   ) ::  sym_pre, sym_vel
  integer                ::  m
  
  
  !--- Koeffizienten bestimmen ---
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S12,N12,BC_1L,BC_1U,0,1,sym_pre,.false.,cp1CL ,cp1CL_LU ,cp1CR )
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S21,N21,BC_2L,BC_2U,0,1,sym_pre,.false.,cp2CL ,cp2CL_LU ,cp2CR )
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S31,N31,BC_3L,BC_3U,0,1,sym_pre,.false.,cp3CL ,cp3CL_LU ,cp3CR )
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,1,1,sym_vel,.false.,cu1CL ,cu1CL_LU ,cu1CR )
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,1,1,sym_vel,.false.,cv2CL ,cv2CL_LU ,cv2CR )
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,1,1,sym_vel,.false.,cw3CL ,cw3CL_LU ,cw3CR )
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S12,N12,BC_1L,BC_1U,0,2,sym_pre,.false.,cp11CL,cp11CL_LU,cp11CR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S21,N21,BC_2L,BC_2U,0,2,sym_pre,.false.,cp22CL,cp22CL_LU,cp22CR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S31,N31,BC_3L,BC_3U,0,2,sym_pre,.false.,cp33CL,cp33CL_LU,cp33CR)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,1,2,sym_vel,.false.,cu11CL,cu11CL_LU,cu11CR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,1,2,sym_vel,.false.,cv22CL,cv22CL_LU,cv22CR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,1,2,sym_vel,.false.,cw33CL,cw33CL_LU,cw33CR)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,1,sym_vel,.false.,cDu1CL,cDu1CL_LU,cDu1CR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,1,sym_vel,.false.,cDv2CL,cDv2CL_LU,cDv2CR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,1,sym_vel,.false.,cDw3CL,cDw3CL_LU,cDw3CR)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,1,sym_vel,.true. ,cDu1CLT,cDu1CLT_LU,cDu1CRT)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,1,sym_vel,.true. ,cDv2CLT,cDv2CLT_LU,cDv2CRT)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,1,sym_vel,.true. ,cDw3CLT,cDw3CLT_LU,cDw3CRT)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,1,sym_pre,.false.,cGp1CL,cGp1CL_LU,cGp1CR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,1,sym_pre,.false.,cGp2CL,cGp2CL_LU,cGp2CR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,1,sym_pre,.false.,cGp3CL,cGp3CL_LU,cGp3CR)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,1,sym_pre,.true. ,cGp1CLT,cGp1CLT_LU,cGp1CRT)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,1,sym_pre,.true. ,cGp2CLT,cGp2CLT_LU,cGp2CRT)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,1,sym_pre,.true. ,cGp3CLT,cGp3CLT_LU,cGp3CRT)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,0,sym_vel,.false.,cIupCL,cIupCL_LU,cIupCR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,0,sym_vel,.false.,cIvpCL,cIvpCL_LU,cIvpCR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,0,sym_vel,.false.,cIwpCL,cIwpCL_LU,cIwpCR)
  
                   call diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,0,sym_pre,.false.,cIpuCL,cIpuCL_LU,cIpuCR)
                   call diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,0,sym_pre,.false.,cIpvCL,cIpvCL_LU,cIpvCR)
  if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,0,sym_pre,.false.,cIpwCL,cIpwCL_LU,cIpwCR)
  
  
  
  do m = 1, n_conc
                      call diff_coeffs_compact(N1,NB1,ndL,ndR,S1c(m),N1c(m),BCc_1L(m),BCc_1U(m),0,1,sym_pre,.false.,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m))
                      call diff_coeffs_compact(N2,NB2,ndL,ndR,S2c(m),N2c(m),BCc_2L(m),BCc_2U(m),0,1,sym_pre,.false.,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m))
     if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S3c(m),N3c(m),BCc_3L(m),BCc_3U(m),0,1,sym_pre,.false.,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m))
     
                      call diff_coeffs_compact(N1,NB1,ndL,ndR,S1c(m),N1c(m),BCc_1L(m),BCc_1U(m),0,2,sym_pre,.false.,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m))
                      call diff_coeffs_compact(N2,NB2,ndL,ndR,S2c(m),N2c(m),BCc_2L(m),BCc_2U(m),0,2,sym_pre,.false.,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m))
     if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S3c(m),N3c(m),BCc_3L(m),BCc_3U(m),0,2,sym_pre,.false.,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m))
  
                      call diff_coeffs_compact(N1,NB1,ndL,ndR,S11   ,N11   ,BCc_1L(m),BCc_1U(m),2,0,sym_pre,.false.,cIcuCL(-ndL,0,m),cIcuCL_LU(1,0,m),cIcuCR(-ndR,0,m))
                      call diff_coeffs_compact(N2,NB2,ndL,ndR,S22   ,N22   ,BCc_2L(m),BCc_2U(m),2,0,sym_pre,.false.,cIcvCL(-ndL,0,m),cIcvCL_LU(1,0,m),cIcvCR(-ndR,0,m))
     if (dimens == 3) call diff_coeffs_compact(N3,NB3,ndL,ndR,S33   ,N33   ,BCc_3L(m),BCc_3U(m),2,0,sym_pre,.false.,cIcwCL(-ndL,0,m),cIcwCL_LU(1,0,m),cIcwCR(-ndR,0,m))
  end do
  
  
  !--- Schur-Komplemente und Hilfsvektoren bestimmen ---
  ! TEST!!!
  ! Im Vergleich zu der auskommentierten Version unten wurde nur die IF-Abfrage modifiziert, da nun IMMER das
  ! Schur-Komplement fuer die Behandlung der Block-Grenzen gerechnet wird! (Pivoting-Problem an Raendern)
                   call compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp1CL ,dimS1,BC_1L_global,Sp1 ,Wp1 ,COMM_BAR1)
                   call compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp2CL ,dimS2,BC_2L_global,Sp2 ,Wp2 ,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp3CL ,dimS3,BC_3L_global,Sp3 ,Wp3 ,COMM_BAR3)
  
                   call compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu1CL ,dimS1,BC_1L_global,Su1 ,Wu1 ,COMM_BAR1)
                   call compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv2CL ,dimS2,BC_2L_global,Sv2 ,Wv2 ,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw3CL ,dimS3,BC_3L_global,Sw3 ,Ww3 ,COMM_BAR3)
  
                   call compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp11CL,dimS1,BC_1L_global,Sp11,Wp11,COMM_BAR1)
                   call compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp22CL,dimS2,BC_2L_global,Sp22,Wp22,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp33CL,dimS3,BC_3L_global,Sp33,Wp33,COMM_BAR3)
  
                   call compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu11CL,dimS1,BC_1L_global,Su11,Wu11,COMM_BAR1)
                   call compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv22CL,dimS2,BC_2L_global,Sv22,Wv22,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw33CL,dimS3,BC_3L_global,Sw33,Ww33,COMM_BAR3)
  
                   call compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CL,dimS1,BC_1L_global,SDu1,WDu1,COMM_BAR1)
                   call compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CL,dimS2,BC_2L_global,SDv2,WDv2,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CL,dimS3,BC_3L_global,SDw3,WDw3,COMM_BAR3)
  
                   call compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CLT,dimS1,BC_1L_global,SDu1T,WDu1T,COMM_BAR1)
                   call compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CLT,dimS2,BC_2L_global,SDv2T,WDv2T,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CLT,dimS3,BC_3L_global,SDw3T,WDw3T,COMM_BAR3)
  
                   call compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CL,dimS1,BC_1L_global,SGp1,WGp1,COMM_BAR1)
                   call compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CL,dimS2,BC_2L_global,SGp2,WGp2,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CL,dimS3,BC_3L_global,SGp3,WGp3,COMM_BAR3)
  
                   call compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CLT,dimS1,BC_1L_global,SGp1T,WGp1T,COMM_BAR1)
                   call compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CLT,dimS2,BC_2L_global,SGp2T,WGp2T,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CLT,dimS3,BC_3L_global,SGp3T,WGp3T,COMM_BAR3)
  
                   call compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cIupCL,dimS1,BC_1L_global,SIup,WIup,COMM_BAR1)
                   call compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cIvpCL,dimS2,BC_2L_global,SIvp,WIvp,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cIwpCL,dimS3,BC_3L_global,SIwp,WIwp,COMM_BAR3)
  
                   call compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cIpuCL,dimS1,BC_1L_global,SIpu,WIpu,COMM_BAR1)
                   call compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cIpvCL,dimS2,BC_2L_global,SIpv,WIpv,COMM_BAR2)
  if (dimens == 3) call compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cIpwCL,dimS3,BC_3L_global,SIpw,WIpw,COMM_BAR3)
  
  
  do m = 1, n_conc
                      call compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc1CL (-ndL,0,m),dimS1,BC_1L_global,Sc1 (1,1,m),Wc1 (1,0,m),COMM_BAR1)
                      call compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc2CL (-ndL,0,m),dimS2,BC_2L_global,Sc2 (1,1,m),Wc2 (1,0,m),COMM_BAR2)
     if (dimens == 3) call compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc3CL (-ndL,0,m),dimS3,BC_3L_global,Sc3 (1,1,m),Wc3 (1,0,m),COMM_BAR3)
     
                      call compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc11CL(-ndL,0,m),dimS1,BC_1L_global,Sc11(1,1,m),Wc11(1,0,m),COMM_BAR1)
                      call compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc22CL(-ndL,0,m),dimS2,BC_2L_global,Sc22(1,1,m),Wc22(1,0,m),COMM_BAR2)
     if (dimens == 3) call compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc33CL(-ndL,0,m),dimS3,BC_3L_global,Sc33(1,1,m),Wc33(1,0,m),COMM_BAR3)
     
                      call compact_parallel(S11   ,N11   ,N1,iB(1,1),NB1,ndL,cIcuCL(-ndL,0,m),dimS1,BC_1L_global,SIcu(1,1,m),WIcu(1,0,m),COMM_BAR1)
                      call compact_parallel(S22   ,N22   ,N2,iB(2,1),NB2,ndL,cIcvCL(-ndL,0,m),dimS2,BC_2L_global,SIcv(1,1,m),WIcv(1,0,m),COMM_BAR2)
     if (dimens == 3) call compact_parallel(S33   ,N33   ,N3,iB(3,1),NB3,ndL,cIcwCL(-ndL,0,m),dimS3,BC_3L_global,SIcw(1,1,m),WIcw(1,0,m),COMM_BAR3)
  end do
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp1CL ,dimS1,BC_1L_global,Sp1 ,Wp1 ,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp2CL ,dimS2,BC_2L_global,Sp2 ,Wp2 ,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp3CL ,dimS3,BC_3L_global,Sp3 ,Wp3 ,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu1CL ,dimS1,BC_1L_global,Su1 ,Wu1 ,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv2CL ,dimS2,BC_2L_global,Sv2 ,Wv2 ,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw3CL ,dimS3,BC_3L_global,Sw3 ,Ww3 ,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp11CL,dimS1,BC_1L_global,Sp11,Wp11,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp22CL,dimS2,BC_2L_global,Sp22,Wp22,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp33CL,dimS3,BC_3L_global,Sp33,Wp33,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu11CL,dimS1,BC_1L_global,Su11,Wu11,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv22CL,dimS2,BC_2L_global,Sv22,Wv22,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw33CL,dimS3,BC_3L_global,Sw33,Ww33,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CL,dimS1,BC_1L_global,SDu1,WDu1,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CL,dimS2,BC_2L_global,SDv2,WDv2,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CL,dimS3,BC_3L_global,SDw3,WDw3,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CLT,dimS1,BC_1L_global,SDu1T,WDu1T,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CLT,dimS2,BC_2L_global,SDv2T,WDv2T,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CLT,dimS3,BC_3L_global,SDw3T,WDw3T,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CL,dimS1,BC_1L_global,SGp1,WGp1,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CL,dimS2,BC_2L_global,SGp2,WGp2,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CL,dimS3,BC_3L_global,SGp3,WGp3,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CLT,dimS1,BC_1L_global,SGp1T,WGp1T,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CLT,dimS2,BC_2L_global,SGp2T,WGp2T,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CLT,dimS3,BC_3L_global,SGp3T,WGp3T,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cIupCL,dimS1,BC_1L_global,SIup,WIup,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cIvpCL,dimS2,BC_2L_global,SIvp,WIvp,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cIwpCL,dimS3,BC_3L_global,SIwp,WIwp,COMM_BAR3)
  !
  !IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cIpuCL,dimS1,BC_1L_global,SIpu,WIpu,COMM_BAR1)
  !IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cIpvCL,dimS2,BC_2L_global,SIpv,WIpv,COMM_BAR2)
  !IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cIpwCL,dimS3,BC_3L_global,SIpw,WIpw,COMM_BAR3)
  !
  !
  !DO m = 1, n_conc
  !   ! Anmerkung: BC_... ist so ok, es wird nur auf Periodizitaet getestet.
  !   IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc1CL (-ndL,0,m),dimS1,BC_1L_global,Sc1 (1,1,m),Wc1 (1,0,m),COMM_BAR1)
  !   IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc2CL (-ndL,0,m),dimS2,BC_2L_global,Sc2 (1,1,m),Wc2 (1,0,m),COMM_BAR2)
  !   IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc3CL (-ndL,0,m),dimS3,BC_3L_global,Sc3 (1,1,m),Wc3 (1,0,m),COMM_BAR3)
  !   
  !   IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc11CL(-ndL,0,m),dimS1,BC_1L_global,Sc11(1,1,m),Wc11(1,0,m),COMM_BAR1)
  !   IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc22CL(-ndL,0,m),dimS2,BC_2L_global,Sc22(1,1,m),Wc22(1,0,m),COMM_BAR2)
  !   IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc33CL(-ndL,0,m),dimS3,BC_3L_global,Sc33(1,1,m),Wc33(1,0,m),COMM_BAR3)
  !   
  !   IF ( NB1 .GT. 1 .OR. BC_1L_global == -1                   ) CALL compact_parallel(S11   ,N11   ,N1,iB(1,1),NB1,ndL,cIcuCL(-ndL,0,m),dimS1,BC_1L_global,SIcu(1,1,m),WIcu(1,0,m),COMM_BAR1)
  !   IF ( NB2 .GT. 1 .OR. BC_2L_global == -1                   ) CALL compact_parallel(S22   ,N22   ,N2,iB(2,1),NB2,ndL,cIcvCL(-ndL,0,m),dimS2,BC_2L_global,SIcv(1,1,m),WIcv(1,0,m),COMM_BAR2)
  !   IF ((NB3 .GT. 1 .OR. BC_3L_global == -1) .AND. dimens == 3) CALL compact_parallel(S33   ,N33   ,N3,iB(3,1),NB3,ndL,cIcwCL(-ndL,0,m),dimS3,BC_3L_global,SIcw(1,1,m),WIcw(1,0,m),COMM_BAR3)
  !END DO
  
  
  end subroutine init_compact
  
  
  
  
  
  
  
  
  
  
  
  subroutine diff_coeffs_compact(Nmax,NB,ndL,ndR,SS,NN,BCL,BCU,grid_type,abl,sym,transp_yes,ccL,ccL_LU,ccR)
  
  implicit none
  
  integer, intent(in   ) ::  Nmax, NB, ndL, ndR, BCL, BCU, SS, NN, grid_type, abl, sym
  logical, intent(in   ) ::  transp_yes
  real   , intent(out  ) ::  ccL   (-ndL:ndL   ,0:Nmax)
  real   , intent(out  ) ::  ccR   (-ndR:ndR   ,0:Nmax)
  real   , intent(out  ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  
  real                   ::  ccLT  (-ndL:ndL,-ndL:(Nmax+ndL))
  real                   ::  ccRT  (-ndR:ndR,-ndR:(Nmax+ndR))
  
  integer                ::  pivot(0:Nmax)
  
  integer                ::  i, ii
  real                   ::  dd
  
  
  ! TEST!!! Neumann-RB werden weiterhin separat und explizit gerechnet!
  !=====================================================================================================
  if (ndL == 2 .and. ndR == 3) then
     if ((grid_type == 0) .or. (grid_type == 1)) then
        if (abl == 0) then
           ! Filter ...
           do i = SS, NN
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
           end if
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        else if (abl == 1) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/1., 10., 20., 10., 1./)
              ccR(-ndR:ndR,i) = (/-1., -101., -425., 0., 425., 101., 1./) / 30.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              ccL(-ndL:ndL,SS+0) = (/0.,  0.,  1.,  0., 0./)
              ccL(-ndL:ndL,SS+1) = (/0.,  1.,  4.,  1., 0./)
              ccL(-ndL:ndL,SS+2) = (/1., 16., 36., 16., 1./)
              !ccL(-ndL:ndL,SS+0) = (/0.,  0.,  1.,  2., 0./)
              !ccL(-ndL:ndL,SS+1) = (/0.,  1.,  4.,  1., 0./)
              !ccL(-ndL:ndL,SS+2) = (/1., 16., 36., 16., 1./)
              
              ccR(-ndR:ndR,SS+0) = (/0.,   0.,   -3., -10.,  18., -6., 1./) / 12.
              ccR(-ndR:ndR,SS+1) = (/0.,   0.,   -3.,   0.,   3.,  0., 0./)
              ccR(-ndR:ndR,SS+2) = (/0., -25., -160.,   0., 160., 25., 0./) / 6.
              !ccR(-ndR:ndR,SS+0) = (/0.,   0.,    0., -5.,   4.,  1., 0./) / 2.
              !ccR(-ndR:ndR,SS+1) = (/0.,   0.,   -3.,  0.,   3.,  0., 0./)
              !ccR(-ndR:ndR,SS+2) = (/0., -25., -160.,  0., 160., 25., 0./) / 6.
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              ccL(-ndL:ndL,NN-0) = (/0.,  0.,  1.,  0., 0./)
              ccL(-ndL:ndL,NN-1) = (/0.,  1.,  4.,  1., 0./)
              ccL(-ndL:ndL,NN-2) = (/1., 16., 36., 16., 1./)
              !ccL(-ndL:ndL,NN-0) = (/0.,  2.,  1.,  0., 0./)
              !ccL(-ndL:ndL,NN-1) = (/0.,  1.,  4.,  1., 0./)
              !ccL(-ndL:ndL,NN-2) = (/1., 16., 36., 16., 1./)
              
              ccR(-ndR:ndR,NN-0) = (/-1.,   6.,  -18., 10.,   3.,  0., 0./) / 12.
              ccR(-ndR:ndR,NN-1) = (/ 0.,   0.,   -3.,  0.,   3.,  0., 0./)
              ccR(-ndR:ndR,NN-2) = (/ 0., -25., -160.,  0., 160., 25., 0./) / 6.
              !ccR(-ndR:ndR,NN-0) = (/0.,  -1.,   -4., 5.,   0.,  0., 0./) / 2.
              !ccR(-ndR:ndR,NN-1) = (/0.,   0.,   -3., 0.,   3.,  0., 0./)
              !ccR(-ndR:ndR,NN-2) = (/0., -25., -160., 0., 160., 25., 0./) / 6.
           end if
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        else if (abl == 2) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/43., 668., 1798., 668., 43./)
              ccR(-ndR:ndR,i) = (/79., 4671., 9585., -28670., 9585., 4671., 79./) / 9.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              ccL(-ndL:ndL,SS+0) = (/ 0.,   0.,    1.,   0.,  0./)
              ccL(-ndL:ndL,SS+1) = (/ 0.,   1.,   10.,   1.,  0./)
              ccL(-ndL:ndL,SS+2) = (/23., 688., 2358., 688., 23./)
              !ccL(-ndL:ndL,SS+0) = (/ 0.,   0.,    1.,  11.,  0./)
              !ccL(-ndL:ndL,SS+1) = (/ 0.,   1.,   10.,   1.,  0./)
              !ccL(-ndL:ndL,SS+2) = (/23., 688., 2358., 688., 23./)
              
              ccR(-ndR:ndR,SS+0) = (/0.,   0.,   11.,   -20.,    6.,   4., -1./) / 12.
              ccR(-ndR:ndR,SS+1) = (/0.,   0.,   12.,   -24.,   12.,   0.,  0./)
              ccR(-ndR:ndR,SS+2) = (/0., 465., 1920., -4770., 1920., 465.,  0./)
              !ccR(-ndR:ndR,SS+0) = (/0.,   0.,    0.,    13.,  -27.,  15., -1./)
              !ccR(-ndR:ndR,SS+1) = (/0.,   0.,   12.,   -24.,   12.,   0.,  0./)
              !ccR(-ndR:ndR,SS+2) = (/0., 465., 1920., -4770., 1920., 465.,  0./)
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              ccL(-ndL:ndL,NN-0) = (/ 0.,   0.,    1.,   0.,  0./)
              ccL(-ndL:ndL,NN-1) = (/ 0.,   1.,   10.,   1.,  0./)
              ccL(-ndL:ndL,NN-2) = (/23., 688., 2358., 688., 23./)
              !ccL(-ndL:ndL,NN-0) = (/ 0.,  11.,    1.,   0.,  0./)
              !ccL(-ndL:ndL,NN-1) = (/ 0.,   1.,   10.,   1.,  0./)
              !ccL(-ndL:ndL,NN-2) = (/23., 688., 2358., 688., 23./)
              
              ccR(-ndR:ndR,NN-0) = (/-1.,   4.,    6.,   -20.,   11.,   0., 0./) / 12.
              ccR(-ndR:ndR,NN-1) = (/ 0.,   0.,   12.,   -24.,   12.,   0., 0./)
              ccR(-ndR:ndR,NN-2) = (/ 0., 465., 1920., -4770., 1920., 465., 0./)
              !ccR(-ndR:ndR,NN-0) = (/-1.,  15.,  -27.,    13.,    0.,   0., 0./)
              !ccR(-ndR:ndR,NN-1) = (/ 0.,   0.,   12.,   -24.,   12.,   0., 0./)
              !ccR(-ndR:ndR,NN-2) = (/ 0., 465., 1920., -4770., 1920., 465., 0./)
           end if
        end if
     !**************************************************************************************************
     else if (grid_type == 2) then ! TEST!!! Sollte frueher oder spaeter mit grid_type == 3 verschmolzen werden (ndR verallgemeinern!)
        if (abl == 0) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/5., 60., 126., 60., 5./)
              ccR(-ndR:ndR,i) = (/0., 1., 45., 210., 210., 45., 1./) / 2.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              ccL(-ndL:ndL,SS+0) = (/0.,  0.,  0.,  1., 0./)
              ccL(-ndL:ndL,SS+1) = (/0.,  3., 10.,  3., 0./)
              !ccL(-ndL:ndL,SS+0) = (/0.,  0.,  1.,  0., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+1) = (/0.,  1.,  6.,  1., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+2) = (/1., 28., 70., 28., 1./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,SS+0) = (/0., 0., 0., -1.,  9.,  9., -1./) / 16.
              ccR(-ndR:ndR,SS+1) = (/0., 0., 1., 15., 15.,  1.,  0./) / 2.
              !ccR(-ndR:ndR,SS+0) = (/0., 0., 0.,  5., 15., -5.,  1./) / 16. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+1) = (/0., 0., 0.,  4.,  4.,  0.,  0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+2) = (/0., 0., 8., 56., 56.,  8.,  0./) ! TEST!!! ORIG!
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              ccL(-ndL:ndL,NN-0) = (/0.,  1.,  0.,  0., 0./)
              ccL(-ndL:ndL,NN-1) = (/0.,  3., 10.,  3., 0./)
              !ccL(-ndL:ndL,NN-0) = (/0.,  0.,  1.,  0., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-1) = (/0.,  1.,  6.,  1., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-2) = (/1., 28., 70., 28., 1./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,NN-0) = (/0., -1.,  9.,  9., -1., 0., 0./) / 16.
              ccR(-ndR:ndR,NN-1) = (/0.,  0.,  1., 15., 15., 1., 0./) / 2.
              !ccR(-ndR:ndR,NN-0) = (/0.,  1., -5., 15.,  5., 0., 0./) / 16. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-1) = (/0.,  0.,  0.,  4.,  4., 0., 0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-2) = (/0.,  0.,  8., 56., 56., 8., 0./) ! TEST!!! ORIG!
           end if
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        else if (abl == 1) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/9675., 193700., 577058., 193700., 9675./)
              ccR(-ndR:ndR,i) = (/0., -69049., -2525875., -6834250., 6834250., 2525875., 69049./) / 15.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              ccL(-ndL:ndL,SS+0) = (/   0.,     0.,     0.,     1.,   0./)
              ccL(-ndL:ndL,SS+1) = (/   0.,     9.,    62.,     9.,   0./)
              !ccL(-ndL:ndL,SS+0) = (/   0.,     0.,     1.,     0.,   0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+1) = (/   0.,     1.,    22.,     1.,   0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+2) = (/ 183., 12228., 51338., 12228., 183./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,SS+0) = (/0.,  0.,      0.,       1.,    -27.,    27., -1./) / 24.
              ccR(-ndR:ndR,SS+1) = (/0.,  0.,    -17.,    -189.,    189.,    17.,  0./) / 3.
              !ccR(-ndR:ndR,SS+0) = (/0.,  0.,      0.,     -23.,     21.,     3., -1./) / 24. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+1) = (/0.,  0.,      0.,     -24.,     24.,     0.,  0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+2) = (/0.,  0., -29360., -140400., 140400., 29360.,  0./) / 3. ! TEST!!! ORIG!
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              ccL(-ndL:ndL,NN-0) = (/   0.,      1.,      0.,      0.,    0./)
              ccL(-ndL:ndL,NN-1) = (/   0.,      9.,     62.,      9.,    0./)
              !ccL(-ndL:ndL,NN-0) = (/   0.,      0.,      1.,      0.,    0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-1) = (/   0.,      1.,     22.,      1.,    0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-2) = (/ 183.,  12228.,  51338.,  12228.,  183./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,NN-0) = (/0.,  1.,    -27.,      27.,     -1.,     0.,  0./) / 24.
              ccR(-ndR:ndR,NN-1) = (/0.,  0.,    -17.,    -189.,    189.,    17.,  0./) / 3.
              !ccR(-ndR:ndR,NN-0) = (/0.,  1.,     -3.,     -21.,     23.,     0.,  0./) / 24. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-1) = (/0.,  0.,      0.,     -24.,     24.,     0.,  0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-2) = (/0.,  0., -29360., -140400., 140400., 29360.,  0./) / 3. ! TEST!!! ORIG!
           end if
        end if
     else if (grid_type == 3) then
        if (abl == 0) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/5., 60., 126., 60., 5./)
              ccR(-ndR:ndR,i) = (/1., 45., 210., 210., 45., 1., 0./) / 2.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              ccL(-ndL:ndL,SS+0) = (/0.,  0.,  0.,  1., 0./)
              ccL(-ndL:ndL,SS+1) = (/0.,  3., 10.,  3., 0./)
              !ccL(-ndL:ndL,SS+0) = (/0.,  0.,  1.,  0., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+1) = (/0.,  1.,  6.,  1., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+2) = (/1., 28., 70., 28., 1./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,SS+0) = (/0., 0., -1.,  9.,  9., -1., 0./) / 16.
              ccR(-ndR:ndR,SS+1) = (/0., 1., 15., 15.,  1.,  0., 0./) / 2.
              !ccR(-ndR:ndR,SS+0) = (/0., 0.,  5., 15., -5.,  1., 0./) / 16. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+1) = (/0., 0.,  4.,  4.,  0.,  0., 0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+2) = (/0., 8., 56., 56.,  8.,  0., 0./) ! TEST!!! ORIG!
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              ccL(-ndL:ndL,NN-0) = (/0.,  1.,  0.,  0., 0./)
              ccL(-ndL:ndL,NN-1) = (/0.,  3., 10.,  3., 0./)
              !ccL(-ndL:ndL,NN-0) = (/0.,  0.,  1.,  0., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-1) = (/0.,  1.,  6.,  1., 0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-2) = (/1., 28., 70., 28., 1./) ! TEST!!! ORIG!
              
              ccR(-ndR:ndR,NN-0) = (/-1.,  9.,  9.,  -1., 0., 0., 0./) / 16.
              ccR(-ndR:ndR,NN-1) = (/ 0.,  1., 15.,  15., 1., 0., 0./) / 2.
              !ccR(-ndR:ndR,NN-0) = (/ 1., -5., 15.,   5., 0., 0., 0./) / 16. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-1) = (/ 0.,  0.,  4.,   4., 0., 0., 0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-2) = (/ 0.,  8., 56.,  56., 8., 0., 0./) ! TEST!!! ORIG!
           end if
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        else if (abl == 1) then
           do i = SS, NN
              ccL(-ndL:ndL,i) = (/9675., 193700., 577058., 193700., 9675./)
              ccR(-ndR:ndR,i) = (/-69049., -2525875., -6834250., 6834250., 2525875., 69049., 0./) / 15.
           end do
           !--------------------------------------------------------------------------------------------
           if (BCL >= 1) then
              !ccL(-ndL:ndL,SS+0) = (/   0.,     0.,     0.,     1.,   0./) ! Brueger et. al.
              ccL(-ndL:ndL,SS+0) = (/   0.,     0.,     0.,     1.,   0./)
              ccL(-ndL:ndL,SS+1) = (/   0.,     9.,    62.,     9.,   0./)
              !ccL(-ndL:ndL,SS+0) = (/   0.,     0.,     1.,     0.,   0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+1) = (/   0.,     1.,    22.,     1.,   0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,SS+2) = (/ 183., 12228., 51338., 12228., 183./) ! TEST!!! ORIG!
              
              !ccR(-ndR:ndR,SS+0) = (/ 0.,      0.,       1.,    -27.,    27., -1., 0./) / 24. ! Brueger et. al.
              ccR(-ndR:ndR,SS+0) = (/ 0.,      0.,       1.,    -27.,    27., -1., 0./) / 24.
              ccR(-ndR:ndR,SS+1) = (/ 0.,    -17.,    -189.,    189.,    17.,  0., 0./) / 3.
              !ccR(-ndR:ndR,SS+0) = (/ 0.,      0.,     -23.,     21.,     3., -1., 0./) / 24. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+1) = (/ 0.,      0.,     -24.,     24.,     0.,  0., 0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,SS+2) = (/ 0., -29360., -140400., 140400., 29360.,  0., 0./) / 3. ! TEST!!! ORIG!
           end if
           !--------------------------------------------------------------------------------------------
           if (BCU >= 1) then
              !ccL(-ndL:ndL,NN-0) = (/   0.,      1.,      0.,      0.,    0./) ! Brueger et. al.
              ccL(-ndL:ndL,NN-0) = (/   0.,      1.,      0.,      0.,    0./)
              ccL(-ndL:ndL,NN-1) = (/   0.,      9.,     62.,      9.,    0./)
              !ccL(-ndL:ndL,NN-0) = (/   0.,      0.,      1.,      0.,    0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-1) = (/   0.,      1.,     22.,      1.,    0./) ! TEST!!! ORIG!
              !ccL(-ndL:ndL,NN-2) = (/ 183.,  12228.,  51338.,  12228.,  183./) ! TEST!!! ORIG!
              
              !ccR(-ndR:ndR,NN-0) = (/ 1.,    -27.,      27.,     -1.,     0.,  0., 0./) / 24. ! Brueger et. al.
              ccR(-ndR:ndR,NN-0) = (/ 1.,    -27.,      27.,     -1.,     0.,  0., 0./) / 24.
              ccR(-ndR:ndR,NN-1) = (/ 0.,    -17.,    -189.,    189.,    17.,  0., 0./) / 3.
              !ccR(-ndR:ndR,NN-0) = (/ 1.,     -3.,     -21.,     23.,     0.,  0., 0./) / 24. ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-1) = (/ 0.,      0.,     -24.,     24.,     0.,  0., 0./) ! TEST!!! ORIG!
              !ccR(-ndR:ndR,NN-2) = (/ 0., -29360., -140400., 140400., 29360.,  0., 0./) / 3. ! TEST!!! ORIG!
           end if
        end if
     end if
  !=====================================================================================================
  else if (ndL == 3 .and. ndR == 3) then
     do i = SS, NN
        ccL(-ndL:ndL,i) = (/ 1. ,  36., 225. ,400.,225. ,36. ,1. /) / 400.
        ccR(-ndR:ndR,i) = (/-49.,-924.,-2625.,0.  ,2625.,924.,49./) / 4000.
     end do
     !--------------------------------------------------------------------------------------------------
     if (BCL >= 1) then
        ccL(-ndL:ndL,SS+0) = (/0.,0.,0. ,1. ,2. ,0.,0./)
        ccL(-ndL:ndL,SS+1) = (/0.,0.,1. ,4. ,1. ,0.,0./) / 4.
        ccL(-ndL:ndL,SS+2) = (/0.,1.,16.,36.,16.,1.,0./) / 36.
        
        ccR(-ndR:ndR,SS+0) = (/0., 0. , 0.  ,-5. ,4.  , 1. ,0./) / 2.
        ccR(-ndR:ndR,SS+1) = (/0., 0. ,-3.  , 0. ,3.  , 0. ,0./) / 4.
        ccR(-ndR:ndR,SS+2) = (/0.,-25.,-160., 0. ,160., 25.,0./) / 216.
     end if
     !--------------------------------------------------------------------------------------------------
     if (BCU >= 1) then
        ccL(-ndL:ndL,NN-0) = (/0.,0.,2. ,1. ,0. ,0.,0./)
        ccL(-ndL:ndL,NN-1) = (/0.,0.,1. ,4. ,1. ,0.,0./) / 4.
        ccL(-ndL:ndL,NN-2) = (/0.,1.,16.,36.,16.,1.,0./) / 36.
        
        ccR(-ndR:ndR,NN-0) = (/ 0.,-1. ,-4.  ,5. ,0.  ,0. ,0./) / 2.
        ccR(-ndR:ndR,NN-1) = (/ 0., 0. ,-3.  ,0. ,3.  ,0. ,0./) / 4.
        ccR(-ndR:ndR,NN-2) = (/ 0.,-25.,-160.,0. ,160.,25.,0./) / 216.
     end if
  !=====================================================================================================
  else if (ndL == 1 .and. ndR == 1) then
     do i = SS, NN
        ccL(-ndL:ndL,i) = (/ 1.,4.,1./) / 4.
        ccR(-ndR:ndR,i) = (/-3.,0.,3./) / 4.
     end do
     !--------------------------------------------------------------------------------------------------
     if (BCL >= 1) then
        ccL(-ndL:ndL,SS+0) = (/0., 1.,0./) ! TEST!!! Schlechte Wahl ...
        
        ccR(-ndR:ndR,SS+0) = (/0.,-1.,1./)
     end if
     !--------------------------------------------------------------------------------------------------
     if (BCU >= 1) then
        ccL(-ndL:ndL,NN-0) = (/ 0.,1.,0./)
        
        ccR(-ndR:ndR,NN-0) = (/-1.,1.,0./)
     end if
  !=====================================================================================================
  else
     if (rank == 0) write(*,*) 'ERROR! No coefficients specified for this choice of ndL, ndR.'
     call MPI_FINALIZE(merror)
     stop
  end if
  !=====================================================================================================
  
  
  ! TEST!!! nochmals richtig durchtesten fuer ndR /= 3 und ndL /= 2 ...
  !=====================================================================================================
  !=== Symmetrie =======================================================================================
  !=====================================================================================================
  if (BCL == -2) then
     if (grid_type == 0) then
        do i = 1, ndL
           do ii = 1, 1-i+ndL
              ccL(1-i+ii,i) = ccL(1-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           end do
        end do
        do i = 1, ndR
           do ii = 1, 1-i+ndR
              ccR(1-i+ii,i) = ccR(1-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 1) then
        do i = 1, ndL
           do ii = 1, 1-i+ndL
              ccL(0-i+ii,i) = ccL(0-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           end do
        end do
        do i = 1, ndR
           do ii = 1, 1-i+ndR
              ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 2) then
        do i = 1, ndL
           do ii = 1, 1-i+ndL
              ccL(0-i+ii,i) = ccL(0-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           end do
        end do
        do i = 1, ndR
           do ii = 1, 1-i+ndR
              ccR(1-i+ii,i) = ccR(1-i+ii,i) + ccR(1-i-ii,i)*sym ! TEST!!! Nur wegen aktuellem Shift bei grid_type == 2
              ccR(1-i-ii,i) = 0.
              !ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(0-i-ii,i)*sym
              !ccR(0-i-ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 3) then
        do i = 1, ndL
           do ii = 1, 1-i+ndL
              ccL(1-i+ii,i) = ccL(1-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           end do
        end do
        do i = 1, ndR
           do ii = 1, 1-i+ndR
              ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           end do
        end do
     end if
  end if
  !=====================================================================================================
  if (BCU == -2) then
     if (grid_type == 0) then
        do i = Nmax-ndL+1, Nmax-0
           do ii = 1, i+ndL-(Nmax-0)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i  +ii,i)*sym*(-1.)**abl
              ccL(Nmax-i  +ii,i) = 0.
           end do
        end do
        do i = Nmax-ndR+1, Nmax-0
           do ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i  +ii,i)*sym
              ccR(Nmax-i  +ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 1) then
        do i = Nmax-ndL, Nmax-1
           do ii = 1, i+ndL-(Nmax-1)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i-1+ii,i)*sym*(-1.)**abl
              ccL(Nmax-i-1+ii,i) = 0.
           end do
        end do
        do i = Nmax-ndR+0, Nmax-1
           do ii = 1, i+ndR-(Nmax-1)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i-1+ii,i)*sym
              ccR(Nmax-i-1+ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 2) then
        do i = Nmax-ndL+0, Nmax-1
           do ii = 1, i+ndL-(Nmax-1)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i-1+ii,i)*sym*(-1.)**abl
              ccL(Nmax-i-1+ii,i) = 0.
           end do
        end do
        do i = Nmax-ndR+1, Nmax-0
           do ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i  +ii,i)*sym! TEST!!! Nur wegen aktuellem Shift bei grid_type == 2
              ccR(Nmax-i  +ii,i) = 0.
        !      ccR(Nmax-i-1-ii,i) = ccR(Nmax-i-1-ii,i) + ccR(Nmax-i-1+ii,i)*sym
        !      ccR(Nmax-i-1+ii,i) = 0.
           end do
        end do
     end if
     !--------------------------------------------------------------------------------------------------
     if (grid_type == 3) then
        do i = Nmax-ndL+1, Nmax-0
           do ii = 1, i+ndL-(Nmax-0)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i  +ii,i)*sym*(-1.)**abl
              ccL(Nmax-i  +ii,i) = 0.
           end do
        end do
        do i = Nmax-ndR+1, Nmax-0
           do ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i-1+ii,i)*sym
              ccR(Nmax-i-1+ii,i) = 0.
           end do
        end do
     end if
  end if
  !=====================================================================================================
  
  
  
  !=====================================================================================================
  !=== Transposition ===================================================================================
  !=====================================================================================================
  if (transp_yes) then
     ccLT = 0.
     ccRT = 0.
     ccLT(-ndL:ndL,SS:NN) = ccL(-ndL:ndL,SS:NN)
     ccRT(-ndR:ndR,SS:NN) = ccR(-ndR:ndR,SS:NN)
     
     !*************************************
     ! TEST!!! Dieser Abschnitt ist nur erlaubt, wenn die Stencil-Koeffizienten im
     !         Innern (bzw. genauer gesagt: and den Blockgrenzen) identisch sind!
     !         (Ansonsten m�sste per MPI ausgetauscht werden ...)
     !         (gilt hier auf f�r Periodizit�t!)
     if (BCL == 0 .or. BCL == -1) then
        do ii = -ndL, ndL
           ccLT(ii,-ndL:(SS-1)) = ccL(ii,SS)
        end do
        do ii = -ndR, ndR
           ccRT(ii,-ndR:(SS-1)) = ccR(ii,SS)
        end do
     end if
     if (BCU == 0 .or. BCU == -1) then
        do ii = -ndL, ndL
           ccLT(ii,(NN+1):(NN+ndL)) = ccL(ii,NN)
        end do
        do ii = -ndR, ndR
           ccRT(ii,(NN+1):(NN+ndR)) = ccR(ii,NN)
        end do
     end if
     !*************************************
     ccL = 0.
     ccR = 0.
     do i = SS, NN ! ok ...
        do ii = -ndL, ndL
           ccL(ii,i) = ccLT(-ii,i+ii)
        end do
     end do
     !DO i = SS, NN ! TEST!!! ... Intervall SS:NN nicht ok! 
     do i = 0, Nmax ! TEST!!! R�nder sollten noch auf Null gesetzt werden!
        do ii = -ndR, ndR
           ccR(ii,i) = ccRT(-ii,i+ii)
        end do
     end do
     
  end if
  !=====================================================================================================
  
  
  
  !=====================================================================================================
  !=== LU-Zerlegung ====================================================================================
  !=====================================================================================================
  ccL_LU = 0.
  
  !--- Transponierte Koeffizienten ---
  do i = SS, NN
     do ii = -ndL, ndL
        if ((i+ii >= SS) .and. (i+ii <= NN)) ccL_LU(2*ndL-ii+1,i+ii) = ccL(ii,i)
     end do
  end do
  
  !--- LU-Zerlegung ---
  if (NB == 1 .and. BCL /= -1) then
     call myDGBTRF(NN-SS+1      ,ndL,ccL_LU(1,SS    ),3*ndL+1,pivot(SS    ))
  else
     call myDGBTRF(NN-SS+1-2*ndL,ndL,ccL_LU(1,SS+ndL),3*ndL+1,pivot(SS+ndL))
  end if
  !=====================================================================================================
  
  
  end subroutine diff_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  subroutine myDGBTRF(N,ndL,AB,LDAB,pivot)
  
  implicit none
  
  integer, intent(in   )  ::  ndL, LDAB, N
  integer, intent(out  )  ::  pivot( * )
  real   , intent(inout)  ::  AB( LDAB, * )
  integer, parameter      ::  NBMAX = 64
  integer, parameter      ::  LDWORK = NBMAX+1
  integer                 ::  I, I2, I3, II, IP, J, J2, J3, JJ, JM, JP, JU, K2, KM, KV, NW
  real                    ::  TEMP
  real                    ::  WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX )
  integer                 ::  IDAMAX
  
  external IDAMAX
  external DCOPY, DGEMM, DGER, DLASWP, DSCAL, DSWAP, DTRSM
  
  
  KV = ndL + ndL
  
  do J = ndL + 2, MIN( KV, N )
     do I = KV - J + 2, ndL
        AB( I, J ) = 0.
     end do
  end do
  
  JU = 1
  
  do J = 1, N
     I2 = MIN(ndL-1,N-J)
     I3 = MIN(1,N-J-ndL+1)
     
     if(J+KV <= N) then
        do I = 1, ndL
           AB(I,J+KV) = 0.
        end do
     end if
     
     KM = MIN(ndL,N-J)
     
     !--- Pivoting ---
     !JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
     JP = 1 ! TEST!!!
     
     pivot(J) = JP
     
     if (AB(KV+JP,J) /= 0.) then
        JU = MAX( JU, MIN(J+ndL+JP-1,N))
        if (JP /= 1) then
           if (JP-1 < ndL) then
              call DSWAP( 1, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
           else
              call DSWAP( 0, AB( KV+1, J ), LDAB-1,WORK31( JP-ndL, 1 ), LDWORK )
              call DSWAP( 1, AB( KV+1, J ), LDAB-1,AB( KV+JP, J ), LDAB-1 )
           end if
        end if
        call DSCAL( KM, 1. / AB( KV+1, J ), AB( KV+2, J ),1 )
        JM = MIN( JU, J )
        if (JM > J) call DGER( KM, JM-J, -1., AB( KV+2, J ), 1, AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), LDAB-1 )
     end if
     NW = MIN( 1, I3 )
     if(NW > 0) call DCOPY( NW, AB( KV+ndL+1, J ), 1, WORK31( 1, 1 ), 1 )
     
     
     
     if (J+1 <= N) then
        
        J2 = MIN( JU-J+1, KV ) - 1
        J3 = MAX( 0, JU-J-KV+1 )
        
        call DLASWP( J2, AB( KV, J+1 ), LDAB-1, 1, 1, pivot( J ), 1 )
        
        pivot( J ) = pivot( J ) + J - 1
        
        K2 = J + J2
        do I = 1, J3
           JJ = K2 + I
           do II = J + I - 1, J
              IP = pivot( II )
              if (IP /= II) then
                 TEMP = AB( KV+1+II-JJ, JJ )
                 AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                 AB( KV+1+IP-JJ, JJ ) = TEMP
              end if
           end do
        end do
        
        if (J2 > 0) then
           call DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J2, 1., AB( KV+1, J ), LDAB-1, AB( KV, J+1 ), LDAB-1 )
           if (I2 > 0) call DGEMM( 'No transpose', 'No transpose', I2, J2, 1, -1., AB( KV+2, J ), LDAB-1, AB( KV, J+1 ), LDAB-1, 1., AB( KV+1, J+1 ), LDAB-1 )
           if (I3 > 0) call DGEMM( 'No transpose', 'No transpose', I3, J2, 1, -1., WORK31, LDWORK, AB( KV, J+1 ), LDAB-1, 1., AB( KV+ndL, J+1 ), LDAB-1 )
        end if
        
        if (J3 > 0) then
           do JJ = 1, J3
              do II = JJ, 1
                 WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
              end do
           end do
           
           call DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J3, 1., AB( KV+1, J ), LDAB-1, WORK13, LDWORK )
           if( I2>0 ) call DGEMM( 'No transpose', 'No transpose', I2, J3, 1, -1., AB( KV+2, J ), LDAB-1, WORK13, LDWORK, 1., AB( 2, J+KV ), LDAB-1 )
           if( I3>0 ) call DGEMM( 'No transpose', 'No transpose', I3, J3, 1, -1., WORK31, LDWORK, WORK13, LDWORK, 1., AB( 1+ndL, J+KV ), LDAB-1 )
           
           do JJ = 1, J3
              do II = JJ, 1
                 AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
              end do
           end do
        end if
        
     else
        pivot( J ) = pivot( J ) + J - 1
     end if
     
     JP = pivot( J ) - J + 1
     if (JP /= 1) then
        if (JP-1 < ndL) then
           call DSWAP( 0, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
        else
           call DSWAP( 0, AB( KV+1, J ), LDAB-1, WORK31( JP-ndL, 1 ), LDWORK )
        end if
     end if
     NW = MIN( I3, 1 )
     if( NW>0 ) call DCOPY( NW, WORK31( 1, 1 ), 1, AB( KV+ndL+1, J ), 1 )
     
  end do
  
  
  end subroutine myDGBTRF
  
  
  
  
  
  
  
  
  
  
  
  subroutine compact_parallel(SS,NN,Nmax,iB,NB,ndL,ccL,dimS,BCL,Schur,WW,COMM)
  
  implicit none
  
  integer, intent(in)   ::  SS, NN
  integer, intent(in)   ::  Nmax
  integer, intent(in)   ::  iB, NB
  integer, intent(in)   ::  ndL
  real   , intent(in)   ::  ccL(-ndL:ndL,0:Nmax)
  integer, intent(in)   ::  dimS
  integer, intent(in)   ::  BCL
  real   , intent(out)  ::  Schur(1:dimS,1:dimS)
  real   , intent(out)  ::  WW(1:2*ndL,0:Nmax)
  integer, intent(in)   ::  COMM
  
  real                  ::  VV(1:2*ndL,0:Nmax)
  
  real   , allocatable  ::  VV1(:,:), VV2(:,:)
  real   , allocatable  ::  WW1(:,:), WW2(:,:)
  real                  ::  Schur_global(1:dimS,1:dimS)
  
  real                  ::  ccL_ext(0:3*ndL,SS:NN)
  integer               ::  pivot1(0:Nmax), pivot2(1:dimS)
  integer               ::  dimG, dimA1(1:NB), dimA2(1:NB), dimA, i0
  integer               ::  i, j, ii
  integer               ::  info
  
  
  !-------------------------------------------------------------------------------------------------------------!
  !  | H  G | | u |   | f |          | H  G          | | u |   | f          |                                   !
  !  |      | |   | = |   |   <==>   |               | |   | = |            |                                   !
  !  | D  B | | p |   | b |          | 0  B-DH^(-1)G | | p |   | b-DH^(-1)f |                                   !
  !                                                                                                             !
  !                                                                                                             !
  !  M :=  H^(-1)                                                                                               !
  !  W :=  DM            ("WW")                                                                                 !
  !  V :=  G             ("VV")                                                                                 !
  !  S := (B-WV)^(-1)    ("Schur")                                                                              !
  !                                                                                                             !
  !                                                                                                             !
  !  | u |   | M(f-VS(b-Wf)) |   | M(I-VSW) -MVS | | f |                                                        !
  !  |   | = |               | = |               | |   |                                                        !
  !  | p |   |      S(b-Wf)  |   |     -SW     S | | b |                                                        !
  !                                                                                                             !
  !                                                                                                             !
  !  | M(I-VSW) -MVS |   | M  0 |  | I  0 |   | V  0 | | S   S | | W  0 |                                       !
  !  |               | = |      | (|      | - |      | |       | |      |)                                      !
  !  |     -SW     S |   | 0  I |  | 0  0 |   | 0  I | | S  -S | | 0  I |                                       !
  !                                                                                                             !
  !                                                                                                             !
  !                                                                                                             !
  !  | M(I-VSW) -MVS |t    | I  0 |    | Wt  0 | | St   St | | Vt  0 |  | Mt  0 |   | (I-Wt St Vt)Mt  -Wt St |  !
  !  |               |  = (|      | -  |       | |         | |       |) |       | = |                        |  !
  !  |     -SW     S |     | 0  0 |    | 0   I | | St  -St | | 0   I |  | 0   I |   |      -St Vt Mt      St |  !
  !                                                                                                             !
  !                                                                                                             !
  !-------------------------------------------------------------------------------------------------------------!
  
  
  dimG = NN-SS+1-2*ndL
  
  VV = 0.
  WW = 0.
  
  Schur        = 0.
  Schur_global = 0.
  
  
  i = 2*ndL*(iB-1)
  
  do ii = 1, ndL
     !--- Zentraler Block (global) ---
     Schur(i+ii    ,(i+1    ):(i+1*ndL)) = ccL((1-ii):(ndL-ii),SS+ii-1  )
     Schur(i+ii+ndL,(i+1+ndL):(i+2*ndL)) = ccL((1-ii):(ndL-ii),NN+ii-ndL)
     
     !--- restliche Bloecke (global) ---
     if (iB /= 1                 ) Schur(i+ii,(i   -ndL+ii):i   ) = ccL(-ndL:-ii,SS+ii-1)
     if (iB == 1  .and. BCL == -1) Schur(i+ii,(dimS-ndL+ii):dimS) = ccL(-ndL:-ii,SS+ii-1)
     
     if (iB /= NB                ) Schur(i+ii+ndL,(i+2*ndL+1):(i+2*ndL+ii)) = ccL((ndL-ii+1):ndL,NN+ii-ndL)
     if (iB == NB .and. BCL == -1) Schur(i+ii+ndL,         1 :         ii ) = ccL((ndL-ii+1):ndL,NN+ii-ndL)
     
     !--- restliche Bloecke (lokal) ---
     VV(    ii : ndL    ,(SS+ndL +ii-1)) = ccL(     -ndL :-ii,SS+ii-1+ndL)
     VV((1+ndL):(ii+ndL),(SS+dimG+ii-1)) = ccL((1+ndL-ii):ndL,NN+ii-2*ndL)
     
     WW(ii    ,(SS+ndL      ):(SS+ndL-1+ii  )) = ccL((ndL-ii+1):ndL,SS+ii-1  )
     WW(ii+ndL,(SS+dimG+ii-1):(SS+ndL-1+dimG)) = ccL(  -ndL    :-ii,NN+ii-ndL)
  end do
  
  
  !--- Bandmatrix erweitern fuer LAPACK ---
  ccL_ext = 0.
  do i = SS, NN
     ccL_ext(ndL:3*ndL,i) = ccL(-ndL:ndL,i)
  end do
  
  
  !--- Bandloeser fuer WW ---
  call DGBTRF(dimG,dimG,ndL,ndL,ccL_ext(0,SS+ndL),3*ndL+1,pivot1,info)
  do i = 1, 2*ndL
     call DGBTRS('N',dimG,ndL,ndL,1,ccL_ext(0,SS+ndL),3*ndL+1,pivot1,WW(i,(SS+ndL):(SS+ndL+dimG-1):1),dimG,info)
  end do
  
  
  !--- Globale Dimensionen und Positionen bestimmen ---
  dimA1     = 0
  dimA1(iB) = dimG
  
  call MPI_ALLREDUCE(dimA1,dimA2,NB,MPI_INTEGER,MPI_SUM,COMM,merror)
  
  dimA = 0
  do i = 1, NB
     dimA = dimA + dimA2(i)
  end do
  
  i0   = 0
  do i = 1, iB-1
     i0   = i0   + dimA2(i)
  end do
  
  allocate(VV1(1:dimS,1:dimA))
  allocate(VV2(1:dimS,1:dimA))
  allocate(WW1(1:dimS,1:dimA))
  allocate(WW2(1:dimS,1:dimA))
  
  VV1 = 0.
  VV2 = 0.
  WW1 = 0.
  WW2 = 0.
  
  
  !--- Umspeichern auf lokalen Array ---
  do i = 1, dimG
     VV1((1+2*ndL*(iB-1)):2*ndL*iB,i+i0) = VV(1:2*ndL,(SS+ndL+i-1))
     WW1((1+2*ndL*(iB-1)):2*ndL*iB,i+i0) = WW(1:2*ndL,(SS+ndL+i-1))
  end do
  
  
  !--- Matrizen global zusammenbauen ---
  call MPI_ALLREDUCE(VV1,VV2,dimS*dimA,MPI_REAL8,MPI_SUM,COMM,merror)
  call MPI_ALLREDUCE(WW1,WW2,dimS*dimA,MPI_REAL8,MPI_SUM,COMM,merror)
  call MPI_ALLREDUCE(Schur,Schur_global,dimS*dimS,MPI_REAL8,MPI_SUM,COMM,merror)
  
  
  !--- Schur-Komplement bilden ---
  do j = 1, dimS
     do i = 1, dimS
        do ii = 1, dimA
           Schur_global(i,j) = Schur_global(i,j) - WW2(i,ii)*VV2(j,ii)
        end do
     end do
  end do
  
  deallocate(VV1)
  deallocate(VV2)
  deallocate(WW1)
  deallocate(WW2)
  
  
  !--- Einheitsmatrix als RHS ---
  Schur = 0.
  do i = 1, dimS
     Schur(i,i) = 1.
  end do
  
  
  !--- Inverse des Schur-Komplements ---
  ! TEST!! Hier gibt es in Lapack einen Floating Point exception ..
  !CALL DGESV(dimS,dimS,Schur_global,dimS,pivot2,Schur,dimS,info)
  call matrix_invert(dimS,Schur_global,Schur)
  
  
  end subroutine compact_parallel
  
  
  
  
  
  
  
  
  
  
  
  subroutine test_coeffs_compact()

  implicit none
  
  integer                ::  h, i, j, k
  integer                ::  abl, grid
  real                   ::  dd1, max_diff, max_diff_global
  real                   ::  omega, phase, xref
  real                   ::  pi
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  omega = 4.
  phase = 1.
  
  !===========================================================================================================
  do h = 1, 8
     
     if (h == 1) then
        abl  = 0
        grid = 2
     else if (h == 2) then
        abl  = 0
        grid = 3
     else if (h == 3) then
        abl  = 1
        grid = 0
     else if (h == 4) then
        abl  = 1
        grid = 1
     else if (h == 5) then
        abl  = 1
        grid = 2
     else if (h == 6) then
        abl  = 1
        grid = 3
     else if (h == 7) then
        abl  = 2
        grid = 0
     else if (h == 8) then
        abl  = 2
        grid = 1
     end if
     
     
     xref = phase*0.25*L1/omega
     
     if (BC_1L_global == -2 .or. BC_1U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L1
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (grid == 0 .or. grid == 2) pp(i,j,k) =  COS(dd1*(x1p(i)-xref))
              
              if (grid == 0 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x1p(i)-xref))
              if (grid == 0 .and. abl == 2) rh(i,j,k) = -COS(dd1*(x1p(i)-xref))
              
              if (grid == 3 .and. abl == 0) rh(i,j,k) =  SIN(dd1*(x1p(i)-xref))
              if (grid == 3 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x1p(i)-xref))
           end do
        end do
     end do
     do k = S31B, N31B
        do j = S21B, N21B
           do i = S11B, N11B
              if (grid == 1 .or. grid == 3) pp(i,j,k) =  SIN(dd1*(x1u(i)-xref))
              
              if (grid == 1 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x1u(i)-xref))
              if (grid == 1 .and. abl == 2) rh(i,j,k) = -SIN(dd1*(x1u(i)-xref))
              
              if (grid == 2 .and. abl == 0) rh(i,j,k) =  COS(dd1*(x1u(i)-xref))
              if (grid == 2 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x1u(i)-xref))
           end do
        end do
     end do
     
     pp = pp * dd1**(-abl)
     
     if (grid == 0 .or. grid == 2) call exchange(1,0,pp)
     if (grid == 1 .or. grid == 3) call exchange(1,1,pp)
     
     if (abl == 0 .and. grid == 2) call apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,pp,Ap)
     if (abl == 0 .and. grid == 3) call apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,pp,Ap)
     
     if (abl == 1 .and. grid == 0) call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ap)
     if (abl == 1 .and. grid == 1) call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ap)
     if (abl == 1 .and. grid == 2) call apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,pp,Ap)
     if (abl == 1 .and. grid == 3) call apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,pp,Ap)
     
     if (abl == 2 .and. grid == 0) call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ap)
     if (abl == 2 .and. grid == 0) call apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,pp,rr)
     
     if (abl == 2 .and. grid == 1) call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ap)
     if (abl == 2 .and. grid == 1) call apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,pp,rr)
     
     
     max_diff = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (abl == 1 .and. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx1DM(i)
              
              if (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              if (abl == 1 .and. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx1pM(i)
              if (abl == 2 .and. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx1pM(i)**2 + Ap(i,j,k)*ddx1pM(i)
              
              if (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
           
        end do
     end do
     
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (abl == 1 .and. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx1uM(i)
              if (abl == 1 .and. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx1GM(i)
              if (abl == 2 .and. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx1uM(i)**2 + Ap(i,j,k)*ddx1uM(i)
              
              if (grid == 1 .or. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     call MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     if (rank == 0) write(*,*) max_diff_global
     
     
     if (1 == 1 .and. h == 8) then
     if (grid == 0) then
        j = N2p
        k = N3p
        do i = S12, N12
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(55,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 2 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(56,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 3 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(57,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 4 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(58,*) x1p(i), Ap(i,j,k)
        end do
     else if (grid == 3) then
        j = N2p
        k = N3p
        do i = S1p, N1p
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(55,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 2 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(56,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 3 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(57,*) x1p(i), Ap(i,j,k)
           if (iB(1,1) == 4 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(58,*) x1p(i), Ap(i,j,k)
        end do
     else
        j = N21
        k = N31
        do i = S11, N11
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(55,*) x1u(i), Ap(i,j,k)
           if (iB(1,1) == 2 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(56,*) x1u(i), Ap(i,j,k)
           if (iB(1,1) == 3 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(57,*) x1u(i), Ap(i,j,k)
           if (iB(1,1) == 4 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(58,*) x1u(i), Ap(i,j,k)
        end do
     end if
     end if
  end do
  
  !===========================================================================================================
  if (rank == 0) write(*,*)
  
  do h = 1, 8
     
     if (h == 1) then
        abl  = 0
        grid = 2
     else if (h == 2) then
        abl  = 0
        grid = 3
     else if (h == 3) then
        abl  = 1
        grid = 0
     else if (h == 4) then
        abl  = 1
        grid = 1
     else if (h == 5) then
        abl  = 1
        grid = 2
     else if (h == 6) then
        abl  = 1
        grid = 3
     else if (h == 7) then
        abl  = 2
        grid = 0
     else if (h == 8) then
        abl  = 2
        grid = 1
     end if
     
     
     xref = phase*0.25*L2/omega
     
     if (BC_2L_global == -2 .or. BC_2U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L2
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (grid == 0 .or. grid == 2) pp(i,j,k) =  COS(dd1*(x2p(j)-xref))
              
              if (grid == 0 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x2p(j)-xref))
              if (grid == 0 .and. abl == 2) rh(i,j,k) = -COS(dd1*(x2p(j)-xref))
              
              if (grid == 3 .and. abl == 0) rh(i,j,k) =  SIN(dd1*(x2p(j)-xref))
              if (grid == 3 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x2p(j)-xref))
           end do
        end do
     end do
     do k = S32B, N32B
        do j = S22B, N22B
           do i = S12B, N12B
              if (grid == 1 .or. grid == 3) pp(i,j,k) =  SIN(dd1*(x2v(j)-xref))
              
              if (grid == 1 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x2v(j)-xref))
              if (grid == 1 .and. abl == 2) rh(i,j,k) = -SIN(dd1*(x2v(j)-xref))
              
              if (grid == 2 .and. abl == 0) rh(i,j,k) =  COS(dd1*(x2v(j)-xref))
              if (grid == 2 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x2v(j)-xref))
           end do
        end do
     end do
     
     pp = pp * dd1**(-abl)
     
     if (grid == 0 .or. grid == 2) call exchange(2,0,pp)
     if (grid == 1 .or. grid == 3) call exchange(2,2,pp)
     
     if (abl == 0 .and. grid == 2) call apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,pp,Ap)
     if (abl == 0 .and. grid == 3) call apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,pp,Ap)
     
     if (abl == 1 .and. grid == 0) call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ap)
     if (abl == 1 .and. grid == 1) call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ap)
     if (abl == 1 .and. grid == 2) call apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,pp,Ap)
     if (abl == 1 .and. grid == 3) call apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,pp,Ap)
     
     if (abl == 2 .and. grid == 0) call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ap)
     if (abl == 2 .and. grid == 0) call apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,pp,rr)
     
     if (abl == 2 .and. grid == 1) call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ap)
     if (abl == 2 .and. grid == 1) call apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,pp,rr)
     
     
     max_diff = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (abl == 1 .and. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx2DM(j)
              
              if (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (abl == 1 .and. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx2pM(j)
              if (abl == 2 .and. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx2pM(j)**2 + Ap(i,j,k)*ddx2pM(j)
              
              if (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     do k = S32, N32
        do j = S22, N22
           do i = S12, N12
              if (abl == 1 .and. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx2vM(j)
              if (abl == 1 .and. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx2GM(j)
              if (abl == 2 .and. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx2vM(j)**2 + Ap(i,j,k)*ddx2vM(j)
              
              if (grid == 1 .or. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     call MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     if (rank == 0) write(*,*) max_diff_global
     
     
     if (1 == 1 .and. h == 8) then
     if (grid == 0) then
        i = N1p
        k = N3p
        do j = S21, N21
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(65,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 2 .and. iB(3,1) == 1) write(66,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 3 .and. iB(3,1) == 1) write(67,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 4 .and. iB(3,1) == 1) write(68,*) x2p(j), Ap(i,j,k), rh(i,j,k)
        end do
     else if (grid == 3) then
        i = N1p
        k = N3p
        do j = S2p, N2p
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(65,*) x2p(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 2 .and. iB(3,1) == 1) write(66,*) x2p(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 3 .and. iB(3,1) == 1) write(67,*) x2p(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 4 .and. iB(3,1) == 1) write(68,*) x2p(j), Ap(i,j,k)
        end do
     else
        i = N12
        k = N32
        do j = S22, N22
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(65,*) x2v(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 2 .and. iB(3,1) == 1) write(66,*) x2v(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 3 .and. iB(3,1) == 1) write(67,*) x2v(j), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 4 .and. iB(3,1) == 1) write(68,*) x2v(j), Ap(i,j,k)
        end do
     end if
     end if
  end do
  
  !===========================================================================================================
  if (dimens == 3) then
  
  if (rank == 0) write(*,*)
  
  do h = 1, 8
     
     if (h == 1) then
        abl  = 0
        grid = 2
     else if (h == 2) then
        abl  = 0
        grid = 3
     else if (h == 3) then
        abl  = 1
        grid = 0
     else if (h == 4) then
        abl  = 1
        grid = 1
     else if (h == 5) then
        abl  = 1
        grid = 2
     else if (h == 6) then
        abl  = 1
        grid = 3
     else if (h == 7) then
        abl  = 2
        grid = 0
     else if (h == 8) then
        abl  = 2
        grid = 1
     end if
     
     
     xref = phase*0.25*L3/omega
     
     if (BC_3L_global == -2 .or. BC_3U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L3
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (grid == 0 .or. grid == 2) pp(i,j,k) =  COS(dd1*(x3p(k)-xref))
              
              if (grid == 0 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x3p(k)-xref))
              if (grid == 0 .and. abl == 2) rh(i,j,k) = -COS(dd1*(x3p(k)-xref))
              
              if (grid == 3 .and. abl == 0) rh(i,j,k) =  SIN(dd1*(x3p(k)-xref))
              if (grid == 3 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x3p(k)-xref))
           end do
        end do
     end do
     do k = S33B, N33B
        do j = S23B, N23B
           do i = S13B, N13B
              if (grid == 1 .or. grid == 3) pp(i,j,k) =  SIN(dd1*(x3w(k)-xref))
              
              if (grid == 1 .and. abl == 1) rh(i,j,k) =  COS(dd1*(x3w(k)-xref))
              if (grid == 1 .and. abl == 2) rh(i,j,k) = -SIN(dd1*(x3w(k)-xref))
              
              if (grid == 2 .and. abl == 0) rh(i,j,k) =  COS(dd1*(x3w(k)-xref))
              if (grid == 2 .and. abl == 1) rh(i,j,k) = -SIN(dd1*(x3w(k)-xref))
           end do
        end do
     end do
     
     pp = pp * dd1**(-abl)
     
     if (grid == 0 .or. grid == 2) call exchange(3,0,pp)
     if (grid == 1 .or. grid == 3) call exchange(3,3,pp)
     
     if (abl == 0 .and. grid == 2) call apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,pp,Ap)
     if (abl == 0 .and. grid == 3) call apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,pp,Ap)
     
     if (abl == 1 .and. grid == 0) call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ap)
     if (abl == 1 .and. grid == 1) call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ap)
     if (abl == 1 .and. grid == 2) call apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,pp,Ap)
     if (abl == 1 .and. grid == 3) call apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,pp,Ap)
     
     if (abl == 2 .and. grid == 0) call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ap)
     if (abl == 2 .and. grid == 0) call apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,pp,rr)
     
     if (abl == 2 .and. grid == 1) call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ap)
     if (abl == 2 .and. grid == 1) call apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,pp,rr)
     
     
     max_diff = 0.
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              if (abl == 1 .and. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx3DM(k)
              
              if (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     do k = S31, N31
        do j = S21, N21
           do i = S11, N11
              if (abl == 1 .and. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx3pM(k)
              if (abl == 2 .and. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx3pM(k)**2 + Ap(i,j,k)*ddx3pM(k)
              
              if (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     do k = S33, N33
        do j = S23, N23
           do i = S13, N13
              if (abl == 1 .and. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx3wM(k)
              if (abl == 1 .and. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx3GM(k)
              if (abl == 2 .and. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx3wM(k)**2 + Ap(i,j,k)*ddx3wM(k)
              
              if (grid == 1 .or. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           end do
        end do
     end do
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     call MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     if (rank == 0) write(*,*) max_diff_global
     
     
     if (1 == 1 .and. h == 3) then
     if (grid == 0) then
        i = N1p
        j = N2p
        do k = S31, N31
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(75,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 2) write(76,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 3) write(77,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 4) write(78,*) x3p(k), Ap(i,j,k)
        end do
     else if (grid == 3) then
        i = N1p
        j = N2p
        do k = S3p, N3p
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(75,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 2) write(76,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 3) write(77,*) x3p(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 4) write(78,*) x3p(k), Ap(i,j,k)
        end do
     else
        i = N13
        j = N23
        do k = S33, N33
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 1) write(75,*) x3w(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 2) write(76,*) x3w(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 3) write(77,*) x3w(k), Ap(i,j,k)
           if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. iB(3,1) == 4) write(78,*) x3w(k), Ap(i,j,k)
        end do
     end if
     end if
  end do
  
  end if
  !===========================================================================================================
  call MPI_FINALIZE(merror)
  stop
  
  
  end subroutine test_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  !> \brief initializes all FD coeff array with the help of diff_coeffs
  subroutine FD_coeffs()
  
  implicit none
  
  integer               ::  m
  
  
  !===========================================================================================================
  call diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1D)
  call diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1U)
  call diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1D)
  call diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1U)
  
  call diff_coeffs(1,0, 0,.false.    ,dim_ncb1c,ncb1f,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cFp1 )
  call diff_coeffs(1,0, 0,.false.    ,dim_ncb1c,ncb1f,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cFu1 )
  
  call diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp1  )
  call diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu1  )
  
  call diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp11 )
  call diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu11 )
  
  call diff_coeffs(1,1, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cGp1 )
  call diff_coeffs(1,1, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cDu1 )
  
  call diff_coeffs(1,0, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cIpu )
  call diff_coeffs(1,0, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cIup )
  !-----------------------------------------------------------------------------------------------------------
  call diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2D)
  call diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2U)
  call diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2D)
  call diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2U)
  
  call diff_coeffs(2,0, 0,.false.    ,dim_ncb2c,ncb2f,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cFp2 )
  call diff_coeffs(2,0, 0,.false.    ,dim_ncb2c,ncb2f,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cFv2 )
  
  call diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp2  )
  call diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv2  )
  
  call diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp22 )
  call diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv22 )
  
  call diff_coeffs(2,1, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cGp2 )
  call diff_coeffs(2,1, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cDv2 )
  
  call diff_coeffs(2,0, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cIpv )
  call diff_coeffs(2,0, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cIvp )
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  call diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3D)
  call diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3U)
  call diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3D)
  call diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3U)
  
  call diff_coeffs(3,0, 0,.false.    ,dim_ncb3c,ncb3f,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cFp3 )
  call diff_coeffs(3,0, 0,.false.    ,dim_ncb3c,ncb3f,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cFw3 )
  
  call diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp3  )
  call diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw3  )
  
  call diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp33 )
  call diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw33 )
  
  call diff_coeffs(3,1, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cGp3 )
  call diff_coeffs(3,1, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cDw3 )
  
  call diff_coeffs(3,0, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cIpw )
  call diff_coeffs(3,0, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cIwp )
  end if
  !===========================================================================================================
  
  
  do m = 1, n_conc
     !========================================================================================================
     call diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,n1L,n1U,cNc1D(n1L,0,m))
     call diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,n1L,n1U,cNc1U(n1L,0,m))
     call diff_coeffs(1,0, 0,.false.    ,dim_ncb1c,ncb1f,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cFc1 (b1L,0,m))
     call diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cc1  (b1L,0,m))
     call diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cc11 (b1L,0,m))
     call diff_coeffs(1,0, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BCc_1L(m),BCc_1U(m),2,N1,g1L,g1U,cIcu (g1L,0,m))
     !--------------------------------------------------------------------------------------------------------
     call diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,n2L,n2U,cNc2D(n2L,0,m))
     call diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,n2L,n2U,cNc2U(n2L,0,m))
     call diff_coeffs(2,0, 0,.false.    ,dim_ncb2c,ncb2f,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cFc2 (b2L,0,m))
     call diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cc2  (b2L,0,m))
     call diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cc22 (b2L,0,m))
     call diff_coeffs(2,0, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BCc_2L(m),BCc_2U(m),2,N2,g2L,g2U,cIcv (g2L,0,m))
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     call diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,n3L,n3U,cNc3D(n3L,0,m))
     call diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,n3L,n3U,cNc3U(n3L,0,m))
     call diff_coeffs(3,0, 0,.false.    ,dim_ncb3c,ncb3f,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cFc3 (b3L,0,m))
     call diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cc3  (b3L,0,m))
     call diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cc33 (b3L,0,m))
     call diff_coeffs(3,0, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BCc_3L(m),BCc_3U(m),2,N3,g3L,g3U,cIcw (g3L,0,m))
     end if
     !========================================================================================================
  end do
  
  !===========================================================================================================
                   call diff_coeffs(1,-1,0,mapping_yes,dim_ncb1c,ncb1c,x1p,x1p,BC_1L,BC_1U,0,N1,b1L,b1U,cInt1)
                   call diff_coeffs(2,-1,0,mapping_yes,dim_ncb2c,ncb2c,x2p,x2p,BC_2L,BC_2U,0,N2,b2L,b2U,cInt2)
  if (dimens == 3) call diff_coeffs(3,-1,0,mapping_yes,dim_ncb3c,ncb3c,x3p,x3p,BC_3L,BC_3U,0,N3,b3L,b3U,cInt3)
  !===========================================================================================================
  
  
  end subroutine FD_coeffs
  
  
  
  
  
  
  
  
  
  
  
  subroutine test_coeffs()
  
  implicit none
  
  
  ! ACHTUNG!!! Upwind-Differenzen noch nicht getestet!!!
  !===========================================================================================================
  if (iB(2,1) == 1 .and. iB(3,1) == 1) then
     call test_diff(1,0,x1p              ,x1p              ,N1,b1L,b1U,cFp1 ,BC_1L,0,'cFp1' )
     call test_diff(1,0,x1u              ,x1u              ,N1,b1L,b1U,cFu1 ,BC_1L,1,'cFu1' )
     call test_diff(1,1,x1u              ,x1u              ,N1,b1L,b1U,cu1  ,BC_1L,1,'cu1'  )
     call test_diff(1,1,x1p              ,x1p              ,N1,b1L,b1U,cp1  ,BC_1L,0,'cp1'  )
     call test_diff(1,2,x1u              ,x1u              ,N1,b1L,b1U,cu11 ,BC_1L,1,'cu11' )
     call test_diff(1,2,x1p              ,x1p              ,N1,b1L,b1U,cp11 ,BC_1L,0,'cp11' )
     
     !CALL test_diff(1,1,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cGp1 ,BC_1L,2,'cGp1' ) ! TEST!!! Funktioniert nicht mehr ...
     !CALL test_diff(1,1,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cDu1 ,BC_1L,3,'cDu1' )
     !CALL test_diff(1,0,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cIpu ,BC_1L,2,'cIpu' )
     !CALL test_diff(1,0,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cIup ,BC_1L,3,'cIup' )
     
     call test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1D,BC_1L,1,'cNu1D')
     call test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1U,BC_1L,1,'cNu1U')
     call test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1D,BC_1L,0,'cNp1D')
     call test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1U,BC_1L,0,'cNp1U')
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (iB(1,1) == 1 .and. iB(3,1) == 1) then
     call test_diff(2,0,x2p              ,x2p              ,N2,b2L,b2U,cFp2 ,BC_2L,0,'cFp2' )
     call test_diff(2,0,x2v              ,x2v              ,N2,b2L,b2U,cFv2 ,BC_2L,1,'cFv2' )
     call test_diff(2,1,x2p              ,x2p              ,N2,b2L,b2U,cp2  ,BC_2L,0,'cp2'  )
     call test_diff(2,1,x2v              ,x2v              ,N2,b2L,b2U,cv2  ,BC_2L,1,'cv2'  )
     call test_diff(2,2,x2p              ,x2p              ,N2,b2L,b2U,cp22 ,BC_2L,0,'cp22' )
     call test_diff(2,2,x2v              ,x2v              ,N2,b2L,b2U,cv22 ,BC_2L,1,'cv22' )
     
     !CALL test_diff(2,1,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cGp2 ,BC_2L,2,'cGp2' )
     !CALL test_diff(2,1,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cDv2 ,BC_2L,3,'cDv2' )
     !CALL test_diff(2,0,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cIpv ,BC_2L,2,'cIpv' )
     !CALL test_diff(2,0,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cIvp ,BC_2L,3,'cIvp' )
     
     call test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2D,BC_2L,1,'cNv2D')
     call test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2U,BC_2L,1,'cNv2U')
     call test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2D,BC_2L,0,'cNp2D')
     call test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2U,BC_2L,0,'cNp2U')
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (iB(1,1) == 1 .and. iB(2,1) == 1 .and. dimens == 3) then
     call test_diff(3,0,x3p              ,x3p              ,N3,b3L,b3U,cFp3 ,BC_3L,0,'cFp3' )
     call test_diff(3,0,x3w              ,x3w              ,N3,b3L,b3U,cFw3 ,BC_3L,1,'cFw3' )
     call test_diff(3,1,x3p              ,x3p              ,N3,b3L,b3U,cp3  ,BC_3L,0,'cp3'  )
     call test_diff(3,1,x3w              ,x3w              ,N3,b3L,b3U,cw3  ,BC_3L,1,'cw3'  )
     call test_diff(3,2,x3p              ,x3p              ,N3,b3L,b3U,cp33 ,BC_3L,0,'cp33' )
     call test_diff(3,2,x3w              ,x3w              ,N3,b3L,b3U,cw33 ,BC_3L,1,'cw33' )
     
     !CALL test_diff(3,1,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cGp3 ,BC_3L,2,'cGp3' )
     !CALL test_diff(3,1,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cDw3 ,BC_3L,3,'cDw3' )
     !CALL test_diff(3,0,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cIpw ,BC_3L,2,'cIpw' )
     !CALL test_diff(3,0,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cIwp ,BC_3L,3,'cIwp' )
     
     call test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3D,BC_3L,1,'cNw3D')
     call test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3U,BC_3L,1,'cNw3U')
     call test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3D,BC_3L,0,'cNp3D')
     call test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3U,BC_3L,0,'cNp3U')
  end if
  !===========================================================================================================
  
  
  end subroutine test_coeffs
  
  
  
  
  
  
  
  
  
  
  !> \brief calculates differential coefficients.
  !!
  !! \param[in] dir direction
  !! \param[in] abl degree of derivative
  !! \param[in] upwind (0: central, 1: up, -1: low)
  !! \param[in] mapping_yes ??
  !! \param[in] dim_ncb dimension of ??
  !! \param[in] n_coeff_bound ??
  !! \param[in] xC coordinate in
  !! \param[in] xE coordinate out
  !! \param[in] BCL
  !! \param[in] BCU
  !! \param[in] grid_type (0: p, 1: u, 2: v, 3:w )
  !! \param[in] bl  lower bound
  !! \param[in] bu  upper bound
  !! \param[out] cc coefficients
  subroutine diff_coeffs(dir,abl,upwind,mapping_yes,dim_ncb,n_coeff_bound,xC,xE,BCL,BCU,grid_type,Nmax,bL,bU,cc)
  
  implicit none
  
  integer, intent(in)   ::  dir
  integer, intent(in)   ::  abl
  integer, intent(in)   ::  upwind
  logical, intent(in)   ::  mapping_yes
  integer, intent(in)   ::  dim_ncb
  integer, intent(in)   ::  n_coeff_bound(1:dim_ncb)
  integer, intent(in)   ::  Nmax
  integer, intent(in)   ::  bL
  integer, intent(in)   ::  bU
  real   , intent(in)   ::  xC(bL:(Nmax+bU))
  real   , intent(in)   ::  xE(bL:(Nmax+bU))
  integer, intent(in)   ::  BCL
  integer, intent(in)   ::  BCU
  integer, intent(in)   ::  grid_type
  
  real   , intent(out)  ::  cc(bL:bU,0:Nmax)
  
  integer               ::  n_coeff
  integer               ::  dim_n_coeff_bound
  integer               ::  i, ii, iC, iStart, SShift
  integer               ::  k, kk
  
  integer               ::  left, right
  
  real                  ::  dxi(1:2)
  
  real                  ::  dxL, dxU ! Für Integrationskoeffizienten
  
  real   , allocatable  ::  cc_xi(:,:)
  real   , allocatable  ::  deltaX(:)
  
  logical               ::  filter_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Upwinding wird auf dem Rand unterdrücken, da Differenzenstencils dort ohnehin schief sind.!
  !              - Generell koennte/sollte die Initialisierung der Index-Grenzen auch schon vor dieser       !
  !                Routine ausgefuehrt werden, um hier die Uebersichtlichkeit zu verbessern.                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  if (mapping_yes .and. abl > 2) then
     if (rank == 0) write(*,*) 'ERROR! Can`t handle derivatives > 2 combined with mapping ...'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  if (dir == 1) SShift = iShift ! TEST!!! Besser als Argument übergeben ...
  if (dir == 2) SShift = jShift
  if (dir == 3) SShift = kShift
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  if (grid_type == 0 .or. grid_type == 3) then
     iStart = 1
  else
     iStart = 0
  end if
  !===========================================================================================================
  
  filter_yes = .false.
  if (abl == 0 .and. (grid_type == 0 .or. grid_type == 1)) filter_yes = .true.
  
  
  dim_n_coeff_bound = SIZE(n_coeff_bound)
  
  cc = 0.
  
  do i = iStart, Nmax
     
     !========================================================================================================
     !=== Stencil-Breite auslesen ============================================================================
     !========================================================================================================
     if      (BCL > 0 .and. i <= (iStart - 1 + dim_n_coeff_bound)) then
        n_coeff = n_coeff_bound(i + 1 - iStart)
        !IF (upwind /= 0 .AND. i == iStart) n_coeff = 0
     else if (BCU > 0 .and. i >= (Nmax   + 1 - dim_n_coeff_bound)) then
        n_coeff = n_coeff_bound(Nmax + 1 - i)
        !IF (upwind /= 0 .AND. i == Nmax  ) n_coeff = 0
     else
        n_coeff = n_coeff_bound(dim_n_coeff_bound)
     end if
     
     !========================================================================================================
     
     if (n_coeff > 0 .and. n_coeff > abl) then
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten RECHTS vom Entwicklungspunkt bestimmen =================================
        !=====================================================================================================
        if (grid_type == 0 .or. grid_type == 1) then
           !==================================================================================================
           !=== "normales" Gitter, keine Zwischengitterpunkte ================================================
           !==================================================================================================
           if      (BCL > 0 .and. i <= (iStart - 1 + n_coeff/2)) then
              right = n_coeff - i + iStart - 1
           else if (BCU > 0 .and. i >= (Nmax   + 1 - n_coeff/2)) then
              right = Nmax - i
           else
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= ungerade Anzahl Koeffizienten)!
              if (MOD(n_coeff,2) == 0) then
                 if (rank == 0) then
                    write(*,'(a   )') 'ERROR! Choose odd number of coefficients!'
                    write(*,'(a,i4)') '    direction =', dir
                    write(*,'(a,i4)') '    grid_type =', grid_type
                    write(*,'(a,i4)') '            i =', i+SShift
                 end if
              end if
              right = (n_coeff-1)/2
           end if
        else if (grid_type == 2) then
           !==================================================================================================
           !=== Druck ==> Impuls =============================================================================
           !==================================================================================================
           if      (BCL > 0 .and. i <= n_coeff/2) then
              right = n_coeff - i - iStart
           else if (BCU > 0 .and. i >= (Nmax + 0 - n_coeff/2)) then
              right = Nmax - i
           else
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              if (MOD(n_coeff,2) /= 0) then
                 if (rank == 0) then
                    write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    write(*,'(a,i4)') '    direction =', dir
                    write(*,'(a,i4)') '    grid_type =', grid_type
                    write(*,'(a,i4)') '            i =', i+SShift
                 end if
                 call MPI_FINALIZE(merror)
                 stop
              end if
              right = n_coeff/2
           end if
        else if (grid_type == 3) then
           !==================================================================================================
           !=== Geschwindigkeit ==> Konti ====================================================================
           !==================================================================================================
           if      (BCL > 0 .and. i <= n_coeff/2) then
              right = n_coeff - i - iStart
           else if (BCU > 0 .and. i >= (Nmax + 1 - n_coeff/2)) then
              right = Nmax - i
           else
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              if (MOD(n_coeff,2) /= 0) then
                 if (rank == 0) then
                    write(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    write(*,'(a,i4)') '    direction =', dir
                    write(*,'(a,i4)') '    grid_type =', grid_type
                    write(*,'(a,i4)') '            i =', i+SShift
                 end if
                 call MPI_FINALIZE(merror)
                 stop
              end if
              right = n_coeff/2 - 1
           end if
        end if
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten LINKS vom Entwicklungspunkt bestimmen ==================================
        !=====================================================================================================
        ! (0-ter bzw. 1-ter Koeffizient links vom Entwicklungspunkt
        !          == "zentraler" Koeffizient im Speicher)
        left = right - n_coeff + 1
        !=====================================================================================================
        
        
        !=====================================================================================================
        if (upwind == -1) then
           if (.not. (BCL > 0 .and. i < (iStart + n_coeff/2))) then
              n_coeff = n_coeff - 1
              left    = left    + 1
           end if
        end if
        if (upwind ==  1) then
           if (.not. (BCU > 0 .and. i > (Nmax   - n_coeff/2))) then
              n_coeff = n_coeff - 1
              right   = right   - 1
           end if
        end if
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Stencilanordnung testen =========================================================================
        !=====================================================================================================
        if (right > bU .or. left < bL) then
           if (rank == 0) then
              !WRITE(*,'(a   )') 'WARNING! The FD-Stencil does probably not fit into provided array!'
              write(*,'(a   )') 'ERROR! Stencil doesn`t fit into provided array!'
              write(*,'(a,i4)') '    direction =', dir
              write(*,'(a,i4)') '    grid_type =', grid_type
              write(*,'(a,i4)') '            i =', i+SShift
           end if
           call MPI_FINALIZE(merror)
           stop
        end if
        !=====================================================================================================
        
        
        allocate(deltaX(1:n_coeff))
        
        
        !=====================================================================================================
        !=== räumliche Abstände zum Entwicklungspunkt ========================================================
        !=====================================================================================================
        do ii = 1, n_coeff
           ! Koeffizienten-Punkte ("iC"):
           iC = i + left + (ii-1)
           
           if (mapping_yes) then
              if      (grid_type == 2) then
                 deltaX(ii) = REAL(iC-i)-0.5
              else if (grid_type == 3) then
                 deltaX(ii) = REAL(iC-i)+0.5
              else
                 deltaX(ii) = REAL(iC-i)
              end if
           else
              deltaX(ii) = xC(iC) - xE(i)
           end if
        end do
        !=====================================================================================================
        
        
        
        !=====================================================================================================
        !=== Integrations-Intervall ==========================================================================
        !=====================================================================================================
        ! Anmerkungen: - Kein Mapping aus Genauigkeitsgründen vorgesehen.
        !              - Nur für Druckgitter vorgesehen.
        if (abl == -1) then
           
           if (BCL > 0 .and. i == iStart) then
              dxL = 0.
           else
              dxL = (xC(i-1) - xE(i)) / 2. ! der Einfachheit halber anstelle der Geschwindigkeitspunkte ...
           end if
           if (BCU > 0 .and. i == Nmax  ) then
              dxU = 0.
           else
              dxU = (xC(i+1) - xE(i)) / 2.
           end if
           
        end if
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Bestimmung der Koeffizienten ====================================================================
        !=====================================================================================================
        ! Anmerkung: Filter und Integratoren sollen nicht gemappt werden (siehe z.B. Stolz).
        if (mapping_yes .and. abl /= 0 .and. abl /= -1) then
           
           !--- derivatives (mapped) ---
           
           allocate(cc_xi(left:right,1:abl))
           
           cc_xi = 0.
           dxi   = 0.
           
           do k = 1, abl
              if (n_coeff <= 7) then ! TEST!!!
                 !kk = k
                 !INCLUDE 'FD_coeffs_expl_map.f90'
                 call diff_coeffs_exact(k,n_coeff,deltaX(1),cc_xi(left:right,k)) ! TEST!!!
              else
                 call FD_coeffs_solver(k,filter_yes,n_coeff,deltaX,cc_xi(left:right,k))
              end if
              
              do ii = left, right
                 dxi(k) = dxi(k) + cc_xi(ii,k)*xC(i+ii)
              end do
           end do
           
           !---------------------------------------------------!
           ! Mapping:                                          !
           ! d/dx     = (xi') * d/dxi                          !
           ! d^2/dx^2 = (xi") * d/dxi + (xi')**2 * d^2/dxi^2   !
           !                                                   !
           ! Mapping-Faktoren:                                 !
           ! xi' =   1  /  x'                                  !
           ! xi" = - x" / (x')**3                              !
           !---------------------------------------------------!
           
           if (abl == 1) cc(left:right,i) = cc_xi(left:right,1)/dxi(1)
           if (abl == 2) cc(left:right,i) = cc_xi(left:right,2)/dxi(1)**2 - cc_xi(left:right,1)*dxi(2)/dxi(1)**3
           
           deallocate(cc_xi)
           
        else if (mapping_yes .and. abl == 0 .and. .not. filter_yes) then ! TEST!!!
           
           !--- interpolation (mapped) ---
           
           !k  = 1
           !kk = 0
           !ALLOCATE(cc_xi(left:right,k:k))
           !INCLUDE 'FD_coeffs_expl_map.f90'
           !cc(left:right,i) = cc_xi(left:right,k)
           !DEALLOCATE(cc_xi)
           
           if (n_coeff <= 7) then ! TEST!!!
              call diff_coeffs_exact(abl,n_coeff,deltaX(1),cc(left:right,i))
           else
              call FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           end if
           
        else
           
           if (abl == -1) then
              !--- integration coefficients ---
              call FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc(left:right,i))
           else
              !--- interpolation and derivatives ---
              call FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           end if
           
        end if
        !=====================================================================================================
        
        
        deallocate(deltaX)
        
        
        !=====================================================================================================
        !=== Symmetrie =======================================================================================
        !=====================================================================================================
        if (BCL == -2 .and. abl == -1 .and. i == iStart) then
           cc(0,i) = 0.5*cc(0,i)
           do ii = bL, -1
              cc(ii,i) = 0.
           end do
        end if
        if (BCU == -2 .and. abl == -1 .and. i == Nmax  ) then
           cc(0,i) = 0.5*cc(0,i)
           do ii = 1, bU
              cc(ii,i) = 0.
           end do
        end if
        !-----------------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        ! TEST!!! abl == -1 ???
        if (BCL == -2 .and. (i+bL) < 1) then
           if (grid_type == 0 .or. (grid_type == 2 .and. i >= 1)) then
              do ii = 1, 1-(i+bL)
                 cc(1-i+ii,i) = cc(1-i+ii,i) + cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              end do
           end if
           if ((grid_type == 1 .and. i >= 1) .or. grid_type == 3) then
              do ii = 1, 1-(i+bL)
                 cc(0-i+ii,i) = cc(0-i+ii,i) - cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              end do
           end if
        end if
        !-----------------------------------------------------------------------------------------------------
        if (BCU == -2 .and. (i+bU) > (Nmax+0)) then
           if (grid_type == 0 .or. (grid_type == 2 .and. i <= Nmax-1)) then
              do ii = 1, i+bU-(Nmax-0)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) + cc(Nmax-i  +ii,i)
                 cc(Nmax-i  +ii,i) = 0.
              end do
           end if
        end if
        if (BCU == -2 .and. (i+bU) > (Nmax-1)) then
           if ((grid_type == 1 .and. i <= Nmax-1) .or. grid_type == 3) then
              do ii = 1, i+bU-(Nmax-1)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) - cc(Nmax-i-1+ii,i)
                 cc(Nmax-i-1+ii,i) = 0.
              end do
           end if
        end if
        !=====================================================================================================
        
     end if
     
  end do
  
  
  end subroutine diff_coeffs
  
  
  
  
  
  
  
  
  
  
  
  subroutine diff_coeffs_exact(abl,n_coeff,deltaX,cc) ! TEST!!!
  
  implicit none
  
  integer, intent(in   ) ::  abl
  integer, intent(in   ) ::  n_coeff
  real   , intent(in   ) ::  deltaX
  real   , intent(out  ) ::  cc(1:n_coeff)
  
  ! TEST!!!
  ! - korrekt ausgerichtet?
  ! - Vorzeichen ok?
  
  
  if (n_coeff == 2) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,-1./)/2.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 1., 1./)/2.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/-1., 3./)/2.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     
     if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
  end if
  if (n_coeff == 3) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/15.,-10., 3./)/8.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,  6.,-1./)/8.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/-1.,  6., 3./)/8.
     if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/ 3.,-10.,15./)/8.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-2.,  3.,-1./)/1.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1.,  1., 0./)/1.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/ 0., -1., 1./)/1.
     if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/ 1., -3., 2./)/1.
     
     if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-3.,  4.,-1./)/2.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/-1.,  0., 1./)/2.
     if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/ 1., -4., 3./)/2.
     
     if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
  end if
  if (n_coeff == 4) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/ 35.,-35.,  21., -5./)/16.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/  5., 15.,  -5.,  1./)/16.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/ -1.,  9.,   9., -1./)/16.
     if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/  1., -5.,  15.,  5./)/16.
     if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/ -5., 21., -35., 35./)/16.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-71.,141., -93., 23./)/24.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-23., 21.,   3., -1./)/24.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/  1.,-27.,  27., -1./)/24.
     if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/  1., -3., -21., 23./)/24.
     if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/-23., 93.,-141., 71./)/24.
     
     if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-11., 18.,  -9.,  2./)/6.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/ -2., -3.,   6., -1./)/6.
     if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/  1., -6.,   3.,  2./)/6.
     if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/ -2.,  9., -18., 11./)/6.
     
     if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/  2., -5.,   4., -1./)/1.
     if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/  1., -2.,   1.,  0./)/1.
     if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/  0.,  1.,  -2.,  1./)/1.
     if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,  4.,  -5.,  2./)/1.
  end if
  if (n_coeff == 5) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/315.,-420., 378.,-180., 35./)/128.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/ 35., 140., -70.,  28., -5./)/128.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/ -5.,  60.,  90., -20.,  3./)/128.
     if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/  3., -20.,  90.,  60., -5./)/128.
     if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/ -5.,  28., -70., 140., 35./)/128.
     if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/ 35.,-180., 378.,-420.,315./)/128.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-93., 229.,-225., 111.,-22./)/24.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-22.,  17.,   9.,  -5.,  1./)/24.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/  1., -27.,  27.,  -1.,  0./)/24.
     if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/  0.,   1., -27.,  27., -1./)/24.
     if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/ -1.,   5.,  -9., -17., 22./)/24.
     if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/ 22.,-111., 225.,-229., 93./)/24.
     
     if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/-25.,  48., -36.,  16., -3./)/12.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/ -3., -10.,  18.,  -6.,  1./)/12.
     if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/  1.,  -8.,   0.,   8., -1./)/12.
     if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/ -1.,   6., -18.,  10.,  3./)/12.
     if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/  3., -16.,  36., -48., 25./)/12.
     
     if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/ 35.,-104., 114., -56., 11./)/12.
     if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/ 11., -20.,   6.,   4., -1./)/12.
     if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,  16., -30.,  16., -1./)/12.
     if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/ -1.,   4.,   6., -20., 11./)/12.
     if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/ 11., -56., 114.,-104., 35./)/12.
  end if
  if (n_coeff == 6) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/  693.,-1155.,  1386., -990.,   385., -63./)/256.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/   63.,  315.,  -210.,  126.,   -45.,   7./)/256.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/   -7.,  105.,   210.,  -70.,    21.,  -3./)/256.
     if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/    3.,  -25.,   150.,  150.,   -25.,   3./)/256.
     if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/   -3.,   21.,   -70.,  210.,   105.,  -7./)/256.
     if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/    7.,  -45.,   126., -210.,   315.,  63./)/256.
     if (deltaX == -5.5 .and. abl == 0) cc(1:n_coeff) = (/  -63.,  385.,  -990., 1386., -1155., 693./)/256.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-9129.,26765.,-34890.,25770.,-10205.,1689./)/1920.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/-1689., 1005.,  1430.,-1110.,   435., -71./)/1920.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/   71.,-2115.,  2070.,   10.,   -45.,   9./)/1920.
     if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/   -9.,  125., -2250., 2250.,  -125.,   9./)/1920.
     if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/   -9.,   45.,   -10.,-2070.,  2115., -71./)/1920.
     if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/   71., -435.,  1110.,-1430., -1005.,1689./)/1920.
     if (deltaX == -5.5 .and. abl == 1) cc(1:n_coeff) = (/-1689.,10205.,-25770.,34890.,-26765.,9129./)/1920.
     
     if (deltaX ==  0.0 .and. abl == 1) cc(1:n_coeff) = (/ -137.,  300.,  -300.,  200.,   -75.,  12./)/60.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/  -12.,  -65.,   120.,  -60.,    20.,  -3./)/60.
     if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/    3.,  -30.,   -20.,   60.,   -15.,   2./)/60.
     if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/   -2.,   15.,   -60.,   20.,    30.,  -3./)/60.
     if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/    3.,  -20.,    60., -120.,    65.,  12./)/60.
     if (deltaX == -5.0 .and. abl == 1) cc(1:n_coeff) = (/  -12.,   75.,  -200.,  300.,  -300., 137./)/60.
     
     if (deltaX ==  0.0 .and. abl == 2) cc(1:n_coeff) = (/   45., -154.,   214., -156.,    61., -10./)/12.
     if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/   10.,  -15.,    -4.,   14.,    -6.,   1./)/12.
     if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/   -1.,   16.,   -30.,   16.,    -1.,   0./)/12.
     if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/    0.,   -1.,    16.,  -30.,    16.,  -1./)/12.
     if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/    1.,   -6.,    14.,   -4.,   -15.,  10./)/12.
     if (deltaX == -5.0 .and. abl == 2) cc(1:n_coeff) = (/  -10.,   61.,  -156.,  214.,  -154.,  45./)/12.
  end if
  if (n_coeff == 7) then
     if (deltaX ==  0.5 .and. abl == 0) cc(1:n_coeff) = (/  3003., -6006.,  9009.,-8580.,   5005., -1638.,  231./)/1024.
     if (deltaX == -0.5 .and. abl == 0) cc(1:n_coeff) = (/   231.,  1386., -1155.,  924.,   -495.,   154.,  -21./)/1024.
     if (deltaX == -1.5 .and. abl == 0) cc(1:n_coeff) = (/   -21.,   378.,   945., -420.,    189.,   -54.,    7./)/1024.
     if (deltaX == -2.5 .and. abl == 0) cc(1:n_coeff) = (/     7.,   -70.,   525.,  700.,   -175.,    42.,   -5./)/1024.
     if (deltaX == -3.5 .and. abl == 0) cc(1:n_coeff) = (/    -5.,    42.,  -175.,  700.,    525.,   -70.,    7./)/1024.
     if (deltaX == -4.5 .and. abl == 0) cc(1:n_coeff) = (/     7.,   -54.,   189., -420.,    945.,   378.,  -21./)/1024.
     if (deltaX == -5.5 .and. abl == 0) cc(1:n_coeff) = (/   -21.,   154.,  -495.,  924.,  -1155.,  1386.,  231./)/1024.
     if (deltaX == -6.5 .and. abl == 0) cc(1:n_coeff) = (/   231., -1638.,  5005.,-8580.,   9009., -6006., 3003./)/1024.
     
     if (deltaX ==  0.5 .and. abl == 1) cc(1:n_coeff) = (/-10756., 36527.,-59295., 58310.,-34610., 11451.,-1627./)/1920.
     if (deltaX == -0.5 .and. abl == 1) cc(1:n_coeff) = (/ -1627.,   633.,  2360., -2350.,  1365.,  -443.,   62./)/1920.
     if (deltaX == -1.5 .and. abl == 1) cc(1:n_coeff) = (/    62., -2061.,  1935.,   190.,  -180.,    63.,   -9./)/1920.
     if (deltaX == -2.5 .and. abl == 1) cc(1:n_coeff) = (/    -9.,   125., -2250.,  2250.,  -125.,     9.,    0./)/1920.
     if (deltaX == -3.5 .and. abl == 1) cc(1:n_coeff) = (/     0.,    -9.,   125., -2250.,  2250.,  -125.,    9./)/1920.
     if (deltaX == -4.5 .and. abl == 1) cc(1:n_coeff) = (/     9.,   -63.,   180.,  -190., -1935.,  2061.,  -62./)/1920.
     if (deltaX == -5.5 .and. abl == 1) cc(1:n_coeff) = (/   -62.,   443., -1365.,  2350., -2360.,  -633., 1627./)/1920.
     if (deltaX == -6.5 .and. abl == 1) cc(1:n_coeff) = (/  1627.,-11451., 34610.,-58310., 59295.,-36527.,10756./)/1920.
     
     if (deltaX == -0.0 .and. abl == 1) cc(1:n_coeff) = (/  -147.,   360.,  -450.,   400.,  -225.,    72.,  -10./)/60.
     if (deltaX == -1.0 .and. abl == 1) cc(1:n_coeff) = (/   -10.,   -77.,   150.,  -100.,    50.,   -15.,    2./)/60.
     if (deltaX == -2.0 .and. abl == 1) cc(1:n_coeff) = (/     2.,   -24.,   -35.,    80.,   -30.,     8.,   -1./)/60.
     if (deltaX == -3.0 .and. abl == 1) cc(1:n_coeff) = (/    -1.,     9.,   -45.,     0.,    45.,    -9.,    1./)/60.
     if (deltaX == -4.0 .and. abl == 1) cc(1:n_coeff) = (/     1.,    -8.,    30.,   -80.,    35.,    24.,   -2./)/60.
     if (deltaX == -5.0 .and. abl == 1) cc(1:n_coeff) = (/    -2.,    15.,   -50.,   100.,  -150.,    77.,   10./)/60.
     if (deltaX == -6.0 .and. abl == 1) cc(1:n_coeff) = (/    10.,   -72.,   225.,  -400.,   450.,  -360.,  147./)/60.
     
     if (deltaX == -0.0 .and. abl == 2) cc(1:n_coeff) = (/   812., -3132.,  5265., -5080.,  2970.,  -972.,  137./)/180.
     if (deltaX == -1.0 .and. abl == 2) cc(1:n_coeff) = (/   137.,  -147.,  -225.,   470.,  -285.,    93.,  -13./)/180.
     if (deltaX == -2.0 .and. abl == 2) cc(1:n_coeff) = (/   -13.,   228.,  -420.,   200.,    15.,   -12.,    2./)/180.
     if (deltaX == -3.0 .and. abl == 2) cc(1:n_coeff) = (/     2.,   -27.,   270.,  -490.,   270.,   -27.,    2./)/180.
     if (deltaX == -4.0 .and. abl == 2) cc(1:n_coeff) = (/     2.,   -12.,    15.,   200.,  -420.,   228.,  -13./)/180.
     if (deltaX == -5.0 .and. abl == 2) cc(1:n_coeff) = (/   -13.,    93.,  -285.,   470.,  -225.,  -147.,  137./)/180.
     if (deltaX == -6.0 .and. abl == 2) cc(1:n_coeff) = (/   137.,  -972.,  2970., -5080.,  5265., -3132.,  812./)/180.
  end if
  
  
  end subroutine diff_coeffs_exact
  
  
  
  
  
  
  
  
  
  
  
  subroutine test_diff(dir,abl,xC,xE,Nmax,bL,bU,cc,BCL,grid_type,name)
  
  implicit none
  
  integer, intent(in)   ::  dir
  integer, intent(in)   ::  abl
  integer, intent(in)   ::  Nmax
  integer, intent(in)   ::  bL
  integer, intent(in)   ::  bU
  real   , intent(in)   ::  xC(bL:(Nmax+bU))
  real   , intent(in)   ::  xE(bL:(Nmax+bU))
  real   , intent(in)   ::  cc(bL:bU,0:Nmax)
  integer, intent(in)   ::  BCL
  integer, intent(in)   ::  grid_type
  
  character(*), intent(in) ::  name
  
  integer               ::  i, ii, iStartC, iStartE, SShift
  integer               ::  nw, nw_max, n_dx
  real                  ::  wave, wave_mod(1:2), dx
  
  character(len=1)      ::  part
  real                  ::  pi
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Hier sollen alle Differenzen-Stencils innerhalb eines Blockes getestet werden, analog zu  !
  !                ihrer Berechnung. Daher wird bewusst nicht mittels MPI in eine Datei geschrieben, sondern !
  !                in verschiedene Blöcke, um Fehler bei der MPI-Programmierung auszuschliessen.             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  if (dir == 1) then
     write(part,'(i1.1)') iB(1,1)
     SShift = iShift
  else if (dir == 2) then
     write(part,'(i1.1)') iB(2,1)
     SShift = jShift
  else if (dir == 3) then
     write(part,'(i1.1)') iB(3,1)
     SShift = kShift
  else
     if (rank == 0) write(*,*) 'ERROR! Wrong input at subroutine `test_diff`!'
     call MPI_FINALIZE(merror)
     stop
  end if
  
  
  open(10,file='test_'//name//'_transfer_block'//part//'.txt',status='UNKNOWN')
  open(11,file='test_'//name//'_coeffs_block'  //part//'.txt',status='UNKNOWN')
  
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  if (BCL <= 0 .or. grid_type == 0 .or. grid_type == 3) then
     iStartE = 1
  else
     iStartE = 0
  end if
  
  
  !===========================================================================================================
  !=== Startindex für Koeffizienten ==========================================================================
  !===========================================================================================================
  if (BCL > 0 .and. (grid_type == 1 .or. grid_type == 3)) then
     iStartC = 0
  else
     iStartC = 1
  end if
  
  
  !===========================================================================================================
  !=== Test der Koeffizienten ================================================================================
  !===========================================================================================================
  nw_max = 100
  
  do i = iStartE, Nmax
     
     do nw = 1, nw_max
        
        wave     = pi*REAL(nw)/REAL(nw_max)
        wave_mod = 0.
        
        !=== Referenz-Gitterweite bestimmen ==================================================================
        dx   = 0.
        n_dx = 0
        
        do ii = bL, bU-1
           !IF (i+ii .GE. iStartC .AND. i+ii .LT. Nmax .AND. ABS(xC(i+ii+1)-xC(i+ii)) .GT. dx) THEN
           if (i+ii >= iStartC .and. i+ii < Nmax) then
              dx   = dx + ABS(xC(i+ii+1)-xC(i+ii))
              n_dx = n_dx + 1
           end if
        end do
        
        dx = dx / REAL(n_dx)
        
        
        !=== Transferfunktion plotten ========================================================================
        if (abl == 0) then
           do ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))
              wave_mod(2) = wave_mod(2) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))
           end do
        else if (abl == 1) then
           do ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx
              wave_mod(2) = wave_mod(2) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx/wave**2
           end do
        else if (abl == 2) then
           do ii = bL, bU
              wave_mod(1) = wave_mod(1) - cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx**2
              wave_mod(2) = wave_mod(2) - cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx**2/wave
           end do
        end if
        
        write(10,'(i7,3E26.17)') i+SShift, wave, wave_mod(1:2)
        
     end do
     
     write(10,*)
     write(11,'(i7,100E26.17)') i+SShift, cc(:,i)
     
  end do
  
  close(10)
  close(11)
  
  
  end subroutine test_diff
  
  
  
  
  
  
  
  
  
  
  
  subroutine FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc)
  
  implicit none
  
  integer, intent(in)   ::  abl
  logical, intent(in)   ::  filter_yes
  integer, intent(in)   ::  n_coeff
  real   , intent(in)   ::  deltaX(1:n_coeff)
  real   , intent(out)  ::  cc    (1:n_coeff)
  
  integer               ::  i, j
  
  real                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  real                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  real                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  do i = 1, n_coeff
     do j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     end do
  end do
  
  
  !===========================================================================================================
  !=== Zusatzbedingungen =====================================================================================
  !===========================================================================================================
  if (filter_yes) then
     ! G(pi) = 0.
     do i = 1, n_coeff
        polyn_vals(i,n_coeff) = (-1.)**i
     end do
  end if
  
  
  !===========================================================================================================
  !=== Lösen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  call Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Funktionswert am Entwicklungspunkt (deltaX = 0) =======================================================
  !===========================================================================================================
  const = 1.
  do i = 1, abl-1
     const = const*REAL(1+i)
  end do
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  cc(1:n_coeff) = const*polyn_vals_inv(1+abl,1:n_coeff)
  
  
  end subroutine FD_coeffs_solver
  
  
  
  
  
  
  
  
  
  
  
  subroutine FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc)
  
  implicit none
  
  integer, intent(in)   ::  n_coeff
  real   , intent(in)   ::  deltaX(1:n_coeff)
  real   , intent(out)  ::  cc    (1:n_coeff)
  real   , intent(in)   ::  dxL, dxU
  
  integer               ::  i, j, k
  
  real                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  real                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  real                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  do i = 1, n_coeff
     do j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     end do
  end do
  
  
  !===========================================================================================================
  !=== L�sen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  call Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  !  (Matrix-Vektor-Multiplikation)
  cc = 0.
  do j = 1, n_coeff
     do i = 1, n_coeff
        cc(j) = cc(j) + (dxU**i - dxL**i)/REAL(i)*polyn_vals_inv(i,j)
     end do
  end do
  
  
  end subroutine FD_coeffs_solver_integral
  
  
  
  
  
  
  
  
  
  
  
  subroutine get_stencil()
  
  implicit none
  
  integer                ::  i, imax
  integer                ::  j, jmax
  integer                ::  k, kmax
  
  integer                ::  g
  
  real                   ::  cDu1R(-1:0,1:N1)
  real                   ::  cDv2R(-1:0,1:N2)
  real                   ::  cDw3R(-1:0,1:N3)
  
  real                   ::  cGp1R( 0:1,0:N1)
  real                   ::  cGp2R( 0:1,0:N2)
  real                   ::  cGp3R( 0:1,0:N3)
  
  
  cdg1 = 0.
  cdg2 = 0.
  cdg3 = 0.
  
  do g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     cDu1R = 0.
     cGp1R = 0.
     
     !---------------------------------------------------------------------------!
     ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
     !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
     !---------------------------------------------------------------------------!
     
     do i = 1, imax
        cDu1R(-1,i) = -1./(x1uR(i  ,g) - x1uR(i-1,g))
        cDu1R( 0,i) =  1./(x1uR(i  ,g) - x1uR(i-1,g))
     end do
     
     do i = 1, imax-1
        cGp1R( 0,i) = -1./(x1pR(i+1,g) - x1pR(i  ,g))
        cGp1R( 1,i) =  1./(x1pR(i+1,g) - x1pR(i  ,g))
     end do
     
     if (BC(1,1,g) > 0) then
        cGp1R( :,0   ) =  0.
        cDu1R( :,1   ) =  0.
        !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
        cDu1R( 0,1   ) =  1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     else
        cGp1R( 0,0   ) = -1./(x1pR(1,g) - x1pR(0,g))
        cGp1R( 1,0   ) =  1./(x1pR(1,g) - x1pR(0,g))
     end if
     
     if (BC(2,1,g) > 0) then
        cGp1R( :,imax) =  0.
        cDu1R( :,imax) =  0.
        !cDu1R(-1,imax) = -1./(x1u (N1      ) - x1u (N1-1  )) ! TEST!!! Das geht in die Hose ...
        cDu1R(-1,imax) = -1./(y1u (M1      ) - y1u (M1-1  ))
     else
        cGp1R( 0,imax) = -1./(x1pR(imax+1,g) - x1pR(imax,g))
        cGp1R( 1,imax) =  1./(x1pR(imax+1,g) - x1pR(imax,g))
     end if
     !--------------------------------------------------------------------------------------------------------
     do i = 1, imax
        cdg1(-1,i,g) = cDu1R(-1,i)*cGp1R(0,i-1)
        cdg1( 0,i,g) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
        cdg1( 1,i,g) =                            cDu1R(0,i)*cGp1R(1,i)
     end do
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     if (BC(1,1,g) > 0) cdg1( :,1   ,g) = 2.*cdg1(:,1   ,g) ! TEST!!! Ist das wirklich so optimal?
     if (BC(2,1,g) > 0) cdg1( :,imax,g) = 2.*cdg1(:,imax,g)
     
     if (BC(1,1,g) == -2) then
        cdg1( 1,1   ,g) = cdg1( 1,1   ,g) + cdg1(-1,1   ,g)
        cdg1(-1,1   ,g) = 0.
     end if
     if (BC(2,1,g) == -2) then
        cdg1(-1,imax,g) = cdg1(-1,imax,g) + cdg1( 1,imax,g)
        cdg1( 1,imax,g) = 0.
     end if
     !========================================================================================================
     jmax = NN(2,g)
     
     cDv2R = 0.
     cGp2R = 0.
     
     do j = 1, jmax
        cDv2R(-1,j) = -1./(x2vR(j  ,g) - x2vR(j-1,g))
        cDv2R( 0,j) =  1./(x2vR(j  ,g) - x2vR(j-1,g))
     end do
     
     do j = 1, jmax-1
        cGp2R( 0,j) = -1./(x2pR(j+1,g) - x2pR(j  ,g))
        cGp2R( 1,j) =  1./(x2pR(j+1,g) - x2pR(j  ,g))
     end do
     
     if (BC(1,2,g) > 0) then
        cGp2R( :,0   ) =  0.
        cDv2R( :,1   ) =  0.
        !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
        cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
     else
        cGp2R( 0,0   ) = -1./(x2pR(1,g) - x2pR(0  ,g))
        cGp2R( 1,0   ) =  1./(x2pR(1,g) - x2pR(0  ,g))
     end if
     
     if (BC(2,2,g) > 0) then
        cGp2R( :,jmax) =  0.
        cDv2R( :,jmax) =  0.
        !cDv2R(-1,jmax) = -1./(x2v (N2      ) - x2v (N2-1  )) ! TEST!!! Das geht in die Hose ...
        cDv2R(-1,jmax) = -1./(y2v (M2      ) - y2v (M2-1  ))
     else
        cGp2R( 0,jmax) = -1./(x2pR(jmax+1,g) - x2pR(jmax,g))
        cGp2R( 1,jmax) =  1./(x2pR(jmax+1,g) - x2pR(jmax,g))
     end if
     !--------------------------------------------------------------------------------------------------------
     do j = 1, jmax
        cdg2(-1,j,g) = cDv2R(-1,j)*cGp2R(0,j-1)
        cdg2( 0,j,g) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
        cdg2( 1,j,g) =                            cDv2R(0,j)*cGp2R(1,j)
     end do
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     if (BC(1,2,g) > 0) cdg2( :,1   ,g) = 2.*cdg2(:,1   ,g)
     if (BC(2,2,g) > 0) cdg2( :,jmax,g) = 2.*cdg2(:,jmax,g)
     
     if (BC(1,2,g) == -2) then
        cdg2( 1,1   ,g) = cdg2( 1,1   ,g) + cdg2(-1,1   ,g)
        cdg2(-1,1   ,g) = 0.
     end if
     if (BC(2,2,g) == -2) then
        cdg2(-1,jmax,g) = cdg2(-1,jmax,g) + cdg2( 1,jmax,g)
        cdg2( 1,jmax,g) = 0.
     end if
     !========================================================================================================
     if (dimens == 3) then
     
     kmax = NN(3,g)
     
     cDw3R = 0.
     cGp3R = 0.
     
     do k = 1, kmax
        cDw3R(-1,k) = -1./(x3wR(k  ,g) - x3wR(k-1,g))
        cDw3R( 0,k) =  1./(x3wR(k  ,g) - x3wR(k-1,g))
     end do
     
     do k = 1, kmax-1
        cGp3R( 0,k) = -1./(x3pR(k+1,g) - x3pR(k  ,g))
        cGp3R( 1,k) =  1./(x3pR(k+1,g) - x3pR(k  ,g))
     end do
     
     if (BC(1,3,g) > 0) then
        cGp3R( :,0   ) =  0.
        cDw3R( :,1   ) =  0.
        !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
        cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     else
        cGp3R( 0,0   ) = -1./(x3pR(1,g) - x3pR(0,g))
        cGp3R( 1,0   ) =  1./(x3pR(1,g) - x3pR(0,g))
     end if
     
     if (BC(2,3,g) > 0) then
        cGp3R( :,kmax) =  0.
        cDw3R( :,kmax) =  0.
        !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
        cDw3R(-1,kmax) = -1./(y3w (M3      ) - y3w (M3-1  ))
     else
        cGp3R( 0,kmax) = -1./(x3pR(kmax+1,g) - x3pR(kmax,g))
        cGp3R( 1,kmax) =  1./(x3pR(kmax+1,g) - x3pR(kmax,g))
     end if
     !--------------------------------------------------------------------------------------------------------
     do k = 1, kmax
        cdg3(-1,k,g) = cDw3R(-1,k)*cGp3R(0,k-1)
        cdg3( 0,k,g) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
        cdg3( 1,k,g) =                            cDw3R(0,k)*cGp3R(1,k)
     end do
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     if (BC(1,3,g) > 0) cdg3( :,1   ,g) = 2.*cdg3(:,1   ,g)
     if (BC(2,3,g) > 0) cdg3( :,kmax,g) = 2.*cdg3(:,kmax,g)
     
     if (BC(1,3,g) == -2) then
        cdg3( 1,1   ,g) = cdg3( 1,1   ,g) + cdg3(-1,1   ,g)
        cdg3(-1,1   ,g) = 0.
     end if
     if (BC(2,3,g) == -2) then
        cdg3(-1,kmax,g) = cdg3(-1,kmax,g) + cdg3( 1,kmax,g)
        cdg3( 1,kmax,g) = 0.
     end if
     
     end if
     !========================================================================================================
     
  end do
  
  
  end subroutine get_stencil
  
  
  
  
  
  
  
  
  
  
  
  subroutine get_stencil_transp
  
  implicit none
  
  integer                ::  i, j, k, g
  real   , allocatable   ::  work(:,:)
  
  
  do g = 1, n_grids
     
     !========================================================================================================
     allocate(work(-1:1,0:(NN(1,g)+1)))
     
     work = 0.
     do i = 1, NN(1,g)
        work(-1:1,i) = cdg1(-1:1,i,g)
     end do
     
     cdg1(:,:,g) = 0.
     do i = 1, NN(1,g)
        cdg1(-1,i,g) = work( 1,i-1)
        cdg1( 0,i,g) = work( 0,i  )
        cdg1( 1,i,g) = work(-1,i+1)
     end do
     
     deallocate(work)
     
     if (BC(1,1,g) == 0 .or. BC(1,1,g) == -1) then
        i = 1
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     end if
     if (BC(2,1,g) == 0 .or. BC(2,1,g) == -1) then
        i = NN(1,g)
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     end if
     !========================================================================================================
     allocate(work(-1:1,0:(NN(2,g)+1)))
     
     work = 0.
     do j = 1, NN(2,g)
        work(-1:1,j) = cdg2(-1:1,j,g)
     end do
     
     cdg2(:,:,g) = 0.
     do j = 1, NN(2,g)
        cdg2(-1,j,g) = work( 1,j-1)
        cdg2( 0,j,g) = work( 0,j  )
        cdg2( 1,j,g) = work(-1,j+1)
     end do
     
     deallocate(work)
     
     if (BC(1,2,g) == 0 .or. BC(1,2,g) == -1) then
        j = 1
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     end if
     if (BC(2,2,g) == 0 .or. BC(2,2,g) == -1) then
        j = NN(2,g)
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     end if
     !========================================================================================================
     if (dimens == 3) then
     
     allocate(work(-1:1,0:(NN(3,g)+1)))
     
     work = 0.
     do k = 1, NN(3,g)
        work(-1:1,k) = cdg3(-1:1,k,g)
     end do
     
     cdg3(:,:,g) = 0.
     do k = 1, NN(3,g)
        cdg3(-1,k,g) = work( 1,k-1)
        cdg3( 0,k,g) = work( 0,k  )
        cdg3( 1,k,g) = work(-1,k+1)
     end do
     
     deallocate(work)
     
     if (BC(1,3,g) == 0 .or. BC(1,3,g) == -1) then
        k = 1
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     end if
     if (BC(2,3,g) == 0 .or. BC(2,3,g) == -1) then
        k = NN(3,g)
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     end if
     
     end if
     !========================================================================================================
     
  end do
  
  
  ! sicherheitshalber (cdg3 = 0. sollte eigentlich schon erf�llt sein):
  if (dimens == 2) cdg3 = 0.
  
  
  end subroutine get_stencil_transp
  
  
  
  
  
  
  
  
  
  
  
  subroutine get_stencil_Helm()
  
  implicit none
  
  integer                ::  imax, jmax, kmax
  integer                ::  g, m
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Bei Symmetrie ergeben sich automatisch auf der Symmetrie-Ebene Dirichlet-RB für die       !
  !                Normalkomponente der Geschwindigkeit (betrifft cu11R, cv22R, cw33R und tritt aufgrund des !
  !                Gittertyps nur auf den gröberen Gittern auf). Die zugehörigen Restriktionskoeffizienten   !
  !                liefern automatisch entsprechend die "Randbedingung" vel_n = 0.                           !
  !              - Bei Symmetrie-RB werden in "CALL diff_coeffs" die Koeffizienten am Rand gemäss Gittertyp  !
  !                manipuliert, d.h. bei cp11R, cp22R, cp33R werden Geschwindigkeitskomponenten parallel     !
  !                zur Symmetrieebene (bzw. Druck) angenommen, so dass cu11R, cv22R, cw33R nach dem Kopieren !
  !                von cp11R, cp22R, cp33R entsprechend für Geschwindigkeitskomponenten normal zur           !
  !                Symmetrie-Ebene geändert werden müssen.                                                   !
  !----------------------------------------------------------------------------------------------------------!
  
  ! sicherheitshalber (dimens == 2, ...):
  cp11R = 0.
  cp22R = 0.
  cp33R = 0.
  
  cu11R = 0.
  cv22R = 0.
  cw33R = 0.
  
  cc11R = 0.
  cc22R = 0.
  cc33R = 0.
  
  
  !===========================================================================================================
  !=== grobe & feinstes Gitter ===============================================================================
  !===========================================================================================================
  do g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     call diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1pR(-1:imax+1,g),x1pR(-1:imax+1,g),BC_1L,BC_1U,0,imax,-1,1,cp11R(-1,0,g))
     
     if (BC_1L == 1 .or. BC_1L == 3) then
        cp11R( :,1   ,g) =  0.
        cp11R( 0,1   ,g) =  1.
     end if
     if (BC_1U == 1 .or. BC_1U == 3) then
        cp11R( :,imax,g) =  0.
        cp11R( 0,imax,g) =  1.
     end if
     
     if (BC_1L == 2 .or. BC_1L == 4) then
        cp11R(-1,1   ,g) =  0.
        cp11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cp11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     end if
     if (BC_1U == 2 .or. BC_1U == 4) then
        cp11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 1,imax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     cu11R(:,:,g) = cp11R(:,:,g)
     
     if (BC_1L ==  2) then
        cu11R( :,1   ,g) =  0.
        cu11R( 0,1   ,g) =  1.
     end if
     if (BC_1U ==  2) then
        cu11R( :,imax,g) =  0.
        cu11R( 0,imax,g) =  1.
     end if
     
     if (BC_1L ==  3) then
        cu11R(-1,1   ,g) =  0.
        cu11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cu11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     end if
     if (BC_1U ==  3) then
        cu11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 1,imax,g) =  0.
     end if
     
     if (BC_1L == -2) then
        cu11R(-1,1   ,g) =  0.
        cu11R( 1,1   ,g) =  0.
     end if
     if (BC_1U == -2) then
        cu11R(-1,imax,g) =  0.
        cu11R( 1,imax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     do m = 1, n_conc
        cc11R(:,:,g,m) = cp11R(:,:,g)
        
        if (BCc_1L(m) > 0) then
           cc11R(-1,1   ,g,m) =  0.
           cc11R( 0,1   ,g,m) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
           cc11R( 1,1   ,g,m) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
        end if
        if (BCc_1U(m) > 0) then
           cc11R(-1,imax,g,m) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
           cc11R( 0,imax,g,m) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
           cc11R( 1,imax,g,m) =  0.
        end if
     end do
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     jmax = NN(2,g)
     
     call diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2pR(-1:jmax+1,g),x2pR(-1:jmax+1,g),BC_2L,BC_2U,0,jmax,-1,1,cp22R(-1,0,g))
     
     if (BC_2L == 1 .or. BC_2L == 3) then
        cp22R(:,1    ,g) =  0.
        cp22R(0,1    ,g) =  1.
     end if
     if (BC_2U == 1 .or. BC_2U == 3) then
        cp22R( :,jmax,g) =  0.
        cp22R( 0,jmax,g) =  1.
     end if
     
     if (BC_2L == 2 .or. BC_2L == 4) then
        cp22R(-1,1   ,g) =  0.
        cp22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cp22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     end if
     if (BC_2U == 2 .or. BC_2U == 4) then
        cp22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 1,jmax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     cv22R(:,:,g) = cp22R(:,:,g)
     
     if (BC_2L ==  2) then
        cv22R( :,1   ,g) =  0.
        cv22R( 0,1   ,g) =  1.
     end if
     if (BC_2U ==  2) then
        cv22R( :,jmax,g) =  0.
        cv22R( 0,jmax,g) =  1.
     end if
     
     if (BC_2L ==  3) then
        cv22R(-1,1   ,g) =  0.
        cv22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cv22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     end if
     if (BC_2U ==  3) then
        cv22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 1,jmax,g) =  0.
     end if
     
     if (BC_2L == -2) then
        cv22R(-1,1   ,g) =  0.
        cv22R( 1,1   ,g) =  0.
     end if
     if (BC_2U == -2) then
        cv22R(-1,jmax,g) =  0.
        cv22R( 1,jmax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     do m = 1, n_conc
        cc22R(:,:,g,m) = cp22R(:,:,g)
        
        if (BCc_2L(m) > 0) then
           cc22R(-1,1   ,g,m) =  0.
           cc22R( 0,1   ,g,m) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
           cc22R( 1,1   ,g,m) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
        end if
        if (BCc_2U(m) > 0) then
           cc22R(-1,jmax,g,m) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
           cc22R( 0,jmax,g,m) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
           cc22R( 1,jmax,g,m) =  0.
        end if
     end do
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     
     kmax = NN(3,g)
     
     call diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3pR(-1:kmax+1,g),x3pR(-1:kmax+1,g),BC_3L,BC_3U,0,kmax,-1,1,cp33R(-1,0,g))
     
     if (BC_3L == 1 .or. BC_3L == 3) then
        cp33R( :,1   ,g) =  0.
        cp33R( 0,1   ,g) =  1.
     end if
     if (BC_3U == 1 .or. BC_3U == 3) then
        cp33R( :,kmax,g) =  0.
        cp33R( 0,kmax,g) =  1.
     end if
     
     if (BC_3L == 2 .or. BC_3L == 4) then
        cp33R(-1,1   ,g) =  0.
        cp33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cp33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     end if
     if (BC_3U == 2 .or. BC_3U == 4) then
        cp33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 1,kmax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     cw33R(:,:,g) = cp33R(:,:,g)
     
     if (BC_3L ==  2) then
        cw33R( :,1   ,g) =  0.
        cw33R( 0,1   ,g) =  1.
     end if
     if (BC_3U ==  2) then
        cw33R( :,kmax,g) =  0.
        cw33R( 0,kmax,g) =  1.
     end if
     
     if (BC_3L ==  3) then
        cw33R(-1,1   ,g) =  0.
        cw33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cw33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     end if
     if (BC_3U ==  3) then
        cw33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 1,kmax,g) =  0.
     end if
     
     if (BC_3L == -2) then
        cw33R(-1,1   ,g) =  0.
        cw33R( 1,1   ,g) =  0.
     end if
     if (BC_3U == -2) then
        cw33R(-1,kmax,g) =  0.
        cw33R( 1,kmax,g) =  0.
     end if
     !--------------------------------------------------------------------------------------------------------
     do m = 1, n_conc
        cc33R(:,:,g,m) = cp33R(:,:,g)
        
        if (BCc_3L(m) > 0) then
           cc33R(-1,1   ,g,m) =  0.
           cc33R( 0,1   ,g,m) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
           cc33R( 1,1   ,g,m) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
        end if
        if (BCc_3U(m) > 0) then
           cc33R(-1,kmax,g,m) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
           cc33R( 0,kmax,g,m) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
           cc33R( 1,kmax,g,m) =  0.
        end if
     end do
     
     end if
     !========================================================================================================
     
  end do
  
  g = 1
  
  !===========================================================================================================
  !=== feinstes Gitter (Tangential-Koeffizienten überschreiben) ==============================================
  !===========================================================================================================
  call diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1u(-1:N1+1),x1u(-1:N1+1),BC_1L,BC_1U,1,N1,-1,1,cu11R(-1,0,g))
  
  if (BC_1L == 1 .or. BC_1L == 2) then
     cu11R(-1,0 ,g) = 0.
     cu11R( 0,0 ,g) = 1.- (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =     (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
  end if
  if (BC_1U == 1 .or. BC_1U == 2) then
     cu11R(-1,N1,g) = 1.- (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =     (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) = 0.
  end if
  
  if (BC_1L == 3 .or. BC_1L == 4) then
     cu11R(-1,0 ,g) =  0.
     cu11R( 0,0 ,g) = -1./(x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =  1./(x1u(1 )-x1u(0   ))
  end if
  if (BC_1U == 3 .or. BC_1U == 4) then
     cu11R(-1,N1,g) = -1./(x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =  1./(x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) =  0.
  end if
  !-----------------------------------------------------------------------------------------------------------
  call diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2v(-1:N2+1),x2v(-1:N2+1),BC_2L,BC_2U,1,N2,-1,1,cv22R(-1,0,g))
  
  if (BC_2L == 1 .or. BC_2L == 2) then
     cv22R(-1,0 ,g) = 0.
     cv22R( 0,0 ,g) = 1.- (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =     (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
  end if
  if (BC_2U == 1 .or. BC_2U == 2) then
     cv22R(-1,N2,g) = 1.- (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =     (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) = 0.
  end if
  
  if (BC_2L == 3 .or. BC_2L == 4) then
     cv22R(-1,0 ,g) =  0.
     cv22R( 0,0 ,g) = -1./(x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =  1./(x2v(1 )-x2v(0   ))
  end if
  if (BC_2U == 3 .or. BC_2U == 4) then
     cv22R(-1,N2,g) = -1./(x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =  1./(x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) =  0.
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  
  call diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3w(-1:N3+1),x3w(-1:N3+1),BC_3L,BC_3U,1,N3,-1,1,cw33R(-1,0,g))
  
  if (BC_3L == 1 .or. BC_3L == 2) then
     cw33R(-1,0 ,g) = 0.
     cw33R( 0,0 ,g) = 1.- (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =     (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
  end if
  if (BC_3U == 1 .or. BC_3U == 2) then
     cw33R(-1,N3,g) = 1.- (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =     (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) = 0.
  end if
  
  if (BC_3L == 3 .or. BC_3L == 4) then
     cw33R(-1,0 ,g) =  0.
     cw33R( 0,0 ,g) = -1./(x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =  1./(x3w(1 )-x3w(0   ))
  end if
  if (BC_3U == 3 .or. BC_3U == 4) then
     cw33R(-1,N3,g) = -1./(x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =  1./(x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) =  0.
  end if
  
  end if
  !===========================================================================================================
  
  
  end subroutine get_stencil_Helm
  
  
  
  
  
  
  
  
  
  
  
  subroutine interp_coeffs()
  
  implicit none
  
  integer               ::  i, i0, di
  integer               ::  j, j0, dj
  integer               ::  k, k0, dk
  integer               ::  g
  
  real                  ::  Dx12, Dx10
  
  
  cI1 = 0.
  cI2 = 0.
  cI3 = 0.
  
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  do g = 1, n_grids-1
     
     di = (M1-1)/((NN(1,g)-1)*NB(1,g))
     dj = (M2-1)/((NN(2,g)-1)*NB(2,g))
     dk = (M3-1)/((NN(3,g)-1)*NB(3,g))
     
     !--------------------------------------------------------------------------------------------------------
     do i = 2, NN(1,g)-1, 2
        i0 = 1 + (i-1)*di + (M1-1)*(iB(1,g)-1)/NB(1,g)
        
        Dx10 = y1p(i0   )-y1p(i0-di)
        Dx12 = y1p(i0+di)-y1p(i0-di)
        
        cI1(1,i,g) = 1.- Dx10/Dx12
        cI1(2,i,g) =     Dx10/Dx12
     end do
     !--------------------------------------------------------------------------------------------------------
     do j = 2, NN(2,g)-1, 2
        j0 = 1 + (j-1)*dj + (M2-1)*(iB(2,g)-1)/NB(2,g)
        
        Dx10 = y2p(j0   )-y2p(j0-dj)
        Dx12 = y2p(j0+dj)-y2p(j0-dj)
        
        cI2(1,j,g) = 1.- Dx10/Dx12
        cI2(2,j,g) =     Dx10/Dx12
     end do
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     do k = 2, NN(3,g)-1, 2
        k0 = 1 + (k-1)*dk + (M3-1)*(iB(3,g)-1)/NB(3,g)
        
        Dx10 = y3p(k0   )-y3p(k0-dk)
        Dx12 = y3p(k0+dk)-y3p(k0-dk)
        
        cI3(1,k,g) = 1.- Dx10/Dx12
        cI3(2,k,g) =     Dx10/Dx12
     end do
     end if
     !--------------------------------------------------------------------------------------------------------
  end do
  !===========================================================================================================
  
  
  end subroutine interp_coeffs
  
  
  
  
  
  
  
  
  
  
  
  subroutine interp_coeffs_Helm()
  
  implicit none
  
  integer               ::  i, j, k
  real                  ::  Dx12, Dx1a, Dx1b
  
  
  cIH1 = 0.
  cIH2 = 0.
  cIH3 = 0.
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  do i = 1, N1-2, 2
     Dx1a = x1u(i  )-x1p(i)
     Dx1b = x1u(i+1)-x1p(i)
     Dx12 = x1p(i+2)-x1p(i)
     
     cIH1(1,i  ) = 1.- Dx1a/Dx12
     cIH1(2,i  ) =     Dx1a/Dx12
     
     cIH1(1,i+1) = 1.- Dx1b/Dx12
     cIH1(2,i+1) =     Dx1b/Dx12
  end do
  
  !--- Randbedingungen (Extrapolation) ---
  if (BC_1L > 0) then
     Dx1a = x1u(0)-x1p(1)
     Dx12 = x1p(3)-x1p(1)
     
     cIH1(1,0) = 1.- Dx1a/Dx12
     cIH1(2,0) =     Dx1a/Dx12
  end if
  if (BC_1U > 0) then
     Dx1a = x1u(N1)-x1p(N1-2)
     Dx12 = x1p(N1)-x1p(N1-2)
     
     cIH1(1,N1) = 1.- Dx1a/Dx12
     cIH1(2,N1) =     Dx1a/Dx12
  end if
  !-----------------------------------------------------------------------------------------------------------
  do j = 1, N2-2, 2
     Dx1a = x2v(j  )-x2p(j)
     Dx1b = x2v(j+1)-x2p(j)
     Dx12 = x2p(j+2)-x2p(j)
     
     cIH2(1,j  ) = 1.- Dx1a/Dx12
     cIH2(2,j  ) =     Dx1a/Dx12
     
     cIH2(1,j+1) = 1.- Dx1b/Dx12
     cIH2(2,j+1) =     Dx1b/Dx12
  end do
  
  !--- Randbedingungen (Extrapolation) ---
  if (BC_2L > 0) then
     Dx1a = x2v(0)-x2p(1)
     Dx12 = x2p(3)-x2p(1)
     
     cIH2(1,0) = 1.- Dx1a/Dx12
     cIH2(2,0) =     Dx1a/Dx12
  end if
  if (BC_2U > 0) then
     Dx1a = x2v(N2)-x2p(N2-2)
     Dx12 = x2p(N2)-x2p(N2-2)
     
     cIH2(1,N2) = 1.- Dx1a/Dx12
     cIH2(2,N2) =     Dx1a/Dx12
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  
  do k = 1, N3-2, 2
     Dx1a = x3w(k  )-x3p(k)
     Dx1b = x3w(k+1)-x3p(k)
     Dx12 = x3p(k+2)-x3p(k)
     
     cIH3(1,k  ) = 1.- Dx1a/Dx12
     cIH3(2,k  ) =     Dx1a/Dx12
     
     cIH3(1,k+1) = 1.- Dx1b/Dx12
     cIH3(2,k+1) =     Dx1b/Dx12
  end do
  
  !--- Randbedingungen (Extrapolation) ---
  if (BC_3L > 0) then
     Dx1a = x3w(0)-x3p(1)
     Dx12 = x3p(3)-x3p(1)
     
     cIH3(1,0) = 1.- Dx1a/Dx12
     cIH3(2,0) =     Dx1a/Dx12
  end if
  if (BC_3U > 0) then
     Dx1a = x3w(N3)-x3p(N3-2)
     Dx12 = x3p(N3)-x3p(N3-2)
     
     cIH3(1,N3) = 1.- Dx1a/Dx12
     cIH3(2,N3) =     Dx1a/Dx12
  end if
  
  end if
  !===========================================================================================================
  
  
  end subroutine interp_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  
  subroutine restr_coeffs()
  
  implicit none
  
  integer               ::  g
  integer               ::  i, ii, iimax
  integer               ::  j, jj, jjmax
  integer               ::  k, kk, kkmax
  
  integer               ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U ! TEST!!!
  
  
  cR1 = 0.
  cR2 = 0.
  cR3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  do g = 2, n_grids
     
     iimax = (NN(1,g)-1)/n_gather(1,g)+1
     jjmax = (NN(2,g)-1)/n_gather(2,g)+1
     kkmax = (NN(3,g)-1)/n_gather(3,g)+1
     
     !iimax = (NN(1,g)-1)*NB(1,g)/NB(1,g-1)+1 ! TEST!!! alt ...
     !jjmax = (NN(2,g)-1)*NB(2,g)/NB(2,g-1)+1
     !kkmax = (NN(3,g)-1)*NB(3,g)/NB(3,g-1)+1
     !--------------------------------------------------------------------------------------------------------
     if (n_gather(1,g) > 1) then
        BC_1L = BC(1,1,g-1) ! TEST!!! evtl. wieder auf dem feinen Gitter speichern ...
        BC_1U = BC(2,1,g-1)
     else
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
     end if
     !--------------------------------------------------------------------------------------------------------
     if (n_gather(2,g) > 1) then
        BC_2L = BC(1,2,g-1)
        BC_2U = BC(2,2,g-1)
     else
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
     end if
     !--------------------------------------------------------------------------------------------------------
     if (n_gather(3,g) > 1) then
        BC_3L = BC(1,3,g-1)
        BC_3U = BC(2,3,g-1)
     else
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
     end if
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     do ii = 1, iimax
        cR1(-1,ii,g) = 1./4.
        cR1( 0,ii,g) = 2./4.
        cR1( 1,ii,g) = 1./4.
     end do
     
     if (BC_1L > 0) then
        cR1(-1,1,g) = 0.
        cR1( 0,1,g) = 1.
        cR1( 1,1,g) = 0.
        
        cR1(-1,2,g) = 0. ! TEST!!! Sollte evtl. noch ergaenzt werden ...
        cR1( 0,2,g) = 1.
        cR1( 1,2,g) = 0.
     end if
     if (BC_1L == -2) then
        cR1( 1,1,g) = cR1( 1,1,g) + cR1(-1,1,g)
        cR1(-1,1,g) = 0.
     end if
     
     if (BC_1U > 0) then
        cR1(-1,iimax  ,g) = 0.
        cR1( 0,iimax  ,g) = 1.
        cR1( 1,iimax  ,g) = 0.
        
        cR1(-1,iimax-1,g) = 0. ! TEST!!!
        cR1( 0,iimax-1,g) = 1.
        cR1( 1,iimax-1,g) = 0.
     end if
     if (BC_1U == -2) then
        cR1(-1,iimax,g) = cR1( 1,iimax,g) + cR1(-1,iimax,g)
        cR1( 1,iimax,g) = 0.
     end if
     !--------------------------------------------------------------------------------------------------------
     do jj = 1, jjmax
        cR2(-1,jj,g) = 1./4.
        cR2( 0,jj,g) = 2./4.
        cR2( 1,jj,g) = 1./4.
     end do
     
     if (BC_2L > 0) then
        cR2(-1,1,g) = 0.
        cR2( 0,1,g) = 1.
        cR2( 1,1,g) = 0.
        
        cR2(-1,2,g) = 0. ! TEST!!!
        cR2( 0,2,g) = 1.
        cR2( 1,2,g) = 0.
     end if
     if (BC_2L == -2) then
        cR2( 1,1,g) = cR2( 1,1,g) + cR2(-1,1,g)
        cR2(-1,1,g) = 0.
     end if
     
     if (BC_2U > 0) then
        cR2(-1,jjmax  ,g) = 0.
        cR2( 0,jjmax  ,g) = 1.
        cR2( 1,jjmax  ,g) = 0.
        
        cR2(-1,jjmax-1,g) = 0. ! TEST!!!
        cR2( 0,jjmax-1,g) = 1.
        cR2( 1,jjmax-1,g) = 0.
     end if
     if (BC_2U == -2) then
        cR2(-1,jjmax,g) = cR2( 1,jjmax,g) + cR2(-1,jjmax,g)
        cR2( 1,jjmax,g) = 0.
     end if
     !--------------------------------------------------------------------------------------------------------
     if (dimens == 3) then
     do kk = 1, kkmax
        cR3(-1,kk,g) = 1./4.
        cR3( 0,kk,g) = 2./4.
        cR3( 1,kk,g) = 1./4.
     end do
     
     if (BC_3L > 0) then
        cR3(-1,1,g) = 0.
        cR3( 0,1,g) = 1.
        cR3( 1,1,g) = 0.
        
        cR3(-1,2,g) = 0. ! TEST!!!
        cR3( 0,2,g) = 1.
        cR3( 1,2,g) = 0.
     end if
     if (BC_3L == -2) then
        cR3( 1,1,g) = cR3( 1,1,g) + cR3(-1,1,g)
        cR3(-1,1,g) = 0.
     end if
     
     if (BC_3U > 0) then
        cR3(-1,kkmax  ,g) = 0.
        cR3( 0,kkmax  ,g) = 1.
        cR3( 1,kkmax  ,g) = 0.
        
        cR3(-1,kkmax-1,g) = 0. ! TEST!!!
        cR3( 0,kkmax-1,g) = 1.
        cR3( 1,kkmax-1,g) = 0.
     end if
     if (BC_3U == -2) then
        cR3(-1,kkmax,g) = cR3( 1,kkmax,g) + cR3(-1,kkmax,g)
        cR3( 1,kkmax,g) = 0.
     end if
     
     end if
     !--------------------------------------------------------------------------------------------------------
  end do
  !===========================================================================================================
  
  
  if (1 == 2) then ! TEST!!!
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  cRest1 = 0.
  cRest2 = 0.
  cRest3 = 0.
  
  do g = 1, n_grids-1
     iimax = NN(1,g)
     jjmax = NN(2,g)
     kkmax = NN(3,g)
     
     if (ncb1f(dim_ncb1c) < iimax) then ! TEST!!! Quick and dirty ...
        call diff_coeffs(1,0,0,.false.,dim_ncb1c,ncb1r,x1pR(b1L,g),x1pR(b1L,g),BC_1L,BC_1U,0,iimax,b1L,b1U,cRest1(b1L,0,g))
     else
        do i = 1, iimax
           cRest1(-1:1,i,g) = (/0.25,0.5,0.25/) ! TEST!!! angepasste Koeffizienten? (mapping_yes = .FALSE. waere besser ...)
        end do
        if (BC_1L > 0) cRest1(-1:1,1    ,g) = (/0. ,1. ,0. /)
        if (BC_1L ==  -2) cRest1(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        if (BC_1U > 0) cRest1(-1:1,iimax,g) = (/0. ,1. ,0. /)
        if (BC_1U ==  -2) cRest1(-1:1,iimax,g) = (/0.5,0.5,0. /)
     end if
     
     
     if (ncb2f(dim_ncb2c) < jjmax) then
        call diff_coeffs(2,0,0,.false.,dim_ncb2c,ncb2r,x2pR(b2L,g),x2pR(b2L,g),BC_2L,BC_2U,0,jjmax,b2L,b2U,cRest2(b2L,0,g))
     else
        do j = 1, jjmax
           cRest2(-1:1,j,g) = (/0.25,0.5,0.25/)
        end do
        if (BC_2L > 0) cRest2(-1:1,1    ,g) = (/0. ,1. ,0. /)
        if (BC_2L ==  -2) cRest2(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        if (BC_2U > 0) cRest2(-1:1,jjmax,g) = (/0. ,1. ,0. /)
        if (BC_2U ==  -2) cRest2(-1:1,jjmax,g) = (/0.5,0.5,0. /)
     end if
     
     
     if (dimens == 3) then
     if (ncb3f(dim_ncb3c) < kkmax) then
        call diff_coeffs(3,0,0,.false.,dim_ncb3c,ncb3r,x3pR(b3L,g),x3pR(b3L,g),BC_3L,BC_3U,0,kkmax,b3L,b3U,cRest3(b3L,0,g))
     else
        do k = 1, kkmax
           cRest3(-1:1,k,g) = (/0.25,0.5,0.25/)
        end do
        if (BC_3L > 0) cRest3(-1:1,1    ,g) = (/0. ,1. ,0. /)
        if (BC_3L ==  -2) cRest3(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        if (BC_3U > 0) cRest3(-1:1,kkmax,g) = (/0. ,1. ,0. /)
        if (BC_3U ==  -2) cRest3(-1:1,kkmax,g) = (/0.5,0.5,0. /)
     end if
     end if
     
  end do
  !===========================================================================================================
  end if
  
  
  end subroutine restr_coeffs
  
  
  
  
  
  
  
  
  
  
  
  subroutine restr_coeffs_Helm()
  
  implicit none
  
  integer               ::  i, j, k
  real                  ::  Dx12, Dx1a
  
  
  cRH1 = 0.
  cRH2 = 0.
  cRH3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  do i = 1, N1, 2
     Dx1a = x1p(i)-x1u(i-1)
     Dx12 = x1u(i)-x1u(i-1)
     
     cRH1(1,i) = 1.- Dx1a/Dx12
     cRH1(2,i) =     Dx1a/Dx12
  end do
  
  if (BC_1L > 0) then
     cRH1(1,1 ) = 1.
     cRH1(2,1 ) = 0.
  end if
  if (BC_1L == -2) then
     cRH1(:,1 ) = 0.
  end if
  
  if (BC_1U > 0) then
     cRH1(1,N1) = 0.
     cRH1(2,N1) = 1.
  end if
  if (BC_1U == -2) then
     cRH1(:,N1) = 0.
  end if
  !-----------------------------------------------------------------------------------------------------------
  do j = 1, N2, 2
     Dx1a = x2p(j)-x2v(j-1)
     Dx12 = x2v(j)-x2v(j-1)
     
     cRH2(1,j) = 1.- Dx1a/Dx12
     cRH2(2,j) =     Dx1a/Dx12
  end do
  
  if (BC_2L > 0) then
     cRH2(1,1 ) = 1.
     cRH2(2,1 ) = 0.
  end if
  if (BC_2L == -2) then
     cRH2(:,1 ) = 1.
  end if
  
  if (BC_2U > 0) then
     cRH2(1,N2) = 0.
     cRH2(2,N2) = 1.
  end if
  if (BC_2U == -2) then
     cRH2(:,N2) = 0.
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
  
  do k = 1, N3, 2
     Dx1a = x3p(k)-x3w(k-1)
     Dx12 = x3w(k)-x3w(k-1)
     
     cRH3(1,k) = 1.- Dx1a/Dx12
     cRH3(2,k) =     Dx1a/Dx12
  end do
  
  if (BC_3L > 0) then
     cRH3(1,1 ) = 1.
     cRH3(2,1 ) = 0.
  end if
  if (BC_3L == -2) then
     cRH3(:,1 ) = 1.
  end if
  
  if (BC_3U > 0) then
     cRH3(1,N3) = 0.
     cRH3(2,N3) = 1.
  end if
  if (BC_3U == -2) then
     cRH3(:,N3) = 0.
  end if
  
  end if
  !===========================================================================================================
  
  
  end subroutine restr_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  
  subroutine get_weights() ! TEST!!! bezieht sich bislang nur auf explizite Differenzen ... ODER: rauswerfen!
  
  ! revised: 24.10.07
  
  implicit none
  
  integer                ::  i, ii
  integer                ::  j, jj
  integer                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Diese Abschaetzung ist im Feld etwas konservativer (Faktor 1.5-2.5) im Vergleich zu       !
  !                Differenzen 2. Konvergenzordnung. Am Rand dagegen sehr viel strikter!                     !
  !              - dx_iL, dx_iU dürfen nicht der Ableitung der Mapping-Funktion an der Wand gleichgesetzt    !
  !                werden, da diese im Extremfall auch verschwinden kann.                                    !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
  if (1 == 1) then ! TEST!!! Rand muesste noch extra behandelt werden!
  !===========================================================================================================
  !=== Gewichte für überprüfung der Divergenzfreiheit ========================================================
  !===========================================================================================================
  weight = 0.
  !-----------------------------------------------------------------------------------------------------------
  if (dimens == 3) then
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              do ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              end do
              do jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              end do
              do kk = d3L, d3U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
              end do
              weight(i,j,k) = 1./weight(i,j,k)
           end do
        end do
     end do
     
  !-----------------------------------------------------------------------------------------------------------
  else
     
     do k = S3p, N3p
        do j = S2p, N2p
           do i = S1p, N1p
              do ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              end do
              do jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              end do
              weight(i,j,k) = 1./weight(i,j,k)
           end do
        end do
     end do
     
  end if
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  else
     weight = 1.
  end if
  
  if (1 == 1) then ! TEST!!!
  !-----------------------------------------------------------------------------------------------------------
  if (BC_1L > 0) then
     i = 1
     do k = S3p, N3p
        do j = S2p, N2p
           do ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  if (BC_1U > 0) then
     i = N1
     do k = S3p, N3p
        do j = S2p, N2p
           do ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BC_2L > 0) then
     j = 1
     do k = S3p, N3p
        do i = S1p, N1p
           do jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  if (BC_2U > 0) then
     j = N2
     do k = S3p, N3p
        do i = S1p, N1p
           do jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  if (BC_3L > 0) then
     k = 1
     do j = S2p, N2p
        do i = S1p, N1p
           do kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  if (BC_3U > 0) then
     k = N3
     do j = S2p, N2p
        do i = S1p, N1p
           do kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           end do
           weight(i,j,k) = 1./weight(i,j,k)
        end do
     end do
  end if
  !-----------------------------------------------------------------------------------------------------------
  end if
  
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  if (BC_1L > 0 .and. BC_2L > 0) weight(1 ,1 ,1:N3) = 0.
  if (BC_1L > 0 .and. BC_2U > 0) weight(1 ,N2,1:N3) = 0.
  if (BC_1U > 0 .and. BC_2L > 0) weight(N1,1 ,1:N3) = 0.
  if (BC_1U > 0 .and. BC_2U > 0) weight(N1,N2,1:N3) = 0.
  
  if (BC_1L > 0 .and. BC_3L > 0) weight(1 ,1:N2,1 ) = 0.
  if (BC_1L > 0 .and. BC_3U > 0) weight(1 ,1:N2,N3) = 0.
  if (BC_1U > 0 .and. BC_3L > 0) weight(N1,1:N2,1 ) = 0.
  if (BC_1U > 0 .and. BC_3U > 0) weight(N1,1:N2,N3) = 0.
  
  if (BC_2L > 0 .and. BC_3L > 0) weight(1:N1,1 ,1 ) = 0.
  if (BC_2L > 0 .and. BC_3U > 0) weight(1:N1,1 ,N3) = 0.
  if (BC_2U > 0 .and. BC_3L > 0) weight(1:N1,N2,1 ) = 0.
  if (BC_2U > 0 .and. BC_3U > 0) weight(1:N1,N2,N3) = 0.
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Relaxation der Konzentration ==========================================================================
  !===========================================================================================================
  if (concentration_yes) then
     !---------------------------------------------------------------!
     ! Andere Variante (aus Differenzenquotienten der 1. Ableitung): !
     ! d/dx|_i ~= c1(:,i) ==> dj/dx|_i ~= j*c1(i+j,i)                !
     !      dj := 1       ==> 1./dx|_i ~= j*c1(i+j,i) ==> dx|_i      !
     !---------------------------------------------------------------!
     
     ! ACHTUNG!!! usReSc(m,conc_nu) wird wegen der (m,conc_nu)-Abh�ngigkeit nicht verwendet ...
     ! zum anderen ist jetzt cc11R von conc_nu abh�ngig!!
     !idx_1L = cc11R(0,1 ,1) !- usReSc(m,conc_nu)
     !idx_1U = cc11R(0,N1,1) !- usReSc(m,conc_nu)
     !idx_2L = cc22R(0,1 ,1) !- usReSc(m,conc_nu)
     !idx_2U = cc22R(0,N2,1) !- usReSc(m,conc_nu)
     !idx_3L = cc33R(0,1 ,1) !- usReSc(m,conc_nu)
     !idx_3U = cc33R(0,N3,1) !- usReSc(m,conc_nu)
     
     ! TEST!!! kann jetzt direkt mit cc11R verbunden werden! Skalierung wird aber sicherheitshalber vorerst nicht verwendet:
     idx_1L = 1.
     idx_1U = 1.
     idx_2L = 1.
     idx_2U = 1.
     idx_3L = 1.
     idx_3U = 1.
     
     dx_1L = 1./ idx_1L
     dx_1U = 1./ idx_1U
     dx_2L = 1./ idx_2L
     dx_2U = 1./ idx_2U
     dx_3L = 1./ idx_3L
     dx_3U = 1./ idx_3U
  end if
  !===========================================================================================================
  
  
  end subroutine get_weights
  
  
  
  
  
  
  
  
  
  
  
  subroutine Matrix_invert(N,matrix,matrix_inv)
  
  implicit none
  
  integer, intent(in)  ::  N
  real   , intent(in)  ::  matrix     (1:N,1:N)
  real   , intent(out) ::  matrix_inv (1:N,1:N)
  
  real                 ::  matrix_left(1:N,1:N)
  real                 ::  mult1, mult2
  real                 ::  eps
  integer              ::  i, j, k
  
  real                 ::  store
  integer              ::  maxValue, maxValuePos
  
  
  
  eps = 10.**(-20) ! double precision
  !eps = 10.**(-??) ! single precision
  
  matrix_left = matrix
  
  
  ! Rechte Seite initialisieren (= Einheitsmatrix):
  matrix_inv = 0.
  
  do i = 1, N
     matrix_inv(i,i) = 1.
  end do  
  
  
  !--- Vorwaertsschritt -------------------------------------
  ! linke Seite umformen in obere Dreiecksmatrix
  ! rechte Seite umformen in untere Dreiecksmatrix
  
  do j = 1, N-1
     
     ! Pivoting == Umsortieren der aktuellen Untermatrix (j:N,j:N)
     ! (Diagonalelemente der linken Dreiecksmatrix sind betragsmaessig zu maximieren)
     ! Groesster Wert in Spalte j = zukuenftiges Diagonalelement j,j
     maxValue    = ABS(matrix_left(j,j))
     maxValuePos = j
     do i = j+1, N
        if (ABS(matrix_left(i,j)) > maxValue) then
           maxValue    = ABS(matrix_left(i,j))
           maxValuePos = i
        end if
     end do
     
     ! Zeilen vertauschen:
     if (maxValuePos /= j) then
        do i = 1, N
           store                      = matrix_left(maxValuePos,i)
           matrix_left(maxValuePos,i) = matrix_left(j,i)
           matrix_left(j,i)           = store
           
           store                      = matrix_inv(maxValuePos,i)
           matrix_inv(maxValuePos,i)  = matrix_inv(j,i)
           matrix_inv(j,i)            = store
        end do
     end if
     
     mult1 = matrix_left(j,j)
     if (ABS(mult1) < eps .and. rank == 0) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (1)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', j
     end if
     
     do i = j+1, N
        mult2 = matrix_left(i,j)/mult1
                
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (linke Seite):
        do k = j, N
           ! To avoid underflow:
           if (.not. (ABS(mult2) < eps .and. ABS(matrix_left(j,k)) < eps)) then
              matrix_left(i,k) = matrix_left(i,k) - matrix_left(j,k)*mult2
           end if
        end do
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        do k = 1, N
           ! To avoid underflow:
           if (.not. (ABS(mult2) < eps .and. ABS(matrix_inv(j,k)) < eps)) then
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           end if
        end do
        
        ! Komponente i,j explizit zu Null setzen:
        matrix_left(i,j) = 0.
     end do
     
  end do
  
  
  !--- Rueckwaertsschritt -----------------------------------
  ! linke Seite umformen in Einheitsmatrix
  ! rechte Seite umformen in gesuchte inverse Matrix
  
  do j = N, 2, -1
     
     ! Multiplikator:
     mult1 = matrix_left(j,j)
     if (ABS(mult1) < eps .and. rank == 0) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (2)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', j
     end if
     
     do i = 1, j-1
        mult2 = matrix_left(i,j)/mult1
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        do k = 1, N
           ! To avoid underflow:
           if (.not. (ABS(mult2) < eps .and. ABS(matrix_inv(j,k)) < eps)) then
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           end if
        end do
        
        ! nicht-Diagonalelement explizit zu 0 setzen:
        matrix_left(i,j) = 0.
     end do
     
  end do
  
  
  ! linke Matrix auf Einheitsmatrix bringen:
  do i = 1, N
     mult1 = matrix_left(i,i)
     if (ABS(mult1) < eps .and. rank == 0) then
        write(*,'(a)'      ) 'WARNING! Matrix is probably singular (3)!'
        write(*,'(a,E13.5)') '       ... dividing by', mult1
        write(*,'(a,i2)'   ) '                   i =', i
     end if
     matrix_inv(i,:) = matrix_inv(i,:)/mult1
  end do
  
  
  end subroutine Matrix_invert
  
  
  
  
end module mod_coeffs
