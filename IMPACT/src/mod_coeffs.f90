!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief module providing functions to initiliaze stencil arrays
MODULE mod_coeffs
  
  
  USE mod_dims
  USE mod_vars
  USE mod_exchange
  USE mod_diff ! (apply_compact)
  
  
  PRIVATE
  
  PUBLIC FD_coeffs_compact, init_compact, diff_coeffs_compact, compact_parallel, test_coeffs_compact
  PUBLIC FD_coeffs, test_diff, diff_coeffs, diff_coeffs_exact, test_coeffs
  PUBLIC FD_coeffs_solver, FD_coeffs_solver_integral
  PUBLIC interp_coeffs, interp_coeffs_Helm, restr_coeffs, restr_coeffs_Helm
  PUBLIC get_stencil, get_stencil_Helm, get_stencil_transp
  PUBLIC get_weights
  PUBLIC Matrix_invert
  
  PUBLIC myDGBTRF
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
!pgi$g unroll = n:8
!!pgi$r unroll = n:8
!!pgi$l unroll = n:8
  
  
  
  SUBROUTINE FD_coeffs_compact
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k
  !----------------------------------------------------
  INTEGER             ::  disp1_global(1:NB1), recv1_global(1:NB1), rank_bar1 ! TEST!!!
  INTEGER             ::  disp2_global(1:NB2), recv2_global(1:NB2), rank_bar2
  INTEGER             ::  disp3_global(1:NB3), recv3_global(1:NB3), rank_bar3
  !----------------------------------------------------
  
  
  ! Anmerkung: So ist das sicherlich noch nicht sehr effizient, sollte aber grundsaetzlich die richtige Strategie sein ...
  
  
  ! TEST!!!
  !=====================================================================================================
  disp1 = 0
  disp2 = 0
  disp3 = 0
  
  CALL MPI_COMM_RANK(COMM_BAR1,rank_bar1,merror)
  CALL MPI_COMM_RANK(COMM_BAR2,rank_bar2,merror)
  CALL MPI_COMM_RANK(COMM_BAR3,rank_bar3,merror)
  
  disp1(rank_bar1+1) = 2*ndL*(iB(1,1)-1)*(N2+1)*(N3+1)
  disp2(rank_bar2+1) = 2*ndL*(iB(2,1)-1)*(N1+1)*(N3+1)
  disp3(rank_bar3+1) = 2*ndL*(iB(3,1)-1)*(N1+1)*(N2+1)
  
  CALL MPI_ALLREDUCE(disp1,disp1_global,NB1,MPI_INTEGER,MPI_SUM,COMM_BAR1,merror)
  CALL MPI_ALLREDUCE(disp2,disp2_global,NB2,MPI_INTEGER,MPI_SUM,COMM_BAR2,merror)
  CALL MPI_ALLREDUCE(disp3,disp3_global,NB3,MPI_INTEGER,MPI_SUM,COMM_BAR3,merror)
  
  disp1 = disp1_global
  disp2 = disp2_global
  disp3 = disp3_global
  
  recv1 = 2*ndl*(N2+1)*(N3+1)
  recv2 = 2*ndl*(N1+1)*(N3+1)
  recv3 = 2*ndl*(N1+1)*(N2+1)
  !=====================================================================================================
  
  
  
  
  CALL init_compact(-1,-1)
  
  !=====================================================================================================
  j = 1
  k = 1
  DO i = S1p, N1p
     pp(i,j,k) = x1p(i) - (REAL(i+iShift)-1.0)*L1/REAL(M1-1)
  END DO
  
  CALL exchange(1,0,pp)
  CALL apply_compact(1,0,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,pp,Ap)
  CALL apply_compact(1,0,S12,1,1,N12,1,1,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ar)
  CALL apply_compact(1,0,S12,1,1,N12,1,1,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,pp,rh)
  
  DO i = S11, N11
     dx1GM (i) = 1./(Ap(i,j,k) + L1/REAL(M1-1))
  END DO
  DO i = S12, N12
     dx1pM (i) = 1./(Ar(i,j,k) + L1/REAL(M1-1))
     ddx1pM(i) = -rh(i,j,k)*dx1pM(i)**3
  END DO
  !-----------------------------------------------------------------------------------------------------
  DO i = S11B, N11B
     pp(i,j,k) = x1u(i) - (REAL(i+iShift)-0.5)*L1/REAL(M1-1)
  END DO
  
  CALL exchange(1,1,pp)
  CALL apply_compact(1,1,S1p,1,1,N1p,1,1,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,pp,Ap)
  CALL apply_compact(1,1,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ar)
  CALL apply_compact(1,1,S11,1,1,N11,1,1,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,pp,rh)
  
  DO i = S1p, N1p
     dx1DM (i) = 1./(Ap(i,j,k) + L1/REAL(M1-1))
  END DO
  DO i = S11, N11
     dx1uM (i) = 1./(Ar(i,j,k) + L1/REAL(M1-1))
     ddx1uM(i) = -rh(i,j,k)*dx1uM(i)**3
  END DO
  !=====================================================================================================
  i = 1
  k = 1
  DO j = S2p, N2p
     pp(i,j,k) = x2p(j) - (REAL(j+jShift)-1.0)*L2/REAL(M2-1)
  END DO
  
  CALL exchange(2,0,pp)
  CALL apply_compact(2,0,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,pp,Ap)
  CALL apply_compact(2,0,1,S21,1,1,N21,1,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ar)
  CALL apply_compact(2,0,1,S21,1,1,N21,1,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,pp,rh)
  
  DO j = S22, N22
     dx2GM (j) = 1./(Ap(i,j,k) + L2/REAL(M2-1))
  END DO
  DO j = S21, N21
     dx2pM (j) = 1./(Ar(i,j,k) + L2/REAL(M2-1))
     ddx2pM(j) = -rh(i,j,k)*dx2pM(j)**3
  END DO
  !-----------------------------------------------------------------------------------------------------
  DO j = S22B, N22B
     pp(i,j,k) = x2v(j) - (REAL(j+jShift)-0.5)*L2/REAL(M2-1)
  END DO
  
  CALL exchange(2,2,pp)
  CALL apply_compact(2,2,1,S2p,1,1,N2p,1,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,pp,Ap)
  CALL apply_compact(2,2,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ar)
  CALL apply_compact(2,2,1,S22,1,1,N22,1,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,pp,rh)
  
  DO j = S2p, N2p
     dx2DM (j) = 1./(Ap(i,j,k) + L2/REAL(M2-1))
  END DO
  DO j = S22, N22
     dx2vM (j) = 1./(Ar(i,j,k) + L2/REAL(M2-1))
     ddx2vM(j) = -rh(i,j,k)*dx2vM(j)**3
  END DO
  !=====================================================================================================
  IF (dimens == 3) THEN
  i = 1
  j = 1
  DO k = S3p, N3p
     pp(i,j,k) = x3p(k) - (REAL(k+kShift)-1.0)*L3/REAL(M3-1)
  END DO
  
  CALL exchange(3,0,pp)
  CALL apply_compact(3,0,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,pp,Ap)
  CALL apply_compact(3,0,1,1,S31,1,1,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ar)
  CALL apply_compact(3,0,1,1,S31,1,1,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,pp,rh)
  
  DO k = S33, N33
     dx3GM (k) = 1./(Ap(i,j,k) + L3/REAL(M3-1))
  END DO
  DO k = S31, N31
     dx3pM (k) = 1./(Ar(i,j,k) + L3/REAL(M3-1))
     ddx3pM(k) = -rh(i,j,k)*dx3pM(k)**3
  END DO
  !-----------------------------------------------------------------------------------------------------
  DO k = S33B, N33B
     pp(i,j,k) = x3w(k) - (REAL(k+kShift)-0.5)*L3/REAL(M3-1)
  END DO
  
  CALL exchange(3,3,pp)
  CALL apply_compact(3,3,1,1,S3p,1,1,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,pp,Ap)
  CALL apply_compact(3,3,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ar)
  CALL apply_compact(3,3,1,1,S33,1,1,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,pp,rh)
  
  DO k = S3p, N3p
     dx3DM (k) = 1./(Ap(i,j,k) + L3/REAL(M3-1))
  END DO
  DO k = S33, N33
     dx3wM (k) = 1./(Ar(i,j,k) + L3/REAL(M3-1))
     ddx3wM(k) = -rh(i,j,k)*dx3wM(k)**3
  END DO
  END IF
  !=====================================================================================================
  
  CALL init_compact( 1,-1)
  
  
  END SUBROUTINE FD_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE init_compact(sym_pre,sym_vel)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  sym_pre, sym_vel
  INTEGER                ::  m
  
  
  !--- Koeffizienten bestimmen ---
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S12,N12,BC_1L,BC_1U,0,1,sym_pre,.FALSE.,cp1CL ,cp1CL_LU ,cp1CR )
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S21,N21,BC_2L,BC_2U,0,1,sym_pre,.FALSE.,cp2CL ,cp2CL_LU ,cp2CR )
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S31,N31,BC_3L,BC_3U,0,1,sym_pre,.FALSE.,cp3CL ,cp3CL_LU ,cp3CR )
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,1,1,sym_vel,.FALSE.,cu1CL ,cu1CL_LU ,cu1CR )
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,1,1,sym_vel,.FALSE.,cv2CL ,cv2CL_LU ,cv2CR )
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,1,1,sym_vel,.FALSE.,cw3CL ,cw3CL_LU ,cw3CR )
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S12,N12,BC_1L,BC_1U,0,2,sym_pre,.FALSE.,cp11CL,cp11CL_LU,cp11CR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S21,N21,BC_2L,BC_2U,0,2,sym_pre,.FALSE.,cp22CL,cp22CL_LU,cp22CR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S31,N31,BC_3L,BC_3U,0,2,sym_pre,.FALSE.,cp33CL,cp33CL_LU,cp33CR)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,1,2,sym_vel,.FALSE.,cu11CL,cu11CL_LU,cu11CR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,1,2,sym_vel,.FALSE.,cv22CL,cv22CL_LU,cv22CR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,1,2,sym_vel,.FALSE.,cw33CL,cw33CL_LU,cw33CR)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,1,sym_vel,.FALSE.,cDu1CL,cDu1CL_LU,cDu1CR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,1,sym_vel,.FALSE.,cDv2CL,cDv2CL_LU,cDv2CR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,1,sym_vel,.FALSE.,cDw3CL,cDw3CL_LU,cDw3CR)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,1,sym_vel,.TRUE. ,cDu1CLT,cDu1CLT_LU,cDu1CRT)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,1,sym_vel,.TRUE. ,cDv2CLT,cDv2CLT_LU,cDv2CRT)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,1,sym_vel,.TRUE. ,cDw3CLT,cDw3CLT_LU,cDw3CRT)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,1,sym_pre,.FALSE.,cGp1CL,cGp1CL_LU,cGp1CR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,1,sym_pre,.FALSE.,cGp2CL,cGp2CL_LU,cGp2CR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,1,sym_pre,.FALSE.,cGp3CL,cGp3CL_LU,cGp3CR)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,1,sym_pre,.TRUE. ,cGp1CLT,cGp1CLT_LU,cGp1CRT)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,1,sym_pre,.TRUE. ,cGp2CLT,cGp2CLT_LU,cGp2CRT)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,1,sym_pre,.TRUE. ,cGp3CLT,cGp3CLT_LU,cGp3CRT)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S1p,N1p,BC_1L,BC_1U,3,0,sym_vel,.FALSE.,cIupCL,cIupCL_LU,cIupCR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S2p,N2p,BC_2L,BC_2U,3,0,sym_vel,.FALSE.,cIvpCL,cIvpCL_LU,cIvpCR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S3p,N3p,BC_3L,BC_3U,3,0,sym_vel,.FALSE.,cIwpCL,cIwpCL_LU,cIwpCR)
  
                   CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11,N11,BC_1L,BC_1U,2,0,sym_pre,.FALSE.,cIpuCL,cIpuCL_LU,cIpuCR)
                   CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22,N22,BC_2L,BC_2U,2,0,sym_pre,.FALSE.,cIpvCL,cIpvCL_LU,cIpvCR)
  IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33,N33,BC_3L,BC_3U,2,0,sym_pre,.FALSE.,cIpwCL,cIpwCL_LU,cIpwCR)
  
  
  
  DO m = 1, n_conc
                      CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S1c(m),N1c(m),BCc_1L(m),BCc_1U(m),0,1,sym_pre,.FALSE.,cc1CL (-ndL,0,m),cc1CL_LU (1,0,m),cc1CR (-ndR,0,m))
                      CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S2c(m),N2c(m),BCc_2L(m),BCc_2U(m),0,1,sym_pre,.FALSE.,cc2CL (-ndL,0,m),cc2CL_LU (1,0,m),cc2CR (-ndR,0,m))
     IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S3c(m),N3c(m),BCc_3L(m),BCc_3U(m),0,1,sym_pre,.FALSE.,cc3CL (-ndL,0,m),cc3CL_LU (1,0,m),cc3CR (-ndR,0,m))
     
                      CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S1c(m),N1c(m),BCc_1L(m),BCc_1U(m),0,2,sym_pre,.FALSE.,cc11CL(-ndL,0,m),cc11CL_LU(1,0,m),cc11CR(-ndR,0,m))
                      CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S2c(m),N2c(m),BCc_2L(m),BCc_2U(m),0,2,sym_pre,.FALSE.,cc22CL(-ndL,0,m),cc22CL_LU(1,0,m),cc22CR(-ndR,0,m))
     IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S3c(m),N3c(m),BCc_3L(m),BCc_3U(m),0,2,sym_pre,.FALSE.,cc33CL(-ndL,0,m),cc33CL_LU(1,0,m),cc33CR(-ndR,0,m))
  
                      CALL diff_coeffs_compact(N1,NB1,ndL,ndR,S11   ,N11   ,BCc_1L(m),BCc_1U(m),2,0,sym_pre,.FALSE.,cIcuCL(-ndL,0,m),cIcuCL_LU(1,0,m),cIcuCR(-ndR,0,m))
                      CALL diff_coeffs_compact(N2,NB2,ndL,ndR,S22   ,N22   ,BCc_2L(m),BCc_2U(m),2,0,sym_pre,.FALSE.,cIcvCL(-ndL,0,m),cIcvCL_LU(1,0,m),cIcvCR(-ndR,0,m))
     IF (dimens == 3) CALL diff_coeffs_compact(N3,NB3,ndL,ndR,S33   ,N33   ,BCc_3L(m),BCc_3U(m),2,0,sym_pre,.FALSE.,cIcwCL(-ndL,0,m),cIcwCL_LU(1,0,m),cIcwCR(-ndR,0,m))
  END DO
  
  
  !--- Schur-Komplemente und Hilfsvektoren bestimmen ---
  ! TEST!!!
  ! Im Vergleich zu der auskommentierten Version unten wurde nur die IF-Abfrage modifiziert, da nun IMMER das
  ! Schur-Komplement fuer die Behandlung der Block-Grenzen gerechnet wird! (Pivoting-Problem an Raendern)
                   CALL compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp1CL ,dimS1,BC_1L_global,Sp1 ,Wp1 ,COMM_BAR1)
                   CALL compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp2CL ,dimS2,BC_2L_global,Sp2 ,Wp2 ,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp3CL ,dimS3,BC_3L_global,Sp3 ,Wp3 ,COMM_BAR3)
  
                   CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu1CL ,dimS1,BC_1L_global,Su1 ,Wu1 ,COMM_BAR1)
                   CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv2CL ,dimS2,BC_2L_global,Sv2 ,Wv2 ,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw3CL ,dimS3,BC_3L_global,Sw3 ,Ww3 ,COMM_BAR3)
  
                   CALL compact_parallel(S12,N12,N1,iB(1,1),NB1,ndL,cp11CL,dimS1,BC_1L_global,Sp11,Wp11,COMM_BAR1)
                   CALL compact_parallel(S21,N21,N2,iB(2,1),NB2,ndL,cp22CL,dimS2,BC_2L_global,Sp22,Wp22,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S31,N31,N3,iB(3,1),NB3,ndL,cp33CL,dimS3,BC_3L_global,Sp33,Wp33,COMM_BAR3)
  
                   CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cu11CL,dimS1,BC_1L_global,Su11,Wu11,COMM_BAR1)
                   CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cv22CL,dimS2,BC_2L_global,Sv22,Wv22,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cw33CL,dimS3,BC_3L_global,Sw33,Ww33,COMM_BAR3)
  
                   CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CL,dimS1,BC_1L_global,SDu1,WDu1,COMM_BAR1)
                   CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CL,dimS2,BC_2L_global,SDv2,WDv2,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CL,dimS3,BC_3L_global,SDw3,WDw3,COMM_BAR3)
  
                   CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cDu1CLT,dimS1,BC_1L_global,SDu1T,WDu1T,COMM_BAR1)
                   CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cDv2CLT,dimS2,BC_2L_global,SDv2T,WDv2T,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cDw3CLT,dimS3,BC_3L_global,SDw3T,WDw3T,COMM_BAR3)
  
                   CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CL,dimS1,BC_1L_global,SGp1,WGp1,COMM_BAR1)
                   CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CL,dimS2,BC_2L_global,SGp2,WGp2,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CL,dimS3,BC_3L_global,SGp3,WGp3,COMM_BAR3)
  
                   CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cGp1CLT,dimS1,BC_1L_global,SGp1T,WGp1T,COMM_BAR1)
                   CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cGp2CLT,dimS2,BC_2L_global,SGp2T,WGp2T,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cGp3CLT,dimS3,BC_3L_global,SGp3T,WGp3T,COMM_BAR3)
  
                   CALL compact_parallel(S1p,N1p,N1,iB(1,1),NB1,ndL,cIupCL,dimS1,BC_1L_global,SIup,WIup,COMM_BAR1)
                   CALL compact_parallel(S2p,N2p,N2,iB(2,1),NB2,ndL,cIvpCL,dimS2,BC_2L_global,SIvp,WIvp,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S3p,N3p,N3,iB(3,1),NB3,ndL,cIwpCL,dimS3,BC_3L_global,SIwp,WIwp,COMM_BAR3)
  
                   CALL compact_parallel(S11,N11,N1,iB(1,1),NB1,ndL,cIpuCL,dimS1,BC_1L_global,SIpu,WIpu,COMM_BAR1)
                   CALL compact_parallel(S22,N22,N2,iB(2,1),NB2,ndL,cIpvCL,dimS2,BC_2L_global,SIpv,WIpv,COMM_BAR2)
  IF (dimens == 3) CALL compact_parallel(S33,N33,N3,iB(3,1),NB3,ndL,cIpwCL,dimS3,BC_3L_global,SIpw,WIpw,COMM_BAR3)
  
  
  DO m = 1, n_conc
                      CALL compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc1CL (-ndL,0,m),dimS1,BC_1L_global,Sc1 (1,1,m),Wc1 (1,0,m),COMM_BAR1)
                      CALL compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc2CL (-ndL,0,m),dimS2,BC_2L_global,Sc2 (1,1,m),Wc2 (1,0,m),COMM_BAR2)
     IF (dimens == 3) CALL compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc3CL (-ndL,0,m),dimS3,BC_3L_global,Sc3 (1,1,m),Wc3 (1,0,m),COMM_BAR3)
     
                      CALL compact_parallel(S1c(m),N1c(m),N1,iB(1,1),NB1,ndL,cc11CL(-ndL,0,m),dimS1,BC_1L_global,Sc11(1,1,m),Wc11(1,0,m),COMM_BAR1)
                      CALL compact_parallel(S2c(m),N2c(m),N2,iB(2,1),NB2,ndL,cc22CL(-ndL,0,m),dimS2,BC_2L_global,Sc22(1,1,m),Wc22(1,0,m),COMM_BAR2)
     IF (dimens == 3) CALL compact_parallel(S3c(m),N3c(m),N3,iB(3,1),NB3,ndL,cc33CL(-ndL,0,m),dimS3,BC_3L_global,Sc33(1,1,m),Wc33(1,0,m),COMM_BAR3)
     
                      CALL compact_parallel(S11   ,N11   ,N1,iB(1,1),NB1,ndL,cIcuCL(-ndL,0,m),dimS1,BC_1L_global,SIcu(1,1,m),WIcu(1,0,m),COMM_BAR1)
                      CALL compact_parallel(S22   ,N22   ,N2,iB(2,1),NB2,ndL,cIcvCL(-ndL,0,m),dimS2,BC_2L_global,SIcv(1,1,m),WIcv(1,0,m),COMM_BAR2)
     IF (dimens == 3) CALL compact_parallel(S33   ,N33   ,N3,iB(3,1),NB3,ndL,cIcwCL(-ndL,0,m),dimS3,BC_3L_global,SIcw(1,1,m),WIcw(1,0,m),COMM_BAR3)
  END DO
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
  
  
  END SUBROUTINE init_compact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE diff_coeffs_compact(Nmax,NB,ndL,ndR,SS,NN,BCL,BCU,grid_type,abl,sym,transp_yes,ccL,ccL_LU,ccR)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  Nmax, NB, ndL, ndR, BCL, BCU, SS, NN, grid_type, abl, sym
  LOGICAL, INTENT(IN   ) ::  transp_yes
  REAL   , INTENT(OUT  ) ::  ccL   (-ndL:ndL   ,0:Nmax)
  REAL   , INTENT(OUT  ) ::  ccR   (-ndR:ndR   ,0:Nmax)
  REAL   , INTENT(OUT  ) ::  ccL_LU(1:(3*ndL+1),0:Nmax)
  
  REAL                   ::  ccLT  (-ndL:ndL,-ndL:(Nmax+ndL))
  REAL                   ::  ccRT  (-ndR:ndR,-ndR:(Nmax+ndR))
  
  INTEGER                ::  pivot(0:Nmax)
  
  INTEGER                ::  i, ii
  REAL                   ::  dd
  
  
  ! TEST!!! Neumann-RB werden weiterhin separat und explizit gerechnet!
  !=====================================================================================================
  IF (ndL == 2 .AND. ndR == 3) THEN
     IF ((grid_type == 0) .OR. (grid_type == 1)) THEN
        IF (abl == 0) THEN
           ! Filter ...
           DO i = SS, NN
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
           END IF
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        ELSE IF (abl == 1) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/1., 10., 20., 10., 1./)
              ccR(-ndR:ndR,i) = (/-1., -101., -425., 0., 425., 101., 1./) / 30.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        ELSE IF (abl == 2) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/43., 668., 1798., 668., 43./)
              ccR(-ndR:ndR,i) = (/79., 4671., 9585., -28670., 9585., 4671., 79./) / 9.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        END IF
     !**************************************************************************************************
     ELSE IF (grid_type == 2) THEN ! TEST!!! Sollte frueher oder spaeter mit grid_type == 3 verschmolzen werden (ndR verallgemeinern!)
        IF (abl == 0) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/5., 60., 126., 60., 5./)
              ccR(-ndR:ndR,i) = (/0., 1., 45., 210., 210., 45., 1./) / 2.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        ELSE IF (abl == 1) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/9675., 193700., 577058., 193700., 9675./)
              ccR(-ndR:ndR,i) = (/0., -69049., -2525875., -6834250., 6834250., 2525875., 69049./) / 15.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        END IF
     ELSE IF (grid_type == 3) THEN
        IF (abl == 0) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/5., 60., 126., 60., 5./)
              ccR(-ndR:ndR,i) = (/1., 45., 210., 210., 45., 1., 0./) / 2.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        ELSE IF (abl == 1) THEN
           DO i = SS, NN
              ccL(-ndL:ndL,i) = (/9675., 193700., 577058., 193700., 9675./)
              ccR(-ndR:ndR,i) = (/-69049., -2525875., -6834250., 6834250., 2525875., 69049., 0./) / 15.
           END DO
           !--------------------------------------------------------------------------------------------
           IF (BCL >= 1) THEN
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
           END IF
           !--------------------------------------------------------------------------------------------
           IF (BCU >= 1) THEN
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
           END IF
        END IF
     END IF
  !=====================================================================================================
  ELSE IF (ndL == 3 .AND. ndR == 3) THEN
     DO i = SS, NN
        ccL(-ndL:ndL,i) = (/ 1. ,  36., 225. ,400.,225. ,36. ,1. /) / 400.
        ccR(-ndR:ndR,i) = (/-49.,-924.,-2625.,0.  ,2625.,924.,49./) / 4000.
     END DO
     !--------------------------------------------------------------------------------------------------
     IF (BCL >= 1) THEN
        ccL(-ndL:ndL,SS+0) = (/0.,0.,0. ,1. ,2. ,0.,0./)
        ccL(-ndL:ndL,SS+1) = (/0.,0.,1. ,4. ,1. ,0.,0./) / 4.
        ccL(-ndL:ndL,SS+2) = (/0.,1.,16.,36.,16.,1.,0./) / 36.
        
        ccR(-ndR:ndR,SS+0) = (/0., 0. , 0.  ,-5. ,4.  , 1. ,0./) / 2.
        ccR(-ndR:ndR,SS+1) = (/0., 0. ,-3.  , 0. ,3.  , 0. ,0./) / 4.
        ccR(-ndR:ndR,SS+2) = (/0.,-25.,-160., 0. ,160., 25.,0./) / 216.
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (BCU >= 1) THEN
        ccL(-ndL:ndL,NN-0) = (/0.,0.,2. ,1. ,0. ,0.,0./)
        ccL(-ndL:ndL,NN-1) = (/0.,0.,1. ,4. ,1. ,0.,0./) / 4.
        ccL(-ndL:ndL,NN-2) = (/0.,1.,16.,36.,16.,1.,0./) / 36.
        
        ccR(-ndR:ndR,NN-0) = (/ 0.,-1. ,-4.  ,5. ,0.  ,0. ,0./) / 2.
        ccR(-ndR:ndR,NN-1) = (/ 0., 0. ,-3.  ,0. ,3.  ,0. ,0./) / 4.
        ccR(-ndR:ndR,NN-2) = (/ 0.,-25.,-160.,0. ,160.,25.,0./) / 216.
     END IF
  !=====================================================================================================
  ELSE IF (ndL == 1 .AND. ndR == 1) THEN
     DO i = SS, NN
        ccL(-ndL:ndL,i) = (/ 1.,4.,1./) / 4.
        ccR(-ndR:ndR,i) = (/-3.,0.,3./) / 4.
     END DO
     !--------------------------------------------------------------------------------------------------
     IF (BCL >= 1) THEN
        ccL(-ndL:ndL,SS+0) = (/0., 1.,0./) ! TEST!!! Schlechte Wahl ...
        
        ccR(-ndR:ndR,SS+0) = (/0.,-1.,1./)
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (BCU >= 1) THEN
        ccL(-ndL:ndL,NN-0) = (/ 0.,1.,0./)
        
        ccR(-ndR:ndR,NN-0) = (/-1.,1.,0./)
     END IF
  !=====================================================================================================
  ELSE
     IF (rank == 0) WRITE(*,*) 'ERROR! No coefficients specified for this choice of ndL, ndR.'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  !=====================================================================================================
  
  
  ! TEST!!! nochmals richtig durchtesten fuer ndR /= 3 und ndL /= 2 ...
  !=====================================================================================================
  !=== Symmetrie =======================================================================================
  !=====================================================================================================
  IF (BCL == -2) THEN
     IF (grid_type == 0) THEN
        DO i = 1, ndL
           DO ii = 1, 1-i+ndL
              ccL(1-i+ii,i) = ccL(1-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           END DO
        END DO
        DO i = 1, ndR
           DO ii = 1, 1-i+ndR
              ccR(1-i+ii,i) = ccR(1-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 1) THEN
        DO i = 1, ndL
           DO ii = 1, 1-i+ndL
              ccL(0-i+ii,i) = ccL(0-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           END DO
        END DO
        DO i = 1, ndR
           DO ii = 1, 1-i+ndR
              ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 2) THEN
        DO i = 1, ndL
           DO ii = 1, 1-i+ndL
              ccL(0-i+ii,i) = ccL(0-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           END DO
        END DO
        DO i = 1, ndR
           DO ii = 1, 1-i+ndR
              ccR(1-i+ii,i) = ccR(1-i+ii,i) + ccR(1-i-ii,i)*sym ! TEST!!! Nur wegen aktuellem Shift bei grid_type == 2
              ccR(1-i-ii,i) = 0.
              !ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(0-i-ii,i)*sym
              !ccR(0-i-ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 3) THEN
        DO i = 1, ndL
           DO ii = 1, 1-i+ndL
              ccL(1-i+ii,i) = ccL(1-i+ii,i) + ccL(1-i-ii,i)*sym*(-1.)**abl
              ccL(1-i-ii,i) = 0.
           END DO
        END DO
        DO i = 1, ndR
           DO ii = 1, 1-i+ndR
              ccR(0-i+ii,i) = ccR(0-i+ii,i) + ccR(1-i-ii,i)*sym
              ccR(1-i-ii,i) = 0.
           END DO
        END DO
     END IF
  END IF
  !=====================================================================================================
  IF (BCU == -2) THEN
     IF (grid_type == 0) THEN
        DO i = Nmax-ndL+1, Nmax-0
           DO ii = 1, i+ndL-(Nmax-0)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i  +ii,i)*sym*(-1.)**abl
              ccL(Nmax-i  +ii,i) = 0.
           END DO
        END DO
        DO i = Nmax-ndR+1, Nmax-0
           DO ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i  +ii,i)*sym
              ccR(Nmax-i  +ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 1) THEN
        DO i = Nmax-ndL, Nmax-1
           DO ii = 1, i+ndL-(Nmax-1)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i-1+ii,i)*sym*(-1.)**abl
              ccL(Nmax-i-1+ii,i) = 0.
           END DO
        END DO
        DO i = Nmax-ndR+0, Nmax-1
           DO ii = 1, i+ndR-(Nmax-1)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i-1+ii,i)*sym
              ccR(Nmax-i-1+ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 2) THEN
        DO i = Nmax-ndL+0, Nmax-1
           DO ii = 1, i+ndL-(Nmax-1)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i-1+ii,i)*sym*(-1.)**abl
              ccL(Nmax-i-1+ii,i) = 0.
           END DO
        END DO
        DO i = Nmax-ndR+1, Nmax-0
           DO ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i  +ii,i)*sym! TEST!!! Nur wegen aktuellem Shift bei grid_type == 2
              ccR(Nmax-i  +ii,i) = 0.
        !      ccR(Nmax-i-1-ii,i) = ccR(Nmax-i-1-ii,i) + ccR(Nmax-i-1+ii,i)*sym
        !      ccR(Nmax-i-1+ii,i) = 0.
           END DO
        END DO
     END IF
     !--------------------------------------------------------------------------------------------------
     IF (grid_type == 3) THEN
        DO i = Nmax-ndL+1, Nmax-0
           DO ii = 1, i+ndL-(Nmax-0)
              ccL(Nmax-i  -ii,i) = ccL(Nmax-i  -ii,i) + ccL(Nmax-i  +ii,i)*sym*(-1.)**abl
              ccL(Nmax-i  +ii,i) = 0.
           END DO
        END DO
        DO i = Nmax-ndR+1, Nmax-0
           DO ii = 1, i+ndR-(Nmax-0)
              ccR(Nmax-i  -ii,i) = ccR(Nmax-i  -ii,i) + ccR(Nmax-i-1+ii,i)*sym
              ccR(Nmax-i-1+ii,i) = 0.
           END DO
        END DO
     END IF
  END IF
  !=====================================================================================================
  
  
  
  !=====================================================================================================
  !=== Transposition ===================================================================================
  !=====================================================================================================
  IF (transp_yes) THEN
     ccLT = 0.
     ccRT = 0.
     ccLT(-ndL:ndL,SS:NN) = ccL(-ndL:ndL,SS:NN)
     ccRT(-ndR:ndR,SS:NN) = ccR(-ndR:ndR,SS:NN)
     
     !*************************************
     ! TEST!!! Dieser Abschnitt ist nur erlaubt, wenn die Stencil-Koeffizienten im
     !         Innern (bzw. genauer gesagt: and den Blockgrenzen) identisch sind!
     !         (Ansonsten m�sste per MPI ausgetauscht werden ...)
     !         (gilt hier auf f�r Periodizit�t!)
     IF (BCL == 0 .OR. BCL == -1) THEN
        DO ii = -ndL, ndL
           ccLT(ii,-ndL:(SS-1)) = ccL(ii,SS)
        END DO
        DO ii = -ndR, ndR
           ccRT(ii,-ndR:(SS-1)) = ccR(ii,SS)
        END DO
     END IF
     IF (BCU == 0 .OR. BCU == -1) THEN
        DO ii = -ndL, ndL
           ccLT(ii,(NN+1):(NN+ndL)) = ccL(ii,NN)
        END DO
        DO ii = -ndR, ndR
           ccRT(ii,(NN+1):(NN+ndR)) = ccR(ii,NN)
        END DO
     END IF
     !*************************************
     ccL = 0.
     ccR = 0.
     DO i = SS, NN ! ok ...
        DO ii = -ndL, ndL
           ccL(ii,i) = ccLT(-ii,i+ii)
        END DO
     END DO
     !DO i = SS, NN ! TEST!!! ... Intervall SS:NN nicht ok! 
     DO i = 0, Nmax ! TEST!!! R�nder sollten noch auf Null gesetzt werden!
        DO ii = -ndR, ndR
           ccR(ii,i) = ccRT(-ii,i+ii)
        END DO
     END DO
     
  END IF
  !=====================================================================================================
  
  
  
  !=====================================================================================================
  !=== LU-Zerlegung ====================================================================================
  !=====================================================================================================
  ccL_LU = 0.
  
  !--- Transponierte Koeffizienten ---
  DO i = SS, NN
     DO ii = -ndL, ndL
        IF ((i+ii >= SS) .AND. (i+ii <= NN)) ccL_LU(2*ndL-ii+1,i+ii) = ccL(ii,i)
     END DO
  END DO
  
  !--- LU-Zerlegung ---
  IF (NB == 1 .AND. BCL /= -1) THEN
     CALL myDGBTRF(NN-SS+1      ,ndL,ccL_LU(1,SS    ),3*ndL+1,pivot(SS    ))
  ELSE
     CALL myDGBTRF(NN-SS+1-2*ndL,ndL,ccL_LU(1,SS+ndL),3*ndL+1,pivot(SS+ndL))
  END IF
  !=====================================================================================================
  
  
  END SUBROUTINE diff_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE myDGBTRF(N,ndL,AB,LDAB,pivot)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   )  ::  ndL, LDAB, N
  INTEGER, INTENT(OUT  )  ::  pivot( * )
  REAL   , INTENT(INOUT)  ::  AB( LDAB, * )
  INTEGER, PARAMETER      ::  NBMAX = 64
  INTEGER, PARAMETER      ::  LDWORK = NBMAX+1
  INTEGER                 ::  I, I2, I3, II, IP, J, J2, J3, JJ, JM, JP, JU, K2, KM, KV, NW
  REAL                    ::  TEMP
  REAL                    ::  WORK13( LDWORK, NBMAX ), WORK31( LDWORK, NBMAX )
  INTEGER                 ::  IDAMAX
  
  EXTERNAL IDAMAX
  EXTERNAL DCOPY, DGEMM, DGER, DLASWP, DSCAL, DSWAP, DTRSM
  
  
  KV = ndL + ndL
  
  DO J = ndL + 2, MIN( KV, N )
     DO I = KV - J + 2, ndL
        AB( I, J ) = 0.
     END DO
  END DO
  
  JU = 1
  
  DO J = 1, N
     I2 = MIN(ndL-1,N-J)
     I3 = MIN(1,N-J-ndL+1)
     
     IF(J+KV <= N) THEN
        DO I = 1, ndL
           AB(I,J+KV) = 0.
        END DO
     END IF
     
     KM = MIN(ndL,N-J)
     
     !--- Pivoting ---
     !JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
     JP = 1 ! TEST!!!
     
     pivot(J) = JP
     
     IF (AB(KV+JP,J) /= 0.) THEN
        JU = MAX( JU, MIN(J+ndL+JP-1,N))
        IF (JP /= 1) THEN
           IF (JP-1 < ndL) THEN
              CALL DSWAP( 1, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
           ELSE
              CALL DSWAP( 0, AB( KV+1, J ), LDAB-1,WORK31( JP-ndL, 1 ), LDWORK )
              CALL DSWAP( 1, AB( KV+1, J ), LDAB-1,AB( KV+JP, J ), LDAB-1 )
           END IF
        END IF
        CALL DSCAL( KM, 1. / AB( KV+1, J ), AB( KV+2, J ),1 )
        JM = MIN( JU, J )
        IF (JM > J) CALL DGER( KM, JM-J, -1., AB( KV+2, J ), 1, AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), LDAB-1 )
     END IF
     NW = MIN( 1, I3 )
     IF(NW > 0) CALL DCOPY( NW, AB( KV+ndL+1, J ), 1, WORK31( 1, 1 ), 1 )
     
     
     
     IF (J+1 <= N) THEN
        
        J2 = MIN( JU-J+1, KV ) - 1
        J3 = MAX( 0, JU-J-KV+1 )
        
        CALL DLASWP( J2, AB( KV, J+1 ), LDAB-1, 1, 1, pivot( J ), 1 )
        
        pivot( J ) = pivot( J ) + J - 1
        
        K2 = J + J2
        DO I = 1, J3
           JJ = K2 + I
           DO II = J + I - 1, J
              IP = pivot( II )
              IF (IP /= II) THEN
                 TEMP = AB( KV+1+II-JJ, JJ )
                 AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                 AB( KV+1+IP-JJ, JJ ) = TEMP
              END IF
           END DO
        END DO
        
        IF (J2 > 0) THEN
           CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J2, 1., AB( KV+1, J ), LDAB-1, AB( KV, J+1 ), LDAB-1 )
           IF (I2 > 0) CALL DGEMM( 'No transpose', 'No transpose', I2, J2, 1, -1., AB( KV+2, J ), LDAB-1, AB( KV, J+1 ), LDAB-1, 1., AB( KV+1, J+1 ), LDAB-1 )
           IF (I3 > 0) CALL DGEMM( 'No transpose', 'No transpose', I3, J2, 1, -1., WORK31, LDWORK, AB( KV, J+1 ), LDAB-1, 1., AB( KV+ndL, J+1 ), LDAB-1 )
        END IF
        
        IF (J3 > 0) THEN
           DO JJ = 1, J3
              DO II = JJ, 1
                 WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
              END DO
           END DO
           
           CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', 1, J3, 1., AB( KV+1, J ), LDAB-1, WORK13, LDWORK )
           IF( I2>0 ) CALL DGEMM( 'No transpose', 'No transpose', I2, J3, 1, -1., AB( KV+2, J ), LDAB-1, WORK13, LDWORK, 1., AB( 2, J+KV ), LDAB-1 )
           IF( I3>0 ) CALL DGEMM( 'No transpose', 'No transpose', I3, J3, 1, -1., WORK31, LDWORK, WORK13, LDWORK, 1., AB( 1+ndL, J+KV ), LDAB-1 )
           
           DO JJ = 1, J3
              DO II = JJ, 1
                 AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
              END DO
           END DO
        END IF
        
     ELSE
        pivot( J ) = pivot( J ) + J - 1
     END IF
     
     JP = pivot( J ) - J + 1
     IF (JP /= 1) THEN
        IF (JP-1 < ndL) THEN
           CALL DSWAP( 0, AB( KV+1, J ), LDAB-1, AB( KV+JP, J ), LDAB-1 )
        ELSE
           CALL DSWAP( 0, AB( KV+1, J ), LDAB-1, WORK31( JP-ndL, 1 ), LDWORK )
        END IF
     END IF
     NW = MIN( I3, 1 )
     IF( NW>0 ) CALL DCOPY( NW, WORK31( 1, 1 ), 1, AB( KV+ndL+1, J ), 1 )
     
  END DO
  
  
  END SUBROUTINE myDGBTRF
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE compact_parallel(SS,NN,Nmax,iB,NB,ndL,ccL,dimS,BCL,Schur,WW,COMM)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)   ::  SS, NN
  INTEGER, INTENT(IN)   ::  Nmax
  INTEGER, INTENT(IN)   ::  iB, NB
  INTEGER, INTENT(IN)   ::  ndL
  REAL   , INTENT(IN)   ::  ccL(-ndL:ndL,0:Nmax)
  INTEGER, INTENT(IN)   ::  dimS
  INTEGER, INTENT(IN)   ::  BCL
  REAL   , INTENT(OUT)  ::  Schur(1:dimS,1:dimS)
  REAL   , INTENT(OUT)  ::  WW(1:2*ndL,0:Nmax)
  INTEGER, INTENT(IN)   ::  COMM
  
  REAL                  ::  VV(1:2*ndL,0:Nmax)
  
  REAL   , ALLOCATABLE  ::  VV1(:,:), VV2(:,:)
  REAL   , ALLOCATABLE  ::  WW1(:,:), WW2(:,:)
  REAL                  ::  Schur_global(1:dimS,1:dimS)
  
  REAL                  ::  ccL_ext(0:3*ndL,SS:NN)
  INTEGER               ::  pivot1(0:Nmax), pivot2(1:dimS)
  INTEGER               ::  dimG, dimA1(1:NB), dimA2(1:NB), dimA, i0
  INTEGER               ::  i, j, ii
  INTEGER               ::  info
  
  
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
  
  DO ii = 1, ndL
     !--- Zentraler Block (global) ---
     Schur(i+ii    ,(i+1    ):(i+1*ndL)) = ccL((1-ii):(ndL-ii),SS+ii-1  )
     Schur(i+ii+ndL,(i+1+ndL):(i+2*ndL)) = ccL((1-ii):(ndL-ii),NN+ii-ndL)
     
     !--- restliche Bloecke (global) ---
     IF (iB /= 1                 ) Schur(i+ii,(i   -ndL+ii):i   ) = ccL(-ndL:-ii,SS+ii-1)
     IF (iB == 1  .AND. BCL == -1) Schur(i+ii,(dimS-ndL+ii):dimS) = ccL(-ndL:-ii,SS+ii-1)
     
     IF (iB /= NB                ) Schur(i+ii+ndL,(i+2*ndL+1):(i+2*ndL+ii)) = ccL((ndL-ii+1):ndL,NN+ii-ndL)
     IF (iB == NB .AND. BCL == -1) Schur(i+ii+ndL,         1 :         ii ) = ccL((ndL-ii+1):ndL,NN+ii-ndL)
     
     !--- restliche Bloecke (lokal) ---
     VV(    ii : ndL    ,(SS+ndL +ii-1)) = ccL(     -ndL :-ii,SS+ii-1+ndL)
     VV((1+ndL):(ii+ndL),(SS+dimG+ii-1)) = ccL((1+ndL-ii):ndL,NN+ii-2*ndL)
     
     WW(ii    ,(SS+ndL      ):(SS+ndL-1+ii  )) = ccL((ndL-ii+1):ndL,SS+ii-1  )
     WW(ii+ndL,(SS+dimG+ii-1):(SS+ndL-1+dimG)) = ccL(  -ndL    :-ii,NN+ii-ndL)
  END DO
  
  
  !--- Bandmatrix erweitern fuer LAPACK ---
  ccL_ext = 0.
  DO i = SS, NN
     ccL_ext(ndL:3*ndL,i) = ccL(-ndL:ndL,i)
  END DO
  
  
  !--- Bandloeser fuer WW ---
  CALL DGBTRF(dimG,dimG,ndL,ndL,ccL_ext(0,SS+ndL),3*ndL+1,pivot1,info)
  DO i = 1, 2*ndL
     CALL DGBTRS('N',dimG,ndL,ndL,1,ccL_ext(0,SS+ndL),3*ndL+1,pivot1,WW(i,(SS+ndL):(SS+ndL+dimG-1):1),dimG,info)
  END DO
  
  
  !--- Globale Dimensionen und Positionen bestimmen ---
  dimA1     = 0
  dimA1(iB) = dimG
  
  CALL MPI_ALLREDUCE(dimA1,dimA2,NB,MPI_INTEGER,MPI_SUM,COMM,merror)
  
  dimA = 0
  DO i = 1, NB
     dimA = dimA + dimA2(i)
  END DO
  
  i0   = 0
  DO i = 1, iB-1
     i0   = i0   + dimA2(i)
  END DO
  
  ALLOCATE(VV1(1:dimS,1:dimA))
  ALLOCATE(VV2(1:dimS,1:dimA))
  ALLOCATE(WW1(1:dimS,1:dimA))
  ALLOCATE(WW2(1:dimS,1:dimA))
  
  VV1 = 0.
  VV2 = 0.
  WW1 = 0.
  WW2 = 0.
  
  
  !--- Umspeichern auf lokalen Array ---
  DO i = 1, dimG
     VV1((1+2*ndL*(iB-1)):2*ndL*iB,i+i0) = VV(1:2*ndL,(SS+ndL+i-1))
     WW1((1+2*ndL*(iB-1)):2*ndL*iB,i+i0) = WW(1:2*ndL,(SS+ndL+i-1))
  END DO
  
  
  !--- Matrizen global zusammenbauen ---
  CALL MPI_ALLREDUCE(VV1,VV2,dimS*dimA,MPI_REAL8,MPI_SUM,COMM,merror)
  CALL MPI_ALLREDUCE(WW1,WW2,dimS*dimA,MPI_REAL8,MPI_SUM,COMM,merror)
  CALL MPI_ALLREDUCE(Schur,Schur_global,dimS*dimS,MPI_REAL8,MPI_SUM,COMM,merror)
  
  
  !--- Schur-Komplement bilden ---
  DO j = 1, dimS
     DO i = 1, dimS
        DO ii = 1, dimA
           Schur_global(i,j) = Schur_global(i,j) - WW2(i,ii)*VV2(j,ii)
        END DO
     END DO
  END DO
  
  DEALLOCATE(VV1)
  DEALLOCATE(VV2)
  DEALLOCATE(WW1)
  DEALLOCATE(WW2)
  
  
  !--- Einheitsmatrix als RHS ---
  Schur = 0.
  DO i = 1, dimS
     Schur(i,i) = 1.
  END DO
  
  
  !--- Inverse des Schur-Komplements ---
  ! TEST!! Hier gibt es in Lapack einen Floating Point exception ..
  !CALL DGESV(dimS,dimS,Schur_global,dimS,pivot2,Schur,dimS,info)
  CALL matrix_invert(dimS,Schur_global,Schur)
  
  
  END SUBROUTINE compact_parallel
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_coeffs_compact
  
  IMPLICIT NONE
  
  INTEGER                ::  h, i, j, k
  INTEGER                ::  abl, grid
  REAL                   ::  dd1, max_diff, max_diff_global
  REAL                   ::  omega, phase, xref
  REAL                   ::  pi
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  omega = 4.
  phase = 1.
  
  !===========================================================================================================
  DO h = 1, 8
     
     IF (h == 1) THEN
        abl  = 0
        grid = 2
     ELSE IF (h == 2) THEN
        abl  = 0
        grid = 3
     ELSE IF (h == 3) THEN
        abl  = 1
        grid = 0
     ELSE IF (h == 4) THEN
        abl  = 1
        grid = 1
     ELSE IF (h == 5) THEN
        abl  = 1
        grid = 2
     ELSE IF (h == 6) THEN
        abl  = 1
        grid = 3
     ELSE IF (h == 7) THEN
        abl  = 2
        grid = 0
     ELSE IF (h == 8) THEN
        abl  = 2
        grid = 1
     END IF
     
     
     xref = phase*0.25*L1/omega
     
     IF (BC_1L_global == -2 .OR. BC_1U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L1
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (grid == 0 .OR. grid == 2) pp(i,j,k) =  COS(dd1*(x1p(i)-xref))
              
              IF (grid == 0 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x1p(i)-xref))
              IF (grid == 0 .AND. abl == 2) rh(i,j,k) = -COS(dd1*(x1p(i)-xref))
              
              IF (grid == 3 .AND. abl == 0) rh(i,j,k) =  SIN(dd1*(x1p(i)-xref))
              IF (grid == 3 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x1p(i)-xref))
           END DO
        END DO
     END DO
     DO k = S31B, N31B
        DO j = S21B, N21B
           DO i = S11B, N11B
              IF (grid == 1 .OR. grid == 3) pp(i,j,k) =  SIN(dd1*(x1u(i)-xref))
              
              IF (grid == 1 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x1u(i)-xref))
              IF (grid == 1 .AND. abl == 2) rh(i,j,k) = -SIN(dd1*(x1u(i)-xref))
              
              IF (grid == 2 .AND. abl == 0) rh(i,j,k) =  COS(dd1*(x1u(i)-xref))
              IF (grid == 2 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x1u(i)-xref))
           END DO
        END DO
     END DO
     
     pp = pp * dd1**(-abl)
     
     IF (grid == 0 .OR. grid == 2) CALL exchange(1,0,pp)
     IF (grid == 1 .OR. grid == 3) CALL exchange(1,1,pp)
     
     IF (abl == 0 .AND. grid == 2) CALL apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cIpuCL,cIpuCL_LU,cIpuCR,WIpu,SIpu,pp,Ap)
     IF (abl == 0 .AND. grid == 3) CALL apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cIupCL,cIupCL_LU,cIupCR,WIup,SIup,pp,Ap)
     
     IF (abl == 1 .AND. grid == 0) CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ap)
     IF (abl == 1 .AND. grid == 1) CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ap)
     IF (abl == 1 .AND. grid == 2) CALL apply_compact(1,0,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cGp1CL,cGp1CL_LU,cGp1CR,WGp1,SGp1,pp,Ap)
     IF (abl == 1 .AND. grid == 3) CALL apply_compact(1,1,S1p,S2p,S3p,N1p,N2p,N3p,N1,ndL,ndR,dimS1,cDu1CL,cDu1CL_LU,cDu1CR,WDu1,SDu1,pp,Ap)
     
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp1CL ,cp1CL_LU ,cp1CR ,Wp1 ,Sp1 ,pp,Ap)
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(1,0,S12,S22,S32,N12,N22,N32,N1,ndL,ndR,dimS1,cp11CL,cp11CL_LU,cp11CR,Wp11,Sp11,pp,rr)
     
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu1CL ,cu1CL_LU ,cu1CR ,Wu1 ,Su1 ,pp,Ap)
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(1,1,S11,S21,S31,N11,N21,N31,N1,ndL,ndR,dimS1,cu11CL,cu11CL_LU,cu11CR,Wu11,Su11,pp,rr)
     
     
     max_diff = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (abl == 1 .AND. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx1DM(i)
              
              IF (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (abl == 1 .AND. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx1pM(i)
              IF (abl == 2 .AND. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx1pM(i)**2 + Ap(i,j,k)*ddx1pM(i)
              
              IF (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
           
        END DO
     END DO
     
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (abl == 1 .AND. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx1uM(i)
              IF (abl == 1 .AND. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx1GM(i)
              IF (abl == 2 .AND. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx1uM(i)**2 + Ap(i,j,k)*ddx1uM(i)
              
              IF (grid == 1 .OR. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     IF (rank == 0) WRITE(*,*) max_diff_global
     
     
     IF (1 == 1 .AND. h == 8) THEN
     IF (grid == 0) THEN
        j = N2p
        k = N3p
        DO i = S12, N12
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(55,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 2 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(56,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 3 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(57,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 4 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(58,*) x1p(i), Ap(i,j,k)
        END DO
     ELSE IF (grid == 3) THEN
        j = N2p
        k = N3p
        DO i = S1p, N1p
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(55,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 2 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(56,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 3 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(57,*) x1p(i), Ap(i,j,k)
           IF (iB(1,1) == 4 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(58,*) x1p(i), Ap(i,j,k)
        END DO
     ELSE
        j = N21
        k = N31
        DO i = S11, N11
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(55,*) x1u(i), Ap(i,j,k)
           IF (iB(1,1) == 2 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(56,*) x1u(i), Ap(i,j,k)
           IF (iB(1,1) == 3 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(57,*) x1u(i), Ap(i,j,k)
           IF (iB(1,1) == 4 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(58,*) x1u(i), Ap(i,j,k)
        END DO
     END IF
     END IF
  END DO
  
  !===========================================================================================================
  IF (rank == 0) WRITE(*,*)
  
  DO h = 1, 8
     
     IF (h == 1) THEN
        abl  = 0
        grid = 2
     ELSE IF (h == 2) THEN
        abl  = 0
        grid = 3
     ELSE IF (h == 3) THEN
        abl  = 1
        grid = 0
     ELSE IF (h == 4) THEN
        abl  = 1
        grid = 1
     ELSE IF (h == 5) THEN
        abl  = 1
        grid = 2
     ELSE IF (h == 6) THEN
        abl  = 1
        grid = 3
     ELSE IF (h == 7) THEN
        abl  = 2
        grid = 0
     ELSE IF (h == 8) THEN
        abl  = 2
        grid = 1
     END IF
     
     
     xref = phase*0.25*L2/omega
     
     IF (BC_2L_global == -2 .OR. BC_2U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L2
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (grid == 0 .OR. grid == 2) pp(i,j,k) =  COS(dd1*(x2p(j)-xref))
              
              IF (grid == 0 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x2p(j)-xref))
              IF (grid == 0 .AND. abl == 2) rh(i,j,k) = -COS(dd1*(x2p(j)-xref))
              
              IF (grid == 3 .AND. abl == 0) rh(i,j,k) =  SIN(dd1*(x2p(j)-xref))
              IF (grid == 3 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x2p(j)-xref))
           END DO
        END DO
     END DO
     DO k = S32B, N32B
        DO j = S22B, N22B
           DO i = S12B, N12B
              IF (grid == 1 .OR. grid == 3) pp(i,j,k) =  SIN(dd1*(x2v(j)-xref))
              
              IF (grid == 1 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x2v(j)-xref))
              IF (grid == 1 .AND. abl == 2) rh(i,j,k) = -SIN(dd1*(x2v(j)-xref))
              
              IF (grid == 2 .AND. abl == 0) rh(i,j,k) =  COS(dd1*(x2v(j)-xref))
              IF (grid == 2 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x2v(j)-xref))
           END DO
        END DO
     END DO
     
     pp = pp * dd1**(-abl)
     
     IF (grid == 0 .OR. grid == 2) CALL exchange(2,0,pp)
     IF (grid == 1 .OR. grid == 3) CALL exchange(2,2,pp)
     
     IF (abl == 0 .AND. grid == 2) CALL apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cIpvCL,cIpvCL_LU,cIpvCR,WIpv,SIpv,pp,Ap)
     IF (abl == 0 .AND. grid == 3) CALL apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cIvpCL,cIvpCL_LU,cIvpCR,WIvp,SIvp,pp,Ap)
     
     IF (abl == 1 .AND. grid == 0) CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ap)
     IF (abl == 1 .AND. grid == 1) CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ap)
     IF (abl == 1 .AND. grid == 2) CALL apply_compact(2,0,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cGp2CL,cGp2CL_LU,cGp2CR,WGp2,SGp2,pp,Ap)
     IF (abl == 1 .AND. grid == 3) CALL apply_compact(2,2,S1p,S2p,S3p,N1p,N2p,N3p,N2,ndL,ndR,dimS2,cDv2CL,cDv2CL_LU,cDv2CR,WDv2,SDv2,pp,Ap)
     
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp2CL ,cp2CL_LU ,cp2CR ,Wp2 ,Sp2 ,pp,Ap)
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(2,0,S11,S21,S31,N11,N21,N31,N2,ndL,ndR,dimS2,cp22CL,cp22CL_LU,cp22CR,Wp22,Sp22,pp,rr)
     
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv2CL ,cv2CL_LU ,cv2CR ,Wv2 ,Sv2 ,pp,Ap)
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(2,2,S12,S22,S32,N12,N22,N32,N2,ndL,ndR,dimS2,cv22CL,cv22CL_LU,cv22CR,Wv22,Sv22,pp,rr)
     
     
     max_diff = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (abl == 1 .AND. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx2DM(j)
              
              IF (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (abl == 1 .AND. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx2pM(j)
              IF (abl == 2 .AND. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx2pM(j)**2 + Ap(i,j,k)*ddx2pM(j)
              
              IF (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     DO k = S32, N32
        DO j = S22, N22
           DO i = S12, N12
              IF (abl == 1 .AND. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx2vM(j)
              IF (abl == 1 .AND. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx2GM(j)
              IF (abl == 2 .AND. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx2vM(j)**2 + Ap(i,j,k)*ddx2vM(j)
              
              IF (grid == 1 .OR. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     IF (rank == 0) WRITE(*,*) max_diff_global
     
     
     IF (1 == 1 .AND. h == 8) THEN
     IF (grid == 0) THEN
        i = N1p
        k = N3p
        DO j = S21, N21
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(65,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 2 .AND. iB(3,1) == 1) WRITE(66,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 3 .AND. iB(3,1) == 1) WRITE(67,*) x2p(j), Ap(i,j,k), rh(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 4 .AND. iB(3,1) == 1) WRITE(68,*) x2p(j), Ap(i,j,k), rh(i,j,k)
        END DO
     ELSE IF (grid == 3) THEN
        i = N1p
        k = N3p
        DO j = S2p, N2p
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(65,*) x2p(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 2 .AND. iB(3,1) == 1) WRITE(66,*) x2p(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 3 .AND. iB(3,1) == 1) WRITE(67,*) x2p(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 4 .AND. iB(3,1) == 1) WRITE(68,*) x2p(j), Ap(i,j,k)
        END DO
     ELSE
        i = N12
        k = N32
        DO j = S22, N22
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(65,*) x2v(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 2 .AND. iB(3,1) == 1) WRITE(66,*) x2v(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 3 .AND. iB(3,1) == 1) WRITE(67,*) x2v(j), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 4 .AND. iB(3,1) == 1) WRITE(68,*) x2v(j), Ap(i,j,k)
        END DO
     END IF
     END IF
  END DO
  
  !===========================================================================================================
  IF (dimens == 3) THEN
  
  IF (rank == 0) WRITE(*,*)
  
  DO h = 1, 8
     
     IF (h == 1) THEN
        abl  = 0
        grid = 2
     ELSE IF (h == 2) THEN
        abl  = 0
        grid = 3
     ELSE IF (h == 3) THEN
        abl  = 1
        grid = 0
     ELSE IF (h == 4) THEN
        abl  = 1
        grid = 1
     ELSE IF (h == 5) THEN
        abl  = 1
        grid = 2
     ELSE IF (h == 6) THEN
        abl  = 1
        grid = 3
     ELSE IF (h == 7) THEN
        abl  = 2
        grid = 0
     ELSE IF (h == 8) THEN
        abl  = 2
        grid = 1
     END IF
     
     
     xref = phase*0.25*L3/omega
     
     IF (BC_3L_global == -2 .OR. BC_3U_global == -2) xref = 0.
     
     dd1 = omega*2.*pi/L3
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (grid == 0 .OR. grid == 2) pp(i,j,k) =  COS(dd1*(x3p(k)-xref))
              
              IF (grid == 0 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x3p(k)-xref))
              IF (grid == 0 .AND. abl == 2) rh(i,j,k) = -COS(dd1*(x3p(k)-xref))
              
              IF (grid == 3 .AND. abl == 0) rh(i,j,k) =  SIN(dd1*(x3p(k)-xref))
              IF (grid == 3 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x3p(k)-xref))
           END DO
        END DO
     END DO
     DO k = S33B, N33B
        DO j = S23B, N23B
           DO i = S13B, N13B
              IF (grid == 1 .OR. grid == 3) pp(i,j,k) =  SIN(dd1*(x3w(k)-xref))
              
              IF (grid == 1 .AND. abl == 1) rh(i,j,k) =  COS(dd1*(x3w(k)-xref))
              IF (grid == 1 .AND. abl == 2) rh(i,j,k) = -SIN(dd1*(x3w(k)-xref))
              
              IF (grid == 2 .AND. abl == 0) rh(i,j,k) =  COS(dd1*(x3w(k)-xref))
              IF (grid == 2 .AND. abl == 1) rh(i,j,k) = -SIN(dd1*(x3w(k)-xref))
           END DO
        END DO
     END DO
     
     pp = pp * dd1**(-abl)
     
     IF (grid == 0 .OR. grid == 2) CALL exchange(3,0,pp)
     IF (grid == 1 .OR. grid == 3) CALL exchange(3,3,pp)
     
     IF (abl == 0 .AND. grid == 2) CALL apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cIpwCL,cIpwCL_LU,cIpwCR,WIpw,SIpw,pp,Ap)
     IF (abl == 0 .AND. grid == 3) CALL apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cIwpCL,cIwpCL_LU,cIwpCR,WIwp,SIwp,pp,Ap)
     
     IF (abl == 1 .AND. grid == 0) CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ap)
     IF (abl == 1 .AND. grid == 1) CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ap)
     IF (abl == 1 .AND. grid == 2) CALL apply_compact(3,0,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cGp3CL,cGp3CL_LU,cGp3CR,WGp3,SGp3,pp,Ap)
     IF (abl == 1 .AND. grid == 3) CALL apply_compact(3,3,S1p,S2p,S3p,N1p,N2p,N3p,N3,ndL,ndR,dimS3,cDw3CL,cDw3CL_LU,cDw3CR,WDw3,SDw3,pp,Ap)
     
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp3CL ,cp3CL_LU ,cp3CR ,Wp3 ,Sp3 ,pp,Ap)
     IF (abl == 2 .AND. grid == 0) CALL apply_compact(3,0,S11,S21,S31,N11,N21,N31,N3,ndL,ndR,dimS3,cp33CL,cp33CL_LU,cp33CR,Wp33,Sp33,pp,rr)
     
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw3CL ,cw3CL_LU ,cw3CR ,Ww3 ,Sw3 ,pp,Ap)
     IF (abl == 2 .AND. grid == 1) CALL apply_compact(3,3,S13,S23,S33,N13,N23,N33,N3,ndL,ndR,dimS3,cw33CL,cw33CL_LU,cw33CR,Ww33,Sw33,pp,rr)
     
     
     max_diff = 0.
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              IF (abl == 1 .AND. grid == 3) Ap(i,j,k) = Ap(i,j,k)*dx3DM(k)
              
              IF (grid == 3) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     DO k = S31, N31
        DO j = S21, N21
           DO i = S11, N11
              IF (abl == 1 .AND. grid == 0) Ap(i,j,k) = Ap(i,j,k)*dx3pM(k)
              IF (abl == 2 .AND. grid == 0) Ap(i,j,k) = rr(i,j,k)*dx3pM(k)**2 + Ap(i,j,k)*ddx3pM(k)
              
              IF (grid == 0) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     DO k = S33, N33
        DO j = S23, N23
           DO i = S13, N13
              IF (abl == 1 .AND. grid == 1) Ap(i,j,k) = Ap(i,j,k)*dx3wM(k)
              IF (abl == 1 .AND. grid == 2) Ap(i,j,k) = Ap(i,j,k)*dx3GM(k)
              IF (abl == 2 .AND. grid == 1) Ap(i,j,k) = rr(i,j,k)*dx3wM(k)**2 + Ap(i,j,k)*ddx3wM(k)
              
              IF (grid == 1 .OR. grid == 2) max_diff = MAX(max_diff,ABS(Ap(i,j,k)-rh(i,j,k)))
           END DO
        END DO
     END DO
     
     !CALL MPI_ALLREDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,COMM_CART,merror)
     CALL MPI_REDUCE(max_diff,max_diff_global,1,MPI_REAL8,MPI_MAX,0,COMM_CART,merror) ! TEST!!!
     
     IF (rank == 0) WRITE(*,*) max_diff_global
     
     
     IF (1 == 1 .AND. h == 3) THEN
     IF (grid == 0) THEN
        i = N1p
        j = N2p
        DO k = S31, N31
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(75,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 2) WRITE(76,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 3) WRITE(77,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 4) WRITE(78,*) x3p(k), Ap(i,j,k)
        END DO
     ELSE IF (grid == 3) THEN
        i = N1p
        j = N2p
        DO k = S3p, N3p
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(75,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 2) WRITE(76,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 3) WRITE(77,*) x3p(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 4) WRITE(78,*) x3p(k), Ap(i,j,k)
        END DO
     ELSE
        i = N13
        j = N23
        DO k = S33, N33
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 1) WRITE(75,*) x3w(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 2) WRITE(76,*) x3w(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 3) WRITE(77,*) x3w(k), Ap(i,j,k)
           IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. iB(3,1) == 4) WRITE(78,*) x3w(k), Ap(i,j,k)
        END DO
     END IF
     END IF
  END DO
  
  END IF
  !===========================================================================================================
  CALL MPI_FINALIZE(merror)
  STOP
  
  
  END SUBROUTINE test_coeffs_compact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE FD_coeffs
  
  IMPLICIT NONE
  
  INTEGER               ::  m
  
  
  !===========================================================================================================
  CALL diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1D)
  CALL diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,n1L,n1U,cNp1U)
  CALL diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1D)
  CALL diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,n1L,n1U,cNu1U)
  
  CALL diff_coeffs(1,0, 0,.FALSE.    ,dim_ncb1c,ncb1f,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cFp1 )
  CALL diff_coeffs(1,0, 0,.FALSE.    ,dim_ncb1c,ncb1f,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cFu1 )
  
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp1  )
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu1  )
  
  CALL diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BC_1L,BC_1U,0,N1,b1L,b1U,cp11 )
  CALL diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1u              ,x1u              ,BC_1L,BC_1U,1,N1,b1L,b1U,cu11 )
  
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cGp1 )
  CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cDu1 )
  
  CALL diff_coeffs(1,0, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BC_1L,BC_1U,2,N1,g1L,g1U,cIpu )
  CALL diff_coeffs(1,0, 0,mapping_yes,dim_ncb1d,ncb1d,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),BC_1L,BC_1U,3,N1,d1L,d1U,cIup )
  !-----------------------------------------------------------------------------------------------------------
  CALL diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2D)
  CALL diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,n2L,n2U,cNp2U)
  CALL diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2D)
  CALL diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,n2L,n2U,cNv2U)
  
  CALL diff_coeffs(2,0, 0,.FALSE.    ,dim_ncb2c,ncb2f,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cFp2 )
  CALL diff_coeffs(2,0, 0,.FALSE.    ,dim_ncb2c,ncb2f,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cFv2 )
  
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp2  )
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv2  )
  
  CALL diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BC_2L,BC_2U,0,N2,b2L,b2U,cp22 )
  CALL diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2v              ,x2v              ,BC_2L,BC_2U,1,N2,b2L,b2U,cv22 )
  
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cGp2 )
  CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cDv2 )
  
  CALL diff_coeffs(2,0, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BC_2L,BC_2U,2,N2,g2L,g2U,cIpv )
  CALL diff_coeffs(2,0, 0,mapping_yes,dim_ncb2d,ncb2d,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),BC_2L,BC_2U,3,N2,d2L,d2U,cIvp )
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  CALL diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3D)
  CALL diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,n3L,n3U,cNp3U)
  CALL diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3D)
  CALL diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,n3L,n3U,cNw3U)
  
  CALL diff_coeffs(3,0, 0,.FALSE.    ,dim_ncb3c,ncb3f,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cFp3 )
  CALL diff_coeffs(3,0, 0,.FALSE.    ,dim_ncb3c,ncb3f,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cFw3 )
  
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp3  )
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw3  )
  
  CALL diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BC_3L,BC_3U,0,N3,b3L,b3U,cp33 )
  CALL diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3w              ,x3w              ,BC_3L,BC_3U,1,N3,b3L,b3U,cw33 )
  
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cGp3 )
  CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cDw3 )
  
  CALL diff_coeffs(3,0, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BC_3L,BC_3U,2,N3,g3L,g3U,cIpw )
  CALL diff_coeffs(3,0, 0,mapping_yes,dim_ncb3d,ncb3d,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),BC_3L,BC_3U,3,N3,d3L,d3U,cIwp )
  END IF
  !===========================================================================================================
  
  
  DO m = 1, n_conc
     !========================================================================================================
     CALL diff_coeffs(1,1,-1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,n1L,n1U,cNc1D(n1L,0,m))
     CALL diff_coeffs(1,1, 1,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,n1L,n1U,cNc1U(n1L,0,m))
     CALL diff_coeffs(1,0, 0,.FALSE.    ,dim_ncb1c,ncb1f,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cFc1 (b1L,0,m))
     CALL diff_coeffs(1,1, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cc1  (b1L,0,m))
     CALL diff_coeffs(1,2, 0,mapping_yes,dim_ncb1c,ncb1c,x1p              ,x1p              ,BCc_1L(m),BCc_1U(m),0,N1,b1L,b1U,cc11 (b1L,0,m))
     CALL diff_coeffs(1,0, 0,mapping_yes,dim_ncb1g,ncb1g,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),BCc_1L(m),BCc_1U(m),2,N1,g1L,g1U,cIcu (g1L,0,m))
     !--------------------------------------------------------------------------------------------------------
     CALL diff_coeffs(2,1,-1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,n2L,n2U,cNc2D(n2L,0,m))
     CALL diff_coeffs(2,1, 1,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,n2L,n2U,cNc2U(n2L,0,m))
     CALL diff_coeffs(2,0, 0,.FALSE.    ,dim_ncb2c,ncb2f,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cFc2 (b2L,0,m))
     CALL diff_coeffs(2,1, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cc2  (b2L,0,m))
     CALL diff_coeffs(2,2, 0,mapping_yes,dim_ncb2c,ncb2c,x2p              ,x2p              ,BCc_2L(m),BCc_2U(m),0,N2,b2L,b2U,cc22 (b2L,0,m))
     CALL diff_coeffs(2,0, 0,mapping_yes,dim_ncb2g,ncb2g,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),BCc_2L(m),BCc_2U(m),2,N2,g2L,g2U,cIcv (g2L,0,m))
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     CALL diff_coeffs(3,1,-1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,n3L,n3U,cNc3D(n3L,0,m))
     CALL diff_coeffs(3,1, 1,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,n3L,n3U,cNc3U(n3L,0,m))
     CALL diff_coeffs(3,0, 0,.FALSE.    ,dim_ncb3c,ncb3f,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cFc3 (b3L,0,m))
     CALL diff_coeffs(3,1, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cc3  (b3L,0,m))
     CALL diff_coeffs(3,2, 0,mapping_yes,dim_ncb3c,ncb3c,x3p              ,x3p              ,BCc_3L(m),BCc_3U(m),0,N3,b3L,b3U,cc33 (b3L,0,m))
     CALL diff_coeffs(3,0, 0,mapping_yes,dim_ncb3g,ncb3g,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),BCc_3L(m),BCc_3U(m),2,N3,g3L,g3U,cIcw (g3L,0,m))
     END IF
     !========================================================================================================
  END DO
  
  !===========================================================================================================
                   CALL diff_coeffs(1,-1,0,mapping_yes,dim_ncb1c,ncb1c,x1p,x1p,BC_1L,BC_1U,0,N1,b1L,b1U,cInt1)
                   CALL diff_coeffs(2,-1,0,mapping_yes,dim_ncb2c,ncb2c,x2p,x2p,BC_2L,BC_2U,0,N2,b2L,b2U,cInt2)
  IF (dimens == 3) CALL diff_coeffs(3,-1,0,mapping_yes,dim_ncb3c,ncb3c,x3p,x3p,BC_3L,BC_3U,0,N3,b3L,b3U,cInt3)
  !===========================================================================================================
  
  
  END SUBROUTINE FD_coeffs
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_coeffs
  
  IMPLICIT NONE
  
  
  ! ACHTUNG!!! Upwind-Differenzen noch nicht getestet!!!
  !===========================================================================================================
  IF (iB(2,1) == 1 .AND. iB(3,1) == 1) THEN
     CALL test_diff(1,0,x1p              ,x1p              ,N1,b1L,b1U,cFp1 ,BC_1L,0,'cFp1' )
     CALL test_diff(1,0,x1u              ,x1u              ,N1,b1L,b1U,cFu1 ,BC_1L,1,'cFu1' )
     CALL test_diff(1,1,x1u              ,x1u              ,N1,b1L,b1U,cu1  ,BC_1L,1,'cu1'  )
     CALL test_diff(1,1,x1p              ,x1p              ,N1,b1L,b1U,cp1  ,BC_1L,0,'cp1'  )
     CALL test_diff(1,2,x1u              ,x1u              ,N1,b1L,b1U,cu11 ,BC_1L,1,'cu11' )
     CALL test_diff(1,2,x1p              ,x1p              ,N1,b1L,b1U,cp11 ,BC_1L,0,'cp11' )
     
     !CALL test_diff(1,1,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cGp1 ,BC_1L,2,'cGp1' ) ! TEST!!! Funktioniert nicht mehr ...
     !CALL test_diff(1,1,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cDu1 ,BC_1L,3,'cDu1' )
     !CALL test_diff(1,0,x1p(g1L:(N1+g1U)),x1u(g1L:(N1+g1U)),N1,g1L,g1U,cIpu ,BC_1L,2,'cIpu' )
     !CALL test_diff(1,0,x1u(d1L:(N1+d1U)),x1p(d1L:(N1+d1U)),N1,d1L,d1U,cIup ,BC_1L,3,'cIup' )
     
     CALL test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1D,BC_1L,1,'cNu1D')
     CALL test_diff(1,1,x1u(n1L:(N1+n1U)),x1u(n1L:(N1+n1U)),N1,n1L,n1U,cNu1U,BC_1L,1,'cNu1U')
     CALL test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1D,BC_1L,0,'cNp1D')
     CALL test_diff(1,1,x1p(n1L:(N1+n1U)),x1p(n1L:(N1+n1U)),N1,n1L,n1U,cNp1U,BC_1L,0,'cNp1U')
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (iB(1,1) == 1 .AND. iB(3,1) == 1) THEN
     CALL test_diff(2,0,x2p              ,x2p              ,N2,b2L,b2U,cFp2 ,BC_2L,0,'cFp2' )
     CALL test_diff(2,0,x2v              ,x2v              ,N2,b2L,b2U,cFv2 ,BC_2L,1,'cFv2' )
     CALL test_diff(2,1,x2p              ,x2p              ,N2,b2L,b2U,cp2  ,BC_2L,0,'cp2'  )
     CALL test_diff(2,1,x2v              ,x2v              ,N2,b2L,b2U,cv2  ,BC_2L,1,'cv2'  )
     CALL test_diff(2,2,x2p              ,x2p              ,N2,b2L,b2U,cp22 ,BC_2L,0,'cp22' )
     CALL test_diff(2,2,x2v              ,x2v              ,N2,b2L,b2U,cv22 ,BC_2L,1,'cv22' )
     
     !CALL test_diff(2,1,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cGp2 ,BC_2L,2,'cGp2' )
     !CALL test_diff(2,1,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cDv2 ,BC_2L,3,'cDv2' )
     !CALL test_diff(2,0,x2p(g2L:(N2+g2U)),x2v(g2L:(N2+g2U)),N2,g2L,g2U,cIpv ,BC_2L,2,'cIpv' )
     !CALL test_diff(2,0,x2v(d2L:(N2+d2U)),x2p(d2L:(N2+d2U)),N2,d2L,d2U,cIvp ,BC_2L,3,'cIvp' )
     
     CALL test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2D,BC_2L,1,'cNv2D')
     CALL test_diff(2,1,x2v(n2L:(N2+n2U)),x2v(n2L:(N2+n2U)),N2,n2L,n2U,cNv2U,BC_2L,1,'cNv2U')
     CALL test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2D,BC_2L,0,'cNp2D')
     CALL test_diff(2,1,x2p(n2L:(N2+n2U)),x2p(n2L:(N2+n2U)),N2,n2L,n2U,cNp2U,BC_2L,0,'cNp2U')
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (iB(1,1) == 1 .AND. iB(2,1) == 1 .AND. dimens == 3) THEN
     CALL test_diff(3,0,x3p              ,x3p              ,N3,b3L,b3U,cFp3 ,BC_3L,0,'cFp3' )
     CALL test_diff(3,0,x3w              ,x3w              ,N3,b3L,b3U,cFw3 ,BC_3L,1,'cFw3' )
     CALL test_diff(3,1,x3p              ,x3p              ,N3,b3L,b3U,cp3  ,BC_3L,0,'cp3'  )
     CALL test_diff(3,1,x3w              ,x3w              ,N3,b3L,b3U,cw3  ,BC_3L,1,'cw3'  )
     CALL test_diff(3,2,x3p              ,x3p              ,N3,b3L,b3U,cp33 ,BC_3L,0,'cp33' )
     CALL test_diff(3,2,x3w              ,x3w              ,N3,b3L,b3U,cw33 ,BC_3L,1,'cw33' )
     
     !CALL test_diff(3,1,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cGp3 ,BC_3L,2,'cGp3' )
     !CALL test_diff(3,1,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cDw3 ,BC_3L,3,'cDw3' )
     !CALL test_diff(3,0,x3p(g3L:(N3+g3U)),x3w(g3L:(N3+g3U)),N3,g3L,g3U,cIpw ,BC_3L,2,'cIpw' )
     !CALL test_diff(3,0,x3w(d3L:(N3+d3U)),x3p(d3L:(N3+d3U)),N3,d3L,d3U,cIwp ,BC_3L,3,'cIwp' )
     
     CALL test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3D,BC_3L,1,'cNw3D')
     CALL test_diff(3,1,x3w(n3L:(N3+n3U)),x3w(n3L:(N3+n3U)),N3,n3L,n3U,cNw3U,BC_3L,1,'cNw3U')
     CALL test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3D,BC_3L,0,'cNp3D')
     CALL test_diff(3,1,x3p(n3L:(N3+n3U)),x3p(n3L:(N3+n3U)),N3,n3L,n3U,cNp3U,BC_3L,0,'cNp3U')
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE test_coeffs
  
  
  
  
  
  
  
  
  
  
  !> \brief calculates differential coefficients
  !! \param[in] dir direction
  !! \param[in] abl ??
  !! \param[in] upwind ?
  !! \param[in] mapping_yes ??
  !! \param[in[ dim_ncb dimension of ??
  !! \param[in] n_coeff_bound ??
  !! \param[in] xC ??
  !! \param[in] xE ??
  !! \param[in] BCL ??
  !! \param[in] BCU ??
  !! \param[in] grid_type ??
  !! \param[in] bl ??
  !! \param[in] bu ??
  !! \param[out] cc coefficients
  SUBROUTINE diff_coeffs(dir,abl,upwind,mapping_yes,dim_ncb,n_coeff_bound,xC,xE,BCL,BCU,grid_type,Nmax,bL,bU,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)   ::  dir
  INTEGER, INTENT(IN)   ::  abl
  INTEGER, INTENT(IN)   ::  upwind
  LOGICAL, INTENT(IN)   ::  mapping_yes
  INTEGER, INTENT(IN)   ::  dim_ncb
  INTEGER, INTENT(IN)   ::  n_coeff_bound(1:dim_ncb)
  INTEGER, INTENT(IN)   ::  Nmax
  INTEGER, INTENT(IN)   ::  bL
  INTEGER, INTENT(IN)   ::  bU
  REAL   , INTENT(IN)   ::  xC(bL:(Nmax+bU))
  REAL   , INTENT(IN)   ::  xE(bL:(Nmax+bU))
  INTEGER, INTENT(IN)   ::  BCL
  INTEGER, INTENT(IN)   ::  BCU
  INTEGER, INTENT(IN)   ::  grid_type
  
  REAL   , INTENT(OUT)  ::  cc(bL:bU,0:Nmax)
  
  INTEGER               ::  n_coeff
  INTEGER               ::  dim_n_coeff_bound
  INTEGER               ::  i, ii, iC, iStart, SShift
  INTEGER               ::  k, kk
  
  INTEGER               ::  left, right
  
  REAL                  ::  dxi(1:2)
  
  REAL                  ::  dxL, dxU ! Für Integrationskoeffizienten
  
  REAL   , ALLOCATABLE  ::  cc_xi(:,:)
  REAL   , ALLOCATABLE  ::  deltaX(:)
  
  LOGICAL               ::  filter_yes
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Upwinding wird auf dem Rand unterdrücken, da Differenzenstencils dort ohnehin schief sind.!
  !              - Generell koennte/sollte die Initialisierung der Index-Grenzen auch schon vor dieser       !
  !                Routine ausgefuehrt werden, um hier die Uebersichtlichkeit zu verbessern.                 !
  !----------------------------------------------------------------------------------------------------------!
  
  
  IF (mapping_yes .AND. abl > 2) THEN
     IF (rank == 0) WRITE(*,*) 'ERROR! Can`t handle derivatives > 2 combined with mapping ...'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  IF (dir == 1) SShift = iShift ! TEST!!! Besser als Argument übergeben ...
  IF (dir == 2) SShift = jShift
  IF (dir == 3) SShift = kShift
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  IF (grid_type == 0 .OR. grid_type == 3) THEN
     iStart = 1
  ELSE
     iStart = 0
  END IF
  !===========================================================================================================
  
  filter_yes = .FALSE.
  IF (abl == 0 .AND. (grid_type == 0 .OR. grid_type == 1)) filter_yes = .TRUE.
  
  
  dim_n_coeff_bound = SIZE(n_coeff_bound)
  
  cc = 0.
  
  DO i = iStart, Nmax
     
     !========================================================================================================
     !=== Stencil-Breite auslesen ============================================================================
     !========================================================================================================
     IF      (BCL > 0 .AND. i <= (iStart - 1 + dim_n_coeff_bound)) THEN
        n_coeff = n_coeff_bound(i + 1 - iStart)
        !IF (upwind /= 0 .AND. i == iStart) n_coeff = 0
     ELSE IF (BCU > 0 .AND. i >= (Nmax   + 1 - dim_n_coeff_bound)) THEN
        n_coeff = n_coeff_bound(Nmax + 1 - i)
        !IF (upwind /= 0 .AND. i == Nmax  ) n_coeff = 0
     ELSE
        n_coeff = n_coeff_bound(dim_n_coeff_bound)
     END IF
     
     !========================================================================================================
     
     IF (n_coeff > 0 .AND. n_coeff > abl) THEN
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten RECHTS vom Entwicklungspunkt bestimmen =================================
        !=====================================================================================================
        IF (grid_type == 0 .OR. grid_type == 1) THEN
           !==================================================================================================
           !=== "normales" Gitter, keine Zwischengitterpunkte ================================================
           !==================================================================================================
           IF      (BCL > 0 .AND. i <= (iStart - 1 + n_coeff/2)) THEN
              right = n_coeff - i + iStart - 1
           ELSE IF (BCU > 0 .AND. i >= (Nmax   + 1 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= ungerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) == 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose odd number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
              END IF
              right = (n_coeff-1)/2
           END IF
        ELSE IF (grid_type == 2) THEN
           !==================================================================================================
           !=== Druck ==> Impuls =============================================================================
           !==================================================================================================
           IF      (BCL > 0 .AND. i <= n_coeff/2) THEN
              right = n_coeff - i - iStart
           ELSE IF (BCU > 0 .AND. i >= (Nmax + 0 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) /= 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
                 CALL MPI_FINALIZE(merror)
                 STOP
              END IF
              right = n_coeff/2
           END IF
        ELSE IF (grid_type == 3) THEN
           !==================================================================================================
           !=== Geschwindigkeit ==> Konti ====================================================================
           !==================================================================================================
           IF      (BCL > 0 .AND. i <= n_coeff/2) THEN
              right = n_coeff - i - iStart
           ELSE IF (BCU > 0 .AND. i >= (Nmax + 1 - n_coeff/2)) THEN
              right = Nmax - i
           ELSE
              ! Ausserhalb des Randes werden zentrale Differenzen gefordert (= gerade Anzahl Koeffizienten)!
              IF (MOD(n_coeff,2) /= 0) THEN
                 IF (rank == 0) THEN
                    WRITE(*,'(a   )') 'ERROR! Choose even number of coefficients!'
                    WRITE(*,'(a,i4)') '    direction =', dir
                    WRITE(*,'(a,i4)') '    grid_type =', grid_type
                    WRITE(*,'(a,i4)') '            i =', i+SShift
                 END IF
                 CALL MPI_FINALIZE(merror)
                 STOP
              END IF
              right = n_coeff/2 - 1
           END IF
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Anzahl der Koeffizienten LINKS vom Entwicklungspunkt bestimmen ==================================
        !=====================================================================================================
        ! (0-ter bzw. 1-ter Koeffizient links vom Entwicklungspunkt
        !          == "zentraler" Koeffizient im Speicher)
        left = right - n_coeff + 1
        !=====================================================================================================
        
        
        !=====================================================================================================
        IF (upwind == -1) THEN
           IF (.NOT. (BCL > 0 .AND. i < (iStart + n_coeff/2))) THEN
              n_coeff = n_coeff - 1
              left    = left    + 1
           END IF
        END IF
        IF (upwind ==  1) THEN
           IF (.NOT. (BCU > 0 .AND. i > (Nmax   - n_coeff/2))) THEN
              n_coeff = n_coeff - 1
              right   = right   - 1
           END IF
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Stencilanordnung testen =========================================================================
        !=====================================================================================================
        IF (right > bU .OR. left < bL) THEN
           IF (rank == 0) THEN
              !WRITE(*,'(a   )') 'WARNING! The FD-Stencil does probably not fit into provided array!'
              WRITE(*,'(a   )') 'ERROR! Stencil doesn`t fit into provided array!'
              WRITE(*,'(a,i4)') '    direction =', dir
              WRITE(*,'(a,i4)') '    grid_type =', grid_type
              WRITE(*,'(a,i4)') '            i =', i+SShift
           END IF
           CALL MPI_FINALIZE(merror)
           STOP
        END IF
        !=====================================================================================================
        
        
        ALLOCATE(deltaX(1:n_coeff))
        
        
        !=====================================================================================================
        !=== räumliche Abstände zum Entwicklungspunkt ========================================================
        !=====================================================================================================
        DO ii = 1, n_coeff
           ! Koeffizienten-Punkte ("iC"):
           iC = i + left + (ii-1)
           
           IF (mapping_yes) THEN
              IF      (grid_type == 2) THEN
                 deltaX(ii) = REAL(iC-i)-0.5
              ELSE IF (grid_type == 3) THEN
                 deltaX(ii) = REAL(iC-i)+0.5
              ELSE
                 deltaX(ii) = REAL(iC-i)
              END IF
           ELSE
              deltaX(ii) = xC(iC) - xE(i)
           END IF
        END DO
        !=====================================================================================================
        
        
        
        !=====================================================================================================
        !=== Integrations-Intervall ==========================================================================
        !=====================================================================================================
        ! Anmerkungen: - Kein Mapping aus Genauigkeitsgründen vorgesehen.
        !              - Nur für Druckgitter vorgesehen.
        IF (abl == -1) THEN
           
           IF (BCL > 0 .AND. i == iStart) THEN
              dxL = 0.
           ELSE
              dxL = (xC(i-1) - xE(i)) / 2. ! der Einfachheit halber anstelle der Geschwindigkeitspunkte ...
           END IF
           IF (BCU > 0 .AND. i == Nmax  ) THEN
              dxU = 0.
           ELSE
              dxU = (xC(i+1) - xE(i)) / 2.
           END IF
           
        END IF
        !=====================================================================================================
        
        
        !=====================================================================================================
        !=== Bestimmung der Koeffizienten ====================================================================
        !=====================================================================================================
        ! Anmerkung: Filter und Integratoren sollen nicht gemappt werden (siehe z.B. Stolz).
        IF (mapping_yes .AND. abl /= 0 .AND. abl /= -1) THEN
           
           !--- derivatives (mapped) ---
           
           ALLOCATE(cc_xi(left:right,1:abl))
           
           cc_xi = 0.
           dxi   = 0.
           
           DO k = 1, abl
              IF (n_coeff <= 7) THEN ! TEST!!!
                 !kk = k
                 !INCLUDE 'FD_coeffs_expl_map.f90'
                 CALL diff_coeffs_exact(k,n_coeff,deltaX(1),cc_xi(left:right,k)) ! TEST!!!
              ELSE
                 CALL FD_coeffs_solver(k,filter_yes,n_coeff,deltaX,cc_xi(left:right,k))
              END IF
              
              DO ii = left, right
                 dxi(k) = dxi(k) + cc_xi(ii,k)*xC(i+ii)
              END DO
           END DO
           
           !---------------------------------------------------!
           ! Mapping:                                          !
           ! d/dx     = (xi') * d/dxi                          !
           ! d^2/dx^2 = (xi") * d/dxi + (xi')**2 * d^2/dxi^2   !
           !                                                   !
           ! Mapping-Faktoren:                                 !
           ! xi' =   1  /  x'                                  !
           ! xi" = - x" / (x')**3                              !
           !---------------------------------------------------!
           
           IF (abl == 1) cc(left:right,i) = cc_xi(left:right,1)/dxi(1)
           IF (abl == 2) cc(left:right,i) = cc_xi(left:right,2)/dxi(1)**2 - cc_xi(left:right,1)*dxi(2)/dxi(1)**3
           
           DEALLOCATE(cc_xi)
           
        ELSE IF (mapping_yes .AND. abl == 0 .AND. .NOT. filter_yes) THEN ! TEST!!!
           
           !--- interpolation (mapped) ---
           
           !k  = 1
           !kk = 0
           !ALLOCATE(cc_xi(left:right,k:k))
           !INCLUDE 'FD_coeffs_expl_map.f90'
           !cc(left:right,i) = cc_xi(left:right,k)
           !DEALLOCATE(cc_xi)
           
           IF (n_coeff <= 7) THEN ! TEST!!!
              CALL diff_coeffs_exact(abl,n_coeff,deltaX(1),cc(left:right,i))
           ELSE
              CALL FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           END IF
           
        ELSE
           
           IF (abl == -1) THEN
              !--- integration coefficients ---
              CALL FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc(left:right,i))
           ELSE
              !--- interpolation and derivatives ---
              CALL FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc(left:right,i))
           END IF
           
        END IF
        !=====================================================================================================
        
        
        DEALLOCATE(deltaX)
        
        
        !=====================================================================================================
        !=== Symmetrie =======================================================================================
        !=====================================================================================================
        IF (BCL == -2 .AND. abl == -1 .AND. i == iStart) THEN
           cc(0,i) = 0.5*cc(0,i)
           DO ii = bL, -1
              cc(ii,i) = 0.
           END DO
        END IF
        IF (BCU == -2 .AND. abl == -1 .AND. i == Nmax  ) THEN
           cc(0,i) = 0.5*cc(0,i)
           DO ii = 1, bU
              cc(ii,i) = 0.
           END DO
        END IF
        !-----------------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------------
        ! TEST!!! abl == -1 ???
        IF (BCL == -2 .AND. (i+bL) < 1) THEN
           IF (grid_type == 0 .OR. (grid_type == 2 .AND. i >= 1)) THEN
              DO ii = 1, 1-(i+bL)
                 cc(1-i+ii,i) = cc(1-i+ii,i) + cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              END DO
           END IF
           IF ((grid_type == 1 .AND. i >= 1) .OR. grid_type == 3) THEN
              DO ii = 1, 1-(i+bL)
                 cc(0-i+ii,i) = cc(0-i+ii,i) - cc(1-i-ii,i)
                 cc(1-i-ii,i) = 0.
              END DO
           END IF
        END IF
        !-----------------------------------------------------------------------------------------------------
        IF (BCU == -2 .AND. (i+bU) > (Nmax+0)) THEN
           IF (grid_type == 0 .OR. (grid_type == 2 .AND. i <= Nmax-1)) THEN
              DO ii = 1, i+bU-(Nmax-0)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) + cc(Nmax-i  +ii,i)
                 cc(Nmax-i  +ii,i) = 0.
              END DO
           END IF
        END IF
        IF (BCU == -2 .AND. (i+bU) > (Nmax-1)) THEN
           IF ((grid_type == 1 .AND. i <= Nmax-1) .OR. grid_type == 3) THEN
              DO ii = 1, i+bU-(Nmax-1)
                 cc(Nmax-i  -ii,i) = cc(Nmax-i  -ii,i) - cc(Nmax-i-1+ii,i)
                 cc(Nmax-i-1+ii,i) = 0.
              END DO
           END IF
        END IF
        !=====================================================================================================
        
     END IF
     
  END DO
  
  
  END SUBROUTINE diff_coeffs
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE diff_coeffs_exact(abl,n_coeff,deltaX,cc) ! TEST!!!
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN   ) ::  abl
  INTEGER, INTENT(IN   ) ::  n_coeff
  REAL   , INTENT(IN   ) ::  deltaX
  REAL   , INTENT(OUT  ) ::  cc(1:n_coeff)
  
  ! TEST!!!
  ! - korrekt ausgerichtet?
  ! - Vorzeichen ok?
  
  
  IF (n_coeff == 2) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,-1./)/2.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 1., 1./)/2.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/-1., 3./)/2.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/-1., 1./)/1.
  END IF
  IF (n_coeff == 3) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/15.,-10., 3./)/8.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,  6.,-1./)/8.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/-1.,  6., 3./)/8.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/ 3.,-10.,15./)/8.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-2.,  3.,-1./)/1.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1.,  1., 0./)/1.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/ 0., -1., 1./)/1.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/ 1., -3., 2./)/1.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-3.,  4.,-1./)/2.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/-1.,  0., 1./)/2.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/ 1., -4., 3./)/2.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/ 1., -2., 1./)/1.
  END IF
  IF (n_coeff == 4) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35.,-35.,  21., -5./)/16.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/  5., 15.,  -5.,  1./)/16.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/ -1.,  9.,   9., -1./)/16.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/  1., -5.,  15.,  5./)/16.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5., 21., -35., 35./)/16.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-71.,141., -93., 23./)/24.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-23., 21.,   3., -1./)/24.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/  1.,-27.,  27., -1./)/24.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/  1., -3., -21., 23./)/24.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/-23., 93.,-141., 71./)/24.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-11., 18.,  -9.,  2./)/6.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/ -2., -3.,   6., -1./)/6.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/  1., -6.,   3.,  2./)/6.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/ -2.,  9., -18., 11./)/6.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/  2., -5.,   4., -1./)/1.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/  1., -2.,   1.,  0./)/1.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/  0.,  1.,  -2.,  1./)/1.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,  4.,  -5.,  2./)/1.
  END IF
  IF (n_coeff == 5) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/315.,-420., 378.,-180., 35./)/128.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35., 140., -70.,  28., -5./)/128.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5.,  60.,  90., -20.,  3./)/128.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/  3., -20.,  90.,  60., -5./)/128.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/ -5.,  28., -70., 140., 35./)/128.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/ 35.,-180., 378.,-420.,315./)/128.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-93., 229.,-225., 111.,-22./)/24.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-22.,  17.,   9.,  -5.,  1./)/24.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/  1., -27.,  27.,  -1.,  0./)/24.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/  0.,   1., -27.,  27., -1./)/24.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/ -1.,   5.,  -9., -17., 22./)/24.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/ 22.,-111., 225.,-229., 93./)/24.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/-25.,  48., -36.,  16., -3./)/12.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/ -3., -10.,  18.,  -6.,  1./)/12.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/  1.,  -8.,   0.,   8., -1./)/12.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/ -1.,   6., -18.,  10.,  3./)/12.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/  3., -16.,  36., -48., 25./)/12.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/ 35.,-104., 114., -56., 11./)/12.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/ 11., -20.,   6.,   4., -1./)/12.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,  16., -30.,  16., -1./)/12.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/ -1.,   4.,   6., -20., 11./)/12.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/ 11., -56., 114.,-104., 35./)/12.
  END IF
  IF (n_coeff == 6) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/  693.,-1155.,  1386., -990.,   385., -63./)/256.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/   63.,  315.,  -210.,  126.,   -45.,   7./)/256.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/   -7.,  105.,   210.,  -70.,    21.,  -3./)/256.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/    3.,  -25.,   150.,  150.,   -25.,   3./)/256.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/   -3.,   21.,   -70.,  210.,   105.,  -7./)/256.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/    7.,  -45.,   126., -210.,   315.,  63./)/256.
     IF (deltaX == -5.5 .AND. abl == 0) cc(1:n_coeff) = (/  -63.,  385.,  -990., 1386., -1155., 693./)/256.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-9129.,26765.,-34890.,25770.,-10205.,1689./)/1920.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/-1689., 1005.,  1430.,-1110.,   435., -71./)/1920.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/   71.,-2115.,  2070.,   10.,   -45.,   9./)/1920.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/   -9.,  125., -2250., 2250.,  -125.,   9./)/1920.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/   -9.,   45.,   -10.,-2070.,  2115., -71./)/1920.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/   71., -435.,  1110.,-1430., -1005.,1689./)/1920.
     IF (deltaX == -5.5 .AND. abl == 1) cc(1:n_coeff) = (/-1689.,10205.,-25770.,34890.,-26765.,9129./)/1920.
     
     IF (deltaX ==  0.0 .AND. abl == 1) cc(1:n_coeff) = (/ -137.,  300.,  -300.,  200.,   -75.,  12./)/60.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/  -12.,  -65.,   120.,  -60.,    20.,  -3./)/60.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/    3.,  -30.,   -20.,   60.,   -15.,   2./)/60.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/   -2.,   15.,   -60.,   20.,    30.,  -3./)/60.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/    3.,  -20.,    60., -120.,    65.,  12./)/60.
     IF (deltaX == -5.0 .AND. abl == 1) cc(1:n_coeff) = (/  -12.,   75.,  -200.,  300.,  -300., 137./)/60.
     
     IF (deltaX ==  0.0 .AND. abl == 2) cc(1:n_coeff) = (/   45., -154.,   214., -156.,    61., -10./)/12.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/   10.,  -15.,    -4.,   14.,    -6.,   1./)/12.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/   -1.,   16.,   -30.,   16.,    -1.,   0./)/12.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/    0.,   -1.,    16.,  -30.,    16.,  -1./)/12.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/    1.,   -6.,    14.,   -4.,   -15.,  10./)/12.
     IF (deltaX == -5.0 .AND. abl == 2) cc(1:n_coeff) = (/  -10.,   61.,  -156.,  214.,  -154.,  45./)/12.
  END IF
  IF (n_coeff == 7) THEN
     IF (deltaX ==  0.5 .AND. abl == 0) cc(1:n_coeff) = (/  3003., -6006.,  9009.,-8580.,   5005., -1638.,  231./)/1024.
     IF (deltaX == -0.5 .AND. abl == 0) cc(1:n_coeff) = (/   231.,  1386., -1155.,  924.,   -495.,   154.,  -21./)/1024.
     IF (deltaX == -1.5 .AND. abl == 0) cc(1:n_coeff) = (/   -21.,   378.,   945., -420.,    189.,   -54.,    7./)/1024.
     IF (deltaX == -2.5 .AND. abl == 0) cc(1:n_coeff) = (/     7.,   -70.,   525.,  700.,   -175.,    42.,   -5./)/1024.
     IF (deltaX == -3.5 .AND. abl == 0) cc(1:n_coeff) = (/    -5.,    42.,  -175.,  700.,    525.,   -70.,    7./)/1024.
     IF (deltaX == -4.5 .AND. abl == 0) cc(1:n_coeff) = (/     7.,   -54.,   189., -420.,    945.,   378.,  -21./)/1024.
     IF (deltaX == -5.5 .AND. abl == 0) cc(1:n_coeff) = (/   -21.,   154.,  -495.,  924.,  -1155.,  1386.,  231./)/1024.
     IF (deltaX == -6.5 .AND. abl == 0) cc(1:n_coeff) = (/   231., -1638.,  5005.,-8580.,   9009., -6006., 3003./)/1024.
     
     IF (deltaX ==  0.5 .AND. abl == 1) cc(1:n_coeff) = (/-10756., 36527.,-59295., 58310.,-34610., 11451.,-1627./)/1920.
     IF (deltaX == -0.5 .AND. abl == 1) cc(1:n_coeff) = (/ -1627.,   633.,  2360., -2350.,  1365.,  -443.,   62./)/1920.
     IF (deltaX == -1.5 .AND. abl == 1) cc(1:n_coeff) = (/    62., -2061.,  1935.,   190.,  -180.,    63.,   -9./)/1920.
     IF (deltaX == -2.5 .AND. abl == 1) cc(1:n_coeff) = (/    -9.,   125., -2250.,  2250.,  -125.,     9.,    0./)/1920.
     IF (deltaX == -3.5 .AND. abl == 1) cc(1:n_coeff) = (/     0.,    -9.,   125., -2250.,  2250.,  -125.,    9./)/1920.
     IF (deltaX == -4.5 .AND. abl == 1) cc(1:n_coeff) = (/     9.,   -63.,   180.,  -190., -1935.,  2061.,  -62./)/1920.
     IF (deltaX == -5.5 .AND. abl == 1) cc(1:n_coeff) = (/   -62.,   443., -1365.,  2350., -2360.,  -633., 1627./)/1920.
     IF (deltaX == -6.5 .AND. abl == 1) cc(1:n_coeff) = (/  1627.,-11451., 34610.,-58310., 59295.,-36527.,10756./)/1920.
     
     IF (deltaX == -0.0 .AND. abl == 1) cc(1:n_coeff) = (/  -147.,   360.,  -450.,   400.,  -225.,    72.,  -10./)/60.
     IF (deltaX == -1.0 .AND. abl == 1) cc(1:n_coeff) = (/   -10.,   -77.,   150.,  -100.,    50.,   -15.,    2./)/60.
     IF (deltaX == -2.0 .AND. abl == 1) cc(1:n_coeff) = (/     2.,   -24.,   -35.,    80.,   -30.,     8.,   -1./)/60.
     IF (deltaX == -3.0 .AND. abl == 1) cc(1:n_coeff) = (/    -1.,     9.,   -45.,     0.,    45.,    -9.,    1./)/60.
     IF (deltaX == -4.0 .AND. abl == 1) cc(1:n_coeff) = (/     1.,    -8.,    30.,   -80.,    35.,    24.,   -2./)/60.
     IF (deltaX == -5.0 .AND. abl == 1) cc(1:n_coeff) = (/    -2.,    15.,   -50.,   100.,  -150.,    77.,   10./)/60.
     IF (deltaX == -6.0 .AND. abl == 1) cc(1:n_coeff) = (/    10.,   -72.,   225.,  -400.,   450.,  -360.,  147./)/60.
     
     IF (deltaX == -0.0 .AND. abl == 2) cc(1:n_coeff) = (/   812., -3132.,  5265., -5080.,  2970.,  -972.,  137./)/180.
     IF (deltaX == -1.0 .AND. abl == 2) cc(1:n_coeff) = (/   137.,  -147.,  -225.,   470.,  -285.,    93.,  -13./)/180.
     IF (deltaX == -2.0 .AND. abl == 2) cc(1:n_coeff) = (/   -13.,   228.,  -420.,   200.,    15.,   -12.,    2./)/180.
     IF (deltaX == -3.0 .AND. abl == 2) cc(1:n_coeff) = (/     2.,   -27.,   270.,  -490.,   270.,   -27.,    2./)/180.
     IF (deltaX == -4.0 .AND. abl == 2) cc(1:n_coeff) = (/     2.,   -12.,    15.,   200.,  -420.,   228.,  -13./)/180.
     IF (deltaX == -5.0 .AND. abl == 2) cc(1:n_coeff) = (/   -13.,    93.,  -285.,   470.,  -225.,  -147.,  137./)/180.
     IF (deltaX == -6.0 .AND. abl == 2) cc(1:n_coeff) = (/   137.,  -972.,  2970., -5080.,  5265., -3132.,  812./)/180.
  END IF
  
  
  END SUBROUTINE diff_coeffs_exact
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE test_diff(dir,abl,xC,xE,Nmax,bL,bU,cc,BCL,grid_type,name)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)   ::  dir
  INTEGER, INTENT(IN)   ::  abl
  INTEGER, INTENT(IN)   ::  Nmax
  INTEGER, INTENT(IN)   ::  bL
  INTEGER, INTENT(IN)   ::  bU
  REAL   , INTENT(IN)   ::  xC(bL:(Nmax+bU))
  REAL   , INTENT(IN)   ::  xE(bL:(Nmax+bU))
  REAL   , INTENT(IN)   ::  cc(bL:bU,0:Nmax)
  INTEGER, INTENT(IN)   ::  BCL
  INTEGER, INTENT(IN)   ::  grid_type
  
  CHARACTER(*), INTENT(IN) ::  name
  
  INTEGER               ::  i, ii, iStartC, iStartE, SShift
  INTEGER               ::  nw, nw_max, n_dx
  REAL                  ::  wave, wave_mod(1:2), dx
  
  CHARACTER(LEN=1)      ::  part
  REAL                  ::  pi
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Hier sollen alle Differenzen-Stencils innerhalb eines Blockes getestet werden, analog zu  !
  !                ihrer Berechnung. Daher wird bewusst nicht mittels MPI in eine Datei geschrieben, sondern !
  !                in verschiedene Blöcke, um Fehler bei der MPI-Programmierung auszuschliessen.             !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !--- pi ---
  pi = 2.*ABS(ACOS(0.))
  
  
  IF (dir == 1) THEN
     WRITE(part,'(i1.1)') iB(1,1)
     SShift = iShift
  ELSE IF (dir == 2) THEN
     WRITE(part,'(i1.1)') iB(2,1)
     SShift = jShift
  ELSE IF (dir == 3) THEN
     WRITE(part,'(i1.1)') iB(3,1)
     SShift = kShift
  ELSE
     IF (rank == 0) WRITE(*,*) 'ERROR! Wrong input at subroutine `test_diff`!'
     CALL MPI_FINALIZE(merror)
     STOP
  END IF
  
  
  OPEN(10,FILE='test_'//name//'_transfer_block'//part//'.txt',STATUS='UNKNOWN')
  OPEN(11,FILE='test_'//name//'_coeffs_block'  //part//'.txt',STATUS='UNKNOWN')
  
  
  !===========================================================================================================
  !=== Startindex für Entwicklungspunkte =====================================================================
  !===========================================================================================================
  IF (BCL <= 0 .OR. grid_type == 0 .OR. grid_type == 3) THEN
     iStartE = 1
  ELSE
     iStartE = 0
  END IF
  
  
  !===========================================================================================================
  !=== Startindex für Koeffizienten ==========================================================================
  !===========================================================================================================
  IF (BCL > 0 .AND. (grid_type == 1 .OR. grid_type == 3)) THEN
     iStartC = 0
  ELSE
     iStartC = 1
  END IF
  
  
  !===========================================================================================================
  !=== Test der Koeffizienten ================================================================================
  !===========================================================================================================
  nw_max = 100
  
  DO i = iStartE, Nmax
     
     DO nw = 1, nw_max
        
        wave     = pi*REAL(nw)/REAL(nw_max)
        wave_mod = 0.
        
        !=== Referenz-Gitterweite bestimmen ==================================================================
        dx   = 0.
        n_dx = 0
        
        DO ii = bL, bU-1
           !IF (i+ii .GE. iStartC .AND. i+ii .LT. Nmax .AND. ABS(xC(i+ii+1)-xC(i+ii)) .GT. dx) THEN
           IF (i+ii >= iStartC .AND. i+ii < Nmax) THEN
              dx   = dx + ABS(xC(i+ii+1)-xC(i+ii))
              n_dx = n_dx + 1
           END IF
        END DO
        
        dx = dx / REAL(n_dx)
        
        
        !=== Transferfunktion plotten ========================================================================
        IF (abl == 0) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))
              wave_mod(2) = wave_mod(2) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))
           END DO
        ELSE IF (abl == 1) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) + cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx
              wave_mod(2) = wave_mod(2) + cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx/wave**2
           END DO
        ELSE IF (abl == 2) THEN
           DO ii = bL, bU
              wave_mod(1) = wave_mod(1) - cc(ii,i)*cos(wave/dx*(xC(i+ii)-xE(i)))*dx**2
              wave_mod(2) = wave_mod(2) - cc(ii,i)*sin(wave/dx*(xC(i+ii)-xE(i)))*dx**2/wave
           END DO
        END IF
        
        WRITE(10,'(i7,3E26.17)') i+SShift, wave, wave_mod(1:2)
        
     END DO
     
     WRITE(10,*)
     WRITE(11,'(i7,100E26.17)') i+SShift, cc(:,i)
     
  END DO
  
  CLOSE(10)
  CLOSE(11)
  
  
  END SUBROUTINE test_diff
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE FD_coeffs_solver(abl,filter_yes,n_coeff,deltaX,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)   ::  abl
  LOGICAL, INTENT(IN)   ::  filter_yes
  INTEGER, INTENT(IN)   ::  n_coeff
  REAL   , INTENT(IN)   ::  deltaX(1:n_coeff)
  REAL   , INTENT(OUT)  ::  cc    (1:n_coeff)
  
  INTEGER               ::  i, j
  
  REAL                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  REAL                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  REAL                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  DO i = 1, n_coeff
     DO j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     END DO
  END DO
  
  
  !===========================================================================================================
  !=== Zusatzbedingungen =====================================================================================
  !===========================================================================================================
  IF (filter_yes) THEN
     ! G(pi) = 0.
     DO i = 1, n_coeff
        polyn_vals(i,n_coeff) = (-1.)**i
     END DO
  END IF
  
  
  !===========================================================================================================
  !=== L�sen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  CALL Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Funktionswert am Entwicklungspunkt (deltaX = 0) =======================================================
  !===========================================================================================================
  const = 1.
  DO i = 1, abl-1
     const = const*REAL(1+i)
  END DO
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  cc(1:n_coeff) = const*polyn_vals_inv(1+abl,1:n_coeff)
  
  
  END SUBROUTINE FD_coeffs_solver
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE FD_coeffs_solver_integral(n_coeff,deltaX,dxL,dxU,cc)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)   ::  n_coeff
  REAL   , INTENT(IN)   ::  deltaX(1:n_coeff)
  REAL   , INTENT(OUT)  ::  cc    (1:n_coeff)
  REAL   , INTENT(IN)   ::  dxL, dxU
  
  INTEGER               ::  i, j, k
  
  REAL                  ::  polyn_vals    (1:n_coeff,1:n_coeff)
  REAL                  ::  polyn_vals_inv(1:n_coeff,1:n_coeff)
  
  REAL                  ::  const
  
  
  !===========================================================================================================
  !=== Aufstellen des Gleichungssystems ======================================================================
  !===========================================================================================================
  DO i = 1, n_coeff
     DO j = 1, n_coeff
        polyn_vals(i,j) = deltaX(i)**(j-1)
     END DO
  END DO
  
  
  !===========================================================================================================
  !=== L�sen des Gleichungssystems ===========================================================================
  !===========================================================================================================
  CALL Matrix_invert(n_coeff,polyn_vals,polyn_vals_inv)
  
  
  !===========================================================================================================
  !=== Koeffizienten bestimmen (explizite Differenzen) =======================================================
  !===========================================================================================================
  !  (Matrix-Vektor-Multiplikation)
  cc = 0.
  DO j = 1, n_coeff
     DO i = 1, n_coeff
        cc(j) = cc(j) + (dxU**i - dxL**i)/REAL(i)*polyn_vals_inv(i,j)
     END DO
  END DO
  
  
  END SUBROUTINE FD_coeffs_solver_integral
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE get_stencil
  
  IMPLICIT NONE
  
  INTEGER                ::  i, imax
  INTEGER                ::  j, jmax
  INTEGER                ::  k, kmax
  
  INTEGER                ::  g
  
  REAL                   ::  cDu1R(-1:0,1:N1)
  REAL                   ::  cDv2R(-1:0,1:N2)
  REAL                   ::  cDw3R(-1:0,1:N3)
  
  REAL                   ::  cGp1R( 0:1,0:N1)
  REAL                   ::  cGp2R( 0:1,0:N2)
  REAL                   ::  cGp3R( 0:1,0:N3)
  
  
  cdg1 = 0.
  cdg2 = 0.
  cdg3 = 0.
  
  DO g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     cDu1R = 0.
     cGp1R = 0.
     
     !---------------------------------------------------------------------------!
     ! Achtung: Falls nicht periodisch, wird hier am Rand ein Fehler gemacht,    !
     !          der nur über die Sonderbehandlung weiter unten korrigiert wird!  !
     !---------------------------------------------------------------------------!
     
     DO i = 1, imax
        cDu1R(-1,i) = -1./(x1uR(i  ,g) - x1uR(i-1,g))
        cDu1R( 0,i) =  1./(x1uR(i  ,g) - x1uR(i-1,g))
     END DO
     
     DO i = 1, imax-1
        cGp1R( 0,i) = -1./(x1pR(i+1,g) - x1pR(i  ,g))
        cGp1R( 1,i) =  1./(x1pR(i+1,g) - x1pR(i  ,g))
     END DO
     
     IF (BC(1,1,g) > 0) THEN
        cGp1R( :,0   ) =  0.
        cDu1R( :,1   ) =  0.
        !cDu1R( 0,1   ) =  1./(x1u (1  ) - x1u (0  ))
        cDu1R( 0,1   ) =  1./(y1u (1  ) - y1u (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     ELSE
        cGp1R( 0,0   ) = -1./(x1pR(1,g) - x1pR(0,g))
        cGp1R( 1,0   ) =  1./(x1pR(1,g) - x1pR(0,g))
     END IF
     
     IF (BC(2,1,g) > 0) THEN
        cGp1R( :,imax) =  0.
        cDu1R( :,imax) =  0.
        !cDu1R(-1,imax) = -1./(x1u (N1      ) - x1u (N1-1  )) ! TEST!!! Das geht in die Hose ...
        cDu1R(-1,imax) = -1./(y1u (M1      ) - y1u (M1-1  ))
     ELSE
        cGp1R( 0,imax) = -1./(x1pR(imax+1,g) - x1pR(imax,g))
        cGp1R( 1,imax) =  1./(x1pR(imax+1,g) - x1pR(imax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO i = 1, imax
        cdg1(-1,i,g) = cDu1R(-1,i)*cGp1R(0,i-1)
        cdg1( 0,i,g) = cDu1R(-1,i)*cGp1R(1,i-1) + cDu1R(0,i)*cGp1R(0,i)
        cdg1( 1,i,g) =                            cDu1R(0,i)*cGp1R(1,i)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     IF (BC(1,1,g) > 0) cdg1( :,1   ,g) = 2.*cdg1(:,1   ,g) ! TEST!!! Ist das wirklich so optimal?
     IF (BC(2,1,g) > 0) cdg1( :,imax,g) = 2.*cdg1(:,imax,g)
     
     IF (BC(1,1,g) == -2) THEN
        cdg1( 1,1   ,g) = cdg1( 1,1   ,g) + cdg1(-1,1   ,g)
        cdg1(-1,1   ,g) = 0.
     END IF
     IF (BC(2,1,g) == -2) THEN
        cdg1(-1,imax,g) = cdg1(-1,imax,g) + cdg1( 1,imax,g)
        cdg1( 1,imax,g) = 0.
     END IF
     !========================================================================================================
     jmax = NN(2,g)
     
     cDv2R = 0.
     cGp2R = 0.
     
     DO j = 1, jmax
        cDv2R(-1,j) = -1./(x2vR(j  ,g) - x2vR(j-1,g))
        cDv2R( 0,j) =  1./(x2vR(j  ,g) - x2vR(j-1,g))
     END DO
     
     DO j = 1, jmax-1
        cGp2R( 0,j) = -1./(x2pR(j+1,g) - x2pR(j  ,g))
        cGp2R( 1,j) =  1./(x2pR(j+1,g) - x2pR(j  ,g))
     END DO
     
     IF (BC(1,2,g) > 0) THEN
        cGp2R( :,0   ) =  0.
        cDv2R( :,1   ) =  0.
        !cDv2R( 0,1   ) =  1./(x2v (1  ) - x2v (0    ))
        cDv2R( 0,1   ) =  1./(y2v (1  ) - y2v (0    )) ! TEST!!! Der Schoenheit halber, s.u. ! TEST!!! Das Ergebnis ist nicht 100% korrekt. Warum???
     ELSE
        cGp2R( 0,0   ) = -1./(x2pR(1,g) - x2pR(0  ,g))
        cGp2R( 1,0   ) =  1./(x2pR(1,g) - x2pR(0  ,g))
     END IF
     
     IF (BC(2,2,g) > 0) THEN
        cGp2R( :,jmax) =  0.
        cDv2R( :,jmax) =  0.
        !cDv2R(-1,jmax) = -1./(x2v (N2      ) - x2v (N2-1  )) ! TEST!!! Das geht in die Hose ...
        cDv2R(-1,jmax) = -1./(y2v (M2      ) - y2v (M2-1  ))
     ELSE
        cGp2R( 0,jmax) = -1./(x2pR(jmax+1,g) - x2pR(jmax,g))
        cGp2R( 1,jmax) =  1./(x2pR(jmax+1,g) - x2pR(jmax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO j = 1, jmax
        cdg2(-1,j,g) = cDv2R(-1,j)*cGp2R(0,j-1)
        cdg2( 0,j,g) = cDv2R(-1,j)*cGp2R(1,j-1) + cDv2R(0,j)*cGp2R(0,j)
        cdg2( 1,j,g) =                            cDv2R(0,j)*cGp2R(1,j)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     IF (BC(1,2,g) > 0) cdg2( :,1   ,g) = 2.*cdg2(:,1   ,g)
     IF (BC(2,2,g) > 0) cdg2( :,jmax,g) = 2.*cdg2(:,jmax,g)
     
     IF (BC(1,2,g) == -2) THEN
        cdg2( 1,1   ,g) = cdg2( 1,1   ,g) + cdg2(-1,1   ,g)
        cdg2(-1,1   ,g) = 0.
     END IF
     IF (BC(2,2,g) == -2) THEN
        cdg2(-1,jmax,g) = cdg2(-1,jmax,g) + cdg2( 1,jmax,g)
        cdg2( 1,jmax,g) = 0.
     END IF
     !========================================================================================================
     IF (dimens == 3) THEN
     
     kmax = NN(3,g)
     
     cDw3R = 0.
     cGp3R = 0.
     
     DO k = 1, kmax
        cDw3R(-1,k) = -1./(x3wR(k  ,g) - x3wR(k-1,g))
        cDw3R( 0,k) =  1./(x3wR(k  ,g) - x3wR(k-1,g))
     END DO
     
     DO k = 1, kmax-1
        cGp3R( 0,k) = -1./(x3pR(k+1,g) - x3pR(k  ,g))
        cGp3R( 1,k) =  1./(x3pR(k+1,g) - x3pR(k  ,g))
     END DO
     
     IF (BC(1,3,g) > 0) THEN
        cGp3R( :,0   ) =  0.
        cDw3R( :,1   ) =  0.
        !cDw3R( 0,1   ) =  1./(x3w (1  ) - x3w (0  ))
        cDw3R( 0,1   ) =  1./(y3w (1  ) - y3w (0  )) ! TEST!!! Der Schoenheit halber, s.u.
     ELSE
        cGp3R( 0,0   ) = -1./(x3pR(1,g) - x3pR(0,g))
        cGp3R( 1,0   ) =  1./(x3pR(1,g) - x3pR(0,g))
     END IF
     
     IF (BC(2,3,g) > 0) THEN
        cGp3R( :,kmax) =  0.
        cDw3R( :,kmax) =  0.
        !cDw3R(-1,kmax) = -1./(x3w (N3      ) - x3w (N3-1  )) ! TEST!!! Das geht in die Hose ...
        cDw3R(-1,kmax) = -1./(y3w (M3      ) - y3w (M3-1  ))
     ELSE
        cGp3R( 0,kmax) = -1./(x3pR(kmax+1,g) - x3pR(kmax,g))
        cGp3R( 1,kmax) =  1./(x3pR(kmax+1,g) - x3pR(kmax,g))
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO k = 1, kmax
        cdg3(-1,k,g) = cDw3R(-1,k)*cGp3R(0,k-1)
        cdg3( 0,k,g) = cDw3R(-1,k)*cGp3R(1,k-1) + cDw3R(0,k)*cGp3R(0,k)
        cdg3( 1,k,g) =                            cDw3R(0,k)*cGp3R(1,k)
     END DO
     
     ! Faktor 2 kommt von der Approximation DH-¹D ~= DI-¹G, wobei I hier das Interpolationspolynom am Rand darstellt
     IF (BC(1,3,g) > 0) cdg3( :,1   ,g) = 2.*cdg3(:,1   ,g)
     IF (BC(2,3,g) > 0) cdg3( :,kmax,g) = 2.*cdg3(:,kmax,g)
     
     IF (BC(1,3,g) == -2) THEN
        cdg3( 1,1   ,g) = cdg3( 1,1   ,g) + cdg3(-1,1   ,g)
        cdg3(-1,1   ,g) = 0.
     END IF
     IF (BC(2,3,g) == -2) THEN
        cdg3(-1,kmax,g) = cdg3(-1,kmax,g) + cdg3( 1,kmax,g)
        cdg3( 1,kmax,g) = 0.
     END IF
     
     END IF
     !========================================================================================================
     
  END DO
  
  
  END SUBROUTINE get_stencil
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE get_stencil_transp
  
  IMPLICIT NONE
  
  INTEGER                ::  i, j, k, g
  REAL   , ALLOCATABLE   ::  work(:,:)
  
  
  DO g = 1, n_grids
     
     !========================================================================================================
     ALLOCATE(work(-1:1,0:(NN(1,g)+1)))
     
     work = 0.
     DO i = 1, NN(1,g)
        work(-1:1,i) = cdg1(-1:1,i,g)
     END DO
     
     cdg1(:,:,g) = 0.
     DO i = 1, NN(1,g)
        cdg1(-1,i,g) = work( 1,i-1)
        cdg1( 0,i,g) = work( 0,i  )
        cdg1( 1,i,g) = work(-1,i+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,1,g) == 0 .OR. BC(1,1,g) == -1) THEN
        i = 1
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     END IF
     IF (BC(2,1,g) == 0 .OR. BC(2,1,g) == -1) THEN
        i = NN(1,g)
        cdg1(-1,i,g) = 1./(x1uR(i-1,g) - x1uR(i-2,g))/(x1pR(i  ,g) - x1pR(i-1,g))
        cdg1( 1,i,g) = 1./(x1uR(i+1,g) - x1uR(i  ,g))/(x1pR(i+1,g) - x1pR(i  ,g))
     END IF
     !========================================================================================================
     ALLOCATE(work(-1:1,0:(NN(2,g)+1)))
     
     work = 0.
     DO j = 1, NN(2,g)
        work(-1:1,j) = cdg2(-1:1,j,g)
     END DO
     
     cdg2(:,:,g) = 0.
     DO j = 1, NN(2,g)
        cdg2(-1,j,g) = work( 1,j-1)
        cdg2( 0,j,g) = work( 0,j  )
        cdg2( 1,j,g) = work(-1,j+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,2,g) == 0 .OR. BC(1,2,g) == -1) THEN
        j = 1
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     END IF
     IF (BC(2,2,g) == 0 .OR. BC(2,2,g) == -1) THEN
        j = NN(2,g)
        cdg2(-1,j,g) = 1./(x2vR(j-1,g) - x2vR(j-2,g))/(x2pR(j  ,g) - x2pR(j-1,g))
        cdg2( 1,j,g) = 1./(x2vR(j+1,g) - x2vR(j  ,g))/(x2pR(j+1,g) - x2pR(j  ,g))
     END IF
     !========================================================================================================
     IF (dimens == 3) THEN
     
     ALLOCATE(work(-1:1,0:(NN(3,g)+1)))
     
     work = 0.
     DO k = 1, NN(3,g)
        work(-1:1,k) = cdg3(-1:1,k,g)
     END DO
     
     cdg3(:,:,g) = 0.
     DO k = 1, NN(3,g)
        cdg3(-1,k,g) = work( 1,k-1)
        cdg3( 0,k,g) = work( 0,k  )
        cdg3( 1,k,g) = work(-1,k+1)
     END DO
     
     DEALLOCATE(work)
     
     IF (BC(1,3,g) == 0 .OR. BC(1,3,g) == -1) THEN
        k = 1
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     END IF
     IF (BC(2,3,g) == 0 .OR. BC(2,3,g) == -1) THEN
        k = NN(3,g)
        cdg3(-1,k,g) = 1./(x3wR(k-1,g) - x3wR(k-2,g))/(x3pR(k  ,g) - x3pR(k-1,g))
        cdg3( 1,k,g) = 1./(x3wR(k+1,g) - x3wR(k  ,g))/(x3pR(k+1,g) - x3pR(k  ,g))
     END IF
     
     END IF
     !========================================================================================================
     
  END DO
  
  
  ! sicherheitshalber (cdg3 = 0. sollte eigentlich schon erf�llt sein):
  IF (dimens == 2) cdg3 = 0.
  
  
  END SUBROUTINE get_stencil_transp
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE get_stencil_Helm
  
  IMPLICIT NONE
  
  INTEGER                ::  imax, jmax, kmax
  INTEGER                ::  g, m
  
  
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
  DO g = 1, n_grids
     
     !========================================================================================================
     imax = NN(1,g)
     
     CALL diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1pR(-1:imax+1,g),x1pR(-1:imax+1,g),BC_1L,BC_1U,0,imax,-1,1,cp11R(-1,0,g))
     
     IF (BC_1L == 1 .OR. BC_1L == 3) THEN
        cp11R( :,1   ,g) =  0.
        cp11R( 0,1   ,g) =  1.
     END IF
     IF (BC_1U == 1 .OR. BC_1U == 3) THEN
        cp11R( :,imax,g) =  0.
        cp11R( 0,imax,g) =  1.
     END IF
     
     IF (BC_1L == 2 .OR. BC_1L == 4) THEN
        cp11R(-1,1   ,g) =  0.
        cp11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cp11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     END IF
     IF (BC_1U == 2 .OR. BC_1U == 4) THEN
        cp11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cp11R( 1,imax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cu11R(:,:,g) = cp11R(:,:,g)
     
     IF (BC_1L ==  2) THEN
        cu11R( :,1   ,g) =  0.
        cu11R( 0,1   ,g) =  1.
     END IF
     IF (BC_1U ==  2) THEN
        cu11R( :,imax,g) =  0.
        cu11R( 0,imax,g) =  1.
     END IF
     
     IF (BC_1L ==  3) THEN
        cu11R(-1,1   ,g) =  0.
        cu11R( 0,1   ,g) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
        cu11R( 1,1   ,g) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
     END IF
     IF (BC_1U ==  3) THEN
        cu11R(-1,imax,g) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 0,imax,g) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
        cu11R( 1,imax,g) =  0.
     END IF
     
     IF (BC_1L == -2) THEN
        cu11R(-1,1   ,g) =  0.
        cu11R( 1,1   ,g) =  0.
     END IF
     IF (BC_1U == -2) THEN
        cu11R(-1,imax,g) =  0.
        cu11R( 1,imax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO m = 1, n_conc
        cc11R(:,:,g,m) = cp11R(:,:,g)
        
        IF (BCc_1L(m) > 0) THEN
           cc11R(-1,1   ,g,m) =  0.
           cc11R( 0,1   ,g,m) = -1./(x1pR(2   ,g) - x1pR(1     ,g))
           cc11R( 1,1   ,g,m) =  1./(x1pR(2   ,g) - x1pR(1     ,g))
        END IF
        IF (BCc_1U(m) > 0) THEN
           cc11R(-1,imax,g,m) = -1./(x1pR(imax,g) - x1pR(imax-1,g))
           cc11R( 0,imax,g,m) =  1./(x1pR(imax,g) - x1pR(imax-1,g))
           cc11R( 1,imax,g,m) =  0.
        END IF
     END DO
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     jmax = NN(2,g)
     
     CALL diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2pR(-1:jmax+1,g),x2pR(-1:jmax+1,g),BC_2L,BC_2U,0,jmax,-1,1,cp22R(-1,0,g))
     
     IF (BC_2L == 1 .OR. BC_2L == 3) THEN
        cp22R(:,1    ,g) =  0.
        cp22R(0,1    ,g) =  1.
     END IF
     IF (BC_2U == 1 .OR. BC_2U == 3) THEN
        cp22R( :,jmax,g) =  0.
        cp22R( 0,jmax,g) =  1.
     END IF
     
     IF (BC_2L == 2 .OR. BC_2L == 4) THEN
        cp22R(-1,1   ,g) =  0.
        cp22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cp22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     END IF
     IF (BC_2U == 2 .OR. BC_2U == 4) THEN
        cp22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cp22R( 1,jmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cv22R(:,:,g) = cp22R(:,:,g)
     
     IF (BC_2L ==  2) THEN
        cv22R( :,1   ,g) =  0.
        cv22R( 0,1   ,g) =  1.
     END IF
     IF (BC_2U ==  2) THEN
        cv22R( :,jmax,g) =  0.
        cv22R( 0,jmax,g) =  1.
     END IF
     
     IF (BC_2L ==  3) THEN
        cv22R(-1,1   ,g) =  0.
        cv22R( 0,1   ,g) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
        cv22R( 1,1   ,g) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
     END IF
     IF (BC_2U ==  3) THEN
        cv22R(-1,jmax,g) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 0,jmax,g) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
        cv22R( 1,jmax,g) =  0.
     END IF
     
     IF (BC_2L == -2) THEN
        cv22R(-1,1   ,g) =  0.
        cv22R( 1,1   ,g) =  0.
     END IF
     IF (BC_2U == -2) THEN
        cv22R(-1,jmax,g) =  0.
        cv22R( 1,jmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO m = 1, n_conc
        cc22R(:,:,g,m) = cp22R(:,:,g)
        
        IF (BCc_2L(m) > 0) THEN
           cc22R(-1,1   ,g,m) =  0.
           cc22R( 0,1   ,g,m) = -1./(x2pR(2   ,g) - x2pR(1     ,g))
           cc22R( 1,1   ,g,m) =  1./(x2pR(2   ,g) - x2pR(1     ,g))
        END IF
        IF (BCc_2U(m) > 0) THEN
           cc22R(-1,jmax,g,m) = -1./(x2pR(jmax,g) - x2pR(jmax-1,g))
           cc22R( 0,jmax,g,m) =  1./(x2pR(jmax,g) - x2pR(jmax-1,g))
           cc22R( 1,jmax,g,m) =  0.
        END IF
     END DO
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     
     kmax = NN(3,g)
     
     CALL diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3pR(-1:kmax+1,g),x3pR(-1:kmax+1,g),BC_3L,BC_3U,0,kmax,-1,1,cp33R(-1,0,g))
     
     IF (BC_3L == 1 .OR. BC_3L == 3) THEN
        cp33R( :,1   ,g) =  0.
        cp33R( 0,1   ,g) =  1.
     END IF
     IF (BC_3U == 1 .OR. BC_3U == 3) THEN
        cp33R( :,kmax,g) =  0.
        cp33R( 0,kmax,g) =  1.
     END IF
     
     IF (BC_3L == 2 .OR. BC_3L == 4) THEN
        cp33R(-1,1   ,g) =  0.
        cp33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cp33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     END IF
     IF (BC_3U == 2 .OR. BC_3U == 4) THEN
        cp33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cp33R( 1,kmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     cw33R(:,:,g) = cp33R(:,:,g)
     
     IF (BC_3L ==  2) THEN
        cw33R( :,1   ,g) =  0.
        cw33R( 0,1   ,g) =  1.
     END IF
     IF (BC_3U ==  2) THEN
        cw33R( :,kmax,g) =  0.
        cw33R( 0,kmax,g) =  1.
     END IF
     
     IF (BC_3L ==  3) THEN
        cw33R(-1,1   ,g) =  0.
        cw33R( 0,1   ,g) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
        cw33R( 1,1   ,g) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
     END IF
     IF (BC_3U ==  3) THEN
        cw33R(-1,kmax,g) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 0,kmax,g) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
        cw33R( 1,kmax,g) =  0.
     END IF
     
     IF (BC_3L == -2) THEN
        cw33R(-1,1   ,g) =  0.
        cw33R( 1,1   ,g) =  0.
     END IF
     IF (BC_3U == -2) THEN
        cw33R(-1,kmax,g) =  0.
        cw33R( 1,kmax,g) =  0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO m = 1, n_conc
        cc33R(:,:,g,m) = cp33R(:,:,g)
        
        IF (BCc_3L(m) > 0) THEN
           cc33R(-1,1   ,g,m) =  0.
           cc33R( 0,1   ,g,m) = -1./(x3pR(2   ,g) - x3pR(1     ,g))
           cc33R( 1,1   ,g,m) =  1./(x3pR(2   ,g) - x3pR(1     ,g))
        END IF
        IF (BCc_3U(m) > 0) THEN
           cc33R(-1,kmax,g,m) = -1./(x3pR(kmax,g) - x3pR(kmax-1,g))
           cc33R( 0,kmax,g,m) =  1./(x3pR(kmax,g) - x3pR(kmax-1,g))
           cc33R( 1,kmax,g,m) =  0.
        END IF
     END DO
     
     END IF
     !========================================================================================================
     
  END DO
  
  g = 1
  
  !===========================================================================================================
  !=== feinstes Gitter (Tangential-Koeffizienten überschreiben) ==============================================
  !===========================================================================================================
  CALL diff_coeffs(1,2,0,mapping_yes,2,(/2,3/),x1u(-1:N1+1),x1u(-1:N1+1),BC_1L,BC_1U,1,N1,-1,1,cu11R(-1,0,g))
  
  IF (BC_1L == 1 .OR. BC_1L == 2) THEN
     cu11R(-1,0 ,g) = 0.
     cu11R( 0,0 ,g) = 1.- (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =     (x1p(1 )-x1u(0   )) / (x1u(1 )-x1u(0   ))
  END IF
  IF (BC_1U == 1 .OR. BC_1U == 2) THEN
     cu11R(-1,N1,g) = 1.- (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =     (x1p(N1)-x1u(N1-1)) / (x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) = 0.
  END IF
  
  IF (BC_1L == 3 .OR. BC_1L == 4) THEN
     cu11R(-1,0 ,g) =  0.
     cu11R( 0,0 ,g) = -1./(x1u(1 )-x1u(0   ))
     cu11R( 1,0 ,g) =  1./(x1u(1 )-x1u(0   ))
  END IF
  IF (BC_1U == 3 .OR. BC_1U == 4) THEN
     cu11R(-1,N1,g) = -1./(x1u(N1)-x1u(N1-1))
     cu11R( 0,N1,g) =  1./(x1u(N1)-x1u(N1-1))
     cu11R( 1,N1,g) =  0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  CALL diff_coeffs(2,2,0,mapping_yes,2,(/2,3/),x2v(-1:N2+1),x2v(-1:N2+1),BC_2L,BC_2U,1,N2,-1,1,cv22R(-1,0,g))
  
  IF (BC_2L == 1 .OR. BC_2L == 2) THEN
     cv22R(-1,0 ,g) = 0.
     cv22R( 0,0 ,g) = 1.- (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =     (x2p(1 )-x2v(0   )) / (x2v(1 )-x2v(0   ))
  END IF
  IF (BC_2U == 1 .OR. BC_2U == 2) THEN
     cv22R(-1,N2,g) = 1.- (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =     (x2p(N2)-x2v(N2-1)) / (x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) = 0.
  END IF
  
  IF (BC_2L == 3 .OR. BC_2L == 4) THEN
     cv22R(-1,0 ,g) =  0.
     cv22R( 0,0 ,g) = -1./(x2v(1 )-x2v(0   ))
     cv22R( 1,0 ,g) =  1./(x2v(1 )-x2v(0   ))
  END IF
  IF (BC_2U == 3 .OR. BC_2U == 4) THEN
     cv22R(-1,N2,g) = -1./(x2v(N2)-x2v(N2-1))
     cv22R( 0,N2,g) =  1./(x2v(N2)-x2v(N2-1))
     cv22R( 1,N2,g) =  0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  CALL diff_coeffs(3,2,0,mapping_yes,2,(/2,3/),x3w(-1:N3+1),x3w(-1:N3+1),BC_3L,BC_3U,1,N3,-1,1,cw33R(-1,0,g))
  
  IF (BC_3L == 1 .OR. BC_3L == 2) THEN
     cw33R(-1,0 ,g) = 0.
     cw33R( 0,0 ,g) = 1.- (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =     (x3p(1 )-x3w(0   )) / (x3w(1 )-x3w(0   ))
  END IF
  IF (BC_3U == 1 .OR. BC_3U == 2) THEN
     cw33R(-1,N3,g) = 1.- (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =     (x3p(N3)-x3w(N3-1)) / (x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) = 0.
  END IF
  
  IF (BC_3L == 3 .OR. BC_3L == 4) THEN
     cw33R(-1,0 ,g) =  0.
     cw33R( 0,0 ,g) = -1./(x3w(1 )-x3w(0   ))
     cw33R( 1,0 ,g) =  1./(x3w(1 )-x3w(0   ))
  END IF
  IF (BC_3U == 3 .OR. BC_3U == 4) THEN
     cw33R(-1,N3,g) = -1./(x3w(N3)-x3w(N3-1))
     cw33R( 0,N3,g) =  1./(x3w(N3)-x3w(N3-1))
     cw33R( 1,N3,g) =  0.
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE get_stencil_Helm
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE interp_coeffs
  
  IMPLICIT NONE
  
  INTEGER               ::  i, i0, di
  INTEGER               ::  j, j0, dj
  INTEGER               ::  k, k0, dk
  INTEGER               ::  g
  
  REAL                  ::  Dx12, Dx10
  
  
  cI1 = 0.
  cI2 = 0.
  cI3 = 0.
  
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  DO g = 1, n_grids-1
     
     di = (M1-1)/((NN(1,g)-1)*NB(1,g))
     dj = (M2-1)/((NN(2,g)-1)*NB(2,g))
     dk = (M3-1)/((NN(3,g)-1)*NB(3,g))
     
     !--------------------------------------------------------------------------------------------------------
     DO i = 2, NN(1,g)-1, 2
        i0 = 1 + (i-1)*di + (M1-1)*(iB(1,g)-1)/NB(1,g)
        
        Dx10 = y1p(i0   )-y1p(i0-di)
        Dx12 = y1p(i0+di)-y1p(i0-di)
        
        cI1(1,i,g) = 1.- Dx10/Dx12
        cI1(2,i,g) =     Dx10/Dx12
     END DO
     !--------------------------------------------------------------------------------------------------------
     DO j = 2, NN(2,g)-1, 2
        j0 = 1 + (j-1)*dj + (M2-1)*(iB(2,g)-1)/NB(2,g)
        
        Dx10 = y2p(j0   )-y2p(j0-dj)
        Dx12 = y2p(j0+dj)-y2p(j0-dj)
        
        cI2(1,j,g) = 1.- Dx10/Dx12
        cI2(2,j,g) =     Dx10/Dx12
     END DO
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     DO k = 2, NN(3,g)-1, 2
        k0 = 1 + (k-1)*dk + (M3-1)*(iB(3,g)-1)/NB(3,g)
        
        Dx10 = y3p(k0   )-y3p(k0-dk)
        Dx12 = y3p(k0+dk)-y3p(k0-dk)
        
        cI3(1,k,g) = 1.- Dx10/Dx12
        cI3(2,k,g) =     Dx10/Dx12
     END DO
     END IF
     !--------------------------------------------------------------------------------------------------------
  END DO
  !===========================================================================================================
  
  
  END SUBROUTINE interp_coeffs
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE interp_coeffs_Helm
  
  IMPLICIT NONE
  
  INTEGER               ::  i, j, k
  REAL                  ::  Dx12, Dx1a, Dx1b
  
  
  cIH1 = 0.
  cIH2 = 0.
  cIH3 = 0.
  
  !===========================================================================================================
  !=== Interpolation, linienweise, 1d ========================================================================
  !===========================================================================================================
  DO i = 1, N1-2, 2
     Dx1a = x1u(i  )-x1p(i)
     Dx1b = x1u(i+1)-x1p(i)
     Dx12 = x1p(i+2)-x1p(i)
     
     cIH1(1,i  ) = 1.- Dx1a/Dx12
     cIH1(2,i  ) =     Dx1a/Dx12
     
     cIH1(1,i+1) = 1.- Dx1b/Dx12
     cIH1(2,i+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_1L > 0) THEN
     Dx1a = x1u(0)-x1p(1)
     Dx12 = x1p(3)-x1p(1)
     
     cIH1(1,0) = 1.- Dx1a/Dx12
     cIH1(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_1U > 0) THEN
     Dx1a = x1u(N1)-x1p(N1-2)
     Dx12 = x1p(N1)-x1p(N1-2)
     
     cIH1(1,N1) = 1.- Dx1a/Dx12
     cIH1(2,N1) =     Dx1a/Dx12
  END IF
  !-----------------------------------------------------------------------------------------------------------
  DO j = 1, N2-2, 2
     Dx1a = x2v(j  )-x2p(j)
     Dx1b = x2v(j+1)-x2p(j)
     Dx12 = x2p(j+2)-x2p(j)
     
     cIH2(1,j  ) = 1.- Dx1a/Dx12
     cIH2(2,j  ) =     Dx1a/Dx12
     
     cIH2(1,j+1) = 1.- Dx1b/Dx12
     cIH2(2,j+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_2L > 0) THEN
     Dx1a = x2v(0)-x2p(1)
     Dx12 = x2p(3)-x2p(1)
     
     cIH2(1,0) = 1.- Dx1a/Dx12
     cIH2(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_2U > 0) THEN
     Dx1a = x2v(N2)-x2p(N2-2)
     Dx12 = x2p(N2)-x2p(N2-2)
     
     cIH2(1,N2) = 1.- Dx1a/Dx12
     cIH2(2,N2) =     Dx1a/Dx12
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  DO k = 1, N3-2, 2
     Dx1a = x3w(k  )-x3p(k)
     Dx1b = x3w(k+1)-x3p(k)
     Dx12 = x3p(k+2)-x3p(k)
     
     cIH3(1,k  ) = 1.- Dx1a/Dx12
     cIH3(2,k  ) =     Dx1a/Dx12
     
     cIH3(1,k+1) = 1.- Dx1b/Dx12
     cIH3(2,k+1) =     Dx1b/Dx12
  END DO
  
  !--- Randbedingungen (Extrapolation) ---
  IF (BC_3L > 0) THEN
     Dx1a = x3w(0)-x3p(1)
     Dx12 = x3p(3)-x3p(1)
     
     cIH3(1,0) = 1.- Dx1a/Dx12
     cIH3(2,0) =     Dx1a/Dx12
  END IF
  IF (BC_3U > 0) THEN
     Dx1a = x3w(N3)-x3p(N3-2)
     Dx12 = x3p(N3)-x3p(N3-2)
     
     cIH3(1,N3) = 1.- Dx1a/Dx12
     cIH3(2,N3) =     Dx1a/Dx12
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE interp_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE restr_coeffs ! TEST!!! aufraeumen und Variablen substituieren ...
  
  IMPLICIT NONE
  
  INTEGER               ::  g
  INTEGER               ::  i, ii, iimax
  INTEGER               ::  j, jj, jjmax
  INTEGER               ::  k, kk, kkmax
  
  INTEGER               ::  BC_1L, BC_1U, BC_2L, BC_2U, BC_3L, BC_3U ! TEST!!!
  
  
  cR1 = 0.
  cR2 = 0.
  cR3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  DO g = 2, n_grids
     
     iimax = (NN(1,g)-1)/n_gather(1,g)+1
     jjmax = (NN(2,g)-1)/n_gather(2,g)+1
     kkmax = (NN(3,g)-1)/n_gather(3,g)+1
     
     !iimax = (NN(1,g)-1)*NB(1,g)/NB(1,g-1)+1 ! TEST!!! alt ...
     !jjmax = (NN(2,g)-1)*NB(2,g)/NB(2,g-1)+1
     !kkmax = (NN(3,g)-1)*NB(3,g)/NB(3,g-1)+1
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(1,g) > 1) THEN
        BC_1L = BC(1,1,g-1) ! TEST!!! evtl. wieder auf dem feinen Gitter speichern ...
        BC_1U = BC(2,1,g-1)
     ELSE
        BC_1L = BC(1,1,g)
        BC_1U = BC(2,1,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(2,g) > 1) THEN
        BC_2L = BC(1,2,g-1)
        BC_2U = BC(2,2,g-1)
     ELSE
        BC_2L = BC(1,2,g)
        BC_2U = BC(2,2,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (n_gather(3,g) > 1) THEN
        BC_3L = BC(1,3,g-1)
        BC_3U = BC(2,3,g-1)
     ELSE
        BC_3L = BC(1,3,g)
        BC_3U = BC(2,3,g)
     END IF
     !--------------------------------------------------------------------------------------------------------
     !--------------------------------------------------------------------------------------------------------
     DO ii = 1, iimax
        cR1(-1,ii,g) = 1./4.
        cR1( 0,ii,g) = 2./4.
        cR1( 1,ii,g) = 1./4.
     END DO
     
     IF (BC_1L > 0) THEN
        cR1(-1,1,g) = 0.
        cR1( 0,1,g) = 1.
        cR1( 1,1,g) = 0.
        
        cR1(-1,2,g) = 0. ! TEST!!! Sollte evtl. noch ergaenzt werden ...
        cR1( 0,2,g) = 1.
        cR1( 1,2,g) = 0.
     END IF
     IF (BC_1L == -2) THEN
        cR1( 1,1,g) = cR1( 1,1,g) + cR1(-1,1,g)
        cR1(-1,1,g) = 0.
     END IF
     
     IF (BC_1U > 0) THEN
        cR1(-1,iimax  ,g) = 0.
        cR1( 0,iimax  ,g) = 1.
        cR1( 1,iimax  ,g) = 0.
        
        cR1(-1,iimax-1,g) = 0. ! TEST!!!
        cR1( 0,iimax-1,g) = 1.
        cR1( 1,iimax-1,g) = 0.
     END IF
     IF (BC_1U == -2) THEN
        cR1(-1,iimax,g) = cR1( 1,iimax,g) + cR1(-1,iimax,g)
        cR1( 1,iimax,g) = 0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     DO jj = 1, jjmax
        cR2(-1,jj,g) = 1./4.
        cR2( 0,jj,g) = 2./4.
        cR2( 1,jj,g) = 1./4.
     END DO
     
     IF (BC_2L > 0) THEN
        cR2(-1,1,g) = 0.
        cR2( 0,1,g) = 1.
        cR2( 1,1,g) = 0.
        
        cR2(-1,2,g) = 0. ! TEST!!!
        cR2( 0,2,g) = 1.
        cR2( 1,2,g) = 0.
     END IF
     IF (BC_2L == -2) THEN
        cR2( 1,1,g) = cR2( 1,1,g) + cR2(-1,1,g)
        cR2(-1,1,g) = 0.
     END IF
     
     IF (BC_2U > 0) THEN
        cR2(-1,jjmax  ,g) = 0.
        cR2( 0,jjmax  ,g) = 1.
        cR2( 1,jjmax  ,g) = 0.
        
        cR2(-1,jjmax-1,g) = 0. ! TEST!!!
        cR2( 0,jjmax-1,g) = 1.
        cR2( 1,jjmax-1,g) = 0.
     END IF
     IF (BC_2U == -2) THEN
        cR2(-1,jjmax,g) = cR2( 1,jjmax,g) + cR2(-1,jjmax,g)
        cR2( 1,jjmax,g) = 0.
     END IF
     !--------------------------------------------------------------------------------------------------------
     IF (dimens == 3) THEN
     DO kk = 1, kkmax
        cR3(-1,kk,g) = 1./4.
        cR3( 0,kk,g) = 2./4.
        cR3( 1,kk,g) = 1./4.
     END DO
     
     IF (BC_3L > 0) THEN
        cR3(-1,1,g) = 0.
        cR3( 0,1,g) = 1.
        cR3( 1,1,g) = 0.
        
        cR3(-1,2,g) = 0. ! TEST!!!
        cR3( 0,2,g) = 1.
        cR3( 1,2,g) = 0.
     END IF
     IF (BC_3L == -2) THEN
        cR3( 1,1,g) = cR3( 1,1,g) + cR3(-1,1,g)
        cR3(-1,1,g) = 0.
     END IF
     
     IF (BC_3U > 0) THEN
        cR3(-1,kkmax  ,g) = 0.
        cR3( 0,kkmax  ,g) = 1.
        cR3( 1,kkmax  ,g) = 0.
        
        cR3(-1,kkmax-1,g) = 0. ! TEST!!!
        cR3( 0,kkmax-1,g) = 1.
        cR3( 1,kkmax-1,g) = 0.
     END IF
     IF (BC_3U == -2) THEN
        cR3(-1,kkmax,g) = cR3( 1,kkmax,g) + cR3(-1,kkmax,g)
        cR3( 1,kkmax,g) = 0.
     END IF
     
     END IF
     !--------------------------------------------------------------------------------------------------------
  END DO
  !===========================================================================================================
  
  
  IF (1 == 2) THEN ! TEST!!!
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  cRest1 = 0.
  cRest2 = 0.
  cRest3 = 0.
  
  DO g = 1, n_grids-1
     iimax = NN(1,g)
     jjmax = NN(2,g)
     kkmax = NN(3,g)
     
     IF (ncb1f(dim_ncb1c) < iimax) THEN ! TEST!!! Quick and dirty ...
        CALL diff_coeffs(1,0,0,.FALSE.,dim_ncb1c,ncb1r,x1pR(b1L,g),x1pR(b1L,g),BC_1L,BC_1U,0,iimax,b1L,b1U,cRest1(b1L,0,g))
     ELSE
        DO i = 1, iimax
           cRest1(-1:1,i,g) = (/0.25,0.5,0.25/) ! TEST!!! angepasste Koeffizienten? (mapping_yes = .FALSE. waere besser ...)
        END DO
        IF (BC_1L > 0) cRest1(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_1L ==  -2) cRest1(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_1U > 0) cRest1(-1:1,iimax,g) = (/0. ,1. ,0. /)
        IF (BC_1U ==  -2) cRest1(-1:1,iimax,g) = (/0.5,0.5,0. /)
     END IF
     
     
     IF (ncb2f(dim_ncb2c) < jjmax) THEN
        CALL diff_coeffs(2,0,0,.FALSE.,dim_ncb2c,ncb2r,x2pR(b2L,g),x2pR(b2L,g),BC_2L,BC_2U,0,jjmax,b2L,b2U,cRest2(b2L,0,g))
     ELSE
        DO j = 1, jjmax
           cRest2(-1:1,j,g) = (/0.25,0.5,0.25/)
        END DO
        IF (BC_2L > 0) cRest2(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_2L ==  -2) cRest2(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_2U > 0) cRest2(-1:1,jjmax,g) = (/0. ,1. ,0. /)
        IF (BC_2U ==  -2) cRest2(-1:1,jjmax,g) = (/0.5,0.5,0. /)
     END IF
     
     
     IF (dimens == 3) THEN
     IF (ncb3f(dim_ncb3c) < kkmax) THEN
        CALL diff_coeffs(3,0,0,.FALSE.,dim_ncb3c,ncb3r,x3pR(b3L,g),x3pR(b3L,g),BC_3L,BC_3U,0,kkmax,b3L,b3U,cRest3(b3L,0,g))
     ELSE
        DO k = 1, kkmax
           cRest3(-1:1,k,g) = (/0.25,0.5,0.25/)
        END DO
        IF (BC_3L > 0) cRest3(-1:1,1    ,g) = (/0. ,1. ,0. /)
        IF (BC_3L ==  -2) cRest3(-1:1,1    ,g) = (/0. ,0.5,0.5/)
        IF (BC_3U > 0) cRest3(-1:1,kkmax,g) = (/0. ,1. ,0. /)
        IF (BC_3U ==  -2) cRest3(-1:1,kkmax,g) = (/0.5,0.5,0. /)
     END IF
     END IF
     
  END DO
  !===========================================================================================================
  END IF
  
  
  END SUBROUTINE restr_coeffs
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE restr_coeffs_Helm
  
  IMPLICIT NONE
  
  INTEGER               ::  i, j, k
  REAL                  ::  Dx12, Dx1a
  
  
  cRH1 = 0.
  cRH2 = 0.
  cRH3 = 0.
  
  !===========================================================================================================
  !=== Restriktion, linienweise, 1d ==========================================================================
  !===========================================================================================================
  DO i = 1, N1, 2
     Dx1a = x1p(i)-x1u(i-1)
     Dx12 = x1u(i)-x1u(i-1)
     
     cRH1(1,i) = 1.- Dx1a/Dx12
     cRH1(2,i) =     Dx1a/Dx12
  END DO
  
  IF (BC_1L > 0) THEN
     cRH1(1,1 ) = 1.
     cRH1(2,1 ) = 0.
  END IF
  IF (BC_1L == -2) THEN
     cRH1(:,1 ) = 0.
  END IF
  
  IF (BC_1U > 0) THEN
     cRH1(1,N1) = 0.
     cRH1(2,N1) = 1.
  END IF
  IF (BC_1U == -2) THEN
     cRH1(:,N1) = 0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  DO j = 1, N2, 2
     Dx1a = x2p(j)-x2v(j-1)
     Dx12 = x2v(j)-x2v(j-1)
     
     cRH2(1,j) = 1.- Dx1a/Dx12
     cRH2(2,j) =     Dx1a/Dx12
  END DO
  
  IF (BC_2L > 0) THEN
     cRH2(1,1 ) = 1.
     cRH2(2,1 ) = 0.
  END IF
  IF (BC_2L == -2) THEN
     cRH2(:,1 ) = 1.
  END IF
  
  IF (BC_2U > 0) THEN
     cRH2(1,N2) = 0.
     cRH2(2,N2) = 1.
  END IF
  IF (BC_2U == -2) THEN
     cRH2(:,N2) = 0.
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
  
  DO k = 1, N3, 2
     Dx1a = x3p(k)-x3w(k-1)
     Dx12 = x3w(k)-x3w(k-1)
     
     cRH3(1,k) = 1.- Dx1a/Dx12
     cRH3(2,k) =     Dx1a/Dx12
  END DO
  
  IF (BC_3L > 0) THEN
     cRH3(1,1 ) = 1.
     cRH3(2,1 ) = 0.
  END IF
  IF (BC_3L == -2) THEN
     cRH3(:,1 ) = 1.
  END IF
  
  IF (BC_3U > 0) THEN
     cRH3(1,N3) = 0.
     cRH3(2,N3) = 1.
  END IF
  IF (BC_3U == -2) THEN
     cRH3(:,N3) = 0.
  END IF
  
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE restr_coeffs_Helm
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE get_weights ! TEST!!! bezieht sich bislang nur auf explizite Differenzen ... ODER: rauswerfen!
  
  ! revised: 24.10.07
  
  IMPLICIT NONE
  
  INTEGER                ::  i, ii
  INTEGER                ::  j, jj
  INTEGER                ::  k, kk
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Anmerkungen: - Diese Abschaetzung ist im Feld etwas konservativer (Faktor 1.5-2.5) im Vergleich zu       !
  !                Differenzen 2. Konvergenzordnung. Am Rand dagegen sehr viel strikter!                     !
  !              - dx_iL, dx_iU d�rfen nicht der Ableitung der Mapping-Funktion an der Wand gleichgesetzt    !
  !                werden, da diese im Extremfall auch verschwinden kann.                                    !
  !----------------------------------------------------------------------------------------------------------!
  
  
  
  IF (1 == 1) THEN ! TEST!!! Rand muesste noch extra behandelt werden!
  !===========================================================================================================
  !=== Gewichte f�r �berpr�fung der Divergenzfreiheit ========================================================
  !===========================================================================================================
  weight = 0.
  !-----------------------------------------------------------------------------------------------------------
  IF (dimens == 3) THEN
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              DO ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              END DO
              DO jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              END DO
              DO kk = d3L, d3U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
              END DO
              weight(i,j,k) = 1./weight(i,j,k)
           END DO
        END DO
     END DO
     
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO i = S1p, N1p
              DO ii = d1L, d1U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
              END DO
              DO jj = d2L, d2U
                 weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
              END DO
              weight(i,j,k) = 1./weight(i,j,k)
           END DO
        END DO
     END DO
     
  END IF
  !-----------------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------------
  ELSE
     weight = 1.
  END IF
  
  IF (1 == 1) THEN ! TEST!!!
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_1L > 0) THEN
     i = 1
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_1U > 0) THEN
     i = N1
     DO k = S3p, N3p
        DO j = S2p, N2p
           DO ii = d1L, d1U
              weight(i,j,k) = weight(i,j,k) + ABS(cDu1(ii,i))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_2L > 0) THEN
     j = 1
     DO k = S3p, N3p
        DO i = S1p, N1p
           DO jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_2U > 0) THEN
     j = N2
     DO k = S3p, N3p
        DO i = S1p, N1p
           DO jj = d2L, d2U
              weight(i,j,k) = weight(i,j,k) + ABS(cDv2(jj,j))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  IF (BC_3L > 0) THEN
     k = 1
     DO j = S2p, N2p
        DO i = S1p, N1p
           DO kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  IF (BC_3U > 0) THEN
     k = N3
     DO j = S2p, N2p
        DO i = S1p, N1p
           DO kk = d3L, d3U
              weight(i,j,k) = weight(i,j,k) + ABS(cDw3(kk,k))
           END DO
           weight(i,j,k) = 1./weight(i,j,k)
        END DO
     END DO
  END IF
  !-----------------------------------------------------------------------------------------------------------
  END IF
  
  
  ! Kanten/Ecken sind durch die Randbedingungen immer Divergenz-frei:
  IF (BC_1L > 0 .AND. BC_2L > 0) weight(1 ,1 ,1:N3) = 0.
  IF (BC_1L > 0 .AND. BC_2U > 0) weight(1 ,N2,1:N3) = 0.
  IF (BC_1U > 0 .AND. BC_2L > 0) weight(N1,1 ,1:N3) = 0.
  IF (BC_1U > 0 .AND. BC_2U > 0) weight(N1,N2,1:N3) = 0.
  
  IF (BC_1L > 0 .AND. BC_3L > 0) weight(1 ,1:N2,1 ) = 0.
  IF (BC_1L > 0 .AND. BC_3U > 0) weight(1 ,1:N2,N3) = 0.
  IF (BC_1U > 0 .AND. BC_3L > 0) weight(N1,1:N2,1 ) = 0.
  IF (BC_1U > 0 .AND. BC_3U > 0) weight(N1,1:N2,N3) = 0.
  
  IF (BC_2L > 0 .AND. BC_3L > 0) weight(1:N1,1 ,1 ) = 0.
  IF (BC_2L > 0 .AND. BC_3U > 0) weight(1:N1,1 ,N3) = 0.
  IF (BC_2U > 0 .AND. BC_3L > 0) weight(1:N1,N2,1 ) = 0.
  IF (BC_2U > 0 .AND. BC_3U > 0) weight(1:N1,N2,N3) = 0.
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Relaxation der Konzentration ==========================================================================
  !===========================================================================================================
  IF (concentration_yes) THEN
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
  END IF
  !===========================================================================================================
  
  
  END SUBROUTINE get_weights
  
  
  
  
  
  
  
  
  
  
  
  SUBROUTINE Matrix_invert(N,matrix,matrix_inv)
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  ::  N
  REAL   , INTENT(IN)  ::  matrix     (1:N,1:N)
  REAL   , INTENT(OUT) ::  matrix_inv (1:N,1:N)
  
  REAL                 ::  matrix_left(1:N,1:N)
  REAL                 ::  mult1, mult2
  REAL                 ::  eps
  INTEGER              ::  i, j, k
  
  REAL                 ::  store
  INTEGER              ::  maxValue, maxValuePos
  
  
  
  eps = 10.**(-20) ! double precision
  !eps = 10.**(-??) ! single precision
  
  matrix_left = matrix
  
  
  ! Rechte Seite initialisieren (= Einheitsmatrix):
  matrix_inv = 0.
  
  DO i = 1, N
     matrix_inv(i,i) = 1.
  END DO  
  
  
  !--- Vorwaertsschritt -------------------------------------
  ! linke Seite umformen in obere Dreiecksmatrix
  ! rechte Seite umformen in untere Dreiecksmatrix
  
  DO j = 1, N-1
     
     ! Pivoting == Umsortieren der aktuellen Untermatrix (j:N,j:N)
     ! (Diagonalelemente der linken Dreiecksmatrix sind betragsmaessig zu maximieren)
     ! Groesster Wert in Spalte j = zukuenftiges Diagonalelement j,j
     maxValue    = ABS(matrix_left(j,j))
     maxValuePos = j
     DO i = j+1, N
        IF (ABS(matrix_left(i,j)) > maxValue) THEN
           maxValue    = ABS(matrix_left(i,j))
           maxValuePos = i
        END IF
     END DO
     
     ! Zeilen vertauschen:
     IF (maxValuePos /= j) THEN
        DO i = 1, N
           store                      = matrix_left(maxValuePos,i)
           matrix_left(maxValuePos,i) = matrix_left(j,i)
           matrix_left(j,i)           = store
           
           store                      = matrix_inv(maxValuePos,i)
           matrix_inv(maxValuePos,i)  = matrix_inv(j,i)
           matrix_inv(j,i)            = store
        END DO
     END IF
     
     mult1 = matrix_left(j,j)
     IF (ABS(mult1) < eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (1)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', j
     END IF
     
     DO i = j+1, N
        mult2 = matrix_left(i,j)/mult1
                
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (linke Seite):
        DO k = j, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) < eps .AND. ABS(matrix_left(j,k)) < eps)) THEN
              matrix_left(i,k) = matrix_left(i,k) - matrix_left(j,k)*mult2
           END IF
        END DO
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        DO k = 1, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) < eps .AND. ABS(matrix_inv(j,k)) < eps)) THEN
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           END IF
        END DO
        
        ! Komponente i,j explizit zu Null setzen:
        matrix_left(i,j) = 0.
     END DO
     
  END DO
  
  
  !--- Rueckwaertsschritt -----------------------------------
  ! linke Seite umformen in Einheitsmatrix
  ! rechte Seite umformen in gesuchte inverse Matrix
  
  DO j = N, 2, -1
     
     ! Multiplikator:
     mult1 = matrix_left(j,j)
     IF (ABS(mult1) < eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (2)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', j
     END IF
     
     DO i = 1, j-1
        mult2 = matrix_left(i,j)/mult1
        
        ! Aufaddieren von Zeile j auf aktuelle Zeile i (rechte Seite):
        DO k = 1, N
           ! To avoid underflow:
           IF (.NOT. (ABS(mult2) < eps .AND. ABS(matrix_inv(j,k)) < eps)) THEN
              matrix_inv(i,k) = matrix_inv(i,k) - matrix_inv(j,k)*mult2
           END IF
        END DO
        
        ! nicht-Diagonalelement explizit zu 0 setzen:
        matrix_left(i,j) = 0.
     END DO
     
  END DO
  
  
  ! linke Matrix auf Einheitsmatrix bringen:
  DO i = 1, N
     mult1 = matrix_left(i,i)
     IF (ABS(mult1) < eps .AND. rank == 0) THEN
        WRITE(*,'(a)'      ) 'WARNING! Matrix is probably singular (3)!'
        WRITE(*,'(a,E13.5)') '       ... dividing by', mult1
        WRITE(*,'(a,i2)'   ) '                   i =', i
     END IF
     matrix_inv(i,:) = matrix_inv(i,:)/mult1
  END DO
  
  
  END SUBROUTINE Matrix_invert
  
  
  
  
END MODULE mod_coeffs
