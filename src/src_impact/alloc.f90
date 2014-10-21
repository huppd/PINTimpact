!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************
  
!#ifdef ALLOC
  !===========================================================================================================
  !=== Differenzen-Koeffizienten-Arrays ======================================================================
  !===========================================================================================================
  !--- 1. Ableitung (zentral) --------------------------------------------------------------------------------
  ALLOCATE(cp1  (b1L:b1U,0:N1))
  ALLOCATE(cp2  (b2L:b2U,0:N2))
  ALLOCATE(cp3  (b3L:b3U,0:N3))
  
  ALLOCATE(cu1  (b1L:b1U,0:N1))
  ALLOCATE(cv2  (b2L:b2U,0:N2))
  ALLOCATE(cw3  (b3L:b3U,0:N3))
  
  ALLOCATE(cc1  (b1L:b1U,0:N1,1:n_conc))
  ALLOCATE(cc2  (b2L:b2U,0:N2,1:n_conc))
  ALLOCATE(cc3  (b3L:b3U,0:N3,1:n_conc))
  
  !--- 1. Ableitung (upwind) ---------------------------------------------------------------------------------
  ALLOCATE(cNp1D(n1L:n1U,0:N1))
  ALLOCATE(cNp2D(n2L:n2U,0:N2))
  ALLOCATE(cNp3D(n3L:n3U,0:N3))
  
  ALLOCATE(cNp1U(n1L:n1U,0:N1))
  ALLOCATE(cNp2U(n2L:n2U,0:N2))
  ALLOCATE(cNp3U(n3L:n3U,0:N3))
  
  ALLOCATE(cNu1D(n1L:n1U,0:N1))
  ALLOCATE(cNv2D(n2L:n2U,0:N2))
  ALLOCATE(cNw3D(n3L:n3U,0:N3))
  
  ALLOCATE(cNu1U(n1L:n1U,0:N1))
  ALLOCATE(cNv2U(n2L:n2U,0:N2))
  ALLOCATE(cNw3U(n3L:n3U,0:N3))
  
  ALLOCATE(cNc1D(n1L:n1U,0:N1,1:n_conc))
  ALLOCATE(cNc2D(n2L:n2U,0:N2,1:n_conc))
  ALLOCATE(cNc3D(n3L:n3U,0:N3,1:n_conc))
  
  ALLOCATE(cNc1U(n1L:n1U,0:N1,1:n_conc))
  ALLOCATE(cNc2U(n2L:n2U,0:N2,1:n_conc))
  ALLOCATE(cNc3U(n3L:n3U,0:N3,1:n_conc))
  
  !--- Divergenz ---------------------------------------------------------------------------------------------
  ALLOCATE(cDu1 (d1L:d1U,0:N1))
  ALLOCATE(cDv2 (d2L:d2U,0:N2))
  ALLOCATE(cDw3 (d3L:d3U,0:N3))
  
  !--- Divergenz (transponiert) ------------------------------------------------------------------------------
  ALLOCATE(cDu1T(g1L:g1U,0:N1))
  ALLOCATE(cDv2T(g2L:g2U,0:N2))
  ALLOCATE(cDw3T(g3L:g3U,0:N3))
  
  !--- Gradient ----------------------------------------------------------------------------------------------
  ALLOCATE(cGp1 (g1L:g1U,0:N1))
  ALLOCATE(cGp2 (g2L:g2U,0:N2))
  ALLOCATE(cGp3 (g3L:g3U,0:N3))
  
  !--- Gradient (transponiert) -------------------------------------------------------------------------------
  ALLOCATE(cGp1T(d1L:d1U,0:N1))
  ALLOCATE(cGp2T(d2L:d2U,0:N2))
  ALLOCATE(cGp3T(d3L:d3U,0:N3))
  
  !--- 2. Ableitung (zentral) --------------------------------------------------------------------------------
  ALLOCATE(cp11 (b1L:b1U,0:N1))
  ALLOCATE(cp22 (b2L:b2U,0:N2))
  ALLOCATE(cp33 (b3L:b3U,0:N3))
  
  ALLOCATE(cu11 (b1L:b1U,0:N1))
  ALLOCATE(cv22 (b2L:b2U,0:N2))
  ALLOCATE(cw33 (b3L:b3U,0:N3))
  
  ALLOCATE(cc11 (b1L:b1U,0:N1,1:n_conc))
  ALLOCATE(cc22 (b2L:b2U,0:N2,1:n_conc))
  ALLOCATE(cc33 (b3L:b3U,0:N3,1:n_conc))
  
  !--- Interpolation ----------------------------------------------------------------------------------------- 
  ALLOCATE(cIpu(g1L:g1U,0:N1))
  ALLOCATE(cIpv(g2L:g2U,0:N2))
  ALLOCATE(cIpw(g3L:g3U,0:N3))
  
  ALLOCATE(cIup(d1L:d1U,0:N1))
  ALLOCATE(cIvp(d2L:d2U,0:N2))
  ALLOCATE(cIwp(d3L:d3U,0:N3))
  
  ALLOCATE(cIcu(g1L:g1U,0:N1,1:n_conc))
  ALLOCATE(cIcv(g2L:g2U,0:N2,1:n_conc))
  ALLOCATE(cIcw(g3L:g3U,0:N3,1:n_conc))
  
  !--- Filter ------------------------------------------------------------------------------------------------
  ALLOCATE(cFp1(b1L:b1U,0:N1))
  ALLOCATE(cFp2(b2L:b2U,0:N2))
  ALLOCATE(cFp3(b3L:b3U,0:N3))
  
  ALLOCATE(cFu1(b1L:b1U,0:N1))
  ALLOCATE(cFv2(b2L:b2U,0:N2))
  ALLOCATE(cFw3(b3L:b3U,0:N3))
  
  ALLOCATE(cFc1(b1L:b1U,0:N1,1:n_conc))
  ALLOCATE(cFc2(b2L:b2U,0:N2,1:n_conc))
  ALLOCATE(cFc3(b3L:b3U,0:N3,1:n_conc))
  
  !--- Integrator (nur für Druckgitter) ----------------------------------------------------------------------
  ALLOCATE(cInt1(b1L:b1U,0:N1))
  ALLOCATE(cInt2(b2L:b2U,0:N2))
  ALLOCATE(cInt3(b3L:b3U,0:N3))
  
  
  !*** Compact *************************************************
  ALLOCATE(buffer1A(0:N2,0:N3,1:2*ndL))
  ALLOCATE(buffer2A(0:N1,0:N3,1:2*ndL))
  ALLOCATE(buffer3A(0:N1,0:N2,1:2*ndL))
  
  ALLOCATE(buffer1B(0:N2,0:N3,1:dimS1))
  ALLOCATE(buffer2B(0:N1,0:N3,1:dimS2))
  ALLOCATE(buffer3B(0:N1,0:N2,1:dimS3))
  
  ALLOCATE(disp1(1:NB1)) ! TEST!!!
  ALLOCATE(disp2(1:NB2))
  ALLOCATE(disp3(1:NB3))
  
  ALLOCATE(recv1(1:NB1)) ! TEST!!!
  ALLOCATE(recv2(1:NB2))
  ALLOCATE(recv3(1:NB3))
  
  !--- Schur complement (square) matrices ---
  ALLOCATE(Sp1(1:dimS1,1:dimS1))
  ALLOCATE(Sp2(1:dimS2,1:dimS2))
  ALLOCATE(Sp3(1:dimS3,1:dimS3))
  
  ALLOCATE(Su1(1:dimS1,1:dimS1))
  ALLOCATE(Sv2(1:dimS2,1:dimS2))
  ALLOCATE(Sw3(1:dimS3,1:dimS3))
  
  ALLOCATE(Sc1(1:dimS1,1:dimS1,1:n_conc))
  ALLOCATE(Sc2(1:dimS2,1:dimS2,1:n_conc))
  ALLOCATE(Sc3(1:dimS3,1:dimS3,1:n_conc))
  
  ALLOCATE(Sp11(1:dimS1,1:dimS1))
  ALLOCATE(Sp22(1:dimS2,1:dimS2))
  ALLOCATE(Sp33(1:dimS3,1:dimS3))
  
  ALLOCATE(Su11(1:dimS1,1:dimS1))
  ALLOCATE(Sv22(1:dimS2,1:dimS2))
  ALLOCATE(Sw33(1:dimS3,1:dimS3))
  
  ALLOCATE(Sc11(1:dimS1,1:dimS1,1:n_conc))
  ALLOCATE(Sc22(1:dimS2,1:dimS2,1:n_conc))
  ALLOCATE(Sc33(1:dimS3,1:dimS3,1:n_conc))
  
  ALLOCATE(SDu1(1:dimS1,1:dimS1))
  ALLOCATE(SDv2(1:dimS2,1:dimS2))
  ALLOCATE(SDw3(1:dimS3,1:dimS3))
  
  ALLOCATE(SDu1T(1:dimS1,1:dimS1))
  ALLOCATE(SDv2T(1:dimS2,1:dimS2))
  ALLOCATE(SDw3T(1:dimS3,1:dimS3))
  
  ALLOCATE(SGp1(1:dimS1,1:dimS1))
  ALLOCATE(SGp2(1:dimS2,1:dimS2))
  ALLOCATE(SGp3(1:dimS3,1:dimS3))
  
  ALLOCATE(SGp1T(1:dimS1,1:dimS1))
  ALLOCATE(SGp2T(1:dimS2,1:dimS2))
  ALLOCATE(SGp3T(1:dimS3,1:dimS3))
  
  ALLOCATE(SIpu(1:dimS1,1:dimS1))
  ALLOCATE(SIpv(1:dimS2,1:dimS2))
  ALLOCATE(SIpw(1:dimS3,1:dimS3))
  
  ALLOCATE(SIup(1:dimS1,1:dimS1))
  ALLOCATE(SIvp(1:dimS2,1:dimS2))
  ALLOCATE(SIwp(1:dimS3,1:dimS3))
  
  ALLOCATE(SIcu(1:dimS1,1:dimS1,1:n_conc))
  ALLOCATE(SIcv(1:dimS2,1:dimS2,1:n_conc))
  ALLOCATE(SIcw(1:dimS3,1:dimS3,1:n_conc))
  
  !--- bottom (horizontally oriented) matrices ---
  ALLOCATE(Wp1(1:2*ndL,0:N1))
  ALLOCATE(Wp2(1:2*ndL,0:N2))
  ALLOCATE(Wp3(1:2*ndL,0:N3))
  
  ALLOCATE(Wu1(1:2*ndL,0:N1))
  ALLOCATE(Wv2(1:2*ndL,0:N2))
  ALLOCATE(Ww3(1:2*ndL,0:N3))
  
  ALLOCATE(Wc1(1:2*ndL,0:N1,1:n_conc))
  ALLOCATE(Wc2(1:2*ndL,0:N2,1:n_conc))
  ALLOCATE(Wc3(1:2*ndL,0:N3,1:n_conc))
  
  ALLOCATE(Wp11(1:2*ndL,0:N1))
  ALLOCATE(Wp22(1:2*ndL,0:N2))
  ALLOCATE(Wp33(1:2*ndL,0:N3))
  
  ALLOCATE(Wu11(1:2*ndL,0:N1))
  ALLOCATE(Wv22(1:2*ndL,0:N2))
  ALLOCATE(Ww33(1:2*ndL,0:N3))
  
  ALLOCATE(Wc11(1:2*ndL,0:N1,1:n_conc))
  ALLOCATE(Wc22(1:2*ndL,0:N2,1:n_conc))
  ALLOCATE(Wc33(1:2*ndL,0:N3,1:n_conc))
  
  ALLOCATE(WDu1(1:2*ndL,0:N1))
  ALLOCATE(WDv2(1:2*ndL,0:N2))
  ALLOCATE(WDw3(1:2*ndL,0:N3))
  
  ALLOCATE(WDu1T(1:2*ndL,0:N1))
  ALLOCATE(WDv2T(1:2*ndL,0:N2))
  ALLOCATE(WDw3T(1:2*ndL,0:N3))
  
  ALLOCATE(WGp1(1:2*ndL,0:N1))
  ALLOCATE(WGp2(1:2*ndL,0:N2))
  ALLOCATE(WGp3(1:2*ndL,0:N3))
  
  ALLOCATE(WGp1T(1:2*ndL,0:N1))
  ALLOCATE(WGp2T(1:2*ndL,0:N2))
  ALLOCATE(WGp3T(1:2*ndL,0:N3))
  
  ALLOCATE(WIpu(1:2*ndL,0:N1))
  ALLOCATE(WIpv(1:2*ndL,0:N2))
  ALLOCATE(WIpw(1:2*ndL,0:N3))
  
  ALLOCATE(WIup(1:2*ndL,0:N1))
  ALLOCATE(WIvp(1:2*ndL,0:N2))
  ALLOCATE(WIwp(1:2*ndL,0:N3))
  
  ALLOCATE(WIcu(1:2*ndL,0:N1,1:n_conc))
  ALLOCATE(WIcv(1:2*ndL,0:N2,1:n_conc))
  ALLOCATE(WIcw(1:2*ndL,0:N3,1:n_conc))
  
  !--- implicit finite difference coefficients ("L"=left-hand side) ---
  ALLOCATE(cp1CL(-ndL:ndL,0:N1))
  ALLOCATE(cp2CL(-ndL:ndL,0:N2))
  ALLOCATE(cp3CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cu1CL(-ndL:ndL,0:N1))
  ALLOCATE(cv2CL(-ndL:ndL,0:N2))
  ALLOCATE(cw3CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cc1CL(-ndL:ndL,0:N1,1:n_conc))
  ALLOCATE(cc2CL(-ndL:ndL,0:N2,1:n_conc))
  ALLOCATE(cc3CL(-ndL:ndL,0:N3,1:n_conc))
  
  ALLOCATE(cp11CL(-ndL:ndL,0:N1))
  ALLOCATE(cp22CL(-ndL:ndL,0:N2))
  ALLOCATE(cp33CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cu11CL(-ndL:ndL,0:N1))
  ALLOCATE(cv22CL(-ndL:ndL,0:N2))
  ALLOCATE(cw33CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cc11CL(-ndL:ndL,0:N1,1:n_conc))
  ALLOCATE(cc22CL(-ndL:ndL,0:N2,1:n_conc))
  ALLOCATE(cc33CL(-ndL:ndL,0:N3,1:n_conc))
  
  ALLOCATE(cDu1CL(-ndL:ndL,0:N1))
  ALLOCATE(cDv2CL(-ndL:ndL,0:N2))
  ALLOCATE(cDw3CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cDu1CLT(-ndL:ndL,0:N1))
  ALLOCATE(cDv2CLT(-ndL:ndL,0:N2))
  ALLOCATE(cDw3CLT(-ndL:ndL,0:N3))
  
  ALLOCATE(cGp1CL(-ndL:ndL,0:N1))
  ALLOCATE(cGp2CL(-ndL:ndL,0:N2))
  ALLOCATE(cGp3CL(-ndL:ndL,0:N3))
  
  ALLOCATE(cGp1CLT(-ndL:ndL,0:N1))
  ALLOCATE(cGp2CLT(-ndL:ndL,0:N2))
  ALLOCATE(cGp3CLT(-ndL:ndL,0:N3))
  
  ALLOCATE(cIpuCL(-ndL:ndL,0:N1))
  ALLOCATE(cIpvCL(-ndL:ndL,0:N2))
  ALLOCATE(cIpwCL(-ndL:ndL,0:N3))
  
  ALLOCATE(cIupCL(-ndL:ndL,0:N1))
  ALLOCATE(cIvpCL(-ndL:ndL,0:N2))
  ALLOCATE(cIwpCL(-ndL:ndL,0:N3))
  
  ALLOCATE(cIcuCL(-ndL:ndL,0:N1,1:n_conc))
  ALLOCATE(cIcvCL(-ndL:ndL,0:N2,1:n_conc))
  ALLOCATE(cIcwCL(-ndL:ndL,0:N3,1:n_conc))
  
  ALLOCATE(cp1CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cp2CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cp3CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cu1CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cv2CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cw3CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cc1CL_LU(1:(3*ndL+1),0:N1,1:n_conc))
  ALLOCATE(cc2CL_LU(1:(3*ndL+1),0:N2,1:n_conc))
  ALLOCATE(cc3CL_LU(1:(3*ndL+1),0:N3,1:n_conc))
  
  ALLOCATE(cp11CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cp22CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cp33CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cu11CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cv22CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cw33CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cc11CL_LU(1:(3*ndL+1),0:N1,1:n_conc))
  ALLOCATE(cc22CL_LU(1:(3*ndL+1),0:N2,1:n_conc))
  ALLOCATE(cc33CL_LU(1:(3*ndL+1),0:N3,1:n_conc))
  
  ALLOCATE(cDu1CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cDv2CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cDw3CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cDu1CLT_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cDv2CLT_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cDw3CLT_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cGp1CL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cGp2CL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cGp3CL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cGp1CLT_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cGp2CLT_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cGp3CLT_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cIpuCL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cIpvCL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cIpwCL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cIupCL_LU(1:(3*ndL+1),0:N1))
  ALLOCATE(cIvpCL_LU(1:(3*ndL+1),0:N2))
  ALLOCATE(cIwpCL_LU(1:(3*ndL+1),0:N3))
  
  ALLOCATE(cIcuCL_LU(1:(3*ndL+1),0:N1,1:n_conc))
  ALLOCATE(cIcvCL_LU(1:(3*ndL+1),0:N2,1:n_conc))
  ALLOCATE(cIcwCL_LU(1:(3*ndL+1),0:N3,1:n_conc))
  
  !--- explicit finite difference coefficients ("R"=right-hand side) ---
  ALLOCATE(cp1CR(-ndR:ndR,0:N1))
  ALLOCATE(cp2CR(-ndR:ndR,0:N2))
  ALLOCATE(cp3CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cu1CR(-ndR:ndR,0:N1))
  ALLOCATE(cv2CR(-ndR:ndR,0:N2))
  ALLOCATE(cw3CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cc1CR(-ndR:ndR,0:N1,1:n_conc))
  ALLOCATE(cc2CR(-ndR:ndR,0:N2,1:n_conc))
  ALLOCATE(cc3CR(-ndR:ndR,0:N3,1:n_conc))
  
  ALLOCATE(cp11CR(-ndR:ndR,0:N1))
  ALLOCATE(cp22CR(-ndR:ndR,0:N2))
  ALLOCATE(cp33CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cu11CR(-ndR:ndR,0:N1))
  ALLOCATE(cv22CR(-ndR:ndR,0:N2))
  ALLOCATE(cw33CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cc11CR(-ndR:ndR,0:N1,1:n_conc))
  ALLOCATE(cc22CR(-ndR:ndR,0:N2,1:n_conc))
  ALLOCATE(cc33CR(-ndR:ndR,0:N3,1:n_conc))
  
  ALLOCATE(cDu1CR(-ndR:ndR,0:N1))
  ALLOCATE(cDv2CR(-ndR:ndR,0:N2))
  ALLOCATE(cDw3CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cDu1CRT(-ndR:ndR,0:N1))
  ALLOCATE(cDv2CRT(-ndR:ndR,0:N2))
  ALLOCATE(cDw3CRT(-ndR:ndR,0:N3))
  
  ALLOCATE(cGp1CR(-ndR:ndR,0:N1))
  ALLOCATE(cGp2CR(-ndR:ndR,0:N2))
  ALLOCATE(cGp3CR(-ndR:ndR,0:N3))
  
  ALLOCATE(cGp1CRT(-ndR:ndR,0:N1))
  ALLOCATE(cGp2CRT(-ndR:ndR,0:N2))
  ALLOCATE(cGp3CRT(-ndR:ndR,0:N3))
  
  ALLOCATE(cIpuCR(-ndR:ndR,0:N1))
  ALLOCATE(cIpvCR(-ndR:ndR,0:N2))
  ALLOCATE(cIpwCR(-ndR:ndR,0:N3))
  
  ALLOCATE(cIupCR(-ndR:ndR,0:N1))
  ALLOCATE(cIvpCR(-ndR:ndR,0:N2))
  ALLOCATE(cIwpCR(-ndR:ndR,0:N3))
  
  ALLOCATE(cIcuCR(-ndR:ndR,0:N1,1:n_conc))
  ALLOCATE(cIcvCR(-ndR:ndR,0:N2,1:n_conc))
  ALLOCATE(cIcwCR(-ndR:ndR,0:N3,1:n_conc))
  !****************************************************
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ALLOCATE(cp11R(-1:1,0:N1,1:n_grids_max))
  ALLOCATE(cp22R(-1:1,0:N2,1:n_grids_max))
  ALLOCATE(cp33R(-1:1,0:N3,1:n_grids_max))
  
  ALLOCATE(cu11R(-1:1,0:N1,1:n_grids_max))
  ALLOCATE(cv22R(-1:1,0:N2,1:n_grids_max))
  ALLOCATE(cw33R(-1:1,0:N3,1:n_grids_max))
  
  ALLOCATE(cc11R(-1:1,0:N1,1:n_grids_max,1:n_conc))
  ALLOCATE(cc22R(-1:1,0:N2,1:n_grids_max,1:n_conc))
  ALLOCATE(cc33R(-1:1,0:N3,1:n_grids_max,1:n_conc))
  
  ALLOCATE(cdg1 (-1:1,1:N1,1:n_grids_max))
  ALLOCATE(cdg2 (-1:1,1:N2,1:n_grids_max))
  ALLOCATE(cdg3 (-1:1,1:N3,1:n_grids_max))
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  ALLOCATE(cI1(1:2,1:N1,1:n_grids_max))
  ALLOCATE(cI2(1:2,1:N2,1:n_grids_max))
  ALLOCATE(cI3(1:2,1:N3,1:n_grids_max))
  
  ALLOCATE(cIH1(1:2,0:N1))
  ALLOCATE(cIH2(1:2,0:N2))
  ALLOCATE(cIH3(1:2,0:N3))
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  ALLOCATE(cR1 (-1:1,1:N1,2:n_grids_max))
  ALLOCATE(cR2 (-1:1,1:N2,2:n_grids_max))
  ALLOCATE(cR3 (-1:1,1:N3,2:n_grids_max))
  
  ALLOCATE(cRest1(b1L:b1U,0:N1,1:n_grids_max-1)) ! TEST!!!
  ALLOCATE(cRest2(b2L:b2U,0:N2,1:n_grids_max-1))
  ALLOCATE(cRest3(b3L:b3U,0:N3,1:n_grids_max-1))
  
  ALLOCATE(cRH1(1:2,1:N1))
  ALLOCATE(cRH2(1:2,1:N2))
  ALLOCATE(cRH3(1:2,1:N3))
  
  
  !===========================================================================================================
  !=== Gitterspezifikationen =================================================================================
  !===========================================================================================================
  !--- physiklische Koordinaten (global) ---------------------------------------------------------------------
  ALLOCATE(y1p(1:M1))
  ALLOCATE(y2p(1:M2))
  ALLOCATE(y3p(1:M3))
  
  ALLOCATE(y1u(0:M1))
  ALLOCATE(y2v(0:M2))
  ALLOCATE(y3w(0:M3))
  
  !--- physiklische Koordinaten (Block) ----------------------------------------------------------------------
  ALLOCATE(x1p(b1L:(N1+b1U)))
  ALLOCATE(x2p(b2L:(N2+b2U)))
  ALLOCATE(x3p(b3L:(N3+b3U)))
  
  ALLOCATE(x1u(b1L:(N1+b1U)))
  ALLOCATE(x2v(b2L:(N2+b2U)))
  ALLOCATE(x3w(b3L:(N3+b3U)))
  
  !--- physiklische Koordinaten (Block, Multigrid) -----------------------------------------------------------
  ALLOCATE(x1pR(b1L:(N1+b1U),1:n_grids_max))
  ALLOCATE(x2pR(b2L:(N2+b2U),1:n_grids_max))
  ALLOCATE(x3pR(b3L:(N3+b3U),1:n_grids_max))
  
  ALLOCATE(x1uR(b1L:(N1+b1U),1:n_grids_max))
  ALLOCATE(x2vR(b2L:(N2+b2U),1:n_grids_max))
  ALLOCATE(x3wR(b3L:(N3+b3U),1:n_grids_max))
  
  !--- Gitterweiten (global) ---------------------------------------------------------------------------------
  ALLOCATE(dy1p(1:M1))
  ALLOCATE(dy2p(1:M2))
  ALLOCATE(dy3p(1:M3))
  
  ALLOCATE(dy1u(0:M1))
  ALLOCATE(dy2v(0:M2))
  ALLOCATE(dy3w(0:M3))
  
  !--- Gitterweiten (Block) ----------------------------------------------------------------------------------
  ALLOCATE(dx1p(1:N1))
  ALLOCATE(dx2p(1:N2))
  ALLOCATE(dx3p(1:N3))
  
  ALLOCATE(dx1u(0:N1))
  ALLOCATE(dx2v(0:N2))
  ALLOCATE(dx3w(0:N3))
  
  ALLOCATE(dx1DM(1:N1))
  ALLOCATE(dx2DM(1:N2))
  ALLOCATE(dx3DM(1:N3))
  
  ALLOCATE(dx1pM(1:N1))
  ALLOCATE(dx2pM(1:N2))
  ALLOCATE(dx3pM(1:N3))
  
  ALLOCATE(ddx1pM(1:N1))
  ALLOCATE(ddx2pM(1:N2))
  ALLOCATE(ddx3pM(1:N3))
  
  ALLOCATE(dx1GM(0:N1))
  ALLOCATE(dx2GM(0:N2))
  ALLOCATE(dx3GM(0:N3))
  
  ALLOCATE(dx1uM(0:N1))
  ALLOCATE(dx2vM(0:N2))
  ALLOCATE(dx3wM(0:N3))
  
  ALLOCATE(ddx1uM(0:N1))
  ALLOCATE(ddx2vM(0:N2))
  ALLOCATE(ddx3wM(0:N3))
  
  !--- Smagorinsky-Modell ---
  ALLOCATE(dx1pS(1:N1))
  ALLOCATE(dx2pS(1:N2))
  ALLOCATE(dx3pS(1:N3))
  
  ALLOCATE(dx1uS(0:N1))
  ALLOCATE(dx2vS(0:N2))
  ALLOCATE(dx3wS(0:N3))
  
  
  !===========================================================================================================
  !=== Arbeitsfelder =========================================================================================
  !===========================================================================================================
  !--- Geschwindigkeiten -------------------------------------------------------------------------------------
  ALLOCATE(vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  
  !--- nicht-linearer Term -----------------------------------------------------------------------------------
  ALLOCATE(nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  
  !--- Recht-Hand-Seite --------------------------------------------------------------------------------------
  ALLOCATE(rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  
  !--- Druck -------------------------------------------------------------------------------------------------
  ALLOCATE(pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Konzentrationsfelder ----------------------------------------------------------------------------------
  ALLOCATE(conc(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc))
  
  !--- nicht-linearer Term (Konzentration) -------------------------------------------------------------------
  ALLOCATE(nlco(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc))
  
  !--- Ausfluss-RB (Geschwindigkeitsfeld) --------------------------------------------------------------------
  ALLOCATE(bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2))
  ALLOCATE(bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2))
  ALLOCATE(bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2))
  
  ALLOCATE(bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2))
  ALLOCATE(bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2))
  ALLOCATE(bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2))
  
  ALLOCATE(bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2))
  ALLOCATE(bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2))
  ALLOCATE(bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2))
  
  ALLOCATE(drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2))
  ALLOCATE(drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2))
  ALLOCATE(drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2))
  
  
  !--- Sedimentations-RB (Konzentration) ---------------------------------------------------------------------
  ALLOCATE(sed1L (1:N2,1:N3,1:n_conc))
  ALLOCATE(sed1U (1:N2,1:N3,1:n_conc))
  
  ALLOCATE(sed2L (1:N1,1:N3,1:n_conc))
  ALLOCATE(sed2U (1:N1,1:N3,1:n_conc))
  
  ALLOCATE(sed3L (1:N1,1:N2,1:n_conc))
  ALLOCATE(sed3U (1:N1,1:N2,1:n_conc))
  
  ALLOCATE(conc1L(1:N2,1:N3,1:n_conc))
  ALLOCATE(conc1U(1:N2,1:N3,1:n_conc))
  
  ALLOCATE(conc2L(1:N1,1:N3,1:n_conc))
  ALLOCATE(conc2U(1:N1,1:N3,1:n_conc))
  
  ALLOCATE(conc3L(1:N1,1:N2,1:n_conc))
  ALLOCATE(conc3U(1:N1,1:N2,1:n_conc))
  
  !--- Sedimentation (Konzentration, Ausschrieb) -------------------------------------------------------------
  ALLOCATE(dep1L_conc(1:N2,1:N3,1:n_conc))
  ALLOCATE(dep1U_conc(1:N2,1:N3,1:n_conc))
  
  ALLOCATE(dep2L_conc(1:N1,1:N3,1:n_conc))
  ALLOCATE(dep2U_conc(1:N1,1:N3,1:n_conc))
  
  ALLOCATE(dep3L_conc(1:N1,1:N2,1:n_conc))
  ALLOCATE(dep3U_conc(1:N1,1:N2,1:n_conc))
  
  !--- Sedimentation (Partikel, Ausschrieb) ------------------------------------------------------------------
  ALLOCATE(dep1L_part(1:N2,1:N3))
  ALLOCATE(dep1U_part(1:N2,1:N3))
  
  ALLOCATE(dep2L_part(1:N1,1:N3))
  ALLOCATE(dep2U_part(1:N1,1:N3))
  
  ALLOCATE(dep3L_part(1:N1,1:N2))
  ALLOCATE(dep3U_part(1:N1,1:N2))
  
  !--- Residuum ----------------------------------------------------------------------------------------------
  ALLOCATE(res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Druckgradient (eine Komponente) -----------------------------------------------------------------------
  ALLOCATE(gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Gewichte für Divergenzfreiheit ------------------------------------------------------------------------
  ALLOCATE(weight(1:N1,1:N2,1:N3))
  
  !--- fluid density on velocity grid ------------------------------------------------------------------------
!#ifdef NONBOUSSINESQ
!  ALLOCATE(dens(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)) ! Checkt der Compiler nicht, wurde herausgezogen ...
!#endif
  
  !--- Null-Raum-Vektor --------------------------------------------------------------------------------------
  ALLOCATE(psi     (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(psi_vel (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))
  ALLOCATE(psi_rel1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  ALLOCATE(th11(1:N2,1:N3,1:2))
  ALLOCATE(th12(0:N1,1:N3,1:2))
  ALLOCATE(th13(0:N1,1:N2,1:2))
  
  ALLOCATE(th21(0:N2,1:N3,1:2))
  ALLOCATE(th22(1:N1,1:N3,1:2))
  ALLOCATE(th23(1:N1,0:N2,1:2))
  
  ALLOCATE(th31(1:N2,0:N3,1:2))
  ALLOCATE(th32(1:N1,0:N3,1:2))
  ALLOCATE(th33(1:N1,1:N2,1:2))
  
  !--- Multigrid ---------------------------------------------------------------------------------------------
  ALLOCATE(vec1C(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  
  !--- BiCGstab / Richardson ---------------------------------------------------------------------------------
  ALLOCATE(pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- product_div_grad --------------------------------------------------------------------------------------
  ALLOCATE(dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Hilfsfeld (kompakte Differenzen) ----------------------------------------------------------------------
  ALLOCATE(com(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Hilfsfelder (Druckiteration) --------------------------------------------------------------------------
  ALLOCATE(work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  ALLOCATE(work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U)))
  
  !--- Lagrange-Partikel -------------------------------------------------------------------------------------
  ALLOCATE(particles(1:n_args,1:n_part_max))
  
  !--- Linienrelaxation --------------------------------------------------------------------------------------
  ALLOCATE(vec1(1:N1))
  ALLOCATE(vec2(1:N2))
  ALLOCATE(vec3(1:N3))
  
  ALLOCATE(dia1(1:N1))
  ALLOCATE(dia2(1:N2))
  ALLOCATE(dia3(1:N3))
  
  ALLOCATE(SOR1(1:N1))
  ALLOCATE(SOR2(1:N2))
  ALLOCATE(SOR3(1:N3))
  
  ALLOCATE(band1(1:2,1:N1))
  ALLOCATE(band2(1:2,1:N2))
  ALLOCATE(band3(1:2,1:N3))
  
  !===========================================================================================================
  !=== Indizierung (Intervallgrenzen, Verschiebungen) ========================================================
  !===========================================================================================================
  !--- Konzentrationen (exklusive Rand) ----------------------------------------------------------------------
  ALLOCATE(S1c(1:n_conc))
  ALLOCATE(S2c(1:n_conc))
  ALLOCATE(S3c(1:n_conc))
  
  ALLOCATE(N1c(1:n_conc))
  ALLOCATE(N2c(1:n_conc))
  ALLOCATE(N3c(1:n_conc))
  
  
  !===========================================================================================================
  !=== Randbedingungen =======================================================================================
  !===========================================================================================================
  !--- lokal (Block) -----------------------------------------------------------------------------------------
  ALLOCATE(BCc_1L(1:n_conc))
  ALLOCATE(BCc_2L(1:n_conc))
  ALLOCATE(BCc_3L(1:n_conc))
  
  ALLOCATE(BCc_1U(1:n_conc))
  ALLOCATE(BCc_2U(1:n_conc))
  ALLOCATE(BCc_3U(1:n_conc))
  
  !--- field properties --------------------------------------------------------------------------------------
  ALLOCATE(recvR(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(recvI(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(dispR(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(dispI(    1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  ALLOCATE(sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max))
  
  
  !===========================================================================================================
  !=== physikalische Parameter ===============================================================================
  !===========================================================================================================
  !ALLOCATE(Sc (1:n_conc)) ! TEST!!! nach usr_config verschoben
  !ALLOCATE(Ric(1:n_conc))
  !ALLOCATE(Rip(1:n_conc))
  !ALLOCATE(usc(1:n_conc))
  !ALLOCATE(usp(1:n_conc))
  !ALLOCATE(Stp(1:n_conc))
  
  ALLOCATE(velReSc1(1:N2,1:N3,1:2,1:n_grids_max))
  ALLOCATE(velReSc2(1:N1,1:N3,1:2,1:n_grids_max))
  ALLOCATE(velReSc3(1:N1,1:N2,1:2,1:n_grids_max))
  
  ALLOCATE(usReSc(1:3,1:n_conc))
  ALLOCATE(us_vec(1:3,1:n_conc))
  

  !===========================================================================================================
  !=== MPI ===================================================================================================
  !===========================================================================================================
  !--- Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes) -----------------------
  ALLOCATE(bar1_size(1:NB1))
  ALLOCATE(bar2_size(1:NB2))
  ALLOCATE(bar3_size(1:NB3))
  
  ALLOCATE(bar1_offset(1:NB1))
  ALLOCATE(bar2_offset(1:NB2))
  ALLOCATE(bar3_offset(1:NB3))
  
  ! --- stored velocity for periodicity criterion
  ALLOCATE(velp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)) !< added
! ALLOCATE(vel (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3))


!#endif
