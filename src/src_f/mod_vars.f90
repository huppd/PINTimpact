!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> \brief contains many variables, like wokring stuff and stencils
module mod_vars
  
  use iso_c_binding

  use mod_dims
  
  implicit none
  
  
#ifdef ALLOC
  
  !===========================================================================================================
  !> \brief Raemliche Dimensionen
  !!
  !! 2D wird ueber M3==2 eingeschaltet
  integer                ::  dimens
  
  
  !===========================================================================================================
  !> \{ \name Domain- und Blockspezifikationen
  !===========================================================================================================
  !--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  ! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !--------------------------------------------------------------
  ! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  ! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  ! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  ! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  ! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  ! ...
  !> \{
  !! \brief Domain grid points
  integer                ::  N1
  integer                ::  N2
  integer                ::  N3
  !> \}
  
  !> \{
  !! \brief --- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  integer, parameter     ::  n_grids_max = 15
  integer                ::  n_grids, n_grids_limit
  !>\}
  
  !> \{
  !! \brief --- Dimensionen -------------------------------------------------------------------------------------------
  integer, parameter     ::  dim_ncb1c = SIZE(ncb1c)
  integer, parameter     ::  dim_ncb1g = SIZE(ncb1g)
  integer, parameter     ::  dim_ncb1d = SIZE(ncb1d)
  
  integer, parameter     ::  dim_ncb2c = SIZE(ncb2c)
  integer, parameter     ::  dim_ncb2g = SIZE(ncb2g)
  integer, parameter     ::  dim_ncb2d = SIZE(ncb2d)
  
  integer, parameter     ::  dim_ncb3c = SIZE(ncb3c)
  integer, parameter     ::  dim_ncb3g = SIZE(ncb3g)
  integer, parameter     ::  dim_ncb3d = SIZE(ncb3d)
  !> \}
  
  !> \{
  !! \brief Anzahl Stencil-Koeffizienten (Feld)
  !!
  !! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  integer, parameter     ::  nc1c = ncb1c(dim_ncb1c)
  integer, parameter     ::  nc1s = ncb1g(dim_ncb1g)
  
  integer, parameter     ::  nc2c = ncb2c(dim_ncb2c)
  integer, parameter     ::  nc2s = ncb2g(dim_ncb2g)
  
  integer, parameter     ::  nc3c = ncb3c(dim_ncb3c)
  integer, parameter     ::  nc3s = ncb3g(dim_ncb3g)
  
  !> \} \}
  !===========================================================================================================
  !> \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  !> \{ \brief zentral
  integer, parameter     ::  b1U = nc1s/2
  integer, parameter     ::  b2U = nc2s/2
  integer, parameter     ::  b3U = nc3s/2
  
  integer, parameter     ::  b1L = -b1U
  integer, parameter     ::  b2L = -b2U
  integer, parameter     ::  b3L = -b3U
  !> \}
  
  !> \{ \brief upwind (nicht-linear)
  !!
  !! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !!  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  integer, parameter     ::  n1L = b1L
  integer, parameter     ::  n2L = b2L
  integer, parameter     ::  n3L = b3L
  
  integer, parameter     ::  n1U = b1U
  integer, parameter     ::  n2U = b2U
  integer, parameter     ::  n3U = b3U
  !> \}
  
  !> \{ \brief Divergenz
  integer, parameter     ::  d1L = b1L
  integer, parameter     ::  d2L = b2L
  integer, parameter     ::  d3L = b3L
  
  integer, parameter     ::  d1U = b1U-1
  integer, parameter     ::  d2U = b2U-1
  integer, parameter     ::  d3U = b3U-1
  !> \}
  
  !> \{ \brief Gradient
  integer, parameter     ::  g1L = b1L+1
  integer, parameter     ::  g2L = b2L+1
  integer, parameter     ::  g3L = b3L+1
  
  integer, parameter     ::  g1U = b1U
  integer, parameter     ::  g2U = b2U
  integer, parameter     ::  g3U = b3U
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief 1. Ableitung (zentral)
  real   , allocatable   ::  cp1  (:,:)
  real   , allocatable   ::  cp2  (:,:)
  real   , allocatable   ::  cp3  (:,:)
  
  real   , allocatable   ::  cu1  (:,:)
  real   , allocatable   ::  cv2  (:,:)
  real   , allocatable   ::  cw3  (:,:)
  
  ! Eigene Stencils notwendig, da
  !   - Symmetrie des Geschwindigkeitsfeldes nicht unbedingt Symmetrie bei Konzentrationen bedeutet,
  !   - Randbedinungen darauf gespeichert werden koennen.
  real   , allocatable   ::  cc1  (:,:,:)
  real   , allocatable   ::  cc2  (:,:,:)
  real   , allocatable   ::  cc3  (:,:,:)
  
  !> \}
  !> \{ \brief 1. Ableitung (upwind)
  real   , allocatable   ::  cNp1D(:,:)
  real   , allocatable   ::  cNp2D(:,:)
  real   , allocatable   ::  cNp3D(:,:)
  
  real   , allocatable   ::  cNp1U(:,:)
  real   , allocatable   ::  cNp2U(:,:)
  real   , allocatable   ::  cNp3U(:,:)
  
  real   , allocatable   ::  cNu1D(:,:)
  real   , allocatable   ::  cNv2D(:,:)
  real   , allocatable   ::  cNw3D(:,:)
  
  real   , allocatable   ::  cNu1U(:,:)
  real   , allocatable   ::  cNv2U(:,:)
  real   , allocatable   ::  cNw3U(:,:)
  
  real   , allocatable   ::  cNc1D(:,:,:)
  real   , allocatable   ::  cNc2D(:,:,:)
  real   , allocatable   ::  cNc3D(:,:,:)
  
  real   , allocatable   ::  cNc1U(:,:,:)
  real   , allocatable   ::  cNc2U(:,:,:)
  real   , allocatable   ::  cNc3U(:,:,:)
  !> \}
  
  !> \{ \brief Divergenz
  real   , allocatable   ::  cDu1 (:,:)
  real   , allocatable   ::  cDv2 (:,:)
  real   , allocatable   ::  cDw3 (:,:)
  !> \}
  
  !> \{ \brief Divergenz (transponiert)
  real   , allocatable   ::  cDu1T(:,:)
  real   , allocatable   ::  cDv2T(:,:)
  real   , allocatable   ::  cDw3T(:,:)
  !> \}
  
  !> \{ \brief Gradient
  real   , allocatable   ::  cGp1 (:,:)
  real   , allocatable   ::  cGp2 (:,:)
  real   , allocatable   ::  cGp3 (:,:)
  !> \}
  
  !> \{ \brief Gradient (transponiert)
  real   , allocatable   ::  cGp1T(:,:)
  real   , allocatable   ::  cGp2T(:,:)
  real   , allocatable   ::  cGp3T(:,:)
  !> \}
  
  !> \{ \brief 2. Ableitung (zentral)
  real   , allocatable   ::  cp11 (:,:)
  real   , allocatable   ::  cp22 (:,:)
  real   , allocatable   ::  cp33 (:,:)
  
  real   , allocatable   ::  cu11 (:,:)
  real   , allocatable   ::  cv22 (:,:)
  real   , allocatable   ::  cw33 (:,:)
  
  real   , allocatable   ::  cc11 (:,:,:)
  real   , allocatable   ::  cc22 (:,:,:)
  real   , allocatable   ::  cc33 (:,:,:)
  !> \}
  
  !> \{ \brief Interpolation
  real   , allocatable   ::  cIpu(:,:)
  real   , allocatable   ::  cIpv(:,:)
  real   , allocatable   ::  cIpw(:,:)
  
  real   , allocatable   ::  cIup(:,:)
  real   , allocatable   ::  cIvp(:,:)
  real   , allocatable   ::  cIwp(:,:)
  
  real   , allocatable   ::  cIcu(:,:,:)
  real   , allocatable   ::  cIcv(:,:,:)
  real   , allocatable   ::  cIcw(:,:,:)
  !> \}
  
  !> \{ \brief Filter
  real   , allocatable   ::  cFp1(:,:)
  real   , allocatable   ::  cFp2(:,:)
  real   , allocatable   ::  cFp3(:,:)
  
  real   , allocatable   ::  cFu1(:,:)
  real   , allocatable   ::  cFv2(:,:)
  real   , allocatable   ::  cFw3(:,:)
  
  real   , allocatable   ::  cFc1(:,:,:)
  real   , allocatable   ::  cFc2(:,:,:)
  real   , allocatable   ::  cFc3(:,:,:)
  !> \}
  
  !> \{ \brief Integrator (nur für Druckgitter)
  real   , allocatable   ::  cInt1(:,:)
  real   , allocatable   ::  cInt2(:,:)
  real   , allocatable   ::  cInt3(:,:)
  !> \}
  
  
  !> \{ \brief Compact
  integer                ::  dimS1
  integer                ::  dimS2
  integer                ::  dimS3
  
  real   , allocatable   ::  buffer1A(:,:,:), buffer1B(:,:,:)
  real   , allocatable   ::  buffer2A(:,:,:), buffer2B(:,:,:)
  real   , allocatable   ::  buffer3A(:,:,:), buffer3B(:,:,:)
  !----------------------------------------------------
  integer, allocatable   ::  disp1(:), recv1(:)
  integer, allocatable   ::  disp2(:), recv2(:)
  integer, allocatable   ::  disp3(:), recv3(:)
  !----------------------------------------------------
  !> \}
  
  !> \{ \brief Schur complement (square) matrices ---
  real   , allocatable   ::  Sp1(:,:), Su1(:,:), Sc1(:,:,:), Sp11(:,:), Su11(:,:), Sc11(:,:,:)
  real   , allocatable   ::  Sp2(:,:), Sv2(:,:), Sc2(:,:,:), Sp22(:,:), Sv22(:,:), Sc22(:,:,:)
  real   , allocatable   ::  Sp3(:,:), Sw3(:,:), Sc3(:,:,:), Sp33(:,:), Sw33(:,:), Sc33(:,:,:)
  
  real   , allocatable   ::  SDu1(:,:), SGp1(:,:), SIpu(:,:), SIup(:,:), SIcu(:,:,:)
  real   , allocatable   ::  SDv2(:,:), SGp2(:,:), SIpv(:,:), SIvp(:,:), SIcv(:,:,:)
  real   , allocatable   ::  SDw3(:,:), SGp3(:,:), SIpw(:,:), SIwp(:,:), SIcw(:,:,:)
  
  real   , allocatable   ::  SDu1T(:,:), SGp1T(:,:)
  real   , allocatable   ::  SDv2T(:,:), SGp2T(:,:)
  real   , allocatable   ::  SDw3T(:,:), SGp3T(:,:)
  !> \}
  
  !> \{ \brief bottom (horizontally oriented) matrices
  real   , allocatable   ::  Wp1(:,:), Wu1(:,:), Wc1(:,:,:), Wp11(:,:), Wu11(:,:), Wc11(:,:,:)
  real   , allocatable   ::  Wp2(:,:), Wv2(:,:), Wc2(:,:,:), Wp22(:,:), Wv22(:,:), Wc22(:,:,:)
  real   , allocatable   ::  Wp3(:,:), Ww3(:,:), Wc3(:,:,:), Wp33(:,:), Ww33(:,:), Wc33(:,:,:)
  
  real   , allocatable   ::  WDu1(:,:), WGp1(:,:), WIpu(:,:), WIup(:,:), WIcu(:,:,:)
  real   , allocatable   ::  WDv2(:,:), WGp2(:,:), WIpv(:,:), WIvp(:,:), WIcv(:,:,:)
  real   , allocatable   ::  WDw3(:,:), WGp3(:,:), WIpw(:,:), WIwp(:,:), WIcw(:,:,:)
  
  real   , allocatable   ::  WDu1T(:,:), WGp1T(:,:)
  real   , allocatable   ::  WDv2T(:,:), WGp2T(:,:)
  real   , allocatable   ::  WDw3T(:,:), WGp3T(:,:)
  !> \}
  
  !> \{ \brief implicit finite difference coefficients ("L"=left-hand side)
  real   , allocatable   ::  cp1CL(:,:), cu1CL(:,:), cc1CL(:,:,:), cp11CL(:,:), cu11CL(:,:), cc11CL(:,:,:)
  real   , allocatable   ::  cp2CL(:,:), cv2CL(:,:), cc2CL(:,:,:), cp22CL(:,:), cv22CL(:,:), cc22CL(:,:,:)
  real   , allocatable   ::  cp3CL(:,:), cw3CL(:,:), cc3CL(:,:,:), cp33CL(:,:), cw33CL(:,:), cc33CL(:,:,:)
  
  real   , allocatable   ::  cDu1CL(:,:), cGp1CL(:,:), cIpuCL(:,:), cIupCL(:,:), cIcuCL(:,:,:)
  real   , allocatable   ::  cDv2CL(:,:), cGp2CL(:,:), cIpvCL(:,:), cIvpCL(:,:), cIcvCL(:,:,:)
  real   , allocatable   ::  cDw3CL(:,:), cGp3CL(:,:), cIpwCL(:,:), cIwpCL(:,:), cIcwCL(:,:,:)
  
  real   , allocatable   ::  cDu1CLT(:,:), cGp1CLT(:,:)
  real   , allocatable   ::  cDv2CLT(:,:), cGp2CLT(:,:)
  real   , allocatable   ::  cDw3CLT(:,:), cGp3CLT(:,:)
  
  real   , allocatable   ::  cp1CL_LU(:,:), cu1CL_LU(:,:), cc1CL_LU(:,:,:), cp11CL_LU(:,:), cu11CL_LU(:,:), cc11CL_LU(:,:,:)
  real   , allocatable   ::  cp2CL_LU(:,:), cv2CL_LU(:,:), cc2CL_LU(:,:,:), cp22CL_LU(:,:), cv22CL_LU(:,:), cc22CL_LU(:,:,:)
  real   , allocatable   ::  cp3CL_LU(:,:), cw3CL_LU(:,:), cc3CL_LU(:,:,:), cp33CL_LU(:,:), cw33CL_LU(:,:), cc33CL_LU(:,:,:)
  
  real   , allocatable   ::  cDu1CL_LU(:,:), cGp1CL_LU(:,:), cIpuCL_LU(:,:), cIupCL_LU(:,:), cIcuCL_LU(:,:,:)
  real   , allocatable   ::  cDv2CL_LU(:,:), cGp2CL_LU(:,:), cIpvCL_LU(:,:), cIvpCL_LU(:,:), cIcvCL_LU(:,:,:)
  real   , allocatable   ::  cDw3CL_LU(:,:), cGp3CL_LU(:,:), cIpwCL_LU(:,:), cIwpCL_LU(:,:), cIcwCL_LU(:,:,:)
  
  real   , allocatable   ::  cDu1CLT_LU(:,:), cGp1CLT_LU(:,:)
  real   , allocatable   ::  cDv2CLT_LU(:,:), cGp2CLT_LU(:,:)
  real   , allocatable   ::  cDw3CLT_LU(:,:), cGp3CLT_LU(:,:)
  !> \}
  
  !> \{ \brief explicit finite difference coefficients ("R"=right-hand side)
  real   , allocatable   ::  cp1CR(:,:), cu1CR(:,:), cc1CR(:,:,:), cp11CR(:,:), cu11CR(:,:), cc11CR(:,:,:)
  real   , allocatable   ::  cp2CR(:,:), cv2CR(:,:), cc2CR(:,:,:), cp22CR(:,:), cv22CR(:,:), cc22CR(:,:,:)
  real   , allocatable   ::  cp3CR(:,:), cw3CR(:,:), cc3CR(:,:,:), cp33CR(:,:), cw33CR(:,:), cc33CR(:,:,:)
  
  real   , allocatable   ::  cDu1CR(:,:), cGp1CR(:,:), cIpuCR(:,:), cIupCR(:,:), cIcuCR(:,:,:)
  real   , allocatable   ::  cDv2CR(:,:), cGp2CR(:,:), cIpvCR(:,:), cIvpCR(:,:), cIcvCR(:,:,:)
  real   , allocatable   ::  cDw3CR(:,:), cGp3CR(:,:), cIpwCR(:,:), cIwpCR(:,:), cIcwCR(:,:,:)
  
  real   , allocatable   ::  cDu1CRT(:,:), cGp1CRT(:,:)
  real   , allocatable   ::  cDv2CRT(:,:), cGp2CRT(:,:)
  real   , allocatable   ::  cDw3CRT(:,:), cGp3CRT(:,:)
  !> \}
  !****************************************************
  
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  real   , allocatable   ::  cp11R(:,:,:)
  real   , allocatable   ::  cp22R(:,:,:)
  real   , allocatable   ::  cp33R(:,:,:)
  
  real   , allocatable   ::  cu11R(:,:,:)
  real   , allocatable   ::  cv22R(:,:,:)
  real   , allocatable   ::  cw33R(:,:,:)
  
  real   , allocatable   ::  cc11R(:,:,:,:)
  real   , allocatable   ::  cc22R(:,:,:,:)
  real   , allocatable   ::  cc33R(:,:,:,:)
  
  real   , allocatable   ::  cdg1 (:,:,:)
  real   , allocatable   ::  cdg2 (:,:,:)
  real   , allocatable   ::  cdg3 (:,:,:)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  real   , allocatable   ::  cI1(:,:,:)
  real   , allocatable   ::  cI2(:,:,:)
  real   , allocatable   ::  cI3(:,:,:)
  
  real   , allocatable   ::  cIH1(:,:)
  real   , allocatable   ::  cIH2(:,:)
  real   , allocatable   ::  cIH3(:,:)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  real   , allocatable   ::  cR1 (:,:,:)
  real   , allocatable   ::  cR2 (:,:,:)
  real   , allocatable   ::  cR3 (:,:,:)
  
  real   , allocatable   ::  cRest1(:,:,:) ! TEST!!!
  real   , allocatable   ::  cRest2(:,:,:)
  real   , allocatable   ::  cRest3(:,:,:)
  
  real   , allocatable   ::  cRH1(:,:)
  real   , allocatable   ::  cRH2(:,:)
  real   , allocatable   ::  cRH3(:,:)
  
  
  !===========================================================================================================
  !> \{ \name Gitterspezifikationen
  !===========================================================================================================
  !> \{
  !! \brief physiklische Koordinaten (global)
  real   , allocatable   ::  y1p(:), y1u(:)
  real   , allocatable   ::  y2p(:), y2v(:)
  real   , allocatable   ::  y3p(:), y3w(:)
  !> \}
  
  !> \{
  !! \brief physiklische Koordinaten (Block)
  real   , allocatable   ::  x1p(:), x1u(:)
  real   , allocatable   ::  x2p(:), x2v(:)
  real   , allocatable   ::  x3p(:), x3w(:)
  !> \}
  
  !> \{ \brief physiklische Koordinaten (Block, Multigrid)
  real   , allocatable   ::  x1pR(:,:), x1uR(:,:)
  real   , allocatable   ::  x2pR(:,:), x2vR(:,:)
  real   , allocatable   ::  x3pR(:,:), x3wR(:,:)
  !> \}
  
  !> \{ \brief Gitterweiten (global)
  real   , allocatable   ::  dy1p(:), dy1u(:)
  real   , allocatable   ::  dy2p(:), dy2v(:)
  real   , allocatable   ::  dy3p(:), dy3w(:)
  !> \}
  
  !> \{ \brief Gitterweiten (Block)
  real   , allocatable   ::  dx1p(:), dx1u(:)
  real   , allocatable   ::  dx2p(:), dx2v(:)
  real   , allocatable   ::  dx3p(:), dx3w(:)
  !> \}
  
  real   , allocatable   ::  dx1DM(:), dx1pM(:), ddx1pM(:)
  real   , allocatable   ::  dx2DM(:), dx2pM(:), ddx2pM(:)
  real   , allocatable   ::  dx3DM(:), dx3pM(:), ddx3pM(:)
  
  real   , allocatable   ::  dx1GM(:), dx1uM(:), ddx1uM(:)
  real   , allocatable   ::  dx2GM(:), dx2vM(:), ddx2vM(:)
  real   , allocatable   ::  dx3GM(:), dx3wM(:), ddx3wM(:)
  
  !> \{ \brief Smagorinsky-Modell
  real   , allocatable   ::  dx1pS(:), dx1uS(:) ! TEST!!!
  real   , allocatable   ::  dx2pS(:), dx2vS(:)
  real   , allocatable   ::  dx3pS(:), dx3wS(:)
  !> \}
  
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Arbeitsfelder
  !===========================================================================================================
  !> \brief Geschwindigkeiten
  real   , allocatable   ::  vel(:,:,:,:)
  
  !> \brief nicht-linearer Term
  real   , allocatable   ::  nl (:,:,:,:)
  
  !> \brief Recht-Hand-Seite
  real   , allocatable   ::  rhs(:,:,:,:)
  
  !> \brief Druck
  real   , allocatable   ::  pre(:,:,:)
  
  !> \brief Konzentrationsfelder
  real   , allocatable   ::  conc(:,:,:,:)
  
  !> \brief nicht-linearer Term (Konzentration)
  real   , allocatable   ::  nlco(:,:,:,:)
  
  !> \{ \brief Ausfluss-RB (Geschwindigkeitsfeld)
  !!
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  real   , allocatable   ::  bc11(:,:,:), nlbc11(:,:,:)
  real   , allocatable   ::  bc12(:,:,:), nlbc12(:,:,:)
  real   , allocatable   ::  bc13(:,:,:), nlbc13(:,:,:)
  
  real   , allocatable   ::  bc21(:,:,:), nlbc21(:,:,:)
  real   , allocatable   ::  bc22(:,:,:), nlbc22(:,:,:)
  real   , allocatable   ::  bc23(:,:,:), nlbc23(:,:,:)
  
  real   , allocatable   ::  bc31(:,:,:), nlbc31(:,:,:)
  real   , allocatable   ::  bc32(:,:,:), nlbc32(:,:,:)
  real   , allocatable   ::  bc33(:,:,:), nlbc33(:,:,:)
  
  real   , allocatable   ::  drift1(:,:,:)
  real   , allocatable   ::  drift2(:,:,:)
  real   , allocatable   ::  drift3(:,:,:)
  !> \}
  
  
  !> \{ \brief Sedimentations-RB (Konzentration)
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehörigen Randbedingungen gespeichert werden.
  real   , allocatable   ::  sed1L (:,:,:)
  real   , allocatable   ::  sed1U (:,:,:)
  
  real   , allocatable   ::  sed2L (:,:,:)
  real   , allocatable   ::  sed2U (:,:,:)
  
  real   , allocatable   ::  sed3L (:,:,:)
  real   , allocatable   ::  sed3U (:,:,:)
  
  real   , allocatable   ::  conc1L(:,:,:)
  real   , allocatable   ::  conc1U(:,:,:)
  
  real   , allocatable   ::  conc2L(:,:,:)
  real   , allocatable   ::  conc2U(:,:,:)
  
  real   , allocatable   ::  conc3L(:,:,:)
  real   , allocatable   ::  conc3U(:,:,:)
  !> \}
  
  !> \{ \brief Sedimentation (Konzentration, Ausschrieb)
  real   , allocatable   ::  dep1L_conc(:,:,:)
  real   , allocatable   ::  dep1U_conc(:,:,:)
  
  real   , allocatable   ::  dep2L_conc(:,:,:)
  real   , allocatable   ::  dep2U_conc(:,:,:)
  
  real   , allocatable   ::  dep3L_conc(:,:,:)
  real   , allocatable   ::  dep3U_conc(:,:,:)
  !> \}
  
  !> \< \brief Sedimentation (Partikel, Ausschrieb)
  ! TEST!!! 1:n_spec???
  real   , allocatable   ::  dep1L_part(:,:)
  real   , allocatable   ::  dep1U_part(:,:)
  
  real   , allocatable   ::  dep2L_part(:,:)
  real   , allocatable   ::  dep2U_part(:,:)
  
  real   , allocatable   ::  dep3L_part(:,:)
  real   , allocatable   ::  dep3U_part(:,:)
  !> \}
  
  !> \brief Residuum
  real   , allocatable   ::  res (:,:,:)
  
  !> \brief Druckgradient (eine Komponente)
  real   , allocatable   ::  gpre(:,:,:)
  
  !> \brief Gewichte f�r Divergenzfreiheit
  real   , allocatable   ::  weight(:,:,:)
  
  !> \brief fluid density on velocity grid
#ifdef NONBOUSSINESQ
  real   , allocatable   ::  dens(:,:,:,:)
#endif
  
  !> \{ \brief Null-Raum-Vektor
  real   , allocatable   ::  psi    (:,:,:)
  real   , allocatable   ::  psi_vel(:,:,:,:)
  
  real   , allocatable   ::  psi_rel1 (:,:,:)
  real   , allocatable   ::  psi_rel2 (:,:,:)
  real   , allocatable   ::  psi_rel3 (:,:,:)
  real   , allocatable   ::  psi_rel4 (:,:,:)
  real   , allocatable   ::  psi_rel5 (:,:,:)
  real   , allocatable   ::  psi_rel6 (:,:,:)
  real   , allocatable   ::  psi_rel7 (:,:,:)
  real   , allocatable   ::  psi_rel8 (:,:,:)
  real   , allocatable   ::  psi_rel9 (:,:,:)
  real   , allocatable   ::  psi_rel10(:,:,:)
  real   , allocatable   ::  psi_rel11(:,:,:)
  real   , allocatable   ::  psi_rel12(:,:,:)
  real   , allocatable   ::  psi_rel13(:,:,:)
  real   , allocatable   ::  psi_rel14(:,:,:)
  real   , allocatable   ::  psi_rel15(:,:,:)
  
  real   , allocatable   ::  th11(:,:,:)
  real   , allocatable   ::  th12(:,:,:)
  real   , allocatable   ::  th13(:,:,:)
  
  real   , allocatable   ::  th21(:,:,:)
  real   , allocatable   ::  th22(:,:,:)
  real   , allocatable   ::  th23(:,:,:)
  
  real   , allocatable   ::  th31(:,:,:)
  real   , allocatable   ::  th32(:,:,:)
  real   , allocatable   ::  th33(:,:,:)
  !> \}
  
  !> \{ \brief Multigrid
  real   , allocatable   ::  vec1C (:,:,:)
  
  real   , allocatable   ::  vec2A (:,:,:)
  real   , allocatable   ::  vec2B (:,:,:)
  real   , allocatable   ::  vec2C (:,:,:)
  
  real   , allocatable   ::  vec3A (:,:,:)
  real   , allocatable   ::  vec3B (:,:,:)
  real   , allocatable   ::  vec3C (:,:,:)
  
  real   , allocatable   ::  vec4A (:,:,:)
  real   , allocatable   ::  vec4B (:,:,:)
  real   , allocatable   ::  vec4C (:,:,:)
  
  real   , allocatable   ::  vec5A (:,:,:)
  real   , allocatable   ::  vec5B (:,:,:)
  real   , allocatable   ::  vec5C (:,:,:)
  
  real   , allocatable   ::  vec6A (:,:,:)
  real   , allocatable   ::  vec6B (:,:,:)
  real   , allocatable   ::  vec6C (:,:,:)
  
  real   , allocatable   ::  vec7A (:,:,:)
  real   , allocatable   ::  vec7B (:,:,:)
  real   , allocatable   ::  vec7C (:,:,:)
  
  real   , allocatable   ::  vec8A (:,:,:)
  real   , allocatable   ::  vec8B (:,:,:)
  real   , allocatable   ::  vec8C (:,:,:)
  
  real   , allocatable   ::  vec9A (:,:,:)
  real   , allocatable   ::  vec9B (:,:,:)
  real   , allocatable   ::  vec9C (:,:,:)
  
  real   , allocatable   ::  vec10A(:,:,:)
  real   , allocatable   ::  vec10B(:,:,:)
  real   , allocatable   ::  vec10C(:,:,:)
  
  real   , allocatable   ::  vec11A(:,:,:)
  real   , allocatable   ::  vec11B(:,:,:)
  real   , allocatable   ::  vec11C(:,:,:)
  
  real   , allocatable   ::  vec12A(:,:,:)
  real   , allocatable   ::  vec12B(:,:,:)
  real   , allocatable   ::  vec12C(:,:,:)
  
  real   , allocatable   ::  vec13A(:,:,:)
  real   , allocatable   ::  vec13B(:,:,:)
  real   , allocatable   ::  vec13C(:,:,:)
  
  real   , allocatable   ::  vec14A(:,:,:)
  real   , allocatable   ::  vec14B(:,:,:)
  real   , allocatable   ::  vec14C(:,:,:)
  
  real   , allocatable   ::  vec15A(:,:,:)
  real   , allocatable   ::  vec15B(:,:,:)
  !> \}
  
  
  !> \{ \brief BiCGstab / Richardson
  real   , allocatable   ::  pp(:,:,:)
  real   , allocatable   ::  Ap(:,:,:)
  real   , allocatable   ::  rr(:,:,:)
  real   , allocatable   ::  rh(:,:,:)
  real   , allocatable   ::  Ar(:,:,:)
  real   , allocatable   ::  z1(:,:,:)
  real   , allocatable   ::  z2(:,:,:)
  !> \}
  
  !> \brief product_div_grad
  real   , allocatable   ::  dig(:,:,:)
  
  !> \brief Hilfsfeld (kompakte Differenzen)
  real   , allocatable   ::  com(:,:,:)
  
  !> \{ \brief Hilfsfelder (Druckiteration)
  !!
  !! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  !! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  real   , allocatable   ::  work1(:,:,:)
  real   , allocatable   ::  work2(:,:,:)
  real   , allocatable   ::  work3(:,:,:)
  !> \}
  
  !> \{ \brief Lagrange-Partikel
  integer                ::  n_part, n_part_tot, n_groups
  real   , allocatable   ::  particles(:,:)
  !> \}
  
  !> \brief Smagorinsky LES-Modell
  real   , allocatable   ::  smag_diffus(:,:,:,:)
  
  !> \{ \brief Linienrelaxation
  real   , allocatable   ::  vec1(:), dia1(:), SOR1(:), band1(:,:) ! TEST!!! siehe unten ...
  real   , allocatable   ::  vec2(:), dia2(:), SOR2(:), band2(:,:)
  real   , allocatable   ::  vec3(:), dia3(:), SOR3(:), band3(:,:)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Indizierung (Intervallgrenzen, Verschiebungen)
  !!
  !!===========================================================================================================
  !!--- Block-Index -------------------------------------------------------------------------------------------
  !!INTEGER               ::  iB, jB, kB ! TEST!!! iB(1:3,1:n_grids_max) hierher verschieben ...
  
  !> \{ \brief Indexverschiebung (Block --> global)
  integer                ::  iShift, jShift, kShift
  !> \}
  
  !--- Domaingrösse (Periodizität-bereinigt) -----------------------------------------------------------------
  integer                ::  dim1, dim2, dim3
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  integer                ::  S1p, S2p, S3p
  integer                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  integer                ::  S11B, S21B, S31B
  integer                ::  S12B, S22B, S32B
  integer                ::  S13B, S23B, S33B
  
  integer                ::  N11B, N21B, N31B
  integer                ::  N12B, N22B, N32B
  integer                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  integer                ::  S11, S21, S31
  integer                ::  S12, S22, S32
  integer                ::  S13, S23, S33
  
  integer                ::  N11, N21, N31
  integer                ::  N12, N22, N32
  integer                ::  N13, N23, N33
  
  !--- Konzentrationen (exklusive Rand) ----------------------------------------------------------------------
  integer, allocatable   ::  S1c(:), S2c(:), S3c(:)
  integer, allocatable   ::  N1c(:), N2c(:), N3c(:)
  
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  integer                ::  S1R, S2R, S3R
  integer                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  integer                ::  S11R, S22R, S33R
  integer                ::  d11R, d22R, d33R
  
  !> \{ \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)
  integer, parameter     ::  ls1 = -1
  integer, parameter     ::  ls2 = -1
  integer, parameter     ::  ls3 = -1
  !> \}
  
  !> \{ \brief Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  !!
  !! ex = -1: unten <--  oben
  !! ex =  0: unten <--> oben
  !! ex =  1: unten  --> oben
  integer                ::  ex1, ex2, ex3
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Boundary conditions
  !!
  !!    symmetrische, zentrale Stencils:
  !!     - Symmetrie-RB:   BC = -2
  !!     - periodische RB: BC = -1
  !!
  !!    schiefe, nicht-zentrale Stencils:
  !!     - Nachbar-Block:  BC =  0
  !!     - Dirichlet-RB:   BC =  1
  !!     - Neumann-RB:     BC =  2
  !!     - Robin-RB:       BC =  3

  !> \{ \brief global boundary conditions
  logical                ::  outlet(1:3,1:2,1:3)
  logical, allocatable   ::  isopycnal(:,:,:)
  
  integer                ::  BC_1L_global, BC_1U_global
  integer                ::  BC_2L_global, BC_2U_global
  integer                ::  BC_3L_global, BC_3U_global
  !> \}
  
  !> \{ \brief lokal (Block)
  integer                ::  BC_1L, BC_1U
  integer                ::  BC_2L, BC_2U
  integer                ::  BC_3L, BC_3U
  
  integer, allocatable   ::  BCc_1L(:), BCc_1U(:)
  integer, allocatable   ::  BCc_2L(:), BCc_2U(:)
  integer, allocatable   ::  BCc_3L(:), BCc_3U(:)
  !> \}
  
  !> \{ \brief field properties
  integer                ::  n_gather(1:3,1:n_grids_max)
  integer                ::  NN (1:3,1:n_grids_max)
  integer                ::  NB (1:3,1:n_grids_max)
  integer                ::  iB (1:3,1:n_grids_max)
  integer                ::  SNF(1:2,1:3,1:n_grids_max)
  integer                ::  SNB(1:2,1:3,1:n_grids_max)
  !> first index lower upper/ second index direction, third index grid level
  integer                ::  BC (1:2,1:3,1:n_grids_max)
  integer                ::  ngb(1:2,1:3,1:n_grids_max)
  integer                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  integer                ::  rankc2(1:n_grids_max)
  logical                ::  participate_yes(1:n_grids_max)
  integer, allocatable   ::  recvR(  :,:), recvI(  :,:)
  integer, allocatable   ::  dispR(  :,:), dispI(  :,:)
  integer, allocatable   ::  offsR(:,:,:), offsI(:,:,:)
  integer, allocatable   ::  sizsR(:,:,:), sizsI(:,:,:)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ name physikalische Parameter
  !===========================================================================================================
  real                   ::  L1, L2, L3
  real                   ::  Re, Fro
  real   , allocatable   ::  Sc (:)
  real   , allocatable   ::  Ric(:)
  real   , allocatable   ::  Rip(:) ! TEST!!!
  real   , allocatable   ::  usc(:) ! TEST!!! ebenfalls umbenennen!
  real   , allocatable   ::  usp(:) ! TEST!!!
  real   , allocatable   ::  Stp(:) ! TEST!!!
  
  real   , allocatable   ::  velReSc1(:,:,:,:)
  real   , allocatable   ::  velReSc2(:,:,:,:)
  real   , allocatable   ::  velReSc3(:,:,:,:)
  
  real   , allocatable   ::  usReSc (:,:)
  real   , allocatable   ::  us_vec (:,:)
  real                   ::  gravity(1:3)
  
  real                   ::  dx_1L, dx_1U
  real                   ::  dx_2L, dx_2U
  real                   ::  dx_3L, dx_3U
  
  real                   ::  idx_1L, idx_1U
  real                   ::  idx_2L, idx_2U
  real                   ::  idx_3L, idx_3U
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name numerische Parameter
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  real                   ::  CFL
  real                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  integer                ::  timestep, timestep_old, substep, n_timesteps
  logical                ::  mapping_yes, upwind_yes, upwind_conc_yes
  logical                ::  Euler_yes, nonBoussinesq_yes, Stokes_yes, twostep_yes
  logical                ::  comp_visc_yes, comp_conv_yes, comp_inter_yes, comp_div_yes, comp_grad_yes
  logical, parameter     ::  filter_BC_yes = .true. ! TEST!!!
  integer                ::  timeint_mode, LES_mode, forcing_mode
  integer                ::  bulkflow_dir
  integer                ::  n_lp_vel    , n_hp_vel
  integer, allocatable   ::  n_lp_conc(:), n_hp_conc(:)
  real                   ::  chi_vel
  real   , allocatable   ::  chi_conc (:)
  real                   ::  vel_bulk ! TEST!!!
  
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  real   , parameter     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  real   , parameter     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  integer, parameter     ::  RK_steps = 3
  
  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  real   , parameter     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)
  integer, parameter     ::  n_stab = SIZE(stabilitylimit)
  
  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  real                   ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  integer                ::  Int_dtime, Int_lev_pre
  
  integer                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  logical                ::  write_large, write_med, write_small
  real                   ::  time_out_scal, dtime_out_scal
  real                   ::  time_out_vect, dtime_out_vect
  
  logical                ::  write_out_scal, write_out_vect
  logical                ::  new_dtime, finish_yes
    
  integer                ::  write_count
  integer                ::  restart
  character(len=3)       ::  restart_char
  
  integer                ::  n_conc_old

  
  !> \}
  !===========================================================================================================
  !> \{ \name weitere Steuerungsoptionen
  !===========================================================================================================
  integer                ::  task
  logical                ::  read_nullspace_yes
  logical                ::  concentration_yes
  logical                ::  particles_yes
  logical                ::  nullspace_yes, nullspace_coarse_yes
  logical                ::  nullspace_ortho_yes
  
  logical                ::  write_stout_yes
  logical                ::  log_iteration_yes
  logical                ::  write_restart_yes
  logical                ::  write_lambda2_yes
  logical                ::  write_test_yes
  
  !> \{ \brief globale Laufindizes
  integer                ::  direction
  integer                ::  conc_nu
  !> \}
  
  !> \brief explizite Behandlung von Ecken bei Dirichlet-Randbedingungen
  !!
  !! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort k�nstlich zu Null gesetzt wird.)
  logical, parameter     ::  corner_yes = .true.
  
  !> \{ \brief Systemzeit
  integer                ::  elatime, ctime(1:8)
  integer                ::  day, hour, minu, sec, msec
  !> \}
  
  !> \brief Partikel-Energie
  real                   ::  Ekin_part
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsparameter
  !===========================================================================================================
  !> \{ Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten
  real                   ::  epsU, epsU0
  !> \}
  
  !> \brief Glaetter
  logical                ::  Jacobi_yes
  
  !> \{ \brief max. Anzahl Iterationen
  integer                ::  n_it_outer
  integer                ::  n_it_Poisson
  integer                ::  n_it_Helmh_vel
  integer                ::  n_it_Helmh_conc
  !> \}
  
  !> \{ \brief erwartete Konvergenzrate (äussere Iteration)
  real                   ::  precRatio0 (1:RK_steps)
  real                   ::  precOffset0(1:RK_steps)
  real   , allocatable   ::  precRatio  (:,:)
  real   , allocatable   ::  precOffset (:,:)
  !> \}
  
  !> \{ \brief Null-Initialisierung (äussere Iteration)
  logical                ::  init_pre(1:RK_steps), init_vel(1:RK_steps), init_conc(1:RK_steps)
  !> \}
  
  !> \{ \name Vorkonditionierung (Multigrid)
  integer                ::  precond_outer
  integer                ::  precond_Poisson
  integer                ::  precond_Helmh_vel
  integer                ::  precond_Helmh_conc
  !> \}
  
  !> \{ \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  integer                ::  n_relax_down, n_relax_up, n_relax_bottom
  !> \}
  
  !> \brief implizite Richtungen bei Linienrelaxation (Multigrid)
  integer                ::  impl_dir(1:3)
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  logical                ::  weighting_yes
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsstatistik
  !===========================================================================================================
  real                   ::  dtime_average
  real                   ::  max_div_init(1:2)
  integer                ::  number_poisson
  
  !> \{ \brief Zähler
  integer                ::  countO(1:RK_steps)
  integer                ::  countP(1:RK_steps,1:2)
  integer                ::  countH(1:RK_steps,1:3)
  !> \}
  
  !> \{ \brief Konvergenzrate
  real                   ::  ratioO(1:RK_steps)
  real                   ::  ratioH(1:RK_steps,1:3)
  real                   ::  ratioP(1:RK_steps,1:2)
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name MPI
  !===========================================================================================================
  !> \{ \brief Kommunikatoren
  integer                ::  COMM_CART
  
  integer                ::  COMM_SLICE1, COMM_BAR1
  integer                ::  COMM_SLICE2, COMM_BAR2
  integer                ::  COMM_SLICE3, COMM_BAR3
  !> \}
  
  !> \{ \brief Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes)
  !!
  !! (for MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  integer, allocatable   ::  bar1_size(:), bar1_offset(:)
  integer, allocatable   ::  bar2_size(:), bar2_offset(:)
  integer, allocatable   ::  bar3_size(:), bar3_offset(:)
  !> \}
  
  !> \{ \brief Ränge der Prozesse
  integer                ::  rank
  integer                ::  rank_bar1, rank_slice1
  integer                ::  rank_bar2, rank_slice2
  integer                ::  rank_bar3, rank_slice3
  !> \}
  
  !> \{ \brief Ränge der Nachbarprozesse (in kartesischem Gitter)
  integer                ::  rank1L, rank1U
  integer                ::  rank2L, rank2U
  integer                ::  rank3L, rank3U
  !> \}
  
  !> \brief Error-Handle
  integer                ::  merror
  
  !> \{ \brief Request-Handles
  !!
  !! (müssen offenbar global angelegt werden)
  integer                ::  req1L, req1U
  integer                ::  req2L, req2U
  integer                ::  req3L, req3U
  !> \}
  
  !> \}
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  integer                ::  herror
  
  !> \{ \name periodicity criterions
  real   , allocatable   ::  velp(:,:,:,:)       !< \brief stored velocity, used as reference to previous period
  real                   ::  tp
  real                   ::  freq               !< \brief frequency
  real                   ::  periodic_tol       !< \brief tolerance which, we use as stopping criterion. \f[ ||u(t)-u(t-T)|| < tol \f]
  !> \}
  
  
#else
  
  
  
  
  !===========================================================================================================
  !> \brief Raemliche Dimensionen
  !===========================================================================================================
  !> 2D wird ueber M3==2 eingeschaltet
  integer                ::  dimens
  
  
  !===========================================================================================================
  !> \{ \name Domain- und Blockspezifikationen
  !!===========================================================================================================
  !!--- zulaessige Blockgroessen ------------------------------------------------------------------------------
  !!--------------------------------------------------------------
  !! n  x 1   2   4   8   16   32   64   128   256  512   ...
  !!--------------------------------------------------------------
  !! n2 = 2,  4,  8, 16,  32,  64, 128,  256,  512, 1024, ...
  !! n3 = 3,  6, 12, 24,  48,  96, 192,  384,  768, 1536, ...
  !! n5 = 5, 10, 20, 40,  80, 160, 320,  640, 1280, 2560, ...
  !! n7 = 7, 14, 28, 56, 112, 224, 448,  896, 1792, 3584, ...
  !! n9 = 9, 18, 36, 72, 144, 288, 576, 1152, 2304, 4608, ...
  !! ...
  integer, parameter     ::  N1 = 1+(M1-1)/NB1
  integer, parameter     ::  N2 = 1+(M2-1)/NB2
  integer, parameter     ::  N3 = 1+(M3-1)/NB3
  
  !> \{ \brief Anzahl grobe Gitter (Multigrid)
  integer, parameter     ::  n_grids_max = 15
  integer                ::  n_grids, n_grids_limit
  !> \}
  
  !> \{ \brief Dimensionen
  integer, parameter     ::  dim_ncb1c = SIZE(ncb1c)
  integer, parameter     ::  dim_ncb1g = SIZE(ncb1g)
  integer, parameter     ::  dim_ncb1d = SIZE(ncb1d)
  
  integer, parameter     ::  dim_ncb2c = SIZE(ncb2c)
  integer, parameter     ::  dim_ncb2g = SIZE(ncb2g)
  integer, parameter     ::  dim_ncb2d = SIZE(ncb2d)
                         
  integer, parameter     ::  dim_ncb3c = SIZE(ncb3c)
  integer, parameter     ::  dim_ncb3g = SIZE(ncb3g)
  integer, parameter     ::  dim_ncb3d = SIZE(ncb3d)
  !> \}
  
  !> \{ \brief Anzahl Stencil-Koeffizienten (Feld)
  !!
  !! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  integer, parameter     ::  nc1c = ncb1c(dim_ncb1c)
  integer, parameter     ::  nc1s = ncb1g(dim_ncb1g)
                         
  integer, parameter     ::  nc2c = ncb2c(dim_ncb2c)
  integer, parameter     ::  nc2s = ncb2g(dim_ncb2g)
                         
  integer, parameter     ::  nc3c = ncb3c(dim_ncb3c)
  integer, parameter     ::  nc3s = ncb3g(dim_ncb3g)
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief zentral
  integer, parameter     ::  b1U = nc1s/2
  integer, parameter     ::  b2U = nc2s/2
  integer, parameter     ::  b3U = nc3s/2
  
  integer, parameter     ::  b1L = -b1U
  integer, parameter     ::  b2L = -b2U
  integer, parameter     ::  b3L = -b3U
  !> \}
  
  !> \{ \brief upwind (nicht-linear).
  !! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !!  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  integer, parameter     ::  n1L = b1L
  integer, parameter     ::  n2L = b2L
  integer, parameter     ::  n3L = b3L
  
  integer, parameter     ::  n1U = b1U
  integer, parameter     ::  n2U = b2U
  integer, parameter     ::  n3U = b3U
  !> \}
  
  !> \{ \brief Divergenz
  integer, parameter     ::  d1L = b1L
  integer, parameter     ::  d2L = b2L
  integer, parameter     ::  d3L = b3L
  
  integer, parameter     ::  d1U = b1U-1
  integer, parameter     ::  d2U = b2U-1
  integer, parameter     ::  d3U = b3U-1
  !> \}
  
  !> \{ \brief Gradient
  integer, parameter     ::  g1L = b1L+1
  integer, parameter     ::  g2L = b2L+1
  integer, parameter     ::  g3L = b3L+1
  
  integer, parameter     ::  g1U = b1U
  integer, parameter     ::  g2U = b2U
  integer, parameter     ::  g3U = b3U
  !> \}
  
  
  !===========================================================================================================
  !> \{ \name Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief 1. Ableitung (zentral)
  real                   ::  cp1  (b1L:b1U,0:N1)
  real                   ::  cp2  (b2L:b2U,0:N2)
  real                   ::  cp3  (b3L:b3U,0:N3)
  
  real                   ::  cu1  (b1L:b1U,0:N1)
  real                   ::  cv2  (b2L:b2U,0:N2)
  real                   ::  cw3  (b3L:b3U,0:N3)
  !> \}
  
  !> \{ Eigene Stencils. notwendig, da
  !!   - Symmetrie des Geschwindigkeitsfeldes nicht unbedingt Symmetrie bei Konzentrationen bedeutet,
  !!   - Randbedinungen darauf gespeichert werden koennen.
  real                   ::  cc1  (b1L:b1U,0:N1,1:n_conc)
  real                   ::  cc2  (b2L:b2U,0:N2,1:n_conc)
  real                   ::  cc3  (b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief 1. Ableitung (upwind)
  real                   ::  cNp1D(n1L:n1U,0:N1)
  real                   ::  cNp2D(n2L:n2U,0:N2)
  real                   ::  cNp3D(n3L:n3U,0:N3)
  
  real                   ::  cNp1U(n1L:n1U,0:N1)
  real                   ::  cNp2U(n2L:n2U,0:N2)
  real                   ::  cNp3U(n3L:n3U,0:N3)
  
  real                   ::  cNu1D(n1L:n1U,0:N1)
  real                   ::  cNv2D(n2L:n2U,0:N2)
  real                   ::  cNw3D(n3L:n3U,0:N3)
  
  real                   ::  cNu1U(n1L:n1U,0:N1)
  real                   ::  cNv2U(n2L:n2U,0:N2)
  real                   ::  cNw3U(n3L:n3U,0:N3)
  
  real                   ::  cNc1D(n1L:n1U,0:N1,1:n_conc)
  real                   ::  cNc2D(n2L:n2U,0:N2,1:n_conc)
  real                   ::  cNc3D(n3L:n3U,0:N3,1:n_conc)
  
  real                   ::  cNc1U(n1L:n1U,0:N1,1:n_conc)
  real                   ::  cNc2U(n2L:n2U,0:N2,1:n_conc)
  real                   ::  cNc3U(n3L:n3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Divergenz
  real                   ::  cDu1 (d1L:d1U,0:N1)
  real                   ::  cDv2 (d2L:d2U,0:N2)
  real                   ::  cDw3 (d3L:d3U,0:N3)
  !> \}
  
  !> \{ \brief Divergenz (transponiert)
  real                   ::  cDu1T(g1L:g1U,0:N1)
  real                   ::  cDv2T(g2L:g2U,0:N2)
  real                   ::  cDw3T(g3L:g3U,0:N3)
  !> \}
  
  !> \{ \brief Gradient
  real                   ::  cGp1 (g1L:g1U,0:N1)
  real                   ::  cGp2 (g2L:g2U,0:N2)
  real                   ::  cGp3 (g3L:g3U,0:N3)
  !> \}
  
  !> \{ \brief Gradient (transponiert)
  real                   ::  cGp1T(d1L:d1U,0:N1)
  real                   ::  cGp2T(d2L:d2U,0:N2)
  real                   ::  cGp3T(d3L:d3U,0:N3)
  !> \}
  
  !> \{ \brief 2. Ableitung (zentral)
  real                   ::  cp11 (b1L:b1U,0:N1)
  real                   ::  cp22 (b2L:b2U,0:N2)
  real                   ::  cp33 (b3L:b3U,0:N3)
  
  real                   ::  cu11 (b1L:b1U,0:N1)
  real                   ::  cv22 (b2L:b2U,0:N2)
  real                   ::  cw33 (b3L:b3U,0:N3)
  
  real                   ::  cc11 (b1L:b1U,0:N1,1:n_conc)
  real                   ::  cc22 (b2L:b2U,0:N2,1:n_conc)
  real                   ::  cc33 (b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Interpolation
  real                   ::  cIpu(g1L:g1U,0:N1)
  real                   ::  cIpv(g2L:g2U,0:N2)
  real                   ::  cIpw(g3L:g3U,0:N3)
  
  real                   ::  cIup(d1L:d1U,0:N1)
  real                   ::  cIvp(d2L:d2U,0:N2)
  real                   ::  cIwp(d3L:d3U,0:N3)
  
  real                   ::  cIcu(g1L:g1U,0:N1,1:n_conc)
  real                   ::  cIcv(g2L:g2U,0:N2,1:n_conc)
  real                   ::  cIcw(g3L:g3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Filter
  real                   ::  cFp1(b1L:b1U,0:N1)
  real                   ::  cFp2(b2L:b2U,0:N2)
  real                   ::  cFp3(b3L:b3U,0:N3)
  
  real                   ::  cFu1(b1L:b1U,0:N1)
  real                   ::  cFv2(b2L:b2U,0:N2)
  real                   ::  cFw3(b3L:b3U,0:N3)
  
  real                   ::  cFc1(b1L:b1U,0:N1,1:n_conc)
  real                   ::  cFc2(b2L:b2U,0:N2,1:n_conc)
  real                   ::  cFc3(b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief --- Integrator (only for Druckgitter)
  real                   ::  cInt1(b1L:b1U,0:N1)
  real                   ::  cInt2(b2L:b2U,0:N2)
  real                   ::  cInt3(b3L:b3U,0:N3)
  !> \}
  
  
  !> \{ \brief Compact
  integer, parameter     ::  dimS1 = 2*ndL*NB1
  integer, parameter     ::  dimS2 = 2*ndL*NB2
  integer, parameter     ::  dimS3 = 2*ndL*NB3
  
  real                   ::  buffer1A(0:N2,0:N3,1:2*ndL), buffer1B(0:N2,0:N3,1:dimS1)
  real                   ::  buffer2A(0:N1,0:N3,1:2*ndL), buffer2B(0:N1,0:N3,1:dimS2)
  real                   ::  buffer3A(0:N1,0:N2,1:2*ndL), buffer3B(0:N1,0:N2,1:dimS3)
  !----------------------------------------------------
  integer                ::  disp1(1:NB1), recv1(1:NB1)
  integer                ::  disp2(1:NB2), recv2(1:NB2)
  integer                ::  disp3(1:NB3), recv3(1:NB3)
  !----------------------------------------------------
  !> \}
  
  !> \{ \brief Schur complement (square) matrices
  real                   ::  Sp1(1:dimS1,1:dimS1), Su1(1:dimS1,1:dimS1), Sc1(1:dimS1,1:dimS1,1:n_conc), Sp11(1:dimS1,1:dimS1), Su11(1:dimS1,1:dimS1), Sc11(1:dimS1,1:dimS1,1:n_conc)
  real                   ::  Sp2(1:dimS2,1:dimS2), Sv2(1:dimS2,1:dimS2), Sc2(1:dimS2,1:dimS2,1:n_conc), Sp22(1:dimS2,1:dimS2), Sv22(1:dimS2,1:dimS2), Sc22(1:dimS2,1:dimS2,1:n_conc)
  real                   ::  Sp3(1:dimS3,1:dimS3), Sw3(1:dimS3,1:dimS3), Sc3(1:dimS3,1:dimS3,1:n_conc), Sp33(1:dimS3,1:dimS3), Sw33(1:dimS3,1:dimS3), Sc33(1:dimS3,1:dimS3,1:n_conc)
  
  real                   ::  SDu1(1:dimS1,1:dimS1), SGp1(1:dimS1,1:dimS1), SIpu(1:dimS1,1:dimS1), SIup(1:dimS1,1:dimS1), SIcu(1:dimS1,1:dimS1,1:n_conc)
  real                   ::  SDv2(1:dimS2,1:dimS2), SGp2(1:dimS2,1:dimS2), SIpv(1:dimS2,1:dimS2), SIvp(1:dimS2,1:dimS2), SIcv(1:dimS2,1:dimS2,1:n_conc)
  real                   ::  SDw3(1:dimS3,1:dimS3), SGp3(1:dimS3,1:dimS3), SIpw(1:dimS3,1:dimS3), SIwp(1:dimS3,1:dimS3), SIcw(1:dimS3,1:dimS3,1:n_conc)
  
  real                   ::  SDu1T(1:dimS1,1:dimS1), SGp1T(1:dimS1,1:dimS1)
  real                   ::  SDv2T(1:dimS2,1:dimS2), SGp2T(1:dimS2,1:dimS2)
  real                   ::  SDw3T(1:dimS3,1:dimS3), SGp3T(1:dimS3,1:dimS3)
  !> \}
  
  !> \{ \brief bottom (horizontally oriented) matrices
  real                   ::  Wp1(1:2*ndL,0:N1), Wu1(1:2*ndL,0:N1), Wc1(1:2*ndL,0:N1,1:n_conc), Wp11(1:2*ndL,0:N1), Wu11(1:2*ndL,0:N1), Wc11(1:2*ndL,0:N1,1:n_conc)
  real                   ::  Wp2(1:2*ndL,0:N2), Wv2(1:2*ndL,0:N2), Wc2(1:2*ndL,0:N2,1:n_conc), Wp22(1:2*ndL,0:N2), Wv22(1:2*ndL,0:N2), Wc22(1:2*ndL,0:N2,1:n_conc)
  real                   ::  Wp3(1:2*ndL,0:N3), Ww3(1:2*ndL,0:N3), Wc3(1:2*ndL,0:N3,1:n_conc), Wp33(1:2*ndL,0:N3), Ww33(1:2*ndL,0:N3), Wc33(1:2*ndL,0:N3,1:n_conc)
  
  real                   ::  WDu1(1:2*ndL,0:N1), WGp1(1:2*ndL,0:N1), WIpu(1:2*ndL,0:N1), WIup(1:2*ndL,0:N1), WIcu(1:2*ndL,0:N1,1:n_conc)
  real                   ::  WDv2(1:2*ndL,0:N2), WGp2(1:2*ndL,0:N2), WIpv(1:2*ndL,0:N2), WIvp(1:2*ndL,0:N2), WIcv(1:2*ndL,0:N2,1:n_conc)
  real                   ::  WDw3(1:2*ndL,0:N3), WGp3(1:2*ndL,0:N3), WIpw(1:2*ndL,0:N3), WIwp(1:2*ndL,0:N3), WIcw(1:2*ndL,0:N3,1:n_conc)
  
  real                   ::  WDu1T(1:2*ndL,0:N1), WGp1T(1:2*ndL,0:N1)
  real                   ::  WDv2T(1:2*ndL,0:N2), WGp2T(1:2*ndL,0:N2)
  real                   ::  WDw3T(1:2*ndL,0:N3), WGp3T(1:2*ndL,0:N3)
  !> \}
  
  !> \{ \brief implicit finite difference coefficients ("L"=left-hand side)
  real                   ::  cp1CL (-ndL:ndL,0:N1), cu1CL(-ndL:ndL,0:N1), cc1CL(-ndL:ndL,0:N1,1:n_conc), cp11CL(-ndL:ndL,0:N1), cu11CL(-ndL:ndL,0:N1), cc11CL(-ndL:ndL,0:N1,1:n_conc)
  real                   ::  cp2CL (-ndL:ndL,0:N2), cv2CL(-ndL:ndL,0:N2), cc2CL(-ndL:ndL,0:N2,1:n_conc), cp22CL(-ndL:ndL,0:N2), cv22CL(-ndL:ndL,0:N2), cc22CL(-ndL:ndL,0:N2,1:n_conc)
  real                   ::  cp3CL (-ndL:ndL,0:N3), cw3CL(-ndL:ndL,0:N3), cc3CL(-ndL:ndL,0:N3,1:n_conc), cp33CL(-ndL:ndL,0:N3), cw33CL(-ndL:ndL,0:N3), cc33CL(-ndL:ndL,0:N3,1:n_conc)
  
  real                   ::  cDu1CL(-ndL:ndL,0:N1), cGp1CL(-ndL:ndL,0:N1), cIpuCL(-ndL:ndL,0:N1), cIupCL(-ndL:ndL,0:N1), cIcuCL(-ndL:ndL,0:N1,1:n_conc)
  real                   ::  cDv2CL(-ndL:ndL,0:N2), cGp2CL(-ndL:ndL,0:N2), cIpvCL(-ndL:ndL,0:N2), cIvpCL(-ndL:ndL,0:N2), cIcvCL(-ndL:ndL,0:N2,1:n_conc)
  real                   ::  cDw3CL(-ndL:ndL,0:N3), cGp3CL(-ndL:ndL,0:N3), cIpwCL(-ndL:ndL,0:N3), cIwpCL(-ndL:ndL,0:N3), cIcwCL(-ndL:ndL,0:N3,1:n_conc)
  
  real                   ::  cDu1CLT(-ndL:ndL,0:N1), cGp1CLT(-ndL:ndL,0:N1)
  real                   ::  cDv2CLT(-ndL:ndL,0:N2), cGp2CLT(-ndL:ndL,0:N2)
  real                   ::  cDw3CLT(-ndL:ndL,0:N3), cGp3CLT(-ndL:ndL,0:N3)
  
  real                   ::  cp1CL_LU(1:(3*ndL+1),0:N1), cu1CL_LU(1:(3*ndL+1),0:N1), cc1CL_LU(1:(3*ndL+1),0:N1,1:n_conc), cp11CL_LU(1:(3*ndL+1),0:N1), cu11CL_LU(1:(3*ndL+1),0:N1), cc11CL_LU(1:(3*ndL+1),0:N1,1:n_conc)
  real                   ::  cp2CL_LU(1:(3*ndL+1),0:N2), cv2CL_LU(1:(3*ndL+1),0:N2), cc2CL_LU(1:(3*ndL+1),0:N2,1:n_conc), cp22CL_LU(1:(3*ndL+1),0:N2), cv22CL_LU(1:(3*ndL+1),0:N2), cc22CL_LU(1:(3*ndL+1),0:N2,1:n_conc)
  real                   ::  cp3CL_LU(1:(3*ndL+1),0:N3), cw3CL_LU(1:(3*ndL+1),0:N3), cc3CL_LU(1:(3*ndL+1),0:N3,1:n_conc), cp33CL_LU(1:(3*ndL+1),0:N3), cw33CL_LU(1:(3*ndL+1),0:N3), cc33CL_LU(1:(3*ndL+1),0:N3,1:n_conc)
  
  real                   ::  cDu1CL_LU(1:(3*ndL+1),0:N1), cGp1CL_LU(1:(3*ndL+1),0:N1), cIpuCL_LU(1:(3*ndL+1),0:N1), cIupCL_LU(1:(3*ndL+1),0:N1), cIcuCL_LU(1:(3*ndL+1),0:N1,1:n_conc)
  real                   ::  cDv2CL_LU(1:(3*ndL+1),0:N2), cGp2CL_LU(1:(3*ndL+1),0:N2), cIpvCL_LU(1:(3*ndL+1),0:N2), cIvpCL_LU(1:(3*ndL+1),0:N2), cIcvCL_LU(1:(3*ndL+1),0:N2,1:n_conc)
  real                   ::  cDw3CL_LU(1:(3*ndL+1),0:N3), cGp3CL_LU(1:(3*ndL+1),0:N3), cIpwCL_LU(1:(3*ndL+1),0:N3), cIwpCL_LU(1:(3*ndL+1),0:N3), cIcwCL_LU(1:(3*ndL+1),0:N3,1:n_conc)
  
  real                   ::  cDu1CLT_LU(1:(3*ndL+1),0:N1), cGp1CLT_LU(1:(3*ndL+1),0:N1)
  real                   ::  cDv2CLT_LU(1:(3*ndL+1),0:N2), cGp2CLT_LU(1:(3*ndL+1),0:N2)
  real                   ::  cDw3CLT_LU(1:(3*ndL+1),0:N3), cGp3CLT_LU(1:(3*ndL+1),0:N3)
  !> \}
  
  !> \{ \brief explicit finite difference coefficients ("R"=right-hand side)
  real                   ::  cp1CR(-ndR:ndR,0:N1), cu1CR(-ndR:ndR,0:N1), cc1CR(-ndR:ndR,0:N1,1:n_conc), cp11CR(-ndR:ndR,0:N1), cu11CR(-ndR:ndR,0:N1), cc11CR(-ndR:ndR,0:N1,1:n_conc)
  real                   ::  cp2CR(-ndR:ndR,0:N2), cv2CR(-ndR:ndR,0:N2), cc2CR(-ndR:ndR,0:N2,1:n_conc), cp22CR(-ndR:ndR,0:N2), cv22CR(-ndR:ndR,0:N2), cc22CR(-ndR:ndR,0:N2,1:n_conc)
  real                   ::  cp3CR(-ndR:ndR,0:N3), cw3CR(-ndR:ndR,0:N3), cc3CR(-ndR:ndR,0:N3,1:n_conc), cp33CR(-ndR:ndR,0:N3), cw33CR(-ndR:ndR,0:N3), cc33CR(-ndR:ndR,0:N3,1:n_conc)
  
  real                   ::  cDu1CR(-ndR:ndR,0:N1), cGp1CR(-ndR:ndR,0:N1), cIpuCR(-ndR:ndR,0:N1), cIupCR(-ndR:ndR,0:N1), cIcuCR(-ndR:ndR,0:N1,1:n_conc)
  real                   ::  cDv2CR(-ndR:ndR,0:N2), cGp2CR(-ndR:ndR,0:N2), cIpvCR(-ndR:ndR,0:N2), cIvpCR(-ndR:ndR,0:N2), cIcvCR(-ndR:ndR,0:N2,1:n_conc)
  real                   ::  cDw3CR(-ndR:ndR,0:N3), cGp3CR(-ndR:ndR,0:N3), cIpwCR(-ndR:ndR,0:N3), cIwpCR(-ndR:ndR,0:N3), cIcwCR(-ndR:ndR,0:N3,1:n_conc)
  
  real                   ::  cDu1CRT(-ndR:ndR,0:N1), cGp1CRT(-ndR:ndR,0:N1)
  real                   ::  cDv2CRT(-ndR:ndR,0:N2), cGp2CRT(-ndR:ndR,0:N2)
  real                   ::  cDw3CRT(-ndR:ndR,0:N3), cGp3CRT(-ndR:ndR,0:N3)
  !> \}
  !****************************************************
  
  
  !> \{ \brief 2. Ableitung (Multigrid).
  !! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  real                   ::  cp11R(-1:1,0:N1,1:n_grids_max)
  real                   ::  cp22R(-1:1,0:N2,1:n_grids_max)
  real                   ::  cp33R(-1:1,0:N3,1:n_grids_max)
  
  real                   ::  cu11R(-1:1,0:N1,1:n_grids_max)
  real                   ::  cv22R(-1:1,0:N2,1:n_grids_max)
  real                   ::  cw33R(-1:1,0:N3,1:n_grids_max)
  
  real                   ::  cc11R(-1:1,0:N1,1:n_grids_max,1:n_conc)
  real                   ::  cc22R(-1:1,0:N2,1:n_grids_max,1:n_conc)
  real                   ::  cc33R(-1:1,0:N3,1:n_grids_max,1:n_conc)
  
  real                   ::  cdg1 (-1:1,1:N1,1:n_grids_max)
  real                   ::  cdg2 (-1:1,1:N2,1:n_grids_max)
  real                   ::  cdg3 (-1:1,1:N3,1:n_grids_max)
  !> \}
  
  !> \{ \brief Interpolation (Multigrid)
  real                   ::  cI1(1:2,1:N1,1:n_grids_max)
  real                   ::  cI2(1:2,1:N2,1:n_grids_max)
  real                   ::  cI3(1:2,1:N3,1:n_grids_max)
  
  real                   ::  cIH1(1:2,0:N1)
  real                   ::  cIH2(1:2,0:N2)
  real                   ::  cIH3(1:2,0:N3)
  !> \}
  
  !> \{ \brief Restriktion (Multigrid)
  real                   ::  cR1 (-1:1,1:N1,2:n_grids_max)
  real                   ::  cR2 (-1:1,1:N2,2:n_grids_max)
  real                   ::  cR3 (-1:1,1:N3,2:n_grids_max)
  
  real                   ::  cRest1(b1L:b1U,0:N1,1:n_grids_max-1) ! TEST!!!
  real                   ::  cRest2(b2L:b2U,0:N2,1:n_grids_max-1)
  real                   ::  cRest3(b3L:b3U,0:N3,1:n_grids_max-1)
  
  real                   ::  cRH1(1:2,1:N1)
  real                   ::  cRH2(1:2,1:N2)
  real                   ::  cRH3(1:2,1:N3)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Gitterspezifikationen
  !===========================================================================================================
  !> \{ \brief physiklische Koordinaten (global)
  real                   ::  y1p(1:M1), y1u(0:M1)
  real                   ::  y2p(1:M2), y2v(0:M2)
  real                   ::  y3p(1:M3), y3w(0:M3)
  !> \}
  
  !> \{ \brief physiklische Koordinaten (Block)
  real                   ::  x1p(b1L:(N1+b1U)), x1u(b1L:(N1+b1U))
  real                   ::  x2p(b2L:(N2+b2U)), x2v(b2L:(N2+b2U))
  real                   ::  x3p(b3L:(N3+b3U)), x3w(b3L:(N3+b3U))
  !> \}
  
  !> \{ \brief --- physiklische Koordinaten (Block, Multigrid)
  real                   ::  x1pR(b1L:(N1+b1U),1:n_grids_max), x1uR(b1L:(N1+b1U),1:n_grids_max)
  real                   ::  x2pR(b2L:(N2+b2U),1:n_grids_max), x2vR(b2L:(N2+b2U),1:n_grids_max)
  real                   ::  x3pR(b3L:(N3+b3U),1:n_grids_max), x3wR(b3L:(N3+b3U),1:n_grids_max)
  !> \}
  
  !> \{ \brief Gitterweiten (global)
  real                   ::  dy1p(1:M1), dy1u(0:M1)
  real                   ::  dy2p(1:M2), dy2v(0:M2)
  real                   ::  dy3p(1:M3), dy3w(0:M3)
  !> \}
  
  !> \{ \brief Gitterweiten (Block)
  real                   ::  dx1p(1:N1), dx1u(0:N1)
  real                   ::  dx2p(1:N2), dx2v(0:N2)
  real                   ::  dx3p(1:N3), dx3w(0:N3)
  
  real                   ::  dx1DM(1:N1), dx1pM(1:N1), ddx1pM(1:N1)
  real                   ::  dx2DM(1:N2), dx2pM(1:N2), ddx2pM(1:N2)
  real                   ::  dx3DM(1:N3), dx3pM(1:N3), ddx3pM(1:N3)
  
  real                   ::  dx1GM(0:N1), dx1uM(0:N1), ddx1uM(0:N1)
  real                   ::  dx2GM(0:N2), dx2vM(0:N2), ddx2vM(0:N2)
  real                   ::  dx3GM(0:N3), dx3wM(0:N3), ddx3wM(0:N3)
  !> \}
  
  !> \{ \brief Smagorinsky-Modell
  real                   ::  dx1pS(1:N1), dx1uS(0:N1) ! TEST!!!
  real                   ::  dx2pS(1:N2), dx2vS(0:N2)
  real                   ::  dx3pS(1:N3), dx3wS(0:N3)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Arbeitsfelder
  !===========================================================================================================
  !> \brief Geschwindigkeiten -------------------------------------------------------------------------------------
  real                   ::  vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief nicht-linearer Term -----------------------------------------------------------------------------------
  real                   ::  nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief Recht-Hand-Seite --------------------------------------------------------------------------------------
  real                   ::  rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief Druck -------------------------------------------------------------------------------------------------
  real                   ::  pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Konzentrationsfelder ----------------------------------------------------------------------------------
  real                   ::  conc(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc)
  
  !> \brief nicht-linearer Term (Konzentration) -------------------------------------------------------------------
  real                   ::  nlco(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc)
  
  !> \{ \brief Ausfluss-RB (Geschwindigkeitsfeld)
  !!
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehürigen Randbedingungen gespeichert werden.
  real                   ::  bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2)
  real                   ::  bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2)
  real                   ::  bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2)
  
  real                   ::  bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2)
  real                   ::  bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2)
  real                   ::  bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2)
  
  real                   ::  bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2)
  real                   ::  bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2)
  real                   ::  bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2)
  
  real                   ::  drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2)
  real                   ::  drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2)
  real                   ::  drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2)
  !>\}

  !> \{ \brief Sedimentations-RB (Konzentration)
  !!
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehürigen Randbedingungen gespeichert werden.
  real                   ::  sed1L (1:N2,1:N3,1:n_conc)
  real                   ::  sed1U (1:N2,1:N3,1:n_conc)
  
  real                   ::  sed2L (1:N1,1:N3,1:n_conc)
  real                   ::  sed2U (1:N1,1:N3,1:n_conc)
  
  real                   ::  sed3L (1:N1,1:N2,1:n_conc)
  real                   ::  sed3U (1:N1,1:N2,1:n_conc)
  
  real                   ::  conc1L(1:N2,1:N3,1:n_conc)
  real                   ::  conc1U(1:N2,1:N3,1:n_conc)
  
  real                   ::  conc2L(1:N1,1:N3,1:n_conc)
  real                   ::  conc2U(1:N1,1:N3,1:n_conc)
  
  real                   ::  conc3L(1:N1,1:N2,1:n_conc)
  real                   ::  conc3U(1:N1,1:N2,1:n_conc)
  !> \}
  
  !> \{ \brief Sedimentation (Konzentration, Ausschrieb)
  real                   ::  dep1L_conc(1:N2,1:N3,1:n_conc)
  real                   ::  dep1U_conc(1:N2,1:N3,1:n_conc)
  
  real                   ::  dep2L_conc(1:N1,1:N3,1:n_conc)
  real                   ::  dep2U_conc(1:N1,1:N3,1:n_conc)
  
  real                   ::  dep3L_conc(1:N1,1:N2,1:n_conc)
  real                   ::  dep3U_conc(1:N1,1:N2,1:n_conc)
  !> \}
  
  !> \{ \brief Sedimentation (Partikel, Ausschrieb).
  !! TEST!!! 1:n_spec???
  real                   ::  dep1L_part(1:N2,1:N3)
  real                   ::  dep1U_part(1:N2,1:N3)
  
  real                   ::  dep2L_part(1:N1,1:N3)
  real                   ::  dep2U_part(1:N1,1:N3)
  
  real                   ::  dep3L_part(1:N1,1:N2)
  real                   ::  dep3U_part(1:N1,1:N2)
  !> \}

  !> \brief Residuum ----------------------------------------------------------------------------------------------
  real                   ::  res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Druckgradient (eine Komponente) -----------------------------------------------------------------------
  real                   ::  gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Gewichte for Divergenzfreiheit ------------------------------------------------------------------------
  real                   ::  weight(1:N1,1:N2,1:N3)
  
#ifdef NONBOUSSINESQ
  !> \brief fluid density on velocity grid ------------------------------------------------------------------------
  real                   ::  dens(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
#endif
  
  !> \{ \brief Null-Raum-Vektor
  real                   ::  psi      (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  psi_vel  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  real                   ::  psi_rel1 (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real   , allocatable   ::  psi_rel2 (:,:,:)
  real   , allocatable   ::  psi_rel3 (:,:,:)
  real   , allocatable   ::  psi_rel4 (:,:,:)
  real   , allocatable   ::  psi_rel5 (:,:,:)
  real   , allocatable   ::  psi_rel6 (:,:,:)
  real   , allocatable   ::  psi_rel7 (:,:,:)
  real   , allocatable   ::  psi_rel8 (:,:,:)
  real   , allocatable   ::  psi_rel9 (:,:,:)
  real   , allocatable   ::  psi_rel10(:,:,:)
  real   , allocatable   ::  psi_rel11(:,:,:)
  real   , allocatable   ::  psi_rel12(:,:,:)
  real   , allocatable   ::  psi_rel13(:,:,:)
  real   , allocatable   ::  psi_rel14(:,:,:)
  real   , allocatable   ::  psi_rel15(:,:,:)
  
  real                   ::  th11(1:N2,1:N3,1:2)
  real                   ::  th12(0:N1,1:N3,1:2)
  real                   ::  th13(0:N1,1:N2,1:2)
  
  real                   ::  th21(0:N2,1:N3,1:2)
  real                   ::  th22(1:N1,1:N3,1:2)
  real                   ::  th23(1:N1,0:N2,1:2)
  
  real                   ::  th31(1:N2,0:N3,1:2)
  real                   ::  th32(1:N1,0:N3,1:2)
  real                   ::  th33(1:N1,1:N2,1:2)
  !> \}
  
  !> \{ \brief Multigrid
  real                   ::  vec1C (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  real   , allocatable   ::  vec2A (:,:,:)
  real   , allocatable   ::  vec2B (:,:,:)
  real   , allocatable   ::  vec2C (:,:,:)
  
  real   , allocatable   ::  vec3A (:,:,:)
  real   , allocatable   ::  vec3B (:,:,:)
  real   , allocatable   ::  vec3C (:,:,:)
  
  real   , allocatable   ::  vec4A (:,:,:)
  real   , allocatable   ::  vec4B (:,:,:)
  real   , allocatable   ::  vec4C (:,:,:)
  
  real   , allocatable   ::  vec5A (:,:,:)
  real   , allocatable   ::  vec5B (:,:,:)
  real   , allocatable   ::  vec5C (:,:,:)
  
  real   , allocatable   ::  vec6A (:,:,:)
  real   , allocatable   ::  vec6B (:,:,:)
  real   , allocatable   ::  vec6C (:,:,:)
  
  real   , allocatable   ::  vec7A (:,:,:)
  real   , allocatable   ::  vec7B (:,:,:)
  real   , allocatable   ::  vec7C (:,:,:)
  
  real   , allocatable   ::  vec8A (:,:,:)
  real   , allocatable   ::  vec8B (:,:,:)
  real   , allocatable   ::  vec8C (:,:,:)
  
  real   , allocatable   ::  vec9A (:,:,:)
  real   , allocatable   ::  vec9B (:,:,:)
  real   , allocatable   ::  vec9C (:,:,:)
  
  real   , allocatable   ::  vec10A(:,:,:)
  real   , allocatable   ::  vec10B(:,:,:)
  real   , allocatable   ::  vec10C(:,:,:)
  
  real   , allocatable   ::  vec11A(:,:,:)
  real   , allocatable   ::  vec11B(:,:,:)
  real   , allocatable   ::  vec11C(:,:,:)
  
  real   , allocatable   ::  vec12A(:,:,:)
  real   , allocatable   ::  vec12B(:,:,:)
  real   , allocatable   ::  vec12C(:,:,:)
  
  real   , allocatable   ::  vec13A(:,:,:)
  real   , allocatable   ::  vec13B(:,:,:)
  real   , allocatable   ::  vec13C(:,:,:)
  
  real   , allocatable   ::  vec14A(:,:,:)
  real   , allocatable   ::  vec14B(:,:,:)
  real   , allocatable   ::  vec14C(:,:,:)
  
  real   , allocatable   ::  vec15A(:,:,:)
  real   , allocatable   ::  vec15B(:,:,:)
  !> \}
  
  
  !> \{ \brief BiCGstab / Richardson
  real                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !> \}
  
  !> \brief product_div_grad --------------------------------------------------------------------------------------
  real                   ::  dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Hilfsfeld (kompakte Differenzen) ----------------------------------------------------------------------
  real                   ::  com(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \{ \brief Hilfsfelder (Druckiteration).
  !! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  !! - wird auch fuer interpolierte Geschwindigkeiten in rhs_vel und rhs_conc verwendet.
  real                   ::  work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  real                   ::  work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !> \}
  
  !> \{ \brief Lagrange-Partikel
  integer                ::  n_part, n_part_tot, n_groups
  real                   ::  particles(1:n_args,1:n_part_max)
  !> \}
  
  !> \{ \brief Smagorinsky LES-Modell
  real   , allocatable   ::  smag_diffus(:,:,:,:)
  !> \}
  
  !> \{ \brief Linienrelaxation
  real                   ::  vec1(1:N1), dia1(1:N1), SOR1(1:N1), band1(1:2,1:N1) ! TEST!!! vec1, dia muessen auch in mod_Helmholtz ausgetauscht werden!
  real                   ::  vec2(1:N2), dia2(1:N2), SOR2(1:N2), band2(1:2,1:N2) ! TEST!!! SOR muesste man idealerweise auch in band eingliedern!
  real                   ::  vec3(1:N3), dia3(1:N3), SOR3(1:N3), band3(1:2,1:N3) ! dia3(:) --> band3(1,:), vec3(:) --> band3(2,:), etc. ...
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \brief Indizierung (Intervallgrenzen, Verschiebungen)
  !===========================================================================================================
  !> \{ \brief Block-Index
  !INTEGER                ::  iB, jB, kB ! TEST!!!
  
  !> \} \brief Indexverschiebung (Block --> global)
  integer                ::  iShift, jShift, kShift
  !> \}
  
  !> \{ \brief Domaingrösse (Periodizität-bereinigt)
  integer                ::  dim1, dim2, dim3
  !> \}
  
  !> \{ \brief Druck / Konzentrationen (inklusive Rand)
  integer                ::  S1p, S2p, S3p
  integer                ::  N1p, N2p, N3p
  !> \}
  
  !> \{ \brief Geschwindigkeiten (inklusive Rand)
  !!
  !! start values for the velocity including the boundary value
  !! first int indicates the spatial direction, the second in denotes the velocity direction
  integer                ::  S11B, S21B, S31B
  integer                ::  S12B, S22B, S32B
  integer                ::  S13B, S23B, S33B
  !> \}
  
  !> \{ \brief end values for the velocity including the boundary value
  integer                ::  N11B, N21B, N31B
  integer                ::  N12B, N22B, N32B
  integer                ::  N13B, N23B, N33B
  !> \}
  
  !> \{ \brief Geschwindigkeiten (exklusive Rand)
  !! start and end indices for the velocity without the boundary nodes
  integer                ::  S11, S21, S31
  integer                ::  S12, S22, S32
  integer                ::  S13, S23, S33
  
  integer                ::  N11, N21, N31
  integer                ::  N12, N22, N32
  integer                ::  N13, N23, N33
  !> \}
  
  !> \{ \brief Konzentrationen (exklusive Rand)
  integer                ::  S1c(1:n_conc), S2c(1:n_conc), S3c(1:n_conc)
  integer                ::  N1c(1:n_conc), N2c(1:n_conc), N3c(1:n_conc)
  !> \}
  
  !> \{ \brief grobe Gitter (Multigrid, INklusive Rand)
  integer                ::  S1R, S2R, S3R
  integer                ::  d1R, d2R, d3R
  !> \}
  
  !> \{ \brief grobe Gitter (Multigrid, EXklusive Rand)
  integer                ::  S11R, S22R, S33R
  integer                ::  d11R, d22R, d33R
  !> \}
  
  !> \{ \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)
  integer, parameter     ::  ls1 = -1
  integer, parameter     ::  ls2 = -1
  integer, parameter     ::  ls3 = -1
  !> \}
  
  !> \{ \brief Austauschrichtung (Multigrid)
  !!
  !! ex = -1: unten <--  oben
  !! ex =  0: unten <--> oben
  !! ex =  1: unten  --> oben
  integer                ::  ex1, ex2, ex3
  !\}
  
  
  !\}

  !===========================================================================================================

  !> \{ \name Boundary conditions
  !!    symmetrische, zentrale Stencils:
  !!     - Symmetrie-RB:   BC = -2
  !!     - periodische RB: BC = -1
  !!
  !!    schiefe, nicht-zentrale Stencils:
  !!     - Nachbar-Block:  BC =  0
  !!     - Dirichlet-RB:   BC =  1
  !!     - Neumann-RB:     BC =  2
  !!     - Robin-RB:       BC =  3

  !> \{ \name global boundary conditions
  logical                ::  outlet   (1:3,1:2,1:3)
  logical                ::  isopycnal(1:3,1:2,1:n_conc)
  
  integer                ::  BC_1L_global, BC_1U_global
  integer                ::  BC_2L_global, BC_2U_global
  integer                ::  BC_3L_global, BC_3U_global
  !> \}
  
  !> \{ \name lokal (Block) boundary conditions
  integer                ::  BC_1L, BC_1U
  integer                ::  BC_2L, BC_2U
  integer                ::  BC_3L, BC_3U
  
  integer                ::  BCc_1L(1:n_conc), BCc_1U(1:n_conc)
  integer                ::  BCc_2L(1:n_conc), BCc_2U(1:n_conc)
  integer                ::  BCc_3L(1:n_conc), BCc_3U(1:n_conc)
  !> \}
  
  !> \{ \name field properties
  integer                ::  n_gather(1:3,1:n_grids_max)
  integer                ::  NN (1:3,1:n_grids_max)
  integer                ::  NB (1:3,1:n_grids_max)
  integer                ::  iB (1:3,1:n_grids_max)
  integer                ::  SNF(1:2,1:3,1:n_grids_max)
  integer                ::  SNB(1:2,1:3,1:n_grids_max)
  integer                ::  BC (1:2,1:3,1:n_grids_max)
  integer                ::  ngb(1:2,1:3,1:n_grids_max)
  integer                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  integer                ::  rankc2(1:n_grids_max)
  logical                ::  participate_yes(1:n_grids_max)
  integer                ::  recvR(    1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  recvI(    1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  dispR(    1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  dispI(    1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  integer                ::  sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name physikalische Parameter
  !===========================================================================================================
  real                   ::  L1, L2, L3
  real                   ::  Re, Fro
  real                   ::  Sc (1:n_conc)
  real                   ::  Ric(1:n_conc)
  real                   ::  Rip(1:n_spec) ! TEST!!!
  real                   ::  usc(1:n_conc)
  real                   ::  usp(1:n_spec) ! TEST!!!
  real                   ::  Stp(1:n_spec) ! TEST!!!
  
  real                   ::  velReSc1(1:N2,1:N3,1:2,1:n_grids_max)
  real                   ::  velReSc2(1:N1,1:N3,1:2,1:n_grids_max)
  real                   ::  velReSc3(1:N1,1:N2,1:2,1:n_grids_max)
  
  real                   ::  usReSc (1:3,1:n_conc)
  real                   ::  us_vec (1:3,1:n_conc)
  real                   ::  gravity(1:3)
  
  real                   ::  dx_1L, dx_1U
  real                   ::  dx_2L, dx_2U
  real                   ::  dx_3L, dx_3U
  
  real                   ::  idx_1L, idx_1U
  real                   ::  idx_2L, idx_2U
  real                   ::  idx_3L, idx_3U
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name numerische Parameter
  !===========================================================================================================
  !> \{ \brief allgemein
  real                   ::  CFL
  real                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  integer                ::  timestep, timestep_old, substep, n_timesteps
  logical                ::  mapping_yes, upwind_yes, upwind_conc_yes
  logical                ::  Euler_yes, nonBoussinesq_yes, Stokes_yes, twostep_yes
  logical                ::  comp_visc_yes, comp_conv_yes, comp_inter_yes, comp_div_yes, comp_grad_yes
  logical, parameter     ::  filter_BC_yes = .true. ! TEST!!!
  !> timeint_mode = 0 (CN-RK3)
  !! timeint_mode = 1 (RK3-O3)
  integer                ::  timeint_mode, LES_mode, forcing_mode
  integer                ::  bulkflow_dir
  integer                ::  n_lp_vel           , n_hp_vel
  integer                ::  n_lp_conc(1:n_conc), n_hp_conc(1:n_conc)
  real                   ::  chi_vel
  real                   ::  chi_conc (1:n_conc)
  real                   ::  vel_bulk ! TEST!!!
  !> \}
  
  
  !> \{ \brief Runge-Kutta-Koeffizienten
  real   , parameter     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  real   , parameter     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  integer, parameter     ::  RK_steps = 3
  !> \}
  
  !> \{ \brief look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi)
  real   , parameter     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
                                      &               2.290031261, 2.361554127, 2.418567407, 2.462989697,  &
                                      &               2.496169963, 2.519146008, 2.532795254, 2.537935070,  &
                                      &               2.535397854, 2.526091466, 2.511046932, 2.491448818,  &
                                      &               2.468639045, 2.444084180, 2.419302172, 2.395757241,  &
                                      &               2.374745783, 2.357302135, 2.344145473, 2.335672458,  &
                                      &               2.331985072, 2.332936948, 2.338183901, 2.347230689,  &
                                      &               2.359471631, 2.374225928, 2.390769340, 2.408363261,  &
                                      &               2.426281290, 2.443832601, 2.460381269, 2.475360992,  &
                                      &               2.488285197, 2.498753090, 2.506452564, 2.511161051,  &
                                      &               2.512745327 /)
  integer, parameter     ::  n_stab = SIZE(stabilitylimit)
  !> \}
  
  !> \brief Helmholtz-Vorfaktor
  !!
  !! \f[0 \le \theta_L \le 1 \f]
  !! 0 means fully explicit and 1 means fully implicit
  real                   ::  thetaL
  !> \brief Helmholtz-Vorfaktor
  !!
  !! \f[ \mathrm{multL = -(1-thetaL)(aRK(substep)+bRK(substep))\frac{\Delta t}{Re} }\f]
  !! computed before mod_diff::Helmholtz call
  real                   ::  multL
  
  !> \{ \brief zeitliche Steuerung
  integer                ::  Int_dtime, Int_lev_pre
  
  integer                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  logical                ::  write_large, write_med, write_small
  real                   ::  time_out_scal, dtime_out_scal
  real                   ::  time_out_vect, dtime_out_vect
  
  logical                ::  write_out_scal, write_out_vect  
  logical                ::  new_dtime, finish_yes
    
  integer                ::  write_count
  integer                ::  restart
  character(len=3)       ::  restart_char
  
  integer                ::  n_conc_old
  !> \}

  !===========================================================================================================
  !> \{ \brief weitere Steuerungsoptionen
  !===========================================================================================================
  integer                ::  task
  logical                ::  read_nullspace_yes
  logical                ::  concentration_yes
  logical                ::  particles_yes
  logical                ::  nullspace_yes, nullspace_coarse_yes
  logical                ::  nullspace_ortho_yes
  
  logical                ::  write_stout_yes
  logical                ::  log_iteration_yes
  logical                ::  write_restart_yes
  logical                ::  write_lambda2_yes
  logical                ::  write_test_yes
  !> \}
  
  !> \{ \brief globale Laufindizes
  integer                ::  direction
  integer                ::  conc_nu
  !> \}
  
  !> \brief explizite Behandlung von Ecken bei Dirichlet-Randbedingungen
  !!
  !! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  logical, parameter     ::  corner_yes = .true.
  
  !> \{ \brief Systemzeit
  integer                ::  elatime, ctime(1:8)
  integer                ::  day, hour, minu, sec, msec
  !> \}
  
  !> \brief Partikel-Energie
  real                   ::  Ekin_part
  
  
  !===========================================================================================================
  !> \{ \name Iterationsparameter
  !===========================================================================================================
  !> \{ \brief Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten
  real                   ::  epsU, epsU0
  !> \}
  
  !> \brief Glaetter
  logical                ::  Jacobi_yes
  
  !> \{ \brief max. Anzahl Iterationen
  integer                ::  n_it_outer
  integer                ::  n_it_Poisson
  integer                ::  n_it_Helmh_vel
  integer                ::  n_it_Helmh_conc
  !> \}
  
  !> \{ \brief erwartete Konvergenzrate (äussere Iteration)
  real                   ::  precRatio0 (1:RK_steps)
  real                   ::  precOffset0(1:RK_steps)
  real, allocatable      ::  precRatio  (:,:)
  real, allocatable      ::  precOffset (:,:)
  !> \}
  
  !> \{ \brief Null-Initialisierung (äussere Iteration)
  logical                ::  init_pre(1:RK_steps), init_vel(1:RK_steps), init_conc(1:RK_steps)
  !> \}
  
  !> \{ \brief Vorkonditionierung (Multigrid)
  integer                ::  precond_outer
  integer                ::  precond_Poisson
  integer                ::  precond_Helmh_vel
  integer                ::  precond_Helmh_conc
  !> \}
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  integer                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !> \brief implizite Richtungen bei Linienrelaxation (Multigrid)
  integer                ::  impl_dir(1:3)
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  logical                ::  weighting_yes
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsstatistik
  !===========================================================================================================
  real                   ::  dtime_average
  real                   ::  max_div_init(1:2)
  integer                ::  number_poisson
  
  !> \{ \brief Zähler
  integer                ::  countO(1:RK_steps)
  integer                ::  countP(1:RK_steps,1:2)
  integer                ::  countH(1:RK_steps,1:3)
  !> \}
  
  !> \{ \brief Konvergenzrate
  real                   ::  ratioO(1:RK_steps)
  real                   ::  ratioH(1:RK_steps,1:3)
  real                   ::  ratioP(1:RK_steps,1:2)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name MPI
  !===========================================================================================================
  !> \{ \brief Kommunikatoren
  integer                ::  COMM_CART
  
  integer                ::  COMM_SLICE1, COMM_BAR1
  integer                ::  COMM_SLICE2, COMM_BAR2
  integer                ::  COMM_SLICE3, COMM_BAR3
  !> \}
  
  !> \{ \brief Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes)
  !!
  !! (for MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  integer                ::  bar1_size(1:NB1), bar1_offset(1:NB1)
  integer                ::  bar2_size(1:NB2), bar2_offset(1:NB2)
  integer                ::  bar3_size(1:NB3), bar3_offset(1:NB3)
  !> \{
  
  !> \{ \brief Ränge der Prozesse
  integer                ::  rank
  integer                ::  rank_bar1, rank_slice1
  integer                ::  rank_bar2, rank_slice2
  integer                ::  rank_bar3, rank_slice3
  !> \}
  
  !> \{ \brief Ränge der Nachbarprozesse (in kartesischem Gitter)
  integer                ::  rank1L, rank1U
  integer                ::  rank2L, rank2U
  integer                ::  rank3L, rank3U
  !> \}
  
  !> \brief Error-Handle
  integer                ::  merror
  
  !> \{ \brief Request-Handles
  !!
  !! (müssen offenbar global angelegt werden)
  integer                ::  req1L, req1U
  integer                ::  req2L, req2U
  integer                ::  req3L, req3U
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name HDF5
  !===========================================================================================================
  integer                ::  herror
  !> \}
  
  !> \{ \name periodicity criterions
  real   , allocatable   ::  velp(:,:,:,:)       !< \brief stored velocity, used as reference to previous period
  real                   ::  tp
  real                   ::  freq               !< \brief frequency
  real                   ::  periodic_tol       !< \brief tolerance which, we use as stopping criterion. \f[ ||u(t)-u(t-T)|| < tol \f]
  !> \}
#endif

contains
  function fwrite_test_yes()bind(c,name='fwrite_test_yes') result(write_test_yes_out)
  implicit none
  logical(c_bool)                ::  write_test_yes_out
  write_test_yes_out = write_test_yes

  end function fwrite_test_yes
  

  function get_task()bind(c,name='fget_task') result(taskout)
  implicit none
  integer(c_int)                ::  taskout
  taskout = task

  end function get_task


  
end module mod_vars
