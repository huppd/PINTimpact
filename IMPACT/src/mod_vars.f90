!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> @brief contains many variables, like wokring stuff and stencils
MODULE mod_vars
  
  USE mod_dims
  
  IMPLICIT NONE
  
  
#ifdef ALLOC
  
  !===========================================================================================================
  !> \brief Raemliche Dimensionen
  !!
  !! 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
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
  INTEGER                ::  N1
  INTEGER                ::  N2
  INTEGER                ::  N3
  !> \}
  
  !> \{
  !! \brief --- Anzahl grobe Gitter (Multigrid) -----------------------------------------------------------------------
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
  !>\}
  
  !> \{
  !! \brief --- Dimensionen -------------------------------------------------------------------------------------------
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
  
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)
  !> \}
  
  !> \{
  !! \brief Anzahl Stencil-Koeffizienten (Feld)
  !!
  !! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
  
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
  
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)
  
  !> \} \}
  !===========================================================================================================
  !> \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !--- zentral -----------------------------------------------------------------------------------------------
  !> \{ \brief zentral
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U
  !> \}
  
  !> \{ \brief upwind (nicht-linear)
  !!
  !! (aktuell wird nicht zwischen auf- und abwärtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !!  zu können, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U
  !> \}
  
  !> \{ \brief Divergenz
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1
  !> \}
  
  !> \{ \brief Gradient
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief 1. Ableitung (zentral)
  REAL   , ALLOCATABLE   ::  cp1  (:,:)
  REAL   , ALLOCATABLE   ::  cp2  (:,:)
  REAL   , ALLOCATABLE   ::  cp3  (:,:)
  
  REAL   , ALLOCATABLE   ::  cu1  (:,:)
  REAL   , ALLOCATABLE   ::  cv2  (:,:)
  REAL   , ALLOCATABLE   ::  cw3  (:,:)
  
  ! Eigene Stencils notwendig, da
  !   - Symmetrie des Geschwindigkeitsfeldes nicht unbedingt Symmetrie bei Konzentrationen bedeutet,
  !   - Randbedinungen darauf gespeichert werden koennen.
  REAL   , ALLOCATABLE   ::  cc1  (:,:,:)
  REAL   , ALLOCATABLE   ::  cc2  (:,:,:)
  REAL   , ALLOCATABLE   ::  cc3  (:,:,:)
  
  !> \}
  !> \{ \brief 1. Ableitung (upwind)
  REAL   , ALLOCATABLE   ::  cNp1D(:,:)
  REAL   , ALLOCATABLE   ::  cNp2D(:,:)
  REAL   , ALLOCATABLE   ::  cNp3D(:,:)
  
  REAL   , ALLOCATABLE   ::  cNp1U(:,:)
  REAL   , ALLOCATABLE   ::  cNp2U(:,:)
  REAL   , ALLOCATABLE   ::  cNp3U(:,:)
  
  REAL   , ALLOCATABLE   ::  cNu1D(:,:)
  REAL   , ALLOCATABLE   ::  cNv2D(:,:)
  REAL   , ALLOCATABLE   ::  cNw3D(:,:)
  
  REAL   , ALLOCATABLE   ::  cNu1U(:,:)
  REAL   , ALLOCATABLE   ::  cNv2U(:,:)
  REAL   , ALLOCATABLE   ::  cNw3U(:,:)
  
  REAL   , ALLOCATABLE   ::  cNc1D(:,:,:)
  REAL   , ALLOCATABLE   ::  cNc2D(:,:,:)
  REAL   , ALLOCATABLE   ::  cNc3D(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cNc1U(:,:,:)
  REAL   , ALLOCATABLE   ::  cNc2U(:,:,:)
  REAL   , ALLOCATABLE   ::  cNc3U(:,:,:)
  !> \}
  
  !> \{ \brief Divergenz
  REAL   , ALLOCATABLE   ::  cDu1 (:,:)
  REAL   , ALLOCATABLE   ::  cDv2 (:,:)
  REAL   , ALLOCATABLE   ::  cDw3 (:,:)
  !> \}
  
  !> \{ \brief Divergenz (transponiert)
  REAL   , ALLOCATABLE   ::  cDu1T(:,:)
  REAL   , ALLOCATABLE   ::  cDv2T(:,:)
  REAL   , ALLOCATABLE   ::  cDw3T(:,:)
  !> \}
  
  !> \{ \brief Gradient
  REAL   , ALLOCATABLE   ::  cGp1 (:,:)
  REAL   , ALLOCATABLE   ::  cGp2 (:,:)
  REAL   , ALLOCATABLE   ::  cGp3 (:,:)
  !> \}
  
  !> \{ \brief Gradient (transponiert)
  REAL   , ALLOCATABLE   ::  cGp1T(:,:)
  REAL   , ALLOCATABLE   ::  cGp2T(:,:)
  REAL   , ALLOCATABLE   ::  cGp3T(:,:)
  !> \}
  
  !> \{ \brief 2. Ableitung (zentral)
  REAL   , ALLOCATABLE   ::  cp11 (:,:)
  REAL   , ALLOCATABLE   ::  cp22 (:,:)
  REAL   , ALLOCATABLE   ::  cp33 (:,:)
  
  REAL   , ALLOCATABLE   ::  cu11 (:,:)
  REAL   , ALLOCATABLE   ::  cv22 (:,:)
  REAL   , ALLOCATABLE   ::  cw33 (:,:)
  
  REAL   , ALLOCATABLE   ::  cc11 (:,:,:)
  REAL   , ALLOCATABLE   ::  cc22 (:,:,:)
  REAL   , ALLOCATABLE   ::  cc33 (:,:,:)
  !> \}
  
  !> \{ \brief Interpolation
  REAL   , ALLOCATABLE   ::  cIpu(:,:)
  REAL   , ALLOCATABLE   ::  cIpv(:,:)
  REAL   , ALLOCATABLE   ::  cIpw(:,:)
  
  REAL   , ALLOCATABLE   ::  cIup(:,:)
  REAL   , ALLOCATABLE   ::  cIvp(:,:)
  REAL   , ALLOCATABLE   ::  cIwp(:,:)
  
  REAL   , ALLOCATABLE   ::  cIcu(:,:,:)
  REAL   , ALLOCATABLE   ::  cIcv(:,:,:)
  REAL   , ALLOCATABLE   ::  cIcw(:,:,:)
  !> \}
  
  !> \{ \brief Filter
  REAL   , ALLOCATABLE   ::  cFp1(:,:)
  REAL   , ALLOCATABLE   ::  cFp2(:,:)
  REAL   , ALLOCATABLE   ::  cFp3(:,:)
  
  REAL   , ALLOCATABLE   ::  cFu1(:,:)
  REAL   , ALLOCATABLE   ::  cFv2(:,:)
  REAL   , ALLOCATABLE   ::  cFw3(:,:)
  
  REAL   , ALLOCATABLE   ::  cFc1(:,:,:)
  REAL   , ALLOCATABLE   ::  cFc2(:,:,:)
  REAL   , ALLOCATABLE   ::  cFc3(:,:,:)
  !> \}
  
  !> \{ \brief Integrator (nur f�r Druckgitter)
  REAL   , ALLOCATABLE   ::  cInt1(:,:)
  REAL   , ALLOCATABLE   ::  cInt2(:,:)
  REAL   , ALLOCATABLE   ::  cInt3(:,:)
  !> \}
  
  
  !> \{ \brief Compact
  INTEGER                ::  dimS1
  INTEGER                ::  dimS2
  INTEGER                ::  dimS3
  
  REAL   , ALLOCATABLE   ::  buffer1A(:,:,:), buffer1B(:,:,:)
  REAL   , ALLOCATABLE   ::  buffer2A(:,:,:), buffer2B(:,:,:)
  REAL   , ALLOCATABLE   ::  buffer3A(:,:,:), buffer3B(:,:,:)
  !----------------------------------------------------
  INTEGER, ALLOCATABLE   ::  disp1(:), recv1(:)
  INTEGER, ALLOCATABLE   ::  disp2(:), recv2(:)
  INTEGER, ALLOCATABLE   ::  disp3(:), recv3(:)
  !----------------------------------------------------
  !> \}
  
  !> \{ \brief Schur complement (square) matrices ---
  REAL   , ALLOCATABLE   ::  Sp1(:,:), Su1(:,:), Sc1(:,:,:), Sp11(:,:), Su11(:,:), Sc11(:,:,:)
  REAL   , ALLOCATABLE   ::  Sp2(:,:), Sv2(:,:), Sc2(:,:,:), Sp22(:,:), Sv22(:,:), Sc22(:,:,:)
  REAL   , ALLOCATABLE   ::  Sp3(:,:), Sw3(:,:), Sc3(:,:,:), Sp33(:,:), Sw33(:,:), Sc33(:,:,:)
  
  REAL   , ALLOCATABLE   ::  SDu1(:,:), SGp1(:,:), SIpu(:,:), SIup(:,:), SIcu(:,:,:)
  REAL   , ALLOCATABLE   ::  SDv2(:,:), SGp2(:,:), SIpv(:,:), SIvp(:,:), SIcv(:,:,:)
  REAL   , ALLOCATABLE   ::  SDw3(:,:), SGp3(:,:), SIpw(:,:), SIwp(:,:), SIcw(:,:,:)
  
  REAL   , ALLOCATABLE   ::  SDu1T(:,:), SGp1T(:,:)
  REAL   , ALLOCATABLE   ::  SDv2T(:,:), SGp2T(:,:)
  REAL   , ALLOCATABLE   ::  SDw3T(:,:), SGp3T(:,:)
  !> \}
  
  !> \{ \brief bottom (horizontally oriented) matrices
  REAL   , ALLOCATABLE   ::  Wp1(:,:), Wu1(:,:), Wc1(:,:,:), Wp11(:,:), Wu11(:,:), Wc11(:,:,:)
  REAL   , ALLOCATABLE   ::  Wp2(:,:), Wv2(:,:), Wc2(:,:,:), Wp22(:,:), Wv22(:,:), Wc22(:,:,:)
  REAL   , ALLOCATABLE   ::  Wp3(:,:), Ww3(:,:), Wc3(:,:,:), Wp33(:,:), Ww33(:,:), Wc33(:,:,:)
  
  REAL   , ALLOCATABLE   ::  WDu1(:,:), WGp1(:,:), WIpu(:,:), WIup(:,:), WIcu(:,:,:)
  REAL   , ALLOCATABLE   ::  WDv2(:,:), WGp2(:,:), WIpv(:,:), WIvp(:,:), WIcv(:,:,:)
  REAL   , ALLOCATABLE   ::  WDw3(:,:), WGp3(:,:), WIpw(:,:), WIwp(:,:), WIcw(:,:,:)
  
  REAL   , ALLOCATABLE   ::  WDu1T(:,:), WGp1T(:,:)
  REAL   , ALLOCATABLE   ::  WDv2T(:,:), WGp2T(:,:)
  REAL   , ALLOCATABLE   ::  WDw3T(:,:), WGp3T(:,:)
  !> \}
  
  !> \{ \brief implicit finite difference coefficients ("L"=left-hand side)
  REAL   , ALLOCATABLE   ::  cp1CL(:,:), cu1CL(:,:), cc1CL(:,:,:), cp11CL(:,:), cu11CL(:,:), cc11CL(:,:,:)
  REAL   , ALLOCATABLE   ::  cp2CL(:,:), cv2CL(:,:), cc2CL(:,:,:), cp22CL(:,:), cv22CL(:,:), cc22CL(:,:,:)
  REAL   , ALLOCATABLE   ::  cp3CL(:,:), cw3CL(:,:), cc3CL(:,:,:), cp33CL(:,:), cw33CL(:,:), cc33CL(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CL(:,:), cGp1CL(:,:), cIpuCL(:,:), cIupCL(:,:), cIcuCL(:,:,:)
  REAL   , ALLOCATABLE   ::  cDv2CL(:,:), cGp2CL(:,:), cIpvCL(:,:), cIvpCL(:,:), cIcvCL(:,:,:)
  REAL   , ALLOCATABLE   ::  cDw3CL(:,:), cGp3CL(:,:), cIpwCL(:,:), cIwpCL(:,:), cIcwCL(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CLT(:,:), cGp1CLT(:,:)
  REAL   , ALLOCATABLE   ::  cDv2CLT(:,:), cGp2CLT(:,:)
  REAL   , ALLOCATABLE   ::  cDw3CLT(:,:), cGp3CLT(:,:)
  
  REAL   , ALLOCATABLE   ::  cp1CL_LU(:,:), cu1CL_LU(:,:), cc1CL_LU(:,:,:), cp11CL_LU(:,:), cu11CL_LU(:,:), cc11CL_LU(:,:,:)
  REAL   , ALLOCATABLE   ::  cp2CL_LU(:,:), cv2CL_LU(:,:), cc2CL_LU(:,:,:), cp22CL_LU(:,:), cv22CL_LU(:,:), cc22CL_LU(:,:,:)
  REAL   , ALLOCATABLE   ::  cp3CL_LU(:,:), cw3CL_LU(:,:), cc3CL_LU(:,:,:), cp33CL_LU(:,:), cw33CL_LU(:,:), cc33CL_LU(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CL_LU(:,:), cGp1CL_LU(:,:), cIpuCL_LU(:,:), cIupCL_LU(:,:), cIcuCL_LU(:,:,:)
  REAL   , ALLOCATABLE   ::  cDv2CL_LU(:,:), cGp2CL_LU(:,:), cIpvCL_LU(:,:), cIvpCL_LU(:,:), cIcvCL_LU(:,:,:)
  REAL   , ALLOCATABLE   ::  cDw3CL_LU(:,:), cGp3CL_LU(:,:), cIpwCL_LU(:,:), cIwpCL_LU(:,:), cIcwCL_LU(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CLT_LU(:,:), cGp1CLT_LU(:,:)
  REAL   , ALLOCATABLE   ::  cDv2CLT_LU(:,:), cGp2CLT_LU(:,:)
  REAL   , ALLOCATABLE   ::  cDw3CLT_LU(:,:), cGp3CLT_LU(:,:)
  !> \}
  
  !> \{ \brief explicit finite difference coefficients ("R"=right-hand side)
  REAL   , ALLOCATABLE   ::  cp1CR(:,:), cu1CR(:,:), cc1CR(:,:,:), cp11CR(:,:), cu11CR(:,:), cc11CR(:,:,:)
  REAL   , ALLOCATABLE   ::  cp2CR(:,:), cv2CR(:,:), cc2CR(:,:,:), cp22CR(:,:), cv22CR(:,:), cc22CR(:,:,:)
  REAL   , ALLOCATABLE   ::  cp3CR(:,:), cw3CR(:,:), cc3CR(:,:,:), cp33CR(:,:), cw33CR(:,:), cc33CR(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CR(:,:), cGp1CR(:,:), cIpuCR(:,:), cIupCR(:,:), cIcuCR(:,:,:)
  REAL   , ALLOCATABLE   ::  cDv2CR(:,:), cGp2CR(:,:), cIpvCR(:,:), cIvpCR(:,:), cIcvCR(:,:,:)
  REAL   , ALLOCATABLE   ::  cDw3CR(:,:), cGp3CR(:,:), cIpwCR(:,:), cIwpCR(:,:), cIcwCR(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cDu1CRT(:,:), cGp1CRT(:,:)
  REAL   , ALLOCATABLE   ::  cDv2CRT(:,:), cGp2CRT(:,:)
  REAL   , ALLOCATABLE   ::  cDw3CRT(:,:), cGp3CRT(:,:)
  !> \}
  !****************************************************
  
  
  !--- 2. Ableitung (Multigrid) ------------------------------------------------------------------------------ 
  ! Anmerkung: Die Koeffizientensätze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL   , ALLOCATABLE   ::  cp11R(:,:,:)
  REAL   , ALLOCATABLE   ::  cp22R(:,:,:)
  REAL   , ALLOCATABLE   ::  cp33R(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cu11R(:,:,:)
  REAL   , ALLOCATABLE   ::  cv22R(:,:,:)
  REAL   , ALLOCATABLE   ::  cw33R(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cc11R(:,:,:,:)
  REAL   , ALLOCATABLE   ::  cc22R(:,:,:,:)
  REAL   , ALLOCATABLE   ::  cc33R(:,:,:,:)
  
  REAL   , ALLOCATABLE   ::  cdg1 (:,:,:)
  REAL   , ALLOCATABLE   ::  cdg2 (:,:,:)
  REAL   , ALLOCATABLE   ::  cdg3 (:,:,:)
  
  !--- Interpolation (Multigrid) ----------------------------------------------------------------------------- 
  REAL   , ALLOCATABLE   ::  cI1(:,:,:)
  REAL   , ALLOCATABLE   ::  cI2(:,:,:)
  REAL   , ALLOCATABLE   ::  cI3(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cIH1(:,:)
  REAL   , ALLOCATABLE   ::  cIH2(:,:)
  REAL   , ALLOCATABLE   ::  cIH3(:,:)
  
  !--- Restriktion (Multigrid) ------------------------------------------------------------------------------- 
  REAL   , ALLOCATABLE   ::  cR1 (:,:,:)
  REAL   , ALLOCATABLE   ::  cR2 (:,:,:)
  REAL   , ALLOCATABLE   ::  cR3 (:,:,:)
  
  REAL   , ALLOCATABLE   ::  cRest1(:,:,:) ! TEST!!!
  REAL   , ALLOCATABLE   ::  cRest2(:,:,:)
  REAL   , ALLOCATABLE   ::  cRest3(:,:,:)
  
  REAL   , ALLOCATABLE   ::  cRH1(:,:)
  REAL   , ALLOCATABLE   ::  cRH2(:,:)
  REAL   , ALLOCATABLE   ::  cRH3(:,:)
  
  
  !===========================================================================================================
  !> \{ \name Gitterspezifikationen
  !===========================================================================================================
  !> \{
  !! \brief physiklische Koordinaten (global)
  REAL   , ALLOCATABLE   ::  y1p(:), y1u(:)
  REAL   , ALLOCATABLE   ::  y2p(:), y2v(:)
  REAL   , ALLOCATABLE   ::  y3p(:), y3w(:)
  !> \}
  
  !> \{
  !! \brief physiklische Koordinaten (Block)
  REAL   , ALLOCATABLE   ::  x1p(:), x1u(:)
  REAL   , ALLOCATABLE   ::  x2p(:), x2v(:)
  REAL   , ALLOCATABLE   ::  x3p(:), x3w(:)
  !> \}
  
  !> \{ \brief physiklische Koordinaten (Block, Multigrid)
  REAL   , ALLOCATABLE   ::  x1pR(:,:), x1uR(:,:)
  REAL   , ALLOCATABLE   ::  x2pR(:,:), x2vR(:,:)
  REAL   , ALLOCATABLE   ::  x3pR(:,:), x3wR(:,:)
  !> \}
  
  !> \{ \brief Gitterweiten (global)
  REAL   , ALLOCATABLE   ::  dy1p(:), dy1u(:)
  REAL   , ALLOCATABLE   ::  dy2p(:), dy2v(:)
  REAL   , ALLOCATABLE   ::  dy3p(:), dy3w(:)
  !> \}
  
  !> \{ \brief Gitterweiten (Block)
  REAL   , ALLOCATABLE   ::  dx1p(:), dx1u(:)
  REAL   , ALLOCATABLE   ::  dx2p(:), dx2v(:)
  REAL   , ALLOCATABLE   ::  dx3p(:), dx3w(:)
  !> \}
  
  REAL   , ALLOCATABLE   ::  dx1DM(:), dx1pM(:), ddx1pM(:)
  REAL   , ALLOCATABLE   ::  dx2DM(:), dx2pM(:), ddx2pM(:)
  REAL   , ALLOCATABLE   ::  dx3DM(:), dx3pM(:), ddx3pM(:)
  
  REAL   , ALLOCATABLE   ::  dx1GM(:), dx1uM(:), ddx1uM(:)
  REAL   , ALLOCATABLE   ::  dx2GM(:), dx2vM(:), ddx2vM(:)
  REAL   , ALLOCATABLE   ::  dx3GM(:), dx3wM(:), ddx3wM(:)
  
  !> \{ \brief Smagorinsky-Modell
  REAL   , ALLOCATABLE   ::  dx1pS(:), dx1uS(:) ! TEST!!!
  REAL   , ALLOCATABLE   ::  dx2pS(:), dx2vS(:)
  REAL   , ALLOCATABLE   ::  dx3pS(:), dx3wS(:)
  !> \}
  
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Arbeitsfelder
  !===========================================================================================================
  !> \brief Geschwindigkeiten
  REAL   , ALLOCATABLE   ::  vel(:,:,:,:)
  
  !> \brief nicht-linearer Term
  REAL   , ALLOCATABLE   ::  nl (:,:,:,:)
  
  !> \brief Recht-Hand-Seite
  REAL   , ALLOCATABLE   ::  rhs(:,:,:,:)
  
  !> \brief Druck
  REAL   , ALLOCATABLE   ::  pre(:,:,:)
  
  !> \brief Konzentrationsfelder
  REAL   , ALLOCATABLE   ::  conc(:,:,:,:)
  
  !> \brief nicht-linearer Term (Konzentration)
  REAL   , ALLOCATABLE   ::  nlco(:,:,:,:)
  
  !> \{ \brief Ausfluss-RB (Geschwindigkeitsfeld)
  !!
  !! Da die RHS f�r die Konzentrationsfelder nicht �ber die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, m�ssen mindestens die zugeh�rigen Randbedingungen gespeichert werden.
  REAL   , ALLOCATABLE   ::  bc11(:,:,:), nlbc11(:,:,:)
  REAL   , ALLOCATABLE   ::  bc12(:,:,:), nlbc12(:,:,:)
  REAL   , ALLOCATABLE   ::  bc13(:,:,:), nlbc13(:,:,:)
  
  REAL   , ALLOCATABLE   ::  bc21(:,:,:), nlbc21(:,:,:)
  REAL   , ALLOCATABLE   ::  bc22(:,:,:), nlbc22(:,:,:)
  REAL   , ALLOCATABLE   ::  bc23(:,:,:), nlbc23(:,:,:)
  
  REAL   , ALLOCATABLE   ::  bc31(:,:,:), nlbc31(:,:,:)
  REAL   , ALLOCATABLE   ::  bc32(:,:,:), nlbc32(:,:,:)
  REAL   , ALLOCATABLE   ::  bc33(:,:,:), nlbc33(:,:,:)
  
  REAL   , ALLOCATABLE   ::  drift1(:,:,:)
  REAL   , ALLOCATABLE   ::  drift2(:,:,:)
  REAL   , ALLOCATABLE   ::  drift3(:,:,:)
  !> \}
  
  
  !> \{ \brief Sedimentations-RB (Konzentration)
  !! Da die RHS f�r die Konzentrationsfelder nicht �ber die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, m�ssen mindestens die zugeh�rigen Randbedingungen gespeichert werden.
  REAL   , ALLOCATABLE   ::  sed1L (:,:,:)
  REAL   , ALLOCATABLE   ::  sed1U (:,:,:)
  
  REAL   , ALLOCATABLE   ::  sed2L (:,:,:)
  REAL   , ALLOCATABLE   ::  sed2U (:,:,:)
  
  REAL   , ALLOCATABLE   ::  sed3L (:,:,:)
  REAL   , ALLOCATABLE   ::  sed3U (:,:,:)
  
  REAL   , ALLOCATABLE   ::  conc1L(:,:,:)
  REAL   , ALLOCATABLE   ::  conc1U(:,:,:)
  
  REAL   , ALLOCATABLE   ::  conc2L(:,:,:)
  REAL   , ALLOCATABLE   ::  conc2U(:,:,:)
  
  REAL   , ALLOCATABLE   ::  conc3L(:,:,:)
  REAL   , ALLOCATABLE   ::  conc3U(:,:,:)
  !> \}
  
  !> \{ \brief Sedimentation (Konzentration, Ausschrieb)
  REAL   , ALLOCATABLE   ::  dep1L_conc(:,:,:)
  REAL   , ALLOCATABLE   ::  dep1U_conc(:,:,:)
  
  REAL   , ALLOCATABLE   ::  dep2L_conc(:,:,:)
  REAL   , ALLOCATABLE   ::  dep2U_conc(:,:,:)
  
  REAL   , ALLOCATABLE   ::  dep3L_conc(:,:,:)
  REAL   , ALLOCATABLE   ::  dep3U_conc(:,:,:)
  !> \}
  
  !> \< \brief Sedimentation (Partikel, Ausschrieb)
  ! TEST!!! 1:n_spec???
  REAL   , ALLOCATABLE   ::  dep1L_part(:,:)
  REAL   , ALLOCATABLE   ::  dep1U_part(:,:)
  
  REAL   , ALLOCATABLE   ::  dep2L_part(:,:)
  REAL   , ALLOCATABLE   ::  dep2U_part(:,:)
  
  REAL   , ALLOCATABLE   ::  dep3L_part(:,:)
  REAL   , ALLOCATABLE   ::  dep3U_part(:,:)
  !> \}
  
  !> \brief Residuum
  REAL   , ALLOCATABLE   ::  res (:,:,:)
  
  !> \brief Druckgradient (eine Komponente)
  REAL   , ALLOCATABLE   ::  gpre(:,:,:)
  
  !> \brief Gewichte f�r Divergenzfreiheit
  REAL   , ALLOCATABLE   ::  weight(:,:,:)
  
  !> \brief fluid density on velocity grid
#ifdef NONBOUSSINESQ
  REAL   , ALLOCATABLE   ::  dens(:,:,:,:)
#endif
  
  !> \{ \brief Null-Raum-Vektor
  REAL   , ALLOCATABLE   ::  psi    (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_vel(:,:,:,:)
  
  REAL   , ALLOCATABLE   ::  psi_rel1 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th11(:,:,:)
  REAL   , ALLOCATABLE   ::  th12(:,:,:)
  REAL   , ALLOCATABLE   ::  th13(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th21(:,:,:)
  REAL   , ALLOCATABLE   ::  th22(:,:,:)
  REAL   , ALLOCATABLE   ::  th23(:,:,:)
  
  REAL   , ALLOCATABLE   ::  th31(:,:,:)
  REAL   , ALLOCATABLE   ::  th32(:,:,:)
  REAL   , ALLOCATABLE   ::  th33(:,:,:)
  !> \}
  
  !> \{ \brief Multigrid
  REAL   , ALLOCATABLE   ::  vec1C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec15B(:,:,:)
  !> \}
  
  
  !> \{ \brief BiCGstab / Richardson
  REAL   , ALLOCATABLE   ::  pp(:,:,:)
  REAL   , ALLOCATABLE   ::  Ap(:,:,:)
  REAL   , ALLOCATABLE   ::  rr(:,:,:)
  REAL   , ALLOCATABLE   ::  rh(:,:,:)
  REAL   , ALLOCATABLE   ::  Ar(:,:,:)
  REAL   , ALLOCATABLE   ::  z1(:,:,:)
  REAL   , ALLOCATABLE   ::  z2(:,:,:)
  !> \}
  
  !> \brief product_div_grad
  REAL   , ALLOCATABLE   ::  dig(:,:,:)
  
  !> \brief Hilfsfeld (kompakte Differenzen)
  REAL   , ALLOCATABLE   ::  com(:,:,:)
  
  !> \{ \brief Hilfsfelder (Druckiteration)
  !!
  !! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  !! - wird auch fuer interpolierte Geschwindigkeiten in rhs_NS und rhs_conc verwendet.
  REAL   , ALLOCATABLE   ::  work1(:,:,:)
  REAL   , ALLOCATABLE   ::  work2(:,:,:)
  REAL   , ALLOCATABLE   ::  work3(:,:,:)
  !> \}
  
  !> \{ \brief Lagrange-Partikel
  INTEGER                ::  n_part, n_part_tot, n_groups
  REAL   , ALLOCATABLE   ::  particles(:,:)
  !> \}
  
  !> \brief Smagorinsky LES-Modell
  REAL   , ALLOCATABLE   ::  smag_diffus(:,:,:,:)
  
  !> \{ \brief Linienrelaxation
  REAL   , ALLOCATABLE   ::  vec1(:), dia1(:), SOR1(:), band1(:,:) ! TEST!!! siehe unten ...
  REAL   , ALLOCATABLE   ::  vec2(:), dia2(:), SOR2(:), band2(:,:)
  REAL   , ALLOCATABLE   ::  vec3(:), dia3(:), SOR3(:), band3(:,:)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Indizierung (Intervallgrenzen, Verschiebungen)
  !!
  !!===========================================================================================================
  !!--- Block-Index -------------------------------------------------------------------------------------------
  !!INTEGER               ::  iB, jB, kB ! TEST!!! iB(1:3,1:n_grids_max) hierher verschieben ...
  
  !> \{ \brief Indexverschiebung (Block --> global)
  INTEGER                ::  iShift, jShift, kShift
  !> \}
  
  !--- Domaingr�sse (Periodizit�t-bereinigt) -----------------------------------------------------------------
  INTEGER                ::  dim1, dim2, dim3
  
  !--- Druck / Konzentrationen (inklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  
  !--- Geschwindigkeiten (inklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  
  !--- Geschwindigkeiten (exklusive Rand) --------------------------------------------------------------------
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
  
  !--- Konzentrationen (exklusive Rand) ----------------------------------------------------------------------
  INTEGER, ALLOCATABLE   ::  S1c(:), S2c(:), S3c(:)
  INTEGER, ALLOCATABLE   ::  N1c(:), N2c(:), N3c(:)
  
  !--- grobe Gitter (Multigrid, INklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  
  !--- grobe Gitter (Multigrid, EXklusive Rand) --------------------------------------------------------------
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  
  !> \{ \brief Ueberlappungskonvention der Bl�cke (Multigrid, siehe mod_setup)
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1
  !> \}
  
  !> \{ \brief Austauschrichtung (Multigrid) -------------------------------------------------------------------------
  !!
  !! ex = -1: unten <--  oben
  !! ex =  0: unten <--> oben
  !! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Randbedingungen
  !!===========================================================================================================
  !!                              _
  !!    Symmetrie-RB:   BC = -2    |
  !!    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !!    Nachbar-Block:  BC =  0   _|
  !!    Dirichlet-RB:   BC =  1    |
  !!    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !!    Robin-RB:       BC =  3   _|

  !> \{ \brief global
  LOGICAL                ::  outlet(1:3,1:2,1:3)
  LOGICAL, ALLOCATABLE   ::  isopycnal(:,:,:)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  !> \}
  
  !> \{ \brief lokal (Block)
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U
  
  INTEGER, ALLOCATABLE   ::  BCc_1L(:), BCc_1U(:)
  INTEGER, ALLOCATABLE   ::  BCc_2L(:), BCc_2U(:)
  INTEGER, ALLOCATABLE   ::  BCc_3L(:), BCc_3U(:)
  !> \}
  
  !> \{ \brief field properties
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER                ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER, ALLOCATABLE   ::  recvR(  :,:), recvI(  :,:)
  INTEGER, ALLOCATABLE   ::  dispR(  :,:), dispI(  :,:)
  INTEGER, ALLOCATABLE   ::  offsR(:,:,:), offsI(:,:,:)
  INTEGER, ALLOCATABLE   ::  sizsR(:,:,:), sizsI(:,:,:)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ name physikalische Parameter
  !===========================================================================================================
  REAL                   ::  L1, L2, L3
  REAL                   ::  Re, Fro
  REAL   , ALLOCATABLE   ::  Sc (:)
  REAL   , ALLOCATABLE   ::  Ric(:)
  REAL   , ALLOCATABLE   ::  Rip(:) ! TEST!!!
  REAL   , ALLOCATABLE   ::  usc(:) ! TEST!!! ebenfalls umbenennen!
  REAL   , ALLOCATABLE   ::  usp(:) ! TEST!!!
  REAL   , ALLOCATABLE   ::  Stp(:) ! TEST!!!
  
  REAL   , ALLOCATABLE   ::  velReSc1(:,:,:,:)
  REAL   , ALLOCATABLE   ::  velReSc2(:,:,:,:)
  REAL   , ALLOCATABLE   ::  velReSc3(:,:,:,:)
  
  REAL   , ALLOCATABLE   ::  usReSc (:,:)
  REAL   , ALLOCATABLE   ::  us_vec (:,:)
  REAL                   ::  gravity(1:3)
  
  REAL                   ::  dx_1L, dx_1U
  REAL                   ::  dx_2L, dx_2U
  REAL                   ::  dx_3L, dx_3U
  
  REAL                   ::  idx_1L, idx_1U
  REAL                   ::  idx_2L, idx_2U
  REAL                   ::  idx_3L, idx_3U
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name numerische Parameter
  !===========================================================================================================
  !--- allgemein ---------------------------------------------------------------------------------------------
  REAL                   ::  CFL
  REAL                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes, upwind_conc_yes
  LOGICAL                ::  Euler_yes, nonBoussinesq_yes, Stokes_yes, twostep_yes
  LOGICAL                ::  comp_visc_yes, comp_conv_yes, comp_inter_yes, comp_div_yes, comp_grad_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  INTEGER                ::  timeint_mode, LES_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel    , n_hp_vel
  INTEGER, ALLOCATABLE   ::  n_lp_conc(:), n_hp_conc(:)
  REAL                   ::  chi_vel
  REAL   , ALLOCATABLE   ::  chi_conc (:)
  REAL                   ::  vel_bulk ! TEST!!!
  
  
  !--- Runge-Kutta-Koeffizienten -----------------------------------------------------------------------------
  REAL   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3
  
  !--- look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi) ---------------------------
  REAL   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
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
  INTEGER, PARAMETER     ::  n_stab = SIZE(stabilitylimit)
  
  !--- Helmholtz-Vorfaktoren ---------------------------------------------------------------------------------
  REAL                   ::  thetaL, multL
  
  !--- zeitliche Steuerung -----------------------------------------------------------------------------------
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
  REAL                   ::  time_out_scal, dtime_out_scal
  REAL                   ::  time_out_vect, dtime_out_vect
  
  LOGICAL                ::  write_out_scal, write_out_vect
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old

  
  !> \}
  !===========================================================================================================
  !> \{ \name weitere Steuerungsoptionen
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  concentration_yes
  LOGICAL                ::  particles_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  
  !> \{ \brief globale Laufindizes
  INTEGER                ::  direction
  INTEGER                ::  conc_nu
  !> \}
  
  !> \brief explizite Behandlung von Ecken bei Dirichlet-Randbedingungen
  !!
  !! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort k�nstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.
  
  !> \{ \brief Systemzeit
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec
  !> \}
  
  !> \brief Partikel-Energie
  REAL                   ::  Ekin_part
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsparameter
  !===========================================================================================================
  !> \{ Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten
  REAL                   ::  epsU, epsU0
  !> \}
  
  !> \brief Glaetter
  LOGICAL                ::  Jacobi_yes
  
  !> \{ \brief max. Anzahl Iterationen
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  INTEGER                ::  n_it_Helmh_conc
  !> \}
  
  !> \{ \brief erwartete Konvergenzrate (�ussere Iteration)
  REAL                   ::  precRatio0 (1:RK_steps)
  REAL                   ::  precOffset0(1:RK_steps)
  REAL   , ALLOCATABLE   ::  precRatio  (:,:)
  REAL   , ALLOCATABLE   ::  precOffset (:,:)
  !> \}
  
  !> \{ \brief Null-Initialisierung (äussere Iteration)
  LOGICAL                ::  init_pre(1:RK_steps), init_vel(1:RK_steps), init_conc(1:RK_steps)
  !> \}
  
  !> \{ \name Vorkonditionierung (Multigrid)
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  INTEGER                ::  precond_Helmh_conc
  !> \}
  
  !> \{ \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  !> \}
  
  !> \brief implizite Richtungen bei Linienrelaxation (Multigrid)
  INTEGER                ::  impl_dir(1:3)
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  LOGICAL                ::  weighting_yes
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsstatistik
  !===========================================================================================================
  REAL                   ::  dtime_average
  REAL                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !> \{ \brief Zähler
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  !> \}
  
  !> \{ \brief Konvergenzrate
  REAL                   ::  ratioO(1:RK_steps)
  REAL                   ::  ratioH(1:RK_steps,1:3)
  REAL                   ::  ratioP(1:RK_steps,1:2)
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name MPI
  !===========================================================================================================
  !> \{ \brief Kommunikatoren
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  !> \}
  
  !> \{ \brief Dimension und Position der Bl�cke innerhalb der Kommunikatoren (Gitter-Indizes)
  !!
  !! (for MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER, ALLOCATABLE   ::  bar1_size(:), bar1_offset(:)
  INTEGER, ALLOCATABLE   ::  bar2_size(:), bar2_offset(:)
  INTEGER, ALLOCATABLE   ::  bar3_size(:), bar3_offset(:)
  !> \}
  
  !> \{ \brief Ränge der Prozesse
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  !> \}
  
  !> \{ \brief Ränge der Nachbarprozesse (in kartesischem Gitter)
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  !> \}
  
  !> \brief Error-Handle
  INTEGER                ::  merror
  
  !> \{ \brief Request-Handles
  !!
  !! (müssen offenbar global angelegt werden)
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  !> \}
  
  !> \}
  !===========================================================================================================
  !=== HDF5 ==================================================================================================
  !===========================================================================================================
  INTEGER                ::  herror
  
  !> \{ \name periodicity criterions
  REAL   , ALLOCATABLE   ::  velp(:,:,:,:)       !< \brief stored velocity, used as reference to previous period
  REAL                   ::  tp
  REAL                   ::  freq               !< \brief frequency
  REAL                   ::  periodic_tol       !< \brief tolerance which, we use as stopping criterion. \f[ ||u(t)-u(t-T)|| < tol \f]
  !> \}
  
  
#else
  
  
  
  
  !===========================================================================================================
  !> \brief Raemliche Dimensionen
  !===========================================================================================================
  !> 2D wird ueber M3==2 eingeschaltet
  INTEGER                ::  dimens
  
  
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
  INTEGER, PARAMETER     ::  N1 = 1+(M1-1)/NB1
  INTEGER, PARAMETER     ::  N2 = 1+(M2-1)/NB2
  INTEGER, PARAMETER     ::  N3 = 1+(M3-1)/NB3
  
  !> \{ \brief Anzahl grobe Gitter (Multigrid)
  INTEGER, PARAMETER     ::  n_grids_max = 15
  INTEGER                ::  n_grids, n_grids_limit
  !> \}
  
  !> \{ \brief Dimensionen
  INTEGER, PARAMETER     ::  dim_ncb1c = SIZE(ncb1c)
  INTEGER, PARAMETER     ::  dim_ncb1g = SIZE(ncb1g)
  INTEGER, PARAMETER     ::  dim_ncb1d = SIZE(ncb1d)
  
  INTEGER, PARAMETER     ::  dim_ncb2c = SIZE(ncb2c)
  INTEGER, PARAMETER     ::  dim_ncb2g = SIZE(ncb2g)
  INTEGER, PARAMETER     ::  dim_ncb2d = SIZE(ncb2d)
                         
  INTEGER, PARAMETER     ::  dim_ncb3c = SIZE(ncb3c)
  INTEGER, PARAMETER     ::  dim_ncb3g = SIZE(ncb3g)
  INTEGER, PARAMETER     ::  dim_ncb3d = SIZE(ncb3d)
  !> \}
  
  !> \{ \brief Anzahl Stencil-Koeffizienten (Feld)
  !!
  !! Anzahl der Koeffizienten im Feld (zentrale Differenzen angenommen):
  INTEGER, PARAMETER     ::  nc1c = ncb1c(dim_ncb1c)
  INTEGER, PARAMETER     ::  nc1s = ncb1g(dim_ncb1g)
                         
  INTEGER, PARAMETER     ::  nc2c = ncb2c(dim_ncb2c)
  INTEGER, PARAMETER     ::  nc2s = ncb2g(dim_ncb2g)
                         
  INTEGER, PARAMETER     ::  nc3c = ncb3c(dim_ncb3c)
  INTEGER, PARAMETER     ::  nc3s = ncb3g(dim_ncb3g)
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name Intervallgrenzen der Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief zentral
  INTEGER, PARAMETER     ::  b1U = nc1s/2
  INTEGER, PARAMETER     ::  b2U = nc2s/2
  INTEGER, PARAMETER     ::  b3U = nc3s/2
  
  INTEGER, PARAMETER     ::  b1L = -b1U
  INTEGER, PARAMETER     ::  b2L = -b2U
  INTEGER, PARAMETER     ::  b3L = -b3U
  !> \}
  
  !> \{ \brief upwind (nicht-linear).
  !! (aktuell wird nicht zwischen auf- und abw�rtsgerichteten Stencils unterschieden, um auch am Rand arbeiten
  !!  zu k�nnen, wo KEINE upwind-Differenzierung verwendet wird)
  INTEGER, PARAMETER     ::  n1L = b1L
  INTEGER, PARAMETER     ::  n2L = b2L
  INTEGER, PARAMETER     ::  n3L = b3L
  
  INTEGER, PARAMETER     ::  n1U = b1U
  INTEGER, PARAMETER     ::  n2U = b2U
  INTEGER, PARAMETER     ::  n3U = b3U
  !> \}
  
  !> \{ \brief Divergenz
  INTEGER, PARAMETER     ::  d1L = b1L
  INTEGER, PARAMETER     ::  d2L = b2L
  INTEGER, PARAMETER     ::  d3L = b3L
  
  INTEGER, PARAMETER     ::  d1U = b1U-1
  INTEGER, PARAMETER     ::  d2U = b2U-1
  INTEGER, PARAMETER     ::  d3U = b3U-1
  !> \}
  
  !> \{ \brief Gradient
  INTEGER, PARAMETER     ::  g1L = b1L+1
  INTEGER, PARAMETER     ::  g2L = b2L+1
  INTEGER, PARAMETER     ::  g3L = b3L+1
  
  INTEGER, PARAMETER     ::  g1U = b1U
  INTEGER, PARAMETER     ::  g2U = b2U
  INTEGER, PARAMETER     ::  g3U = b3U
  !> \}
  
  
  !===========================================================================================================
  !> \{ \name Differenzen-Koeffizienten-Arrays
  !===========================================================================================================
  !> \{ \brief 1. Ableitung (zentral)
  REAL                   ::  cp1  (b1L:b1U,0:N1)
  REAL                   ::  cp2  (b2L:b2U,0:N2)
  REAL                   ::  cp3  (b3L:b3U,0:N3)
  
  REAL                   ::  cu1  (b1L:b1U,0:N1)
  REAL                   ::  cv2  (b2L:b2U,0:N2)
  REAL                   ::  cw3  (b3L:b3U,0:N3)
  !> \}
  
  !> \{ Eigene Stencils. notwendig, da
  !!   - Symmetrie des Geschwindigkeitsfeldes nicht unbedingt Symmetrie bei Konzentrationen bedeutet,
  !!   - Randbedinungen darauf gespeichert werden koennen.
  REAL                   ::  cc1  (b1L:b1U,0:N1,1:n_conc)
  REAL                   ::  cc2  (b2L:b2U,0:N2,1:n_conc)
  REAL                   ::  cc3  (b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief 1. Ableitung (upwind)
  REAL                   ::  cNp1D(n1L:n1U,0:N1)
  REAL                   ::  cNp2D(n2L:n2U,0:N2)
  REAL                   ::  cNp3D(n3L:n3U,0:N3)
  
  REAL                   ::  cNp1U(n1L:n1U,0:N1)
  REAL                   ::  cNp2U(n2L:n2U,0:N2)
  REAL                   ::  cNp3U(n3L:n3U,0:N3)
  
  REAL                   ::  cNu1D(n1L:n1U,0:N1)
  REAL                   ::  cNv2D(n2L:n2U,0:N2)
  REAL                   ::  cNw3D(n3L:n3U,0:N3)
  
  REAL                   ::  cNu1U(n1L:n1U,0:N1)
  REAL                   ::  cNv2U(n2L:n2U,0:N2)
  REAL                   ::  cNw3U(n3L:n3U,0:N3)
  
  REAL                   ::  cNc1D(n1L:n1U,0:N1,1:n_conc)
  REAL                   ::  cNc2D(n2L:n2U,0:N2,1:n_conc)
  REAL                   ::  cNc3D(n3L:n3U,0:N3,1:n_conc)
  
  REAL                   ::  cNc1U(n1L:n1U,0:N1,1:n_conc)
  REAL                   ::  cNc2U(n2L:n2U,0:N2,1:n_conc)
  REAL                   ::  cNc3U(n3L:n3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Divergenz
  REAL                   ::  cDu1 (d1L:d1U,0:N1)
  REAL                   ::  cDv2 (d2L:d2U,0:N2)
  REAL                   ::  cDw3 (d3L:d3U,0:N3)
  !> \}
  
  !> \{ \brief Divergenz (transponiert)
  REAL                   ::  cDu1T(g1L:g1U,0:N1)
  REAL                   ::  cDv2T(g2L:g2U,0:N2)
  REAL                   ::  cDw3T(g3L:g3U,0:N3)
  !> \}
  
  !> \{ \brief Gradient
  REAL                   ::  cGp1 (g1L:g1U,0:N1)
  REAL                   ::  cGp2 (g2L:g2U,0:N2)
  REAL                   ::  cGp3 (g3L:g3U,0:N3)
  !> \}
  
  !> \{ \brief Gradient (transponiert)
  REAL                   ::  cGp1T(d1L:d1U,0:N1)
  REAL                   ::  cGp2T(d2L:d2U,0:N2)
  REAL                   ::  cGp3T(d3L:d3U,0:N3)
  !> \}
  
  !> \{ \brief 2. Ableitung (zentral)
  REAL                   ::  cp11 (b1L:b1U,0:N1)
  REAL                   ::  cp22 (b2L:b2U,0:N2)
  REAL                   ::  cp33 (b3L:b3U,0:N3)
  
  REAL                   ::  cu11 (b1L:b1U,0:N1)
  REAL                   ::  cv22 (b2L:b2U,0:N2)
  REAL                   ::  cw33 (b3L:b3U,0:N3)
  
  REAL                   ::  cc11 (b1L:b1U,0:N1,1:n_conc)
  REAL                   ::  cc22 (b2L:b2U,0:N2,1:n_conc)
  REAL                   ::  cc33 (b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Interpolation
  REAL                   ::  cIpu(g1L:g1U,0:N1)
  REAL                   ::  cIpv(g2L:g2U,0:N2)
  REAL                   ::  cIpw(g3L:g3U,0:N3)
  
  REAL                   ::  cIup(d1L:d1U,0:N1)
  REAL                   ::  cIvp(d2L:d2U,0:N2)
  REAL                   ::  cIwp(d3L:d3U,0:N3)
  
  REAL                   ::  cIcu(g1L:g1U,0:N1,1:n_conc)
  REAL                   ::  cIcv(g2L:g2U,0:N2,1:n_conc)
  REAL                   ::  cIcw(g3L:g3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief Filter
  REAL                   ::  cFp1(b1L:b1U,0:N1)
  REAL                   ::  cFp2(b2L:b2U,0:N2)
  REAL                   ::  cFp3(b3L:b3U,0:N3)
  
  REAL                   ::  cFu1(b1L:b1U,0:N1)
  REAL                   ::  cFv2(b2L:b2U,0:N2)
  REAL                   ::  cFw3(b3L:b3U,0:N3)
  
  REAL                   ::  cFc1(b1L:b1U,0:N1,1:n_conc)
  REAL                   ::  cFc2(b2L:b2U,0:N2,1:n_conc)
  REAL                   ::  cFc3(b3L:b3U,0:N3,1:n_conc)
  !> \}
  
  !> \{ \brief --- Integrator (only for Druckgitter)
  REAL                   ::  cInt1(b1L:b1U,0:N1)
  REAL                   ::  cInt2(b2L:b2U,0:N2)
  REAL                   ::  cInt3(b3L:b3U,0:N3)
  !> \}
  
  
  !> \{ \brief Compact
  INTEGER, PARAMETER     ::  dimS1 = 2*ndL*NB1
  INTEGER, PARAMETER     ::  dimS2 = 2*ndL*NB2
  INTEGER, PARAMETER     ::  dimS3 = 2*ndL*NB3
  
  REAL                   ::  buffer1A(0:N2,0:N3,1:2*ndL), buffer1B(0:N2,0:N3,1:dimS1)
  REAL                   ::  buffer2A(0:N1,0:N3,1:2*ndL), buffer2B(0:N1,0:N3,1:dimS2)
  REAL                   ::  buffer3A(0:N1,0:N2,1:2*ndL), buffer3B(0:N1,0:N2,1:dimS3)
  !----------------------------------------------------
  INTEGER                ::  disp1(1:NB1), recv1(1:NB1)
  INTEGER                ::  disp2(1:NB2), recv2(1:NB2)
  INTEGER                ::  disp3(1:NB3), recv3(1:NB3)
  !----------------------------------------------------
  !> \}
  
  !> \{ \brief Schur complement (square) matrices
  REAL                   ::  Sp1(1:dimS1,1:dimS1), Su1(1:dimS1,1:dimS1), Sc1(1:dimS1,1:dimS1,1:n_conc), Sp11(1:dimS1,1:dimS1), Su11(1:dimS1,1:dimS1), Sc11(1:dimS1,1:dimS1,1:n_conc)
  REAL                   ::  Sp2(1:dimS2,1:dimS2), Sv2(1:dimS2,1:dimS2), Sc2(1:dimS2,1:dimS2,1:n_conc), Sp22(1:dimS2,1:dimS2), Sv22(1:dimS2,1:dimS2), Sc22(1:dimS2,1:dimS2,1:n_conc)
  REAL                   ::  Sp3(1:dimS3,1:dimS3), Sw3(1:dimS3,1:dimS3), Sc3(1:dimS3,1:dimS3,1:n_conc), Sp33(1:dimS3,1:dimS3), Sw33(1:dimS3,1:dimS3), Sc33(1:dimS3,1:dimS3,1:n_conc)
  
  REAL                   ::  SDu1(1:dimS1,1:dimS1), SGp1(1:dimS1,1:dimS1), SIpu(1:dimS1,1:dimS1), SIup(1:dimS1,1:dimS1), SIcu(1:dimS1,1:dimS1,1:n_conc)
  REAL                   ::  SDv2(1:dimS2,1:dimS2), SGp2(1:dimS2,1:dimS2), SIpv(1:dimS2,1:dimS2), SIvp(1:dimS2,1:dimS2), SIcv(1:dimS2,1:dimS2,1:n_conc)
  REAL                   ::  SDw3(1:dimS3,1:dimS3), SGp3(1:dimS3,1:dimS3), SIpw(1:dimS3,1:dimS3), SIwp(1:dimS3,1:dimS3), SIcw(1:dimS3,1:dimS3,1:n_conc)
  
  REAL                   ::  SDu1T(1:dimS1,1:dimS1), SGp1T(1:dimS1,1:dimS1)
  REAL                   ::  SDv2T(1:dimS2,1:dimS2), SGp2T(1:dimS2,1:dimS2)
  REAL                   ::  SDw3T(1:dimS3,1:dimS3), SGp3T(1:dimS3,1:dimS3)
  !> \}
  
  !> \{ \brief bottom (horizontally oriented) matrices
  REAL                   ::  Wp1(1:2*ndL,0:N1), Wu1(1:2*ndL,0:N1), Wc1(1:2*ndL,0:N1,1:n_conc), Wp11(1:2*ndL,0:N1), Wu11(1:2*ndL,0:N1), Wc11(1:2*ndL,0:N1,1:n_conc)
  REAL                   ::  Wp2(1:2*ndL,0:N2), Wv2(1:2*ndL,0:N2), Wc2(1:2*ndL,0:N2,1:n_conc), Wp22(1:2*ndL,0:N2), Wv22(1:2*ndL,0:N2), Wc22(1:2*ndL,0:N2,1:n_conc)
  REAL                   ::  Wp3(1:2*ndL,0:N3), Ww3(1:2*ndL,0:N3), Wc3(1:2*ndL,0:N3,1:n_conc), Wp33(1:2*ndL,0:N3), Ww33(1:2*ndL,0:N3), Wc33(1:2*ndL,0:N3,1:n_conc)
  
  REAL                   ::  WDu1(1:2*ndL,0:N1), WGp1(1:2*ndL,0:N1), WIpu(1:2*ndL,0:N1), WIup(1:2*ndL,0:N1), WIcu(1:2*ndL,0:N1,1:n_conc)
  REAL                   ::  WDv2(1:2*ndL,0:N2), WGp2(1:2*ndL,0:N2), WIpv(1:2*ndL,0:N2), WIvp(1:2*ndL,0:N2), WIcv(1:2*ndL,0:N2,1:n_conc)
  REAL                   ::  WDw3(1:2*ndL,0:N3), WGp3(1:2*ndL,0:N3), WIpw(1:2*ndL,0:N3), WIwp(1:2*ndL,0:N3), WIcw(1:2*ndL,0:N3,1:n_conc)
  
  REAL                   ::  WDu1T(1:2*ndL,0:N1), WGp1T(1:2*ndL,0:N1)
  REAL                   ::  WDv2T(1:2*ndL,0:N2), WGp2T(1:2*ndL,0:N2)
  REAL                   ::  WDw3T(1:2*ndL,0:N3), WGp3T(1:2*ndL,0:N3)
  !> \}
  
  !> \{ \brief implicit finite difference coefficients ("L"=left-hand side)
  REAL                   ::  cp1CL (-ndL:ndL,0:N1), cu1CL(-ndL:ndL,0:N1), cc1CL(-ndL:ndL,0:N1,1:n_conc), cp11CL(-ndL:ndL,0:N1), cu11CL(-ndL:ndL,0:N1), cc11CL(-ndL:ndL,0:N1,1:n_conc)
  REAL                   ::  cp2CL (-ndL:ndL,0:N2), cv2CL(-ndL:ndL,0:N2), cc2CL(-ndL:ndL,0:N2,1:n_conc), cp22CL(-ndL:ndL,0:N2), cv22CL(-ndL:ndL,0:N2), cc22CL(-ndL:ndL,0:N2,1:n_conc)
  REAL                   ::  cp3CL (-ndL:ndL,0:N3), cw3CL(-ndL:ndL,0:N3), cc3CL(-ndL:ndL,0:N3,1:n_conc), cp33CL(-ndL:ndL,0:N3), cw33CL(-ndL:ndL,0:N3), cc33CL(-ndL:ndL,0:N3,1:n_conc)
  
  REAL                   ::  cDu1CL(-ndL:ndL,0:N1), cGp1CL(-ndL:ndL,0:N1), cIpuCL(-ndL:ndL,0:N1), cIupCL(-ndL:ndL,0:N1), cIcuCL(-ndL:ndL,0:N1,1:n_conc)
  REAL                   ::  cDv2CL(-ndL:ndL,0:N2), cGp2CL(-ndL:ndL,0:N2), cIpvCL(-ndL:ndL,0:N2), cIvpCL(-ndL:ndL,0:N2), cIcvCL(-ndL:ndL,0:N2,1:n_conc)
  REAL                   ::  cDw3CL(-ndL:ndL,0:N3), cGp3CL(-ndL:ndL,0:N3), cIpwCL(-ndL:ndL,0:N3), cIwpCL(-ndL:ndL,0:N3), cIcwCL(-ndL:ndL,0:N3,1:n_conc)
  
  REAL                   ::  cDu1CLT(-ndL:ndL,0:N1), cGp1CLT(-ndL:ndL,0:N1)
  REAL                   ::  cDv2CLT(-ndL:ndL,0:N2), cGp2CLT(-ndL:ndL,0:N2)
  REAL                   ::  cDw3CLT(-ndL:ndL,0:N3), cGp3CLT(-ndL:ndL,0:N3)
  
  REAL                   ::  cp1CL_LU(1:(3*ndL+1),0:N1), cu1CL_LU(1:(3*ndL+1),0:N1), cc1CL_LU(1:(3*ndL+1),0:N1,1:n_conc), cp11CL_LU(1:(3*ndL+1),0:N1), cu11CL_LU(1:(3*ndL+1),0:N1), cc11CL_LU(1:(3*ndL+1),0:N1,1:n_conc)
  REAL                   ::  cp2CL_LU(1:(3*ndL+1),0:N2), cv2CL_LU(1:(3*ndL+1),0:N2), cc2CL_LU(1:(3*ndL+1),0:N2,1:n_conc), cp22CL_LU(1:(3*ndL+1),0:N2), cv22CL_LU(1:(3*ndL+1),0:N2), cc22CL_LU(1:(3*ndL+1),0:N2,1:n_conc)
  REAL                   ::  cp3CL_LU(1:(3*ndL+1),0:N3), cw3CL_LU(1:(3*ndL+1),0:N3), cc3CL_LU(1:(3*ndL+1),0:N3,1:n_conc), cp33CL_LU(1:(3*ndL+1),0:N3), cw33CL_LU(1:(3*ndL+1),0:N3), cc33CL_LU(1:(3*ndL+1),0:N3,1:n_conc)
  
  REAL                   ::  cDu1CL_LU(1:(3*ndL+1),0:N1), cGp1CL_LU(1:(3*ndL+1),0:N1), cIpuCL_LU(1:(3*ndL+1),0:N1), cIupCL_LU(1:(3*ndL+1),0:N1), cIcuCL_LU(1:(3*ndL+1),0:N1,1:n_conc)
  REAL                   ::  cDv2CL_LU(1:(3*ndL+1),0:N2), cGp2CL_LU(1:(3*ndL+1),0:N2), cIpvCL_LU(1:(3*ndL+1),0:N2), cIvpCL_LU(1:(3*ndL+1),0:N2), cIcvCL_LU(1:(3*ndL+1),0:N2,1:n_conc)
  REAL                   ::  cDw3CL_LU(1:(3*ndL+1),0:N3), cGp3CL_LU(1:(3*ndL+1),0:N3), cIpwCL_LU(1:(3*ndL+1),0:N3), cIwpCL_LU(1:(3*ndL+1),0:N3), cIcwCL_LU(1:(3*ndL+1),0:N3,1:n_conc)
  
  REAL                   ::  cDu1CLT_LU(1:(3*ndL+1),0:N1), cGp1CLT_LU(1:(3*ndL+1),0:N1)
  REAL                   ::  cDv2CLT_LU(1:(3*ndL+1),0:N2), cGp2CLT_LU(1:(3*ndL+1),0:N2)
  REAL                   ::  cDw3CLT_LU(1:(3*ndL+1),0:N3), cGp3CLT_LU(1:(3*ndL+1),0:N3)
  !> \}
  
  !> \{ \brief explicit finite difference coefficients ("R"=right-hand side)
  REAL                   ::  cp1CR(-ndR:ndR,0:N1), cu1CR(-ndR:ndR,0:N1), cc1CR(-ndR:ndR,0:N1,1:n_conc), cp11CR(-ndR:ndR,0:N1), cu11CR(-ndR:ndR,0:N1), cc11CR(-ndR:ndR,0:N1,1:n_conc)
  REAL                   ::  cp2CR(-ndR:ndR,0:N2), cv2CR(-ndR:ndR,0:N2), cc2CR(-ndR:ndR,0:N2,1:n_conc), cp22CR(-ndR:ndR,0:N2), cv22CR(-ndR:ndR,0:N2), cc22CR(-ndR:ndR,0:N2,1:n_conc)
  REAL                   ::  cp3CR(-ndR:ndR,0:N3), cw3CR(-ndR:ndR,0:N3), cc3CR(-ndR:ndR,0:N3,1:n_conc), cp33CR(-ndR:ndR,0:N3), cw33CR(-ndR:ndR,0:N3), cc33CR(-ndR:ndR,0:N3,1:n_conc)
  
  REAL                   ::  cDu1CR(-ndR:ndR,0:N1), cGp1CR(-ndR:ndR,0:N1), cIpuCR(-ndR:ndR,0:N1), cIupCR(-ndR:ndR,0:N1), cIcuCR(-ndR:ndR,0:N1,1:n_conc)
  REAL                   ::  cDv2CR(-ndR:ndR,0:N2), cGp2CR(-ndR:ndR,0:N2), cIpvCR(-ndR:ndR,0:N2), cIvpCR(-ndR:ndR,0:N2), cIcvCR(-ndR:ndR,0:N2,1:n_conc)
  REAL                   ::  cDw3CR(-ndR:ndR,0:N3), cGp3CR(-ndR:ndR,0:N3), cIpwCR(-ndR:ndR,0:N3), cIwpCR(-ndR:ndR,0:N3), cIcwCR(-ndR:ndR,0:N3,1:n_conc)
  
  REAL                   ::  cDu1CRT(-ndR:ndR,0:N1), cGp1CRT(-ndR:ndR,0:N1)
  REAL                   ::  cDv2CRT(-ndR:ndR,0:N2), cGp2CRT(-ndR:ndR,0:N2)
  REAL                   ::  cDw3CRT(-ndR:ndR,0:N3), cGp3CRT(-ndR:ndR,0:N3)
  !> \}
  !****************************************************
  
  
  !> \{ \brief 2. Ableitung (Multigrid).
  !! Anmerkung: Die Koeffizientens�tze unterscheiden sich z.T. lediglich durch die Randbedingungen.
  REAL                   ::  cp11R(-1:1,0:N1,1:n_grids_max)
  REAL                   ::  cp22R(-1:1,0:N2,1:n_grids_max)
  REAL                   ::  cp33R(-1:1,0:N3,1:n_grids_max)
  
  REAL                   ::  cu11R(-1:1,0:N1,1:n_grids_max)
  REAL                   ::  cv22R(-1:1,0:N2,1:n_grids_max)
  REAL                   ::  cw33R(-1:1,0:N3,1:n_grids_max)
  
  REAL                   ::  cc11R(-1:1,0:N1,1:n_grids_max,1:n_conc)
  REAL                   ::  cc22R(-1:1,0:N2,1:n_grids_max,1:n_conc)
  REAL                   ::  cc33R(-1:1,0:N3,1:n_grids_max,1:n_conc)
  
  REAL                   ::  cdg1 (-1:1,1:N1,1:n_grids_max)
  REAL                   ::  cdg2 (-1:1,1:N2,1:n_grids_max)
  REAL                   ::  cdg3 (-1:1,1:N3,1:n_grids_max)
  !> \}
  
  !> \{ \brief Interpolation (Multigrid)
  REAL                   ::  cI1(1:2,1:N1,1:n_grids_max)
  REAL                   ::  cI2(1:2,1:N2,1:n_grids_max)
  REAL                   ::  cI3(1:2,1:N3,1:n_grids_max)
  
  REAL                   ::  cIH1(1:2,0:N1)
  REAL                   ::  cIH2(1:2,0:N2)
  REAL                   ::  cIH3(1:2,0:N3)
  !> \}
  
  !> \{ \brief Restriktion (Multigrid)
  REAL                   ::  cR1 (-1:1,1:N1,2:n_grids_max)
  REAL                   ::  cR2 (-1:1,1:N2,2:n_grids_max)
  REAL                   ::  cR3 (-1:1,1:N3,2:n_grids_max)
  
  REAL                   ::  cRest1(b1L:b1U,0:N1,1:n_grids_max-1) ! TEST!!!
  REAL                   ::  cRest2(b2L:b2U,0:N2,1:n_grids_max-1)
  REAL                   ::  cRest3(b3L:b3U,0:N3,1:n_grids_max-1)
  
  REAL                   ::  cRH1(1:2,1:N1)
  REAL                   ::  cRH2(1:2,1:N2)
  REAL                   ::  cRH3(1:2,1:N3)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Gitterspezifikationen
  !===========================================================================================================
  !> \{ \brief physiklische Koordinaten (global)
  REAL                   ::  y1p(1:M1), y1u(0:M1)
  REAL                   ::  y2p(1:M2), y2v(0:M2)
  REAL                   ::  y3p(1:M3), y3w(0:M3)
  !> \}
  
  !> \{ \brief physiklische Koordinaten (Block)
  REAL                   ::  x1p(b1L:(N1+b1U)), x1u(b1L:(N1+b1U))
  REAL                   ::  x2p(b2L:(N2+b2U)), x2v(b2L:(N2+b2U))
  REAL                   ::  x3p(b3L:(N3+b3U)), x3w(b3L:(N3+b3U))
  !> \}
  
  !> \{ \brief --- physiklische Koordinaten (Block, Multigrid)
  REAL                   ::  x1pR(b1L:(N1+b1U),1:n_grids_max), x1uR(b1L:(N1+b1U),1:n_grids_max)
  REAL                   ::  x2pR(b2L:(N2+b2U),1:n_grids_max), x2vR(b2L:(N2+b2U),1:n_grids_max)
  REAL                   ::  x3pR(b3L:(N3+b3U),1:n_grids_max), x3wR(b3L:(N3+b3U),1:n_grids_max)
  !> \}
  
  !> \{ \brief Gitterweiten (global)
  REAL                   ::  dy1p(1:M1), dy1u(0:M1)
  REAL                   ::  dy2p(1:M2), dy2v(0:M2)
  REAL                   ::  dy3p(1:M3), dy3w(0:M3)
  !> \}
  
  !> \{ \brief Gitterweiten (Block)
  REAL                   ::  dx1p(1:N1), dx1u(0:N1)
  REAL                   ::  dx2p(1:N2), dx2v(0:N2)
  REAL                   ::  dx3p(1:N3), dx3w(0:N3)
  
  REAL                   ::  dx1DM(1:N1), dx1pM(1:N1), ddx1pM(1:N1)
  REAL                   ::  dx2DM(1:N2), dx2pM(1:N2), ddx2pM(1:N2)
  REAL                   ::  dx3DM(1:N3), dx3pM(1:N3), ddx3pM(1:N3)
  
  REAL                   ::  dx1GM(0:N1), dx1uM(0:N1), ddx1uM(0:N1)
  REAL                   ::  dx2GM(0:N2), dx2vM(0:N2), ddx2vM(0:N2)
  REAL                   ::  dx3GM(0:N3), dx3wM(0:N3), ddx3wM(0:N3)
  !> \}
  
  !> \{ \brief Smagorinsky-Modell
  REAL                   ::  dx1pS(1:N1), dx1uS(0:N1) ! TEST!!!
  REAL                   ::  dx2pS(1:N2), dx2vS(0:N2)
  REAL                   ::  dx3pS(1:N3), dx3wS(0:N3)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name Arbeitsfelder
  !===========================================================================================================
  !> \brief Geschwindigkeiten -------------------------------------------------------------------------------------
  REAL                   ::  vel(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief nicht-linearer Term -----------------------------------------------------------------------------------
  REAL                   ::  nl (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief Recht-Hand-Seite --------------------------------------------------------------------------------------
  REAL                   ::  rhs(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  !> \brief Druck -------------------------------------------------------------------------------------------------
  REAL                   ::  pre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Konzentrationsfelder ----------------------------------------------------------------------------------
  REAL                   ::  conc(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc)
  
  !> \brief nicht-linearer Term (Konzentration) -------------------------------------------------------------------
  REAL                   ::  nlco(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:n_conc)
  
  !> \{ \brief Ausfluss-RB (Geschwindigkeitsfeld)
  !!
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehürigen Randbedingungen gespeichert werden.
  REAL                   ::  bc11(1:N2,1:N3,1:2), nlbc11(1:N2,1:N3,1:2)
  REAL                   ::  bc12(0:N1,1:N3,1:2), nlbc12(0:N1,1:N3,1:2)
  REAL                   ::  bc13(0:N1,1:N2,1:2), nlbc13(0:N1,1:N2,1:2)
  
  REAL                   ::  bc21(0:N2,1:N3,1:2), nlbc21(0:N2,1:N3,1:2)
  REAL                   ::  bc22(1:N1,1:N3,1:2), nlbc22(1:N1,1:N3,1:2)
  REAL                   ::  bc23(1:N1,0:N2,1:2), nlbc23(1:N1,0:N2,1:2)
  
  REAL                   ::  bc31(1:N2,0:N3,1:2), nlbc31(1:N2,0:N3,1:2)
  REAL                   ::  bc32(1:N1,0:N3,1:2), nlbc32(1:N1,0:N3,1:2)
  REAL                   ::  bc33(1:N1,1:N2,1:2), nlbc33(1:N1,1:N2,1:2)
  
  REAL                   ::  drift1(b2L:(N2+b2U),b3L:(N3+b3U),1:2)
  REAL                   ::  drift2(b1L:(N1+b1U),b3L:(N3+b3U),1:2)
  REAL                   ::  drift3(b1L:(N1+b1U),b2L:(N2+b2U),1:2)
  !>\}

  !> \{ \brief Sedimentations-RB (Konzentration)
  !!
  !! Da die RHS für die Konzentrationsfelder nicht über die Runge-Kutta-Zwischenschritte hinweg gespeichert
  !! werden, müssen mindestens die zugehürigen Randbedingungen gespeichert werden.
  REAL                   ::  sed1L (1:N2,1:N3,1:n_conc)
  REAL                   ::  sed1U (1:N2,1:N3,1:n_conc)
  
  REAL                   ::  sed2L (1:N1,1:N3,1:n_conc)
  REAL                   ::  sed2U (1:N1,1:N3,1:n_conc)
  
  REAL                   ::  sed3L (1:N1,1:N2,1:n_conc)
  REAL                   ::  sed3U (1:N1,1:N2,1:n_conc)
  
  REAL                   ::  conc1L(1:N2,1:N3,1:n_conc)
  REAL                   ::  conc1U(1:N2,1:N3,1:n_conc)
  
  REAL                   ::  conc2L(1:N1,1:N3,1:n_conc)
  REAL                   ::  conc2U(1:N1,1:N3,1:n_conc)
  
  REAL                   ::  conc3L(1:N1,1:N2,1:n_conc)
  REAL                   ::  conc3U(1:N1,1:N2,1:n_conc)
  !> \}
  
  !> \{ \brief Sedimentation (Konzentration, Ausschrieb)
  REAL                   ::  dep1L_conc(1:N2,1:N3,1:n_conc)
  REAL                   ::  dep1U_conc(1:N2,1:N3,1:n_conc)
  
  REAL                   ::  dep2L_conc(1:N1,1:N3,1:n_conc)
  REAL                   ::  dep2U_conc(1:N1,1:N3,1:n_conc)
  
  REAL                   ::  dep3L_conc(1:N1,1:N2,1:n_conc)
  REAL                   ::  dep3U_conc(1:N1,1:N2,1:n_conc)
  !> \}
  
  !> \{ \brief Sedimentation (Partikel, Ausschrieb).
  !! TEST!!! 1:n_spec???
  REAL                   ::  dep1L_part(1:N2,1:N3)
  REAL                   ::  dep1U_part(1:N2,1:N3)
  
  REAL                   ::  dep2L_part(1:N1,1:N3)
  REAL                   ::  dep2U_part(1:N1,1:N3)
  
  REAL                   ::  dep3L_part(1:N1,1:N2)
  REAL                   ::  dep3U_part(1:N1,1:N2)
  !> \}

  !> \brief Residuum ----------------------------------------------------------------------------------------------
  REAL                   ::  res (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Druckgradient (eine Komponente) -----------------------------------------------------------------------
  REAL                   ::  gpre(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Gewichte for Divergenzfreiheit ------------------------------------------------------------------------
  REAL                   ::  weight(1:N1,1:N2,1:N3)
  
#ifdef NONBOUSSINESQ
  !> \brief fluid density on velocity grid ------------------------------------------------------------------------
  REAL                   ::  dens(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
#endif
  
  !> \{ \brief Null-Raum-Vektor
  REAL                   ::  psi      (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  psi_vel  (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U),1:3)
  
  REAL                   ::  psi_rel1 (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL   , ALLOCATABLE   ::  psi_rel2 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel3 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel4 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel5 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel6 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel7 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel8 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel9 (:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel10(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel11(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel12(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel13(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel14(:,:,:)
  REAL   , ALLOCATABLE   ::  psi_rel15(:,:,:)
  
  REAL                   ::  th11(1:N2,1:N3,1:2)
  REAL                   ::  th12(0:N1,1:N3,1:2)
  REAL                   ::  th13(0:N1,1:N2,1:2)
  
  REAL                   ::  th21(0:N2,1:N3,1:2)
  REAL                   ::  th22(1:N1,1:N3,1:2)
  REAL                   ::  th23(1:N1,0:N2,1:2)
  
  REAL                   ::  th31(1:N2,0:N3,1:2)
  REAL                   ::  th32(1:N1,0:N3,1:2)
  REAL                   ::  th33(1:N1,1:N2,1:2)
  !> \}
  
  !> \{ \brief Multigrid
  REAL                   ::  vec1C (b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  REAL   , ALLOCATABLE   ::  vec2A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec2C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec3A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec3C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec4A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec4C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec5A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec5C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec6A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec6C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec7A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec7C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec8A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec8C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec9A (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9B (:,:,:)
  REAL   , ALLOCATABLE   ::  vec9C (:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec10A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec10C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec11A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec11C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec12A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec12C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec13A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec13C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec14A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14B(:,:,:)
  REAL   , ALLOCATABLE   ::  vec14C(:,:,:)
  
  REAL   , ALLOCATABLE   ::  vec15A(:,:,:)
  REAL   , ALLOCATABLE   ::  vec15B(:,:,:)
  !> \}
  
  
  !> \{ \brief BiCGstab / Richardson
  REAL                   ::  pp(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  Ap(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rr(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  rh(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  Ar(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  z2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !> \}
  
  !> \brief product_div_grad --------------------------------------------------------------------------------------
  REAL                   ::  dig(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \brief Hilfsfeld (kompakte Differenzen) ----------------------------------------------------------------------
  REAL                   ::  com(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  
  !> \{ \brief Hilfsfelder (Druckiteration).
  !! - rhs wird auch in test_moment nochmals verwendet und kann daher in outer_iteration nicht belegt werden!
  !! - wird auch fuer interpolierte Geschwindigkeiten in rhs_vel und rhs_conc verwendet.
  REAL                   ::  work1(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  work2(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  REAL                   ::  work3(b1L:(N1+b1U),b2L:(N2+b2U),b3L:(N3+b3U))
  !> \}
  
  !> \{ \brief Lagrange-Partikel
  INTEGER                ::  n_part, n_part_tot, n_groups
  REAL                   ::  particles(1:n_args,1:n_part_max)
  !> \}
  
  !> \{ \brief Smagorinsky LES-Modell
  REAL   , ALLOCATABLE   ::  smag_diffus(:,:,:,:)
  !> \}
  
  !> \{ \brief Linienrelaxation
  REAL                   ::  vec1(1:N1), dia1(1:N1), SOR1(1:N1), band1(1:2,1:N1) ! TEST!!! vec1, dia muessen auch in mod_Helmholtz ausgetauscht werden!
  REAL                   ::  vec2(1:N2), dia2(1:N2), SOR2(1:N2), band2(1:2,1:N2) ! TEST!!! SOR muesste man idealerweise auch in band eingliedern!
  REAL                   ::  vec3(1:N3), dia3(1:N3), SOR3(1:N3), band3(1:2,1:N3) ! dia3(:) --> band3(1,:), vec3(:) --> band3(2,:), etc. ...
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \brief Indizierung (Intervallgrenzen, Verschiebungen)
  !===========================================================================================================
  !> \{ \brief Block-Index
  !INTEGER                ::  iB, jB, kB ! TEST!!!
  
  !> \} \brief Indexverschiebung (Block --> global)
  INTEGER                ::  iShift, jShift, kShift
  !> \}
  
  !> \{ \brief Domaingrösse (Periodizität-bereinigt)
  INTEGER                ::  dim1, dim2, dim3
  !> \}
  
  !> \{ \brief Druck / Konzentrationen (inklusive Rand)
  INTEGER                ::  S1p, S2p, S3p
  INTEGER                ::  N1p, N2p, N3p
  !> \}
  
  !> \{ \brief Geschwindigkeiten (inklusive Rand)
  !!
  !! start values for the velocity including the boundary value
  !! first int indicates the spatial direction, the second in denotes the velocity direction
  INTEGER                ::  S11B, S21B, S31B
  INTEGER                ::  S12B, S22B, S32B
  INTEGER                ::  S13B, S23B, S33B
  !> \}
  
  !> \{ \brief end values for the velocity including the boundary value
  INTEGER                ::  N11B, N21B, N31B
  INTEGER                ::  N12B, N22B, N32B
  INTEGER                ::  N13B, N23B, N33B
  !> \}
  
  !> \{ \brief Geschwindigkeiten (exklusive Rand)
  !! start and end indices for the velocity without the boundary nodes
  INTEGER                ::  S11, S21, S31
  INTEGER                ::  S12, S22, S32
  INTEGER                ::  S13, S23, S33
  
  INTEGER                ::  N11, N21, N31
  INTEGER                ::  N12, N22, N32
  INTEGER                ::  N13, N23, N33
  !> \}
  
  !> \{ \brief Konzentrationen (exklusive Rand)
  INTEGER                ::  S1c(1:n_conc), S2c(1:n_conc), S3c(1:n_conc)
  INTEGER                ::  N1c(1:n_conc), N2c(1:n_conc), N3c(1:n_conc)
  !> \}
  
  !> \{ \brief grobe Gitter (Multigrid, INklusive Rand)
  INTEGER                ::  S1R, S2R, S3R
  INTEGER                ::  d1R, d2R, d3R
  !> \}
  
  !> \{ \brief grobe Gitter (Multigrid, EXklusive Rand)
  INTEGER                ::  S11R, S22R, S33R
  INTEGER                ::  d11R, d22R, d33R
  !> \}
  
  !> \{ \brief Ueberlappungskonvention der Blöcke (Multigrid, siehe mod_setup)
  INTEGER, PARAMETER     ::  ls1 = -1
  INTEGER, PARAMETER     ::  ls2 = -1
  INTEGER, PARAMETER     ::  ls3 = -1
  !> \}
  
  !> \{ \brief Austauschrichtung (Multigrid)
  !!
  !! ex = -1: unten <--  oben
  !! ex =  0: unten <--> oben
  !! ex =  1: unten  --> oben
  INTEGER                ::  ex1, ex2, ex3
  !\}
  
  
  !\}
  !===========================================================================================================
  !> \{ \name Randbedingungen
  !!                              _
  !!    Symmetrie-RB:   BC = -2    |
  !!    periodische RB: BC = -1    |- symmetrische, zentrale Stencils
  !!    Nachbar-Block:  BC =  0   _|
  !!    Dirichlet-RB:   BC =  1    |
  !!    Neumann-RB:     BC =  2    |- schiefe, nicht-zentrale Stencils
  !!    Robin-RB:       BC =  3   _|

  !> \{ \brief global
  LOGICAL                ::  outlet   (1:3,1:2,1:3)
  LOGICAL                ::  isopycnal(1:3,1:2,1:n_conc)
  
  INTEGER                ::  BC_1L_global, BC_1U_global
  INTEGER                ::  BC_2L_global, BC_2U_global
  INTEGER                ::  BC_3L_global, BC_3U_global
  !> \}
  
  !> \{ \brief lokal (Block)
  INTEGER                ::  BC_1L, BC_1U
  INTEGER                ::  BC_2L, BC_2U
  INTEGER                ::  BC_3L, BC_3U
  
  INTEGER                ::  BCc_1L(1:n_conc), BCc_1U(1:n_conc)
  INTEGER                ::  BCc_2L(1:n_conc), BCc_2U(1:n_conc)
  INTEGER                ::  BCc_3L(1:n_conc), BCc_3U(1:n_conc)
  !> \}
  
  !> \{ \brief field properties
  INTEGER                ::  n_gather(1:3,1:n_grids_max)
  INTEGER                ::  NN (1:3,1:n_grids_max)
  INTEGER                ::  NB (1:3,1:n_grids_max)
  INTEGER                ::  iB (1:3,1:n_grids_max)
  INTEGER                ::  SNF(1:2,1:3,1:n_grids_max)
  INTEGER                ::  SNB(1:2,1:3,1:n_grids_max)
  INTEGER                ::  BC (1:2,1:3,1:n_grids_max)
  INTEGER                ::  ngb(1:2,1:3,1:n_grids_max)
  INTEGER                ::  comm1(1:n_grids_max), comm2(1:n_grids_max)
  INTEGER                ::  rankc2(1:n_grids_max)
  LOGICAL                ::  participate_yes(1:n_grids_max)
  INTEGER                ::  recvR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  recvI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispR(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  dispI(    1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  offsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsR(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  INTEGER                ::  sizsI(1:3,1:NB1*NB2*NB3,1:n_grids_max)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name physikalische Parameter
  !===========================================================================================================
  REAL                   ::  L1, L2, L3
  REAL                   ::  Re, Fro
  REAL                   ::  Sc (1:n_conc)
  REAL                   ::  Ric(1:n_conc)
  REAL                   ::  Rip(1:n_spec) ! TEST!!!
  REAL                   ::  usc(1:n_conc)
  REAL                   ::  usp(1:n_spec) ! TEST!!!
  REAL                   ::  Stp(1:n_spec) ! TEST!!!
  
  REAL                   ::  velReSc1(1:N2,1:N3,1:2,1:n_grids_max)
  REAL                   ::  velReSc2(1:N1,1:N3,1:2,1:n_grids_max)
  REAL                   ::  velReSc3(1:N1,1:N2,1:2,1:n_grids_max)
  
  REAL                   ::  usReSc (1:3,1:n_conc)
  REAL                   ::  us_vec (1:3,1:n_conc)
  REAL                   ::  gravity(1:3)
  
  REAL                   ::  dx_1L, dx_1U
  REAL                   ::  dx_2L, dx_2U
  REAL                   ::  dx_3L, dx_3U
  
  REAL                   ::  idx_1L, idx_1U
  REAL                   ::  idx_2L, idx_2U
  REAL                   ::  idx_3L, idx_3U
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name numerische Parameter
  !===========================================================================================================
  !> \{ \brief allgemein
  REAL                   ::  CFL
  REAL                   ::  time, dtime, subtime, time_start, time_end, dtime_max, dtime0, dtime_old
  INTEGER                ::  timestep, timestep_old, substep, n_timesteps
  LOGICAL                ::  mapping_yes, upwind_yes, upwind_conc_yes
  LOGICAL                ::  Euler_yes, nonBoussinesq_yes, Stokes_yes, twostep_yes
  LOGICAL                ::  comp_visc_yes, comp_conv_yes, comp_inter_yes, comp_div_yes, comp_grad_yes
  LOGICAL, PARAMETER     ::  filter_BC_yes = .TRUE. ! TEST!!!
  !> timeint_mode = 0 (CN-RK3)
  !! timeint_mode = 1 (RK3-O3)
  INTEGER                ::  timeint_mode, LES_mode, forcing_mode
  INTEGER                ::  bulkflow_dir
  INTEGER                ::  n_lp_vel           , n_hp_vel
  INTEGER                ::  n_lp_conc(1:n_conc), n_hp_conc(1:n_conc)
  REAL                   ::  chi_vel
  REAL                   ::  chi_conc (1:n_conc)
  REAL                   ::  vel_bulk ! TEST!!!
  !> \}
  
  
  !> \{ \brief Runge-Kutta-Koeffizienten
  REAL   , PARAMETER     ::  aRK(1:3) = (/8./15.,  5./12., 3./ 4./)
  REAL   , PARAMETER     ::  bRK(1:3) = (/  0.  ,-17./60.,-5./12./)
  INTEGER, PARAMETER     ::  RK_steps = 3
  !> \}
  
  !> \{ \brief look-up table fuer Stabilitaetsgebiet der Zeitintegration (angle = pi/2,pi)
  REAL   , PARAMETER     ::  stabilitylimit(0:40) = (/1.732050813, 1.943689093, 2.089210537, 2.201001743,  &
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
  INTEGER, PARAMETER     ::  n_stab = SIZE(stabilitylimit)
  !> \}
  
  !> \brief Helmholtz-Vorfaktor
  !!
  !! \f[0 \le \theta_L \le 1 \f]
  !! 0 means fully explicit and 1 means fully implicit
  REAL                   ::  thetaL
  !> \brief Helmholtz-Vorfaktor
  !!
  !! \f[ \mathrm{multL = -(1-thetaL)(aRK(substep)+bRK(substep))\frac{\Delta t}{Re} }\f]
  !! computed before mod_diff::Helmholtz call
  REAL                   ::  multL
  
  !> \{ \brief zeitliche Steuerung
  INTEGER                ::  Int_dtime, Int_lev_pre
  
  INTEGER                ::  stride_large(1:3), stride_med(1:3), stride_small(1:3)
  LOGICAL                ::  write_large, write_med, write_small
  REAL                   ::  time_out_scal, dtime_out_scal
  REAL                   ::  time_out_vect, dtime_out_vect
  
  LOGICAL                ::  write_out_scal, write_out_vect  
  LOGICAL                ::  new_dtime, finish_yes
    
  INTEGER                ::  write_count
  INTEGER                ::  restart
  CHARACTER(LEN=3)       ::  restart_char
  
  INTEGER                ::  n_conc_old
  !> \}

  !===========================================================================================================
  !> \{ \brief weitere Steuerungsoptionen
  !===========================================================================================================
  INTEGER                ::  task
  LOGICAL                ::  read_nullspace_yes
  LOGICAL                ::  concentration_yes
  LOGICAL                ::  particles_yes
  LOGICAL                ::  nullspace_yes, nullspace_coarse_yes
  LOGICAL                ::  nullspace_ortho_yes
  
  LOGICAL                ::  write_stout_yes
  LOGICAL                ::  log_iteration_yes
  LOGICAL                ::  write_restart_yes
  LOGICAL                ::  write_lambda2_yes
  LOGICAL                ::  write_test_yes
  !> \}
  
  !> \{ \brief globale Laufindizes
  INTEGER                ::  direction
  INTEGER                ::  conc_nu
  !> \}
  
  !> \brief explizite Behandlung von Ecken bei Dirichlet-Randbedingungen
  !!
  !! (Hintergrund: Der Druck ist an diesen Orten unbestimmt, so dass er dort künstlich zu Null gesetzt wird.)
  LOGICAL, PARAMETER     ::  corner_yes = .TRUE.
  
  !> \{ \brief Systemzeit
  INTEGER                ::  elatime, ctime(1:8)
  INTEGER                ::  day, hour, minu, sec, msec
  !> \}
  
  !> \brief Partikel-Energie
  REAL                   ::  Ekin_part
  
  
  !===========================================================================================================
  !> \{ \name Iterationsparameter
  !===========================================================================================================
  !> \{ \brief Abbruchkriterium / Absolute Genauigkeit der Geschwindigkeiten
  REAL                   ::  epsU, epsU0
  !> \}
  
  !> \brief Glaetter
  LOGICAL                ::  Jacobi_yes
  
  !> \{ \brief max. Anzahl Iterationen
  INTEGER                ::  n_it_outer
  INTEGER                ::  n_it_Poisson
  INTEGER                ::  n_it_Helmh_vel
  INTEGER                ::  n_it_Helmh_conc
  !> \}
  
  !> \{ \brief erwartete Konvergenzrate (äussere Iteration)
  REAL                   ::  precRatio0 (1:RK_steps)
  REAL                   ::  precOffset0(1:RK_steps)
  REAL, ALLOCATABLE      ::  precRatio  (:,:)
  REAL, ALLOCATABLE      ::  precOffset (:,:)
  !> \}
  
  !> \{ \brief Null-Initialisierung (äussere Iteration)
  LOGICAL                ::  init_pre(1:RK_steps), init_vel(1:RK_steps), init_conc(1:RK_steps)
  !> \}
  
  !> \{ \brief Vorkonditionierung (Multigrid)
  INTEGER                ::  precond_outer
  INTEGER                ::  precond_Poisson
  INTEGER                ::  precond_Helmh_vel
  INTEGER                ::  precond_Helmh_conc
  !> \}
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  INTEGER                ::  n_relax_down, n_relax_up, n_relax_bottom
  
  !> \brief implizite Richtungen bei Linienrelaxation (Multigrid)
  INTEGER                ::  impl_dir(1:3)
  
  !> \brief Anzahl Glättungen pro Gitterlevel (Multigrid)
  LOGICAL                ::  weighting_yes
  
  !> \}
  !===========================================================================================================
  !> \{ \name Iterationsstatistik
  !===========================================================================================================
  REAL                   ::  dtime_average
  REAL                   ::  max_div_init(1:2)
  INTEGER                ::  number_poisson
  
  !> \{ \brief Zähler
  INTEGER                ::  countO(1:RK_steps)
  INTEGER                ::  countP(1:RK_steps,1:2)
  INTEGER                ::  countH(1:RK_steps,1:3)
  !> \}
  
  !> \{ \brief Konvergenzrate
  REAL                   ::  ratioO(1:RK_steps)
  REAL                   ::  ratioH(1:RK_steps,1:3)
  REAL                   ::  ratioP(1:RK_steps,1:2)
  !> \}
  
  
  !> \}
  !===========================================================================================================
  !> \{ \name MPI
  !===========================================================================================================
  !> \{ \brief Kommunikatoren
  INTEGER                ::  COMM_CART
  
  INTEGER                ::  COMM_SLICE1, COMM_BAR1
  INTEGER                ::  COMM_SLICE2, COMM_BAR2
  INTEGER                ::  COMM_SLICE3, COMM_BAR3
  !> \}
  
  !> \{ \brief Dimension und Position der Blöcke innerhalb der Kommunikatoren (Gitter-Indizes)
  !!
  !! (for MPI_GATHERv, MPI_ALLGATHERv, vgl. iShift, jShift, kShift)
  INTEGER                ::  bar1_size(1:NB1), bar1_offset(1:NB1)
  INTEGER                ::  bar2_size(1:NB2), bar2_offset(1:NB2)
  INTEGER                ::  bar3_size(1:NB3), bar3_offset(1:NB3)
  !> \{
  
  !> \{ \brief Ränge der Prozesse
  INTEGER                ::  rank
  INTEGER                ::  rank_bar1, rank_slice1
  INTEGER                ::  rank_bar2, rank_slice2
  INTEGER                ::  rank_bar3, rank_slice3
  !> \}
  
  !> \{ \brief Ränge der Nachbarprozesse (in kartesischem Gitter)
  INTEGER                ::  rank1L, rank1U
  INTEGER                ::  rank2L, rank2U
  INTEGER                ::  rank3L, rank3U
  !> \}
  
  !> \brief Error-Handle
  INTEGER                ::  merror
  
  !> \{ \brief Request-Handles
  !!
  !! (müssen offenbar global angelegt werden)
  INTEGER                ::  req1L, req1U
  INTEGER                ::  req2L, req2U
  INTEGER                ::  req3L, req3U
  !> \}
  
  !> \}
  !===========================================================================================================
  !> \{ \name HDF5
  !===========================================================================================================
  INTEGER                ::  herror
  !> \}
  
  !> \{ \name periodicity criterions
  REAL   , ALLOCATABLE   ::  velp(:,:,:,:)       !< \brief stored velocity, used as reference to previous period
  REAL,                  ::  tp
  !> \}
#endif
  
  
END MODULE mod_vars
