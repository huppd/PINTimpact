!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* Mai 2005 - Dec 2011                                                                                       *
!*************************************************************************************************************

!> @brief pure module for global variables giving global acces to alot of constants concerning the dimensions
module mod_dims
  
  implicit none
  
  !===========================================================================================================
  !> \{ \name verschiedene Dimensionen
  !===========================================================================================================
#ifdef ALLOC
  integer                ::  n_conc !< amount of different concentrations
  integer                ::  n_spec ! TEST!!! muss nicht unbedingt hier stehen, kann auch nach mod_vars ...
  integer                ::  n_part_max
  integer                ::  n_args
#else
  integer, parameter     ::  n_conc = 1 ! TEST!!!
  integer, parameter     ::  n_spec = 1 ! TEST!!! Tests bzgl. n_conc auch auf n_spec anwenden!
  integer, parameter     ::  n_part_max = 100000 ! >= 1 ! TEST!!! Test schreiben
  integer, parameter     ::  n_args   = 14 ! >= 14
#endif
  !> \}
  
  !===========================================================================================================
  !> \{ \name Domain- und Blockspezifikationen
  !===========================================================================================================
  !> \{ \brief Anzahl Blöcke
#ifdef ALLOC
  integer             ::  NB1
  integer             ::  NB2
  integer             ::  NB3
#else
  integer, parameter  ::  NB1 = 2
  integer, parameter  ::  NB2 = 8
  integer, parameter  ::  NB3 = 4
#endif
  !> \}
  
  !> \{ \brief Gesamte Domain (über alle Blöcke)
#ifdef ALLOC
  integer             ::  M1
  integer             ::  M2
  integer             ::  M3
#else
  integer, parameter  ::  M1 = 3*2**4+1
  integer, parameter  ::  M2 = 7*2**6+1
  integer, parameter  ::  M3 = 5*2**5+1
#endif
  !> \}
  !! \}
  
  
  !===========================================================================================================
  !> \{ \name Konvergenzordnung der Differenzenkoeffizienten (Anzahl Koeffizienten)
  !===========================================================================================================
  !*** Compact *************************************************
  !> \{ \brief Compact
  integer, parameter  ::  ndL = 2 ! 3 ! 1
  integer, parameter  ::  ndR = 3 ! 3 ! 1
  !> \}
  !*************************************************************
  
  !--- Anzahl Stencil-Koeffizienten (Rand) -------------------------------------------------------------------
  !INTEGER, PARAMETER  ::  ncb1c(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1f(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1r(2) = (/2,3/)
  !INTEGER, PARAMETER  ::  ncb1g(2) = (/0,2/)
  !INTEGER, PARAMETER  ::  ncb1d(2) = (/2,2/)
  
  !INTEGER, PARAMETER  ::  ncb1c(3) = (/3,4,5/)
  !INTEGER, PARAMETER  ::  ncb1f(3) = (/3,4,5/)
  !INTEGER, PARAMETER  ::  ncb1r(3) = (/2,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(3) = (/2,3,4/)
  !INTEGER, PARAMETER  ::  ncb1d(3) = (/3,4,4/)
  
  ! Stabil   (xi >= 2, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(4) = (/4,5,5,7/)
  !INTEGER, PARAMETER  ::  ncb1f(4) = (/4,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(4) = (/2,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(4) = (/3,4,4,6/)
  !INTEGER, PARAMETER  ::  ncb1d(4) = (/4,4,6,6/)
  
  !> \{ \brief Anzahl Stencil-Koeffizienten (Rand)
  integer, parameter  ::  ncb1c(4) = (/4,5,5,7/)
  integer, parameter  ::  ncb1f(4) = (/4,5,5,5/) ! TEST!!!dim_ncb1c wird als Array-Laenge angenommen ...
  integer, parameter  ::  ncb1r(4) = (/2,3,3,3/) ! TEST!!!dim_ncb1c wird als Array-Laenge angenommen ...
  integer, parameter  ::  ncb1g(4) = (/3,4,4,6/)
  integer, parameter  ::  ncb1d(4) = (/4,4,6,6/)

  integer, parameter  ::  ncb2c(4) = (/4,5,5,7/)
  integer, parameter  ::  ncb2f(4) = (/4,5,5,5/)
  integer, parameter  ::  ncb2r(4) = (/2,3,3,3/)
  integer, parameter  ::  ncb2g(4) = (/3,4,4,6/)
  integer, parameter  ::  ncb2d(4) = (/4,4,6,6/)

  integer, parameter  ::  ncb3c(4) = (/4,5,5,7/)
  integer, parameter  ::  ncb3f(4) = (/4,5,5,5/)
  integer, parameter  ::  ncb3r(4) = (/2,3,3,3/)
  integer, parameter  ::  ncb3g(4) = (/3,4,4,6/)
  integer, parameter  ::  ncb3d(4) = (/4,4,6,6/)
  !> \}
  
  ! Instabil (äquidistant, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(5) = (/5,6,6,7,9/)
  !INTEGER, PARAMETER  ::  ncb1f(5) = (/4,5,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(5) = (/2,3,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(5) = (/0,5,4,6,8/)
  !INTEGER, PARAMETER  ::  ncb1d(5) = (/5,4,6,8,8/)
  
  ! Instabil (äquidistant, Re=10000, N=17)
  !INTEGER, PARAMETER  ::  ncb1c(6) = (/6,7,7,7,9 ,11/)
  !INTEGER, PARAMETER  ::  ncb1f(6) = (/4,5,5,5,5 ,5 /)
  !INTEGER, PARAMETER  ::  ncb1r(6) = (/2,3,3,3,3 ,3 /)
  !INTEGER, PARAMETER  ::  ncb1g(6) = (/5,6,6,6,8 ,10/)
  !INTEGER, PARAMETER  ::  ncb1d(6) = (/6,6,6,8,10,10/)
  
  ! Stabil  (Re=10000, N=65, leicht gestreckt, explizites Forcing)
  !INTEGER, PARAMETER  ::  ncb1c(5) = (/3,7,7,7,9/)
  !INTEGER, PARAMETER  ::  ncb1f(5) = (/4,5,5,5,5/)
  !INTEGER, PARAMETER  ::  ncb1r(5) = (/2,3,3,3,3/)
  !INTEGER, PARAMETER  ::  ncb1g(5) = (/0,6,6,6,8/)
  !INTEGER, PARAMETER  ::  ncb1d(5) = (/6,6,6,8,8/)
  !> \}
  
  
end module mod_dims
