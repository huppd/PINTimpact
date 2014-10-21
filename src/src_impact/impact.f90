!       10        20        30        40        50        60        70        80        90       100       110
!        |         |         |         |         |         |         |         |         |         |         |
!2345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!=============================================================================================================
  !===========================================================================================================
     !========================================================================================================
        !=====================================================================================================
           !==================================================================================================


!*************************************************************************************************************
!* IMPACT                                                                                                    *
!* by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)                     *
!* May 2005 - Dec 2011                                                                                       *
!* by Michael John, Institute of Fluid Dynamics, ETH Zurich (john@ifd.mavt.ethz.ch)                          *
!* May 2013 - May 2013                                                                                       *
!*************************************************************************************************************


program impact
  
!use iso_c_binding

  use mod_dims
  use mod_vars
  use mod_setup
  use mod_coeffs
  use mod_lib
  use mod_inout
  use mod_test
  use mod_geometry
  use mod_timeint
  use HDF5
  
  
  implicit none
  
  INCLUDE 'mpif.h'
  
  
  !----------------------------------------------------------------------------------------------------------!
  ! Annotations / TO-DO list:                                                                                !
  !  - henniger 2011: unite the four init_* (general, parallel, boundaries, limits), already intermingled    !
  !  - mjohn 060513: moved task 3 (eigensolver) to a usr file because mod functions depended themselves on   !
  !                  usr functions, i.e. variables / fields of task 3 were not independent! May be recovered !
  !                  as part of a mod file if new fields are provided and dependence on usr functions lifted !
  !----------------------------------------------------------------------------------------------------------!
  
  
  !===========================================================================================================
  !=== Initialization ========================================================================================
  !===========================================================================================================
  !--- Initialize MPI ----------------------------------------------------------------------------------------
  call MPI_INIT(merror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)
  
  !--- Initialize HDF5 ---------------------------------------------------------------------------------------
  call h5open_f(herror)
  
  !--- Set alarm if queue is being used ----------------------------------------------------------------------
  call init_alarm
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  call configuration
  
  !--- Test of input parameters ------------------------------------------------------------------------------
  call test_parameter
  
   !--- General -------------------------------------------------------------------------------------------
   call init_general
   !--- MPI -----------------------------------------------------------------------------------------------
   call init_parallel
   !--- Type of boundary conditions -----------------------------------------------------------------------
   call init_boundaries
   !--- Limits of indices ---------------------------------------------------------------------------------
   call init_limits
  
  !--- Physical coordinates ----------------------------------------------------------------------------------
  call coordinates
  
  !--- Determine differential coefficients -------------------------------------------------------------------
  call FD_coeffs
  call FD_coeffs_compact
  
  !--- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
  call get_stencil
  call get_stencil_Helm
  
  !--- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
  call interp_coeffs
  call interp_coeffs_Helm
  call restr_coeffs
  call restr_coeffs_Helm
  
  !--- Differenzenkoeffizienten testen -----------------------------------------------------------------------
  if (write_test_yes) call test_coeffs
  if (write_test_yes) call test_coeffs_compact
  
  !--- Weights for stop criterion ----------------------------------------------------------------------------
  call get_weights
  
  !--- Beta (for analysis) -----------------------------------------------------------------------------------
  call get_beta
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Main taks =============================================================================================
  !===========================================================================================================
  !--- DNS / time integration --------------------------------------------------------------------------------
  if (task == 1) call timeintegration
  
  !--- read and average fields -------------------------------------------------------------------------------
  if (task == 2) call postprocess

  !--- solve eigenproblem (removed on 060513) ----------------------------------------------------------------
  !   IF (task == 3) CALL solve_eigenproblem(0)
  ! mjohn 060513 - see annotation above
  if (task == 3) then
   if (rank == 0) then
     write(*,*) 'Task 3 is no longer available. Please read information concerning update on 06 May 2013.'
   end if
  end if
  
  !--- analyze matrices --------------------------------------------------------------------------------------
  if (task == 4) call analyze_matrix(0)
  !===========================================================================================================
  
  
  !===========================================================================================================
  if (rank == 0) then
     if (write_stout_yes) write(*,*)
     if (write_stout_yes) write(*,*) 'DONE ...'
     if (write_stout_yes) write(*,*)
  end if
  
  !---- close HDF5 -------------------------------------------------------------------------------------------
  call h5close_f(herror)
  
  !---- close MPI --------------------------------------------------------------------------------------------
  call MPI_FINALIZE(merror)
  !===========================================================================================================
  
  
end program impact
