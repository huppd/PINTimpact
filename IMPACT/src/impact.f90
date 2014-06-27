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


PROGRAM impact
  
  
  USE mod_dims
  USE mod_vars
  USE mod_setup
  USE mod_coeffs
  USE mod_lib
  USE mod_inout
  USE mod_test
  USE mod_geometry
  USE mod_timeint
  USE HDF5
  
  
  IMPLICIT NONE
  
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
  CALL MPI_INIT(merror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,merror)
  
  !--- Initialize HDF5 ---------------------------------------------------------------------------------------
  CALL h5open_f(herror)
  
  !--- Set alarm if queue is being used ----------------------------------------------------------------------
  CALL init_alarm
  
  !--- Set configuration / topology --------------------------------------------------------------------------
  CALL configuration
  
  !--- Test of input parameters ------------------------------------------------------------------------------
  CALL test_parameter
  
   !--- General -------------------------------------------------------------------------------------------
   CALL init_general
   !--- MPI -----------------------------------------------------------------------------------------------
   CALL init_parallel
   !--- Type of boundary conditions -----------------------------------------------------------------------
   CALL init_boundaries
   !--- Limits of indices ---------------------------------------------------------------------------------
   CALL init_limits
  
  !--- Physical coordinates ----------------------------------------------------------------------------------
  CALL coordinates
  
  !--- Determine differential coefficients -------------------------------------------------------------------
  CALL FD_coeffs
  CALL FD_coeffs_compact
  
  !--- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
  CALL get_stencil
  CALL get_stencil_Helm
  
  !--- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
  CALL interp_coeffs
  CALL interp_coeffs_Helm
  CALL restr_coeffs
  CALL restr_coeffs_Helm
  
  !--- Differenzenkoeffizienten testen -----------------------------------------------------------------------
  IF (write_test_yes) CALL test_coeffs
  IF (write_test_yes) CALL test_coeffs_compact
  
  !--- Weights for stop criterion ----------------------------------------------------------------------------
  CALL get_weights
  
  !--- Beta (for analysis) -----------------------------------------------------------------------------------
  CALL get_beta
  !===========================================================================================================
  
  
  !===========================================================================================================
  !=== Main taks =============================================================================================
  !===========================================================================================================
  !--- DNS / time integration --------------------------------------------------------------------------------
  IF (task == 1) CALL timeintegration
  
  !--- read and average fields -------------------------------------------------------------------------------
  IF (task == 2) CALL postprocess

  !--- solve eigenproblem (removed on 060513) ----------------------------------------------------------------
  !   IF (task == 3) CALL solve_eigenproblem(0)
  ! mjohn 060513 - see annotation above
  IF (task == 3) THEN
   IF (rank == 0) THEN
     WRITE(*,*) 'Task 3 is no longer available. Please read information concerning update on 06 May 2013.'
   END IF
  END IF
  
  !--- analyze matrices --------------------------------------------------------------------------------------
  IF (task == 4) CALL analyze_matrix(0)
  !===========================================================================================================
  
  
  !===========================================================================================================
  IF (rank == 0) THEN
     IF (write_stout_yes) WRITE(*,*)
     IF (write_stout_yes) WRITE(*,*) 'DONE ...'
     IF (write_stout_yes) WRITE(*,*)
  END IF
  
  !---- close HDF5 -------------------------------------------------------------------------------------------
  CALL h5close_f(herror)
  
  !---- close MPI --------------------------------------------------------------------------------------------
  CALL MPI_FINALIZE(merror)
  !===========================================================================================================
  
  
END PROGRAM impact
