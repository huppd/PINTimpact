#include <iostream>
#include "mpi.h"
#include "Belos_config.h"

#include "pimpact.hpp"


int main(int argi, char** argv ) {

	////  --- Initialize MPI ----------------------------------------------------------------------------------------
	MPI_Init(&argi,&argv);
	int rank;
	int ierr =  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	std::cout << "rank?\t" << rank << "\n";
	fset_rank(rank);
	//
	////  --- Set alarm if queue is being used ----------------------------------------------------------------------
	if(0==rank) std::cout << "finit_alarm?\n";
	finit_alarm();
	//
	////  --- Set configuration / topology --------------------------------------------------------------------------
	configuration_();
	//
	////  --- Test of input parameters ------------------------------------------------------------------------------
	ftest_parameter();
	//
	////   --- General -------------------------------------------------------------------------------------------
	finit_general();
	//
	////   --- MPI -----------------------------------------------------------------------------------------------
	finit_parallel();
	//
	////   --- Type of boundary conditions -----------------------------------------------------------------------
	finit_boundaries();
	//
	////   --- Limits of indices ---------------------------------------------------------------------------------
	finit_limits();
	//
	////  --- Physical coordinates ----------------------------------------------------------------------------------
	fcoordinates();
	//
	////  --- Determine differential coefficients -------------------------------------------------------------------
	ffd_coeffs();
	ffd_coeffs_compat();
	//
	////    --- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
	fget_stencil();
	fget_stencil_Helm();
	//
	////  --- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
	finterp_coeffs();
	finterp_coeffs_Helm();
	frestr_coeffs();
	frestr_coeffs_Helm();
	//	int rank = init_impact( argi, argv, Teuchos::null );

	int task;
	task = fget_task();
	if(0==rank) std::cout << "task?\t"<< task << "\n";
	//  --- Differenzenkoeffizienten testen -----------------------------------------------------------------------
	if( fwrite_test_yes() ) {
		if(0==rank) std::cout << "write test ja?\n";
		ftest_coeffs();
		ftest_coeffs_compact();
	}

	//  --- Weights for stop criterion ----------------------------------------------------------------------------
	fget_weights();

	//  --- Beta (for analysis) -----------------------------------------------------------------------------------
	fget_beta();
	//  ===========================================================================================================


	//  ===========================================================================================================
	//  === Main taks =============================================================================================
	//  ===========================================================================================================
	//  --- DNS / time integration --------------------------------------------------------------------------------
	if(task == 1) ftimeintegration();

	//  --- read and average fields -------------------------------------------------------------------------------
	if( task==2 ) postprocess_();

	//  --- solve eigenproblem (removed on 060513) ----------------------------------------------------------------
	//     IF (task == 3) CALL solve_eigenproblem(0)
	//   mjohn 060513 - see annotation above
	if( task==3 )
		if(rank==0)
			std::cout << "Task 3 is no longer available. Please read information concerning update on 06 May 2013.\n";

	//  --- analyze matrices --------------------------------------------------------------------------------------
	//  IF (task == 4) CALL analyze_matrix(0)
	int gridtype = 0;
	if( task==4 ) fanalyze_matrix(gridtype);
	//  ===========================================================================================================


	//  ===========================================================================================================
	if( rank==0 ) {
		bool write_stout_yes = true;
		if( write_stout_yes )
			std::cout << "\nDone ...\n";
	}

	//  ---- close HDF5 -------------------------------------------------------------------------------------------
	//CALL h5close_f(herror)
	//  h5close();


	//  ---- close MPI --------------------------------------------------------------------------------------------
	MPI_Finalize();
	//  ===============================================
}
