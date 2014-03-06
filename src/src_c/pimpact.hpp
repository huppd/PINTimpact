#pragma once
#ifndef PIMPACT_HPP
#define PIMPACT_HPP

#include <iostream>
#include <mpi.h>

#include "Pimpact_DomainSize.hpp"
#include "Pimpact_BoundaryConditions.hpp"

extern "C" {
  void finit_alarm();
  void configuration_();
  void ftest_parameter();
  void finit_general();
  void finit_parallel();
  void finit_boundaries();
  void finit_limits();
  void fcoordinates();
  void ffd_coeffs();
  void ffd_coeffs_compat();
  void fget_stencil();
  void fget_stencil_Helm();
  void finterp_coeffs();
  void finterp_coeffs_Helm();
  void frestr_coeffs();
  void frestr_coeffs_Helm();
  bool fwrite_test_yes();
  void ftest_coeffs();
  void ftest_coeffs_compact();
  void fget_weights();
  void fget_beta();
  int fget_task();
  void fset_rank(int&);
  void ftimeintegration();
  void postprocess_();
  void fanalyze_matrix(int&);
}


/**
 * \brief initializes impact mainly using \c usr_config function
 * \debrecaited
 * @param argi
 * @param argv
 * @param bc
 * @return rank
 */
int init_impact( int argi=0, char** argv=0, Teuchos::RCP<Pimpact::BoundaryConditions> bc=Teuchos::null ){
	//  --- Initialize MPI ----------------------------------------------------------------------------------------
	//		MPI_Init(&argi,&argv);
	int rank;
	int ierr =  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	if( ierr!=0 )
		std::cout << "ierr:\t" << ierr << "rank:\t" << rank << "\n";

	fset_rank(rank);

	//  --- Set alarm if queue is being used ----------------------------------------------------------------------
	//		if(0==rank) std::cout << "finit_alarm?\n";
	finit_alarm();

	//  --- Set configuration / topology --------------------------------------------------------------------------
	configuration_();
	if( bc!=Teuchos::null ){
		bc->set_Impact();
	}

	//  --- Test of input parameters ------------------------------------------------------------------------------
	ftest_parameter();

	//   --- General -------------------------------------------------------------------------------------------
	finit_general();

	//   --- MPI -----------------------------------------------------------------------------------------------
	finit_parallel();

	//   --- Type of boundary conditions -----------------------------------------------------------------------
	finit_boundaries();

	//   --- Limits of indices ---------------------------------------------------------------------------------
	finit_limits();

	//  --- Physical coordinates ----------------------------------------------------------------------------------
	fcoordinates();

	//  --- Determine differential coefficients -------------------------------------------------------------------
	ffd_coeffs();
	ffd_coeffs_compat();

	//    --- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
	fget_stencil();
	fget_stencil_Helm();

	//  --- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
	finterp_coeffs();
	finterp_coeffs_Helm();
	frestr_coeffs();
	frestr_coeffs_Helm();

	return rank;
}


namespace Pimpact{

/// \brief initializes impact, transition for more and more pimpact initializing
/// \return rank
int init_impact_pre(){
	//  --- Initialize MPI ----------------------------------------------------------------------------------------
	//		MPI_Init(&argi,&argv);
	int rank;
	int ierr =  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	if( ierr!=0 )
		std::cout << "ierr:\t" << ierr << "rank:\t" << rank << "\n";

	fset_rank(rank);

	//  --- Set alarm if queue is being used ----------------------------------------------------------------------
	//		if(0==rank) std::cout << "finit_alarm?\n";
	finit_alarm();

	//  --- Set configuration / topology --------------------------------------------------------------------------
	configuration_();
	return( rank );
}


void init_impact_post(){
	//  --- Test of input parameters ------------------------------------------------------------------------------
	ftest_parameter();

	//   --- General -------------------------------------------------------------------------------------------
	finit_general();

	//   --- MPI -----------------------------------------------------------------------------------------------
	finit_parallel();

	//   --- Type of boundary conditions -----------------------------------------------------------------------
	finit_boundaries();

	//   --- Limits of indices ---------------------------------------------------------------------------------
	finit_limits();

	//  --- Physical coordinates ----------------------------------------------------------------------------------
	fcoordinates();

	//  --- Determine differential coefficients -------------------------------------------------------------------
	ffd_coeffs();
	ffd_coeffs_compat();

	//    --- Get stencil of operator div(grad( )), order of convergence is 2 ---------------------------------------
	fget_stencil();
	fget_stencil_Helm();

	//  --- Get interpolation coefficients (multigrid, order of convergence is 2) ---------------------------------
	finterp_coeffs();
	finterp_coeffs_Helm();
	frestr_coeffs();
	frestr_coeffs_Helm();

}

} // endo of namespace pimpact

#endif // PIMPACT_HPP
