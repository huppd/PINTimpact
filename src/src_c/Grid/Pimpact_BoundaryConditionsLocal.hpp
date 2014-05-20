#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONS_HPP
#define PIMPACT_BOUNDARYCONDITIONS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
//	void fsetBC(const int& BC_1L_, const int& BC_1U_, const int& BC_2L_, const int& BC_2U_, const int& BC_3L_, const int& BC_3U_ );
	void fgetBCloc(
	    int& BC_1L_,
	    int& BC_1U_,
	    int& BC_2L_,
	    int& BC_2U_,
	    int& BC_3L_,
	    int& BC_3U_ );
}



namespace Pimpact{

class BoundaryConditionsLocal {

public:

	typedef const Teuchos::Tuple<EBCType,3> TBC3;

//protected:

	TBC3 BCL_local_;
	TBC3 BCU_local_;

//public:

	BoundaryConditionsLocal(
	    EBCType BC1L=DirichletBC,
	    EBCType BC1U=DirichletBC,
	    EBCType BC2L=DirichletBC,
	    EBCType BC2U=DirichletBC,
	    EBCType BC3L=DirichletBC,
	    EBCType BC3U=DirichletBC ):
		BCL_local_( Teuchos::tuple(BC1L, BC2L, BC3L) ),
		BCU_local_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {};

	BoundaryConditionsLocal( TBC3 BCL_local, TBC3 BCU_local ):
		BCL_local_( BCL_local ),
		BCU_local_( BCU_local ) {};

//	void set_Impact(){ fsetBC( BCL_global_[0], BCU_global_[0], BCL_global_[1], BCU_global_[1], BCL_global_[2], BCU_global_[2]);
	};
}; // end of class BoundaryConditions





/// \relates BoundaryConditions
Teuchos::RCP<BoundaryConditionsLocal> createBoudaryConditionsGlobal(
    EDomainType dtype = Dirichelt2DChannel ) {

	switch( dtype ) {
	case AllDirichlet:
		return(
		    Teuchos::rcp(
		        new BoundaryConditionsGlobal() ) );
		break;
	case Dirichelt2DChannel:
	  return(
	      Teuchos::rcp(
	          new BoundaryConditionsGlobal(
	              DirichletBC,
	              DirichletBC,
	              DirichletBC,
	              DirichletBC,
	              PeriodicBC,
	              PeriodicBC ) ) );
		break;
	case Periodic2DChannel:
	  return(
	      Teuchos::rcp(
	          new BoundaryConditionsGlobal(
	              PeriodicBC,
	              PeriodicBC,
	              DirichletBC,
	              DirichletBC,
	              PeriodicBC,
	              PeriodicBC ) ) );
		break;
	case AllNeumann2D:
	  return(
	      Teuchos::rcp(
	          new BoundaryConditionsGlobal(
	              NeumannBC,
	              NeumannBC,
	              NeumannBC,
	              NeumannBC,
	              PeriodicBC,
	              PeriodicBC ) ) );
		break;
	case AllPeriodic:
	  return(
	      Teuchos::rcp(
	          new BoundaryConditionsGlobal(
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC ) ) );
		break;
	case Neumann1Periodic2:
	  return(
	      Teuchos::rcp(
	          new BoundaryConditionsGlobal(
	              NeumannBC,
	              NeumannBC,
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC,
	              PeriodicBC ) ) );
		break;
	default:
		std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
		return( Teuchos::null );
	}
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONS_HPP
