#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
#define PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
	void fsetBC(
	    const int& BC_1L_,
	    const int& BC_1U_,
	    const int& BC_2L_,
	    const int& BC_2U_,
	    const int& BC_3L_,
	    const int& BC_3U_ );
}



namespace Pimpact{

class BoundaryConditionsGlobal {

public:

	typedef const Teuchos::Tuple<EBCType,3> TBC3;

protected:

	TBC3 BCL_global_;
	TBC3 BCU_global_;

public:

	BoundaryConditionsGlobal(
	    EBCType BC1L=DirichletBC,
	    EBCType BC1U=DirichletBC,
	    EBCType BC2L=DirichletBC,
	    EBCType BC2U=DirichletBC,
	    EBCType BC3L=DirichletBC,
	    EBCType BC3U=DirichletBC ):
	      BCL_global_( Teuchos::tuple(BC1L, BC2L, BC3L) ),
	      BCU_global_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {

	  set_Impact();

	};


	BoundaryConditionsGlobal( TBC3 BCL_global, TBC3 BCU_global ):
		BCL_global_( BCL_global ), BCU_global_( BCU_global ) {

	  set_Impact();

	};

	EBCType getBCL( int dim ) {
	  return( BCL_global_[dim] );
	}

	EBCType getBCU( int dim ) {
	  return( BCU_global_[dim] );
	}

	void set_Impact(){
		fsetBC(
		    BCL_global_[0],
		    BCU_global_[0],
		    BCL_global_[1],
		    BCU_global_[1],
		    BCL_global_[2],
		    BCU_global_[2] );
	};

}; // end of class BoundaryConditionsGlobal





/// \relates BoundaryConditionsGlobal
Teuchos::RCP<BoundaryConditionsGlobal> createBoudaryConditionsGlobal(
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


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONSGLOBAL_HPP
