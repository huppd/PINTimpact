#pragma once
#ifndef PIMPACT_BOUNDARYCONDITIONS_HPP
#define PIMPACT_BOUNDARYCONDITIONS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

#include "Pimpact_Types.hpp"



extern "C" {
	void fsetBC(const int& BC_1L_, const int& BC_1U_, const int& BC_2L_, const int& BC_2U_, const int& BC_3L_, const int& BC_3U_ );
}



namespace Pimpact{

class BoundaryConditions {

public:

	typedef const Teuchos::Tuple<EBCType,3> TBC3;

protected:

	TBC3 BCL_global_;
	TBC3 BCU_global_;

public:

	BoundaryConditions( EBCType BC1L=DirichletBC, EBCType BC1U=DirichletBC, EBCType BC2L=DirichletBC, EBCType BC2U=DirichletBC, EBCType BC3L=DirichletBC, EBCType BC3U=DirichletBC ):
		BCL_global_( Teuchos::tuple(BC1L, BC2L, BC3L) ), BCU_global_( Teuchos::tuple(BC1U, BC2U, BC3U) ) {};

	BoundaryConditions( TBC3 BCL_global, TBC3 BCU_global ):
		BCL_global_( BCL_global ), BCU_global_( BCU_global ) {};

	void set_Impact(){
		fsetBC( BCL_global_[0], BCU_global_[0], BCL_global_[1], BCU_global_[1], BCL_global_[2], BCU_global_[2]);
	};
}; // end of class BoundaryConditions



/// \relates BoundaryConditions
Teuchos::RCP<BoundaryConditions> createAllDirichletBC2D() {
	return Teuchos::rcp( new BoundaryConditions( DirichletBC, DirichletBC, DirichletBC, DirichletBC, PeriodicBC, PeriodicBC ) );
}


/// \relates BoundaryConditions
Teuchos::RCP<BoundaryConditions> createPeriodicChannelBC2D() {
	return Teuchos::rcp( new BoundaryConditions( PeriodicBC, PeriodicBC, DirichletBC, DirichletBC, PeriodicBC, PeriodicBC ) );
}



/// \relates BoundaryConditions
Teuchos::RCP<BoundaryConditions> createBC( EDomainType dtype = Dirichelt2DChannel ) {
	switch( dtype ) {
	case AllDirichlet:
		return Teuchos::rcp( new BoundaryConditions() );
		break;
	case Dirichelt2DChannel:
		return createAllDirichletBC2D();
		break;
	case Periodic2DChannel:
		return createPeriodicChannelBC2D();
		break;
	case AllNeumann2D:
	  return Teuchos::rcp( new BoundaryConditions( NeumannBC, NeumannBC, NeumannBC, NeumannBC, PeriodicBC, PeriodicBC ) );
		break;
	case AllPeriodic:
	  return Teuchos::rcp( new BoundaryConditions( PeriodicBC, PeriodicBC, PeriodicBC, PeriodicBC, PeriodicBC, PeriodicBC ) );
		break;
	case Neumann1Periodic2:
	  return Teuchos::rcp( new BoundaryConditions( NeumannBC, NeumannBC, PeriodicBC, PeriodicBC, PeriodicBC, PeriodicBC ) );
		break;
	default:
		std::cout << "!!!Warning: unkown EDomainType:\t" <<dtype<<"\t!!!\n";
		return Teuchos::null;
	}
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_BOUNDARYCONDITIONS_HPP
