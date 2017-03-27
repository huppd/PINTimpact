#pragma once
#ifndef PIMPACT_ANALYSISTOOLS_HPP
#define PIMPACT_ANALYSISTOOLS_HPP


#include <cmath>
#include <vector>

#include "Teuchos_ScalarTraits.hpp"

#include "Pimpact_Space.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact {


template<class SpaceT>
typename SpaceT::Scalar computeEnergy( const VectorField<SpaceT>& vel ) {

	using ST = typename SpaceT::Scalar;
	using OT = typename SpaceT::Ordinal;

	auto space = vel.space();
	auto coord = space->getCoordinatesLocal();

	ScalarField<SpaceT> temp(space);

	ST energy = Teuchos::ScalarTraits<ST>::zero();

	for( F f=F::U; f<SpaceT::sdim; ++f ) {

		space->getInterpolateV2S()->apply( vel(f), temp );

		for( OT k=space->begin(F::S,Z); k<=space->end(F::S,Z); ++k )
			for( OT j=space->begin(F::S,Y); j<=space->end(F::S,Y); ++j )
				for( OT i=space->begin(F::S,X); i<=space->end(F::S,X); ++i ) {
					ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Y,j) * coord->dx(F::S,Z,k);
					energy += volume * std::pow( temp(i,j,k), 2 );
				}

	}
	return( vel.allReduce( energy ) );
}


/// - compute energy local
/// - reduce to rank_dim=0 
/// - having two vectors for y and e(z)
template<class SpaceT>
void computeEnergyY( const VectorField<SpaceT>& vel, std::ostream& out=std::cout ) {

	using ST = typename SpaceT::Scalar;
	using OT = typename SpaceT::Ordinal;


	auto space = vel.space();
	auto coord = space->getCoordinatesLocal();

	ScalarField<SpaceT> temp(space);

	std::vector<ST> energyY( space->nLoc(Y), Teuchos::ScalarTraits<ST>::zero() );

	for( F f=F::U; f<SpaceT::sdim; ++f ) {

		space->getInterpolateV2S()->apply( vel(f), temp );

		for( OT k=space->begin(F::S,Z); k<=space->end(F::S,Z); ++k )
			for( OT j=space->begin(F::S,Y); j<=space->end(F::S,Y); ++j )
				for( OT i=space->begin(F::S,X); i<=space->end(F::S,X); ++i ) {
					ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Z,k);
					energyY[j-space->begin(F::S,Y)] += volume * std::pow( temp(i,j,k), 2 );
				}

	}

	MPI_Reduce(
			(0==space->getProcGrid()->getRankSlice(Y))?
				MPI_IN_PLACE:energyY.data(),	           // void* send_data,
			energyY.data(),                            // void* recv_data,
			space->nLoc(Y),                            // int count,
			MPI_REAL8,                                 // MPI_Datatype datatype,
			MPI_SUM,                                   // MPI_Op op,
			0,                                         // int root,
			space->getProcGrid()->getCommSlice(Y) );   // MPI_Comm communicator);

	std::vector<ST> energyYglobal( space->nGlo(Y)+1, Teuchos::ScalarTraits<ST>::zero() );

	if( 0==space->getProcGrid()->getRankSlice(Y) ) {
		MPI_Gather(
				energyY.data(),                        // void* send_data,
				space->nLoc(Y),                        // int send_count,
				MPI_REAL8,                             // MPI_Datatype send_datatype,
				energyYglobal.data(),                  // void* recv_data,
				space->nLoc(Y),                        // int recv_count,
				MPI_REAL8,                             // MPI_Datatype recv_datatype,
				0,                                     // int root,
				space->getProcGrid()->getCommBar(Y) ); // MPI_Comm communicator);

		if( 0==space->getProcGrid()->getRankBar(Y) )
			for( OT j=0; j<space->nGlo(Y); ++j )
				out << space->getCoordinatesGlobal()->getX(F::S,Y,j+1) << "\t" << energyYglobal[j] << "\n";
	}
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ANALYSISTOOLS_HPP
