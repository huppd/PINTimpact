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

    for( OT k=space->si(F::S,Z); k<=space->ei(F::S,Z); ++k )
      for( OT j=space->si(F::S,Y); j<=space->ei(F::S,Y); ++j )
        for( OT i=space->si(F::S,X); i<=space->ei(F::S,X); ++i ) {
          ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Y,j) * coord->dx(F::S,Z,k);
          energy += volume * std::pow( temp(i,j,k), 2 );
        }

  }
  return( vel.allReduce( energy ) );
}


/// - compute energy local
/// - reduce to rank_dim=0
/// - having two vectors for y and e(z)
/// \f[ e(y) = \int \int (u^2 + v^2 + z^2 ) \exp( -\frac{x}{\gamma}^2/2 ) \f]
template<class SpaceT>
void computeEnergyY( const VectorField<SpaceT>& vel, std::ostream& out=std::cout, typename SpaceT::Scalar gamma=10. ) {

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  auto space = vel.space();
  auto coord = space->getCoordinatesLocal();

  ScalarField<SpaceT> temp(space);

  std::vector<ST> energyY( space->nLoc(Y), Teuchos::ScalarTraits<ST>::zero() );

  for( F f=F::U; f<SpaceT::sdim; ++f ) {

    space->getInterpolateV2S()->apply( vel(f), temp );

    for( OT k=space->si(F::S,Z); k<=space->ei(F::S,Z); ++k )
      for( OT j=space->si(F::S,Y); j<=space->ei(F::S,Y); ++j )
        for( OT i=space->si(F::S,X); i<=space->ei(F::S,X); ++i ) {
          ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Z,k);
          energyY[j-space->si(F::S,Y)] += volume * std::pow( temp(i,j,k), 2 )*std::exp( -0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2 ) );
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

  //std::cout << space->nGlo(Y) << "\n";
  //std::cout << space->nLoc(Y) << "\n";
  //std::cout << space->getProcGrid()->getNP(Y) << "\n";

  std::vector<ST> energyYglobal( space->nGlo(Y), Teuchos::ScalarTraits<ST>::zero() );

  if( 0==space->getProcGrid()->getRankSlice(Y) ) {
    MPI_Gather(
      energyY.data(),                        // void* send_data,
      space->nLoc(Y)-1,                        // int send_count,
      MPI_REAL8,                             // MPI_Datatype send_datatype,
      energyYglobal.data(),                  // void* recv_data,
      space->nLoc(Y)-1,                        // int recv_count,
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
