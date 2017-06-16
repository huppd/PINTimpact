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
/// \f[ e(dir) = \int \int (u^2 + v^2 + z^2 ) \exp( -\frac{x}{\gamma}^2/2 ) \f]
template<class SpaceT>
void computeEnergyDir( const VectorField<SpaceT>& vel, std::ostream& out=std::cout,
    const ECoord dir=ECoord::Y, const typename SpaceT::Scalar gamma=0. ) {

  using ST = typename SpaceT::Scalar;
  using OT = typename SpaceT::Ordinal;

  auto space = vel.space();
  auto coord = space->getCoordinatesLocal();

  ScalarField<SpaceT> temp(space);

  std::vector<ST> energy( space->nLoc(dir), Teuchos::ScalarTraits<ST>::zero() );

  for( F f=F::U; f<SpaceT::sdim; ++f ) {

    space->getInterpolateV2S()->apply( vel(f), temp );

    for( OT k=space->si(F::S,Z); k<=space->ei(F::S,Z); ++k )
      for( OT j=space->si(F::S,Y); j<=space->ei(F::S,Y); ++j )
        for( OT i=space->si(F::S,X); i<=space->ei(F::S,X); ++i ) {
          ST vel_2 = std::pow( temp(i,j,k), 2 );
          switch( dir ) {
            case X: {
              ST volume = coord->dx(F::S,Y,j) * coord->dx(F::S,Z,k);
              energy[i-space->si(F::S,X)] +=
                volume * vel_2 * ((0.==gamma)?1.:std::exp(-0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2)));
              break;
            }
            case Y: {
              ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Z,k);
              energy[j-space->si(F::S,Y)] +=
                volume * vel_2 * ((0.==gamma)?1.:std::exp(-0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2)));
              break;
            }
            case Z: {
              ST volume = coord->dx(F::S,X,i) * coord->dx(F::S,Y,j);
              energy[k-space->si(F::S,Z)] +=
                volume * vel_2 * ((0.==gamma)?1.:std::exp(-0.5*std::pow( coord->getX(F::S,Z,k)/gamma, 2)));
              break;
            }
            case T: {
              std::cout << "Warning: not implemented!!!\n";
              break;
            }
          }
        }
  }

  MPI_Reduce(
    (0==space->getProcGrid()->getRankSlice(dir))?
    MPI_IN_PLACE:energy.data(),	                // void* send_data,
    energy.data(),                              // void* recv_data,
    space->nLoc(dir),                           // int count,
    MPI_REAL8,                                  // MPI_Datatype datatype,
    MPI_SUM,                                    // MPI_Op op,
    0,                                          // int root,
    space->getProcGrid()->getCommSlice(dir) );  // MPI_Comm communicator);

  //std::cout << space->nGlo(dir) << "\n";
  //std::cout << space->nLoc(dir) << "\n";
  //std::cout << space->getProcGrid()->getNP(dir) << "\n";

  std::vector<ST> energyGlobal( space->nGlo(dir), Teuchos::ScalarTraits<ST>::zero() );

  if( 0==space->getProcGrid()->getRankSlice(dir) ) {
    MPI_Gather(
      energy.data(),                           // void* send_data,
      space->nLoc(dir)-1,                      // int send_count,
      MPI_REAL8,                               // MPI_Datatype send_datatype,
      energyGlobal.data(),                     // void* recv_data,
      space->nLoc(dir)-1,                      // int recv_count,
      MPI_REAL8,                               // MPI_Datatype recv_datatype,
      0,                                       // int root,
      space->getProcGrid()->getCommBar(dir) ); // MPI_Comm communicator);

    if( 0==space->getProcGrid()->getRankBar(dir) )
      for( OT j=0; j<space->nGlo(dir); ++j )
        out << space->getCoordinatesGlobal()->getX(F::S,dir,j+1) << "\t" << energyGlobal[j] << "\n";
  }
}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ANALYSISTOOLS_HPP
