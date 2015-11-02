#pragma once
#ifndef PIMPACT_ABSTRACTFIELD_HPP
#define PIMPACT_ABSTRACTFIELD_HPP


#include "mpi.h"

#include "Teuchos_RCP.hpp"

#include "BelosTypes.hpp"

#include "Pimpact_Space.hpp"




namespace Pimpact {

template<class SpaceT>
class AbstractField {

public:

//  typedef const SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  static const int dimension = SpaceT::dimension;

  AbstractField( const Teuchos::RCP<const SpaceT>& space ):space_(space) {};

protected:

  Teuchos::RCP<const SpaceT> space_;

  void reduceNorm( const MPI_Comm& comm, double& norm, Belos::NormType type = Belos::OneNorm ) const {

    Scalar normGlob;

    switch(type) {
    case Belos::OneNorm:
      MPI_Allreduce( &norm, &normGlob, 1, MPI_REAL8, MPI_SUM, comm );
      norm = normGlob;
      break;
    case Belos::TwoNorm:
      MPI_Allreduce( &norm, &normGlob, 1, MPI_REAL8, MPI_SUM, comm );
      norm = std::sqrt( normGlob );
      break;
    case Belos::InfNorm:
      MPI_Allreduce( &norm, &normGlob, 1, MPI_REAL8, MPI_MAX, comm );
      norm = normGlob;
      break;
    }
  }

}; // end of class AbstractField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ABSTRACTFIELD_HPP
