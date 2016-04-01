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

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

  static const int dimension = SpaceT::dimension;

  AbstractField( const Teuchos::RCP<const SpaceT>& space ):space_(space) {};

protected:

  Teuchos::RCP<const SpaceT> space_;

	constexpr Scalar reduce( const MPI_Comm& comm, const double& normLocal, const MPI_Op& op=MPI_SUM ) const {

		Scalar normGlob;
		MPI_Allreduce( &normLocal, &normGlob, 1, MPI_REAL8, op, comm );
		return( normGlob );

	}

}; // end of class AbstractField


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_ABSTRACTFIELD_HPP
