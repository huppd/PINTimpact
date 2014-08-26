#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP


#include "Teuchos_RCP.hpp"
//#include "Teuchos_Array.hpp"
#include <vector>

#include "Pimpact_Space.hpp"

#include "Pimpact_CoarsenStrategy.hpp"




namespace Pimpact {


template<class Field>
class MultiGrid {

  typedef typename Field::Scalar  Scalar;
  typedef typename Field::Ordinal Ordinal;

  static const int dimension = Field::dimension;

  typedef typename Field::SpaceT SpaceT;

  template<class FieldT, class CoarsenStrategy>
  friend
  Teuchos::RCP< MultiGrid<FieldT> >
  createMultiGrid( const Teuchos::RCP<const typename FieldT::SpaceT>& space, int nGridsMax=10 );

protected:

  int nGrids_;

  std::vector< Teuchos::RCP<const SpaceT> > multiSpace_;


  MultiGrid( const std::vector<Teuchos::RCP<const SpaceT> >& multiSpace ):
    nGrids_( multiSpace.size() ),
    multiSpace_(nGrids_) {}

public:

  int getNGrids() const {
    return( nGrids_ );
  }


}; // end of class MultiGrid


/// \relates MultiGrid
template<class Field, class CoarsenStrategy>
Teuchos::RCP< MultiGrid<Field> >
createMultiGrid( const Teuchos::RCP<const typename Field::SpaceT>& space, int maxGrids=10 ) {

  std::vector<Teuchos::RCP<const typename Field::SpaceT> > spaces =
      CoarsenStrategy::getMultiSpace( space, maxGrids );
  return( Teuchos::rcp( new MultiGrid<Field>( spaces ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
