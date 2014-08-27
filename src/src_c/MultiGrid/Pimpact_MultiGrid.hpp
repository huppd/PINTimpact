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

  std::vector< Teuchos::RCP<const SpaceT> > multiSpace_;


  MultiGrid( const std::vector<Teuchos::RCP<const SpaceT> >& multiSpace ):
    multiSpace_(multiSpace) {}

public:

  int getNGrids() const {
    return( multiSpace_.size() );
  }

  void print(  std::ostream& out=std::cout ) const {

    for( int i = 0; i<getNGrids(); ++i ) {
      if( multiSpace_[i]->rankST()==0 ) {
        out << "-------- level: "<< i << "--------\n";
        multiSpace_[i]->print(out);
      }
//      MPI_Barrier( (*i)->commST() );
    }
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
