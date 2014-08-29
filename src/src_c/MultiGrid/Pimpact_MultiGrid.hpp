#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP


#include "Teuchos_RCP.hpp"
//#include "Teuchos_Array.hpp"
#include <vector>

#include "Pimpact_Space.hpp"

#include "Pimpact_CoarsenStrategy.hpp"

#include "Pimpact_RestrictionOp.hpp"




namespace Pimpact {


template<class Field>
class MultiGrid {

  typedef typename Field::Scalar  Scalar;
  typedef typename Field::Ordinal Ordinal;

  static const int dimension = Field::dimension;

  typedef typename Field::SpaceT SpaceT;
  typedef RestrictionOp<Scalar,Ordinal,dimension> RestrictionOpT;

  template<class FieldT, class CoarsenStrategy>
  friend
  Teuchos::RCP< MultiGrid<FieldT> >
  createMultiGrid( const Teuchos::RCP<const typename FieldT::SpaceT>& space, int nGridsMax=10 );

protected:

  std::vector< Teuchos::RCP<const SpaceT> > multiSpace_;
  std::vector< Teuchos::RCP<Field> >        multiField_;
  std::vector< Teuchos::RCP<RestrictionOpT> > restrictionOps_;


  MultiGrid(
      const std::vector<Teuchos::RCP<const SpaceT> >& multiSpace,
      const std::vector<Teuchos::RCP<Field> >& multiField,
      const std::vector<Teuchos::RCP<RestrictionOpT> >& restrictionOps ):
    multiSpace_(multiSpace),
    multiField_(multiField),
    restrictionOps_(restrictionOps) {}

public:

  int getNGrids() const {
    return( multiSpace_.size() );
  }
  Teuchos::RCP<const SpaceT>    getSpace        ( int i ) const { return( multiSpace_[i] ); }
  Teuchos::RCP<Field>           getField        ( int i ) const { return( multiField_[i] ); }
  Teuchos::RCP<RestrictionOpT>  getRestrictionOp( int i ) const { return( restrictionOps_[i] ); }

  void print(  std::ostream& out=std::cout ) const {

    if( multiSpace_[0]->rankST()==0 ) {
      for( int i = 0; i<multiSpace_.size(); ++i ) {
          out << "-------- space: "<< i << "--------\n";
          multiSpace_[i]->print(out);
      }
      for( int i = 0; i<restrictionOps_.size(); ++i ) {
          out << "-------- restrictor: "<< i << "--------\n";
          restrictionOps_[i]->print(out);
      }
    }
  }

}; // end of class MultiGrid


/// \relates MultiGrid
template<class Field, class CoarsenStrategy>
Teuchos::RCP< MultiGrid<Field> >
createMultiGrid( const Teuchos::RCP<const typename Field::SpaceT>& space, int maxGrids=10 ) {

  typedef typename Field::Scalar  S;
  typedef typename Field::Ordinal O;

  static const int d = Field::dimension;

  std::vector<Teuchos::RCP<const typename Field::SpaceT> > spaces =
      CoarsenStrategy::getMultiSpace( space, maxGrids );

  std::vector< Teuchos::RCP<Field> > fields;

  for( unsigned i=0; i<spaces.size(); ++i )
    fields.push_back( Teuchos::rcp( new Field(spaces[i]) ) );

  std::vector< Teuchos::RCP< typename MultiGrid<Field>::RestrictionOpT> > restrictOp;

  for( unsigned i=0; i<spaces.size()-1; ++i )
    restrictOp.push_back( createRestrictionOp<S,O,d>( spaces[i], spaces[i+1] ) );

  return( Teuchos::rcp( new MultiGrid<Field>( spaces, fields, restrictOp ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
