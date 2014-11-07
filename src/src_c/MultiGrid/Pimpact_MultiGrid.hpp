#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Space.hpp"

#include "Pimpact_CoarsenStrategy.hpp"
#include "Pimpact_RestrictionOp.hpp"
#include "Pimpact_InterpolationOp.hpp"




namespace Pimpact {


template<class FField, class CField>
class MultiGrid {

  typedef typename FField::SpaceT FSpaceT;
  typedef typename CField::SpaceT CSpaceT;

  typedef typename FSpaceT::Scalar  Scalar;
  typedef typename FSpaceT::Ordinal Ordinal;

  static const int dimension = FSpaceT::dimension;

  static const int dimNCF = FSpaceT::dimNC;
  static const int dimNCC = CSpaceT::dimNC;


  typedef RestrictionOp<CSpaceT> RestrictionOpT;
  typedef InterpolationOp<CSpaceT> InterpolationOpT;

  template<class FFieldT, class CFieldT, class CoarsenStrategyT>
  friend
  Teuchos::RCP< MultiGrid<FFieldT,CFieldT> >
  createMultiGrid(
      const Teuchos::RCP<const typename FFieldT::SpaceT>& space,
      int nGridsMax=10,
      EField type=EField::S );

protected:

  Teuchos::RCP<const FSpaceT> space_;
  std::vector< Teuchos::RCP<const CSpaceT> > multiSpace_;
  std::vector< Teuchos::RCP<CField> >        multiField_;
  std::vector< Teuchos::RCP<RestrictionOpT> > restrictionOps_;
  std::vector< Teuchos::RCP<InterpolationOpT> > interpolationOps_;


  MultiGrid(
      const Teuchos::RCP<const FSpaceT>& space,
      const std::vector<Teuchos::RCP<const CSpaceT> >& multiSpace,
      const std::vector<Teuchos::RCP<CField> >& multiField,
      const std::vector<Teuchos::RCP<RestrictionOpT> >& restrictionOps,
      const std::vector<Teuchos::RCP<InterpolationOpT> >& interpolationOps ):
    space_(space),
    multiSpace_(multiSpace),
    multiField_(multiField),
    restrictionOps_(restrictionOps),
    interpolationOps_(interpolationOps)
    {}

public:

  int getNGrids() const {
    return( multiSpace_.size() );
  }
  Teuchos::RCP<const FSpaceT>     getSpace          ()        const { return( space_ ); }
  Teuchos::RCP<const CSpaceT>     getSpace          ( int i ) const { return( multiSpace_[i] ); }
  Teuchos::RCP<CField>            getField          ( int i ) const { return( multiField_[i] ); }
  Teuchos::RCP<RestrictionOpT>    getRestrictionOp  ( int i ) const { return( restrictionOps_[i] ); }
  Teuchos::RCP<InterpolationOpT>  getInterpolationOp( int i ) const { return( interpolationOps_[i] ); }

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
      for( int i = 0; i<interpolationOps_.size(); ++i ) {
          out << "-------- interpolator: "<< i << "--------\n";
          interpolationOps_[i]->print(out);
      }
    }
  }

}; // end of class MultiGrid



/// \relates MultiGrid
template<class FField, class CField, class CoarsenStrategy>
Teuchos::RCP< MultiGrid<FField,CField> >
createMultiGrid(
    const Teuchos::RCP<const typename FField::SpaceT>& space,
    int maxGrids=10,
    EField type=EField::S ) {

//  typedef typename FField::SpaceT FSpaceT;
  typedef typename CField::SpaceT CSpaceT;


  std::vector<Teuchos::RCP<const CSpaceT> > spaces =
      CoarsenStrategy::getMultiSpace( space, maxGrids );

  std::vector< Teuchos::RCP<CField> > fields;

  for( unsigned i=0; i<spaces.size(); ++i )
    fields.push_back( Teuchos::rcp( new CField(spaces[i],true,type) ) );

  std::vector< Teuchos::RCP< typename MultiGrid<FField,CField>::RestrictionOpT> > restrictOp;
  std::vector< Teuchos::RCP< typename MultiGrid<FField,CField>::InterpolationOpT> > interpolationOp;

  for( unsigned i=0; i<spaces.size()-1; ++i ) {
    restrictOp.push_back( createRestrictionOp( spaces[i], spaces[i+1] ) );
    interpolationOp.push_back( createInterpolationOp( spaces[i+1], spaces[i] ) );
  }

  return(
      Teuchos::rcp(
          new MultiGrid<FField,CField>(
              space,
              spaces,
              fields,
              restrictOp,
              interpolationOp ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
