#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP



#include "Pimpact_MGSpaces.hpp"




namespace Pimpact {



template<
class FSpaceT,
class CSpaceT,
template<class> class FieldT,
template<class> class FOperatorT,
template<class> class COperatorT,
class SmootherT,
class CGridSolverT >
class MultiGrid {

  typedef MGSpaces<FSpaceT,CSpaceT> MGSpacesT;

  typedef MGFields<MGSpacesT,FieldT> MGFieldsT;

  typedef MGOperators<MGSpacesT, template<class> class FOT, template<class> class COT> MGOperatorsT;


//  template<class FFieldT, class CFieldT, class CoarsenStrategyT>
//  friend Teuchos::RCP< MultiGrid<FFieldT,CFieldT> >
//  createMultiGrid(
//      const Teuchos::RCP<const typename FFieldT::SpaceT>& space,
//      int nGridsMax=10,
//      EField type=EField::S );

protected:

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP<MGFieldsT> x_;
  Teuchos::RCP<MGFieldsT> b_;

  Teuchos::RCP<MGOperatorsT> mgOps_;

  Teuchos::RCP<const CGridSolverT> cGridSolver_;

  MultiGrid(
      const Teuchos::RCP<const FSpaceT>& space,
      const std::vector<Teuchos::RCP<const CSpaceT> >& multiSpace,
      const std::vector<Teuchos::RCP<CField> >& xFields,
      const std::vector<Teuchos::RCP<CField> >& bFields,
      const std::vector<Teuchos::RCP<RestrictionOpT> >& restrictionOps,
      const std::vector<Teuchos::RCP<InterpolationOpT> >& interpolationOps ):
    space_(space),
    xFields_(xFields),
    bFields_(bFields),
    restrictionOps_(restrictionOps),
    interpolationOps_(interpolationOps)
    {}

public:

  int getNGrids() const { return( multiSpace_.size() ); }

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
