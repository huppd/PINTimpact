#pragma once
#ifndef PIMPACT_MULTIGRID_HPP
#define PIMPACT_MULTIGRID_HPP




#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_MGFields.hpp"
#include "Pimpact_MGOperators.hpp"
#include "Pimpact_MGTransfers.hpp"
#include "Pimpact_MGSmoothers.hpp"


/// \defgroup MG MultiGrid
///
/// Multi Grid




namespace Pimpact {


/// \ingroup MG
/// \todo think about templeting everything MGSpaces,  MGOperators,...
template<
class FSpaceT,
class CSpaceT,
template<class> class FieldT,
template<class> class FOperatorT,
template<class> class COperatorT,
template<class> class SmootherT,
class CGridSolverT >
class MultiGrid {

  typedef MGSpaces<FSpaceT,CSpaceT> MGSpacesT;

  typedef MGTransfers<MGSpacesT> MGTransfersT;

  typedef MGFields<MGSpacesT,FieldT> MGFieldsT;

  typedef MGOperators<MGSpacesT,FOperatorT,COperatorT>   MGOperatorsT;

  typedef MGSmoothers<MGSpacesT, SmootherT> MGSmoothersT;

protected:

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP<const MGTransfersT> mgTrans_;

  Teuchos::RCP<const MGOperatorsT> mgOps_;

  Teuchos::RCP<const MGSmoothersT> mgSms_;

  Teuchos::RCP<const CGridSolverT> cGridSolver_;

  Teuchos::RCP<MGFieldsT> x_;
  Teuchos::RCP<MGFieldsT> b_;

  MultiGrid(
      const Teuchos::RCP<const MGSpacesT>& mgSpaces,
      const Teuchos::RCP<const CGridSolverT>& cGridSolver,
      EField type = EField::S ):
        mgSpaces_(mgSpaces),
        mgTrans_( createMGTransfers(mgSpaces) ),
        mgOps_( Pimpact::createMGOperators<FOperatorT,COperatorT>(mgSpaces) ),
        mgSms_( Pimpact::createMGOperators<FOperatorT,COperatorT>(mgSpaces) ),
        x_( createMGField(mgSpaces, type ) ),
        b_( createMGField(mgSpaces, type ) ),
        cGridSolver_(cGridSolver) {}

public:

  //  void print(  std::ostream& out=std::cout ) const {}

  //  }

}; // end of class MultiGrid



/// \relates MultiGrid
template<
template<class> class FieldT,
template<class> class FOperatorT,
template<class> class COperatorT,
template<class> class SmootherT,
class CGridSolverT,
class MGSpacesT >
Teuchos::RCP< MultiGrid<typename MGSpacesT::FSpaceT,typename MGSpacesT::CSpaceT,FieldT,FOperatorT,COperatorT,SmootherT,CGridSolverT> >
createMultiGrid(
    const Teuchos::RCP<const MGSpacesT>& mgSpaces,
    const Teuchos::RCP<const CGridSolverT>& cGridSolver,
    EField type = EField::S  ) {

  return(
      Teuchos::rcp(
          new MultiGrid<typename MGSpacesT::FSpaceT,typename MGSpacesT::CSpaceT,FieldT,FOperatorT,COperatorT,SmootherT,CGridSolverT>(
              mgSpaces, cGridSolver,type)
      )
  );

}
//
////  typedef typename FField::SpaceT FSpaceT;
//  typedef typename CField::SpaceT CSpaceT;
//
//
//  std::vector<Teuchos::RCP<const CSpaceT> > spaces =
//      CoarsenStrategy::getMultiSpace( space, maxGrids );
//
//  std::vector< Teuchos::RCP<CField> > fields;
//
//  for( unsigned i=0; i<spaces.size(); ++i )
//    fields.push_back( Teuchos::rcp( new CField(spaces[i],true,type) ) );
//
//  std::vector< Teuchos::RCP< typename MultiGrid<FField,CField>::RestrictionOpT> > restrictOp;
//  std::vector< Teuchos::RCP< typename MultiGrid<FField,CField>::InterpolationOpT> > interpolationOp;
//
//  for( unsigned i=0; i<spaces.size()-1; ++i ) {
//    restrictOp.push_back( createRestrictionOp( spaces[i], spaces[i+1] ) );
//    interpolationOp.push_back( createInterpolationOp( spaces[i+1], spaces[i] ) );
//  }
//
//  return(
//      Teuchos::rcp(
//          new MultiGrid<FField,CField>(
//              space,
//              spaces,
//              fields,
//              restrictOp,
//              interpolationOp ) ) );
//
//}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MULTIGRID_HPP
