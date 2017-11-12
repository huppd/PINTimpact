#pragma once
#ifndef PIMPACT_MGOPERATORS_HPP
#define PIMPACT_MGOPERATORS_HPP


#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {


/// \ingroup MG
/// \note todo add constructror from FOperator
/// \note todo generalize MGContainer, specify to Fields and Operators
template<class MGST, template<class> class FOT, template<class> class COT>
class MGOperators {

public:

  using MGSpacesT = MGST;

  using FSpaceT = typename MGSpacesT::FSpaceT;
  using CSpaceT = typename MGSpacesT::CSpaceT;

  using FOperatorT = FOT<FSpaceT>;
  using COperatorT = COT<CSpaceT>;

protected:

  // why?
  //  template< template<class> class FOTT, template<class> class COTT, class MGSpacesTT >
  //  friend
  //  Teuchos::RCP<const MGOperators<MGSpacesTT,FOTT,COTT> >
  //  createMGOperators( const Teuchos::RCP<const MGSpacesTT>& space );
  //  template<template<class> class FOTT, template<class> class COTT, class MGSTT >
  //  friend Teuchos::RCP<const MGOperators<MGSTT,FOTT,COTT> >
  //  createMGOperators(
  //      const Teuchos::RCP<const MGSTT>& mgSpaces );

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP<FOperatorT>                fOperator_;
  std::vector<Teuchos::RCP<COperatorT> >  cOperator_;

public:

  MGOperators( const Teuchos::RCP<const MGSpacesT>& mgSpaces,
      const Teuchos::RCP<FOperatorT>& fOperator ):
    mgSpaces_(mgSpaces),
    fOperator_( fOperator ),
    cOperator_( mgSpaces_->getNGrids() ) {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      if( mgSpaces_->participating(i) )
        cOperator_[i] = Teuchos::rcp( new COperatorT( mgSpaces_->get(i) ) );

    // not working on brutus
    //cOperator_.shrink_to_fit();
  }

  constexpr const Teuchos::RCP<const MGSpacesT>& getMGSpaces() const {
    return mgSpaces_;
  }

  constexpr const Teuchos::RCP<FOperatorT>& get() const {
    return fOperator_;
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses space
  constexpr const Teuchos::RCP<COperatorT>& get( int i ) const {
    if( i<0 )
      return cOperator_[mgSpaces_->getNGrids()+i];
    else
      return cOperator_[i];
  }

  void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) const {

    fOperator_->setParameter( para );

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      if( mgSpaces_->participating(i) )
        cOperator_[i]->setParameter( para );
  }

  void print( std::ostream& out=std::cout ) const {

    fOperator_->print();
    for( int i=0; i<mgSpaces_->getNGrids(); ++i ) {
      std::cout << "\n\n --- level: " << i << " ---\n\n";
      get(i)->print();
    }
  }

}; // end of class MGOperators



/// \relates MGOperators
template<template<class> class FOperatorT, template<class> class COperatorT=FOperatorT, class MGSpacesT >
Teuchos::RCP<const MGOperators<MGSpacesT,FOperatorT,COperatorT> >
createMGOperators(
  const Teuchos::RCP<const MGSpacesT>& mgSpaces,
  const Teuchos::RCP<FOperatorT<typename MGSpacesT::FSpaceT> >& fOperator ) {

  return Teuchos::rcp(
      new MGOperators<MGSpacesT,FOperatorT,COperatorT>( mgSpaces, fOperator ) );
}



} // end of namespace Pimpact



#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_DivGradOp.hpp"
#include "Pimpact_DivGradO2Op.hpp"
#include "Pimpact_NonlinearOp.hpp"
#include "Pimpact_ConvectionDiffusionSOp.hpp"
template<class T> using ConvDiffOpT = Pimpact::NonlinearOp<Pimpact::ConvectionDiffusionSOp<T> >;


#endif // end of #ifndef PIMPACT_MGOPERATORS_HPP
