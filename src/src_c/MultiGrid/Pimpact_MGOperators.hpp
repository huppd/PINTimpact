#pragma once
#ifndef PIMPACT_MGOPERATORS_HPP
#define PIMPACT_MGOPERATORS_HPP

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact {


/// \ingroup MG
/// \todo add constructror from FOperator
/// \todo generalize MGContainer, specify to Fields and Operators
template<class MGST, template<class> class FOT, template<class> class COT>
class MGOperators {

public:

  typedef MGST MGSpacesT;

  typedef typename MGSpacesT::FSpaceT FSpaceT;
  typedef typename MGSpacesT::CSpaceT CSpaceT;

  typedef FOT<FSpaceT> FOperatorT;
  typedef COT<CSpaceT> COperatorT;

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

  Teuchos::RCP< FOperatorT >               fOperator_;
  std::vector< Teuchos::RCP<COperatorT> >  cOperator_;

public:

  MGOperators( const Teuchos::RCP<const MGSpacesT>& mgSpaces  ):
    mgSpaces_(mgSpaces),
    fOperator_( Teuchos::rcp( new FOperatorT( mgSpaces_->get() ) ) ),
    cOperator_( mgSpaces_->getNGrids() ) {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      cOperator_[i] = Teuchos::rcp( new COperatorT( mgSpaces_->get(i) ) );

	// not working on brutus
    //cOperator_.shrink_to_fit();

  }

public:

  Teuchos::RCP<const MGSpacesT> getMGSpaces() const { return( mgSpaces_ ); }

  Teuchos::RCP<FOperatorT> get() const { return( fOperator_ ); }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<COperatorT> get( int i ) const {
    if( i<0 )
      return( cOperator_[mgSpaces_->getNGrids()+i] );
    else
      return( cOperator_[i] );
  }

  //  void print(  std::ostream& out=std::cout ) const {
  //
  //  }

}; // end of class MGOperators



/// \relates MGOperators
template<template<class> class FOperatorT, template<class> class COperatorT=FOperatorT, class MGSpacesT >
Teuchos::RCP<const MGOperators<MGSpacesT,FOperatorT,COperatorT> >
createMGOperators(
    const Teuchos::RCP<const MGSpacesT>& mgSpaces ) {

  return(
      Teuchos::rcp(
          new MGOperators<MGSpacesT,FOperatorT,COperatorT>(
              mgSpaces ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MGOPERATORS_HPP
