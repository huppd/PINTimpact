#pragma once
#ifndef PIMPACT_MGFIELDS_HPP
#define PIMPACT_MGFIELDS_HPP

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Types.hpp"




namespace Pimpact {



/// \ingroup MG
template<class MGSpacesT, template<class> class FieldT>
class MGFields {

public:

  typedef typename MGSpacesT::FSpaceT FSpaceT;
  typedef typename MGSpacesT::CSpaceT CSpaceT;

  typedef FieldT<FSpaceT> FFieldT;
  typedef FieldT<CSpaceT> CFieldT;

protected:

//  template< template<class> class FieldTT, class MGSpacesTT >
//  friend
//  Teuchos::RCP< MGFields<MGSpacesTT,FieldTT> >
//  createMGFields(
//      const Teuchos::RCP<const MGSpacesTT>& space,
//      EField type=EField::S  );

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP< FFieldT >                fField_;
  std::vector< Teuchos::RCP<CFieldT> >  cFields_;

public:

  MGFields( const Teuchos::RCP<const MGSpacesT>& mgSpaces,
      EField type=EField::S  ):
    mgSpaces_(mgSpaces),
    fField_( Teuchos::rcp( new FFieldT( mgSpaces_->get(), type ) ) ),
    cFields_( mgSpaces_->getNGrids() ) {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      cFields_[i] = Teuchos::rcp( new CFieldT( mgSpaces_->get(i), type ) );

	// not working on brutus
    //cFields_.shrink_to_fit();

  }

public:

  Teuchos::RCP<FFieldT> get(       ) const { return( fField_ ); }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses space
  Teuchos::RCP<CFieldT> get( int i ) const {
    if( i<0 )
      return( cFields_[mgSpaces_->getNGrids()+i] );
    else
      return( cFields_[i] );
  }

  //  void print(  std::ostream& out=std::cout ) const {
  //
  //  }

}; // end of class MGFields



/// \relates MGFields
template< template<class> class FieldT, class MGSpacesT >
Teuchos::RCP< MGFields<MGSpacesT,FieldT> >
createMGFields(
    const Teuchos::RCP<const MGSpacesT>& mgSpaces,
    EField type=EField::S ) {

  return(
      Teuchos::rcp(
          new MGFields<MGSpacesT,FieldT>(
              mgSpaces,
              type ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_MGFIELDS_HPP
