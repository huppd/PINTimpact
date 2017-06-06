#pragma once
#ifndef PIMPACT_MGFIELDS_HPP
#define PIMPACT_MGFIELDS_HPP


#include <vector>

#include "Teuchos_RCP.hpp"

#include "Pimpact_Utils.hpp"




namespace Pimpact {



/// \ingroup MG
template<class MGSpacesT, template<class> class FieldT>
class MGFields {

public:

  using FSpaceT = typename MGSpacesT::FSpaceT;
  using CSpaceT = typename MGSpacesT::CSpaceT;

  using FFieldT = FieldT<FSpaceT>;
  using CFieldT = FieldT<CSpaceT>;

protected:

//  template< template<class> class FieldTT, class MGSpacesTT >
//  friend
//  Teuchos::RCP< MGFields<MGSpacesTT,FieldTT> >
//  createMGFields(
//      const Teuchos::RCP<const MGSpacesTT>& space,
//      F type=F::S  );

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  FFieldT              fField_;
  std::vector<CFieldT> cFields_;

public:

  MGFields( const Teuchos::RCP<const MGSpacesT>& mgSpaces ):
    mgSpaces_(mgSpaces),
    fField_( mgSpaces_->get() ),
    cFields_() {

    for( int i=0; i<mgSpaces_->getNGrids(); ++i )
      //			if( 0==i || mgSpaces_->participating(i-1) )
      cFields_.push_back( CFieldT( mgSpaces_->get(i) ) );

    // not working on brutus
    //cFields_.shrink_to_fit();
  }

public:

  constexpr const FFieldT& get() {
    return( fField_ );
  }
  FFieldT& get() {
    return( fField_ );
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses space
  constexpr const CFieldT& get( int i ) const {
    if( i<0 )
      return( cFields_[mgSpaces_->getNGrids()+i] );
    else
      return( cFields_[i] );
  }

  /// \brief gets ith operator, similar to python i=-1 is gets you the coarses space
  CFieldT& get( int i )  {
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
createMGFields( const Teuchos::RCP<const MGSpacesT>& mgSpaces ) {

  return( Teuchos::rcp( new MGFields<MGSpacesT,FieldT>( mgSpaces ) ) );
}



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MGFIELDS_HPP
