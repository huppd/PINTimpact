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
//      EField type=EField::S  );

  Teuchos::RCP<const MGSpacesT> mgSpaces_;

  Teuchos::RCP< FFieldT >                fField_;
  std::vector< Teuchos::RCP<CFieldT> >  cFields_;

public:

	MGFields( const Teuchos::RCP<const MGSpacesT>& mgSpaces ):
		mgSpaces_(mgSpaces),
		fField_( Teuchos::rcp( new FFieldT( mgSpaces_->get() ) ) ),
		cFields_( mgSpaces_->getNGrids() ) {

			for( int i=0; i<mgSpaces_->getNGrids(); ++i )
				//			if( 0==i || mgSpaces_->participating(i-1) )
				cFields_[i] = Teuchos::rcp( new CFieldT( mgSpaces_->get(i) ) );

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
createMGFields( const Teuchos::RCP<const MGSpacesT>& mgSpaces ) {

	return( Teuchos::rcp( new MGFields<MGSpacesT,FieldT>( mgSpaces ) ) );

}



} // end of namespace Pimpact


#ifdef COMPILE_ETI
#include "Pimpact_MGSpaces.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"
extern template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::ScalarField >;
extern template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::ScalarField >;
extern template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,3,4>, Pimpact::Space<double,int,3,2> >, Pimpact::VectorField >;
extern template class Pimpact::MGFields< Pimpact::MGSpaces< Pimpact::Space<double,int,4,4>, Pimpact::Space<double,int,4,2> >, Pimpact::VectorField >;
#endif


#endif // end of #ifndef PIMPACT_MGFIELDS_HPP
