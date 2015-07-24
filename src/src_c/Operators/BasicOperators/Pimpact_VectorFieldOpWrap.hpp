#pragma once
#ifndef PIMPACT_VECTORFIELDOPWRAP_HPP
#define PIMPACT_VECTORFIELDOPWRAP_HPP


#include "Pimpact_Types.hpp"
#include "Pimpact_ScalarField.hpp"
#include "Pimpact_VectorField.hpp"



namespace Pimpact{



/// \brief Transfers fields from "coarse" to "fine" spaces, necessary when \c Space::dimNC is  different.
///
/// Goes in both direction. If this is used a lot, it could be beneficial, to
/// seperate StencilWidths with data layout, so having same datalayout for all
/// stencile. Coping could be beneficial because Cash effects are bether
///
/// \tparam FSpaceT fine space type in the sense of stencil order.
/// \tparam CSpaceT coase space type
/// \ingroup BaseOperator
template<class SOpT>
class VectorFieldOpWrap {

public:
	

  typedef typename SOpT::FSpaceT FSpaceT;
  typedef typename SOpT::CSpaceT CSpaceT;

	typedef FSpaceT SpaceT;

//  typedef typename FSpaceT::Scalar Scalar;
//  typedef typename FSpaceT::Ordinal Ordinal;

  typedef ScalarField<typename SOpT::DomainFieldT>  DomainFieldT;
  typedef ScalarField<typename SOpT::RangeFieldT>  RangeFieldT;


protected:

  Teuchos::RCP<const SOpT> sop_;

public:

  template< class SP1T, class SP2T>
  VectorFieldOpWrap(
      const Teuchos::RCP<const SP1T>& fSpace,
      const Teuchos::RCP<const SP2T>& cSpace ):
		sop_( create<SOpT>( fSpace, cSpace ) ) {}

  template< class SP1T, class SP2T>
	VectorFieldOpWrap(
			const Teuchos::RCP<const SP1T>& fSpace,
      const Teuchos::RCP<const SP2T>& cSpace,
			Teuchos::Tuple<int,3> nb ):
		sop_( Teuchos::rcp( new SOpT( fSpace, cSpace, nb )  ) ) {}


  template< class SP1T, class SP2T>
  void apply( const VectorField<SP1T>& x, VectorField<SP2T>& y ) const {

    for( int i=0; i<x.space()->dim(); ++i ) {
      sop_->apply( x.getConstField(i), y.getField(i) );
    }

  }

//  void assignField( const RangeFieldT& mv ) {};

//  bool hasApplyTranspose() const { return( false ); }
//
	Teuchos::RCP<const SpaceT> space() const { return(sop_->space()); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		sop_->setParameter( para );
	}

  void print(  std::ostream& out=std::cout ) const {
    sop_->print();
  }

	const std::string getLabel() const { return( "VectorFieldOpWrap( "+sop_->getLabel()+" ) " ); };

}; // end of class VectorFieldOpWrap



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_VECTORFIELDOPWRAP_HPP
