#pragma once
#ifndef PIMPACT_DTTIMEOP_HPP
#define PIMPACT_DTTIMEOP_HPP


#include <cmath>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Types.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


/// \ingroup TimeOperator
template<class ST>
class DtTimeOp {

public:

  typedef ST SpaceT;

  typedef typename SpaceT::Scalar Scalar;

  typedef typename SpaceT::Ordinal Ordinal;

  typedef TimeField< VectorField<SpaceT> >  DomainFieldT;
  typedef TimeField< VectorField<SpaceT> >  RangeFieldT;

protected:

  Teuchos::RCP<const SpaceT> space_;

  Scalar mulI_;

public:

	DtTimeOp( const Teuchos::RCP<const SpaceT>& space ):
		space_(space) {

			Scalar pi = 4.*std::atan(1.);
			Scalar idt = ((Scalar)space_->nGlo()[3])/2./pi;
			mulI_ =
				space_->getDomainSize()->getAlpha2()*idt/space_->getDomainSize()->getRe();

		};


  void apply( const DomainFieldT& x, RangeFieldT& y ) const {

    x.exchange();

    for( Ordinal i=space_->sInd(S,3); i<space_->eInd(S,3); ++i ) {
       y.getFieldPtr(i)->add( mulI_, x.getConstField(i), -mulI_, x.getConstField(i-1) );
    }

    y.changed();

  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(space_); };

	const std::string getLabel() const { return( "DtTimeOp " ); };

}; // end of class DtTimeOp



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
