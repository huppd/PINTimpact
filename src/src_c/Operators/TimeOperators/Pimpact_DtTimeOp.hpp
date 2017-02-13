#pragma once
#ifndef PIMPACT_DTTIMEOP_HPP
#define PIMPACT_DTTIMEOP_HPP


#include <cmath>

#include "Pimpact_TimeField.hpp"
#include "Pimpact_Utils.hpp"
#include "Pimpact_VectorField.hpp"




namespace Pimpact{


/// \ingroup TimeOperator
template<class ST>
class DtTimeOp {

public:

  using SpaceT = ST;

  using DomainFieldT = TimeField< VectorField<SpaceT> >;
  using RangeFieldT = TimeField< VectorField<SpaceT> >;

protected:

  using Scalar = typename SpaceT::Scalar;
  using Ordinal = typename SpaceT::Ordinal;

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

    for( Ordinal i=space_->begin(F::S,3); i<space_->end(F::S,3); ++i ) {
       y.getField(i).add( mulI_, x.getField(i), -mulI_, x.getField(i-1) );
    }

    y.changed();
  }

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(space_); };

	const std::string getLabel() const { return( "DtTimeOp " ); };

}; // end of class DtTimeOp



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_DTTIMEOP_HPP
