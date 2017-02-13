#pragma once
#ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
#define PIMPACT_INTERPOLATIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Utils.hpp"

#include "Pimpact_TimeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT>
class InterpolationTimeOp  {

public:

	using DomainFieldT = TimeField<typename OperatorT::DomainFieldT>;
	using RangeFieldT = TimeField<typename OperatorT::RangeFieldT>;

	using SpaceT = typename DomainFieldT::SpaceT;

	using FSpaceT = SpaceT;
	using CSpaceT = SpaceT;

	using Ordinal = typename SpaceT::Ordinal;

protected:

  Teuchos::RCP<OperatorT> op_; //  interpolation operator in space
    
public:
    
    InterpolationTimeOp(
                        const Teuchos::RCP<const SpaceT>& spaceC,
                        const Teuchos::RCP<const SpaceT>& spaceF):
    op_(Teuchos::rcp(new OperatorT(spaceC,spaceF)))
    {};

  /// \brief default apply
	void apply( const DomainFieldT& x, RangeFieldT& y) const {
      
		Ordinal d = (spaceF()->nLoc(3)) / (spaceC()->nLoc(3));
	
		x.exchange();

		for( Ordinal i=spaceC()->begin(F::S,3); i <= spaceC()->end(F::S,3); ++i ) { 
			Ordinal iF = d*( i-spaceC()->begin(F::S,3) ) + spaceF()->begin(F::S,3);
			op_->apply( x.getField( i ), y.getField( iF ) );
			if (spaceC()->nLoc(3)==1 && d>1) 
				op_->apply( x.getField(1), y.getField(2) );
		}

		if (d > 1 && spaceC()->nLoc(3)>1) { 

			for( int i=spaceF()->begin(F::S,3) + 1; i < spaceF()->end(F::S,3) ; i=i+d ) {
				y.getField(i).add( 0.5, y.getField(i-1), 0.5, y.getField(i+1) ); 
			}		
		}
	
		y.changed(); 
  }

	Teuchos::RCP<const SpaceT> spaceC() const { return( op_->spaceC() ); };
	Teuchos::RCP<const SpaceT> spaceF() const { return( op_->spaceF() ); };

	Teuchos::RCP<OperatorT> getOperatorPtr() { return( op_ ); }


}; // end of class InterpolationTimeOp

    

/// \relates
template< class OpT >
Teuchos::RCP< InterpolationTimeOp<OpT> > createInterpolationTimeOp(
    const Teuchos::RCP<const typename OpT::CSpaceT>& spaceC, const Teuchos::RCP<const typename OpT::FSpaceT>& spaceF ) {

	return( Teuchos::rcp( new InterpolationTimeOp<OpT>( spaceC,spaceF ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
