#pragma once
#ifndef PIMPACT_INTERPOLATIONTIMEOP_HPP
#define PIMPACT_INTERPOLATIONTIMEOP_HPP

#include "Teuchos_RCP.hpp"

#include <BelosTypes.hpp>

#include "Pimpact_Types.hpp"

#include "Pimpact_TimeField.hpp"



namespace Pimpact {


/// \brief Operator wrapper for \c TimeField's.
/// \ingroup TimeOperator
/// wraps and \c OperatorT and adds the functionality of handling \c TimeField.
template<class OperatorT>
class InterpolationTimeOp  {

public:

	typedef TimeField<typename OperatorT::DomainFieldT> DomainFieldT;
	typedef TimeField<typename OperatorT::RangeFieldT> RangeFieldT;

	typedef typename DomainFieldT::SpaceT SpaceT;

	typedef SpaceT FSpaceT;
	typedef SpaceT CSpaceT;

	typedef typename SpaceT::Ordinal Ordinal;

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

		for( Ordinal i=spaceC()->sInd(S,3); i <= spaceC()->eInd(S,3); ++i ) { 
			Ordinal iF = d*( i-spaceC()->sInd(S,3) ) + spaceF()->sInd(S,3);
			op_->apply( x.getConstField( i ), y.getField( iF ) );
			if (spaceC()->nLoc(3)==1 && d>1) 
				op_->apply( x.getConstField(1), y.getField(2) );
		}

		if (d > 1 && spaceC()->nLoc(3)>1) { 

			for( int i=spaceF()->sInd(S,3) + 1; i < spaceF()->eInd(S,3) ; i=i+d ) {
				y.getFieldPtr(i)->add( 0.5, y.getField(i-1), 0.5, y.getField(i+1) ); 
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
