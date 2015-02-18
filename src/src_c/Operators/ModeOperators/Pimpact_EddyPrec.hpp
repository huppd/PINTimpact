#pragma once
#ifndef PIMPACT_EDDYPREC_HPP
#define PIMPACT_EDDYPREC_HPP

#include "Teuchos_RCP.hpp"

#include "Pimpact_VectorField.hpp"
#include "Pimpact_MultiField.hpp"

#include "Pimpact_Operator.hpp"
#include "Pimpact_OperatorBase.hpp"




namespace Pimpact {



/// \ingroup ModeOperator
template<class OpT>
class EddyPrec {

public:

	typedef typename OpT::SpaceT SpaceT;

	typedef ModeField<typename OpT::DomainFieldT> DomainFieldT;
	typedef ModeField<typename OpT::RangeFieldT > RangeFieldT;

protected:

 	typedef typename SpaceT::Scalar Scalar;

 	Teuchos::RCP<DomainFieldT> temp_;
 	Teuchos::RCP<OpT> op_;

public:

 	EddyPrec( const Teuchos::RCP<OpT>& op ):
 		temp_( create<DomainFieldT>(op->space()) ),
 		op_(op) {};

 	void apply(const DomainFieldT& x, RangeFieldT& y) {

 		temp_->getCFieldPtr()->add( 1., x.getConstCField(),  1., x.getConstSField() );
 		temp_->getSFieldPtr()->add( 1., x.getConstCField(), -1., x.getConstSField() );

 		op_->apply( temp_->getConstCField(), y.getCField() );
 		op_->apply( temp_->getConstSField(), y.getSField() );

 		y.scale( (Scalar)0.5 );

 	}

 	void assignField( const DomainFieldT& mv ) {
//		 op_->assignField( mv- );
 	};


	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	Teuchos::RCP<OpT> getOperator() const { return(op_); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {
		op_->setParameter( para );
	}

 	bool hasApplyTranspose() const { return( false ); }

}; // end of class EddyPrec



/// \relates EddyPrec
template<class OpT>
Teuchos::RCP<EddyPrec<OpT> > createEddyPrec(
		const Teuchos::RCP<OpT>& op ) {

	return( Teuchos::rcp( new EddyPrec<OpT>( op ) ) );

}



} // end of namespace Pimpact

#endif // end of #ifndef PIMPACT_EDDYPREC_HPP
