#pragma once
#ifndef PIMPACT_EDDYPREC_HPP
#define PIMPACT_EDDYPREC_HPP


#include "Teuchos_RCP.hpp"

#include "Pimpact_ModeField.hpp"




namespace Pimpact {



/// \ingroup ModeOperator
template<class OpT>
class EddyPrec {

public:

	using SpaceT = typename OpT::SpaceT;

	using DomainFieldT = ModeField<typename OpT::DomainFieldT>;
	using RangeFieldT = ModeField<typename OpT::RangeFieldT >;

protected:

 	using Scalar = typename SpaceT::Scalar;

	Scalar mulI_;
	Scalar mulC_;
	Scalar mulL_;

 	Teuchos::RCP<OpT> op_;

public:

 	EddyPrec( const Teuchos::RCP<OpT>& op ):
		mulI_(0.),
		mulC_(1.),
		mulL_( 1./op->space()->getDomainSize()->getRe() ),
 		op_(op) {};


 	void apply(const DomainFieldT& x, RangeFieldT& y) {

		Teuchos::RCP<DomainFieldT> temp = create<DomainFieldT>( space() );
	
		// left
		// where is the minus best
		temp->getCFieldPtr()->add( 0.5, x.getConstCField(), -0.5, x.getConstSField() );
		temp->getSFieldPtr()->add( 0.5, x.getConstCField(),  0.5, x.getConstSField() );

		// set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", mulI_ );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );

		op_->setParameter( pl );

		op_->apply( temp->getConstCField(), y.getCField() );
		op_->apply( temp->getConstSField(), y.getSField() );

// right
//		op_->apply( x.getConstCField(), temp->getCField() );
//		op_->apply( x.getConstSField(), temp->getSField() );
//
//		// where is the minus best
//		y.getCFieldPtr()->add( 0.5, temp->getConstCField(), -0.5, temp->getConstSField() );
//		y.getSFieldPtr()->add( 0.5, temp->getConstCField(),  0.5, temp->getConstSField() );
//

//		// just block
//		// set paramters
//		auto pl = Teuchos::parameterList();
//		pl->set<Scalar>( "mulI", 0. );
//		pl->set<Scalar>( "mulC", mulC_ );
//		pl->set<Scalar>( "mulL", mulL_ );
//		op_->setParameter( pl );
//
//		op_->apply( x.getConstCField(), y.getCField() );
//		op_->apply( x.getConstSField(), y.getSField() );
//
//		pl->set<Scalar>( "mulI", mulI_ );
//		pl->set<Scalar>( "mulC", mulC_ );
//		pl->set<Scalar>( "mulL", mulL_ );
//		op_->setParameter( pl );
 	}


 	void assignField( const DomainFieldT& mv ) {
//		 op_->assignField( mv- );
 	};


	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	Teuchos::RCP<OpT> getOperator() const { return(op_); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {

		mulI_ = para->get<Scalar>( "mulI", 0. );
		mulC_ = para->get<Scalar>( "mulC", 1. );
		mulL_ = para->get<Scalar>( "mulL", 1./space()->getDomainSize()->getRe() );

	}

 	bool hasApplyTranspose() const { return( false ); }

	const std::string getLabel() const { return( "EddyPrec" ); };

  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }

}; // end of class EddyPrec



} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EDDYPREC_HPP
