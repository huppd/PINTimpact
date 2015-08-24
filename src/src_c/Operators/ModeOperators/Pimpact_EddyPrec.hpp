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

	typedef typename OpT::SpaceT SpaceT;

	typedef ModeField<typename OpT::DomainFieldT> DomainFieldT;
	typedef ModeField<typename OpT::RangeFieldT > RangeFieldT;

protected:

 	typedef typename SpaceT::Scalar Scalar;

	Scalar mulI_;
	Scalar mulC_;
	Scalar mulL_;

 	Teuchos::RCP<DomainFieldT> temp_;
 	Teuchos::RCP<OpT> op_;

public:

 	EddyPrec( const Teuchos::RCP<OpT>& op ):
		mulI_(0.),
		mulC_(1.),
		mulL_( 1./op->space()->getDomain()->getDomainSize()->getRe() ),
 		temp_( create<DomainFieldT>(op->space()) ),
 		op_(op) {};


 	void apply(const DomainFieldT& x, RangeFieldT& y) {

		// left
		// where is the minus best
		temp_->getCFieldPtr()->add( 0.5, x.getConstCField(), -0.5, x.getConstSField() );
		temp_->getSFieldPtr()->add( 0.5, x.getConstCField(),  0.5, x.getConstSField() );

		// set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", mulI_ );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );

		op_->setParameter( pl );

		op_->apply( temp_->getConstCField(), y.getCField() );
		op_->apply( temp_->getConstSField(), y.getSField() );

// right
//		op_->apply( x.getConstCField(), temp_->getCField() );
//		op_->apply( x.getConstSField(), temp_->getSField() );
//
//		// where is the minus best
//		y.getCFieldPtr()->add( 0.5, temp_->getConstCField(), -0.5, temp_->getConstSField() );
//		y.getSFieldPtr()->add( 0.5, temp_->getConstCField(),  0.5, temp_->getConstSField() );
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
		mulL_ = para->get<Scalar>( "mulL", 1./space()->getDomain()->getDomainSize()->getRe() );

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
