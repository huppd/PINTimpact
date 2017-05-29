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
	using RangeFieldT = ModeField<typename OpT::RangeFieldT>;

protected:

 	using Scalar = typename SpaceT::Scalar;

	Scalar mulI_;
	Scalar mulC_;
	Scalar mulL_;

	int type_;

 	Teuchos::RCP<OpT> op_;

public:

	EddyPrec(
			const Teuchos::RCP<OpT>& op,
			const Teuchos::RCP<Teuchos::ParameterList>& pl=Teuchos::parameterList() ):

		mulI_(0.),
		mulC_(1.),
		mulL_( 1./op->space()->getDomainSize()->getRe() ),
		type_( pl->get<int>("type", 0) ),
 		op_(op) {};

	void apply(const DomainFieldT& x, RangeFieldT& y) const {

		//std:: cout << type_ << "\n";

		switch(type_) {
			case 0: {
								y = x;
								break;
							}
			case 1: {
								applyDTinv( x, y );
								break;
							}
			case 2: {
								applyCDinv( x, y );
								break;
							}
			case 3: {
								applyELinv( x, y );
								break;
							}
			case 4: {
								applyERinv( x, y );
								break;
							}
			default: {
								 applyERinv( x, y );
								 break;
							 }
		}
 	}

	/// left/right same
	/// \todo fix BC
	void applyDTinv( const DomainFieldT& x, RangeFieldT& y ) const {
		//std::cout << "applyDTinv\n";
		y.getCField().add( 0.0,      x.getCField(), -1.0/mulI_, x.getSField(), B::N );
		y.getSField().add( 1.0/mulI_, x.getCField(), 0.0,      x.getSField(), B::N );
	}

	void applyCDinv( const DomainFieldT& x, RangeFieldT& y ) const {
		//std::cout << "applyCDinv\n";

		//set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", 0. );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );
		op_->setParameter( pl );

		op_->apply( x.getCField(), y.getCField() );
		op_->apply( x.getSField(), y.getSField() );
	}

	void applyELinv( const DomainFieldT& x, RangeFieldT& y ) const {

		//std::cout << "applyELinv\n";
		DomainFieldT temp( space() );

		// left
		temp = x;
		temp.getCField().add( 1.0, x.getCField(),  1.0, x.getSField(), B::N );
		temp.getSField().add( 1.0, x.getCField(), -1.0, x.getSField(), B::N );

		//// set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", mulI_ );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );

		op_->setParameter( pl );

		op_->apply( temp.getCField(), y.getCField() );
		op_->apply( temp.getSField(), y.getSField() );

		y.scale( 0.5, B::N );
	}

	void applyERinv( const DomainFieldT& x, RangeFieldT& y ) const {

		//std::cout << "applyERinv\n";
		DomainFieldT temp( space() );

		// right
		// set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", mulI_ );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );

		op_->setParameter( pl );

		op_->apply( x.getCField(), y.getCField() );
		op_->apply( x.getSField(), y.getSField() );

		temp = x;
		temp.getCField().add( 1.0, x.getCField(),  1.0, x.getSField(), B::N );
		temp.getSField().add( 1.0, x.getCField(), -1.0, x.getSField(), B::N );

		y.scale( 0.5, B::N );
	}


	void assignField(const DomainFieldT& mv) {};

	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };

	Teuchos::RCP<OpT> getOperator() const { return(op_); };

	void setParameter( const Teuchos::RCP<Teuchos::ParameterList>& para ) {

		if( para->name()!="Linear Solver" ) {
			mulI_ = para->get<Scalar>( "mulI" );
			mulC_ = para->get<Scalar>( "mulC" );
			mulL_ = para->get<Scalar>( "mulL" );
		}
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
