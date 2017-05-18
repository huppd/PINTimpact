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

 	Teuchos::RCP<OpT> op_;

public:

 	EddyPrec( const Teuchos::RCP<OpT>& op ):
		mulI_(0.),
		mulC_(1.),
		mulL_( 1./op->space()->getDomainSize()->getRe() ),
 		op_(op) {};


	void apply(const DomainFieldT& x, RangeFieldT& y) {

		//if( true || mulI_<std::max( mulC_, mulL_ ) ) {

			////set paramters
			//auto pl = Teuchos::parameterList();
			//pl->set<Scalar>( "mulI", 0. );
			//pl->set<Scalar>( "mulC", mulC_ );
			//pl->set<Scalar>( "mulL", mulL_ );
			//op_->setParameter( pl );

			//op_->apply( x.getCField(), y.getCField() );
			//op_->apply( x.getSField(), y.getSField() );
		//}
		//else
		{

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

			//MultiField<typename OpT::DomainFieldT> mx( space(), 0 );
			//MultiField<typename OpT::RangeFieldT> my( space(), 0 );

			//mx.push_back( Teuchos::rcpFromRef(temp.getCField()) );
			//mx.push_back( Teuchos::rcpFromRef(temp.getSField()) );

			//my.push_back( Teuchos::rcpFromRef(y.getCField()) );
			//my.push_back( Teuchos::rcpFromRef(y.getSField()) );
			//op_->apply( mx, my );
		}

 	}


 	void assignField( const DomainFieldT& mv ) {
//		 op_->assignField( mv- );
 	};


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
