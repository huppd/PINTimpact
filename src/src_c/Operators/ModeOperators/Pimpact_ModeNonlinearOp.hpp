#pragma once
#ifndef PIMPACT_MODENONLINEAROP_HPP
#define PIMPACT_MODENONLINEAROP_HPP


#include "Pimpact_ModeField.hpp"




namespace Pimpact{



/// \ingroup ModeOperator
template<class OpT>
class ModeNonlinearOp {

public:

  using SpaceT = typename OpT::SpaceT;

  using DomainFieldT = ModeField<typename OpT::DomainFieldT >;
  using RangeFieldT = ModeField<typename OpT::RangeFieldT  >;

protected:

  using Scalar = typename SpaceT::Scalar;

	Scalar mulI_;
	Scalar mulC_;
	Scalar mulL_;

	Teuchos::RCP<OpT> op_;

public:

	ModeNonlinearOp( const Teuchos::RCP<const SpaceT>& space ):
		mulI_(0.),
		mulC_(1.),
		mulL_( 1./space->getDomainSize()->getRe() ),
		op_( create<OpT>( space ) ) { };

	ModeNonlinearOp( const Teuchos::RCP<OpT>& op ):
		mulI_(0.),
		mulC_(1.),
		mulL_( 1./op->space()->getDomainSize()->getRe() ),
		op_( op ) { };


  void apply(const DomainFieldT& x, RangeFieldT& y ) const {

		// set paramters
		auto pl = Teuchos::parameterList();
		pl->set<Scalar>( "mulI", 0.    );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );
		op_->setParameter( pl );

		//std::cout << *pl;
		//typename OpT::RangeFieldT temp( space() );

		op_->apply( x.getCField(), y.getCField() );
		y.getCField().add( 1., y.getCField(),  mulI_, x.getSField(), B::N );

		op_->apply( x.getSField(), y.getSField() );
		y.getSField().add( 1., y.getSField(), -mulI_, x.getCField(), B::N );

//		pl->set<Scalar>( "mulI", mulI_ );
//		pl->set<Scalar>( "mulC", mulC_ );
//		pl->set<Scalar>( "mulL", mulL_ );
//		op_->setParameter( pl );

  }


	void assignField( const DomainFieldT& mv ) { };


	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };


	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {
		if( para->name()!="Linear Solver" ) {
			mulI_ = para->get<Scalar>( "mulI" );
			mulC_ = para->get<Scalar>( "mulC" );
			mulL_ = para->get<Scalar>( "mulL" );
		}
	}


	Teuchos::RCP< OpT > getInnerOpPtr() { return( op_ ); }


  bool hasApplyTranspose() const { return( false ); }


	const std::string getLabel() const { return( "ModeNonlinearOp_"+op_->getLabel() ); };


  void print( std::ostream& out=std::cout ) const {
		out << getLabel() << ":\n";
		op_->print( out );
  }


}; // end of class ModeNonlinearOp




/// \relates ModeNonlinearOp
template< class OpT >
Teuchos::RCP< ModeNonlinearOp<OpT> >
createModeNonlinearOp( const Teuchos::RCP<OpT>& op ) {

  return( Teuchos::rcp( new ModeNonlinearOp<OpT>( op ) ) );

}


} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_MODENONLINEAROP_HPP
