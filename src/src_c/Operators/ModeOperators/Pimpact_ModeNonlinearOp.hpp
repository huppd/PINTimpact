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
		pl->set<Scalar>( "mulI", 0. );
		pl->set<Scalar>( "mulC", mulC_ );
		pl->set<Scalar>( "mulL", mulL_ );
		op_->setParameter( pl );

		auto temp = create<typename OpT::RangeFieldT>( space() );

		op_->apply( x.getConstCField(), *temp );
		y.getCField().add( 1., *temp,  mulI_, x.getConstSField() );

		op_->apply( x.getConstSField(), *temp );
		y.getSField().add( 1., *temp, -mulI_, x.getConstCField() );

//		pl->set<Scalar>( "mulI", mulI_ );
//		pl->set<Scalar>( "mulC", mulC_ );
//		pl->set<Scalar>( "mulL", mulL_ );
//		op_->setParameter( pl );

  }


	void assignField( const DomainFieldT& mv ) { };


	constexpr const Teuchos::RCP<const SpaceT>& space() const { return(op_->space()); };


	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {

		mulI_ = para->get<Scalar>( "mulI", 0. );
		mulC_ = para->get<Scalar>( "mulC", 1. );
		mulL_ = para->get<Scalar>( "mulL", 1./space()->getDomainSize()->getRe() );
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
