#pragma once
#ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
#define PIMPACT_DIVGRADO2JSMOOTHER_HPP


#include "Pimpact_DivGradO2Op.hpp"




namespace Pimpact{


extern "C" {

void OP_DivGradO2JSmoother(
    const int& dimens,
    const int* const N,
    const int* const BL,
    const int* const BU,
    const int* const BCL,
    const int* const BCU,
    const double* const cdg1,
    const double* const cdg2,
    const double* const cdg3,
    const double& omega,
    const double* const b,
    const double* const x,
          double* const temp );

}



/// \brief \f$\omega\f$-Jacobian smoother for second Order DivGradOp.
///
///
/// \relates DivGradO2Op
/// \ingroup BaseOperator
/// \todo instead of hardcode 2nd Order it would be pretty to use new space with StencilWidth<3,2>
/// \todo add ParameterList
/// \todo handle corner
template<class OperatorT>
class DivGradO2JSmoother {

public:

  typedef typename OperatorT::SpaceT SpaceT;

  typedef typename SpaceT::Scalar Scalar;
  typedef typename SpaceT::Ordinal Ordinal;

  typedef ScalarField<SpaceT>  DomainFieldT;
  typedef ScalarField<SpaceT>  RangeFieldT;

protected:

  Scalar omega_;
  int nIter_;
	bool levelYes_;

  Teuchos::RCP<DomainFieldT> temp_;

  const Teuchos::RCP<const OperatorT> op_;

public:

  DivGradO2JSmoother(
      const Teuchos::RCP<const SpaceT>& space):
    omega_( (2==space->dim())?0.8:6./7. ),
    nIter_( 2 ),
    levelYes_( false ),
    temp_( createScalarField<SpaceT>(space) ),
    op_( Teuchos::rcp( new OperatorT(space) ) ) {}

  DivGradO2JSmoother(
      const Teuchos::RCP<const OperatorT>& op,
      Teuchos::RCP<Teuchos::ParameterList> pl=Teuchos::parameterList() ):
    omega_( pl->get<Scalar>("omega", (2==op->space()->dim())?0.8:6./7. ) ),
    nIter_( pl->get<int>( "numIters", 4 ) ),
    levelYes_( pl->get<bool>( "level", false ) ),
    temp_( createScalarField<SpaceT>( op->space() ) ),
    op_(op) {}


  /// \f[ y_k = (1-\omega) y_k + \omega D^{-1}( x - A y_k ) \f]
	void apply(const DomainFieldT& x, RangeFieldT& y,
			Belos::ETrans trans=Belos::NOTRANS ) const {

		for( int i=0; i<nIter_; ++i) {
      y.exchange();

			OP_DivGradO2JSmoother(
					space()->dim(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					op_->getC(X),
					op_->getC(Y),
					op_->getC(Z),
					omega_,
					x.getConstRawPtr(),
					y.getConstRawPtr(),
					temp_->getRawPtr() );

			SF_handle_corner(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					temp_->getRawPtr() );

			temp_->changed();
			temp_->exchange();

			OP_DivGradO2JSmoother(
					space()->dim(),
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					op_->getC(X),
					op_->getC(Y),
					op_->getC(Z),
					omega_,
					x.getConstRawPtr(),
					temp_->getConstRawPtr(),
					y.getRawPtr() );

			SF_handle_corner(
					space()->nLoc(),
					space()->bl(),
					space()->bu(),
					space()->getDomain()->getBCLocal()->getBCL(),
					space()->getDomain()->getBCLocal()->getBCU(),
					y.getRawPtr() );

			y.changed();
		}
		if( levelYes_ )
			y.level();

	}

  void assignField( const DomainFieldT& mv ) {};

  bool hasApplyTranspose() const { return( false ); }

	Teuchos::RCP<const SpaceT> space() const { return(op_->space()); };

	void setParameter( Teuchos::RCP<Teuchos::ParameterList> para ) {}

  void print( std::ostream& out=std::cout ) const {
    out << "--- " << getLabel() << " ---\n";
    out << "\t omega: " << omega_ << "\n";
    out << "\t numIter: " << nIter_ << "\n";
    op_->print( out );
  }

	const std::string getLabel() const { return( "DivGradO2JSmoother" ); };

}; // end of class DivGradO2JSmoother



/// \todo move somewhere better
template<template<class> class SmootherT, class OperatorT>
Teuchos::RCP< SmootherT<OperatorT> >
create(
    const Teuchos::RCP<OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl ) {

  return(
      Teuchos::rcp( new SmootherT<OperatorT>( op, pl ) ) );

}


/// \todo move somewhere better
template<class SmootherT, class OperatorT>
Teuchos::RCP< SmootherT >
create(
    const Teuchos::RCP< OperatorT>& op,
    Teuchos::RCP<Teuchos::ParameterList> pl ) {

  return(
      Teuchos::rcp( new SmootherT( op, pl ) ) );

}



} // end of namespace Pimpact


#ifdef COMPILE_ETI
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,2> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,3,4> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,2> > >;
extern template class Pimpact::DivGradO2JSmoother< Pimpact::DivGradO2Op< Pimpact::Space<double,int,4,4> > >;
#endif


#endif // end of #ifndef PIMPACT_DIVGRADO2JSMOOTHER_HPP
